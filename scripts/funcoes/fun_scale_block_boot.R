
# FUNÇÃO ------------------------------------------------------------------

# Versão da função fun_scale_invariance, mas que estima o expoente de escala 'H'
# e inclui a possibilidade de estimar o parâmetro de offset 'theta' e os
# intervalos de confiança utilizando bootstrap por blocos

# Esta é a versão mais atual da função: 25/02/2026

# Utilizada p/ gerar os resultados:
## df_scale_block_boot.pqt; e 
## df_scale_block_boot_offset.pqt

# Argumentos 'piecewise' e 'break.duration' não implementados

fun_scale_block_boot <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                 na.accept,                 # limite de falhas por ano [0,1]
                                 min.years,                 # mínimo de anos em uma série
                                 offset = FALSE,            # incluir e otimizar parametro theta
                                 which.moment = 1:3,        # quais momentos calcular
                                 which.duration = c(1, 24), # intervalo de durações
                                 d.target = NULL,           # sequencia de durações "alvo"
                                 piecewise = FALSE,         # realizar invariância com quebra?
                                 break.duration = NULL,     # duração p/ quebra
                                 min.duration = 3,          # mínimo de durações p/ cálculo
                                 ci.boot.method = "perc",   # método p/ intervalos de confiança
                                 R.boot = 1e4,              # réplicas bootstrap
                                 block.length = 1,          # tamanho do bloco
                                 return.boot = FALSE        # retornar objeto boot
                                 ){
  
  # Pacotes
  if(!require(pacman)){
    install.packages("pacman")
    message("Instalando gerenciador de pacotes 'pacman'...")
  }
  pacman::p_load(pacman, dplyr, tidyr, purrr, pbapply, boot, parallelly)
  
  # Mensagens
  message("\nCalculando intervalos de confiança bootstrap com `boot::tsboot()` e `broom::tidy.boot()`")
  message("Parâmetro de deslocamento (theta): ", offset)
  message("Durações: ", which.duration[1], " a ", which.duration[2], " horas")
  message("Tamanho do bloco: ", block.length)
  message("Nro. réplicas bootstrap: ", R.boot)
  message("Método IC: ", ci.boot.method)
  if(!is.null(d.target)) message("Durações alvo para coeficiente de desagregação: ", paste(d.target, collapse = ", "),  " h")
  if(isTRUE(return.boot)) message("Argumento `return.boot = TRUE`: use `attr(<output name>`, 'boot_obj')` para extrair")
  
  
  # FILTRAR ESTAÇÕES -------------------------------------------------------
  
  # Selecionar estações adequadas
  data <- df.imax %>%
    filter(if_else(year == 2024, TRUE, na_prct <= na.accept),    # filtrar anos por porcentagem de falhas
           between(d, which.duration[1], which.duration[2])) %>% # filtrar durações dentro do intervalo
    group_by(gauge_code) %>%                  # agrupar por estação
    filter(n_distinct(d) >= min.duration,     # só estações com pelo menos 'min.duration'
           n_distinct(year) >= min.years) %>% # so estações com pelo menos 'min.years'
    ungroup() %>% 
    select(gauge_code, year, d, imax)
  

  # DEFINIR FUNÇÕES --------------------------------------------------------

  # Função p/ otimização (minimizar model.error.root)
  # O primeiro argumento da função é o 'parâmetro' a ser otimizado
  fun.obj <- function(offset, data.matrix, d, mom){
    
    # if(offset < 0 | offset > 10) return(1e9)
    
    # Calcular momentos
    moments <- apply(X = data.matrix, MARGIN = 2, function(imax) mean(imax^mom, na.rm = TRUE))
    
    # Adicionar offset theta
    offset.durations <- d + offset
    
    # Regressão linear
    fit <- lm(log(moments) ~ log(offset.durations))
    scale <- -coef(fit)[[2]]/mom
    
    var.model <- (sigma(fit)/mom)^2
    
    return(var.model)
    
  }
  
  # Função que ajusta o modelo linear c/ um offset e
  # calcula os coeficientes de desagregação
  fun.lm <- function(offset = 0, data.matrix, d, mom){
    
    # Ajustar modelo c/ theta otimizado
    moments <- apply(X = data.matrix, MARGIN = 2, function(imax) mean(imax^mom, na.rm = TRUE))
    
    # Adicionar offset theta
    offset.durations <- d + offset
    
    # Regressão linear - c/ offset
    fit <- lm(log(moments) ~ log(offset.durations))
    scale <- -coef(fit)[[2]]/mom
    sigma <- sigma(fit)/mom
    
    # Calcular coeficientes caso definidas durações alvo
    if(is.null(d.target)){
      
      res <- c(scale, offset, sigma)
      names(res) <- c(paste0("mom", mom, "_", c("scale", "offset", "sigma")))
      
      return(res)
      
    } else{
      
      # Coeficientes de desagregação
      D <- max(d)
      disag.coef <- ((d.target + offset)/(D + offset))^(-scale)*(d.target/D)
      res <- c(scale, offset, sigma, disag.coef)
      names(res) <- c(paste0("mom", mom, "_", c("scale", "offset", "sigma",
                                                paste0("coef_d_", d.target))))
      return(res)
      
    }
    
  }
  
  # Cálculo do expoente de escala e coeficientes de desagregação
  fun.lm.boot <- function(offset = TRUE, data.matrix, d, which.moment){
    
    # Separar por duração e ajustar o modeleo linear individualmente por momentos
    regression <- lapply(X = which.moment, FUN = function(mom){
      
      # Conferir se calcula com ou sem offset
      if(isTRUE(offset)){
        
        # Otimização
        optimal <- optimize(f = fun.obj,
                            interval = c(0, max(d)),
                            data.matrix = data.matrix,
                            d = d,
                            mom = mom)
        
        optimal.offset <- optimal$minimum
        
        # Regressão linear
        res <- fun.lm(offset = optimal.offset, data.matrix = data.matrix, d = d, mom = mom)
        
        # Testar expoente de escala
        if(res[[1]] < 0 | res[[1]] > 1) res <- fun.lm(offset = 0, data.matrix = data.matrix, d = d, mom = mom)
        
        return(res)
        
      } else{
        
        # Rodar função s/ offset
        return(fun.lm(offset = 0, data.matrix = data.matrix, d = d, mom = mom))
        
      }
      
    })
    
    return(unlist(regression))
    
  }
  
  # Rodar o bootstrap e extrair intervalos de confiança
  fun.run.boot <- function(data, regional = FALSE){
    
    # Análise local ou regional
    if(isFALSE(regional)){
      
      # Transformar em matriz wide
      data.matrix <- data %>% 
        pivot_wider(names_from = d,          # cada duração em uma coluna
                    values_from = imax,      # contendo valores de imax
                    names_prefix = "d_") %>% # prefixo p/ colunas
        arrange(year) %>% 
        select(-c(gauge_code, year)) %>% 
        as.matrix()
      
      # Extrair durações
      durations <- as.numeric(stringr::str_remove_all(colnames(data.matrix), "d_"))
      
    } else{
      
      # Transformar em matriz wide
      data.matrix <- data %>% 
        unite("col_id", gauge_code, d, sep = "__d__") %>%        # juntar gauge_code e duração em uma única coluna
        pivot_wider(names_from = col_id, values_from = imax) %>% # transformar em formato longo
        arrange(year) %>% 
        select(-year) %>% 
        as.matrix()
      
      # Extrair durações
      durations <- as.numeric(stringr::str_extract(colnames(data.matrix), "(?<=__d__).*"))
      
    }
    
    # Processamento paralelo
    which.parallel.mode <- if(.Platform$OS.type == "windows") "snow" else "multicore"
    n.cores <- parallelly::availableCores(omit = 1)
    
    # Calcular bootstrap p/ série temporal
    boot.obj <- boot::tsboot(tseries = data.matrix,
                             statistic = function(x){
                               fun.lm.boot(data.matrix = x,
                                           d = durations,
                                           which.moment = which.moment,
                                           offset = offset)},
                             R = R.boot,
                             l = block.length,
                             sim = "fixed",
                             parallel = which.parallel.mode,
                             ncpus = n.cores,
                             cl = NULL)
    
    invisible(gc())
    
    # Estimar intervalos de confiança
    names.statistic <- names(boot.obj$t0)
    ls.ci <- lapply(X = seq_along(boot.obj$t0), FUN = function(i){

      # Extrair réplicas
      replicates <- boot.obj$t[,i]

      # Conferir se há alguma NA
      if(sd(replicates, na.rm = TRUE) == 0 || all(is.na(replicates))){

        ci.lower <- NA_real_
        ci.upper <- NA_real_

      }else{

        ci <- tryCatch(boot::boot.ci(boot.obj, conf = 0.95, index = i, type = ci.boot.method),
                       error = function(e) message("Erro ao processar estação", data$gauge_code, ":\n", e))

        # Conferir nome do método e corrigir se necessário
        ci.boot.method <- if_else(ci.boot.method == "perc", "percent", ci.boot.method)

        # Extrair resultados se não forem nulos
        if(!is.null(ci) && !is.null(ci[[ci.boot.method]])){

          ci.values <- ci[[ci.boot.method]]
          ci.lower <- ci.values[4]
          ci.upper <- ci.values[5]

        }else{

          ci.lower <- NA_real_
          ci.upper <- NA_real_

        }

      }

      res <- data.frame(term = names.statistic[i],
                        statistic = boot.obj$t0[i],
                        bias = mean(boot.obj$t[,i] - boot.obj$t0[i]),
                        std_error = sd(boot.obj$t[,i]),
                        ci_lower = ci.lower,
                        ci_upper = ci.upper)

      return(res)

    })
    
    # Organizar resultados
    boot.ci <- bind_rows(ls.ci) %>% 
      mutate(d = paste(which.duration, collapse = "-"),       # intervalo de duração
             n_year = n_distinct(data$year),                  # nro. anos da estação
             gauge_code = if_else(isFALSE(regional),          # se análise regional ou não
                                  as.character(first(data$gauge_code)), # extrair nome estação
                                  "Regional")) %>%            # atribui gauge_code == "Regional"
      mutate(which_moment = as.numeric(stringr::str_extract(term, "(?<=mom)\\d+")),       # extrair numérico
             which_variable = stringr::str_match(term, "^mom\\d+_([^_]+)")[,2],           # extrair scale' / 'coef' / 'offset' / 'sigma'
             d_target = if_else(which_variable == "coef",                                 # se 'coef'
                                as.numeric(stringr::str_extract(term, "(?<=_d_)[0-9.]+")), # extrair duração alvo
                                NA_real_)) %>% 
      select(gauge_code, n_year, d, which_moment, which_variable, d_target,
             statistic, bias, std_error, ci_lower, ci_upper)
    
    rownames(boot.ci) <- NULL
    
    if(isTRUE(return.boot)) attr(boot.ci, "boot_obj") <- boot.obj
    
    return(boot.ci)
    
  }
  
  # APLICAR FUNÇÕES --------------------------------------------------------
  
  # Separar por estação e calcular localmente
  ls.data <- split(x = data, f = data$gauge_code)
  ls.boot <- pbapply::pblapply(X = ls.data, FUN = fun.run.boot)
  df.boot <- bind_rows(ls.boot)
  
  # Calcular modelo regional
  message("\nSimulando bootstrap regional: ", n_distinct(data$gauge_code), " estações")
  df.boot.regional <- fun.run.boot(data, regional = TRUE)
  
  # Transformar em tbl_df
  res <- bind_rows(df.boot, df.boot.regional)
  
  message("Concluído!\n")
  return(res)
  
}