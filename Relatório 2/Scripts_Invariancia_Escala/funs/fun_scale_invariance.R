
# FUNCAO SIMPLE SCALING DIRETA -------------------------------------------

# Essa função avalia invariância de escala para os 1º e 2º momentos e
# gera gráficos para análise
fit.simple.scaling <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                               which.moment = 1:3,        # quais momentos calcular
                               which.duration = c(1, 24), # intervalo de durações
                               min.duration = 3           # nro. mínimo durações p/ regressão
){
  
  
  # Mensagens
  message("Avaliacao empirica da hipotese de invariancia de escala: sentido amplo")
  message(">> Estacoes: ", n_distinct(df.imax$gauge_code))
  message(">> Duracoes: ", which.duration[1], " a ", which.duration[2], " horas")
  message(">> Momentos nao centrais: ", paste0(which.moment, sep = " "))
  
  # Função de regressão
  fun.lm.regression <- function(gg, regional.name = NULL){
    
    which.moment <- gg$which_moment[1] # extrair qual momento
    gg$log_mom <- log(gg$mom)          # logaritmo do momento
    gg$log_d <- log(gg$d)              # logaritmo da duração
    model <- lm(log_mom ~ log_d, data = gg) # regressão linear
    coef <- coef(model)[["log_d"]]     # slope
    scale <- -coef/which.moment        # expoente de escala
    summary <- summary(model)          # resumo
    rsquared <- summary$r.squared      # R²
    se.scale <- summary$coefficients["log_d","Std. Error"]/which.moment # erro padrão expoente de escala
    model.error.root <- summary$sigma  # raiz da variância do erro do modelo
    
    res <- data.frame("gauge_code" = ifelse(is.null(regional.name), gg$gauge_code[1], regional.name),
                      "n_year" = gg$n_year,
                      "d" = gg$d,
                      "which_moment" = factor(which.moment),
                      "mom" = gg$mom,
                      "coef_reg" = coef,
                      "scale" = scale,
                      "rsquared" = rsquared,
                      "std_error" = se.scale,
                      "model_error" = model.error.root)
    
  }
  
  # MOMENTOS --------------------------------------------------------------
  
  # Selecionar estações adequadas
  data <- df.imax |> filter(between(d, which.duration[1], which.duration[2]))
  
  # Calcular 'which.moment' momentos centrais locais (p/ cada estação)
  message(rep("-", getOption("width")))
  message("Estimando momentos nao centrais...")
  df.mom <- data |> 
    group_by(gauge_code, d) |> 
    summarise(
      n_year = n(), 
      purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE)))
    ) |> 
    ungroup() |> 
    pivot_longer(
      cols = -c(gauge_code, d, n_year), 
      names_to = c(".value", "which_moment"), 
      names_pattern = "(mom)(\\d+)"
    ) |> 
    mutate(which_moment = as.numeric(which_moment))
  
  # REGRESSÃO -------------------------------------------------------------
  
  # Regressão entre momentos e duração de cada estação
  ls.gauge <- split(x = df.mom, f = df.mom$gauge_code) # dividir em lista

  message(rep("-", getOption("width")))
  message("Ajustando modelo de regressao linear...")
  regression.local <- pbapply::pblapply(X = ls.gauge, FUN = function(gauge){
    
    # Aplicar para cada 'which.moment' individualmente
    ls.mom <- split(x = gauge, f = gauge$which_moment)
    regression <- lapply(X = ls.mom, FUN = fun.lm.regression)
    regression <- bind_rows(regression)
    
  })
  regression.local <- bind_rows(regression.local)
  
  # Regressão entre momentos e duração regional
  ls.gauge.regional <- split(x = df.mom, f = df.mom$which_moment)
  regression.regional <- lapply(X = ls.gauge.regional, FUN = fun.lm.regression, regional.name = "Regional")
  regression.regional <- bind_rows(regression.regional)
  regional.summary <- regression.regional |> 
    group_by(gauge_code, which_moment) |> 
    reframe(across(-c(n_year, mom), first), # todas observações são iguais, manter só 1ª
            n_year = n())           # nro. observações usadas p/ ajustar modelo lm
  
  
  # RESULTADOS ------------------------------------------------------------
  
  # Dataframe com coeficientes de escala
  scale.coefficient <- bind_rows(regression.local, regression.regional) |> 
    mutate(ci_lower = scale - 1.96*std_error,      # criar intervalos de confiança aproximados
           ci_upper = scale + 1.96*std_error) |>  # baseados no erro padrão do expoente de escala
    group_by(gauge_code, which_moment) |> 
    reframe(n_year = first(n_year),
            scale = first(scale),
            rsquared = first(rsquared),
            std_error = first(std_error),
            ci_lower = first(ci_lower),
            ci_upper = first(ci_upper)) |> 
    mutate(d = paste(which.duration, collapse = "-"))
  
  # Organizar resultados
  res <- list("scale.moments" = regression.local[c("gauge_code", "n_year", "d", "which_moment", "mom")], # contém só os momentos de cada duração
              "scale.coefficient" = scale.coefficient) # contém coeficiente de escala e informações sobre a regressão
  
  return(res)
  
}



# FUNCAO SIMPLE SCALING BOOTSTRAP ----------------------------------------


# Versão da função fun_scale_invariance, mas que estima o expoente de escala 'H'
# e inclui a possibilidade de estimar o parâmetro de offset 'theta' e os
# intervalos de confiança utilizando bootstrap por blocos

# Esta é a versão mais atual da função: 25/02/2026

# Utilizada p/ gerar os resultados:
## df_scale_block_boot.pqt; e 
## df_scale_block_boot_offset.pqt

# Argumentos 'piecewise' e 'break.duration' não implementados
simple.scaling.boot <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                offset = FALSE,            # incluir e otimizar parametro theta
                                which.moment = 1:3,        # quais momentos calcular
                                which.duration = c(1, 24), # intervalo de durações
                                d.target = NULL,           # sequencia de durações "alvo"
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
  message("Calculando intervalos de confiança bootstrap por blocos")
  message(">> Parâmetro de deslocamento (theta): ", offset)
  message(">> Durações: ", which.duration[1], " a ", which.duration[2], " horas")
  message(">> Tamanho do bloco: ", block.length)
  message(">> Nro. réplicas bootstrap: ", R.boot)
  message(">> Método IC: ", ci.boot.method)
  if(!is.null(d.target)) message(">> Durações alvo para coeficiente de desagregação: ", paste(d.target, collapse = ", "),  " h")
  if(isTRUE(return.boot)) message(">> Argumento `return.boot = TRUE`: use `attr(<output name>`, 'boot_obj')` para extrair")
  
  
  # FILTRAR ESTAÇÕES -------------------------------------------------------
  
  # Selecionar estações adequadas
  data <- df.imax |> filter(between(d, which.duration[1], which.duration[2]))
  
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
      data.matrix <- data |> 
        pivot_wider(names_from = d,          # cada duração em uma coluna
                    values_from = imax,      # contendo valores de imax
                    names_prefix = "d_") |> # prefixo p/ colunas
        arrange(wateryear) |> 
        select(-c(gauge_code, wateryear)) |> 
        as.matrix()
      
      # Extrair durações
      durations <- as.numeric(stringr::str_remove_all(colnames(data.matrix), "d_"))
      
    } else{
      
      # Transformar em matriz wide
      data.matrix <- data |> 
        unite("col_id", gauge_code, d, sep = "__d__") |>        # juntar gauge_code e duração em uma única coluna
        pivot_wider(names_from = col_id, values_from = imax) |> # transformar em formato longo
        arrange(wateryear) |> 
        select(-wateryear) |> 
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
    boot.ci <- bind_rows(ls.ci) |> 
      mutate(d = paste(which.duration, collapse = "-"),       # intervalo de duração
             n_year = n_distinct(data$wateryear),                  # nro. anos da estação
             gauge_code = if_else(isFALSE(regional),          # se análise regional ou não
                                  as.character(first(data$gauge_code)), # extrair nome estação
                                  "Regional")) |>            # atribui gauge_code == "Regional"
      mutate(which_moment = as.numeric(stringr::str_extract(term, "(?<=mom)\\d+")),       # extrair numérico
             which_variable = stringr::str_match(term, "^mom\\d+_([^_]+)")[,2],           # extrair scale' / 'coef' / 'offset' / 'sigma'
             d_target = if_else(which_variable == "coef",                                 # se 'coef'
                                as.numeric(stringr::str_extract(term, "(?<=_d_)[0-9.]+")), # extrair duração alvo
                                NA_real_)) |> 
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
