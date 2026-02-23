
# FUNÇÃO ------------------------------------------------------------------

# Versão da função fun_scale_invariance, mas que calcula os intervalos
# de confiança do expoente de escala 'H' através do bootstrap por blocos

# Posteriormente incluir argumentos 'piecewise' e 'break.duration'
fun_scale_block_boot <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                 na.accept,                 # limite de falhas por ano [0,1]
                                 min.years,                 # mínimo de anos em uma série
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
  pacman::p_load(pacman, dplyr, tidyr, purrr, pbapply, broom, boot)
  
  # Mensagens
  message("\nCalculando intervalos de confiança bootstrap com `boot::tsboot()` e `broom::tidy.boot()`")
  message("Durações: ", which.duration[1], " a ", which.duration[2], " horas")
  message("Tamanho do bloco: ", block.length)
  message("Nro. réplicas bootstrap: ", R.boot)
  message("Método IC: ", ci.boot.method)
  if(!is.null(d.target)) message("Durações alvo para coeficiente de desagregação: ", paste(d.target, collapse = ", "),  " h")
  if((isTRUE(return.boot))) message("Argumento `return.boot = TRUE`: use `attr(<output name>`, 'boot_obj')` para extrair")
  
  
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
  
  # Cálculo do expoente de escala e coeficientes de desagregação
  fun.lm.boot <- function(data.matrix, durations){
    
    # Separar por duração e ajustar o modelo linear individualmente
    regression <- lapply(X = which.moment, FUN = function(mom){
      
      # Calcular momentos não centrais p/ cada duração (coluna)
      moments <- apply(X = data.matrix, # matriz c/ linhas (anos) reamostradas
                       MARGIN = 2,      # 2: aplicar função às colunas
                       function(imax) mean(imax^mom, na.rm = TRUE))
      
      # Regressão linear
      fit <- lm(log(moments) ~ log(durations))
      slope <- coef(fit)[2]
      scale <- -slope/mom
      
      if(is.null(d.target)){
        
        res <- scale
        names(res) <- paste0("mom", mom, "_scale")
        
      } else{
        
        # Coeficientes de desagregação
        D <- max(durations)
        disag.coef <- (d.target/D)^(1 - scale)
        res <- c(scale, disag.coef)
        names(res) <- c(paste0("mom", mom, "_scale"),
                        paste0("mom", mom, "_coef_d", d.target))
        
      }
      
      return(res)
      
    })
    
    return(unlist(regression)) # vetor nomeado
    
  } # fim 'fun.lm.boot'
  
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
    
    # Calcular bootstrap p/ série temporal
    boot.obj <- boot::tsboot(tseries = data.matrix,
                             statistic = function(x) fun.lm.boot(x, durations),
                             R = R.boot,
                             l = block.length,
                             sim = "fixed")
    
    # Estimar intervalos de confiança
    boot.ci <- broom::tidy(boot.obj,                          # objeto 'boot'
                           conf.int = TRUE,                   # calcular ICs
                           conf.method = ci.boot.method) %>%  # "perc", "norm", "basic"
      mutate(d = paste(which.duration, collapse = "-"),       # intervalo de duração
             n_year = n_distinct(data$year),                  # nro. anos da estação
             gauge_code = if_else(isFALSE(regional),          # se análise regional ou não
                                  as.character(first(data$gauge_code)), # extrair nome estação
                                  "Regional")) %>%            # atribui gauge_code == "Regional"
      rename(ci_lower = conf.low, ci_upper =  conf.high, std_error = std.error) %>% 
      mutate(which_moment = as.numeric(stringr::str_extract(term, "(?<=mom)\\d+")),         # extrair numérico
             which_variable = if_else(stringr::str_detect(term, "scale"), "scale", "coef"), # extrair 'scale' ou 'coef'
             d_target = if_else(which_variable == "coef",                                   # se 'coef'
                                as.numeric(stringr::str_extract(term, "(?<=_d)[0-9.]+")),   # extrair duração alvo
                                NA_real_)) %>% 
      select(gauge_code, n_year, d, which_moment, which_variable, d_target,
             statistic, bias, std_error, ci_lower, ci_upper)
    
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
  
  message("Concluído!")
  return(res)
  
}