# Versão simplificada da função fun_scale_invariance, mas que calcula os intervalos
# de confiança do expoente de escala 'H' através do bootstrap

fun_scale_invariance_boot <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                      na.accept,                 # limite de falhas por ano [0,1]
                                      min.years,                 # mínimo de anos em uma série
                                      which.moment = 1:3,        # quais momentos calcular
                                      which.duration = c(1, 24), # intervalo de durações
                                      min.duration = 3,          # mínimo de durações p/ cálculo
                                      R.boot = 1e4,              # réplicas bootstrap
                                      cl = NULL){
  
  # Pacotes
  if(!require(pacman)){
    install.packages("pacman")
    message("Instalando gerenciador de pacotes 'pacman'...")
  }
  pacman::p_load(pacman, dplyr, tidyr, purrr, pbapply, broom, boot)
  
  # MOMENTOS --------------------------------------------------------------
  
  # Selecionar estações adequadas
  data <- df.imax %>%
    filter(na_prct <= na.accept,              # filtrar anos por porcentagem de falhas
           between(d, which.duration[1], which.duration[2])) %>% # filtrar durações dentro do intervalo
    group_by(gauge_code) %>%                  # agrupar por estação
    filter(n_distinct(d) >= min.duration,     # só estações com pelo menos 'min.duration'
           n_distinct(year) >= min.years) %>% # so estações com pelo menos 'min.years'
    ungroup()                                 # desagrupar
  
  # REGRESSÃO -------------------------------------------------------------
  
  # Função de regressão c/ bootstrap
  fun.lm.boot <- function(data, idx){
    
    resample <- data[idx,]
    
    # Calcular 'which.moments' momentos centrais locais (p/ cada estação)
    df.mom <- resample %>%
      group_by(d) %>% 
      reframe(n_years = n(),
              purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE)))) %>% 
      pivot_longer(cols = -c(d, n_years),      # manter colunas como estão
                   names_to = c(".value", "which_moment"), # transformar colunas separadas dos momentos
                   names_pattern = "(mom)(\\d+)") %>%      # em duas colunas contendo o valor e qual a ordem
      mutate(which_moment = as.numeric(which_moment),      # transformar ordem em numérico
             log_mom = log(mom),                           # calcular logaritmo do momento
             log_d = log(d))                               # calcular logaritmo da duração [h]
    
    # Separar estação por momentos e ajustar o modelo linear individualmente
    gauge.mom <- split(x = df.mom, f = df.mom$which_moment)
    # gg <- gauge.mom[[2]]
    regression <- lapply(X = gauge.mom, FUN = function(gg){
      
      model <- lm(log_mom~log_d, data = gg)  # ajustar modelo
      coef <- coef(model)[["log_d"]]         # extrair 'slope'
      scale <- -coef/unique(gg$which_moment) # extrair coeficiente de escala
      
    }) # fim 'regression'
    
    # Extrair vetor com coeficientes calculados
    scales <- unlist(regression)
    
  } # fim 'fun.lm.boot'
  
  # APLICAR FUNÇÃO --------------------------------------------------------
  
  # Separar por estação e calcular intervalos de confiança bootstrap
  ls.data <- split(x = data, f = data$gauge_code)
  ls.boot <- pbapply::pblapply(X = ls.data, FUN = function(gauge){
    
    boot.obj <- boot::boot(data = gauge,            # série de 1 estação
                           statistic = fun.lm.boot, # função a ser calculada
                           R = R.boot,              # nro. replicas
                           strata = gauge$d)        # reamostrar por durações
    
    boot.ci <- broom::tidy(boot.obj, conf.int = TRUE) %>% # calcular IC95 'perc'
      mutate(d = first(gauge$d),                   # criar coluna duração
             n_year = n_distinct(gauge$year),      # criar coluna nro. anos
             which_moment = row_number(),          # criar coluna momentos
             gauge_code = first(gauge$gauge_code)) # criar coluna nome estação
    
  }) # fim 'ls.boot'
  
  # Transformar em tbl_df
  df.scale.boot <- bind_rows(ls.boot)
  
  return(df.scale.boot)
  
}