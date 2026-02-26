
rm(list = ls()); invisible(gc())
set.seed(1234)

# FUNÇÃO ------------------------------------------------------------------

# Versão simplificada da função fun_scale_invariance, mas que calcula os intervalos
# de confiança do expoente de escala 'H' através do bootstrap

# A lógica dessa abordagem pode estar errada, pois não preserva a correlação temporal
# anual entre os máximos, já que realiza uma amostragem idependente p/ cada duração

fun_scale_strat_boot <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                 na.accept,                 # limite de falhas por ano [0,1]
                                 min.years,                 # mínimo de anos em uma série
                                 which.moment = 1:3,        # quais momentos calcular
                                 which.duration = c(1, 24), # intervalo de durações
                                 disagg.seq = NULL,         # sequencia de durações "alvo"
                                 piecewise = FALSE,         # realizar invariância com quebra?
                                 break.duration = NULL,     # duração p/ quebra
                                 regional = TRUE,           # calcular modelo regional
                                 plot = FALSE,
                                 min.duration = 3,          # mínimo de durações p/ cálculo
                                 ci.boot.method = "bca",    # método p/ intervalos de confiança
                                 R.boot = 1e4               # réplicas bootstrap
                                 ){
  
  # Pacotes
  if(!require(pacman)){
    install.packages("pacman")
    message("Instalando gerenciador de pacotes 'pacman'...")
  }
  pacman::p_load(pacman, dplyr, tidyr, purrr, pbapply, broom, boot)
  
  # Mensagens
  message("Calculando intervalos de confiança bootstrap...\n")
  message("Durações: ", which.duration[1], " a ", which.duration[2], " horas")
  message("Nro. réplicas bootstrap: ", R.boot)
  message("Método IC: ", ci.boot.method)
  message("Coeficientes de desagregação: ", ifelse(is.null(disagg.seq), FALSE, disagg.seq),  "\n")
  
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
      reframe(n_years = n(), purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE)))) %>% 
      pivot_longer(cols = -c(d, n_years),      # manter colunas como estão
                   names_to = c(".value", "which_moment"), # transformar colunas separadas dos momentos
                   names_pattern = "(mom)(\\d+)") %>%      # em duas colunas contendo o valor e qual a ordem
      mutate(which_moment = as.numeric(which_moment),      # transformar ordem em numérico
             log_mom = log(mom),                           # calcular logaritmo do momento
             log_d = log(d))                               # calcular logaritmo da duração [h]
    
    # Separar estação por momentos e ajustar o modelo linear individualmente
    ls.mom <- split(x = df.mom, f = df.mom$which_moment)
    regression <- lapply(X = ls.mom, FUN = function(gg){
      
      # Expoente de escala
      model <- lm(log_mom ~ log_d, data = gg) # ajustar modelo
      coef <- coef(model)[["log_d"]]          # extrair 'slope'
      scale <- -coef/unique(gg$which_moment)  # expoente de escala
      
      if(is.null(disagg.seq)){
        
        res <- scale
        
      } else{
        
        # Coeficientes de desagregação
        d <- disagg.seq
        D <- max(gg$d) # duração base p/ desagregação
        disagg.coef <- (d/D)^(1 - scale)
        res <- c(scale, disagg.coef)
        names(res) <- c("scale", paste0("coef.d", d))
        
      }
      
      return(res)
      
    }) # fim 'regression'
    
    # Extrair vetor com coeficientes calculados
    scales <- unlist(regression)
    
  } # fim 'fun.lm.boot'
  
  # APLICAR FUNÇÃO --------------------------------------------------------
  
  # Rodar bootstrap
  fun.run.boot <- function(gauge, regional = NULL){
    
    d.groups <- factor(gauge$d)                     # bootstrap estratificado
    boot.obj <- boot::boot(data = gauge,            # série de 1 estação
                           statistic = fun.lm.boot, # função a ser calculada
                           R = R.boot,              # nro. replicas
                           strata = d.groups)       # reamostrar por durações
    
    # Calcular alguma estatística de qualidade do bootstrap
    if(is.null(disagg.seq)){
      
      boot.ci <- broom::tidy(boot.obj, conf.int = TRUE,
                             conf.method = ci.boot.method) %>% # calcular IC95 'bca'
        mutate(d = deparse(which.duration),                # criar coluna intervalo duração
               n_year = n_distinct(gauge$year),            # criar coluna nro. anos
               which_moment = factor(row_number()),        # criar coluna momentos
               gauge_code = ifelse(is.null(regional),
                                   first(gauge$gauge_code),
                                   regional)) %>%   # criar coluna nome estação
        rename(ci_lower = conf.low, ci_upper =  conf.high, std_error = std.error) %>% 
        mutate(which_variable = "scale", target_disagg_d = NA_real_) %>% 
        select(gauge_code, n_year, d, which_moment, which_variable, target_disagg_d,
               statistic, bias, std_error, ci_lower, ci_upper)
      
    } else{
      
      boot.ci <- broom::tidy(boot.obj,
                             conf.int = TRUE,
                             conf.method = ci.boot.method) %>% 
        mutate(d = paste(which.duration, collapse = "-"),
               n_year = n_distinct(gauge$year),
               gauge_code = ifelse(is.null(regional),
                                   first(gauge$gauge_code),
                                   regional)) %>%   # criar coluna nome estação
        rename(ci_lower = conf.low, ci_upper =  conf.high, std_error = std.error) %>% 
        mutate(which_moment = factor(sub("\\..*", "", term)),           # extrair antes primeiro ponto
               which_variable = sub("^[^.]+\\.([^.]+).*", "\\1", term), # extrair scale ou coef
               target_disagg_d = if_else(which_variable == "scale", # somente se for coef (if_else avalia somente TRUE)
                                        NA_real_,             # extrair tudo depois do 'd'
                                        as.numeric(sub(".*d", "", term)))) %>% # NA se não for 'coef'
        select(gauge_code, n_year, d, which_moment, which_variable, target_disagg_d,
               statistic, bias, std_error, ci_lower, ci_upper)
      
    }
    
    return(boot.ci)
    
  } # fim 'fun.run.boot'
  
  # Separar por estação e calcular intervalos de confiança bootstrap
  ls.data <- split(x = data, f = data$gauge_code)
  ls.boot <- pbapply::pblapply(X = ls.data, FUN = fun.run.boot)
  df.boot <- bind_rows(ls.boot)
  
  # Calcular modelo regional
  message("Simulando bootstrap regional: ", n_distinct(data$gauge_code), " estações")
  df.boot.regional <- fun.run.boot(gauge = data, regional = "Regional")
  
  # Transformar em tbl_df
  df.scale.boot <- bind_rows(df.boot, df.boot.regional)
  
  return(df.scale.boot)
  
}