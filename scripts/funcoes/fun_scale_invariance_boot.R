rm(list = ls()); invisible(gc())
set.seed(123)

# Dados subdiários
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")
df.imax <- df.imax[df.imax$gauge_code != "SBPA",] # estação dando problema

# Argumentos
na.accept <- 0.2            # percentual de falhas (0.2 -> 20%)
min.years <- 8              # mínimo de anos na série
which.moment <- 1:3         # qual momento usar de base p/ gerar figuras individuais
which.duration <- c(1, 24)  # qual subconjunto de durações
min.duration <- 3           # quantas durações uma estação deve ter no mínimo

# Versão simplificada da função fun_scale_invariance, mas que calcula os intervalos
# de confiança do expoente de escala 'H' através do bootstrap

fun_scale_invariance_boot <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                      na.accept,                 # limite de falhas por ano [0,1]
                                      min.years,                 # mínimo de anos em uma série
                                      which.moment = 1:3,        # quais momentos calcular
                                      which.duration = c(1, 24), # intervalo de durações
                                      min.duration = 3,          # mínimo de durações p/ cálculo
                                      R.boot = 100000            # réplicas bootstrap
){
  
  # Pacotes
  if(!require(pacman)){
    install.packages("pacman")
    message("Instalando gerenciador de pacotes 'pacman'...")
  }
  pacman::p_load(pacman, dplyr, tidyr, boot, pbapply)
  
  # MOMENTOS --------------------------------------------------------------
  
  # Selecionar estações adequadas
  data <- df.imax %>%
    filter(na_prct <= na.accept,                # filtrar anos por porcentagem de falhas
           between(d, which.duration[1], which.duration[2])) %>% # filtrar durações dentro do intervalo
    group_by(gauge_code) %>%                    # agrupar por estação
    filter(n_distinct(d) >= min.duration) %>%   # so estações com pelo menos 'min.duration' durações
    ungroup() %>%                               # desagrupar   
    group_by(gauge_code) %>%                    # agrupar por estação
    filter(n_distinct(year) >= min.years) %>%   # so estações com pelo menos 'min.year' anos
    ungroup()                                   # desagrupar
  
  # REGRESSÃO -------------------------------------------------------------
  
  gauges <- unique(data$gauge_code)
  data <- data[data$gauge_code == gauges[1],]
  years <- unique(data$year)
  
  sample.years <- sample(years, replace = TRUE)
  idx <- match(sample.years, years)
  resample <- years[idx]
  n.durations <- length(unique(data$d))
  
  for(i in 2:n.durations){
    idx.data <- c(idx, idx*n.durations[i])
  }
  
  data.resample <- data[,]
  
  # Função de regressão c/ bootstrap
  boot.coef <- function(data, idx){
    
    resample <- data[idx,]                        # amostra bootstrap
    
    # Calcular 'which.moments' momentos centrais locais (p/ cada estação)
    df.mom.local <- data %>% 
      group_by(gauge_code, d) %>% 
      reframe(n_years = n(),
              purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE))))
    
    df.mom.regional <- data %>% 
      group_by(d) %>% 
      reframe(n_years = n(),
              purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE))))
    
    df.mom.regional$gauge_code <- "Regional" # adicionar coluna
    df.mom <- 
      bind_rows(df.mom.local, df.mom.regional) %>%         # juntar
      pivot_longer(cols = -c(gauge_code, d, n_years),      # manter colunas como estão
                   names_to = c(".value", "which_moment"), # transformar colunas separadas dos momentos
                   names_pattern = "(mom)(\\d+)") %>%      # em duas colunas contendo o valor e qual a ordem
      mutate(which_moment = as.numeric(which_moment),      # transformar ordem em numérico
             log_mom = log(mom),                           # calcular logaritmo do momento
             log_d = log(d))                               # calcular logaritmo da duração [h]
    
    # Modelo de regressão linear
    model <- lm(log_mom ~ log_d, data = df.mom) # regressão linear
    scale <- -coef(model)[["log_d"]]/which.moment # extrair coeficiente de escala
    
  }
  
}

na.accept <- 0.2            # percentual de falhas (0.2 -> 20%)
min.years <- 8              # mínimo de anos na série
which.duration <- c(1, 24)  # intervalo de durações
which.moment <- 1:3         # qual momento usar de base p/ gerar figuras individuais
min.duration <- 3
R.boot <- 10000

model <- lm(log_mom ~ log_d, data = gg)

install.packages("car")
library(car)
influencePlot(model)

hist(replicates, breaks = 30)
abline(v = boot.scale$t0, col = "red")
plot(boot.scale)
