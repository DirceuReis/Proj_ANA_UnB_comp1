
pacman::p_load(pacman, tidyverse, boot, tidymodels)

# Dados subdiários
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")
df.imax <- df.imax[df.imax$gauge_code != "SBPA",] # estação dando problema

# Argumentos
na.accept <- 0.2            # percentual de falhas (0.2 -> 20%)
min.years <- 8              # mínimo de anos na série
which.moment <- 1:3         # qual momento usar de base p/ gerar figuras individuais
which.duration <- c(1/6, 1)  # qual subconjunto de durações
min.duration <- 3           # quantas durações uma estação deve ter no mínimo
R.boot <- 1e4               # quantas amostras bootstrap


data <- df.imax %>%
  filter(na_prct <= na.accept,                # filtrar anos por porcentagem de falhas
         between(d, which.duration[1], which.duration[2])) %>% # filtrar durações dentro do intervalo
  group_by(gauge_code) %>%                    # agrupar por estação
  filter(n_distinct(d) >= min.duration) %>%   # so estações com pelo menos 'min.duration' durações
  ungroup() %>%                               # desagrupar   
  group_by(gauge_code) %>%                    # agrupar por estação
  filter(n_distinct(year) >= min.years) %>%   # so estações com pelo menos 'min.year' anos
  ungroup()


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
  
data.boot <- data %>%
  select(gauge_code, d, imax, year) %>% 
  # filter(gauge_code == gauges[1]) %>%
  group_by(gauge_code) %>% 
  nest(data = c(d, year, imax)) %>%                                   # criar coluna com as séries
  mutate(boot = purrr::map(data, ~boot::boot(data = .x, statistic = fun.lm.boot, R = R.boot, strata = .x$d) %>% # calcular bootstrap
                                 broom::tidy(conf.int = TRUE) %>%     # formatar como tibble
                             mutate(which_moment = row_number()) %>%  # criar coluna indicando momentos
                             rename(scale = statistic))) %>%
  unnest(boot) %>%                                                    # tirar do nest
  rename_with(~gsub(".", "_", .x, fixed = TRUE)) %>%                  # renomear
  rename(ci_lower = conf_low, ci_upper = conf_high) %>% 
  select(-data)
