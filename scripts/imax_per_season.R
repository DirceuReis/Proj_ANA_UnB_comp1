
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow, ggforce)


# FUNÇÕES -----------------------------------------------------------------

source("scripts/funcoes/fun_imax_agg.R")
source("scripts/funcoes/fun_imax_agg_wateryear.R")
source("scripts/funcoes/fun_filter_set.R")
source("scripts/funcoes/fun_group_ts.R")


# EXTRAIR AMS -------------------------------------------------------------

# Rodar a seção abaixosomente uma vez, possivelmente os arquivos já estão salvos
if(!(file.exists("base/gerados/df_imax_wet.pqt") & file.exists("base/gerados/df_imax_dry.pqt"))){
  
  df.data <- arrow::read_parquet("base/gerados/df_subdaily_data.pqt")
  
  df.data <- fun_filter_set(data = df.data,
                            daily = FALSE,
                            col_names = names(df.data),
                            filter = FALSE); invisible(gc())
  
  ls.data <- split(df.data, f = df.data$gauge_code)
  data.by.time <- fun_group_ts(ls.data, ts_name = "time_step")
  
  d.subhour <- c(10, 15, 20, 30, 40, 45, 50)
  d.subdaily <- seq(1, 23)
  d.daily <- 1:10
  durations <- c(d.subhour/60, d.subdaily, d.daily*24)
  time.steps <- names(data.by.time)
  
  ls.seasons <- list(months.wet = c(10, 11, 12, 1, 2, 3),
                     months.dry = c(4, 5, 6, 7, 8, 9))
  
  imax.seasons <- lapply(X = ls.seasons, FUN = function(s){
    
    message("Processando estações com os meses: ", paste(s, collapse = ", "))
    
    imax.ls <- lapply(X = time.steps, FUN = function(ts){
      
      current.ts <- data.by.time[[ts]]
      ds <- as.numeric(ts)/3600
      valid.durations <- durations[(durations/ds) %% 1 < 1e-8]
      
      message(paste0("Processando conjunto de estações com resolução de ", ds*60, " minutos..."))
      
      out <- fun_imax_agg(data = current.ts,
                          durations = valid.durations,
                          which.mon = s,
                          names = c("datetime", "rain_mm")) 
      
    })
    
  })
  
  df.imax.wet <- arrow::write_parquet(dplyr::bind_rows(imax.seasons$months.wet), "base/gerados/df_imax_wet.pqt")
  df.imax.dry <- arrow::write_parquet(dplyr::bind_rows(imax.seasons$months.dry), "base/gerados/df_imax_dry.pqt")
  
}


# AMS POR ESTAÇÃO DO ANO --------------------------------------------------

min.years <- 8
durations <- c(1, 6, 12, 24, 72)
ds <- 1

# Informações
df.info <- arrow::read_parquet("base/gerados/df_subdaily_info.pqt")

# Série de intensidades acumuladas
df.imax.s <- bind_rows(
  arrow::read_parquet("base/gerados/df_imax_wet.pqt") %>% mutate(season = "wet"),
  arrow::read_parquet("base/gerados/df_imax_dry.pqt") %>% mutate(season = "dry")
) %>% 
  mutate(responsible = df.info$responsible[match(gauge_code, df.info$gauge_code)]) %>% 
  filter(d %in% durations) %>%
  group_by(gauge_code, year) %>% 
  filter(na_prct < 1) %>%             # remover so anos c/ 100% de falhas
  ungroup() %>% 
  group_by(gauge_code, d, season) %>% 
  filter(n_distinct(year) >= min.years) %>% # manter estações c/ + 8 anos
  mutate(n_years = n_distinct(year),  # calcular nro. anos em cada estação
         imax_mean = mean(imax),      # calcular intensidades médias
         imax_cum = cumsum(imax)) %>% # calcular intensidades acumuladas
  ungroup()

labels <- df.imax.s %>% 
  filter(d == ds) %>% 
  group_by(gauge_code) %>% 
  summarise(responsible = first(responsible),
            n_years = first(n_years))

# Gerar figuras imax_cum
for(ds in durations){
  
  message("Gerando gráficos para duração de ", ds, " h...")
  message("Usando estações com pelo menos ", min.years, " anos\n")
  
  df.imax.s %>% 
    # filter(gauge_code == "A801") %>% 
    # filter(responsible != "SGB") %>% 
    filter(d == ds) %>% 
    group_by(gauge_code, season, d) %>% 
    mutate(label_x = min(year), label_y = max(imax_cum)) %>% 
    ungroup() %>% 
    ggplot(aes(x = year, y = imax_cum)) +
    ggforce::facet_wrap_paginate(~gauge_code, scales = "free", ncol = 3, nrow = 3) +
    geom_ribbon(aes(ymin = 0, ymax = imax_cum, color = season, fill = season), alpha = 0.3) +
    geom_line(aes(color = season)) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = paste(responsible, "|", n_years, "anos")), hjust = -0.1, vjust = 1.5,
              family = "Consolas", size = 2, color = "black") +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_color_manual(values = c("dry" = "darkgoldenrod1", wet = "cadetblue2"), aesthetics = c("fill", "color")) +
    labs(x = "", y = "Intensidades acumuladas [mm/h]", color = "", fill = "", caption = paste("Duração:", ds, "h")) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black"),
          strip.text = element_text(size = 8, margin = margin(t = 0.2, b = 0.5)),
          text = element_text(color = "black", family = "Times New Roman", size = 10),
          legend.position = "none")
  
  ggsave(filename = paste0("figuras/sazonalidade/plot_imax_cum_por_estacao_", ds, "h.png"),
         width = 410, height = 297, units = "mm", dpi = 300)
  
}

# Gerar figuras densidade
min.years <- 10
p.dim <- c(3,4)
dir.name <- "figuras/sazonalidade/densidades/"

for(ds in durations){
  
  message("\nGerando gráficos de densidade para durações de ", ds, " h")
  message("Usando estações com pelo menos ", min.years, " anos")
  
  # Calcular distribuições densidade
  p.density <- df.imax.s %>% 
    # filter(responsible != "SGB") %>%
    filter(n_years >= min.years, d == ds) %>% 
    ggplot(aes(x = imax)) +
    ggforce::facet_wrap_paginate(~gauge_code, scales = "free", nrow = p.dim[1], ncol = p.dim[2]) +
    geom_density(aes(color = season, fill = season), alpha = 0.3) +
    geom_vline(aes(xintercept = imax_mean, color = season), linetype = "dashed") +
    geom_text(data = labels %>% filter(n_years >= min.years), aes(x = -Inf, y = Inf, label = paste(responsible, "|", n_years, "anos")),
              hjust = -0.1, vjust = 1.6, family = "Consolas", size = 3, color = "black") +
    scale_color_manual(values = c("dry" = "darkgoldenrod1", wet = "cadetblue2"), aesthetics = c("fill", "color")) +
    labs(x = "Intensidade máxima anual [mm/h]", y = "Densidade", color = "", fill = "", caption = paste("Duração:", ds, "h")) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black"),
          strip.text = element_text(size = 8, margin = margin(t = 0.5, b = 1)),
          text = element_text(color = "black", family = "Times New Roman", size = 10),
          legend.position = "none")
  
  p.seq <- 1:n_pages(p.density) # número de páginas necessárias para plotar o gráfico completo
  
  # Dividir o gráfico acima em 'n_pages' e salvar cada um em uma lista
  ls.p.density <- lapply(X = p.seq, FUN = function(group){
    p.density + ggforce::facet_wrap_paginate(~gauge_code, nrow = p.dim[1], ncol = p.dim[2], page = group, scales = "free")
  }) # fim 'ls.p.density'
  
  if(!dir.exists(dir.name)) dir.create(dir.name)
  for(group in seq_along(ls.p.density)){
    
    message(group)
    ggsave(plot = ls.p.density[[group]], filename = paste0(dir.name, "plot_densidades_", ds, "h_", group,".png"),
           width = 297, height = 210, units = "mm", dpi = 300)
    
  }
  
}

# Selecionar estações que os totais acumulados do período seco superam
# os totais acumulados do período úmido
season.change <- df.imax.s %>% 
  group_by(gauge_code, season, d) %>% 
  slice_tail(n = 1) %>% # selecionar última linha (total de imax_cum)
  pivot_wider(id_cols = c(gauge_code, d, year, responsible),
              names_from = season,
              values_from = imax_cum) %>% 
  mutate(change = if_else(dry > wet, -1, 1),
         rel_dif = (wet - dry)/wet*100) %>% 
  ungroup() %>% 
  pivot_longer(cols = -c("gauge_code", "d", "responsible", "year", "change", "rel_dif"),
               names_to = "season",
               values_to = "imax_cum")
  
season.change %>% 
  filter(d == 1) %>% 
  ggplot(aes(gauge_code, ))
  


# FIT GEV -----------------------------------------------------------------

rm(list = setdiff(ls(), "df.imax.s")); invisible(gc())

# Funções
if(!file.exists("scripts/funcoes/fit_gev/gev_aux.R")){
  
  stop("RScript 'gev_aux.R' with auxiliary functions not found in:\n", getwd())
  
} else{
  
  # Ler funções de ajuste, quantis e funções auxiliarees
  invisible(lapply(X = list.files("scripts/funcoes/fit_gev", pattern = ".R", full.names = TRUE), FUN = source))
  
}

# Pacotes
pacman::p_load(pacman, lmom)

min.years <- 10
ds <- 1
labels <- df.imax.s %>% 
  filter(d == ds) %>% 
  group_by(gauge_code) %>% 
  summarise(responsible = first(responsible),
            n_years = first(n_years))

# Separar dados por estação e por duração
ls.imax.s <- lapply(X = split(df.imax.s[df.imax.s$n_years >= min.years,], df.imax.s[df.imax.s$n_years >= min.years, "season"]),
                    FUN = function(s) split(s, s$d))

# Ajustar GEV
ls.fit <- lapply(ls.imax.s, function(ls.s){
  
  message("Fitting AMS from '", ls.s[[1]][["season"]][1], "' season with LMOM.")
  
  ls.duration <- pbapply::pblapply(ls.s, function(data.imax){
      
    ls.df <- split(data.imax, data.imax$gauge_code)
    ls.fitted <- lapply(ls.df, function(df){
      
      imax <- df$imax
      imax.pos <- imax[imax > 0]
      par.lmom.gev <- fun.lmom(xp = imax.pos, dist = "gev")
      
      res <- tibble(gauge_code = df$gauge_code[1],
                    d = df$d[1],
                    season = df$season[1],
                    par_lmom_gev = list(tibble(par.lmom.gev)))
        
      }) # fim 'ls.fitted'
    
    bind_rows(ls.fitted)
    
  }) # fim 'ls.duration'
  
  bind_rows(ls.duration)
    
}) # fim 'ls.fit'

df.fit <- bind_rows(ls.fit)

# Calcular quantis GEV
# Tempos de retorno
tr <- c(1.1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)

# Abordagem c/ purrr
df.q <- lapply(ls.fit, function(df.s){
  
  df.s %>% 
    mutate(q = map(par_lmom_gev, ~ fun.q.gev(p = 1 - 1/tr, param = unlist(.x)))) %>% 
    unnest_longer(q) %>% 
    mutate(tr = rep(tr, times = n()/length(tr)))
  
}) %>% bind_rows()

# Gerar figuras densidade
durations <- c(1, 24, 72)
min.years <- 10
p.dim <- c(3,4)
dir.name <- "figuras/sazonalidade/frequencia/"

# Atualizar ano hidrológico
shift.wateryear <- function(dates, start.month){
  
  dates <- as.Date(dates)
  shift <- dplyr::if_else(lubridate::month(dates) < start.month, 0, 1)
  wateryear <- lubridate::year(dates) + shift
  
  return(wateryear)
  
}

wateryear.breaks <- c(9,2)
df.imax.s <- df.imax.s %>% mutate(wateryear = shift.wateryear(date, wateryear.breaks[1]))

# Razão entre períodos de maior chuva
# Alguns anos hidrológicos têm dois valores de imax
ratios <- df.imax.s %>%
  group_by(gauge_code) %>% 
  mutate(n_years = n_distinct(year)) %>% 
  group_by(gauge_code, d, wateryear, n_years, season) %>% 
  summarise(imax = max(imax, na.rm = TRUE), .groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(gauge_code, d, wateryear, n_years),
              names_from = season,
              values_from = imax) %>% 
  mutate(which_imax = if_else(wet > dry, "wet", "dry")) %>% 
  group_by(gauge_code, d, n_years) %>% 
  summarise(count_wet = sum(which_imax == "wet", na.rm = TRUE),
            count_dry = sum(which_imax == "dry", na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(ratio = (count_wet/n_years)/(count_dry/n_years),
         label = paste0(round(count_wet/n_years*100,1), "%/", round(count_dry/n_years*100,1), "%"))

for(ds in durations){
  
  message("\nGerando gráficos de frequência para durações de ", ds, " h")
  message("Usando estações com pelo menos ", min.years, " anos")
  
  # Calcular distribuições densidade
  p.frequency <- df.q %>% 
    filter(d == ds) %>% 
    ggplot(aes(x = tr, y = q, color = season)) +
    ggforce::facet_wrap_paginate(~gauge_code, scales = "free", nrow = p.dim[1], ncol = p.dim[2]) +
    geom_line() +
    geom_text(data = labels %>% filter(n_years >= min.years), aes(x = 1, y = Inf, label = paste(responsible, "|", n_years, "anos")),
              hjust = 0, vjust = 2, family = "Consolas", size = 3, color = "black") +
    geom_text(data = ratios %>% filter(d == ds), aes(x = 1, y = Inf, label = label),
              hjust = 0, vjust = 4.2, family = "Consolas", size = 3, color = "black") +
    scale_color_manual(values = c("dry" = "darkgoldenrod1", wet = "cadetblue2")) +
    scale_x_log10(breaks = c(1, 10, 20, 30, 50, 100, 200, 500),
                  labels = c("1", "10", "20", "30", "50", "100", "200", "500")) +
    labs(x = "Return period [years]", y = "Estimated intensity [mm/h]", color = "", caption = paste("Duration:", ds, "h | Method: L-moments")) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black"),
          strip.text = element_text(size = 8, margin = margin(t = 0.5, b = 1)),
          text = element_text(color = "black", family = "Times New Roman", size = 10),
          legend.position = "none"); p.frequency
  
  p.seq <- 1:n_pages(p.frequency) # número de páginas necessárias para plotar o gráfico completo
  
  # Dividir o gráfico acima em 'n_pages' e salvar cada um em uma lista
  ls.p.frequency <- lapply(X = p.seq, FUN = function(group){
    p.frequency + ggforce::facet_wrap_paginate(~gauge_code, nrow = p.dim[1], ncol = p.dim[2], page = group, scales = "free")
  }) # fim 'ls.p.density'
  
  if(!dir.exists(dir.name)) dir.create(dir.name)
  for(group in seq_along(ls.p.frequency)){
    
    message(group)
    ggsave(plot = ls.p.frequency[[group]], filename = paste0(dir.name, "plot_frequencias_", ds, "h_", group,".png"),
           width = 297, height = 210, units = "mm", dpi = 300)
    
  }
  
}

# FREQUÊNCIA IMAX ANO CIVIL -----------------------------------------------

rm(list = ls()); invisible(gc())

# Ler dados imax ano civil e remover estações SGB
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")

# 