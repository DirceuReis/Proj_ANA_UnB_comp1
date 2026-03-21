# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, ggh4x, ggimage, ggsubplot, sf, geobr)


# FUNÇÕES -----------------------------------------------------------------

# Função p/ gráfico rosa dos ventos
plot.imax.windrose <- function(data,
                               filter.d = NULL,
                               filter.gauge = NULL,
                               seasons,
                               legend.position = "none"){
  
  # Filtros
  if(!is.null(filter.d)) data <- data[data$d %in% filter.d,]
  if(!is.null(filter.d)) data <- data[data$gauge_code %in% filter.gauge,]
  
  # Transformar 'date' em que ocorreram os máximos em radianos
  date.to.angle <- function(date){
    
    year.start <- lubridate::floor_date(date, "year")
    year.end <- lubridate::ceiling_date(date, "year")
    year.length <- as.numeric(difftime(year.end, year.start, units = "sec"))
    time.passed <- as.numeric(difftime(date, year.start, units = "sec"))
    
    date.angle <- 2*pi*time.passed/year.length
    return(date.angle)
    
  }
  
  start.month <- lubridate::yday(as.Date(paste0("1800/", 1:12, "/01")))
  start.month.rad <- 2*pi*start.month/365.25
  
  # Calcular angulos
  data$date_rad <- date.to.angle(data$date)
  
  # Plot
  ggplot() +
    coord_polar() +
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = season), seasons, alpha = 0.4) +
    geom_segment(aes(x = date_rad, y = 0, xend = date_rad, yend = imax), data, color = "black", linewidth = 1) +
    scale_x_continuous(breaks = start.month.rad, labels = month(1:12, label = FALSE), limits = c(0, 2*pi)) +
    scale_fill_manual(values = c("Inverno" = "darkgoldenrod1", "Verão" = "steelblue")) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Mês de ocorrência", y = "Intensidades máximas anuais", fill = "") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", linewidth = 0.3),
          panel.background = element_rect(color = "white"),
          legend.position = legend.position,
          aspect.ratio = 1,
          panel.grid.minor = element_blank(),
          plot.background = element_rect(color = NA),
          # plot.title = element_text(size = 8, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          # axis.text.x = element_text(vjust = -0.9, color = "black", size = 8),
          axis.text.x = element_text(size = 8, margin = margin(b = 1)),
          text = element_text(family = "serif", color = "black", size = 10))
  
}

# PLOT REGIONAL ROSA DOS VENTOS -------------------------------------------

# Analisar a "direção" das rosa dos ventos plotando no local de cada estação
# e avaliando a variação espacial

# Ler dados
df.info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.pqt")
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.pqt") %>% 
  mutate(lat = df.info$lat[match(gauge_code, df.info$gauge_code)],   # adicionar latitude
         long = df.info$long[match(gauge_code, df.info$gauge_code)]) # adicionar longitude

# Atribuir informação geográfica
sf.imax <- sf::st_as_sf(df.imax, coords = c("long", "lat"))

# Rio Grande do Sul
sf.rs <- geobr::read_state(code_state = "RS")

# Definir meses de início dos periodos
umido1 <- as.Date(c("1800-01-01", "1800-04-01"))
seco <- as.Date(c("1800-04-01", "1800-09-01"))
umido2 <- as.Date(c("1800-09-01", "1800-12-31"))

seasons <- rbind(date.to.angle(umido1),
                 date.to.angle(seco),
                 date.to.angle(umido2) + c(0, + 2*pi/365.25)) %>% 
  as.data.frame() %>%
  reframe(start = V1, end = V2) %>% 
  mutate(season = c("Verão", "Inverno", "Verão"))

# Avaliar quais são as estações que apresentam boa cobertura do estado
# Duração: 1 h
ls.imax <- split(df.imax, df.imax$gauge_code)
plots <- pbapply::pblapply(X = ls.imax, FUN = function(gauge){
  
  plot.imax.windrose(data = gauge, filter.d = 1, seasons = seasons, legend.position = "none") +
    theme(axis.title = element_blank(),
          panel.border = element_rect(color = NA))
  
})

basemap <- ggplot() +
  geom_sf(data = sf.rs, fill = "grey90", color = "black") +
  coord_sf() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank())

# Pacote `ggimage()` 
plot.regional <- basemap +
  annotation_custom(grob = ggplotGrob(plots[[1]]),
                    xmin = df.imax$long[1] - 0.3,
                    xmax = df.imax$long[1] + 0.3,
                    ymin = df.imax$lat[1] + -0.3,
                    ymax = df.imax$lat[1] + 0.3); plot.regional
