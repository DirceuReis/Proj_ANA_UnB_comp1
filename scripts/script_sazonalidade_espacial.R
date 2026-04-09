rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, ggh4x, ggimage, sf, geobr, ggspatial, ggrepel)

# FUNÇÕES -----------------------------------------------------------------

# Transformar 'date' em que ocorreram os máximos em radianos
date.to.angle <- function(date){
  
  year.start <- lubridate::floor_date(date, "year")
  year.end <- lubridate::ceiling_date(date, "year")
  year.length <- as.numeric(difftime(year.end, year.start, units = "sec"))
  time.passed <- as.numeric(difftime(date, year.start, units = "sec"))
  
  date.angle <- 2*pi*time.passed/year.length
  return(date.angle)
  
}

# Função p/ gráfico rosa dos ventos
plot.imax.windrose <- function(data,
                               filter.d = NULL,
                               filter.gauge = NULL,
                               seasons,
                               legend.position = "none"){
  
  # Filtros
  if(!is.null(filter.d)) data <- data[data$d %in% filter.d,]
  if(!is.null(filter.gauge)) data <- data[data$gauge_code %in% filter.gauge,]
  
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
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = season), seasons, alpha = 1) +
    geom_segment(aes(x = date_rad, y = 0, xend = date_rad, yend = imax), data, color = "black", linewidth = 0.2) +
    scale_x_continuous(breaks = start.month.rad, labels = month(1:12, label = FALSE), limits = c(0, 2*pi)) +
    scale_fill_manual(values = c("Inverno" = "darkgoldenrod1", "Verão" = "steelblue1")) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0, add = 0)) +
    labs(x = "Mês de ocorrência", y = "Intensidades máximas anuais", fill = "") +
    theme_minimal() +
    theme(panel.border = element_blank(),
          plot.background = element_blank(),
          legend.position = legend.position,
          aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(family = "serif", color = "black", size = 10))
  
}

# Função p/ selecionar subconjunto de estações que melhor cubra o RS
source("scripts/funcoes/scatter2.R")

# PLOT REGIONAL ROSA DOS VENTOS -------------------------------------------

# Analisar a "direção" das rosa dos ventos plotando no local de cada estação
# e avaliando a variação espacial

# Filtros
na.accept <- 0.2
min.years <- 8

# Ler dados
df.info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.pqt")
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.pqt") %>% 
  mutate(lat = df.info$lat[match(gauge_code, df.info$gauge_code)],       # adicionar latitude
         long = df.info$long[match(gauge_code, df.info$gauge_code)]) %>% # adicionar longitude
  filter(na_prct <= na.accept) %>%  
  group_by(gauge_code) %>% 
  mutate(n_years = n_distinct(year)) %>% 
  filter(n_years >= min.years) %>%
  ungroup()

# Rio Grande do Sul
sf.rs <- geobr::read_state(code_state = "RS")

# Atribuir informação geográfica
sf.data <- sf::st_as_sf(df.info, coords = c("long", "lat"))
st_crs(sf.data) <- st_crs(sf.rs) # sirgas 2000

# Regiões hidrogeomorfológicas Rio Grande do Sul
sf.geomorf <- sf::st_read("base/fit_stan/geomorfologia.gpkg")
sf.geomorf <- st_transform(sf.geomorf, crs = st_crs(sf.rs)) # sirgas 2000

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

# ALTERAR -----------------------------------------------------------------

# Duração de interesse [h]
durations <- c(1, 6, 12, 24, 72)

for(d in durations){
  
  message("Gerando mapa para duração de ", d, " h...")
  
  ls.imax <- split(df.imax, df.imax$gauge_code)
  plots <- pbapply::pblapply(X = ls.imax, FUN = function(gauge){
    
    plot.imax.windrose(data = gauge, filter.d = d, seasons = seasons, legend.position = "none") +
      theme(axis.title = element_blank(),
            panel.border = element_rect(color = NA))
    
  }); names(ls.imax) <- plots
  
  df.info$plot <- plots[match(df.info$gauge_code, names(plots))]
  
  # Mapa base 
  basemap <- ggplot() +
    # geom_sf(data = sf.rs, fill = "grey90", color = "black") +
    geom_sf(data = sf.geomorf, fill = "grey90", color = "black") +
    coord_sf()
  
  # Selecionar estações que melhor cobrem o estado
  scatter.gauges <- scatter(pol = sf.rs,
                            pts = sf.data[sf.data$gauge_code %in% gauges,],
                            cell.size = c(0.8, 0.8))
  
  # Preparar dados p/ plotagem
  data.plot <- df.info %>% 
    filter(gauge_code %in% scatter.gauges) %>% 
    mutate(n_years = df.imax$n_years[match(gauge_code, df.imax$gauge_code)]) %>% 
    select(gauge_code, lat, long, n_years, plot)
  
  # Usando o pacote 'ggimage'
  basemap +
    ggimage::geom_subview(data = data.plot, aes(x = long, y = lat, subview = plot), width = 0.7, height = 0.7) +
    ggrepel::geom_text_repel(data = data.plot, aes(x = long, y = lat, label = n_years), family = "Consolas",
                             color = "black", size = 2, vjust = -3, segment.color = NA) +
    ggrepel::geom_text_repel(data = data.plot, aes(x = long, y = lat, label = gauge_code), family = "Consolas",
                             color = "black", size = 2, vjust = 3.5, segment.color = NA) +
    labs(x = "", y = "", subtitle = paste("Duração:", d ,"h |", n_distinct(scatter.gauges), "estações")) +
    ggspatial::annotation_scale(location = "br", pad_y = unit(4, "mm"), height = unit(2, "mm"),
                                bar_cols = c("black", "white"),
                                text_family = "Times New Roman") +
    ggspatial::annotation_north_arrow(location = "tr",pad_x = unit(2, "mm"), pad_y = unit(2, "mm"),
                                      which_north = "true",style = ggspatial::north_arrow_nautical(
                                        fill = c("black", "white"),
                                        line_col = "black",
                                        text_family = "Times New Roman")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          text = element_text(family = "Times New Roman", size = 8))
  
  path <- paste0("figuras/sazonalidade/map_polar_", d, "h.png")
  message("Salvando figura em: ", path)
  
  ggsave(filename = path, height = 16, width = 16, units = "cm", dpi = 300)
  
}
