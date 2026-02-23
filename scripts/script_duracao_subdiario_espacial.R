# Este script contém somente a análise dos resultados 'horários' da análise de
# invariância de escala -> 'scale.invariance$hourly' gerado em
# 'script_duracao_subdiarios.R'

# PACOTES -----------------------------------------------------------------

rm(list = ls()); invisible(gc())
pacman::p_load(pacman, arrow, tidyverse, sf, ggplot2, geobr, OpenStreetMap, tidyterra)

font.family <- "serif"
font.size <- 12

# LER DADOS ---------------------------------------------------------------

# Dados exportados de 'df.regression'
df_scale_coefficient <- arrow::read_parquet(file = "base/gerados/df_scale_coefficient.parquet")

# Informações sobre estações subdiárias
df_subdaily_info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.parquet")


# ANÁLISE ESPACIAL ESTAÇÕES -----------------------------------------------

# Extrair coeficiente regional
scale.regional <- df_scale_coefficient %>% 
  filter(which_moment == 1,
         gauge_code == "Regional") %>% 
  summarise(scale = first(scale)) %>% 
  pull(scale)

# Calcular diferença regional entre valores do coeficiente de escala
df_scale_coefficient <- df_scale_coefficient %>% 
  filter(which_moment == 1,             # filtrar somente 1º momento
         gauge_code != "Regional") %>%  # remover valor regional
  mutate(relative_diff = (scale - scale.regional)/scale.regional*100) # calcular diferença relativa entre valores de 'H'

# Juntar coordenadas
df_scale_coefficient$long <- df_subdaily_info$long[match(df_scale_coefficient$gauge_code, df_subdaily_info$gauge_code)]
df_scale_coefficient$lat <- df_subdaily_info$lat[match(df_scale_coefficient$gauge_code, df_subdaily_info$gauge_code)]


# HISTOGRAMA DE SCALE -----------------------------------------------------

# Classes de 'relative_diff' conforme os quantis empíricos
relative.diff.quantiles <- round(quantile(df_scale_coefficient$relative_diff)[c(2,3,4)], 2) # só Q1, Q2 e Q3
names(relative.diff.quantiles) <- c("Q1", "Q2", "Q3")

# Estabelecer um limite p/ seleção de estações
limit <- c(-10, 10) #%

plot.scale.hist <- ({
  
  ggplot(df_scale_coefficient, aes(x = relative_diff)) +
    # fact_wrap(~which_moment, ncol = 2) +
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 20, fill = "yellow2", color = "darkorange3") +
    geom_vline(xintercept = relative.diff.quantiles,
               color = "red", linetype = "dashed", linewidth = 0.4) +
    annotate(geom = "text", label = names(relative.diff.quantiles), x = relative.diff.quantiles, y = Inf,
             vjust = 1.5, hjust = 1.2, family = font.family, size = 4) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 0.11), breaks = seq(0, 0.1, 0.02)) +
    labs(x = "(H<sub>i</sub> - H<sub>r</sub>)/H<sub>r</sub> × 100", y = "") +
    theme_minimal() +
    theme(legend.position =  "none",
          axis.title.x = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
          strip.placement = "none",                         # remover título do 'facet' 
          strip.text = element_blank(),                     # 'strip.' altera configurações do facet_wrap
          plot.background = element_rect(color = "white"),
          aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = font.family, color = "black", size = font.size))
  
}); plot.scale.hist


# OBJETOS 'SF' ------------------------------------------------------------

# Transformar em objeto 'sf'
crs <- 4674 # SIRGAS 2000
sf.scale <- sf::st_as_sf(x = df_scale_coefficient, coords = c("long", "lat"), crs = crs)

# Rio Grande do Sul e municipalidades
sf.uf <- geobr::read_state()
sf.rs <- sf.uf[sf.uf$abbrev_state == "RS",]
sf.mun <- geobr::read_municipality(code_muni = "RS")

# Mapa
basemap <- maptiles::get_tiles(sf.rs, provider = "OpenStreetMap", zoom = 10)


# MAPA --------------------------------------------------------------------

# Durações usuadas
duration.interval <- as.numeric(strsplit(gsub("[c()]", "", df_scale_coefficient$d[1]), ", ")[[1]])

# Plotar variação do 
map <- ({
  
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = basemap, alpha = 0.8) +
    geom_sf(data = sf.rs, color = "black", fill = NA) +
    geom_sf(data = sf.scale, aes(color = relative_diff), alpha = 0.9) +
    scale_color_stepsn(n.breaks = 5, colors = c("#EB3200", "#EB9383", "#E1EBE9", "#519DF5", "#0009F0"),
                       labels = relative.diff.quantiles) +
    labs(x = "Latitude (°S)", y = "Longitude (°W)", color = "DR [%]",
         caption = paste0("Durações: ", paste0(duration.interval, collapse = "-"), "h\nFonte de dados: subdiários")) +
    # ggtitle("Diferença relativa entre\ncoeficientes de escala") +
    theme_minimal() +
    theme(legend.position = c(.14, .221),
          legend.background = element_rect(fill = "white", linewidth = 0.15),
          text = element_text(family = "serif"),
          panel.grid = element_blank())
  
}); map

ggsave(filename = "figuras/scale_invariance/hourly_duration/map_scale_coef_relative_difference.png", plot = map,
       width = 16, height = 14, units = "cm", dpi = 200)
