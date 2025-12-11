
# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, arrow, tidyverse, sf, ggplot2, geobr, OpenStreetMap, tidyterra)


# LER DADOS ---------------------------------------------------------------

# Dados exportados de 'df.regression'
df_scale <- arrow::read_parquet(file = "base/gerados/df_scale_invariance.parquet")

# Informações sobre estações subdiárias
df_subdaily_info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.parquet")


# ANÁLISE ESPACIAL ESTAÇÕES -----------------------------------------------

# Extrair coeficiente regional
scale.regional <- df_scale %>% 
  filter(which_moment == 1,
         gauge_code == "Regional") %>% 
  summarise(scale = mean(scale)) %>% 
  pull(scale)

# Juntar coordenadas
df_scale <- df_scale %>% 
  filter(which_moment == 1,
         gauge_code != "Regional") %>% 
  group_by(gauge_code) %>%
  summarise(n_years = mean(n_years),
            rsquared = mean(rsquared),
            scale_local = mean(scale)) %>% 
  mutate(scale_regional = scale.regional,
         relative_diff = (scale_local - scale_regional)/scale_regional*100)

df_scale$long <- df_subdaily_info$long[match(df_scale$gauge_code, df_subdaily_info$gauge_code)]
df_scale$lat <- df_subdaily_info$lat[match(df_scale$gauge_code, df_subdaily_info$gauge_code)]

# Transformar em objeto 'sf'
crs <- 4674 # SIRGAS 2000
sf.scale <- sf::st_as_sf(x = df_scale, coords = c("long", "lat"), crs = crs)

# Rio Grande do Sul e municipalidades
sf.uf <- geobr::read_state()
sf.rs <- sf.uf[sf.uf$abbrev_state == "RS",]
sf.mun <- geobr::read_municipality(code_muni = "RS")

# Mapa
basemap <- maptiles::get_tiles(sf.rs, provider = "OpenStreetMap", zoom = 10)
map <- ({
  
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = basemap) +
    geom_sf(data = sf.rs, color = "black", fill = NA) +
    geom_sf(data = sf.scale, aes(color = scale), alpha = 0.8) +
    scale_color_binned() +
    labs(x = "Latitude (°S)", y = "Longitude (°W)", title = "Coeficiente de escala: horário", color = "H") +
    theme_minimal() +
    theme(legend.position = "right",
          text = element_text(family = "serif"),
          panel.grid = element_blank())
  
}); map

ggsave(filename = "figuras/scale_invariance/hourly_duration/map_scale_coef.png", plot = map,
       width = 16, height = 14, units = "cm", dpi = 200)
