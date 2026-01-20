# Este script está destinado a calcular covariáveis (variáveis preditoras) regionais para 
# construir um modelo regional de mínimos quadrados para o parâmetro de forma da GEV

# Os arquivos foram baixados usando o código do GEE localizado em 'scripts/gee/extract_anadem_rs.txt'
# que exporta

# PACOTES -----------------------------------------------------------------

rm(list = ls()); invisible(gc())

if(!require(pacman)){
  message("Instalando gerenciador de pacotes 'pacman'...")
  install.packages("pacman")
}

pacman::p_load(arrow, pacman, tidyverse, tidyterra, sf, geobr, ggnewscale)


# LER DADOS ---------------------------------------------------------------

# Configurações p/ lidar c/ rasters usando pacote 'terra'
temp.dir <- tempfile() # paste temporária
dir.create(temp.dir)   # criar pasta temporária
terra::terraOptions(tempdir = temp.dir, memfrac = 0.7) # alocar memória

# Modelo Digital de Elevação (ANADEM) do Rio Grande do Sul
# path.dem <- "C:/Users/tomas/Downloads/anadem_rs"                          # caminho dos arquivos
path.dem <- "C:/Users/tomas/Downloads/anadem_rs_buffer"                          # caminho dos arquivos (c/ buffer)
dem.files <- list.files(path.dem, pattern = "\\.tif$", full.names = TRUE) # listar arquivos
dem.virtual <- terra::vrt(x = dem.files,                                  # cria uma camada virtual
                          filename = "dem_rs_virtual.vrt",
                          overwrite = TRUE)

# Camada espacial Rio Grande do Sul
sf.rs <- geobr::read_state(code_state = "RS")
# Estações
sf.gauges <- st_read(dsn = "base/gerados/sig/sf_rs_estacoes.gpkg")


# VISUALIZAR ESTAÇÕES E DEM -----------------------------------------------

# Recortar DEM com sf do Rio Grande do Sul
dem.clip <- terra::mask(x = dem.virtual, mask = sf.rs) # manter somente RS (remover lagoa dos patos)

# Gráfico
fills <- c("Pluviômetros" = "red", "Pluviógrafos" = "blue")
colors <- c("Pluviômetros" = "red4", "Pluviógrafos" = "blue4")
shapes <- c("Pluviômetros" = 21, "Pluviógrafos" = 25)

# Visualizar
ggplot() +
  geom_spatraster(data = dem.clip, aes(fill = after_stat(value)), alpha = 0.8) + # ou fill = names(dem.clip)
  scale_fill_gradientn(name = "Cota [m]", colors = hcl.colors(100, "viridis"), na.value = NA,
                       guide = guide_colorbar(display = "rectangles", position = "right", barheight = unit(5, "cm"))) +
  geom_sf(data = sf.rs, color = "black", linewidth = 0.5, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = sf.gauges, aes(fill = res, color = res, shape = res)) +
  # coord_sf(crs = 4674, expand = TRUE) +
  scale_fill_manual(values = fills) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  labs(color = "", fill = "", shape = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(family = "serif", size = 12),
        text = element_text(family = "serif", size = 12),
        panel.grid = element_blank())

# ggsave(filename = "figuras/cotas.png", width = 16, height = 14, units = "cm", dpi = 300)


# COVARIÁVEIS -------------------------------------------------------------

## ASPECT -----------------------------------------------------------------

# Calcular slope com pacote 'terra'
dem.aspect <- terra::terrain(x = dem.virtual, v = "aspect", unit = "degrees", neighbors = 8)

# Extrair valores para estações
aspect <- terra::extract(x = dem.aspect,             # fonte do dado a ser extraído
                         y = terra::vect(sf.gauges)) # localização de onde extrair
sf.gauges$aspect <- aspect$aspect                    # atribuir aspect (algumas estações estão fora)

# Visualizar 'aspect'
ggplot() +
  geom_spatraster(data = dem.aspect, aes(fill = after_stat(value)), alpha = 0.8) + # ou fill = names(dem.clip)
  scale_fill_gradientn(name = "Aspect [°]", colors = hcl.colors(100, "viridis"), na.value = NA,
                       guide = guide_colorbar(display = "rectangles", position = "right", barheight = unit(5, "cm"))) +
  geom_sf(data = sf.rs, color = "black", linewidth = 0.5, fill = NA) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(family = "serif", size = 12),
        text = element_text(family = "serif", size = 12),
        panel.grid = element_blank())


## SLOPE ------------------------------------------------------------------

# Calcular aspect com pacote 'terra'
dem.slope <- terra::terrain(x = dem.virtual, v = "slope", unit = "degree", neighbors = 8)

# Extrair valores para estações
slope <- terra::extract(x = dem.slope,              # fonte do dado a ser extraído
                        y = terra::vect(sf.gauges)) # localização de onde extrair
sf.gauges$slope <- slope$slope                      # atribuir slope (algumas estações estão fora)

# Visualizar 'aspect'
ggplot() +
  geom_spatraster(data = dem.slope, aes(fill = after_stat(value)), alpha = 0.8) + # ou fill = names(dem.clip)
  scale_fill_gradientn(name = "Slope [°]", colors = hcl.colors(100, "viridis"), na.value = NA,
                       guide = guide_colorbar(display = "rectangles", position = "right", barheight = unit(5, "cm"))) +
  geom_sf(data = sf.rs, color = "black", linewidth = 0.5, fill = NA) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(family = "serif", size = 12),
        text = element_text(family = "serif", size = 12),
        panel.grid = element_blank())

# SALVAR ------------------------------------------------------------------

df.covariaveis <- sf.gauges %>% 
  mutate(lat = sf::st_coordinates(.)[,2],
         long = sf::st_coordinates(.)[,1]) %>% 
  as.data.frame() %>% 
  select(-geom)

arrow::write_parquet(x = df.covariaveis, sink = "base/gerados/df_covariaveis.parquet")
