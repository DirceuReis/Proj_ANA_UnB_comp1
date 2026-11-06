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

pacman::p_load(arrow, pacman, tidyverse, tidyterra, sf, geobr, ggnewscale, extrafont, fastDummies)
if(!require(rnaturalearthhires)) remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
suppressMessages(extrafont::loadfonts())


# LER DADOS ---------------------------------------------------------------

# Configurações p/ lidar c/ rasters usando pacote 'terra'
temp.dir <- tempfile() # paste temporária
dir.create(temp.dir)   # criar pasta temporária
terra::terraOptions(tempdir = temp.dir, memfrac = 0.7) # alocar memória

# Modelo Digital de Elevação (ANADEM) do Rio Grande do Sul
path.dem <- "base/fonte/anadem_rs/"                         # caminho dos arquivos (c/ buffer)
dem.files <- list.files(path.dem, pattern = "\\.tif$", full.names = TRUE) # listar arquivos
dem.virtual <- terra::vrt(x = dem.files,                    # cria uma camada virtual
                          filename = "dem_rs_virtual.vrt",
                          overwrite = TRUE)

# Camada espacial Rio Grande do Sul
# sf.rs <- geobr::read_state(code_state = "RS")
sf.rs <- sf::st_read("base/sig/sf_rs.gpkg")

# Estações
sf.gauges <- sf::st_read(dsn = "base/sig/sf_daily.gpkg")
sf.gauges <- sf.gauges %>% 
  slice_tail(n = 1, by = gauge_code) %>% 
  mutate(gauge_code = sub("^0", "", gauge_code)) %>% 
  filter(!gauge_code %in% c("SBPA", "SBUG", "SBSM")) %>% 
  arrange(gauge_code)

# Catchments
sf.hidro <- sf::st_read("base/sig/sf_hidro_rs.gpkg")

# Séries
remove.gauges <- c("SBUG", "SBPA", "SBSM")
df.pmax <- arrow::read_parquet("base/gerados/df_daily_pmax.pqt") %>%
  filter(!gauge_code %in% remove.gauges) %>% 
  arrange(gauge_code)
df.data <- arrow::read_parquet("base/gerados/df_daily_data.pqt") %>%
  filter(!gauge_code %in% remove.gauges) %>% 
  arrange(gauge_code, datetime)


# VISUALIZAR ESTAÇÕES E DEM -----------------------------------------------

font <- "Consolas"
fontsize <- 8
textsize <- (8 - 2)/.pt

# Recortar DEM com sf do Rio Grande do Sul
dem.clip <- terra::mask(x = dem.virtual, mask = sf.rs) # manter somente RS (remover lagoa dos patos)

# Visualizar
ggplot() +
  geom_spatraster(data = dem.clip, aes(fill = after_stat(value)), alpha = 0.8) + # ou fill = names(dem.clip)
  scale_fill_gradientn(name = "Elevation [m]", colors = hcl.colors(100, "viridis"), na.value = NA,
                       guide = guide_colorbar(display = "rectangles", position = "right", barheight = unit(5, "cm"))) +
  geom_sf(data = sf.rs, color = "black", linewidth = 0.6, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = sf.gauges, color = "red2") +
  # coord_sf(crs = 4674, expand = TRUE) +
  theme_minimal() +
  theme(text = element_text(family = "Consolas", size = 12),
        panel.grid = element_blank())

# ggsave(filename = "figuras/cotas.png", width = 16, height = 14, units = "cm", dpi = 300)


# COVARIÁVEIS -------------------------------------------------------------


## DUMMY REGIONS ----------------------------------------------------------

# When examining the regional mean model, we noticed that residuals could follow
# a delimitation loosely based on regional macro-catchments.
# Rio Grange do Sul has three main catchments: Rio Uruguai, Rio Guaíba and Coastal
# The last were combined into a Coastal-Guaíba catchment
# Gauges in Santa Catarina are inside Rio Uruguai's catchment
regions <- sf.hidro %>% 
  mutate(region = ifelse(macro_name == "Uruguai", macro_name, "Guaíba-Coastal"))

dummy.roi <- sf.gauges %>%
  st_join(regions[,c("region")]) %>% 
  mutate(region = ifelse(is.na(region), "Uruguai", region),
         z1 = ifelse(region == "Uruguai", 1, 0),
         # z2 = ifelse(region != "Uruguai", 1, 0) # remove z2 if working with two regions only
         ) %>% 
  as.data.frame()

# old.cov <- arrow::read_parquet("base/gerados/covariates/v20260512_sem_cov_hidrografia_rs/df_covariates.pqt")
#
# old.cov <- old.cov %>%
#   mutate(gauge_code = sub("^0", "", gauge_code)) %>%
#   left_join(dummy.roi[,c("gauge_code", "z1")], by = "gauge_code")
#   # left_join(dummy.roi[,c("gauge_code", "z1", "z2")], by = "gauge_code") # remove z2 if working with two regions only
# 
# arrow::write_parquet(old.cov, "base/gerados/covariates/df_covariates.pqt")


## ASPECT -----------------------------------------------------------------

# Calcular slope com pacote 'terra'
# Aspect is factored in the wrong way, it cannot be a cyclical variable
dem.aspect <- terra::terrain(x = dem.virtual, v = "aspect", unit = "degrees", neighbors = 8)

# Extrair valores para estações
aspect <- terra::extract(x = dem.aspect,             # fonte do dado a ser extraído
                         y = terra::vect(sf.gauges)) # localização de onde extrair

sf.gauges$aspect <- aspect$aspect                    # atribuir aspect (algumas estações estão fora)

# Visualizar 'aspect'
# ggplot() +
#   geom_spatraster(data = dem.aspect, aes(fill = after_stat(value)), alpha = 0.8) + # ou fill = names(dem.clip)
#   scale_fill_gradientn(name = "Aspect [°]", colors = hcl.colors(100, "viiridis"), na.value = NA,
#                        guide = guide_colorbar(display = "rectangles", position = "right", barheight = unit(5, "cm"))) +
#   geom_sf(data = sf.rs, color = "black", linewidth = 0.5, fill = NA) +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         legend.text = element_text(family = "serif", size = 12),
#         text = element_text(family = "serif", size = 12),
#         panel.grid = element_blank())




## SLOPE ------------------------------------------------------------------

# Not a covariate anymore

# Calcular aspect com pacote 'terra'
# dem.slope <- terra::terrain(x = dem.virtual, v = "slope", unit = "degree", neighbors = 8)
# 
# # Extrair valores para estações
# slope <- terra::extract(x = dem.slope,              # fonte do dado a ser extraído
#                         y = terra::vect(sf.gauges)) # localização de onde extrair
# 
# # Removendo 'slope' já que há uma limitação quanto à declividade e o local
# # de instalação de estações pluviométricas, fazendo com que grande parte
# # das estações estejam entre 0-5° de declividade
# # sf.gauges$slope <- slope$slope                      # atribuir slope (algumas estações estão fora)
# 
# # Visualizar 'aspect'
# ggplot() +
#   geom_spatraster(data = dem.slope, aes(fill = after_stat(value)), alpha = 0.8) + # ou fill = names(dem.clip)
#   scale_fill_gradientn(name = "Slope [°]", colors = hcl.colors(100, "viridis"), na.value = NA,
#                        guide = guide_colorbar(display = "rectangles", position = "right", barheight = unit(5, "cm"))) +
#   geom_sf(data = sf.rs, color = "black", linewidth = 0.5, fill = NA) +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         legend.text = element_text(family = "serif", size = 12),
#         text = element_text(family = "serif", size = 12),
#         panel.grid = element_blank())


# MEAN ANNUAL RAINFALL ----------------------------------------------------

# Mean annual rainfall of 1 and 3 days
ds <- c(1, 3, 5, 10)
which.mon <- 1:12

source("scripts/funcoes/fun_pmax_agg_wateryear.R")
pmax <- pmax.wateryear(data = split(df.data, f = df.data$gauge_code),
                       durations = ds[1:2],
                       which.mon = which.mon,
                       names = c("datetime", "rain_mm"))

pmax.summary <- pmax %>% 
  filter(!gauge_code %in% remove.gauges) %>% 
  group_by(gauge_code, d) %>% 
  reframe(mean_annual_rainfall = mean(pmax, na.rm = TRUE))

# old.cov <- arrow::read_parquet("base/gerados/covariates/20260508_base_antiga_ate_2023/df_covariates.pqt")
# old.mean.pmax <- old.cov %>%
#   select(gauge_code, mean_1day, mean_3day) %>%
#   mutate(gauge_code = sub("^0", "", gauge_code))
# 
# new.mean.pmax <-
#   pmax.summary %>%
#   filter(d %in% c(1, 3), !gauge_code %in% remove.gauges) %>%
#   pivot_wider(id_cols = gauge_code,
#               names_from = d,
#               names_glue = "mean_{d}day",
#               values_from = mean_annual_rainfall)
# 
# old.mean.pmax %>%
#   left_join(new.mean.pmax, by = "gauge_code", suffix = c("_old", "_new")) %>%
#   mutate(reldif_1day = (mean_1day_old - mean_1day_new)/mean_1day_old*100,
#          reldif_3day = (mean_3day_old - mean_3day_new)/mean_3day_old*100) %>%
#   print(n = 182)
 
 
# old.cov <- old.cov %>%
#   mutate(gauge_code = sub("^0", "", gauge_code),
#          mean_1day = pmax.summary$mean_annual_rainfall[pmax.summary$d == 1 & match(gauge_code, pmax.summary$gauge_code)],
#          mean_3day = pmax.summary$mean_annual_rainfall[pmax.summary$d == 3 & match(gauge_code, pmax.summary$gauge_code)])
# 
# arrow::write_parquet(old.cov, "base/gerados/covariates/df_covariates.pqt")

# DISTANCE FROM COAST -----------------------------------------------------

# Extract coastline data
sf.coast <- rnaturalearthhires::coastline10 %>% sf::st_transform(crs = st_crs(sf.rs))
bbox.rs <- st_bbox(sf.rs) %>% st_as_sfc() %>% st_buffer(dist = 0.01)
sf.coast <- st_crop(sf.coast, bbox.rs) %>% st_transform(crs = st_crs(sf.gauges))

# Distance to coast
coast.dist <- st_distance(sf.gauges, sf.coast)
sf.gauges$coast_dist <- apply(coast.dist, 1, min)/1000 # in kilometers
sf.segments <- st_nearest_points(sf.gauges, sf.coast) %>% 
  st_sf(gauge_code = sf.gauges$gauge_code,
        coast_dist_m = sf.gauges$coast_dist,
        geometry = .)

ggplot() +
  geom_sf(data = sf.rs, fill = "grey90", color = "black") +
  geom_sf(data = sf.segments, color = "firebrick", linetype = "dashed", alpha = 0.5) +
  geom_sf(data = sf.coast, linewidth = 0.8, color = "dodgerblue3") +
  geom_sf(data = sf.gauges, color = "grey10") +
  labs(subtitle = "Distance to the coastline") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(family = font, size = 12),
        text = element_text(family = font, size = 12, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black"))

ggsave("figuras/covariates/plot_coast_dist.png", width = 16, height = 15, units = "cm", dpi = 300)


# MIN AND MAX TEMPERATURES ------------------------------------------------

# Read data from Xavier (BR-DWGD)
temperatures <- arrow::read_parquet("base/gerados/df_daily_temp.pqt")

tmean.year <- temperatures %>%
  mutate(year = lubridate::year(date),
         gauge_code = sub("^0", "", gauge_code)) %>% 
  group_by(gauge_code, year) %>% 
  reframe(tmean = mean((tmean)))

tmean.summary <- tmean.year %>% 
  group_by(gauge_code) %>% 
  summarise(mean_tmean = mean(tmean))


# SALVAR ------------------------------------------------------------------

df.covariates <- sf.gauges %>% 
  filter(!gauge_code %in% remove.gauges) %>% 
  mutate(elevation = df.pmax$elevation[match(gauge_code, df.pmax$gauge_code)],
         long = -sf::st_coordinates(.)[,1], # RS doesn't cross the equator or greenwich
         lat = -sf::st_coordinates(.)[,2],  # changing signs so all coordinates stay positive
         mean_1day = pmax.summary$mean_annual_rainfall[pmax.summary$d == 1 & match(gauge_code, pmax.summary$gauge_code)],
         mean_3day = pmax.summary$mean_annual_rainfall[pmax.summary$d == 3 & match(gauge_code, pmax.summary$gauge_code)],
         mean_tmean = tmean.summary$mean_tmean[match(gauge_code, tmean.summary$gauge_code)],
         z1 = dummy.roi$z1[match(gauge_code, dummy.roi$gauge_code)]) %>% 
  as.data.frame() %>% 
  select(-c(geom)) %>% 
  arrange(gauge_code)

arrow::write_parquet(x = df.covariates, sink = "base/gerados/covariates/df_covariates.pqt")


# EXPLORING COVARIATES ----------------------------------------------------

# Sampling covariance matrix --- returns LMOM shape estimates
source("scripts/funcoes/fun_scov_matrix.R")

covariates <- arrow::read_parquet("base/gerados/covariates/df_covariates.pqt")
data <- arrow::read_parquet("base/df_daily_pmax.pqt")
names(data)

scov.data <- build.scovmatrix(df = data, col = c("gauge_code", "peak_24h"), which.corr = 2)
kappas.lmom <- scov.data[["df_summary"]][c("gauge_code", "kappa")] %>% rename("kappa_lmom" = "kappa")
covariates$kappa_lmom <- kappas.lmom$kappa_lmom[match(covariates$gauge_code, kappas.lmom$gauge_code)]

kappas <- covariates$kappa_lmom
covariate.matrix <- covariates %>% 
  select(-c(gauge_code, kappa_lmom)) %>% 
  as.matrix()

covariate.matrix <- apply(covariate.matrix, 2, scale, center = TRUE, scale = FALSE)
covariate.matrix <- apply(covariate.matrix, 2, log)

names <- colnames(covariate.matrix)

# Create plots
purrr::walk(.x = 1:ncol(covariate.matrix), .f = function(p){

  ggplot() +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_point(aes(x = kappas, y = covariate.matrix[,p]), alpha = 0.8, color = "black") +
    labs(x = "κ", y = names[p]) +
    theme_minimal() +
    theme(legend.position = "bottom",
          text = element_text(family = font, size = 12),
          panel.border = element_rect(color = "black"))

  ggsave(paste0("figuras/covariates/kappa_", names[p], "_scale.png"), width = 16, height = 12, units = "cm", dpi = 300)

})

covariates %>%
  pivot_longer(cols = -c(gauge_code, kappa_lmom),
               names_to = "covariate",
               values_to = "value") %>%
  ggplot(aes(x = value)) +
  facet_wrap(~covariate, nrow = 2, scales = "free_x") +
  geom_histogram(color = "black", fill = "grey", linewidth = 0.3, bins = 30) +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(family = font, size = 12, color = "black"),
        panel.border = element_rect(color = "black"))
# 
# ggsave("figuras/covariates/covariate_histogram.png", width = 20, height = 16, units = "cm", dpi = 300)
