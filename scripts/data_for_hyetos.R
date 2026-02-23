# Este script está destinado a seleção de estações diárias e subdiárias, próximas
# (preferencialmente instaladas no mesmo local) e com séries longas para análise
# com o gerador estocástico de chuvascd.

# Queremos analisar se o gerador mantém a estrutura temporal das séries.

rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

if(!require(pacman)){
  message("Instalando gerenciador de pacotes 'pacman'...")
  install.packages("pacman")
}

# Pacotes
pacman::p_load(pacman, tidyverse, sf, readxl, arrow, geobr, patchwork, tidyterra, ggnewscale, rgdal)


# FUNÇÕES -----------------------------------------------------------------

# Funções
source("scripts/funcoes/fun_close_gauges.R") # fun_close_gauges (calcular distância entre estações)
source("scripts/funcoes/fun_filter_set.R")   # fun_filter_set   (preencher séries subdiárias)


# LER DADOS ---------------------------------------------------------------

# Ler dados diários
path.daily.max <- "base/fonte/consolidado/diario/Máximos anuais para IDF.xlsx"
path.daily <- "base/gerados/df_daily_filled.parquet" # dados já com as datas preenchidas
path.subdaily.info <- "base/gerados/df_subdaily_info.parquet"
path.subdaily.data <- "base/gerados/df_subdaily_data.parquet"
path.imax <- "base/gerados/df_imax.parquet"
path.anadem <- "base/gerados/sig/anadem_rs.tif"
df.daily.max <- readxl::read_excel(path.daily.max)
df.daily.data <- arrow::read_parquet(path.daily)
df.subdaily.info <- arrow::read_parquet(path.subdaily.info)
df.subdaily.data <- arrow::read_parquet(path.subdaily.data)
df.imax <- arrow::read_parquet(path.imax)

# Ler raster RS
raster.anadem <- terra::rast(path.anadem)
anadem <- raster::raster(path.anadem)
plot(anadem)
image(anadem)

terra::plot(raster.anadem)

# Reorganizar 'df.daiy'
df.daily.max <-
  df.daily.max %>% 
  group_by(gauge_code) %>%    # agrupar por estações
  reframe(n_years = n(),      # nro. anos da estação
          long = first(long), # longitude
          lat = first(lat))   # latitude

# Reorganizar 'df.subdaily'
na.prct <- 0.2 # 20% falhas anuais
df.subdaily <- 
  df.imax %>% 
  filter(na_prct <= na.prct, d == 1) %>% # filtrar falhas anuais e somente uma duração
  group_by(gauge_code) %>% # agrupar por estações
  reframe(n_years = n())   # nro. anos da estação

df.subdaily <- merge(df.subdaily, df.subdaily.info)                   # juntar c/ informações de 'df.subdaily.info'
df.subdaily <- df.subdaily[c("gauge_code", "n_years", "long", "lat")] # remover colunas

# Municipalidades e Rio Grande do Sul
sf.rs <- geobr::read_state(code_state = "RS")


# OBJETOS 'SF' ------------------------------------------------------------

# Transformar em objeto 'sf'
# Objetos estão originalmente em coordenada geográfica
crs <- 4674 # SIRGAS 2000
sf.daily <- sf::st_as_sf(x = df.daily.max, coords = c("long", "lat"), crs = crs)
sf.subdaily <- sf::st_as_sf(x = df.subdaily, coords = c("long", "lat"), crs = crs)

# Reprojetar em coordenada projetada
crs.proj <- 31982 # SIRGAS 2000 UTM Zona 22S
sf.daily <- sf::st_transform(x = sf.daily, crs = crs.proj)
sf.subdaily <- sf::st_transform(x = sf.subdaily, crs = crs.proj)

# Juntar em uma única tabela
sf.gauges <- bind_rows(list(daily = sf.daily, subdaily = sf.subdaily), .id = "which")[c(2,1,3,4)]


# DISTÂNCIA ENTRE ESTAÇÕES ------------------------------------------------

# Calcular distância entre estações
ls.gauges <- split(sf.gauges, f = sf.gauges$which)
max.distance <- 30 # filtrar somente estações c/ até 30 m de distância
distance <- fun_close_gauges(sf.x = ls.gauges[["subdaily"]],
                             sf.y = ls.gauges[["daily"]],
                             max.distance = max.distance,
                             gauge.name = "gauge_code")

close.gauges.subdaily <- distance$gauge_subdaily
close.gauges.daily <- distance$gauge_daily

sf.close.daily <- sf.gauges %>% 
  filter(gauge_code %in% close.gauges.daily,
         which == "daily")

sf.close.subdaily <- sf.gauges %>% 
  filter(gauge_code %in% close.gauges.subdaily,
         which == "subdaily")

df.close.gauges <-
  merge(distance, sf.close.daily, by.x = "gauge_daily", by.y = "gauge_code") %>% 
  rename("n_years_daily" = "n_years") %>%  
  select(-c("which", "geometry"))

df.close.gauges <-
  merge(df.close.gauges, sf.close.subdaily, by.x = "gauge_subdaily", by.y = "gauge_code") %>% 
  rename("n_years_subdaily" = "n_years") %>% 
  select(-c("which", "geometry"))

# Extrair 2 observações c/ maior nro. de anos
selected.gauges <- unique(rbind(slice_max(df.close.gauges, order_by = n_years_daily, n = 2),     # diárias
                                slice_max(df.close.gauges, order_by = n_years_subdaily, n = 2))) # subdiárias

selected.gauges <- selected.gauges[selected.gauges$gauge_daily != "SBPA",] # estação do inferno

# Extrair as séries das estações subdiárias selecionadas
df.select.subdaily <- df.subdaily.data[df.subdaily.data$gauge_code %in% selected.gauges$gauge_subdaily,]
df.select.subdaily <- fun_filter_set(data = df.select.subdaily, filter = FALSE)

# Extrair as séries das estações diárias selecionadas
df.select.daily <- df.daily.data[df.daily.data$gauge_code %in% selected.gauges$gauge_daily,] # remover colunsa __index__

# Limpar 'Global Environment'
rm(list = setdiff(ls(), c("sf.rs", "sf.gauges", "df.select.subdaily", "df.select.daily", "selected.gauges", "fun_close_gauges", "fun_filter_set", "raster.anadem")))
invisible(gc())


# EXPORTAR ----------------------------------------------------------------

# Exportar dados como 'PARQUET'
# write.table(selected.gauges, row.names = FALSE, file = "base/gerados/data_for_hyetos/relacao.txt")
# arrow::write_parquet(df.select.daily, sink = "base/gerados/data_for_hyetos/df_select_daily.pqt")
# arrow::write_parquet(df.select.subdaily, sink = "base/gerados/data_for_hyetos/df_select_subdaily.pqt")

# VISUALIZAÇÃO DAS ESTAÇÕES -----------------------------------------------

# Limites de plotagem
lim.y <- c(min(c(df.select.daily$rain_mm, df.select.subdaily$rain_mm), na.rm = TRUE),
           max(c(df.select.daily$rain_mm, df.select.subdaily$rain_mm), na.rm = TRUE))

lim.x <- 
  lubridate::as_datetime(
    c(min(c(df.select.daily$rain_mm, df.select.subdaily$datetime), na.rm = TRUE),
      max(c(df.select.daily$rain_mm, df.select.subdaily$datetime), na.rm = TRUE)))

# Plotar estações subdiárias e diárias lado a lado
plot.select.gauges <- ({

  plot.subdaily <- df.select.subdaily %>% 
    mutate(gauge_code = factor(gauge_code, levels = c("A809", "A844", "A801"))) %>% 
    ggplot(aes(datetime, rain_mm)) +
    facet_wrap(~gauge_code, nrow = length(unique(df.select.subdaily$gauge_code))) +
    # geom_vline(data = df.select.subdaily %>% filter(is.na(rain_mm)), aes(xintercept = datetime), color = "yellow", alpha = 0.5, linewidth = 0.2) +
    geom_line(linewidth = 0.3) +
    lims(x = lim.x, y = lim.y) +
    theme_minimal() +
    labs(x = "Years", y = "Depth [mm]", subtitle = "Subdaily") +
    theme(legend.position =  "none",
          strip.placement = "outside",                      # fora dos rotulos do eixo-y
          strip.text = element_text(size = 12),             # 'strip.' altera configurações do facet_wrap
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),); plot.subdaily

  plot.daily <- df.select.daily %>% 
    mutate(gauge_code = factor(gauge_code, levels = c("83927", "83916", "83967"))) %>% 
    ggplot(aes(datetime, rain_mm)) +
    facet_wrap(~gauge_code, nrow = length(unique(df.select.daily$gauge_code))) +
    # geom_vline(data = df.select.daily %>% filter(is.na(rain_mm)), aes(xintercept = datetime), color = "yellow", alpha = 0.5, linewidth = 0.2) +
    geom_line(linewidth = 0.3) +
    lims(x = lim.x, y = lim.y) +
    labs(x = "Years", y = "Depth [mm]", subtitle = "Daily") +
    theme_minimal() +
    theme(legend.position =  "none",
          strip.placement = "outside",                      # fora dos rotulos do eixo-y
          strip.text = element_text(size = 12),             # 'strip.' altera configurações do facet_wrap
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA)); plot.daily

  plot.both <- 
    plot.daily +
    plot.subdaily + 
    plot_layout(guides = "collect", axes = "collect") & 
    theme(strip.text = element_text(size = 8),       # nome das estações
          plot.subtitle = element_text(hjust = 0.5), # centralizar subtítulo
          text = element_text(color = "black", family = "mono", face = "bold",, size = 10)); plot.both

}) # fim 'plot.select.gauges'

plot.select.gauges

ggsave(plot = plot.select.gauges, filename = "figuras/data_for_hyetos/plot_select_gauges.png",
       width = 16, height = 14, units = "cm", dpi = 300)

# Plotar mapa com localização das estações
# basemap <- maptiles::get_tiles(sf.rs, provider = "OpenStreetMap", zoom = 10)
# basemap <- maptiles::get_tiles(sf.rs, provider = "OpenTopoMap", zoom = 8)

fills <- c("Daily" = "red", "Subdaily" = "blue")
colors <- c("Daily" = "red4", "Subdaily" = "blue4")
shapes <- c("Daily" = 21, "Subdaily" = 25)

data.subdaily <- sf.gauges %>% filter(gauge_code %in% selected.gauges$gauge_subdaily)
data.daily <- sf.gauges %>% filter(gauge_code %in% selected.gauges$gauge_daily)


map <- ({
  
  ggplot() +
    # tidyterra::geom_spatraster_rgb(data = basemap, alpha = 0.7) +
    tidyterra::geom_spatraster(data = raster.anadem) +
    scale_fill_gradient(low = "grey10", high = "grey90", na.value = "transparent") +
    ggnewscale::new_scale_fill() +
    geom_sf(data = sf.rs, color = "black", fill = NA, linewidth = 0.5) +
    geom_sf(data = data.subdaily, aes(fill = "Subdaily", color = "Subdaily", shape = "Subdaily"), size = 3.5, alpha = 0.8) + 
    geom_sf(data = data.daily, aes(fill = "Daily", color = "Daily", shape = "Daily"), size = 3, alpha = 0.8) + 
    geom_sf_text(data = data.subdaily, aes(label = gauge_code), vjust = 0.6, hjust = -0.4,
                 family = "mono", fontface = "bold", color = "blue", size = 3) + 
    geom_sf_text(data = data.daily, aes(label = gauge_code), vjust = -0.7, hjust = -0.4,
                 family = "mono", fontface = "bold", color = "red", size = 3) + 
    scale_fill_manual(values = fills) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(fill = "", shape = "", color = "", x = "", y = "") +
    coord_sf(expand = FALSE) + # preenche espaço entre a borda do plot e o basemap
    theme_minimal() +
    theme(legend.position = c(0.89, 0.09),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_blank(),
          text = element_text(color = "black", family = "mono", face = "bold",, size = 10))
  
}); map

# ggsave(plot = map, filename = "figuras/data_for_hyetos/plot_map_close_gauges.png",
#        width = 16, height = 16, units = "cm", dpi = 300)
