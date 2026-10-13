
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(geobr)
library(terra)
library(tidyterra)
library(gstat)
library(stars)


# TEMA 'GGPLOT' -----------------------------------------------------------

font <- "sans"
fontsize <- 8
textsize = fontsize/.pt
set_theme(
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = c(0, 0),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key.height = unit(fontsize, "pt"),
      strip.background = element_blank(),
      text = element_text(size = 8),
      panel.spacing.y = unit(0.5, "pt"),
    )
)


# CAMINHOS ----------------------------------------------------------------

root <- "Relatório 2"
path.inv.scale <- file.path(root, "Scripts_Invariancia_Escala", "resultados", "subdiarios")
path.imax <- file.path(root, "Scripts_Maximas", "resultados", "BR", "dataframes")
path.plot <- file.path(root, "Scripts_Invariancia_Escala", "figuras")


# LER DADOS ---------------------------------------------------------------

mom.df <- arrow::read_parquet(file.path(path.inv.scale, "df_momentos_subdiario.pqt"))               # momentos nao centrais p/ intensidades maximas anuais subdiarias
scale.df <- arrow::read_parquet(file.path(path.inv.scale, "df_expoente_escala_subdiario.pqt"))      # expoentes de escala 'H'
imax.preqc.df <- arrow::read_parquet(file.path(path.imax, "imax_sem_iqr_subdiario_Uniplu.parquet")) # intensidades maximas anuaisi
info.df <- readRDS(file.path(path.imax, "metadados_estacoes_subdiario_Uniplu.rds"))                 # informacoes


# PREPARAR DADOS ----------------------------------------------------------

na.accept <- 0.2
min.years <- 8
min.durations <- 3
epsg <- 4674

# Intensidades maximas anuais
imax.df <- imax.preqc.df |> 
  filter(na_prct_yr <= na.accept) |>
  group_by(gauge_code) |> 
  mutate(n_years = n_distinct(wateryear)) |> 
  filter(
    n_distinct(d) >= min.durations,
    n_years >= min.years
  ) |> 
  ungroup() |> 
  rename(imax = imax_mm_h) |> 
  select(gauge_code, wateryear, d, imax, n_years, state) |> 
  left_join(
    info.df |> 
      select(gauge_code, long, lat),
    by = "gauge_code"
  )

gauges <- unique(imax.df$gauge_code)

info.df <- info.df |> filter(gauge_code %in% gauges)

# Momentos estimados p/ cada estacao
mom.plot.df <- mom.df |> 
  rename(n_years = n_year) |> 
  left_join(
    imax.df |> 
      select(-c(wateryear, imax)) |> 
      slice_head(n = 1, by = c(gauge_code, d)),
    by = c("gauge_code", "d")
  )

# Dados geoespaciais
uf.sf <- geobr::read_state(year = 2025)
br.sf <- geobr::read_country(year = 2025)
world.sf <- rnaturalearth::ne_countries(continent = "south america", scale = "large")
subdaily.sf <- info.df |> 
  select(-policy) |> 
  mutate(n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)]) |> 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(epsg))


# VISUALIZACAO ------------------------------------------------------------

## VARIACAO ESPACIAL 'H' --------------------------------------------------

# Funcao p/ gerar quebras equivalentes baseadas nos quantis
get.breaks <- function(x, n){
  probs <- seq(0, 1, length.out = n + 1)
  quantile(x, probs, type = 1, names = FALSE)
}

# Funcao p/ rodar IDW
run.idw <- function(points, grid, scale, idp = 2){
  form <- as.formula(paste(scale, "~ 1"))
  gstat::idw(formula = form, locations = points, newdata = grid, idp = idp)
}

bbox.br <- st_bbox(br.sf)
bbox.br[3] <- -35

scale.sf <- scale.df |> 
  filter(gauge_code != "Regional") |> 
  left_join(
    info.df |> select(gauge_code, state, long, lat),
    by = "gauge_code"
  ) |> 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(epsg))

scale.sp <- as_Spatial(scale.sf)

# Criar camada 'SpatRaster'
grid <- st_make_grid(br.sf, cellsize = 0.1, what = "centers", square = TRUE) |> 
  st_intersection(br.sf) |> 
  st_as_sf()

# Interpolação IDW
scale.breaks <- get.breaks(scale.sf$scale, 5)
scale.idw <- scale.sf |> 
  group_split(which_moment) |> 
  map_dfr(function(df){
    q <- unique(df$which_moment)
    run.idw(df, grid, "scale", idp = 2) |> 
      mutate(
        which_moment = q, 
        interval = cut(var1.pred, scale.breaks, right = FALSE, include.lowest = TRUE)
      )
  }); invisible(gc())

coords <- st_coordinates(scale.idw)
scale.r <- scale.idw |> 
  st_drop_geometry() |> 
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  )

scale.colors <- RColorBrewer::brewer.pal(length(scale.breaks), "YlOrRd")
switch.legend <- c(FALSE, FALSE, TRUE)
subtitle <- c("a", "b", "c")

### 'H' por ordem (q = 1, 1.5, 2, 2.5, 3) ---------------------------------

# Interpolação
scale.ls <- split(scale.r, scale.r$which_moment)
scale.plot.ls <- lapply(seq_along(scale.ls), function(q){
  df <- scale.ls[[q]]
  ggplot() +
    # geom_sf(data = df, aes(color = interval), size = 0.05) +
    geom_tile(data = df, aes(x = x, y = y, fill = interval), show.legend = switch.legend[q]) +
    geom_sf(data = subdaily.sf, color = "grey10", shape = 3, size = 0.3, stroke = 0.08, inherit.aes = FALSE) +
    geom_sf(data = uf.sf, color = "grey30", fill = NA, linewidth = 0.1) +
    geom_sf(data = br.sf, color = "black", fill = NA, linewidth = 0.4) +
    scale_fill_manual(name = expression(H), values = scale.colors) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    annotation_scale(
      location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
      bar_cols = c("black", "aliceblue"), line_width = 0.4
    ) +
    annotation_north_arrow(
      location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
      style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "white"))
    ) +
    coord_sf(crs = st_crs(epsg), xlim = bbox.br[c(1,3)], ylim = bbox.br[c(2,4)]) +
    annotate(
      geom = "text", x = -Inf, y = Inf, label = sprintf("(%s) q = %g", subtitle[q], q),
      family = font, size = textsize, hjust = -0.2, vjust = 2
    ) +
    labs(x = NULL, y = NULL) +
    theme(
      plot.margin = margin(t = 10, r = 1, l = 1),
      panel.background = element_rect(fill = "white"),
      legend.key.spacing.y = unit(0, "pt"),
      legend.key.size = unit(fontsize*0.9, "pt"),
      legend.margin = margin(t = -5, l = -2),
      legend.position = c(0.03, 0.02),
      legend.justification = c(0, 0)
    )
})

layout <- "
  1111
  #22#
"

plot.scale.q <- wrap_plots(scale.plot.ls[1:2]) / wrap_plots(scale.plot.ls[3]) + 
  plot_layout(guides = "keep", design = layout); plot.scale.q

ggsave(
  file.path(path.plot, "plot_scale_q_subdiario.png"), plot.scale.q, 
  width = 15, height = 16, units = "cm", dpi = 600
)

