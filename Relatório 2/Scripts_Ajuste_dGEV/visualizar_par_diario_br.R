
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(arrow)
library(forcats)
library(ggplot2)
library(ggtext)
library(RColorBrewer)
library(patchwork)
library(geobr)
library(ggspatial)
library(sf)
library(terra)
library(gstat)
library(rnaturalearth)
library(tidyterra)


# CAMINHOS ----------------------------------------------------------------

path.root <- "Relatório 2/Scripts_Ajuste_dGEV"
path.imax <- file.path("Relatório 2", "Scripts_Maximas", "resultados", "BR", "dataframes", "df_imax_uniplu_hidro.parquet")
path.fit <- file.path(path.root, "resultados/diarios/fit/profile_likelihood")
path.plot <- file.path(path.root, "figuras")


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais
# imax.df <- arrow::read_parquet(path.imax)
imax.df <- arrow::read_parquet(path.imax) |> 
  filter(na_prct_yr <= 0.2) |> 
  group_by(gauge_code) |> 
  mutate(n_years = n_distinct(wateryear)) |> 
  filter(n_years >= 30) |> 
  ungroup()

gauges <- unique(imax.df$gauge_code)

# Ajuste perfil de verossimilhanca
par.df <- arrow::read_parquet(file.path(path.fit, "df_diario_dgev_par.pqt"))

# Nomear modelos
which.prior <- unique(par.df$prior)
prior.names <- c(gmle = "B(0.1, 0.122)") # versão c/ sinais trocados da média
prior.names <- setNames("gmle", which.prior)
prior.labels <- setNames("GMLE", prior.names)

# Atribuir nomes dos modelos a coluna prior
par.df <- par.df |> 
  filter(gauge_code %in% gauges) |> 
  mutate(
    prior = prior.names[prior],
    prior_label = prior.labels[prior],
    par = ifelse(parname == "xi", -par, par) # trocar o sinal do parâmetro de forma p/ manter consistência
  )                                          # ao longo do relatório

# Camadas geoespaciais
sa.sf <- rnaturalearth::ne_countries(continent = "south america", scale = "large")
br.sf <- geobr::read_country(year = 2025)
uf.sf <- geobr::read_state(year = 2025)
daily.sf <- imax.df |> 
  select(gauge_code, long, lat) |> 
  slice_head(n = 1, by = gauge_code) |> 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(4674)) |> 
  select(gauge_code, geometry)


# VISUALIZACAO ------------------------------------------------------------

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



## PARAMETRO DE FORMA -----------------------------------------------------

# Transformar em objeto sf
par.sf <- par.df |> 
  left_join(daily.sf, by = "gauge_code") |> 
  select(gauge_code, parname, par, ci_lower, ci_upper, prior, prior_label, geometry) |> 
  st_as_sf()

# Criar uma camada SpatRaster
bbox.br <- st_bbox(br.sf)
buffer <- 5
grid <- st_make_grid(br.sf, cellsize = 0.05, what = "centers", square = TRUE) |> 
  st_as_sf(crs = st_crs(br.sf)) |> 
  st_intersection(br.sf); gc()

# Função para interpolação IDW
run.idw <- function(points, grid, parname, idp = 2.0){
  formula <- as.formula(paste(parname, "~ 1"))
  gstat::idw(formula = formula, locations = points, newdata = grid, idp = idp)
}

# Função para definir intervalos
get.breaks <- function(x, n){
  probs <- seq(0, 1, length.out = n + 1)
  quantile(x, probs, type = 1, names = FALSE)
}

# Definir intervalos prévios
par.breaks <- split(par.sf$par, par.sf$parname) |> 
  lapply(function(x) get.breaks(x, 4)) |> 
  print()

# Redefinir intervalos para xi
par.breaks[["xi"]] <- c(-0.37, -0.10, -0.05, 0.05, 0.10, 0.46) # GMLE

# Interpolar
par.grid.sf <- par.sf |> 
  filter(parname %in% c("xi", "H")) |> 
  group_split(parname) |> 
  purrr::map_dfr(function(df){
    par.name <- unique(df$parname)
    run.idw(df, grid, "par", idp = 2) |> 
      mutate(
        parname = par.name,
        interval = cut(var1.pred, par.breaks[[par.name]], right = FALSE, include.lowest = TRUE)
      )
  }); invisible(gc())

coords <- st_coordinates(par.grid.sf)
par.r.df <- par.grid.sf |> 
  st_drop_geometry() |> 
  mutate(long = coords[,1], lat = coords[,2])

# Cores
colors.ls <- list(
  xi = c("#2166AC", "#67A9CF", "#F7F7F7", "#FDDBC7", "#B2182B"),
  H = RColorBrewer::brewer.pal(length(par.breaks$H), "YlOrRd"),
  theta = RColorBrewer::brewer.pal(length(par.breaks$theta), "GnBu")
)

# subtitle <- c("a", "b", "c")
par.order <- c("xi", "H")
par.labels <- c("(a) ξ<sub>j</sub>", "(b) H<sub>j</sub>")
par.ls <- split(par.r.df, par.r.df$parname)[par.order]

plot.ls <- lapply(seq_along(par.ls), function(i){
  
  df <- par.ls[[i]]
  par.name <- unique(df$parname)
  ggplot() +
    geom_tile(aes(x = long, y = lat, fill = interval), df) +
    geom_sf(data = uf.sf, color = "grey50", fill = NA, linewidth = 0.2, inherit.aes = FALSE) +
    geom_sf(data = br.sf, color = "black", fill = NA, linewidth = 0.3) +
    # geom_sf(data = daily.sf, color = "grey10", shape = 3, size = 1, inherit.aes = FALSE) +
    scale_fill_manual(name = NULL, values = colors.ls[[par.name]]) +
    coord_sf(crs = st_crs(4674), xlim = bbox.br[c(1,3)], ylim = bbox.br[c(2,4)], expand = TRUE, clip = "off") +
    annotation_scale(
      location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
      bar_cols = c("black", "white"), line_width = 0.4
    ) +
    annotation_north_arrow(
      location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
      pad_x = unit(0.1, "pt"), pad_y = unit(4, "pt"),
      style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "white"))
    ) +
    annotate(
      geom = "richtext", x = -Inf, y = Inf, label = sprintf("%s", par.labels[i]),
      family = font, size = textsize*1.2, hjust = -0.01, vjust = 1.1, fill = NA, label.color = NA
    ) +
    labs(x = "Longitude", y = "Latitude") +
    theme(
      legend.position = c(0, 0),
      legend.justification = c(0, 0),
      legend.key.spacing.y = unit(2, "pt"),
      legend.key.size = unit(10, "pt"),
      panel.border = element_rect(linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
  
}) # fim 'plot.ls'

# plot.dgev.par <- 
#   wrap_plots(plot.ls[[1]]) / wrap_plots(plot.ls[[2]]) +
#   plot_layout(axis_titles = "collect"); plot.dgev.par

# ggsave(
#   file.path(path.plot, "mapa_dgev_par_diario_br.png"), plot.dgev.par, 
#   width = 15, height = 20, units = "cm", dpi = 600
# )

ggsave(
  file.path(path.plot, "mapa_dgev_par_xi_diario_br.png"), plot.ls[[1]],
  width = 15, height = 13, units = "cm", dpi = 600
)

ggsave(
  file.path(path.plot, "mapa_dgev_par_H_diario_br.png"), plot.ls[[2]],
  width = 15, height = 13, units = "cm", dpi = 600
)
