
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
path.imax <- file.path("Relatório 2", "Scripts_Ajuste_dGEV", "resultados", "subdiarios", "dataframes")
path.fit <- file.path(path.root, "resultados/subdiarios/fit/profile_likelihood")
path.plot <- file.path(path.root, "figuras")


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais
imax.df <- arrow::read_parquet(file.path(path.imax, "df_imax_subdiario_ano_civil_rs.pqt")) |> 
  filter(na_prct_yr <= 0.2) |> 
  group_by(gauge_code) |> 
  mutate(n_years = n_distinct(wateryear)) |> 
  filter(n_years >= 8) |> 
  ungroup()

# Informações sobre as estações subdiárias do RS
info.df <- arrow::read_parquet(file.path(path.imax, "df_subdiario_info.pqt")) |> 
  filter(gauge_code %in% unique(imax.df$gauge_code))

# Ajuste perfil de verossimilhanca
par.df <- bind_rows(
  arrow::read_parquet(file.path(path.fit, "df_subdiario_dgev_rs_par_mle.pqt")),
  arrow::read_parquet(file.path(path.fit, "df_subdiario_dgev_rs_par_rgmle.pqt"))
)

# Nomear modelos
which.prior <- unique(par.df$prior)
prior.names <- c(
  mle = NA_character_,         # s/ distribuicao a priori
  gmle = "B(0.1, 0.122)",      # distribuicao a priori de Martins e Stedinger (2000)
  rgmle = "B(-0.0178, 0.0627)" # distribuicao a priori regional c/ modelo GLS
) # versão c/ sinais trocados da média
prior.names <- setNames(prior.names[c("mle", "rgmle")], which.prior)
prior.labels <- setNames(c("MLE", "RGMLE"), prior.names)

# Atribuir nomes dos modelos a coluna prior
par.df <- par.df |> 
  mutate(
    prior = prior.names[prior],
    prior_label = prior.labels[prior],
    par = ifelse(parname == "xi", -par, par)
  )

# Camadas geoespaciais
sa.sf <- rnaturalearth::ne_countries(continent = "south america", scale = "large")
br.sf <- geobr::read_country(year = 2025)
rs.sf <- geobr::read_state(year = 2025, code_state = "RS")
subdaily.sf <- info.df |> 
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
  filter(prior_label == "RGMLE") |> 
  left_join(subdaily.sf, by = "gauge_code") |> 
  select(gauge_code, parname, par, ci_lower, ci_upper, prior, prior_label, geometry) |> 
  st_as_sf()

# Criar uma camada SpatRaster
bbox.rs <- st_bbox(rs.sf)
buffer <- 5
grid <- st_make_grid(rs.sf, cellsize = 0.01, what = "centers", square = TRUE) |> 
  st_as_sf(crs = st_crs(rs.sf)) |> 
  st_intersection(rs.sf)

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
  lapply(function(x) get.breaks(x, 4))

# Redefinir intervalos para xi
# par.breaks[["xi"]] <- c(-0.46, -0.2, -0.05, 0.05, 0.2, 0.49) # MLE
par.breaks[["xi"]] <- c(-0.24, -0.1, -0.05, 0.05, 0.1, 0.25) # RGMLE

# Interpolar
par.grid.sf <- par.sf |> 
  filter(parname %in% c("xi", "H", "theta")) |> 
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

subtitle <- c("a", "b", "c")
par.order <- c("xi", "H", "theta")
par.labels <- c("ξ<sub>j</sub>", "H<sub>j</sub>", "θ")
par.ls <- split(par.r.df, par.r.df$parname)[par.order]

plot.ls <- lapply(seq_along(par.ls), function(i){
  
  df <- par.ls[[i]]
  par.name <- unique(df$parname)
  ggplot() +
    geom_tile(aes(x = long, y = lat, fill = interval), df) +
    geom_sf(data = subdaily.sf, color = "grey10", shape = 3, size = 1, inherit.aes = FALSE) +
    geom_sf(data = rs.sf, color = "black", fill = NA, linewidth = 0.3, inherit.aes = FALSE) +
    scale_fill_manual(name = NULL, values = colors.ls[[par.name]]) +
    coord_sf(crs = st_crs(4674), xlim = bbox.rs[c(1,3)], ylim = bbox.rs[c(2,4)], expand = TRUE, clip = "off") +
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

layout <- "
  1111
  #22#
"

plot.dgev.par <- 
  wrap_plots(plot.ls[1:2]) / wrap_plots(plot.ls[3]) + 
  plot_layout(design = layout); plot.dgev.par

ggsave(
  file.path(path.plot, "mapa_dgev_par.png"), plot.dgev.par, 
  width = 15, height = 15, units = "cm", dpi = 600
)
