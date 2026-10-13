
# Análise de falhas do conjunto de estações subdiárias usando os dados em
# "Relatório 2\Scripts_Ajuste_dGEV\resultados\subdiarios\dataframes\df_imax_subdiario.pqt"


# PACOTES ----------------------------------------------------------------

pkg <- c(
  "tidyverse", "arrow",
  "sf", "geobr",
  "ggplot2", "ggspatial", "rnaturalearth", "patchwork"
)

missing.pkg <- pkg[!pkg %in% installed.packages()[,"Package"]]
if(!require(pak)) install.packages("pak")
if(length(missing.pkg)) pak::pkg_install(missing.pkg)

library(tidyverse)
library(patchwork)
library(arrow)
library(sf)
library(geobr)
library(ggspatial)
library(rnaturalearth)


# VISUALIZAÇÃO -----------------------------------------------------------

fontsize <- 8   # pts
textsize <- (fontsize - 2)/.pt # p/ geom_text, annotate, etc
fig.width <- 15 # cm
ggplot2::set_theme(
  theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(-5.5, 2, 2, 2),
      legend.key.size = unit(fontsize, "pt"),
      text = element_text(color = "black", size = fontsize),
      axis.text = element_text(color = "black"),
      strip.text = element_text(hjust = 0, size = fontsize - 2, margin = margin(2, 2, 2, 2), color = "black"),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5)
    )
)


# CAMINHOS ---------------------------------------------------------------

path.imax <- "Relatório 2/Scripts_Ajuste_dGEV/resultados/subdiarios/dataframes/df_imax_subdiario.pqt"
path.plot <- "Relatório 2/Scripts_Ajuste_dGEV/figuras"


# LER DADOS --------------------------------------------------------------

# Intensidades máximas anuais
if(!file.exists(file.path(path.imax))) stop("Intensidades maximas anuais nao encontradas")
imax.df <- arrow::read_parquet(file.path(path.imax))

# Dados geoespaciais
world.sf <- rnaturalearth::ne_countries(continent = "south america", scale = "large")
br.sf <- geobr::read_country(year = 2025)
uf.sf <- geobr::read_state(year = 2025)


# ANÁLISE QC -------------------------------------------------------------

quality.df <- imax.df |> 
  select(gauge_code, wateryear, na_prct_yr) |> 
  group_by(gauge_code, wateryear) |> 
  summarise(na_prct_yr = first(na_prct_yr)*100, .groups = "drop")

limites.na <- seq(10, 40, 10)
limites.min.years <- seq_len(10)

grid <- expand.grid(na_accept = limites.na, min_years = limites.min.years) |> 
  mutate(n_gauges = purrr::map2_int(.x = na_accept, .y = min_years, .f = ~{
    quality.df |> 
      filter(na_prct_yr <= .x) |> # manter estações c/ no máximo `na.accept` falhas anuais
      count(gauge_code) |>     # contar qtas estações
      filter(n >= .y) |>       # filtrar estações c/ pelo menos `min_years` anos
      nrow()                   # contar qtas estações
  }))

# Limiares escolhidos
min.years <- 8
na.accept <- 20
plot.anos.estacao <- grid |> 
  mutate(
    na_accept = factor(na_accept, levels = limites.na),
    fill = if_else(min_years == min.years & na_accept == na.accept, "escolhido", NA_character_)
  ) |> 
  ggplot(aes(x = min_years, y = n_gauges, fill = fill)) +
  facet_wrap(~na_accept, nrow = length(limites.na), labeller = labeller(na_accept = function(x) paste0(x, "% falhas"))) +
  geom_col(color = "black", show.legend = FALSE, linewidth = 0.2, width = 0.8) +
  geom_text(aes(label = if_else(!is.na(fill), "▼", "")), vjust = -1.2, size = 3.5, color = "grey10") +
  geom_text(aes(label = n_gauges), vjust = -0.5, size = textsize*1.05) +
  scale_x_continuous(breaks = seq(min(grid$min_years), max(grid$min_years), by = 1)) +
  scale_y_continuous(limits = c(0, max(grid$n_gauges)*1.1)) +
  scale_fill_manual(values = "grey10", na.value = "grey80") +
  labs(x = "Número de anos por estação", y = "Número de estações"); plot.anos.estacao

ggsave(
  file.path(path.plot, "analise_qc_subdiario.png"), plot.anos.estacao, 
  width = fig.width, height = 15, units = "cm", dpi = 300
)


# MAPA QC ----------------------------------------------------------------

valid.gauges <- imax.df |> 
  filter(na_prct_yr <= na.accept/100) |> 
  group_by(gauge_code) |> 
  filter(n_distinct(wateryear) >= min.years) |> 
  ungroup() |> 
  select(gauge_code, wateryear, d, imax) |> 
  mutate(valid_qc = "Válidas") |> 
  slice_head(n = 1, by = gauge_code) |> 
  pull(valid_qc, name = gauge_code)

length(valid.gauges)

# Pontos estações subdiárias
gauges.sf <- imax.df |> 
  slice_head(n = 1, by = gauge_code) |> 
  select(gauge_code, d, wateryear, na_prct_yr, estado, start_month, long, lat) |> 
  mutate(
    n_years = n_distinct(wateryear),
    valid_qc = valid.gauges[gauge_code],
    valid_qc = factor(if_else(is.na(valid_qc), "Inválidas", valid_qc), levels = c("Inválidas", "Válidas")),
    .by = gauge_code
  ) |> 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(4764))

gauges.sf |> 
  group_by(valid_qc) |> 
  count()

invisible(gc())


## BRASIL E ESTAÇÕES SUBDIÁRIAS ------------------------------------------

color.qc <- setNames(c("grey70", "dodgerblue3"), unique(gauges.sf$valid_qc))
size.qc <- setNames(c(0.05, 0.1), unique(gauges.sf$valid_qc))

bbox.br <- st_bbox(br.sf)

# Mapa
plot.map.br <- ggplot() +
  geom_sf(data = world.sf, color = "lightblue", fill = "transparent", linewidth = 3) +
  geom_sf(data = world.sf, color = "grey30", fill = "grey90") +
  geom_sf(data = uf.sf, color = "grey10", fill = "lightgreen", linewidth = 0.2) +
  geom_sf(aes(color = "Inválidas"), gauges.sf[gauges.sf$valid_qc == "Inválidas",], size = 0.05, alpha = 0.8) +
  geom_sf(aes(color = "Válidas"), gauges.sf[gauges.sf$valid_qc == "Válidas",], size = 0.1, alpha = 0.8) +
  guides(color = guide_legend(label.hjust = 0, override.aes = list(size = 1.2))) +
  coord_sf(crs = st_crs(4764), xlim = bbox.br[c(1,3)], ylim = bbox.br[c(2,4)]) +
  scale_color_manual(name = "Controle de qualidade", values = color.qc) +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_cex = 0.6,
    bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),  which_north = "true",
    style = north_arrow_nautical(text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  theme(
    # plot.margin = margin(t = 10),
    panel.background = element_rect(fill = "aliceblue"),
    legend.text = element_text(size = fontsize),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    # legend.key.spacing.x = unit(fontsize*2.5, "pt"),
    # legend.margin = margin(t = -5, l = -2),
    # # legend.spacing.y = unit(5, "pt"),
    # # legend.spacing.x = unit(-4, "pt"),
    # # legend.position = c(0.985, 0.12),
    # # legend.justification = c(1, 1)
  ); plot.map.br

ggsave(
  file.path(path.plot, "mapa_analise_qc_subdiario.png"), plot.map.br,
  width = fig.width, height = 15, units = "cm", dpi = 300
)

# Gráfico de barras
uf.df <- gauges.sf |> 
  st_drop_geometry() |> 
  select(gauge_code, estado, valid_qc) |> 
  slice_head(n = 1, by = gauge_code) |> 
  group_by(estado) |> 
  count(valid_qc) |> 
  print(n = Inf)

uf.order <- uf.df |> 
  filter(valid_qc == "Inválidas") |> 
  arrange(desc(n)) |> 
  pull(estado)

plot.gauges.by.uf <- uf.df |> 
  mutate(
    estado = factor(estado, uf.order),
    label_n = sprintf("%-3g", n)
  ) |> 
  arrange(estado) |> 
  ggplot(aes(y = estado, x = n, fill = valid_qc)) +
  geom_col(position = "stack") +
  geom_text(aes(label = label_n), hjust = -0.2, position = "stack", size = textsize) +
  scale_fill_manual(name = "Controle de qualidade", values = color.qc) +
  scale_x_continuous(name = "Número de estações", expand = expansion(mult = c(0.01, 0.08))) +
  labs(y = NULL) +
  theme(
    axis.ticks.y = element_line(linewidth = 0.3),
    legend.position = c(1,1),
    legend.justification = c(1,1)
  ); plot.gauges.by.uf

# plot.qc <- 
#   plot.map.br +
#   plot.gauges.by.uf +
#   plot_layout(ncol = 1, guides = "collect"); plot.qc

ggsave(
  file.path(path.plot, "qc_por_estado.png"), plot.gauges.by.uf,
  width = fig.width, height = 12, units = "cm", dpi = 300
)

