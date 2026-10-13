
rm(list = ls()); invisible(gc())


# PACOTES ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggtext)
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
textsize = (fontsize - 2)/.pt
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
D <- max(mom.df$d)

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

# Selecionar estações mais longas por estado
longest.gauges <- info.df |> 
  mutate(n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)]) |> 
  arrange(state, desc(n_years)) |> 
  slice_head(n = 1, by = "state")


# VISUALIZACAO -----------------------------------------------------------

# Filtrar estações mais longas por estado
mom.gauge.df <- mom.df |> 
  rename(n_years = n_year) |> 
  filter(gauge_code %in% longest.gauges$gauge_code) |> 
  left_join(
    scale.df |> 
      filter(gauge_code %in% longest.gauges$gauge_code) |> 
      select(gauge_code, which_moment, scale, rsquared),
    by = c("gauge_code", "which_moment"),
    relationship = "many-to-many"
  ) |> 
  mutate(
    duration_ratio = d/D,
    state = imax.df$state[match(gauge_code, imax.df$gauge_code)],
    gauge_label = sprintf("%s - %s | %d anos", state, gauge_code, n_years),
    scale_label = sprintf("q = %g | H = %0.4f | R² = %0.4f", which_moment, scale, rsquared),
  ) |> 
  print()

vlines <- tibble(
  x = c(1, 24)/D,
  label = sprintf("d = %s h", c(1, 24))
)

moment.label <- tibble(
  which_moment = unique(mom.df$which_moment),
  label = paste("q =", which_moment)
)

ufs <- sort(unique(mom.gauge.df$state))
walk(ufs, function(uf){

  # Figura principal
  p1 <- mom.gauge.df |> 
    filter(state == uf) |> 
    ggplot(aes(x = duration_ratio, y = mom)) +
    facet_wrap(~gauge_label, ncol = 1, scales = "free_y") +  
    geom_vline(aes(xintercept = x), vlines, linewidth = 0.2, linetype = "dashed", color = "grey80") +
    geom_text(
      aes(label = label, x = x, y = 0), vlines, family = font, size = textsize, inherit.aes = FALSE,
      vjust = -1, hjust = 1.15
    ) +
    geom_smooth(
      aes(group = which_moment), color = "grey20", linetype = "solid", alpha = 0.8,
      method = "lm", formula = y~x, se = FALSE, linewidth = 0.3
    ) +
    geom_point(aes(shape = scale_label), color = "black", fill = "black", size = 1.2, stroke = 0) +
    scale_x_log10(name = "d/D", limits = c(1/6/24, 1)) +
    scale_y_log10(name = "E[I<sub>d</sub><sup>q</sup>]", labels = scales::label_number(big.mark = ".")) +
    scale_shape_manual(name = NULL, values = c(21, 24, 22)) +
    theme(
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      axis.title.y = ggtext::element_markdown()
    ) +
    guides(
      shape = guide_legend(override.aes = list(alpha = 1, color = "black", fill = "black", size = 1.4))
    )
  
  p2 <- scale.df |> 
    mutate(state = imax.df$state[match(gauge_code, imax.df$gauge_code)]) |> 
    filter(state == uf) |> 
    ggplot(aes(x = rsquared)) +
    facet_wrap(~which_moment, ncol = 1) +
    geom_histogram(
      aes(y = after_stat(count/sum(count))),
      bins = 30, fill = "grey90", color = "black", linewidth = 0.2
    ) +
    geom_text(
      aes(label = label, x = Inf, y = Inf), moment.label, family = font, size = textsize,
      inherit.aes = FALSE, hjust = 1.2, vjust = 1.6
    ) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(x = "R²", y = NULL) +
    theme(
      legend.position = "none",
      strip.text = element_blank(),
      axis.title.y = ggtext::element_markdown(),
      plot.margin = margin(0.5, 0.5, 0.5, 5.5)
    )
    
  scale.invariance.plot <- 
    wrap_plots(p1) +
    wrap_plots(p2) +
    plot_layout(widths = c(1, 0.5), guides = "keep")

  ggsave(
    file.path(path.plot, sprintf("plot_invariancia_escala_momentos_%s.png", uf)), scale.invariance.plot,
    width = 15, height = 12, units = "cm", dpi = 600
  )

})
