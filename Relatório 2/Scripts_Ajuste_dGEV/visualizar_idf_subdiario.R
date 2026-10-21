
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(arrow)
library(forcats)
library(ggplot2)
library(patchwork)


# CAMINHOS ----------------------------------------------------------------

path.root <- "Relatório 2/Scripts_Ajuste_dGEV/"
path.imax <- file.path("Relatório 2", "Scripts_Maximas", "resultados", "BR", "dataframes")
path.fit <- file.path(path.root, "resultados/subdiarios/fit/profile_likelihood/")
path.plot <- file.path(path.root, "figuras")


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais
imax.df <- arrow::read_parquet(file.path(path.imax, "imax_sem_iqr_subdiario_Uniplu.parquet"))

# Ajuste perfil de verossimilhanca
par.df <- arrow::read_parquet(file.path(path.fit, "df_subdiario_dgev_par.pqt"))
idf.df <- arrow::read_parquet(file.path(path.fit, "df_subdiario_dgev_idf.pqt"))


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


## IDF --------------------------------------------------------------------

# Tempos de retorno
tr.todos <- c(2, 5, 10, 25, 50, 100, 200)
tr1 <- c(10, 50, 100)
tr2 <- c(5, 25, 200)

# Cores
tr.cores <- setNames(c("#d55e00", "#e69F00", "#009e73", "#56b4e9", "#0072b2", "#cc79a7", "#7b3294"), tr.todos)

# Estações foram escolhidas buscando reunir 3 estações c/ d < 1 h e 3 c/ d > 1 h
# das mais longas disponíveis: 1 c/ xi approx 0, ourta c/ xi > 0 e outra c/ xi < 0
gg <- c("SBAR", "SBCF", "SBRP", "420800501A", "293330704A", "314610703A") # estações usadas pela Larissa

gauge.labels <- imax.df |> 
  filter(gauge_code %in% gg) |> 
  slice_head(n = 1, by = gauge_code) |> 
  mutate(
    gauge_code = factor(gauge_code, gg),
    gauge_label = sprintf("%s - %s", state, gauge_code)
  ) |> 
  arrange(gauge_code) |> 
  pull(gauge_label, name = gauge_code)


### Rio Grande do Sul -----------------------------------------------------

par.gg <- par.df |> 
  filter(gauge_code %in% gg) |>  
  left_join(
    imax.df |> slice_head(n = 1, by = gauge_code) |> select(gauge_code, state),
    by = "gauge_code",
  ) |> 
  mutate(
    gauge_label = sprintf("%s - %s", state, gauge_code),
    gauge_label = factor(gauge_label, gauge.labels)
  )

idf.gg <- idf.df |> 
  filter(gauge_code %in% gg) |> 
  left_join(
    imax.df |> slice_head(n = 1, by = gauge_code) |> select(gauge_code, state),
    by = "gauge_code",
  ) |> 
  mutate(
    gauge_label = sprintf("%s - %s", state, gauge_code),
    gauge_label = factor(gauge_label, gauge.labels)
  )

info.gg <- imax.df |> 
  filter(gauge_code %in% gg) |> 
  slice_head(n = 1, by = gauge_code) |> 
  select(-c(d, imax_mm_h)) |> 
  mutate(
    gauge_label = sprintf("%s - %s", state, gauge_code),
    gauge_label = factor(gauge_label, gauge.labels)
  )
  
textsize.fct <- 0.6

# Tempos de retorno 1
idf.plot.tr1 <- idf.gg |> 
  filter(return_period %in% tr1) |> 
  ggplot(aes(x = d, color = factor(return_period), fill = factor(return_period))) +
  facet_wrap(~gauge_label, ncol = 2, scales = "free_y") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linewidth = 0, show.legend = FALSE) +
  geom_line(aes(y = q), linewidth = 0.3) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype= "dashed", color = "grey") + 
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("ξ = %7.3f", -par)), filter(par.gg, parname == "xi"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", par)), filter(par.gg, parname == "H"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + textsize*textsize.fct, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ = %7.3f", par)), filter(par.gg, parname == "theta"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 2*textsize*textsize.fct, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("N = %g anos", n_valid_years)), info.gg, inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 3*textsize*textsize.fct, family = font, size = textsize
  ) +
  scale_x_continuous(
    name = "Duração (h)",
    trans = "log10"
  ) +
  scale_y_continuous(
    name = "Intensidade (mm/h)",
    trans = "log10"
  ) +
  scale_color_manual(
    name = "Tempo de retorno (anos)", values = tr.cores,
    aesthetics = list("color", "fill")
  ); idf.plot.tr1

ggsave(
  file.path(path.plot, "idf_subdiario_tr1.png"), idf.plot.tr1, 
  width = 15, height = 18, units = "cm", dpi = 600
)

# Tempos de retorno 2
idf.plot.tr2 <- idf.gg |> 
  filter(return_period %in% tr2) |> 
  ggplot(aes(x = d, color = factor(return_period), fill = factor(return_period))) +
  facet_wrap(~gauge_label, ncol = 2, scales = "free_y") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linewidth = 0, show.legend = FALSE) +
  geom_line(aes(y = q), linewidth = 0.3) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype= "dashed", color = "grey") + 
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("ξ = %7.3f", -par)), filter(par.gg, parname == "xi"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", par)), filter(par.gg, parname == "H"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + textsize*textsize.fct, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ = %7.3f", par)), filter(par.gg, parname == "theta"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 2*textsize*textsize.fct, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("N = %g anos", n_valid_years)), info.gg, inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 3*textsize*textsize.fct, family = font, size = textsize
  ) +
  scale_x_continuous(
    name = "Duração (h)",
    trans = "log10"
  ) +
  scale_y_continuous(
    name = "Intensidade (mm/h)",
    trans = "log10"
  ) +
  scale_color_manual(
    name = "Tempo de retorno (anos)", values = tr.cores,
    aesthetics = list("color", "fill")
  ); idf.plot.tr2

ggsave(
  file.path(path.plot, "idf_subdiario_tr2.png"), idf.plot.tr2, 
  width = 15, height = 18, units = "cm", dpi = 600
)
