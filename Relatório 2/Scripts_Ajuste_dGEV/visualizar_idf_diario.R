
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
path.fit <- file.path(path.root, "resultados/diarios/fit/profile_likelihood/")
path.plot <- file.path(path.root, "figuras")


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais
imax.df <- arrow::read_parquet(file.path(path.imax, "df_imax_uniplu_hidro.parquet"))

# Ajuste perfil de verossimilhanca
par.df <- arrow::read_parquet(file.path(path.fit, "df_diario_dgev_par.pqt"))
idf.df <- arrow::read_parquet(file.path(path.fit, "df_diario_dgev_idf.pqt"))


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
tr.cores <- setNames(c("#d55e00", "#e69200", "#e69F00", "#009e73", "#0072b2", "#cc79a7", "#7b3294"), tr.todos)

# Estações
gg.br <- c(
  "739006", # CE
  "1842004", # MG
  "735041",  # PE
  "2549006", # PR
  "535011",  # RN
  "2849002"  # SC
)

gauge.labels <- imax.df |> 
  filter(gauge_code %in% gg.br) |> 
  slice_head(n = 1, by = gauge_code) |> 
  mutate(
    gauge_code = factor(gauge_code, gg.br),
    gauge_label = sprintf("%s - %s", estado, gauge_code)
  ) |> 
  arrange(gauge_code) |> 
  pull(gauge_label, name = gauge_code)

par.rs <- par.df |> 
  filter(gauge_code %in% gg.br) |>  
  left_join(
    imax.df |> slice_head(n = 1, by = gauge_code) |> select(gauge_code, estado),
    by = "gauge_code",
  ) |> 
  mutate(
    gauge_label = sprintf("%s - %s", estado, gauge_code),
    gauge_label = factor(gauge_label, gauge.labels)
  )

idf.rs <- idf.df |> 
  filter(gauge_code %in% gg.br) |> 
  left_join(
    imax.df |> slice_head(n = 1, by = gauge_code) |> select(gauge_code, estado),
    by = "gauge_code",
  ) |> 
  mutate(
    gauge_label = sprintf("%s - %s", estado, gauge_code),
    gauge_label = factor(gauge_label, gauge.labels)
  )
  

# Todos os tempos de retorno, s/ intervalos de confiança
idf.plot.tr.todos <- idf.rs |>
  filter(return_period %in% tr.todos) |> 
  ggplot(aes(x = d, color = factor(return_period))) +
  facet_wrap(~gauge_label, ncol = 2, scales = "free") +
  geom_line(aes(y = q), linewidth = 0.3) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("ξ = %7.3f", -par)), filter(par.rs, parname == "xi"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", par)), filter(par.rs, parname == "H"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + textsize*0.7, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ = %7.3f", par)), filter(par.rs, parname == "theta"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 2*textsize*0.7, family = font, size = textsize
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
    name = "Tempo de retorno (anos)", values = tr.cores
  ); idf.plot.tr.todos

ggsave(
  file.path(path.plot, "idf_diario_tr_todos.png"), idf.plot.tr.todos,
  width = 15, height = 18, units = "cm", dpi = 600
)

# Tempos de retorno 1
idf.plot.tr1 <- idf.rs |> 
  filter(return_period %in% tr1) |> 
  ggplot(aes(x = d, color = factor(return_period), fill = factor(return_period))) +
  facet_wrap(~gauge_label, ncol = 2, scales = "free") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linewidth = 0, show.legend = FALSE) +
  geom_line(aes(y = q), linewidth = 0.3) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("ξ = %7.3f", -par)), filter(par.rs, parname == "xi"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", par)), filter(par.rs, parname == "H"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + textsize*0.7, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ = %7.3f", par)), filter(par.rs, parname == "theta"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 2*textsize*0.7, family = font, size = textsize
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
  file.path(path.plot, "idf_diario_tr1.png"), idf.plot.tr1, 
  width = 15, height = 18, units = "cm", dpi = 600
)

# Tempos de retorno 2
idf.plot.tr2 <- idf.rs |> 
  filter(return_period %in% tr2) |> 
  ggplot(aes(x = d, color = factor(return_period), fill = factor(return_period))) +
  facet_wrap(~gauge_label, ncol = 2, scales = "free") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linewidth = 0, show.legend = FALSE) +
  geom_line(aes(y = q), linewidth = 0.3) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("ξ = %7.3f", -par)), filter(par.rs, parname == "xi"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", par)), filter(par.rs, parname == "H"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + textsize*0.7, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ = %7.3f", par)), filter(par.rs, parname == "theta"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 2*textsize*0.7, family = font, size = textsize
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
  file.path(path.plot, "idf_diario_tr2.png"), idf.plot.tr2, 
  width = 15, height = 18, units = "cm", dpi = 600
)
