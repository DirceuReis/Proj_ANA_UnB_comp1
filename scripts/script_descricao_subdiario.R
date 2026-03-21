# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  pacman,
  tidyverse,
  arrow,
  lmom,
  patchwork, # combinar gráficos
  xtable,    # printa tabelas latex
  boot       # bootstrap
)

font_size <- 10

# LER DADOS ---------------------------------------------------------------

# Funções
source("scripts/funcoes/fun_boot_lratio.R") # estimar intervalos de confiança c/ bootstrap p/ razões de l-momentos

# Dados
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")
df.info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.pqt")
# df_imax <- df_imax[df_imax$gauge_code != "SBPA",]


# QUALIDADE DOS DADOS -----------------------------------------------------

# Manter somente estação, falhas e anos
df.quality <- df.imax %>% 
  mutate(responsible = df.info$responsible[match(gauge_code, df.info$gauge_code)]) %>% 
  select(gauge_code, responsible, year, na_prct) %>% 
  group_by(gauge_code, responsible, year) %>% 
  reframe(na_prct = first(na_prct)*100)

# Limites de falha e mínimo de anos
thresholds <- seq(10, 40, 10) # limites de falha, de 10 a 40%
min.years <- 1:10             # ao menos 1 ano até ao menos 10 anos de dados completos

# Criar tbl_df c/ possibilidades de falha e anos
grid <-
  expand_grid(threshold = thresholds, min_years = min.years) %>% # combinações
  mutate(n_gauges = purrr::map2_int(.x = threshold, .y = min_years, ~{ # analisar qtas estações p/ cada par
    df.quality %>% 
      filter(na_prct <= .x) %>% # filtrar estações c/ no máximo 'threshold' falhas
      count(gauge_code) %>%     # contar qtas estações
      filter(n >= .y) %>%       # filtrar estações c/ no mínimo 'min_years' anos
      nrow()                    # contar qtas estações
  })) # fim 'grid'

resp_names <- paste(unique(df.info$responsible), collapse = ", ")

# Visualização
color <- c("green4", "lightgreen", "orange", "#F10000")
plot.anos.estacao <- ({
  grid %>% 
    mutate(threshold = factor(threshold, levels = thresholds)) %>% # ordenar valores de 'threshold'
    ggplot(aes(x = min_years, y = n_gauges, fill = threshold)) +
    geom_col() +
    geom_text(aes(label = n_gauges), vjust = -0.3, size = 3, family = "serif") +
    facet_wrap(~threshold, nrow = 4, labeller = labeller(threshold = function(x) paste0(x, "% falhas"))) +
    scale_x_continuous(breaks = seq(min(grid$min_years), max(grid$min_years), by = 1)) +
    scale_y_continuous(limits = c(0, max(grid$n_gauges) + 50)) +
    scale_fill_manual(values = color) +
    labs(x = "Número de anos por estação",
         y = "Número de estações",
         caption =  paste0("Responsáveis: ", resp_names)) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"),
          legend.position =  "none",
          plot.caption = element_text(size = 8),
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 12, family = "serif"))
}); plot.anos.estacao # fim 'plot.anos.estacao'

ggsave(filename = "figuras/plot_anos_estacao.png", plot = plot.anos.estacao,
       width = 16, height = 16, units = "cm", dpi = 200)


# FORMA DOS DADOS ---------------------------------------------------------

# Filtros
na.accept <- 0.2
min.years <- 8
which.durations <- c(1/6, 0.25, 1, 6, 12, 24, 48, 72, 120, 240)

# L-momentos
df_est <- df_imax %>% 
  filter(na_prct <= na.accept,
         d %in% which.durations) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(year) >= min.years) %>% 
  ungroup() %>% 
  group_by(gauge_code, d) %>%             # calcular estatísticas: 
  reframe(n_year = n_distinct(year),
          tau3 = lmom::samlmu(imax)[[3]], # L-assimetria
          tau4 = lmom::samlmu(imax)[[4]]) # L-curtose
  
# Tabela resumo p/ tau3
df_summary_tau3 <- df_est %>% 
  na.omit() %>% 
  group_by(d) %>% 
  reframe(min = min(tau3),
          q0.25 = quantile(tau3, 0.25),
          mean = mean(tau3, na.rm = TRUE),
          q0.50 = quantile(tau3, 0.50),
          q0.75 = quantile(tau3, 0.75),
          q0.90 = quantile(tau3, 0.90),
          max = max(tau3),
          sd = sd(tau3))

# Tabela resumo p/ tau4
df_summary_tau4 <- df_est %>% 
  na.omit() %>% 
  group_by(d) %>% 
  reframe(min = min(tau4),
          q0.25 = quantile(tau4, 0.25),
          mean = mean(tau4, na.rm = TRUE),
          q0.50 = quantile(tau4, 0.50),
          q0.75 = quantile(tau4, 0.75),
          q0.90 = quantile(tau4, 0.90),
          max = max(tau4),
          sd = sd(tau4))

# Gerar resultados p/ latex
print(xtable(x = df_summary_tau3, digits = 3), include.rownames = FALSE)
print(xtable(x = df_summary_tau4, digits = 3), include.rownames = FALSE)


# BOOTSTRAP L-CURTOSE -----------------------------------------------------

data <- df_imax %>% 
  filter(na_prct <= na.accept) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(year) >= min.years)

# REFAZER usando o corpo da função de bootstrap por blocos

# Realizar um procedimento bootstrap p/ estimar intervalos de confiança
# p/ estimativas de L-curtose e do coeficiente angular da regressão linear
df_lratio_ci <- fun_boot_lratio(data = data,
                                rep = 2000,
                                which_lratio = 4,
                                signf = 0.05,
                                ci_type = "perc",
                                col_names = c("imax", "d", "gauge_code", "na_prct", "year"))

# Adicionar intervalos de confiança em 'df_est'
df_est <- merge(df_est, df_lratio_ci[, c("gauge_code", "ds", "ci_l", "ci_u")], by = c("gauge_code", "ds"), all.x = TRUE)

# Somente L-curtose
df_slopes <- df_est %>% 
  filter(gauge_code != "SBPA") %>% 
  filter(n() > 1, !is.na(tau4), n_distinct(ds) > 2) %>% 
  select(gauge_code, tau4, ds) %>% 
  group_by(gauge_code) %>% 
  reframe(slope = coef(lm(tau4~ds))[2],
          p_value = summary(lm(tau4~ds))$coefficients[2,4]) %>% 
  mutate(signf = p_value < 0.05)

# Visualização
# Boxplot L-assimetria e L-curtose
plot_tau3 <- df_est %>% 
  mutate(ds = as.factor(ds*60)) %>% 
  ggplot(aes(x = ds, y = tau3)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill = "orange", linewidth = 0.5, outlier.shape = 4, outlier.size = 2) +
  theme_minimal() +
  labs(x = "Durações [min]", y = "L-assimetria") +
  theme(legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "Times New Roman", color = "black", size = font_size)); plot_tau3

plot_tau4 <- df_est %>% 
  mutate(ds = as.factor(ds*60)) %>% 
  ggplot(aes(x = ds, y = tau4)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(fill = "darkorange2", linewidth = 0.5, outlier.shape = 4, outlier.size = 2) +
  theme_minimal() +
  labs(x = "Durações [min]", y = "L-curtose") +
  theme(legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "Times New Roman", color = "black", size = font_size)); plot_tau4

plot_lratios <- (plot_tau3 / plot_tau4); plot_lratios

ggsave(plot = plot_lratios, filename = "figuras/plot_lratios.png",
       width = 160, height = 200, units = "mm", dpi = 200)

# L-curtose (tau4) vs. duração (ds) e regressão linear
plot_tau4_lm <- df_est %>% 
  ggplot(aes(x = ds, y = tau4)) +
  # facet_wrap(~gauge_code, nrow = 4, ncol = 11) +
  facet_wrap(~gauge_code, nrow = 11, ncol = 4) +
  geom_hline(yintercept = 0.4, color = "red", linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.4, color = "grey40", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", linewidth = 0.6) +
  geom_point(pch = 8) +
  coord_cartesian(ylim = c(-1,1)) + # coord_cartesian no lugar do scale_y_continuous não corta os elementos fora do plot
  theme_minimal() +
  labs(x = "Durações [h]", y = "L-curtose") +
  theme(legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "Times New Roman", color = "black", size = font_size)); plot_tau4_lm

# ggsave(plot = plot_tau4_lm, filename = "figuras/plot_tau4_lm.png",
#        width = 210, height = 297, dpi = 300, units = "mm",
#        # device = cairo_pdf,
#        ) # cairo_pdf plota textos c/ fontes que não são padrão


# SAZONALIDADE ------------------------------------------------------------

# Usar um plot de rosa dos ventos
na.accept <- 0.2    # percentual de falhas (0.2 -> 20%)
min.years <- 8      # mínimo de anos na série
which.moment <- 1:3 # qual momento usar de base p/ gerar figuras individuais
min.duration <- 3   # mínimo de durações diferentes necessário
which.duration <- c(1, 24)

# Alice sugeriu fazer a barra usando 'position_stack'
# e a cor indicando "intervalos" de intensidade
# ou um heatmap
plot.month.frequency1 <- {
  df_imax %>%
    filter(na_prct <= na.accept,              # filtrar anos por porcentagem de falhas
           d %in% c(1, 2, 4, 6, 8, 12, 18, 24),
           # d %in% c(24, 48, 72, 140, 240),
    ) %>% # filtrar durações dentro do intervalo
    group_by(gauge_code) %>%                  # agrupar por estação
    filter(n_distinct(d) >= min.duration,     # só estações com pelo menos 'min.duration'
           n_distinct(year) >= min.years) %>% # so estações com pelo menos 'min.years'
    ungroup() %>% 
    mutate(month = factor(lubridate::month(date, label = TRUE))) %>%
    # mutate(month = factor(lubridate::month(date))) %>%
    group_by(month, d) %>% 
    reframe(n_imax = n()) %>% 
    ggplot(aes(x = month, y = n_imax, fill = month)) +
    geom_col(color = "black", linewidth = 0.3, alpha = 0.8, width = 0.8, show.legend = FALSE) +
    facet_wrap(~d, ncol = 4, labeller = labeller(d = ~paste(.x, "h"))) +
    scale_y_continuous(expand = c(0,0)) +
    coord_radial() +
    labs(x = "", y = "Frequência de máximos anuais") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 8, margin = margin(b = 1)),
          axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.background = element_rect(color = "white"),
          # panel.border = element_rect(color = "black", fill = NA),
          text = element_text(size = 10, family = "serif")) +
    guides(theta = guide_axis_theta(angle = 0))
}; plot.month.frequency1

ggsave(plot = plot.month.frequency1, filename = "figuras/month_frequency_imax_days.png",
       width = 16, height = 5, units = "cm", dpi = 300)
