# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); gc()

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  pacman,
  tidyverse,
  arrow,
  lubridate,
  lmom,
  ggplot2,   
  patchwork, # combinar gráficos
  xtable,    # printa tabelas latex
  boot       # bootstrap
)


# LER DADOS ---------------------------------------------------------------

# Funções
source("scripts/funcoes/fun_boot_lratio.R") # estimar intervalos de confiança c/ bootstrap p/ razões de l-momentos

# Dados
df_imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")


# ESTATÍSTICAS DESCRITIVAS ------------------------------------------------

# L-momentos
df_est <- df_imax %>% 
  na.omit() %>% 
  group_by(gauge_code, ds, time_steps) %>%    # calcular estatísticas: 
  reframe(tau3 = lmom::samlmu(imax)[[3]],     # L-assimetria
          tau4 = lmom::samlmu(imax)[[4]]) %>% # L-curtose
  filter(gauge_code != "SBPA")
  
# Tabela resumo p/ tau3
df_summary_tau3 <- df_est %>% 
  na.omit() %>% 
  group_by(ds) %>% 
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
  group_by(ds) %>% 
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

# Realizar um procedimento bootstrap p/ estimar intervalos de confiança
# p/ estimativas de L-curtose e do coeficiente angular da regressão linear
df_lratio_ci <- fun_boot_lratio(data = df_imax,
                                rep = 2000,
                                which_lratio = 4,
                                signf = 0.05,
                                ci_type = "norm",
                                col_names = c("imax", "ds", "gauge_code"))

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


# VISUALIZAÇÃO ------------------------------------------------------------

font_size <- 10

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

ggsave(plot = plot_tau4_lm, filename = "figuras/plot_tau4_lm.png",
       width = 210, height = 297, dpi = 300, units = "mm",
       # device = cairo_pdf,
       ) # cairo_pdf plota textos c/ fontes que não são padrão
