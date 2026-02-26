
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow)

# LER DADOS ---------------------------------------------------------------

df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")
df.scale.boot <- arrow::read_parquet(file = "base/gerados/scale_invariance/df_scale_block_boot.pqt")

# FUNÇÕES -----------------------------------------------------------------

source("scripts/funcoes/fun_scale_invariance.R")


# RODAR FUNÇÃO ------------------------------------------------------------

# Argumentos
na.accept <- 0.2    # percentual de falhas (0.2 -> 20%)
min.years <- 8      # mínimo de anos na série
which.moment <- 1:3 # qual momento usar de base p/ gerar figuras individuais
min.duration <- 3   # mínimo de durações diferentes necessário   

# Avaliar diferentes modelos
duration.intervals <- list(c(10/60, 10*24), # todas as durações
                           c(10, 60)/60,    # sub-horários
                           c(1, 24),        # horários
                           c(1, 10)*24)     # diários

# Regressão + coeficiente de escala + plots
scale.invariance <- lapply(X = duration.intervals, FUN = function(interval){
  
  res <- fun_scale_invariance(df.imax = df.imax,
                              na.accept = na.accept,
                              min.years = min.years,
                              which.moment = which.moment,
                              which.duration = interval,
                              min.duration = 3,
                              font.family = "serif",
                              font.size = 12,
                              plot.dim = c(5,4))
  
}); names(scale.invariance) <- c("all", "subhourly", "hourly", "daily")

# Extrair os momentos
scale.moments <- lapply(X =  scale.invariance, FUN = function(interval){
  interval$scale.moments
}); scale.moments <- bind_rows(scale.moments)

scale.coefficient <- lapply(X = scale.invariance, FUN = function(interval){
  interval$scale.coefficient
});  scale.coefficient <- bind_rows(scale.coefficient)


# MOMENTOS P/ UMA ESTAÇÃO -------------------------------------------------

# Avaliar estações com mais anos e durações e selecionar uma
gauges <- scale.moments %>%
  group_by(gauge_code, which_moment) %>%
  reframe(min_d = min(d),
          max_d = max(d),
          n_durations = n_distinct(d),
          n_year = first(n_year),
          duration_years = n_durations*n_year) %>% 
  arrange(desc(n_durations), desc(duration_years)) %>% 
  pull(gauge_code) %>% 
  unique()

gauge <- gauges[1]

duration.intervals <- unique(scale.coefficient$d)

# Momentos e expoente de escala local
df.mom.gauge <- 
  merge(x = scale.moments %>%
          filter(gauge_code == gauge),
        y = scale.coefficient %>%
          filter(gauge_code == gauge, d == duration.intervals[1]) %>% 
          select(gauge_code, which_moment, scale, rsquared),
        by = c("gauge_code", "which_moment"),
        all.x = FALSE)

# Figura momentos p/ uma estação
plot.scale.gauge <- 
  ggplot(data = df.mom.gauge, aes(x = d, y = mom)) +
  geom_point(color = "black", alpha = 0.9,
             aes(shape = factor(which_moment,
                                labels = paste0("q = ", unique(which_moment), ": H = ", round(unique(scale), 4),
                                                " | R² = ", round(unique(rsquared), 4))))) +
  geom_smooth(aes(group = which_moment),
              color = "grey20", linetype = "solid", alpha = 0.8,
              method = "lm", formula = y~x, se = FALSE, linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_comma()) +
  labs(y = "E[I<sub>d</sub><sup>q</sup>]", x = "Durações [h]", shape = "",
       caption = paste0("Estação: ", gauge, " | ", df.mom.gauge$n_year[1], " anos")) +
  theme_minimal() +
  theme(legend.position =  c(0.98, 0.98),
        legend.justification = c(1,0.7),
        legend.text = element_text(size = 9),
        axis.title.y = ggtext::element_markdown(),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1, color = "black"))); plot.scale.gauge

ggsave(plot = plot.scale.gauge, filename = "figuras/scale_invariance/workshop/plot_scale_gauge.png",
       width = 16, height = 12, units = "cm", dpi = 300)

# Figura histograma R²
plot.r2.hist <- scale.coefficient %>% 
  filter(gauge_code != "Regional",
         d == duration.intervals[1]) %>% 
  ggplot(aes(x = rsquared)) +
  facet_wrap(~which_moment, nrow = length(which.moment), labeller = as_labeller(function(x) paste0("q = ", x))) +
  geom_histogram(aes(y = after_stat(count/sum(count))), bins = 35, fill = "yellow2", color = "darkorange3") +
  geom_vline(xintercept = quantile(scale.coefficient$rsquared[scale.coefficient$gauge_code != "Regional"])[[3]],
             color = "red", linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "R²", y = "") +
  theme_minimal() +
  theme(legend.position =  "none",
        strip.placement = "inside",                      # fora dos rotulos do eixo-y
        strip.text = element_text(size = 10),             # 'strip.' altera configurações do facet_wrap
        axis.title.y = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12)); plot.r2.hist

ggsave(plot = plot.r2.hist, filename = "figuras/scale_invariance/workshop/plot_r2_hist.png",
       width = 8, height = 12, units = "cm", dpi = 300)


# INTERVALOS DE CONFIANÇA LOCAIS ------------------------------------------

groups <- unique(df.scale.boot$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

df.scale.boot <- df.scale.boot %>% 
  mutate(d = recode(d, !!!groups),
         d = factor(d, levels = group.names))

# Filtrar expoentes de escala 'scale'
df.scale <- df.scale.boot %>% 
  filter(which_variable == "scale",
         which_moment == 1,
         d != "Todas") %>% 
  select(gauge_code, n_year, d, statistic, ci_lower, ci_upper)

colors <- c("Todas" = "black", "Sub-horários" = "orange", "Horários" = "steelblue2", "Diários" =  "green3")

dodge.width <- 0.5

order <- df.scale %>% 
  filter(d == "Horários") %>% 
  arrange(statistic) %>% 
  pull(gauge_code) %>% 
  unique()

# Figura intervalos de confiança locais c/ bootstrap por blocos
plot.scale.ci.block.boot <- df.scale %>% 
  filter(gauge_code != "Regional") %>% 
  mutate(gauge_code = factor(gauge_code, levels = order)) %>% 
  ggplot(mapping = aes(y = statistic, ymin = ci_lower, ymax = ci_upper, color = d)) +
  geom_errorbar(aes(x = gauge_code), width = 0, linewidth = 0.8, alpha = 0.5, position = position_dodge(width = dodge.width)) +
  geom_point(aes(x = gauge_code), shape = 16, size = 1.5, alpha = 0.9, position = position_dodge(dodge.width)) +
  labs(x = "", y = "H", color = "Faixa de duração") +
  scale_color_manual(values = colors, aesthetics = c("color", "fill")) + 
  theme_minimal() +
  theme(legend.position = c(0.02,0.02),
        legend.justification = c(0,0),
        legend.key.spacing.y = unit(-2, "mm"),
        legend.margin = margin(l = 2, r = 2, unit = "mm"),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12))

ggsave(plot = plot.scale.ci.block.boot, filename = "figuras/scale_invariance/workshop/plot_scale_ci_boot.png",
       width = 16, height = 12, units = "cm", dpi = 300, bg = "white")

# Figura intervalos de confiança locais c/ bootstrap por blocos + faixas regional
plot.scale.boot.regional <- plot.scale.ci.block.boot +
  geom_rect(data = df.scale %>% filter(gauge_code == "Regional"), alpha = 0.25, color = "transparent", linewidth = 0,
            aes(xmin = -Inf, xmax = Inf, ymin = ci_lower, ymax = ci_upper, fill = d), show.legend = FALSE)

ggsave(plot = plot.scale.boot.regional, filename = "figuras/scale_invariance/workshop/plot_scale_ci_regional.png",
       width = 16, height = 12, units = "cm", dpi = 300, bg = "white")  

# MODELO REGIONAL HORÁRIO -------------------------------------------------

# Momentos e expoente de escala regionais
df.mom.regional <- 
  merge(x = scale.moments %>%
          mutate(gauge_code = "Regional") %>% 
          filter(d >= 1 & d <= 24),
        y = scale.coefficient %>%
          filter(gauge_code == "Regional", d == duration.intervals[3]) %>% 
          select(gauge_code, which_moment, scale, rsquared),
        by = c("gauge_code", "which_moment"),
        all.x = FALSE)

# Figura modelo regional horário
plot.scale.regional <- ggplot(data = df.mom.regional, aes(x = d, y = mom)) +
  geom_point(color = "grey", alpha = 0.5,
             aes(shape = factor(which_moment, labels = paste0("q = ", unique(which_moment),
                                                              ": H = ", round(unique(scale), 4),
                                                              " | R² = ", round(unique(rsquared), 3))))) +
  geom_smooth(aes(group = which_moment), color = "grey20", method = "lm", formula = y~x, se = FALSE, linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_comma()) +
  labs(y = "E[I<sub>d</sub><sup>q</sup>]", x = "Durações [h]", shape = "",
       caption = paste0("Modelo regional | Durações: ", duration.intervals[3], " h")) +
  theme_minimal() +
  theme(legend.position =  c(0.98, 0.98),
        legend.justification = c(1,0.7),
        legend.text = element_text(size = 9),
        axis.title.y = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
        plot.background = element_rect(color = "white"),
        # aspect.ratio = 1,
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1, color = "grey70")))

ggsave(plot = plot.scale.regional, filename = "figuras/scale_invariance/workshop/plot_scale_regional.png",
       width = 16, height = 12, units = "cm", dpi = 300)

# COEFICIENTES DE DESAGREGAÇÃO --------------------------------------------

source("scripts/script_coef_desagregacao.R")

ggsave(plot = plot.coef.disag, filename = "figuras/scale_invariance/workshop/plot_coef_desagregacao_block.png", bg = "white",
       width = 32, height = 16, units = "cm", dpi = 300)
