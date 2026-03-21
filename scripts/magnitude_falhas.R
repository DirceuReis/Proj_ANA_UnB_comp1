# Análise da magnitude dos máximos anuais removidos pelos filtros
# de percentual máximo de falhas e número mínimo de anos

rm(list = ls()); invisible(gc())

pacman::p_load(pacman, tidyverse, ggforce)

df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")
df.info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.pqt")

# Adicionar 'flag' de qualidade referente ao máximo de falhas anual
# 1: qualidade alta
# 0: qualidade baixa
na.accept <- 0.2
min.years <- 8
df.imax <- df.imax %>% mutate(flag = if_else(na_prct <= na.accept, 1, 0))

# Manter somente as observações do 'time_step' padrão da estação
data <- df.imax %>% 
  mutate(responsible = df.info$responsible[match(gauge_code, df.info$gauge_code)]) %>% 
  group_by(gauge_code) %>% 
  filter(d == min(d),
         n_distinct(year) >= min.years) %>% 
  ungroup()

n.gauges <- df.info %>% 
  filter(gauge_code %in% data$gauge_code) %>% 
  group_by(responsible) %>% 
  summarise(n_gauges = n())

# Mediana qualidade alta
median.hq <- data %>%
  group_by(gauge_code) %>% 
  filter(flag == 1) %>% 
  summarise(median = median(imax))

# Visualizar todas as estações
plot.magnitude.falhas <- 
  ggplot(data, aes(x = year, y = imax, color = factor(flag), shape = factor(flag))) +
  ggforce::facet_wrap_paginate(~gauge_code, scales = "free", ncol = 16, nrow = 14) +
  geom_hline(data = median.hq, aes(yintercept = median), linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 1) +
  scale_color_manual(values = c("1" = "steelblue", "0" = "red")) +
  scale_shape_manual(values = c("1" = 16, "0" = 8)) +
  labs(x = "") +
  theme_minimal() +
  theme(legend.position = "none",
        # legend.margin = margin(t = -10),
        # strip.switch.pad.wrap = margin(t = -10, b = -10),
        strip.text = element_text(size = 3.5, margin = margin(t = -0.1, b = 0.2)),
        axis.text.x = element_text(margin = margin(t = -0.5)),
        axis.text.y = element_text(margin = margin(l = -0.5)),
        text = element_text(color = "black", size = 4, family = "serif"),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank()); plot.magnitude.falhas

ggsave(plot = plot.magnitude.falhas, filename = "figuras/plot_magnitude_falhas.pdf",
       height = 297, width = 410, units = "mm", dpi = 300, device = "pdf")

# Calcular diferença relativa entre observações 'flag == 0' (qualidade baixa)
# e a mediana das observações 'flag == 1' (qualidade boa)
data.lq <- data %>% 
  filter(flag == 0) %>% 
  mutate(median_hq = median.hq$median[match(gauge_code, median.hq$gauge_code)],
         rel_dif = (imax - median_hq)/median_hq)

# Visualizar boxplots diferença relativa
plot.falhas.rel.dif <- ggplot(data.lq, aes(x = factor(responsible), y = rel_dif)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(outlier.shape = 8) +
  scale_y_continuous(limits = c(-1, 5), labels = scales::percent)  +
  labs(x = "", y = "Diferença relativa",
       caption = latex2exp::TeX(input = r"($DR = \frac{imax_{>20\%} - Md_{<20\%}}{Md_{<20\%}} \times 100$)", italic = TRUE)) +
  geom_text(aes(label = n_gauges, y = Inf), n.gauges, vjust = 1.6, size = 3, family = "serif") +
  theme_minimal() +
  theme(plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black")); plot.falhas.rel.dif

ggsave(plot = plot.falhas.rel.dif, filename = "figuras/plot_magnitude_falhas_boxplot.png",
       width = 16, height = 12, units = "cm", dpi = 300)
