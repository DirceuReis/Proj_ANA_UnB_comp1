
# O objetivo deste script é investigar variações sazonais nas intensidades
# máximas anuais para diferentes durações, avaliando o mês em que ocorrem
# e se pode haver diferença entre o ano hidrológico e o ano civil.

# Duas figuras serão geradas:
## 1. Boxplots do dia-ano de ocorrência dos máximos para cada estação;
## 2. Gráficos de rosa dos ventos para cada estação e um mapa indicando
##    a localização destas estações no Rio Grande do Sul.

# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow, patchwork, ggimage)

# LER DADOS ---------------------------------------------------------------

df.info <- arrow::read_parquet(file = "base/gerados/df_subdaily_info.pqt")
df.imax <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")


# ANÁLISE SAZONALIDADE ----------------------------------------------------

# Critérios de qualidade
na.accept <- 0.2
min.years <- 8

# Organização do conjunto de dados
data <- df.imax %>% 
  filter(na_prct <= na.accept) %>% 
  group_by(gauge_code, d) %>% 
  filter(n_distinct(year) >= min.years) %>% 
  ungroup() %>% 
  mutate(y_day = lubridate::yday(date),
         responsible = df.info$responsible[match(gauge_code, df.info$gauge_code)],
         lat = df.info$lat[match(gauge_code, df.info$gauge_code)],
         long = df.info$long[match(gauge_code, df.info$gauge_code)]) %>% 
  select(gauge_code, d, imax, date, y_day, year, responsible, lat, long)

start.month <- lubridate::yday(as.Date(paste0("1800/", 1:12, "/01")))
month.label <- data.frame(y_day = start.month,
                          day = as.Date(paste0("1800/", 1:12, "/01")),
                          label = lubridate::month(as.Date(paste0("1800/", 1:12, "/01")), label = FALSE))

# Ordenar pela mediana de imax p/ 1 h; ou
order <- data %>%
  filter(d == 1) %>%
  group_by(gauge_code) %>%
  summarise(y_day_median = median(y_day)) %>%
  arrange(desc(y_day_median)) %>%
  pull(gauge_code)

# Ordenar pela distância interquartílica de imax p/ 1h
order <- data %>%
  filter(d == 1) %>%
  group_by(gauge_code) %>%
  summarise(iqr = IQR(y_day)) %>%
  arrange(desc(iqr)) %>%
  pull(gauge_code)

# Durações
durations <- c(1, 6, 12, 24, 72)

# Boxplots
for(ds in durations){
  
  message("Duração: ", ds, " h")
  
  # order <- data %>% 
  #   filter(d == ds) %>% 
  #   group_by(gauge_code) %>% 
  #   summarise(y_day_median = median(y_day)) %>% 
  #   arrange(desc(y_day_median)) %>% 
  #   pull(gauge_code)
  
  data %>% 
    filter(d == ds) %>% 
    mutate(gauge_code = factor(gauge_code, levels = order)) %>% 
    ggplot(aes(x = gauge_code, y = y_day)) +
    stat_boxplot(geom = "errorbar", width = 0.4) +
    geom_boxplot(outlier.shape = 4, width = 0.5, fill = "grey90") +
    geom_hline(yintercept = start.month, linewidth = 0.4, linetype = "dashed", color = "grey10") +
    geom_text(aes(x =  Inf, y = y_day, label = label), month.label, family = "serif", size = 3, vjust = -0.9) +
    # geom_text(aes(x = gauge_code, y = Inf, label = n_years), n.years, family = "serif", size = 2.5, hjust = 0) +
    theme_minimal() +
    coord_flip(clip = "off") +
    labs(x = "", y = "Dia", caption = paste("Duração: ", ds, " h")) +
    theme(panel.border = element_rect(color = "black", linewidth = 0.3),
          plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"),
          panel.background = element_rect(color = "white"),
          axis.text.x = element_text(angle = 0, hjust = 1),
          text = element_text(family = "serif", color = "black", size = 10))
  
  ggsave(filename = paste0("figuras/sazonalidade/plot_sazonalidade_boxplot_", ds, "h.png"),
         width = 210, height = 297, units = "mm", dpi = 300)
  
}

# Todas as durações no mesmo gráfico - facet_wrap
data %>% 
  filter(d %in% durations) %>% 
  mutate(gauge_code = factor(gauge_code, levels = order)) %>% 
  ggplot(aes(x = gauge_code, y = y_day)) +
  facet_wrap(~factor(d), nrow = 1, strip.position = "bottom", labeller = as_labeller(function(x) paste("Duração:", x, "h"))) +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(outlier.shape = 4, width = 0.5, fill = "grey90") +
  geom_hline(yintercept = start.month, linewidth = 0.4, linetype = "dashed", color = "grey10") +
  geom_text(aes(x =  Inf, y = y_day, label = label), month.label, family = "serif", size = 2.5, vjust = -0.9) +
  theme_minimal() +
  coord_flip(clip = "off") +
  labs(x = "", y = "Dia") +
  theme(panel.border = element_rect(color = "black", linewidth = 0.3),
        plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"),
        panel.background = element_rect(color = "white"),
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(family = "serif", color = "black", size = 10))

ggsave(filename = paste0("figuras/sazonalidade/plot_sazonalidade_boxplot_facet_wrap.png"),
       width = 297, height = 210, units = "mm", dpi = 300)

# Todas as durações no mesmo gráfico - fill
data %>% 
  filter(d %in% durations) %>% 
  mutate(gauge_code = factor(gauge_code, levels = order),
         d = factor(d, labels = paste(unique(d), "h"))
         ) %>% 
  ggplot(aes(x = factor(d), y = y_day)) +
  facet_wrap(~gauge_code) +
  coord_cartesian(clip = "off") +
  stat_boxplot(geom = "errorbar", linewidth = 0.2, width = 0.2) +
  geom_boxplot(aes(group = d), outlier.shape = 4, width = 0.5, linewidth = 0.2, fill = "grey90") +
  theme_minimal() +
  labs(x = "", y = "Dia") +
  theme(panel.border = element_rect(color = "black", linewidth = 0.3),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "lines"),
        panel.background = element_rect(color = "white"),
        legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 1),
        text = element_text(family = "serif", color = "black", size = 10))


ggsave(filename = paste0("figuras/sazonalidade/plot_sazonalidade_boxplot_fill.png"),
       width = 297, height = 210, units = "mm", dpi = 300)


# ROSA DOS VENTOS ---------------------------------------------------------

# Função p/ transformar datas em radianos
date.to.angle <- function(date){
  
  year.start <- lubridate::floor_date(date, "year")
  year.end <- lubridate::ceiling_date(date, "year")
  year.length <- as.numeric(difftime(year.end, year.start, units = "sec"))
  time.passed <- as.numeric(difftime(date, year.start, units = "sec"))
  
  date.angle <- 2*pi*time.passed/year.length
  return(date.angle)
  
}

start.month.rad <- 2*pi*start.month/365.25

# Selecionar 4 estações c/ maior IQR
iqr.large <- data %>% 
  filter(d == 1) %>%
  group_by(gauge_code) %>%
  summarise(iqr = IQR(y_day)) %>%
  arrange(iqr) %>%
  slice_tail(n = 4)

# Selecionar 4 estações c/ menor IQR
iqr.small <- data %>% 
  filter(d == 1) %>%
  group_by(gauge_code) %>%
  summarise(iqr = IQR(y_day)) %>%
  arrange(iqr) %>%
  slice_head(n = 4)

# Visualizar rosa dos ventos da estação com maior distância interquartílica
data %>% 
  filter(d %in% durations,
         gauge_code == iqr.large$gauge_code[1]) %>%
  mutate(date_rad = date.to.angle(date),
         d = factor(d, labels = paste(unique(d), "h"))) %>% 
  ggplot(aes(x = date_rad, y = imax)) +
  facet_wrap(~d, ncol = 2, nrow = 2, scales = "free") +
  coord_polar("x", start = 0, direction = 1) +
  geom_segment(aes(y = 0, xend = date_rad, yend = imax), color = "grey20") +
  # geom_col() +
  scale_x_continuous(name = "Ocorrência", breaks = start.month.rad,
                     labels = month(1:12, label = TRUE), limits = c(0, 2*pi)) +
  scale_y_continuous(name = "Intensidades máximas anuais [mm]") +
  labs(caption = paste("Estação", iqr.large$gauge_code[1])) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", linewidth = 0.3),
        # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        panel.background = element_rect(color = "white"),
        panel.spacing.x = unit(1.8, "lines"),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text.y = element_blank(), 
        axis.text.x = element_text(vjust = -0.9, color = "black", size = 8, margin = margin(b = 1)),
        text = element_text(family = "serif", color = "black", size = 10))

ggsave(filename = "figuras/sazonalidade/plot_polar_maior_iqr.png",
       width = 16, height = 16, units = "cm", dpi = 300)

# Plotar os 4 IQRs mais extremos da duração de 1 h
pacman::p_load(ggh4x)
plot.small.iqr <- data %>% 
  filter(d %in% durations,
         gauge_code %in% iqr.small$gauge_code) %>%
  mutate(date_rad = date.to.angle(date),
         gauge_code = factor(gauge_code, c(iqr.small$gauge_code, iqr.large$gauge_code))) %>% 
  ggplot(aes(x = date_rad, y = imax)) +
  ggh4x::facet_grid2(d ~ gauge_code, scales = "free", independent = "all", labeller = labeller(d = ~paste(.x, "h"))) +
  coord_polar("x", start = 0, direction = 1) +
  geom_segment(aes(y = 0, xend = date_rad, yend = imax), color = "grey20") +
  # geom_col(color = "grey20",) +
  scale_x_continuous(name = "Ocorrência", breaks = start.month.rad, labels = month(1:12, label = FALSE), limits = c(0, 2*pi)) +
  scale_y_continuous(name = "Intensidades máximas anuais [mm]", expand = c(0,0)) +
  labs(title = paste(min(iqr.small$iqr), "< DIQ <", max(iqr.small$iqr))) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", linewidth = 0.3),
        # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        panel.background = element_rect(color = "white"),
        # panel.spacing.x = unit(1.8, "lines"),
        legend.position = "none",
        # aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(vjust = -0.9, color = "black", size = 8, margin = margin(b = 1)),
        text = element_text(family = "serif", color = "black", size = 10)); plot.small.iqr

plot.large.iqr <- data %>% 
  filter(d %in% durations,
         gauge_code %in% iqr.large$gauge_code) %>%
  mutate(date_rad = date.to.angle(date),
         gauge_code = factor(gauge_code, c(iqr.small$gauge_code, iqr.large$gauge_code))) %>% 
  ggplot(aes(x = date_rad, y = imax)) +
  ggh4x::facet_grid2(d ~ gauge_code, scales = "free", independent = "all", labeller = labeller(d = ~paste(.x, "h"))) +
  coord_polar("x", start = 0, direction = 1) +
  geom_segment(aes(y = 0, xend = date_rad, yend = imax), color = "grey20") +
  # geom_col(color = "grey20") +
  scale_x_continuous(name = "Ocorrência", breaks = start.month.rad, labels = month(1:12, label = FALSE), limits = c(0, 2*pi)) +
  scale_y_continuous(name = "Intensidades máximas anuais [mm]", expand = c(0,0)) +
  labs(title = paste(min(iqr.large$iqr), "< DIQ <", max(iqr.large$iqr))) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", linewidth = 0.3),
        # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        panel.background = element_rect(color = "white"),
        # panel.spacing.x = unit(1.8, "lines"),
        legend.position = "none",
        # aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(vjust = -0.9, color = "black", size = 8),
        text = element_text(family = "serif", color = "black", size = 10)); plot.large.iqr

plot.small.iqr + plot.large.iqr + patchwork::plot_layout(ncol = 2, axis_titles = "collect")

ggsave(filename = "figuras/sazonalidade/plot_polar_iqr.png",
       width = 297, height = 170, units = "mm", dpi = 300)


# INVERNO-VERÃO -----------------------------------------------------------

# Visualizar roda dos ventos das maiores estações
n.years <- data %>% 
  group_by(gauge_code) %>% 
  summarise(n_years = n_distinct(year)) %>% 
  arrange(desc(n_years)) %>% 
  head(4)

# Definir meses de início dos periodos
umido1 <- as.Date(c("1800-01-01", "1800-04-01"))
seco <- as.Date(c("1800-04-01", "1800-09-01"))
umido2 <- as.Date(c("1800-09-01", "1800-12-31"))

seasons <- rbind(date.to.angle(umido1),
                 date.to.angle(seco),
                 date.to.angle(umido2) + c(0, + 2*pi/365.25)) %>% 
  as.data.frame() %>%
  reframe(start = V1, end = V2) %>% 
  mutate(season = c("Verão", "Inverno", "Verão"))

plot.polar <- {(
  data %>% 
    filter(d %in% durations, gauge_code %in% n.years$gauge_code) %>%
    mutate(date_rad = date.to.angle(date), gauge_code = factor(gauge_code, labels = paste0(unique(gauge_code), " (", n.years$n_years, " anos)"))) %>%
    ggplot() +
    ggh4x::facet_grid2(d ~ gauge_code, scales = "free", independent = "all", labeller = labeller(d = ~paste(.x, "h"))) +
    geom_rect(aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = season), seasons, alpha = 0.4) +
    geom_segment(aes(x = date_rad, y = 0, xend = date_rad, yend = imax), color = "black") +
    scale_x_continuous(breaks = start.month.rad, labels = month(1:12, label = FALSE), limits = c(0, 2*pi)) +
    scale_fill_manual(values = c("Inverno" = "darkgoldenrod1", "Verão" = "steelblue")) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Mês de ocorrência", y = "Intensidades máximas anuais", fill = "") +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", linewidth = 0.3),
          panel.background = element_rect(color = "white"),
          legend.position = "bottom",
          aspect.ratio = 1,
          panel.grid.minor = element_blank(),
          plot.background = element_rect(color = "white"),
          plot.title = element_text(size = 8, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          # axis.text.x = element_text(vjust = -0.9, color = "black", size = 8),
          axis.text.x = element_text(size = 8, margin = margin(b = 1)),
          text = element_text(family = "serif", color = "black", size = 10)))}; plot.polar

plot.polar <- plot.polar + coord_polar("x", start = 0, direction = 1); plot.polar

ggsave(plot = plot.polar, filename = "figuras/sazonalidade/plot_polar_umido_seco.png",
       height = 21, width = 16, units = "cm", dpi = 300)


# PLOT REGIONAL ROSA DOS VENTOS -------------------------------------------

# Analisar a "direção" das rosa dos ventos plotando no local de cada estação
# e avaliando a variação espacial

# Avaliar quais são as estações que apresentam boa cobertura do estado
n <- 16 # nro. estações



# Pacote `ggimage()`

# INCIDÊNCIA MÁXIMOS INVERNO-VERÃO ----------------------------------------

# Analisar a frequência com que os máximos anuais identificados no ano civil
# ocorrem ou no inverno ou no verão do hemisfério Sul

