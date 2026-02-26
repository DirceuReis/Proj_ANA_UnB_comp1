
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow)


# LER DADOS ---------------------------------------------------------------

# Expoentes de escala e coeficientes de desagregação locais e regionais
# com bootstrap por bloco
df.scale.boot <- arrow::read_parquet(file = "base/gerados/scale_invariance/df_scale_block_boot.pqt")

groups <- unique(df.scale.boot$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

df.scale.boot <- df.scale.boot %>% 
  mutate(d = recode(d, !!!groups),
         d = factor(d, levels = group.names))

# Coeficientes de desagregação da CETESB
df.cetesb <- data.frame(relacao = c("5 min/30 min", "10 min/30 min", "15 min/30 min", "20 min/30 min", "25 min/30 min", 
                                    "30 min/1 h", "1 h/24 h", "6 h/24 h", "8 h/24 h", "10 h/24 h", "12 h/24 h"),
                        d1 = c(seq(5, 30, 5)/60, 1, 6, 8, 10, 12),
                        d2 = c(rep(0.5, 5), 1, rep(24, 5)),
                        d_target = factor(c(seq(5, 30, 5)/60, 1, 6, 8, 10, 12)),
                        coef_pdmax = c(0.34, 0.54, 0.7, 0.81, 0.91, 0.74, 0.42, 0.72, 0.78, 0.82, 0.85))

df.cetesb$relacao_mod <- as.character(NA)
df.cetesb$coef_cetesb <- as.numeric(NA)

for(i in 1:nrow(df.cetesb)){
  
  if(df.cetesb$d2[i] <= 1){
    
    df.cetesb$relacao_mod[i] <- paste0(df.cetesb$d1[i]*60, " min/1 h")
    df.cetesb$coef_cetesb[i] <- df.cetesb$coef_pdmax[i] * df.cetesb$coef_pdmax[6]
    
  } else{
    
    df.cetesb$relacao_mod[i] <- paste0(df.cetesb$d1[i], " h/24 h")
    df.cetesb$coef_cetesb[i] <- df.cetesb$coef_pdmax[i]
    
  }
}


# VISUALIZAÇÃO ------------------------------------------------------------

# Cores
colors <- c("Todas" = "black",
            "Sub-horários" = "orange",
            "Horários" = "steelblue2",
            "Diários" =  "green3",
            "CETESB" = "midnightblue",
            "IC<sub>R</sub> 95%" = "grey")

filter_coef <- c(1, 6, 12)

# Ordem na figura
order <- df.scale.boot %>%
  filter(d == "Horários", which_variable == "scale") %>% 
  arrange(statistic) %>% 
  pull(gauge_code) %>%
  unique()

# Dados locais
df.local <- df.scale.boot %>% 
  filter(gauge_code != "Regional",
         which_variable == "coef",
         d == "Horários",
         which_moment == 1,
         d_target %in% filter_coef) %>% 
  mutate(gauge_code = factor(gauge_code, levels = order))

# Dados regionais
df.regional <- df.scale.boot %>% 
  filter(gauge_code == "Regional",
         which_variable == "coef",
         d == "Horários",
         which_moment == 1,
         d_target %in% filter_coef) %>% 
  select(-gauge_code) # evitar conflitos c/ ggplot

plot.coef.disag <- ggplot(mapping = aes(color = d)) +
  facet_wrap(~d_target, ncol = 3, labeller = as_labeller(function(x) paste0(x, "/24 h"))) +
  geom_rect(data = df.regional, alpha = 0.4, color = "transparent", linewidth = 0,
            aes(xmin = 0, xmax = Inf, ymin = ci_lower, ymax = ci_upper, fill = "IC<sub>R</sub> 95%")) +
  geom_hline(data = df.cetesb %>% filter(d1 %in% filter_coef), linetype = "dashed",
             aes(yintercept = coef_cetesb, color = "CETESB")) +
  geom_errorbar(data = df.local, width = 0,
                aes(x = gauge_code, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(data = df.local, size = 1,
             aes(x = gauge_code, y = statistic)) +
  scale_color_manual(values = colors, aesthetics = c("color", "fill")) +
  labs(x = "", y = "C<sub>P</sub>", color = "", fill = "") +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(legend.position = "bottom",
        legend.text = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA))

# ggsave(filename = "figuras/scale_invariance/coef_desagregacao_block.png", bg = "white",
#        width = 16, height = 20, units = "cm", dpi = 300)
  
# Interessante que quanto mais distante da duração base 'D' (ex. 24 h), maiores ficam
# os intervalos de confiança, oq pode ser reflexo da propagação de erro e que
# próximo do final das durações alvo 'd' (ex. 1 h) o ajuste do modelo linear já
# não seja tão bom


# COMPARAÇÃO COEFICIENTES -------------------------------------------------

# Comparação coeficientes locais vs. regionais
df.local %>% 
  select(gauge_code, d, d_target, statistic) %>% 
  mutate(local = statistic,
         regional = rep(df.regional$statistic, n_distinct(gauge_code)),
         cetesb = rep(df.cetesb[df.cetesb$d_target %in% filter_coef,][["coef_cetesb"]], n_distinct(gauge_code)),
         rel_dif__regional_local = (regional - local)/regional*100,
         rel_dif__cetesb_regional = (cetesb - regional)/cetesb*100,
         rel_dif__cetesb_local = (cetesb - local)/cetesb*100) %>% 
  pivot_longer(cols = starts_with("rel_dif__"),
               names_to = "rel_dif",
               names_prefix = "rel_dif__",
               values_to = "values") %>% 
  mutate(rel_dif = factor(rel_dif, c("regional_local", "cetesb_local", "cetesb_regional"))) %>% 
  ggplot(aes(x = rel_dif, y = values, color = rel_dif)) +
  geom_jitter() +
  stat_boxplot() +
  theme(legend.position = "none")

# Comparação coeficientes locais vs. CETESB

# Comparação coeficientes regionais vs. CETESB
