
# PACOTES -----------------------------------------------------------------

rm(list = ls()); invisible(gc())

if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,
               pbapply,
               arrow,
               tidyverse,
               tidymodels,
               boot,
               patchwork,
               parallel,
               ggtext,
               gghighlight,
               ggforce
               )


# FUNÇÕES -----------------------------------------------------------------

# Regressão linear e análise invariância de escala
# source("scripts/funcoes/fun_scale_invariance.R")
source("scripts/funcoes/fun_scale_invariance_more_moments.R")
source("scripts/funcoes/fun_scale_invariance_boot.R")

# LER DADOS ---------------------------------------------------------------

# Dados subdiários
df_imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")
df_imax <- df_imax[df_imax$gauge_code != "SBPA",] # estação dando problema


# APLICAR FUNÇÃO 'fun_scale_invariance' -----------------------------------

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
scale.invariance <- pbapply::pblapply(X = duration.intervals, FUN = function(interval){
  
  res <- fun_scale_invariance(df.imax = df_imax,
                              na.accept = na.accept,
                              min.years = min.years,
                              which.moment = which.moment,
                              which.duration = interval,
                              min.duration = 3,
                              font.family = "serif",
                              font.size = 12,
                              plot.dim = c(5,4))

}); names(scale.invariance) <- c("all", "subhourly", "hourly", "daily")

# Dimensões gráfico (pixels)
height <- 2160
width <- 3840
units <- "px"

# Salvar gráficos
for(i in seq_along(scale.invariance)){
  
  name <- names(scale.invariance[i])
  interval <- scale.invariance[[name]]
  
  dir.name <- paste0("figuras/scale_invariance/", min.years, " anos ", na.accept*100, " prct/", name, "_duration/")
  if(isFALSE(dir.exists(dir.name))) dir.create(dir.name)
  
  plot.all <- interval$plot.scale.all
  ggsave(plot = plot.all, filename = paste0(dir.name, "plot_scale_invariance_", name, ".png"),
         width = 16, height = 12, units = "cm", dpi = 300)
  
  for(group in seq_along(interval$plot.scale.each)){

    plot.each <- interval$plot.scale.each[[group]]
    ggsave(plot = plot.each, filename = paste0(dir.name, "plot_scale_invariance_", name, "_group", group, ".png"),
           width = width, height = height, units = units, dpi = 300)
    
  }
  
  message("Gerando figuras em: ", dir.name, "...")
  
}

# Salvar resultados de 'scale.invariance'
scale.moments <- lapply(X =  scale.invariance, FUN = function(interval){
  interval$scale.moments
}); scale.moments <- bind_rows(scale.moments)

scale.coefficient <- lapply(X = scale.invariance, FUN = function(interval){
  interval$scale.coefficient
});  scale.coefficient <- bind_rows(scale.coefficient)


# DIFERENÇA ESTATÍSTICA ---------------------------------------------------

# Intervalos normais aproximados
scale.coefficient$ci_lower <- scale.coefficient$scale - 1.96*scale.coefficient$se
scale.coefficient$ci_upper <- scale.coefficient$scale + 1.96*scale.coefficient$se

arrow::write_parquet(x = scale.coefficient, sink = "base/gerados/df_scale_coefficient.parquet")
df_scale <- arrow::read_parquet(file = "base/gerados/df_scale_coefficient.parquet")


# BOOTSTRAP ---------------------------------------------------------------

# Paralelização
n.cores <- detectCores() - 2
cluster <- parallel::makeCluster(n.cores)

# Calcular intervalos de confiança usando bootstrap
R.boot <- 1e4
# df.scale.boot <- fun_scale_invariance_boot(df.imax = df_imax,
#                                            na.accept = na.accept,
#                                            min.years = min.years,
#                                            which.moment = 1:3,
#                                            which.duration = c(1,24),
#                                            min.duration = 3,
#                                            R.boot = R.boot,
#                                            cl = cluster); parallel::stopCluster(cluster)

scale.boot.ci <- lapply(X = duration.intervals, FUN = function(interval){
  
  message("Estimando intervalo...")
  df.scale.boot <- fun_scale_invariance_boot(df.imax = df_imax,
                                             na.accept = na.accept,
                                             min.years = min.years,
                                             which.moment = 1:3,
                                             which.duration = c(1,24),
                                             min.duration = 3,
                                             R.boot = R.boot,
                                             cl = cluster)
  
}); parallel::stopCluster(cluster); names(scale.invariance) <- c("all", "subhourly", "hourly", "daily")

# arrow::write_parquet(df.scale.boot, sink = "base/gerados/df_scale_boot_hour.parquet")
# df_scale_boot <- arrow::read_parquet(file = "base/gerados/df_scale_boot_hour.parquet")


# VISUALIZAÇÃO INTERVALOS DE CONFIANÇA ------------------------------------

groups <- unique(df_scale$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

colors <- c("Todas" = "black", "Sub-horários" = "orange", "Horários" = "steelblue2", "Diários" =  "green")

# Intervalos normais aproximados
# Todas as durações no mesmo gráfico
df_scale %>% 
  mutate(d = recode(d, !!!groups),
         d = factor(d, levels = group.names),
         gauge_code = forcats::fct_reorder(gauge_code, n_years)) %>% 
  group_by(gauge_code, d) %>% 
  summarise_if(is.numeric, ~mean(.x, na.rm = TRUE)) %>% 
  ggplot(aes(x = gauge_code, y = scale, ymin = ci_lower, ymax = ci_upper, color = d)) +
  geom_errorbar(width = 0, linewidth = 0.8, alpha = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(pch = 16, size = 1.5, alpha = 0.9, position = position_dodge(width = 0.5)) +
  # geom_text(data = . %>% filter(gauge_code != "Regional"), aes(y = Inf, label = n_years),
  #            color = "black", vjust = 2, size = 2.5, family = "mono") +
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(0.3, 0.95, 0.1)) +
  labs(x = "", y = "H", color = "") +
  theme_minimal() +
  theme(legend.position = c(0.98,0.02),
        legend.justification = c(1,0),
        legend.key.spacing.y = unit(-2, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 10))

ggsave(filename = "figuras/scale_invariance/scale_confidence_interval.png",
       width = 16, height = 12, units = "cm", bg = "white")

# Com resultados do bootstrap
df_scale_boot %>% 
  group_by(gauge_code) %>% 
  summarise_all(~mean(.x, na.rm = TRUE)) %>% 
  ggplot(aes(x = gauge_code, y = scale, ymin = ci_lower, ymax = ci_upper)) +
  geom_errorbar(width = 0, linewidth = 0.8, alpha = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(pch = 16, size = 1.5, alpha = 0.9, position = position_dodge(width = 0.5)) +
  scale_y_continuous(breaks = seq(0.3, 0.95, 0.1)) +
  labs(x = "", y = "H") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 10))

# Comparar dois ICs gerados
bind_rows(
  df_scale %>% 
    mutate(d = recode(d, !!!groups),
           d = factor(d, levels = group.names)) %>% 
    filter(gauge_code != "Regional",
           d == "Horários") %>% 
    group_by(gauge_code) %>% 
    summarise_if(is.numeric, ~mean(.x, na.rm = TRUE)) %>% 
    select(gauge_code, scale, ci_lower, ci_upper) %>% 
    mutate(which_ci = "Aproximado"),
  df_scale_boot %>% 
    select(gauge_code, scale, ci_lower, ci_upper) %>% 
    group_by(gauge_code) %>% 
    summarise_all(~mean(.x, na.rm = TRUE)) %>% 
    mutate(which_ci = "Bootstrap")) %>% 
  ggplot(aes(x = gauge_code, y = scale, color = which_ci, ymin = ci_lower, ymax = ci_upper)) +
  geom_errorbar(width = 0, linewidth = 0.8, alpha = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(pch = 16, size = 1.5, alpha = 0.9, position = position_dodge(width = 0.5)) +
  scale_y_continuous(breaks = seq(0.3, 0.95, 0.1)) +
  scale_color_manual(values = c("Aproximado" = "red", "Bootstrap" = "blue")) +
  labs(x = "", y = "H") +
  theme_minimal() +
  theme(legend.position = c(0.98,0.02),
        legend.justification = c(1,0),
        legend.key.spacing.y = unit(-1, "mm"),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 10))

ggsave(filename = "figuras/scale_invariance/scale_confidence_intervals_approx_bootstrap_hourly.png",
       width = 16, height = 12, units = "cm", bg = "white")


# COEFICIENTES DE DESAGREGAÇÃO --------------------------------------------

# Calcular média entre os momentos calculados
groups <- unique(df_scale$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

df.scale.mean <- df_scale %>% 
  group_by(gauge_code, d) %>% 
  summarise_if(is.numeric, ~mean(.x, na.rm = TRUE)) %>% # calcular média entre os momentos
  mutate(relative_width = (ci_upper - ci_lower)/scale,  # largura relativa dos intervalos de confiança
         d = recode(d, !!!groups),                      # alterar nome das durações
         d = factor(d, levels = group.names))           # estabelecer ordem dos grupos de duração

# Calcular teste-t pareado p/ avaliar diferença entre as séries do coeficiente de escala
test.pairwise.ttest <- pairwise.t.test(x = df.scale.mean$scale,
                                       g = df.scale.mean$d,
                                       p.adjust.method = "bonferroni")

print(test.pairwise.ttest)

# ANOVA
test.anova <- aov(formula = scale~d, data = df.scale.mean)
summary(test.anova)

# Tukey
TukeyHSD(test.anova)

# Calcular número de estações utilizadas para calcular cada modelo de H
# e o número efetivo de anos utilizados no cálculo
gauges.by.model <- df.scale.mean %>% 
  filter(gauge_code != "Regional") %>% 
  group_by(d) %>% 
  summarise(n_gauges = n(),
            n_years = sum(n_years))


# Boxplot Hggplot(data = df.scale.mean, aes(x = d, y = se)) +
ggplot(data = df.scale.mean, aes(x = d, y = scale)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(aes(fill = d), linewidth = 0.5, outlier.shape = 4, outlier.size = 2, show.legend = FALSE) +
  geom_jitter(alpha = 0.6, size = 0.8) +
  geom_text(data = gauges.by.model, mapping = aes(x = d, y = Inf, label = sprintf("n. estações: %d", n_gauges)),
            family = "mono", fontface = "bold", size = 3, vjust = 2) +
  # scale_fill_manual(values = c("Todas" = "lightblue4"), na.value = "lightblue2") +
  scale_fill_manual(values = c("Todas" = "steelblue"), na.value = "yellow2") +
  scale_y_continuous(limits = c(0.4, 0.85)) +
  labs(x = "", y = "H") +
  theme_minimal() +
  theme(legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12))

# Boxplot da largura relativa do intervalo de confiança aproximado
ggplot(data = df.scale.mean, aes(x = d, y = relative_width)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(aes(fill = d), linewidth = 0.5, outlier.shape = 4, outlier.size = 2, show.legend = FALSE) +
  # scale_fill_manual(values = c("Todas" = "lightblue4"), na.value = "lightblue2") +
  geom_text(data = gauges.by.model, mapping = aes(x = d, y = Inf, label = sprintf("n. estações: %d", n_gauges)),
            family = "mono", fontface = "bold", size = 3, vjust = 2) +
  scale_fill_manual(values = c("Todas" = "steelblue"), na.value = "yellow2") +
  labs(x = "", y = "Largura relativa dos ICs") +
  theme_minimal() +
  theme(legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12))


# COMPARAÇÃO COEFICIENTES DE DESAGREGAÇÃO ---------------------------------

# Fatores de escala regionais
regional <- df.scale.mean %>% 
  filter(gauge_code == "Regional")

# Tabela CETESB coeficiente de desagregação
df.cetesb <- data.frame(relacao = c("5 min/30 min", "10 min/30 min", "15 min/30 min", "20 min/30 min", "25 min/30 min", 
                                    "30 min/1 h", "1 h/24 h", "6 h/24 h", "8 h/24 h", "10 h/24 h", "12 h/24 h"),
                        d1 = c(seq(5, 30, 5)/60, 1, 6, 8, 10, 12),
                        d2 = c(rep(0.5, 5), 1, rep(24, 5)),
                        coef_pdmax = c(0.34, 0.54, 0.7, 0.81, 0.91, 0.74, 0.42, 0.72, 0.78, 0.82, 0.85))

# Passar todos p/ base de 24 h
df.cetesb$relacao_mod <- as.character(NA) # nova coluna p/ nome das relações
df.cetesb$coef_mod <- as.numeric(NA)      # nova coluna p/ valor coeficientes de desagregação

for(i in 1:nrow(df.cetesb)){
  
  if(df.cetesb$d2[i] <= 1){
    
    df.cetesb$relacao_mod[i] <- paste0(df.cetesb$d1[i]*60, " min/24 h")
    df.cetesb$coef_mod[i] <- df.cetesb$coef_pdmax[i] * df.cetesb$coef_pdmax[6] * df.cetesb$coef_pdmax[7]
    df.cetesb$coef_regional_multi[i] <- (df.cetesb$d1[i]/24)^{-regional[regional$d == group.names[2],][["scale"]] + 1} # H sub-horários
    
    if(df.cetesb$d2[i] == 1){
      df.cetesb$coef_mod[i] <- df.cetesb$coef_pdmax[i] * df.cetesb$coef_pdmax[7]
    }
    
  } else{
    
    df.cetesb$relacao_mod[i] <- paste0(df.cetesb$d1[i], " h/24 h")
    df.cetesb$coef_mod[i] <- df.cetesb$coef_pdmax[i]
    df.cetesb$coef_regional_multi[i] <- (df.cetesb$d1[i]/24)^{-regional[regional$d == group.names[3],][["scale"]] + 1} # H horários
    
  }
  
  # Calcular coeficiente de desagregação c/ H regional
  df.cetesb$coef_regional_all[i] <- (df.cetesb$d1[i]/24)^{-regional[regional$d == group.names[1],][["scale"]] + 1}
  
}

# Gerar tabela em LaTeX
print(xtable::xtable(df.cetesb[-c(2,3)]), include.rownames = FALSE, label = "tab.coef-des-cetesb1")
