
# PACOTES -----------------------------------------------------------------

rm(list = ls()); gc()

if(!require(pacman)) install.packates("pacman")
pacman::p_load(pacman,
               arrow,
               tidyverse,
               lubridate,
               patchwork,
               ggplot2,
               ggtext,
               gghighlight,
               ggforce
               )


# FUNÇÕES -----------------------------------------------------------------

# Regressão linear e análise invariância de escala
source("scripts/funcoes/fun_scale_invariance.R")

# LER DADOS ---------------------------------------------------------------

# Dados
df_imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")


# APLICAR FUNÇÃO 'fun_scale_invariance' -----------------------------------

# Argumentos
threshold <- 0.1            # percentual de falhas (0.1 -> 10%)
min.years <- 5              # mínimo de anos na série
which.durations <- c(1, 24) # intervalo de durações
which.moment <- 1           # qual momento usar de base p/ gerar figuras individuais

# Avaliar diferentes modelos
duration.intervals <- list(c(10/60, 10*24), # todas as durações
                           c(10, 60)/60,    # sub-horários
                           c(1, 24),        # horários
                           c(1, 10)*24)     # diários

# Regressão + coeficiente de escala + plots
scale.invariance <- lapply(X = duration.intervals, FUN = function(interval){
  
  res <- fun_scale_invariance(df.imax = df_imax,
                              threshold = threshold,
                              min.years = min.years,
                              which.moment = which.moment,
                              which.duration = interval,
                              font.family = "serif",
                              font.size = 12,
                              plot.dim = c(6,6))
  
})

names(scale.invariance) <- c("all", "subhourly", "hourly", "daily")

# Dimensões gráfico (pixels)
height <- 2160
width <- 3840
units <- "px"

# Salvar gráficos
for(i in seq_along(scale.invariance)){
  
  name <- names(scale.invariance[i])
  interval <- scale.invariance[[name]]
  
  dir.name <- paste0("figuras/scale_invariance/", name, "_duration/")
  if(isFALSE(dir.exists(dir.name))) dir.create(dir.name)
  
  plot.all <- interval$plot.scale.all
  ggsave(plot = plot.all, filename = paste0(dir.name, "plot_scale_invariance_", name, ".png"),
         width = 16, height = 14, units = "cm", dpi = 300)
  
  for(group in seq_along(interval$plot.scale.each)){

    plot.each <- interval$plot.scale.each[[group]]
    ggsave(plot = plot.each, filename = paste0(dir.name, "plot_scale_invariance_", name, "_group", group, ".png"),
           width = width, height = height, units = units, dpi = 300)

  }
  
}

# Salvar 'scale.invariance'
df.regression <- scale.invariance$hourly$regression
arrow::write_parquet(x = df.regression, sink = "base/gerados/df_scale_invariance.parquet")

# HISTOGRAMA DE SCALE -----------------------------------------------------

# Análise horária

# Diferença relativa entre scale locais e scale regional
# Refazer o cálcuulo da diferença relativa
df.scale <- df.regression %>% 
  filter(which_moment == "1") %>% 
  group_by(gauge_code) %>% 
  summarise(scale = mean(scale)) %>% 
  rename("scale_local" = "scale") %>% 
  mutate(scale_regional = scale.regional,
         relative_diff = (scale_local - scale_regional)/scale_regional*100)

# Estabelecer um limite p/ seleção de estações
limit <- c(-10, 10) #%

plot.scale.hist <- ({

  ggplot(df.scale, aes(x = relative_diff)) +
    # fact_wrap(~which_moment, ncol = 2) +
    geom_ribbon(aes(xmin = limit[1], xmax = limit[2])) +
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 20, fill = "yellow2", color = "darkorange3") +
    geom_vline(xintercept = first(df.scale$scale_regional),
               color = "red", linetype = "dashed", linewidth = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "(H<sub>i</sub> - H<sub>R</sub>)/H<sub>r</sub> × 100", y = "") +
    theme_minimal() +
    theme(legend.position =  "none",
          axis.title.x = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
          strip.placement = "none",                         # remover título do 'facet' 
          strip.text = element_blank(),                     # 'strip.' altera configurações do facet_wrap
          plot.background = element_rect(color = "white"),
          aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = font.family, color = "black", size = font.size))

})
