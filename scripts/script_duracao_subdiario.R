
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
               gghighlight
               )

font.size <- 12
font.family <- "serif"


# LER DADOS ---------------------------------------------------------------

# Dados
df_imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")
# df_imax <- df_imax[df_imax$d >= 0.5,]

# Vetor de estações que possuem ao menos N anos
# Manter somente estação, falhas e anos
df.quality <- df_imax %>% 
  select(gauge_code, year, na_prct) %>% 
  group_by(gauge_code, year) %>% 
  reframe(na_prct = first(na_prct)*100)

# Estabelecer limites p/ análise
threshold <- 10 # 10 %
min.years <- 5  # 5 anos  

# Filtrar dados e extrair vetor c/ codigo das estacoes
gauges <- df.quality %>% 
  filter(na_prct <= threshold) %>% # no máximo 'threshold'%
  count(gauge_code) %>%      # contar nro. anos, cria coluna 'n'
  filter(n >= min.years) %>% # pelo menos 'min.years' anos de dados
  pull(gauge_code) # extrair as estações


# RELAÇÕES ENTRE DURAÇÕES -------------------------------------------------

# Calcular 1º e 2º momentos centrais locais (p/ cada estação)
df.mom.local <- df_imax %>% 
  filter(gauge_code %in% gauges) %>% # filtrar por min.years e threshold
  group_by(gauge_code, d) %>%        # agrupar por duração e estação
  reframe(mom1 = mean(imax, na.rm = TRUE),   # valor esperado de Id:  E[Id^1]
          mom2 = mean(imax^2, na.rm = TRUE)) # valor esperado de Id²: E[Id^2]

# Calcular 1º e 2º momentos centrais regionais
df.mom.regional <- df_imax %>% 
  filter(gauge_code %in% gauges) %>% # filtrar por min.years e threshold
  group_by(d) %>%                    # agrupar por duração
  reframe(mom1 = mean(imax, na.rm = TRUE),   # valor esperado de Id:  E[Id^1]
          mom2 = mean(imax^2, na.rm = TRUE)) # valor esperado de Id²: E[Id^2]

df.mom.regional$gauge_code <- "Regional" # adicionar coluna

# Juntar local e regional
df.mom <- bind_rows(df.mom.local, df.mom.regional) %>% 
  pivot_longer(cols = starts_with("mom"),
               names_to = "moment",
               values_to = "value") %>% 
  mutate(moment = recode(moment, "mom1" = "q = 1", "mom2" = "q = 2"))

## REGRESSÃO --------------------------------------------------------------

# Regressão entre momentos e duração de cada estação
ls.mom <- split(x = df.mom.local, f = df.mom.local$gauge_code) # dividir em lista
df.rsquared <- lapply(X = ls.mom, FUN = function(gauge){
  
  name <- gauge$gauge_code[1] # código da estação
  mom1 <- log(gauge$mom1)     # vetor c/ logs do primeiro momento
  mom2 <- log(gauge$mom2)     # vetor c/ logs do primeiro momento
  durations <- log(gauge$d)   # vetor c/ logs das durações
  r.squared1 <- summary(lm(mom1~durations))$r.squared # regressão (pegar $sigma tambem)
  r.squared2 <- summary(lm(mom2~durations))$r.squared # regressão
  
  data.frame("gauge_code" = name,
             "rsquared1" = r.squared1,
             "rsquared2" = r.squared2)
  
}) # fim 'df.squared'

# Transformar em tbl_df
df.rsquared <- bind_rows(df.rsquared) %>% 
  pivot_longer(cols = starts_with("rsquared"), # unir colunas 'r_squared1' e 'r_squared2'
               names_to = "rsquared", # na coluna 'r_squared'             
               values_to = "value")    # valores na colunas 'value'

rownames(df.rsquared) <- NULL

# 


# VISUALIZAÇÃO ------------------------------------------------------------

# Visualização
plot.mom.scale <- ({
  
  # Curvas
  # Colocar H/r regionais no gráfico
  plot1 <- 
    ggplot(df.mom, aes(x = d, y = value, color = gauge_code)) +
    facet_wrap(~moment, nrow = 2, strip.position = "left") +      # colocar título do facet à esquerda
    geom_line(data = subset(df.mom, gauge_code != "Regional"), alpha = 0.7) +
    geom_smooth(data = subset(df.mom, gauge_code == "Regional"), method = lm, formula = y~x,
                aes(color = "Regional"), alpha = 0.6, linewidth = 0.5) +
    geom_point(data = subset(df.mom, gauge_code == "Regional"), aes(shape = "Regional")) +
    scale_x_log10(breaks = c(0.5, 1, 2, 6, 24, 72, 240)) +
    scale_y_log10() +
    scale_color_manual(values = c("Regional" = "grey10"), na.value = "lightblue") +
    scale_shape_manual(values = c("Regional" = 17)) +
    labs(y = "E[I<sub>d</sub><sup>q</sup>]", x = "Durações [h]")+
    theme_minimal() +
    theme(legend.position =  "none",
          axis.title.y = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
          strip.placement = "outside",                      # fora dos rotulos do eixo-y
          strip.text = element_text(size = 12),             # 'strip.' altera configurações do facet_wrap
          plot.background = element_rect(color = "white"),
          aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = font.family, color = "black", size = font.size)); plot1
  
  # Histogramas de R^2
  # Adicionar % de estaçoes acima de R^2 = 0.99
  plot2 <- 
    ggplot(df.rsquared, aes(x = value)) +
    facet_wrap(~rsquared, nrow = 2) +
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 20, fill = "yellow2", color = "darkorange3") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = c(0.98, 0.99, 1)) +
    labs(x = "R²", y = "") +
    theme_minimal() +
    theme(legend.position =  "none",
          axis.title.y = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
          strip.placement = "none",                         # remover título do 'facet' 
          strip.text = element_blank(),                     # 'strip.' altera configurações do facet_wrap
          plot.background = element_rect(color = "white"),
          aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = font.family, color = "black", size = font.size)); plot2
  
  # Combinar c/ 'patchwork'
  plot1 + plot2 #+ patchwork::plot_layout(widths = c(2,1))
  
}) # fim 'plot.mom.scale'
plot.mom.scale

ggsave(filename = "figuras/plot_mom_scale_invariance.png", plot = plot.mom.scale,
       width = 16, height = 14, units = "cm", dpi = 200)
