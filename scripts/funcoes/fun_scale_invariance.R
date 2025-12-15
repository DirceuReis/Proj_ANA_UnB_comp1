# Essa função avalia invariância de escala para os 1º e 2º momentos e
# gera gráficos para análise
fun_scale_invariance <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                 na.accept,                 # limite de falhas por ano [0,1]
                                 min.years,                 # mínimo de anos em uma série
                                 which.moment = 1,          # qual momento usar para plotar 
                                 which.duration = c(1, 24), # intervalo de durações
                                 font.family = "serif",     # fonte gráficos
                                 font.size = 12,            # tamanho fonte
                                 plot.dim = c(6,6)          # nrow e ncol p/ gráfico c/ todas as estações
                                 ){
  
  # Pacotes
  if(!require(pacman)){
    install.packages("pacman")
    message("Instalando gerenciador de pacotes 'pacman'...")
  }
  pacman::p_load(pacman, tidyverse, patchwork, ggplot2, ggtext, ggforce)
  

  # MOMENTOS --------------------------------------------------------------
  
  # Selecionar estações adequadas
  gauges <- df.imax %>% 
    select(gauge_code, year, na_prct) %>% 
    group_by(gauge_code, year) %>% 
    reframe(na_prct = first(na_prct)) %>% 
    filter(na_prct <= na.accept) %>% 
    count(gauge_code) %>% 
    filter(n >= min.years) %>% 
    pull(gauge_code)
  
  data <- df.imax[df.imax$na_prct <= na.accept & df.imax$gauge_code %in% gauges,]
  
  # Calcular 1º e 2º momentos centrais locais (p/ cada estação)
  df.mom.local <- data %>% 
    group_by(gauge_code, d) %>%                # agrupar por duração e estação
    reframe(mom1 = mean(imax, na.rm = TRUE),   # valor esperado de Id:  E[Id^1]
            mom2 = mean(imax^2, na.rm = TRUE), # valor esperado de Id²: E[Id^2]
            n_years = n())
  
  # Calcular 1º e 2º momentos centrais regionais
  df.mom.regional <- data %>% 
    group_by(d) %>%                            # agrupar por duração
    reframe(mom1 = mean(imax, na.rm = TRUE),   # valor esperado de Id:  E[Id^1]
            mom2 = mean(imax^2, na.rm = TRUE), # valor esperado de Id²: E[Id^2]
            n_years = n())
  
  df.mom.regional$gauge_code <- "Regional"           # adicionar coluna
  df.mom <- bind_rows(df.mom.local, df.mom.regional) # juntar
  df.mom <- df.mom[between(df.mom$d, which.duration[1], which.duration[2]),] # filtrar durações especificadas
  

  # REGRESSÃO -------------------------------------------------------------

  # # Regressão entre momentos e duração de cada estação
  ls.mom <- split(x = df.mom, f = df.mom$gauge_code) # dividir em lista
  scale.moments <- lapply(X = ls.mom, FUN = function(gauge){
    
    name <- gauge$gauge_code[1] # código da estação
    mom1 <- gauge$mom1
    mom2 <- gauge$mom2
    n.years <- gauge$n_years
    log.mom1 <- log(mom1)     # vetor c/ logs do primeiro momento
    log.mom2 <- log(mom2)     # vetor c/ logs do primeiro momento
    log.d <- log(gauge$d)     # vetor c/ logs das durações
    regression1 <- lm(log.mom1~log.d)        # regressão linear
    regression2 <- lm(log.mom2~log.d)        # regressão linear
    scale1 <- coef(regression1)[["log.d"]]   # coeficiente
    scale2 <- coef(regression2)[["log.d"]]   # coeficiente
    r.squared1 <- summary(regression1)$r.squared # regressão (pegar $sigma tambem)
    r.squared2 <- summary(regression2)$r.squared # regressão
    
    data.frame("gauge_code" = name,
               "d" = gauge$d,
               "mom1" = mom1,
               "mom2" = mom2,
               "coef1" = scale1,
               "coef2" = scale2,
               "rsquared1" = r.squared1,
               "rsquared2" = r.squared2,
               "n_years" = n.years)
    
  }) # fim 'scale.moments'
  
  # Transformar em tbl_df
  scale.moments <- bind_rows(scale.moments) %>% 
    pivot_longer(cols = -c(gauge_code, d, n_years), # todas as colunas menos essas
                 names_to = c(".value", "which_moment"),
                 names_pattern = "(mom|coef|rsquared)(\\d+)") %>% 
    mutate(scale = -coef/as.numeric(which_moment),  # coeficiente de escala
           which_moment = as.numeric(which_moment)) # qual momento calculado
  

  # VISUALIZAÇÃO ----------------------------------------------------------
  
  # Coeficiente de escala regional de cada momento
  scale.regional <-  scale.moments %>% 
    filter(gauge_code == "Regional") %>% 
    group_by(gauge_code, which_moment) %>% 
    reframe(scale = first(scale))
  
  # R² mínimo
  min.rsquared <- min(scale.moments$rsquared)
  
  # Gráfico omparação geral entre modelos regionai e local
  plot.scale.all <- ({
    
    # Colocar scale/r regionais no gráfico
    plot1 <- 
      ggplot(scale.moments, aes(x = d, y = mom, color = gauge_code)) +
      facet_wrap(~which_moment, nrow = 2, strip.position = "left",          # colocar título do facet à esquerda
                 labeller = as_labeller(c("1" = "q = 1", "2" = "q = 2"))) + # reformatar rótulos dos momentos    
      # Modelos locais
      geom_smooth(data = subset(scale.moments, gauge_code != "Regional"), method = lm, formula = y~x, se = FALSE,
                  alpha = 0.6, linewidth = 0.4) +
      # Modelos regionais
      geom_smooth(data = subset(scale.moments, gauge_code == "Regional"), method = lm, formula = y~x, se = FALSE,
                  aes(color = "Regional"), alpha = 0.6, linewidth = 0.5) +
      geom_point(data = subset(scale.moments, gauge_code == "Regional"), aes(shape = "Regional")) +
      ggrepel::geom_text_repel(data = scale.regional, # ggrepel ajuda no posionamento dinâmico do 'geom_text'
                               aes(x = Inf, y = Inf, label = paste("H =", round(scale, 4))),
                               family = "serif", size = 3.5, hjust = 1, vjust = 1) +
      # scale_x_log10(breaks = c(0.5, 1, 2, 6, 24, 72, 240)) +
      # scale_x_log10(breaks = c(1, 2, 4, 6, 12, 24)) +
      scale_x_log10() +
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
    
    # scaleistogramas de R^2
    # Adicionar % de estaçoes acima de R^2 = 0.99
    plot2 <- 
      ggplot(scale.moments, aes(x = rsquared)) +
      facet_wrap(~which_moment, nrow = 2) +
      geom_histogram(aes(y = after_stat(count/sum(count))), bins = 20, fill = "yellow2", color = "darkorange3") +
      # geom_vline(xintercept = 0.99, color = "red", linetype = "dashed", linewidth = 0.5) +
      scale_y_continuous(labels = scales::percent) +
      # scale_x_continuous(limits = c(min.rsquared, 1)) +
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
    
  }) # fim 'plot.scale.all
  
  # Gráfico avaliação individual dos modelos locais vs. regional
  # Dividir dados em 'local' e 'regional'
  data.local <- scale.moments %>% filter(gauge_code != "Regional", which_moment == which.moment)
  data.regional <- scale.moments %>% filter(gauge_code == "Regional", which_moment == which.moment) %>% select(!gauge_code)
  
  # Rótulos
  labels <- data.local %>% 
    group_by(gauge_code) %>% 
    summarise(scale_local = mean(scale),      # coeficiente de escala da estação
              n_years = first(n_years),       # número de anos da estação
              rsquared = first(rsquared)) %>% # R² resultante da regressão na estação
    mutate(scale_regional = scale.regional[["scale"]][which.moment]) # coeficiente de escala regional
  
  # Avaliar dimensões do gráfico
  plot.group <- ({
    
    ggplot(mapping = aes(x = d, y = mom)) +
      ggforce::facet_wrap_paginate(~gauge_code, nrow = plot.dim[1], ncol = plot.dim[2]) +
      geom_smooth(data = data.regional, aes(linetype = "Regional", color = "Regional"),
                  method = lm, formula = y~x, se = FALSE, alpha = 0.6, linewidth = 0.3) +
      geom_smooth(data = data.local, aes(linetype = "Local", color = "Local"),
                  method = lm, formula = y~x, se = FALSE, alpha = 0.6, linewidth = 0.3) +
      geom_point(data = data.local, aes(shape = "Momentos observados"),
                 color = "steelblue4", alpha = 0.5, size = 1.1) +
      geom_text(data = labels,
                aes(x = Inf, y = Inf, label = sprintf("H = %.4f", scale_regional)),
                hjust = 1.1, vjust = 2, size = 2.3, color = "grey10", family = "mono", fontface = "bold") +
      geom_text(data = labels,
                aes(x = Inf, y = Inf, label = sprintf("H = %.4f", scale_local)),
                hjust = 1.1, vjust = 3.5, size = 2.3, color = "steelblue4", family = "mono", fontface = "bold") +
      geom_text(data = labels,
                aes(x = Inf, y = Inf, label = sprintf("N = %.0f anos", n_years)),
                hjust = 1.1, vjust = 5, size = 2.3, color = "steelblue4", family = "mono", fontface = "bold") +
      geom_text(data = labels,
                aes(x = Inf, y = Inf, label = sprintf("R² = %.3f", rsquared)),
                hjust = 1.1, vjust = 6.5, size = 2.3, color = "steelblue4", family = "mono", fontface = "bold") +
      # scale_x_log10(breaks = c(1, 2, 4, 6, 12, 24)) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(name = "", values = c("Regional" = "grey10", "Local" = "steelblue")) +
      scale_linetype_manual(name = "", values = c("Regional" = "dashed", "Local" = "solid")) +
      scale_shape_manual(name = "", values = c("Momentos observados" = 17)) +
      labs(y = "E[I<sub>d</sub><sup>q</sup>]", x = "Durações [h]") +
      theme_minimal() +
      theme(legend.position =  "bottom",
            legend.margin = margin(t = -10),
            axis.title.y = ggtext::element_markdown(), # formatação da legenda do eixo-y 
            strip.placement = "outside",               # fora dos rotulos do eixo-y
            strip.text = element_text(size = 6),       # 'strip.' altera configurações do facet_wrap
            plot.background = element_rect(color = "white"),
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.minor = element_blank(),
            text = element_text(family = font.family, color = "black", size = 10)) +
      guides(shape = guide_legend(order = 2))
    
  }) # fim 'plot.group'
  
  plot.seq <- 1:n_pages(plot.group) # número de páginas necessárias para plotar o gráfico completo
  
  # Dividir o gráfico acima em 'n_pages' e salvar cada um em uma lista
  plot.scale.each <- lapply(X = plot.seq, FUN = function(group){
    plot.group + ggforce::facet_wrap_paginate(~gauge_code, nrow = plot.dim[1], ncol = plot.dim[2], page = group, scales = "fixed")
  }) # fim 'plot.scale.each'
  
  
  # RESULTADOS ------------------------------------------------------------
  
  # Dataframe com coeficientes de escala
  scale.coefficient <- scale.moments %>% 
    group_by(gauge_code, which_moment) %>% 
    reframe(n_years = first(n_years),
              scale = first(scale),
              rsquared = first(rsquared)) %>% 
    mutate(d = deparse(which.duration))
  
  # Organizar resultados
  res <- list("scale.moments" = scale.moments[c("gauge_code", "d", "which_moment", "mom")], # contém só os momentos de cada duração
              "scale.coefficient" = scale.coefficient, # contém coeficiente de escala e informações sobre a regressão
              "plot.scale.all" = plot.scale.all,       # gráfico para os dois momentos com todos os modelos + histograma de R²
              "plot.scale.each" = plot.scale.each)     # gráficos individuais de um momento escolhido com modelo local + regional
  
  return(res)

}
