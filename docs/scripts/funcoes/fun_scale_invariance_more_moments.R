# Essa função avalia invariância de escala para os 1º e 2º momentos e
# gera gráficos para análise
fun_scale_invariance <- function(df.imax,                   # data.frame com intensidades máximas anuais para diferentes durações
                                 na.accept,                 # limite de falhas por ano [0,1]
                                 min.years,                 # mínimo de anos em uma série
                                 which.moment = 1:3,       # quais momentos calcular
                                 which.duration = c(1, 24), # intervalo de durações
                                 min.duration = 3,          # nro. mínimo durações p/ regressão
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
  data <- df.imax %>%
    filter(na_prct <= na.accept,              # filtrar anos por porcentagem de falhas
           between(d, which.duration[1], which.duration[2])) %>% # filtrar durações dentro do intervalo
    group_by(gauge_code) %>%                  # agrupar por estação
    filter(n_distinct(d) >= min.duration,     # só estações com pelo menos 'min.duration'
           n_distinct(year) >= min.years) %>% # so estações com pelo menos 'min.years'
    ungroup()                                 # desagrupar
  
  # Calcular 'which.moment' momentos centrais locais (p/ cada estação)
  df.mom.local <- data %>% 
    group_by(gauge_code, d) %>% 
    reframe(n_years = n(),
            purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE))))
  
  # Calcular 'which.moment' momentos centrais regionais
  df.mom.regional <- data %>% 
    group_by(d) %>% 
    reframe(n_years = n(),
            purrr::map_dfc(which.moment, ~ tibble("mom{.x}" := mean(imax^.x, na.rm = TRUE))))
  
  df.mom.regional$gauge_code <- "Regional"                # adicionar coluna
  df.mom <- bind_rows(df.mom.local, df.mom.regional) %>%  # juntar
    pivot_longer(cols = -c(gauge_code, d, n_years),
                 names_to = c(".value", "which_moment"),
                 names_pattern = "(mom)(\\d+)") %>% 
    mutate(which_moment = as.numeric(which_moment))

  # REGRESSÃO -------------------------------------------------------------

  # Regressão entre momentos e duração de cada estação
  ls.mom <- split(x = df.mom, f = df.mom$gauge_code) # dividir em lista
  scale.moments <- lapply(X = ls.mom, FUN = function(gauge){
    
    name <- gauge$gauge_code[1] # código da estação
    
    # Aplicar 
    gauge.mom <- split(x = gauge, f = gauge$which_moment)
    
    regression <- lapply(X = gauge.mom, FUN = function(gg){
      
      which.moment <- gg$which_moment[1] # extrair qual momento
      gg$log_mom <- log(gg$mom)          # logaritmo do momento
      gg$log_d <- log(gg$d)              # logaritmo da duração
      model <- lm(log_mom ~ log_d, data = gg) # regressão linear
      coef <- coef(model)[["log_d"]]     # slope
      scale <- -coef/which.moment        # coeficiente de escala
      summary <- summary(model)          # resumo
      rsquared <- summary$r.squared      # R²
      se <- summary$sigma/which.moment   # erro padrão
      
      res <- data.frame("gauge_code" = name,
                        "n_years" = gg$n_years,
                        "d" = gg$d,
                        "which_moment" = factor(which.moment),
                        "mom" = gg$mom,
                        "coef_reg" = coef,
                        "scale" = scale,
                        "rsquared" = rsquared,
                        "se" = se)
      
    }) # fim 'regression'
    
    # Transformar em tbl_df
    regression <- bind_rows(regression)
    
  })
  
  # Transformar em tbl_df
  scale.moments <- bind_rows(scale.moments)
  

  # VISUALIZAÇÃO ----------------------------------------------------------
  
  # Coeficiente de escala regional de cada momento
  scale.regional <- scale.moments %>% 
    filter(gauge_code == "Regional") %>% 
    group_by(gauge_code, which_moment) %>% 
    reframe(scale = first(scale),
            rsquared = first(rsquared)) %>%
    mutate(scale_vjust = (row_number() - 1)*1.5,
           scale_label = sprintf("q = %d: H = %.4f", which_moment, round(scale, 4)),)
  
  # R² mínimo
  min.rsquared <- min(scale.moments$rsquared)
  median.rsquared <- quantile(scale.moments$rsquared[scale.moments$gauge_code != "Regional"])[[3]]
  
  # Gráfico comparação geral entre modelos regional e local
  plot.scale.all <- ({
    
    # Colocar scale/r regionais no gráfico
    plot1 <- 
      ggplot(scale.moments %>% filter(gauge_code == "Regional"), aes(x = d, y = mom)) +
      geom_smooth(aes(group = which_moment), color = "black", method = "lm", formula = y~x, se = FALSE, linewidth = 0.4) +
      geom_point(aes(shape = factor(which_moment,
                                    labels = paste0("q = ", unique(which_moment), ": H = ", round(unique(scale), 4))))) +
      scale_x_log10() +
      scale_y_log10(labels = scales::label_comma()) +
      labs(subtitle = "a)", y = "E[I<sub>d</sub><sup>q</sup>]", x = "Durações [h]", shape = "")+
      theme_minimal() +
      theme(legend.position =  c(1, 1.08),
            legend.justification = c("right" , "top"),
            legend.text = element_text(size = 9),
            axis.title.y = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
            plot.background = element_rect(color = "white"),
            # aspect.ratio = 1,
            panel.border = element_rect(color = "black", fill = NA),
            text = element_text(family = "serif", color = "black", size = 12)); plot1
    
      # Histogramas de R^2
    plot2 <- 
      ggplot(scale.moments, aes(x = rsquared)) +
      facet_wrap(~which_moment, nrow = length(which.moment), labeller = as_labeller(function(x) paste0("q = ", x))) +
      geom_histogram(aes(y = after_stat(count/sum(count))), bins = 30, fill = "yellow2", color = "darkorange3") +
      geom_vline(xintercept = median.rsquared, color = "red", linetype = "dashed", linewidth = 0.5) +
      scale_y_continuous(labels = scales::percent) +
      # scale_x_continuous(limits = c(min.rsquared, 1)) +
      labs(subtitle = "b)", x = "R²", y = "") +
      theme_minimal() +
      theme(legend.position =  "none",
            strip.placement = "inside",                      # fora dos rotulos do eixo-y
            strip.text = element_text(size = 10),             # 'strip.' altera configurações do facet_wrap
            axis.title.y = ggtext::element_markdown(),        # formatação da legenda do eixo-y 
            plot.background = element_rect(color = "white"),
            panel.border = element_rect(color = "black", fill = NA),
            text = element_text(family = "serif", color = "black", size = 12)); plot2
    
    # Combinar c/ 'patchwork'
    plot1 + plot2 + plot_layout(widths = c(1.5,1))
    
}) # fim 'plot.scale.all'
  
  # Gráfico avaliação individual dos modelos locais vs. regional
  # Dividir dados em 'local' e 'regional'
  data.local <- scale.moments %>% filter(gauge_code != "Regional")
  data.regional <- scale.moments %>% filter(gauge_code == "Regional") %>% select(!gauge_code)
  
  # Rótulos
  labels <- data.local %>% 
    group_by(gauge_code) %>% 
    summarise(scale_local = mean(scale),     # coeficiente de escala da estação
              n_years = first(n_years),      # número de anos da estação
              rsquared = mean(rsquared)) %>% # R² resultante da regressão na estação
    mutate(scale_regional = mean(scale.regional[["scale"]])) # coeficiente de escala regional
  
  # Avaliar dimensões do gráfico
  plot.group <- ({
    
    ggplot(mapping = aes(x = d, y = mom)) +
      ggforce::facet_wrap_paginate(~gauge_code, nrow = plot.dim[1], ncol = plot.dim[2]) +
      geom_smooth(data = data.regional, aes(group = which_moment, linetype = "Regional", color = "Regional"),
                  method = lm, formula = y~x, se = FALSE, alpha = 0.6, linewidth = 0.3) +
      geom_smooth(data = data.local, aes(group = which_moment, linetype = "Local", color = "Local"),
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
      scale_x_log10() +
      scale_y_log10(labels = scales::label_comma()) +
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
            text = element_text(family = "serif", color = "black", size = 10)) +
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
            rsquared = first(rsquared),
            se = first(se)) %>% 
    mutate(d = deparse(which.duration))
  
  # Organizar resultados
  res <- list("scale.moments" = scale.moments[c("gauge_code", "d", "which_moment", "mom")], # contém só os momentos de cada duração
              "scale.coefficient" = scale.coefficient, # contém coeficiente de escala e informações sobre a regressão
              "plot.scale.all" = plot.scale.all,       # gráfico para os dois momentos com todos os modelos + histograma de R²
              "plot.scale.each" = plot.scale.each)     # gráficos individuais de um momento escolhido com modelo local + regional
  
  return(res)

}
