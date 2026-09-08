library(lubridate)
library(purrr)
library(arrow)
library(stringr)
library(dplyr)
library(trend)
library(ggplot2)
library(Kendall)
library(readxl)
library(sf)
library(geobr)
library(tidyr)
library(ggspatial)
library(terra)
library(geodata)


#analise da tendencia 
#leitura dos dados
dados_maximo_3d <- read_excel("C:/Users/bruno/OneDrive/Área de Trabalho/Bruno/UnB/Tendências/Código_original_Saulo/Teste_resultados/indices_precipitacao_anuais.xlsx") %>%
  mutate(
    gauge_code = as.character(gauge_code),
    year = as.integer(year)
  )

max_3d <- dados_maximo_3d %>%
  group_by(year) %>%
  summarise(max_chuva_3dias = max(chuva_max_3dias, na.rm = TRUE))

#transformar estações em ponto
estacoes <- dados_maximo_3d %>%
  filter(
    !is.na(lat),
    !is.na(long)
  ) %>%
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )

#associar cada estação a uma RH
estacoes_rh <- st_join(
  estacoes,
  regioes_hidrograficas %>%
    select(
      RHI_SG,
      RHI_CD,
      RHI_NM
    ),
  join = st_within
)

#criar coluna com o nome da RH
estacoes_rh <- estacoes_rh %>%
  rename(
    regiao_hidrografica = RHI_NM
  )

#incorporar regiao no dataframe


dados_maximo_3d <- dados_maximo_3d %>%
  left_join(
    estacoes_rh %>%
      st_drop_geometry() %>%
      select(
        gauge_code,
        regiao_hidrografica
      ) %>%
      distinct(
        gauge_code,
        .keep_all = TRUE
      ),
    by = "gauge_code"
  )


#filtrar 8 anos
#dados_maximo_subdiario <- dados_maximo_subdiario %>%
#filter(!is.na(chuva_max_horaria)) %>%
#group_by(gauge_code) %>%
#filter(n() >= 8) %>%
#ungroup()


#1. MANN-KENDALL CLÁSSICO
# MK - MTFPW - FDR
# Souza & Reis Jr. (2022)
#
# Aplicação:
# dados_maximo_diario
# gauge_code = estação
# year       = ano
# chuva_max  = índice anual
#




# 1. MANN-KENDALL


mann_kendall_java <- function(
    x,
    nivelSignificancia = 5
) {
  
  
  # Preparar série
  
  
  x <- as.numeric(x)
  
  x <- x[
    is.finite(x)
  ]
  
  N <- length(x)
  
  
  if (N < 2) {
    
    return(
      list(
        Z = NA_real_,
        pvalue = NA_real_,
        significativo = NA
      )
    )
    
  }
  
  
  
  # Postos
  
  
  R <- rank(
    x,
    ties.method = "min"
  )
  
  
  
  # Estatística S
  #
  # Equação (2)
  
  
  S <- 0
  
  for (i in 1:(N - 1)) {
    
    for (j in (i + 1):N) {
      
      S <- S +
        sign(
          R[j] - R[i]
        )
      
    }
    
  }
  
  
  
  # Variância de S
  #
  # Equação (4)
  
  
  tab <- table(x)
  
  t <- tab[
    tab > 1
  ]
  
  
  parte_sem_empates <-
    N *
    (N - 1) *
    (2 * N + 5)
  
  
  parte_empates <- sum(
    t *
      (t - 1) *
      (2 * t + 5)
  )
  
  
  varS <- (
    parte_sem_empates -
      parte_empates
  ) / 18
  
  
  sdS <- sqrt(varS)
  
  
  
  # Estatística Z
  #
  # Equação (5)
  
  
  if (S > 0) {
    
    Z <- (
      S - 1
    ) / sdS
    
  } else if (S < 0) {
    
    Z <- (
      S + 1
    ) / sdS
    
  } else {
    
    Z <- 0
    
  }
  
  
  
  # p-valor bilateral
  
  
  p <- 2 *
    (
      1 -
        pnorm(
          abs(Z)
        )
    )
  
  
  alpha <- nivelSignificancia / 100
  
  
  
  # Resultado
  
  
  list(
    Z = Z,
    pvalue = p,
    significativo = p < alpha
  )
  
}




# 2. SEN'S SLOPE
# Para o MTFPW é preciso da tendência estimada para
# retirar a tendência SOMENTE durante a estimativa de r1.


sens_slope <- function(x) {
  
  x <- as.numeric(x)
  
  x <- x[
    is.finite(x)
  ]
  
  n <- length(x)
  
  
  if (n < 2) {
    return(NA_real_)
  }
  
  
  slopes <- c()
  
  
  for (i in 1:(n - 1)) {
    
    for (j in (i + 1):n) {
      
      slopes <- c(
        slopes,
        (
          x[j] - x[i]
        ) / (
          j - i
        )
      )
      
    }
    
  }
  
  
  median(
    slopes,
    na.rm = TRUE
  )
  
}




# 3. MTFPW
#
# Modified Trend-Free Pre-Whitening
#
# Souza & Reis Jr. (2022)
#
# Etapas:
#
# 1. Estimar tendência
# 2. Detrendar a série
# 3. Estimar r1
# 4. Corrigir viés de r1
# 5. Testar autocorrelação
# 6. Se não significativa:
#       retorna série original
# 7. Se significativa:
#       aplica PW na série ORIGINAL
#


mtfpw_java <- function(x) {
  
  
  # Preparar série
  
  
  x <- as.numeric(x)
  
  x <- x[
    is.finite(x)
  ]
  
  n <- length(x)
  
  
  if (n < 5) {
    
    return(
      list(
        serie = x,
        beta = NA_real_,
        r1 = NA_real_,
        r1_corrigido = NA_real_,
        autocorrelacao = FALSE,
        aplicada = FALSE
      )
    )
    
  }
  
  
  
  # 3.1 SEN'S SLOPE
  
  
  beta <- sens_slope(x)
  
  
  
  # 3.2 DETRENDING
  #
  # O artigo apresenta:
  #
  # Xd_t = Xt - beta_hat * X_bar * t / 10
  #
  # Como a beta da Eq. (7) é uma inclinação relativa,
  # beta * X_bar / 10 equivale à inclinação absoluta.
  #
  # Portanto:
  #
  # Xd = X - SenSlope * t
  
  
  t <- seq_len(n)
  
  x_detr <- x -
    beta * t
  
  
  
  
  # 3.3 AUTOCORRELAÇÃO LAG-1
  #
  # Equação (9)
  #
  # Primeiro calcular r1.
  
  
  x_media <- mean(
    x_detr
  )
  
  
  numerador <- sum(
    (
      x_detr[-n] -
        x_media
    ) *
      (
        x_detr[-1] -
          x_media
      )
  )
  
  
  denominador <- sum(
    (
      x_detr -
        x_media
    )^2
  )
  
  
  r1 <- numerador /
    denominador
  
  
  
  if (!is.finite(r1)) {
    
    return(
      list(
        serie = x,
        beta = beta,
        r1 = NA_real_,
        r1_corrigido = NA_real_,
        autocorrelacao = FALSE,
        aplicada = FALSE
      )
    )
    
  }
  
  
  
  # 3.4 CORREÇÃO DO VIÉS
  #
  # Equação (10)
  #
  # r1* = (n*r1 + 2)/(n - 4)
  
  
  r1_corrigido <- (
    n * r1 + 2
  ) / (
    n - 4
  )
  
  
  
  
  # 3.5 INTERVALO DE CONFIANÇA
  #
  # Equação (8)
  #
  # -1 - 1.96 sqrt(n-2)
  # --------------------
  #       n - 1
  #
  # <= rho1 <=
  #
  # -1 + 1.96 sqrt(n-2)
  # --------------------
  #       n - 1
  
  
  lim_inf <- (
    -1 -
      1.96 *
      sqrt(
        n - 2
      )
  ) / (
    n - 1
  )
  
  
  lim_sup <- (
    -1 +
      1.96 *
      sqrt(
        n - 2
      )
  ) / (
    n - 1
  )
  
  
  
  
  # 3.6 VERIFICAR AUTOCORRELAÇÃO
  #
  # Se r1 estiver dentro do intervalo:
  # série considerada independente.
  
  
  autocorrelacao <- (
    r1 < lim_inf ||
      r1 > lim_sup
  )
  
  
  
  
  # 3.7 SEM AUTOCORRELAÇÃO
  #
  # O artigo determina que o MK seja aplicado à série
  # ORIGINAL.
  
  
  if (!autocorrelacao) {
    
    return(
      list(
        serie = x,
        beta = beta,
        r1 = r1,
        r1_corrigido = r1_corrigido,
        autocorrelacao = FALSE,
        aplicada = FALSE
      )
    )
    
  }
  
  
  
  
  # 3.8 MTFPW
  #
  # IMPORTANTE:
  #
  # O PW é aplicado à SÉRIE ORIGINAL.
  #
  # Y_t = X_t - r1* X_(t-1)
  #
  
  
  x_mtfpw <- (
    x[-1] -
      r1_corrigido *
      x[-n]
  )
  
  
  
  
  # 3.9 RETORNO
  
  
  return(
    list(
      serie = x_mtfpw,
      beta = beta,
      r1 = r1,
      r1_corrigido = r1_corrigido,
      autocorrelacao = TRUE,
      aplicada = TRUE
    )
  )
  
}




# 4. BENJAMINI-HOCHBERG
#
# FDR
#
# Implementação seguindo a lógica do artigo:
#
# P(j) < d(j)
#
# d(j) = j * alpha / m
#
# Encontrar o MAIOR j que satisfaz a condição.
# Todas as hipóteses até j são rejeitadas.
#


fdr_java <- function(
    p,
    alpha = 0.05
) {
  
  
  
  # Resultado inicial
  
  
  sig <- rep(
    FALSE,
    length(p)
  )
  
  
  
  # P-valores válidos
  
  
  validos <- is.finite(p)
  
  
  if (!any(validos)) {
    
    return(sig)
    
  }
  
  
  p_validos <- p[
    validos
  ]
  
  
  m <- length(
    p_validos
  )
  
  
  
  # Ordenar p-valores
  
  
  ordem <- order(
    p_validos
  )
  
  
  p_ordenados <- p_validos[
    ordem
  ]
  
  
  
  # Valores críticos
  #
  # d_i = i * alpha / m
  
  
  i <- seq_len(m)
  
  
  d <- (
    i *
      alpha /
      m
  )
  
  
  
  # Encontrar o maior j
  #
  # P(j) < d(j)
  
  
  candidatos <- which(
    p_ordenados < d
  )
  
  
  
  # Nenhum significativo
  
  
  if (length(candidatos) == 0) {
    
    return(sig)
    
  }
  
  
  
  # Maior j
  
  
  j <- max(
    candidatos
  )
  
  
  
  # Rejeitar todas as hipóteses até j
  
  
  sig_ordenado <- rep(
    FALSE,
    m
  )
  
  
  sig_ordenado[
    1:j
  ] <- TRUE
  
  
  
  # Retornar à ordem original
  
  
  sig_validos <- rep(
    FALSE,
    m
  )
  
  
  sig_validos[
    ordem
  ] <- sig_ordenado
  
  
  sig[
    validos
  ] <- sig_validos
  
  
  return(sig)
  
}





# 5. APLICAÇÃO POR REGIÃO HIDROGRÁFICA:
#
# MK
# MTFPW
# MK após MTFPW



# Lista das Regiões Hidrográficas
RHs <- dados_maximo_3d %>%
  
  filter(
    !is.na(regiao_hidrografica)
  ) %>%
  
  distinct(
    regiao_hidrografica
  ) %>%
  
  pull(
    regiao_hidrografica
  )


# Verificar as regiões encontradas

cat(
  "\n============================================\n"
)

cat(
  "REGIÕES HIDROGRÁFICAS - QX3d\n"
)

cat(
  "============================================\n"
)

print(RHs)

cat(
  "Total de regiões:",
  length(RHs),
  "\n"
)

cat(
  "============================================\n"
)




# FUNÇÃO PARA ANALISAR UMA RH POR VEZ



analisar_RH_3d <- function(nome_RH) {
  
  
  cat(
    "\n--------------------------------------------\n"
  )
  
  cat(
    "Analisando:",
    nome_RH,
    "\n"
  )
  
  cat(
    "--------------------------------------------\n"
  )
  
  
  

  # FILTRAR SOMENTE A RH ATUAL

  
  
  dados_RH <- dados_maximo_3d %>%
    
    filter(
      regiao_hidrografica == nome_RH
    )
  
  
  cat(
    "Estações:",
    n_distinct(
      dados_RH$gauge_code
    ),
    "\n"
  )
  
  
  

  # MK + MTFPW

  
  
  resultado_RH <- dados_RH %>%
    
    arrange(
      gauge_code,
      year
    ) %>%
    
    group_by(
      gauge_code
    ) %>%
    
    summarise(
      
      

      # REGIÃO HIDROGRÁFICA

      
      regiao_hidrografica =
        first(
          regiao_hidrografica
        ),
      
      
      

      # TAMANHO DA SÉRIE

      
      n_anos =
        sum(
          is.finite(
            chuva_max_3dias
          )
        ),
      
      
      

      # MK ORIGINAL

      
      mk_original =
        list(
          mann_kendall_java(
            chuva_max_3dias
          )
        ),
      
      
      Z_MK =
        mk_original[[1]]$Z,
      
      
      p_MK =
        mk_original[[1]]$pvalue,
      
      
      sig_MK =
        mk_original[[1]]$significativo,
      
      
      dir_MK =
        case_when(
          
          Z_MK > 0 ~
            "Crescente",
          
          Z_MK < 0 ~
            "Decrescente",
          
          TRUE ~
            "Sem tendência"
          
        ),
      
      
      

      # MTFPW

      
      mtfpw =
        list(
          mtfpw_java(
            chuva_max_3dias
          )
        ),
      
      
      

      # INFORMAÇÕES DA AUTOCORRELAÇÃO

      
      beta_MTFPW =
        mtfpw[[1]]$beta,
      
      
      r1_MTFPW =
        mtfpw[[1]]$r1,
      
      
      r1_corrigido =
        mtfpw[[1]]$r1_corrigido,
      
      
      autocorrelacao =
        mtfpw[[1]]$autocorrelacao,
      
      
      MTFPW_aplicado =
        mtfpw[[1]]$aplicada,
      
      
      

      # MK APÓS MTFPW

      
      mk_mtfpw =
        list(
          mann_kendall_java(
            mtfpw[[1]]$serie
          )
        ),
      
      
      Z_MTFPW =
        mk_mtfpw[[1]]$Z,
      
      
      p_MTFPW =
        mk_mtfpw[[1]]$pvalue,
      
      
      sig_MTFPW =
        mk_mtfpw[[1]]$significativo,
      
      
      dir_MTFPW =
        case_when(
          
          Z_MTFPW > 0 ~
            "Crescente",
          
          Z_MTFPW < 0 ~
            "Decrescente",
          
          TRUE ~
            "Sem tendência"
          
        )
      
    ) %>%
    
    ungroup()
  
  
  
  

  # FDR
  #
  # IMPORTANTE:
  # O FDR é aplicado SOMENTE aos p-valores
  # das estações da RH atual.

  
  
  resultado_RH$sig_FDR <-
    
    fdr_java(
      resultado_RH$p_MTFPW,
      alpha = 0.05
    )
  
  
  
  

  # DIREÇÃO FINAL

  
  
  resultado_RH <- resultado_RH %>%
    
    mutate(
      
      
      sig_final =
        ifelse(
          sig_FDR,
          "S(5)",
          "NS"
        ),
      
      
      dir_final =
        case_when(
          
          sig_FDR &
            Z_MTFPW > 0 ~
            "Crescente",
          
          sig_FDR &
            Z_MTFPW < 0 ~
            "Decrescente",
          
          TRUE ~
            "Sem tendência"
          
        )
      
    )
  
  
  
  

  # RETORNAR RESULTADO DA RH

  
  
  return(
    resultado_RH
  )
  
}




# 6. EXECUTAR UMA RH POR VEZ



resultados_RH_3d <- lapply(
  
  RHs,
  
  analisar_RH_3d
  
)




# 7. JUNTAR OS RESULTADOS DAS 12 RHs



resultado_3d <- bind_rows(
  resultados_RH_3d
)


# Organizar os resultados

resultado_3d <- resultado_3d %>%
  
  arrange(
    regiao_hidrografica,
    gauge_code
  )




# 8. RESUMO DOS RESULTADOS



cat(
  "\n============================================\n"
)

cat(
  "RESULTADOS QX3d - MK - MTFPW - FDR\n"
)

cat(
  "============================================\n"
)

cat(
  "Regiões hidrográficas analisadas:",
  n_distinct(
    resultado_3d$regiao_hidrografica
  ),
  "\n"
)

cat(
  "Estações analisadas:",
  nrow(
    resultado_3d
  ),
  "\n"
)

cat(
  "Estações com MK significativo:",
  sum(
    resultado_3d$sig_MK,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Estações com MTFPW + MK significativo:",
  sum(
    resultado_3d$sig_MTFPW,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Estações significativas após FDR:",
  sum(
    resultado_3d$sig_FDR,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "MTFPW aplicado:",
  sum(
    resultado_3d$MTFPW_aplicado,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Estações sem autocorrelação significativa:",
  sum(
    !resultado_3d$autocorrelacao,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "============================================\n"
)




# 9. RESUMO POR REGIÃO HIDROGRÁFICA



resumo_RH_3d <- resultado_3d %>%
  
  group_by(
    regiao_hidrografica
  ) %>%
  
  summarise(
    
    n_estacoes =
      n(),
    
    MK_significativo =
      sum(
        sig_MK,
        na.rm = TRUE
      ),
    
    MTFPW_significativo =
      sum(
        sig_MTFPW,
        na.rm = TRUE
      ),
    
    FDR_significativo =
      sum(
        sig_FDR,
        na.rm = TRUE
      ),
    
    MTFPW_aplicado =
      sum(
        MTFPW_aplicado,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


print(
  resumo_RH_3d
)




#CRESCENTES E DECRESCENTES APÓS FDR



resultado_3d %>%
  
  filter(
    sig_FDR
  ) %>%
  
  count(
    regiao_hidrografica,
    dir_final
  )

#10.Gráfico contagem

# 10.1 CONTAGEM DE SIGNIFICATIVAS E NÃO SIGNIFICATIVAS


dados_grafico <- tibble(
  
  metodo = c(
    "MK",
    "MK + MTFPW",
    "MK + MTFPW + FDR"
  ),
  
  Significativa = c(
    sum(resultado_3d$sig_MK, na.rm = TRUE),
    sum(resultado_3d$sig_MTFPW, na.rm = TRUE),
    sum(resultado_3d$sig_FDR, na.rm = TRUE)
  )
) %>%
  
  mutate(
    Total = nrow(resultado_3d),
    `Não significativa` = Total - Significativa
  ) %>%
  
  select(
    metodo,
    Significativa,
    `Não significativa`
  ) %>%
  
  pivot_longer(
    cols = c(
      Significativa,
      `Não significativa`
    ),
    names_to = "significancia",
    values_to = "n"
  ) %>%
  
  mutate(
    metodo = factor(
      metodo,
      levels = c(
        "MK",
        "MK + MTFPW",
        "MK + MTFPW + FDR"
      )
    )
  )



#GRÁFICO


grafico_significancia <- ggplot(
  dados_grafico,
  aes(
    x = metodo,
    y = n,
    fill = significancia
  )
) +
  
  geom_col(
    width = 0.65
  ) +
  
  geom_text(
    aes(
      label = n
    ),
    position = position_stack(vjust = 0.5),
    size = 5
  ) +
  
  scale_fill_manual(
    values = c(
      "Significativa" = "red",
      "Não significativa" = "grey80"
    ),
    name = "Significância"
  ) +
  scale_y_continuous(
    breaks = seq(0, 6000, by = 1000),
    limits = c(0, 6000),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = NULL,
    y = "Número de estações",
    title = "Estações com e sem tendência significativa",
    subtitle = "Comparação entre MK, MTFPW e FDR"
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.title.y = element_text(
      size = 16
    ),
    
    axis.text.x = element_text(
      size = 14
    ),
    
    axis.text.y = element_text(
      size = 13
    ),
    
    legend.title = element_text(
      size = 16
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    plot.title = element_text(
      size = 17,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 13
    )
  )

grafico_significancia

#salvar
ggsave(
  filename = "grafico_sign_metodos_3d_RH.png",
  plot = grafico_significancia,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)



#11. mapa sig e nao sig Brasil
# 1. COORDENADAS DAS ESTAÇÕES


coords_estacoes <- dados_maximo_3d %>%
  
  select(
    gauge_code,
    lat,
    long
  ) %>%
  
  filter(
    !is.na(lat),
    !is.na(long)
  ) %>%
  
  distinct(
    gauge_code,
    .keep_all = TRUE
  )



# 2. JUNTAR COORDENADAS AO RESULTADO DA TENDÊNCIA


resultado_mapa_3d <- resultado_3d %>%
  
  select(
    gauge_code,
    sig_MK,
    sig_MTFPW,
    sig_FDR,
    dir_MK,
    dir_MTFPW
  ) %>%
  
  left_join(
    coords_estacoes,
    by = "gauge_code"
  )



# 3. CLASSIFICAÇÃO FINAL


resultado_mapa_3d <- resultado_mapa_3d %>%
  
  mutate(
    
    tendencia = case_when(
      
      sig_FDR & dir_MTFPW == "Crescente" ~
        "Significativa - crescente",
      
      sig_FDR & dir_MTFPW == "Decrescente" ~
        "Significativa - decrescente",
      
      TRUE ~
        "Não significativa"
    )
  )



# 4. CONFERIR NÚMERO DE ESTAÇÕES


cat(
  "Estações no resultado:",
  n_distinct(resultado_mapa_3d$gauge_code),
  "\n"
)

cat(
  "Estações com coordenadas:",
  n_distinct(
    resultado_mapa_3d$gauge_code[
      !is.na(resultado_mapa_3d$lat) &
        !is.na(resultado_mapa_3d$long)
    ]
  ),
  "\n"
)



# 5. CONVERTER ESTAÇÕES PARA SF


estacoes_sf <- resultado_mapa_3d %>%
  
  filter(
    !is.na(lat),
    !is.na(long)
  ) %>%
  
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )



# 6. BAIXAR LIMITES ADMINISTRATIVOS DO BRASIL


cat("Baixando limites dos estados...\n")


estados_brasil <- geodata::gadm(
  country = "BRA",
  level = 1,
  path = tempdir()
)



# 7. CONVERTER PARA SF


estados_brasil <- st_as_sf(
  estados_brasil
)



# 8. TRANSFORMAR PARA WGS84


estados_brasil <- st_transform(
  estados_brasil,
  4326
)



# 9. VERIFICAR


cat(
  "Número de estados:",
  nrow(estados_brasil),
  "\n"
)



# 10. MAPA

# MAPA DE TODAS AS ESTAÇÕES


mapa_tendencia_3d <- ggplot() +
  
  
  
  # ESTADOS DO BRASIL
  
  
  geom_sf(
    data = estados_brasil,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  
  
  
  # ESTAÇÕES
  
  
  geom_sf(
    data = estacoes_sf,
    aes(
      shape = tendencia,
      fill = tendencia
    ),
    size = 2.5,
    color = "black",
    stroke = 0.4
  ) +
  
  
  
  # SÍMBOLOS
  
  
  scale_shape_manual(
    
    values = c(
      "Significativa - crescente" = 24,
      "Significativa - decrescente" = 25,
      "Não significativa" = 21
    ),
    
    name = "Tendência"
  ) +
  
  
  
  # PREENCHIMENTO
  
  
  scale_fill_manual(
    
    values = c(
      "Significativa - crescente" = "red",
      "Significativa - decrescente" = "red",
      "Não significativa" = "blue"
    ),
    
    name = "Tendência"
  ) +
  
  
  
  # NORTE
  
  
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm"),
    pad_x = unit(0.7, "cm"),
    pad_y = unit(1.5, "cm"),
    style = north_arrow_fancy_orienteering(
      fill = c("black", "white"),
      text_col = "black"
    )
  ) +
  
  
  
  # ESCALA
  
  
  annotation_scale(
    location = "bl",
    width_hint = 0.20,
    height = unit(0.25, "cm"),
    text_cex = 0.7,
    line_width = 0.7,
    pad_x = unit(0.7, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  
  
  
  # LIMITES DO MAPA
  
  
  coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE
  ) +
  
  
  
  # TÍTULO
  
  
  labs(
    title = "Tendências nas precipitações máximas anuais",
    subtitle = "MK + MTFPW + FDR (α = 5%)"
  ) +
  
  
  
  # TEMA
  
  
  theme_minimal() +
  
  theme(
    
    
    # GRADE
    
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    
    
    # CONTORNO CINZA
    
    
    
    
    # EIXOS
    
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    axis.text.x = element_text(
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 13
    ),
    
    
    
    # LEGENDA
    
    
    legend.title = element_text(
      size = 16
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    
    
    # TÍTULO
    
    
    plot.title = element_text(
      size = 17,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 13
    )
  )



# MOSTRAR MAPA


mapa_tendencia_3d

#salvar
ggsave(
  filename = "mapa_tendencias_3d_.png",
  plot = mapa_tendencia_3d,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)


#MAPA SÓ COM AS SIG
# FILTRAR SOMENTE AS ESTAÇÕES SIGNIFICATIVAS


estacoes_sig <- estacoes_sf %>%
  filter(
    sig_FDR == TRUE
  ) %>%
  mutate(
    tendencia_sig = case_when(
      dir_MTFPW == "Crescente" ~ "Crescente",
      dir_MTFPW == "Decrescente" ~ "Decrescente",
      TRUE ~ NA_character_
    )
  )



# MAPA


mapa_sig_3d <- ggplot() +
  
  
  
  # ESTADOS DO BRASIL
  
  
  geom_sf(
    data = estados_brasil,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  
  
  
  # ESTAÇÕES SIGNIFICATIVAS
  
  
  geom_sf(
    data = estacoes_sig,
    aes(
      shape = tendencia_sig,
      fill = tendencia_sig
    ),
    size = 3,
    color = "black",
    stroke = 0.4
  ) +
  
  
  
  # SÍMBOLOS
  
  
  scale_shape_manual(
    
    values = c(
      "Crescente" = 24,
      "Decrescente" = 25
    ),
    
    name = "Tendência"
  ) +
  
  
  
  # PREENCHIMENTO
  
  
  scale_fill_manual(
    
    values = c(
      "Crescente" = "red",
      "Decrescente" = "red"
    ),
    
    name = "Tendência"
  ) +
  
  
  
  # NORTE
  
  
  annotation_north_arrow(
    location = "bl",
    which_north = "true",
    height = unit(0.7, "cm"),
    width = unit(0.7, "cm"),
    pad_x = unit(0.7, "cm"),
    pad_y = unit(1.5, "cm"),
    style = north_arrow_fancy_orienteering(
      fill = c("black", "white"),
      text_col = "black"
    )
  ) +
  
  
  
  # ESCALA
  
  
  annotation_scale(
    location = "bl",
    width_hint = 0.20,
    height = unit(0.25, "cm"),
    text_cex = 0.7,
    line_width = 0.7,
    pad_x = unit(0.7, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  
  
  
  # LIMITES DO MAPA
  
  
  coord_sf(
    xlim = c(-74, -34),
    ylim = c(-34, 6),
    expand = FALSE
  ) +
  
  
  
  # TÍTULO
  
  
  labs(
    title = "Estações com tendência significativa",
    subtitle = "MK + MTFPW + FDR (α = 5%)"
  ) +
  
  
  
  # TEMA
  
  
  theme_minimal() +
  
  theme(
    
    
    # retirar grade
    
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    
    
    # CONTORNO CINZA
    
    
    
    
    
    # EIXOS
    
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    axis.text.x = element_text(
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 13
    ),
    
    
    
    # LEGENDA
    
    
    legend.title = element_text(
      size = 16
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    
    
    # TÍTULO
    
    
    plot.title = element_text(
      size = 17,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 13
    )
  )



# MOSTRAR


mapa_sig_3d

#salvar
ggsave(
  filename = "mapa_tendencias_significativas_3d.png",
  plot = mapa_sig_3d,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

