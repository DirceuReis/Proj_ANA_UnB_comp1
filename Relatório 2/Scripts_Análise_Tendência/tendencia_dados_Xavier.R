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
library(ncdf4)




# 1. CAMINHO DO ARQUIVO XAVIER


arquivo_xavier <- paste0(
  "C:/Users/bruno/OneDrive/Área de Trabalho/Bruno/UnB/",
  "Tendências/Código_original_Saulo/Teste_resultados/",
  "Dados_Xavier/pr_1961_2025_BR-DWGD_monthly_v_3.2.4.nc"
)


# 2. ABRIR NETCDF


nc <- nc_open(arquivo_xavier)

print(nc)


# Ver estrutura
print(nc)


# 3. EXTRAIR COORDENADAS E TEMPO


lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
tempo <- ncvar_get(nc, "time")

# Converter tempo para data
datas <- as.POSIXct(
  "1961-01-01 00:00:00",
  tz = "UTC"
) + tempo * 3600

# Verificar
range(datas)
length(lon)
length(lat)
length(datas)


# 4. LER PRECIPITAÇÃO


pr <- ncvar_get(nc, "pr")

dim(pr)
summary(as.vector(pr))

nc_close(nc)


# 5. TRANSFORMAR GRADE EM DATAFRAME


dados_xavier <- expand.grid(
  longitude = lon,
  latitude = lat,
  time = datas
)

dados_xavier$pr <- as.vector(pr)

# Ano
dados_xavier <- dados_xavier %>%
  mutate(
    year = as.integer(format(time, "%Y"))
  )

# Verificar
head(dados_xavier)
dim(dados_xavier)


# 5. PRECIPITAÇÃO ANUAL


anos <- as.integer(format(datas, "%Y"))

anos_unicos <- sort(unique(anos))

length(anos_unicos)


# 6. AGREGAR PRECIPITAÇÃO MENSAL PARA ANUAL


n_lon <- length(lon)
n_lat <- length(lat)
n_anos <- length(anos_unicos)

pr_max_anual <- array(
  NA_real_,
  dim = c(n_lon, n_lat, n_anos)
)

for (i in seq_along(anos_unicos)) {
  
  ano_atual <- anos_unicos[i]
  indices <- which(anos == ano_atual)
  
  dados_ano <- pr[, , indices]
  
  pr_max_anual[, , i] <- apply(
    dados_ano,
    c(1, 2),
    function(x) {
      
      # Se não houver nenhum mês válido
      if (all(is.na(x))) {
        return(NA_real_)
      }
      
      # Maior precipitação mensal daquele ano
      max(x, na.rm = TRUE)
    }
  )
}

dim(pr_max_anual)


# 7. TRANSFORMAR GRADE ANUAL EM DATAFRAME


dados_xavier_max <- expand.grid(
  longitude = lon,
  latitude = lat,
  year = anos_unicos
)

dados_xavier_max$pr_max <- as.vector(pr_max_anual)

# Verificar
head(dados_xavier_max)

dim(dados_xavier_max)

#anos válidos
dados_xavier_tendencia <- dados_xavier_max %>%
  filter(!is.na(pr_max))

#8. Inicio da análise de tedência

# COORDENADAS ÚNICAS DAS CÉLULAS


grade_xavier <- dados_xavier_tendencia %>%
  distinct(longitude, latitude)



# TRANSFORMAR CÉLULAS EM PONTOS


grade_xavier_sf <- st_as_sf(
  grade_xavier,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)



# REGIÕES HIDROGRÁFICAS


regioes_hidrograficas <- st_read(
  "C:/Users/bruno/OneDrive/Área de Trabalho/Bruno/UnB/Tendências/Código_original_Saulo/Teste_resultados/SNIRH_RHI/SNIRH_RegioesHidrograficas.shp"
) %>%
  select(
    RHI_SG,
    RHI_CD,
    RHI_NM,
    geometry
  ) %>%
  st_transform(4326)



# ASSOCIAR CADA CÉLULA À RH

grade_xavier_rh <- st_join(
  grade_xavier_sf,
  regioes_hidrograficas,
  join = st_within
) %>%
  rename(
    regiao_hidrografica = RHI_NM
  )

#confirir
table(
  grade_xavier_rh$regiao_hidrografica,
  useNA = "ifany"
)

#incorporar a RH aos registros:
dados_xavier_tendencia <- dados_xavier_tendencia %>%
  left_join(
    grade_xavier_rh %>%
      st_drop_geometry() %>%
      select(
        longitude,
        latitude,
        regiao_hidrografica
      ),
    by = c("longitude", "latitude")
  )

#confirir
table(
  dados_xavier_tendencia$regiao_hidrografica,
  useNA = "ifany"
)

#denifir RH analisadas
RHs <- dados_xavier_tendencia %>%
  filter(!is.na(regiao_hidrografica)) %>%
  distinct(regiao_hidrografica) %>%
  pull(regiao_hidrografica)

#funções
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
#
# Equação (6)
#
# Para o MTFPW precisamos da tendência estimada para
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
  #
  # Equação (6)
  
  
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
  # Primeiro calculamos r1.
  
  
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
# Implementação seguindo a lógica descrita no artigo:
#
# P(j) < d(j)
#
# d(j) = j * alpha / m
#
# Encontra-se o MAIOR j que satisfaz a condição.
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


#Análise da tendência
analisar_RH_xavier <- function(nome_RH) {
  
  cat("\n============================================\n")
  cat("Analisando:", nome_RH, "\n")
  cat("============================================\n")
  
  dados_RH <- dados_xavier_tendencia %>%
    filter(regiao_hidrografica == nome_RH)
  
  

  # MK + MTFPW POR CÉLULA

  
  resultado_RH <- dados_RH %>%
    arrange(longitude, latitude, year) %>%
    group_by(longitude, latitude) %>%
    summarise(
      
      regiao_hidrografica = first(regiao_hidrografica),
      
      n_anos = sum(is.finite(pr_max)),
      
      # MK original
      mk_original = list(
        mann_kendall_java(pr_max)
      ),
      
      Z_MK = mk_original[[1]]$Z,
      
      p_MK = mk_original[[1]]$pvalue,
      
      sig_MK = mk_original[[1]]$significativo,
      
      dir_MK = case_when(
        Z_MK > 0 ~ "Crescente",
        Z_MK < 0 ~ "Decrescente",
        TRUE ~ "Sem tendência"
      ),
      
      
      # MTFPW
      mtfpw = list(
        mtfpw_java(pr_max)
      ),
      
      beta_MTFPW = mtfpw[[1]]$beta,
      
      r1_MTFPW = mtfpw[[1]]$r1,
      
      r1_corrigido = mtfpw[[1]]$r1_corrigido,
      
      autocorrelacao = mtfpw[[1]]$autocorrelacao,
      
      MTFPW_aplicado = mtfpw[[1]]$aplicada,
      
      
      # MK após MTFPW
      mk_mtfpw = list(
        mann_kendall_java(
          mtfpw[[1]]$serie
        )
      ),
      
      Z_MTFPW = mk_mtfpw[[1]]$Z,
      
      p_MTFPW = mk_mtfpw[[1]]$pvalue,
      
      sig_MTFPW = mk_mtfpw[[1]]$significativo,
      
      dir_MTFPW = case_when(
        Z_MTFPW > 0 ~ "Crescente",
        Z_MTFPW < 0 ~ "Decrescente",
        TRUE ~ "Sem tendência"
      ),
      
      .groups = "drop"
    )
  
  

  # FDR DENTRO DA RH

  
  resultado_RH$sig_FDR <- fdr_java(
    resultado_RH$p_MTFPW,
    alpha = 0.05
  )
  
  

  # RESULTADO FINAL

  
  resultado_RH <- resultado_RH %>%
    mutate(
      
      sig_final = ifelse(
        sig_FDR,
        "S(5)",
        "NS"
      ),
      
      dir_final = case_when(
        sig_FDR & Z_MTFPW > 0 ~ "Crescente",
        sig_FDR & Z_MTFPW < 0 ~ "Decrescente",
        TRUE ~ "Sem tendência"
      )
    )
  
  
  return(resultado_RH)
}

#rodas para as 12 RHs
resultados_RH_xavier <- lapply(
  RHs,
  analisar_RH_xavier
)

#juntar resultados
resultado_xavier <- bind_rows(
  resultados_RH_xavier
) %>%
  arrange(
    regiao_hidrografica,
    longitude,
    latitude
  )

#conferir
dim(resultado_xavier)
head(resultado_xavier)

#resumo por RH
resumo_RH_xavier <- resultado_xavier %>%
  group_by(regiao_hidrografica) %>%
  summarise(
    n_celulas = n(),
    
    MK_significativo = sum(
      sig_MK,
      na.rm = TRUE
    ),
    
    MTFPW_significativo = sum(
      sig_MTFPW,
      na.rm = TRUE
    ),
    
    FDR_significativo = sum(
      sig_FDR,
      na.rm = TRUE
    ),
    
    MTFPW_aplicado = sum(
      MTFPW_aplicado,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(resumo_RH_xavier)

#apos fdr
resultado_xavier %>%
  filter(sig_FDR) %>%
  count(
    regiao_hidrografica,
    dir_final
  )

#salvar resultado
resultado_xavier_csv <- resultado_xavier %>%
  select(
    -mk_original,
    -mtfpw,
    -mk_mtfpw
  )

str(resultado_xavier_csv)

write.csv(
  resultado_xavier_csv,
  "resultado_xavier_MK_MTFPW_FDR2.csv",
  row.names = FALSE,
  na = ""
)

library(writexl)

write_xlsx(
  resultado_xavier_csv,
  "resultado_xavier_MK_MTFPW_FDR2.xlsx"
)

#mapa em grade dos resultados
library(sf)
library(dplyr)
library(ggplot2)

resultado_xavier_sf <- resultado_xavier %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

#projetar o sistema
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

resultado_xavier_sf <- resultado_xavier_sf %>%
  st_transform(5880)

estados_brasil_5880 <- estados_brasil %>%
  st_transform(5880)

#criar classe de tendendia FDR
resultado_xavier_sf <- resultado_xavier_sf %>%
  mutate(
    tendencia = case_when(
      sig_FDR & Z_MTFPW > 0 ~ "Crescente",
      sig_FDR & Z_MTFPW < 0 ~ "Decrescente",
      TRUE ~ "Não significativa"
    )
  )

#mapa
mapa_xavier <- ggplot() +
  
  # Estados
  geom_sf(
    data = estados_brasil_5880,
    fill = "white",
    color = "grey40",
    linewidth = 0.25
  ) +
  
  # Células Xavier
  geom_sf(
    data = resultado_xavier_sf,
    aes(color = tendencia),
    size = 0.7,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "Decrescente" = "red",
      "Não significativa" = "grey75",
      "Crescente" = "blue"
    )
  ) +
  
  labs(
    #title = "Tendência da precipitação anual – Xavier",
    subtitle = "MK + MTFPW+FDR",
    color = "Tendência"
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(.55, "cm"),
    width = unit(.55, "cm"),
    pad_x = unit(.2, "cm"),
    pad_y = unit(.2, "cm"),
    style = north_arrow_minimal(
      line_col = "black",
      text_col = "black"
    )
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = .20,
    height = unit(.25, "cm"),
    text_cex = .7,
    line_width = .7,
    pad_x = unit(.3, "cm"),
    pad_y = unit(.2, "cm")
  ) +
  
  coord_sf(
    expand = FALSE
  ) +
  
  theme_classic()+
  theme(
    text = element_text(size = 10),
    
    axis.title = element_blank(),
    
    legend.position = "bottom",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 10
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    )
  )

mapa_xavier

#salvar
ggsave(
  filename = "mapa_xavier2_FDR.png",
  plot = mapa_xavier,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)


#mapa com mk
resultado_xavier_sf <- resultado_xavier_sf %>%
  mutate(
    tendencia_MK = case_when(
      
      sig_MK & Z_MK > 0 ~ "Crescente",
      
      sig_MK & Z_MK < 0 ~ "Decrescente",
      
      TRUE ~ "Não significativa"
      
    )
  )

#mapa
mapa_xavier_MK <- ggplot() +
  
  # Estados
  geom_sf(
    data = estados_brasil_5880,
    fill = "white",
    color = "grey40",
    linewidth = 0.25
  ) +
  
  # Células Xavier
  geom_sf(
    data = resultado_xavier_sf,
    aes(color = tendencia_MK),
    size = 0.7,
    alpha = 0.8
  ) +
  
  # Cores
  scale_color_manual(
    values = c(
      "Decrescente" = "red",
      "Não significativa" = "grey75",
      "Crescente" = "blue"
    ),
    name = "Tendência (MK)"
  ) +
  
  # Título
  labs(
    #title = "Tendência da precipitação máxima mensal anual – Xavier",
    subtitle = "Mann-Kendall",
    color = "Tendência"
  ) +
  
  # Seta do norte
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(.55, "cm"),
    width = unit(.55, "cm"),
    pad_x = unit(.2, "cm"),
    pad_y = unit(.2, "cm"),
    style = north_arrow_minimal(
      line_col = "black",
      text_col = "black"
    )
  ) +
  
  # Escala
  annotation_scale(
    location = "bl",
    width_hint = .20,
    height = unit(.25, "cm"),
    text_cex = .7,
    line_width = .7,
    pad_x = unit(.3, "cm"),
    pad_y = unit(.2, "cm")
  ) +
  
  # Projeção
  coord_sf(
    expand = FALSE
  ) +
  
  theme_classic() +
  theme(
    text = element_text(size = 10),
    
    axis.title = element_blank(),
    
    legend.position = "bottom",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 10
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    )
  )



mapa_xavier_MK

#salvar
ggsave(
  filename = "mapa_xavier_MK.png",
  plot = mapa_xavier_MK,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

#MAPA MK+MTPW
resultado_xavier_sf <- resultado_xavier_sf %>%
  mutate(
    tendencia_MK_MTFPW = case_when(
      
      sig_MTFPW & Z_MTFPW > 0 ~ "Crescente",
      
      sig_MTFPW & Z_MTFPW < 0 ~ "Decrescente",
      
      TRUE ~ "Não significativa"
      
    )
  )

#MAPA
mapa_xavier_MK_MTFPW <- ggplot() +
  
  # Estados
  geom_sf(
    data = estados_brasil_5880,
    fill = "white",
    color = "grey40",
    linewidth = 0.25
  ) +
  
  # Células Xavier
  geom_sf(
    data = resultado_xavier_sf,
    aes(color = tendencia_MK_MTFPW),
    size = 0.7,
    alpha = 0.8
  ) +
  
  # Cores
  scale_color_manual(
    values = c(
      "Decrescente" = "red",
      "Não significativa" = "grey75",
      "Crescente" = "blue"
    ),
    name = "Tendência"
  ) +
  
  # Título
  labs(
    #title = "Tendência da precipitação máxima mensal anual – Xavier",
    subtitle = "MK + MTFPW",
    color = "Tendência"
  ) +
  
  # Seta do norte
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(.55, "cm"),
    width = unit(.55, "cm"),
    pad_x = unit(.2, "cm"),
    pad_y = unit(.2, "cm"),
    style = north_arrow_minimal(
      line_col = "black",
      text_col = "black"
    )
  ) +
  
  # Escala
  annotation_scale(
    location = "bl",
    width_hint = .20,
    height = unit(.25, "cm"),
    text_cex = .7,
    line_width = .7,
    pad_x = unit(.3, "cm"),
    pad_y = unit(.2, "cm")
  ) +
  
  # Projeção
  coord_sf(
    expand = FALSE
  ) +
  
  theme_classic() +
  theme(
    text = element_text(size = 10),
    
    axis.title = element_blank(),
    
    legend.position = "bottom",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 10
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    )
  )



mapa_xavier_MK_MTFPW

#salvar
ggsave(
  filename = "mapa_xavier_MK_MTFPW.png",
  plot = mapa_xavier_MK_MTFPW,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

mapas_xavier_metodos<-mapa_xavier_MK | mapa_xavier_MK_MTFPW | mapa_xavier 
mapas_xavier_metodos

#salvar
ggsave(
  filename = "mapas_xavier_metodos.png",
  plot = mapas_xavier_metodos,
  width = 14,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

#gráfico significancia
dados_grafico_xavier <- resultado_xavier %>%
  summarise(
    MK_sig = sum(sig_MK, na.rm = TRUE),
    MK_ns = sum(!sig_MK, na.rm = TRUE),
    
    MTFPW_sig = sum(sig_MTFPW, na.rm = TRUE),
    MTFPW_ns = sum(!sig_MTFPW, na.rm = TRUE),
    
    FDR_sig = sum(sig_FDR, na.rm = TRUE),
    FDR_ns = sum(!sig_FDR, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("metodo", "significancia"),
    names_pattern = "(MK|MTFPW|FDR)_(sig|ns)",
    values_to = "n"
  ) %>%
  mutate(
    significancia = case_when(
      significancia == "sig" ~ "Significativa",
      significancia == "ns" ~ "Não significativa"
    ),
    
    metodo = case_when(
      metodo == "MK" ~ "MK",
      metodo == "MTFPW" ~ "MTFPW",
      metodo == "FDR" ~ "MTFPW + FDR"
    )
  )

#grafico
grafico_significancia_xavier <- ggplot(
  dados_grafico_xavier,
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
    breaks = seq(0, max(dados_grafico_xavier$n) + 10000, by = 10000),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = NULL,
    y = "Número de células da grade"
  ) +
  
  theme_minimal() +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.title.y = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 10
    )
  )

grafico_significancia_xavier

#salvar
ggsave(
  filename = "grafico_significancia_xavier.png",
  plot = grafico_significancia_xavier,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

#gráfico tendendia por RH
library(dplyr)
library(ggplot2)


# 1. CLASSIFICAR A TENDÊNCIA


dados_grafico_xavier <- resultado_xavier %>%
  
  mutate(
    
    tendencia = case_when(
      
      sig_FDR & Z_MTFPW < 0 ~ "Decrescente",
      
      sig_FDR & Z_MTFPW > 0 ~ "Crescente",
      
      TRUE ~ "Não significativa"
      
    )
    
  )



# 2. CONTAR AS CÉLULAS POR RH E TENDÊNCIA


dados_grafico_xavier <- dados_grafico_xavier %>%
  
  count(
    regiao_hidrografica,
    tendencia,
    name = "n_celulas"
  )



# 3. CALCULAR A PROPORÇÃO


dados_grafico_xavier <- dados_grafico_xavier %>%
  
  group_by(
    regiao_hidrografica
  ) %>%
  
  mutate(
    
    total_celulas = sum(n_celulas),
    
    proporcao = 100 * n_celulas / total_celulas
    
  ) %>%
  
  ungroup()



# 4. ORDEM DAS CATEGORIAS


dados_grafico_xavier <- dados_grafico_xavier %>%
  
  mutate(
    
    tendencia = factor(
      tendencia,
      levels = c(
        "Decrescente",
        "Não significativa",
        "Crescente"
      )
    )
    
  )



# 5. CRIAR SIGLAS DAS RHs


dados_grafico_xavier <- dados_grafico_xavier %>%
  
  mutate(
    
    sigla_RH = case_when(
      
      regiao_hidrografica == "ATLÂNTICO LESTE" ~ "ALE",
      
      regiao_hidrografica == "AMAZÔNICA" ~ "AMZ",
      
      regiao_hidrografica == "ATLÂNTICO NORDESTE OCIDENTAL" ~ "ANC",
      
      regiao_hidrografica == "ATLÂNTICO NORDESTE ORIENTAL" ~ "ANO",
      
      regiao_hidrografica == "ATLÂNTICO SUDESTE" ~ "ASD",
      
      regiao_hidrografica == "ATLÂNTICO SUL" ~ "ATS",
      
      regiao_hidrografica == "PARNAÍBA" ~ "PNB",
      
      regiao_hidrografica == "PARAGUAI" ~ "PRG",
      
      regiao_hidrografica == "PARANÁ" ~ "PRN",
      
      regiao_hidrografica == "SÃO FRANCISCO" ~ "SFR",
      
      regiao_hidrografica == "TOCANTINS-ARAGUAIA" ~ "TOA",
      
      regiao_hidrografica == "URUGUAI" ~ "URU",
      
      TRUE ~ NA_character_
      
    )
    
  )



# 6. ORDEM DAS RHs


dados_grafico_xavier <- dados_grafico_xavier %>%
  
  mutate(
    
    sigla_RH = factor(
      sigla_RH,
      levels = c(
        "ALE",
        "AMZ",
        "ANC",
        "ANO",
        "ASD",
        "ATS",
        "PNB",
        "PRG",
        "PRN",
        "SFR",
        "TOA",
        "URU"
      )
      
    )
    
  )



# 7. GRÁFICO


prop_sig_xavier <- ggplot(
  
  dados_grafico_xavier,
  
  aes(
    x = sigla_RH,
    y = proporcao,
    fill = tendencia
  )
  
) +
  
  geom_col(
    width = 0.75
  ) +
  
  geom_text(
    
    aes(
      label = n_celulas
    ),
    
    position = position_stack(
      vjust = 0.5
    ),
    
    size = 3.2,
    fontface = "bold"
    
  ) +
  
  scale_fill_manual(
    
    values = c(
      "Decrescente" = "red",
      "Não significativa" = "grey70",
      "Crescente" = "blue"
    )
    
  ) +
  
  scale_y_continuous(
    
    limits = c(0, 100),
    
    breaks = seq(
      0,
      100,
      10
    ),
    
    labels = function(x) paste0(
      x,
      "%"
    ),
    
    expand = c(0, 0)
    
  ) +
  
  labs(
    
    x = "RH",
    
    y = "Proporção de células da grade (%)",
    
    fill = NULL
    
  ) +
  
  theme_classic() +
  
  theme(
    
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 10,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.position = "bottom",
    
    legend.text = element_text(
      size = 10
    )
    
  )


prop_sig_xavier

#salvar
ggsave(
  filename = "prop_sig_xavier.png",
  plot = prop_sig_xavier,
  width = 12,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)



# MAPA DAS REGIÕES HIDROGRÁFICAS DO BRASIL


mapa_RH <- regioes_hidrograficas %>%
  mutate(
    sigla_RH = case_when(
      
      RHI_NM == "ATLÂNTICO LESTE" ~ "ALE",
      RHI_NM == "AMAZÔNICA" ~ "AMZ",
      RHI_NM == "ATLÂNTICO NORDESTE OCIDENTAL" ~ "ANC",
      RHI_NM == "ATLÂNTICO NORDESTE ORIENTAL" ~ "ANO",
      RHI_NM == "ATLÂNTICO SUDESTE" ~ "ASD",
      RHI_NM == "ATLÂNTICO SUL" ~ "ATS",
      RHI_NM == "PARNAÍBA" ~ "PNB",
      RHI_NM == "PARAGUAI" ~ "PRG",
      RHI_NM == "PARANÁ" ~ "PRN",
      RHI_NM == "SÃO FRANCISCO" ~ "SFR",
      RHI_NM == "TOCANTINS-ARAGUAIA" ~ "TOA",
      RHI_NM == "URUGUAI" ~ "URU",
      
      TRUE ~ NA_character_
      
    )
  ) %>%
  st_transform(5880)



# MAPA

grafico_RH <- ggplot() +
  
  geom_sf(
    data = mapa_RH,
    aes(fill = sigla_RH),
    color = "black",
    linewidth = 0.35
  ) +
  
  scale_fill_viridis_d(
    name = "Região Hidrográfica",
    option = "turbo"
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(.55, "cm"),
    width = unit(.55, "cm"),
    pad_x = unit(.2, "cm"),
    pad_y = unit(.2, "cm"),
    style = north_arrow_minimal(
      line_col = "black",
      text_col = "black"
    )
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = .20,
    height = unit(.25, "cm"),
    text_cex = .7,
    line_width = .7,
    pad_x = unit(.3, "cm"),
    pad_y = unit(.2, "cm")
  ) +
  
  coord_sf(
    crs = 5880,
    default_crs = st_crs(4326),
    expand = FALSE,
    datum = st_crs(4326)
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    
    axis.line.x = element_line(
      color = "black",
      linewidth = 0.6
    ),
    
    axis.line.y = element_line(
      color = "black",
      linewidth = 0.6
    ),
    
    axis.ticks.x = element_line(
      color = "black",
      linewidth = 0.4
    ),
    
    axis.ticks.y = element_line(
      color = "black",
      linewidth = 0.4
    ),
    
    axis.ticks.length = unit(0.15, "cm"),
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(size = 10),
    
    legend.position = "right"
  )
grafico_RH


#salvar

ggsave(
  filename = "mapa_RH.png",
  plot = grafico_RH,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

#MAPA CHUVAS XAVIER

# 1. DADOS PARA O MAPA


dados_mapa_xavier <- dados_xavier_max %>%
  filter(!is.na(pr_max))



# 2. ESTADOS DO BRASIL


estados_mapa <- estados_brasil %>%
  st_transform(5880)



# 3. MAPA


mapa_xavier_max <- ggplot() +
  
  # Grade do Xavier
  geom_tile(
    data = dados_mapa_xavier,
    aes(
      x = longitude,
      y = latitude,
      fill = pr_max
    )
  ) +
  
  # Limites dos estados
  geom_sf(
    data = estados_mapa,
    fill = NA,
    color = "black",
    linewidth = 0.25
  ) +
  
  scale_fill_viridis_c(
    name = "Precipitação\n(mm)",
    option = "turbo",
    na.value = "transparent"
  ) +
  
  coord_sf(
    crs = 5880,
    default_crs = st_crs(4326),
    expand = FALSE
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(.55, "cm"),
    width = unit(.55, "cm"),
    pad_x = unit(.2, "cm"),
    pad_y = unit(.2, "cm"),
    style = north_arrow_minimal(
      line_col = "black",
      text_col = "black"
    )
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = .20,
    height = unit(.25, "cm"),
    text_cex = .7,
    line_width = .7,
    pad_x = unit(.3, "cm"),
    pad_y = unit(.2, "cm")
  ) +
  
  labs(
    title = "Precipitação máxima mensal anual — Xavier",
    subtitle = "1961–2025"
  ) +
  
  theme_void() +
  
  theme(
    
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.position = "right"
    
  )

mapa_xavier_max

ggsave(
  filename = "mapa_xavier_max.png",
  plot = mapa_xavier_max,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

