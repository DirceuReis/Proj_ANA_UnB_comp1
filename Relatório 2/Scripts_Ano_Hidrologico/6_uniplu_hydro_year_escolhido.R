# Ano hidrológico Uniplu — método HÍBRIDO rbar (critério final adotado)
#
# Regra:
#   se rbar >= 0,30 → θ+π (mês oposto aos extremos)
#   senão           → último mês do maior bloco seco (relatório)
#
# Último mês seco (recalculado aqui):
#   1) mediana mensal da precipitação (anos bons)
#   2) mês seco se mediana ≤ mín + 15% da amplitude anual
#   3) início = último mês do maior bloco seco circular
#
# Entrada:
#   df_hydroyear_uniplu.rds  (theta, rbar, coords)
#   stations_analyzed/<codigo>_analyzed.rds  (séries diárias)
# Saída:
#   dataframes/df_hidro_rbar_hibrido.csv / .rds
#   pictures/ano_hidro_rbar_hibrido/

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(geobr)
library(readr)
library(lubridate)
library(patchwork)

# --- caminhos ---
DIR_FLUXO <- "Relatório 2/Scripts_Ano_Hidrologico"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_PIC      <- file.path(DIR_OUT, "pictures")
DIR_OUT_PIC  <- file.path(DIR_PIC, "ano_hidro_rbar_hibrido")
DIR_ANALYZED <- file.path(DIR_OUT, "stations_analyzed")

dir.create(DIR_DF, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT_PIC, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parâmetros
# ---------------------------------------------------------------------------
RBAR_HI <- 0.30
SECO_FRAC <- 0.15     # limiar: mês seco se ≤ mín + 15% da amplitude
FRAC_NA_MES <- 0.20   # mês anual válido se falhas ≤ 20%

meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")

cores_mes <- hsv(
  h = ((seq_len(12) - 1L) / 12 + 0.02) %% 1,
  s = 0.82, v = 0.93
)

names(cores_mes) <- as.character(seq_len(12))

UF_SUL <- c("RS", "SC", "PR")
BB_SUL <- c(xmin = -58.5, xmax = -47.0, ymin = -34.2, ymax = -21.8)

achar_rds <- function(nome) {
  p <- file.path(DIR_DF, nome)
  if (file.exists(p)) return(p)
  stop("Arquivo não encontrado: ", nome)
}

dir_analyzed <- function() {
  if (dir.exists(DIR_ANALYZED) &&
      length(list.files(DIR_ANALYZED, pattern = "_analyzed\\.rds$")) > 0L) {
    return(DIR_ANALYZED)
  }
  
  stop("Pasta stations_analyzed não encontrada.")
}

# ---------------------------------------------------------------------------
# Funções: climatologia + último mês seco
# ---------------------------------------------------------------------------
maior_vao_circular <- function(quiet) {
  
  n <- length(quiet)
  
  if (!any(quiet)) return(NULL)
  if (all(quiet)) return(list(start = 1L, end = n, len = n))
  
  runs <- list()
  i <- 1L
  
  while (i <= n) {
    
    if (!quiet[i]) {
      i <- i + 1L
      next
    }
    
    j <- i
    while (j <= n && quiet[j]) j <- j + 1L
    
    runs[[length(runs) + 1L]] <- c(i, j - 1L)
    i <- j
  }
  
  best <- list(start = NA_integer_, end = NA_integer_, len = -1L)
  
  add_run <- function(a, b, len) {
    if (len > best$len) {
      best$start <<- as.integer(a)
      best$end   <<- as.integer(b)
      best$len   <<- as.integer(len)
    }
  }
  
  wrap <- quiet[1] && quiet[n] && length(runs) >= 2L
  
  if (wrap) {
    
    last <- runs[[length(runs)]]
    first <- runs[[1]]
    
    add_run(last[1], first[2], (n - last[1] + 1L) + first[2])
    
    if (length(runs) > 2L) {
      runs <- runs[seq(2L, length(runs) - 1L)]
    } else {
      runs <- list()
    }
  }
  
  for (r in runs) {
    add_run(r[1], r[2], r[2] - r[1] + 1L)
  }
  
  best
}

# Último mês do maior bloco seco (critério do relatório)
ultimo_mes_seco <- function(p_med, seco_frac = SECO_FRAC) {
  
  p_med <- as.numeric(p_med)
  
  if (length(p_med) != 12L || !any(is.finite(p_med))) {
    return(list(mes = NA_integer_, n_secos = NA_integer_))
  }
  
  pmin <- min(p_med, na.rm = TRUE)
  pmax <- max(p_med, na.rm = TRUE)
  amp <- pmax - pmin
  
  if (!is.finite(amp) || amp <= 0) {
    m <- as.integer(which.min(p_med))
    return(list(mes = m, n_secos = 12L))
  }
  
  lim <- pmin + seco_frac * amp
  seco <- is.finite(p_med) & (p_med <= lim)
  n_secos <- as.integer(sum(seco))
  
  if (n_secos == 0L) {
    return(list(mes = as.integer(which.min(p_med)), n_secos = 0L))
  }
  
  if (n_secos == 12L) {
    return(list(mes = as.integer(which.min(p_med)), n_secos = 12L))
  }
  
  vao <- maior_vao_circular(seco)
  
  list(
    mes = as.integer(vao$end),
    n_secos = n_secos
  )
}

# ---------------------------------------------------------------------------
# Medianas mensais a partir da série diária
# ---------------------------------------------------------------------------
climatologia_mensal <- function(sta, frac_na = FRAC_NA_MES) {
  
  sta <- sta %>%
    mutate(
      dt = as.Date(dt),
      p = as.numeric(p),
      ano = year(dt),
      mes = month(dt)
    ) %>%
    filter(!is.na(dt))
  
  mensal <- sta %>%
    group_by(ano, mes) %>%
    summarise(
      n_obs = sum(!is.na(p)),
      p = sum(p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n_esperado = days_in_month(make_date(ano, mes, 1)),
      n_na = n_esperado - n_obs,
      frac_falha = n_na / n_esperado
    ) %>%
    filter(frac_falha <= frac_na)
  
  if (nrow(mensal) < 12) return(NULL)
  
  clima <- mensal %>%
    group_by(mes) %>%
    summarise(
      p_mediana = median(p, na.rm = TRUE),
      n_anos_mes = n(),
      .groups = "drop"
    )
  
  if (nrow(clima) < 12) return(NULL)
  
  p_vec <- setNames(rep(NA_real_, 12), sprintf("p_med_%02d", 1:12))
  p_vec[sprintf("p_med_%02d", clima$mes)] <- clima$p_mediana
  
  list(
    p_med = as.numeric(p_vec),
    mes_mais_seco = as.integer(clima$mes[which.min(clima$p_mediana)]),
    n_anos_clima = as.integer(max(clima$n_anos_mes))
  )
}

mes_oposto_de_theta <- function(theta, mes_inicio_oposto = NULL) {
  
  n <- length(theta)
  out <- rep(NA_integer_, n)
  
  if (!is.null(mes_inicio_oposto)) {
    m <- as.integer(mes_inicio_oposto)
    ok <- is.finite(m) & m >= 1L & m <= 12L
    out[ok] <- m[ok]
  }
  
  miss <- !is.finite(out)
  
  if (any(miss)) {
    
    th <- as.numeric(theta)
    ok_th <- miss & is.finite(th)
    
    if (any(ok_th)) {
      out[ok_th] <- as.integer(
        month(
          as.Date("1999-01-01") +
            round(((th[ok_th] + pi) %% (2 * pi)) * 365 / (2 * pi)) - 1
        )
      )
    }
  }
  
  out
}

# ---------------------------------------------------------------------------
# Dados base (hydroyear)
# ---------------------------------------------------------------------------
F_HY <- achar_rds("df_hydroyear_uniplu.rds")
DIR_STA <- dir_analyzed()

message("hydroyear: ", F_HY)
message("séries:    ", DIR_STA)

# ---------------------------------------------------------------------------
# Contagem inicial
# ---------------------------------------------------------------------------
hy0 <- readRDS(F_HY) %>%
  mutate(codigo = as.character(codigo))

n_inicio <- n_distinct(hy0$codigo)

message("Estações no início do df_hydroyear_uniplu.rds: ", n_inicio)

# ---------------------------------------------------------------------------
# Filtros básicos
# ---------------------------------------------------------------------------
hy <- hy0 %>%
  filter(
    is.finite(lat),
    is.finite(long),
    is.finite(rbar),
    estado != "BO"
  )

if (!is.null(TEST_UF) && nzchar(TEST_UF)) {
  hy <- hy %>% filter(estado == TEST_UF)
}

n_apos_base <- n_distinct(hy$codigo)

message(
  "Após lat/long/rbar válidos + estado != BO",
  if (!is.null(TEST_UF)) paste0(" + TEST_UF=", TEST_UF) else "",
  ": ", n_apos_base
)

cods <- hy$codigo
n <- length(cods)

message("Calculando climatologia + último mês seco em ", n, " postos...")

clima_list <- vector("list", n)
ok <- 0L

for (i in seq_len(n)) {
  
  f_sta <- file.path(DIR_STA, paste0(cods[i], "_analyzed.rds"))
  if (!file.exists(f_sta)) next
  
  serie <- tryCatch(readRDS(f_sta), error = function(e) NULL)
  
  if (is.null(serie) || !all(c("dt", "p") %in% names(serie))) next
  
  cl <- tryCatch(
    climatologia_mensal(serie[, c("dt", "p")]),
    error = function(e) NULL
  )
  
  if (is.null(cl)) next
  
  u <- ultimo_mes_seco(cl$p_med, SECO_FRAC)
  
  ok <- ok + 1L
  
  clima_list[[ok]] <- tibble(
    codigo = cods[i],
    mes_relatorio = u$mes,
    n_meses_secos = u$n_secos,
    mes_artigo = cl$mes_mais_seco,
    n_anos_clima = cl$n_anos_clima,
    p_med_01 = cl$p_med[1],  p_med_02 = cl$p_med[2],
    p_med_03 = cl$p_med[3],  p_med_04 = cl$p_med[4],
    p_med_05 = cl$p_med[5],  p_med_06 = cl$p_med[6],
    p_med_07 = cl$p_med[7],  p_med_08 = cl$p_med[8],
    p_med_09 = cl$p_med[9],  p_med_10 = cl$p_med[10],
    p_med_11 = cl$p_med[11], p_med_12 = cl$p_med[12]
  )
  
  if (i %% 200L == 0L) {
    message("  climatologia ", i, "/", n, " | válidas: ", ok)
  }
}

clima <- bind_rows(clima_list)

if (nrow(clima) == 0L) {
  stop("Nenhum posto com climatologia válida.")
}

n_clima <- n_distinct(clima$codigo)

message("Climatologia OK: ", n_clima, " / ", n_apos_base)

# ---------------------------------------------------------------------------
# Híbrido
# ---------------------------------------------------------------------------
sta0 <- hy %>%
  select(any_of(c(
    "codigo", "nome", "estado", "lat", "long", "network", "year.size",
    "rbar", "theta", "mes_pico", "n_picos", "mes_inicio_oposto"
  ))) %>%
  inner_join(clima, by = "codigo")

mio <- if ("mes_inicio_oposto" %in% names(sta0)) {
  sta0$mes_inicio_oposto
} else {
  NULL
}

sta <- sta0 %>%
  mutate(
    mes_oposto = mes_oposto_de_theta(theta, mio),
    mes_relatorio = as.integer(mes_relatorio),
    mes_artigo = as.integer(mes_artigo),
    
    classe_rbar = case_when(
      rbar >= RBAR_HI ~ paste0("alta (>=", RBAR_HI, ")"),
      rbar >= 0.20 ~ paste0("media (0.20-", RBAR_HI, ")"),
      TRUE ~ "baixa (<0.20)"
    ),
    
    confia_theta = is.finite(rbar) & rbar >= RBAR_HI,
    
    mes_hibrido_rbar = if_else(
      confia_theta,
      as.integer(mes_oposto),
      as.integer(mes_relatorio)
    ),
    
    fonte_hibrido_rbar = if_else(
      confia_theta,
      "theta_pi",
      "ultimo_mes_seco"
    ),
    
    mes_hibrido_nome = meses_pt[mes_hibrido_rbar],
    mes_relatorio_nome = meses_pt[mes_relatorio],
    mes_oposto_nome = meses_pt[mes_oposto],
    seco_frac = SECO_FRAC,
    is_sul = estado %in% UF_SUL
  ) %>%
  filter(is.finite(mes_hibrido_rbar))

rm(sta0)

# ---------------------------------------------------------------------------
# Contagem final
# ---------------------------------------------------------------------------
n_final <- n_distinct(sta$codigo)

message("\n==================================================")
message("FLUXO DE ESTAÇÕES")
message("==================================================")
message("Início df_hydroyear:             ", n_inicio)
message("Após filtros básicos:            ", n_apos_base)
message("Com climatologia mensal válida:  ", n_clima)
message("Resultado híbrido final:         ", n_final)
message(
  "Perda total:                     ",
  n_inicio - n_final, " (",
  round(100 * (n_inicio - n_final) / n_inicio, 1),
  "%)"
)
message("==================================================\n")

message(
  "Postos: ", nrow(sta),
  " | Sul: ", sum(sta$is_sul),
  " | rbar≥", RBAR_HI, ": ",
  sum(sta$confia_theta),
  " (", round(100 * mean(sta$confia_theta), 1), "%) → θ+π"
)

message(
  "  rbar<", RBAR_HI, ": ",
  sum(!sta$confia_theta),
  " (", round(100 * mean(!sta$confia_theta), 1),
  "%) → último mês seco (limiar ", 100 * SECO_FRAC, "%)"
)

# ---------------------------------------------------------------------------
# Tabelas
# ---------------------------------------------------------------------------
resumo <- sta %>%
  summarise(
    n = n(),
    rbar_med = median(rbar, na.rm = TRUE),
    pct_rbar_ge_limiar = round(100 * mean(rbar >= RBAR_HI), 1),
    pct_fonte_theta_pi = round(100 * mean(fonte_hibrido_rbar == "theta_pi"), 1),
    pct_fonte_ultimo_seco = round(100 * mean(fonte_hibrido_rbar == "ultimo_mes_seco"), 1)
  )

resumo_sul <- sta %>%
  filter(is_sul) %>%
  summarise(
    n = n(),
    rbar_med = median(rbar, na.rm = TRUE),
    pct_rbar_ge_limiar = round(100 * mean(rbar >= RBAR_HI), 1),
    pct_fonte_theta_pi = round(100 * mean(fonte_hibrido_rbar == "theta_pi"), 1),
    pct_fonte_ultimo_seco = round(100 * mean(fonte_hibrido_rbar == "ultimo_mes_seco"), 1)
  )

message("Brasil:")
print(resumo)

message("Sul:")
print(resumo_sul)

write_excel_csv(
  bind_rows(
    resumo %>% mutate(recorte = "Brasil"),
    resumo_sul %>% mutate(recorte = "Sul_RS_SC_PR")
  ),
  file.path(DIR_DF, "df_hidro_rbar_hibrido_resumo.csv")
)

planilha <- sta %>%
  transmute(
    codigo, nome, estado, lat, long, network, year.size, is_sul,
    rbar, classe_rbar, confia_theta, theta,
    mes_oposto, mes_oposto_nome,
    mes_relatorio, mes_relatorio_nome, n_meses_secos, seco_frac,
    mes_artigo,
    mes_hibrido_rbar, fonte_hibrido_rbar, mes_hibrido_nome
  ) %>%
  arrange(estado, codigo)

# Colunas opcionais (podem faltar no hydroyear do script 3)
for (col in c("n_picos", "mes_pico")) {
  
  if (col %in% names(sta)) {
    planilha[[col]] <- sta[[col]][match(planilha$codigo, sta$codigo)]
  }
}

write_excel_csv(
  planilha,
  file.path(DIR_DF, "df_hidro_rbar_hibrido.csv")
)

saveRDS(
  sta,
  file.path(DIR_DF, "df_hidro_rbar_hibrido.rds")
)

tab_uf <- sta %>%
  group_by(estado) %>%
  summarise(
    n = n(),
    rbar_med = round(median(rbar, na.rm = TRUE), 3),
    pct_rbar_ge_limiar = round(100 * mean(rbar >= RBAR_HI), 1),
    pct_fonte_theta_pi = round(100 * mean(fonte_hibrido_rbar == "theta_pi"), 1),
    .groups = "drop"
  ) %>%
  arrange(estado)

write_excel_csv(
  tab_uf,
  file.path(DIR_DF, "df_hidro_rbar_hibrido_uf.csv")
)

# ---------------------------------------------------------------------------
# Mapas
# ---------------------------------------------------------------------------
library(ggspatial)
library(grid)

font <- "sans"
fontsize <- 8
textsize <- fontsize / .pt

FIG_WIDTH <- 15

meses_lab <- c(
  "Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
  "Jul", "Ago", "Set", "Out", "Nov", "Dez"
)

cores_mes <- c(
  "1"  = "#3B4CC0",
  "2"  = "#6F58C9",
  "3"  = "#9C4DC4",
  "4"  = "#C43A9A",
  "5"  = "#E34A6F",
  "6"  = "#F26B38",
  "7"  = "#F4A62A",
  "8"  = "#D8C63A",
  "9"  = "#8FCB3C",
  "10" = "#35B779",
  "11" = "#1FA4A9",
  "12" = "#2A78B8"
)

ufs <- read_state(year = 2020, showProgress = FALSE) %>%
  st_transform(4326) %>%
  mutate(estado = abbrev_state)

br <- read_country(year = 2020, showProgress = FALSE) %>%
  st_transform(4326)

sul_ufs <- ufs %>%
  filter(estado %in% UF_SUL)

tema_mapa <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    
    text = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.background = element_blank(),
    legend.key = element_blank(),
    
    legend.title = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.text = element_text(
      family = font,
      size = fontsize
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = font,
      size = fontsize
    ),
    
    plot.margin = margin( t = 5, r = 5, b = 5, l = 5)
  )

legenda_circular <- function() {
  
  d <- tibble(mes = 1:12, lab = meses_lab, xmin = 0:11, xmax = 1:12 )
  
  ggplot(d) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0.58, ymax = 1.0, fill = factor(mes)),
      color = "white", linewidth = 0.35
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = 1.22, label = lab),
      size = BASE_SIZE / .pt, fontface = "bold", color = "grey15"
    ) +
    scale_fill_manual(values = unname(cores_mes), guide = "none") +
    coord_polar(theta = "x", start = -pi / 2 - pi / 12, direction = 1) +
    ylim(0, 1.45) +
    theme_void(base_size = BASE_SIZE) +
    theme(plot.margin = margin(0, 0, 2, 0))
}

mapa_mes <- function( df, col_mes, bbox = NULL, pt_size = 1.2 ) {
  
  d <- df %>%
    filter(
      is.finite(.data[[col_mes]]),
      is.finite(long),
      is.finite(lat)
    ) %>%
    mutate(
      mes_f = factor(
        as.character(
          as.integer(
            .data[[col_mes]]
          )
        ),
        levels = as.character(
          1:12
        )
      )
    )
  
  g <- ggplot() +
    
    # -------------------------------------------------------------------------
  # Brasil / Sul
  # -------------------------------------------------------------------------
  
  geom_sf( data = if (is.null(bbox)) br else sul_ufs, fill = "grey98", color = "grey50", linewidth = 0.40) +
    
    # -------------------------------------------------------------------------
  # Estados
  # -------------------------------------------------------------------------
  
  geom_sf( data = if (is.null(bbox)) ufs else sul_ufs, fill = NA, color = "grey65", linewidth = 0.18 ) +
    
    # -------------------------------------------------------------------------
  # Estações
  # -------------------------------------------------------------------------
  
  geom_point(data = d, aes( x = long, y = lat, color = mes_f ),
    size = pt_size, alpha = 0.90, stroke = 0 ) +
    
    # -------------------------------------------------------------------------
  # Cores
  # -------------------------------------------------------------------------
  
  scale_color_manual( name = "Mês", values = cores_mes, breaks = as.character( 1:12 ),
    labels = meses_lab, drop = FALSE ) +
    
    # -------------------------------------------------------------------------
  # Barra de escala
  # -------------------------------------------------------------------------
  
  annotation_scale( location = "br", width_hint = 0.20, height = unit( 1.5, "mm" ),
    text_family = font, text_cex = 0.6, bar_cols = c( "black", "aliceblue" ),
    line_width = 0.4 ) +
    
    # -------------------------------------------------------------------------
  # Rosa dos ventos
  # -------------------------------------------------------------------------
  
  annotation_north_arrow( location = "tr", width = unit( 1, "cm" ),
    height = unit( 1, "cm" ),
    which_north = "true",
    style = north_arrow_nautical( text_family = font, text_size = fontsize,
      fill = c( "black", "aliceblue" ) ) ) +
    
    labs( x = NULL, y = NULL ) +
    
    tema_mapa +
    
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      
      legend.key.width = unit( 12, "pt"),
      legend.key.height = unit( 8, "pt" ),
      legend.key.spacing.x = unit( 2, "pt" ),
      legend.spacing.x = unit( 1, "pt" ) ) +
    
    guides( color = guide_legend( nrow = 2, byrow = TRUE,
        title.position = "top", override.aes = list( size = 2.5, alpha = 1 )
      )
    )
  
  if (!is.null(bbox)) {
    
    g <- g + coord_sf( xlim = bbox[ c( "xmin", "xmax" ) ],
        ylim = bbox[ c( "ymin", "ymax" ) ], expand = FALSE )
  }
  g
}

g_hibrido_br <- mapa_mes( sta, "mes_hibrido_rbar", pt_size = 1.25 )

ggsave( file.path( DIR_OUT_PIC, "fig_brasil_hibrido_rbar.png" ),
  g_hibrido_br, width = FIG_WIDTH, height = 12, units = "cm", dpi = 600,
  bg = "white"
)


# =============================================================================
# MAPA DO SUL — ANO HIDROLÓGICO HÍBRIDO
# =============================================================================

g_hibrido_sul <- mapa_mes(
  sta %>%
    filter(
      is_sul
    ),
  "mes_hibrido_rbar",
  bbox = BB_SUL,
  pt_size = 2.3
)

ggsave(
  file.path(
    DIR_OUT_PIC,
    "fig_sul_hibrido_rbar.png"
  ),
  g_hibrido_sul,
  width = FIG_WIDTH,
  height = 11,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# =============================================================================
# MAPA DA FONTE DO MÉTODO — BRASIL
# =============================================================================

cores_fonte <- c(
  "theta_pi" = "#0072B2",
  "ultimo_mes_seco" = "#D55E00"
)

labels_fonte <- c(
  "theta_pi" = expression(theta + pi),
  "ultimo_mes_seco" = "Último mês seco"
)

g_fonte_br <- ggplot() +
  
  geom_sf(
    data = br,
    fill = "grey98",
    color = "grey50",
    linewidth = 0.40
  ) +
  
  geom_sf(
    data = ufs,
    fill = NA,
    color = "grey65",
    linewidth = 0.18
  ) +
  
  geom_point(
    data = sta,
    aes(
      x = long,
      y = lat,
      color = fonte_hibrido_rbar
    ),
    size = 1.2,
    alpha = 0.90,
    stroke = 0
  ) +
  
  scale_color_manual(
    name = "Fonte",
    values = cores_fonte,
    labels = labels_fonte
  ) +
  
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit(
      1.5,
      "mm"
    ),
    text_family = font,
    text_cex = 0.6,
    bar_cols = c(
      "black",
      "aliceblue"
    ),
    line_width = 0.4
  ) +
  
  annotation_north_arrow(
    location = "tr",
    width = unit(
      1,
      "cm"
    ),
    height = unit(
      1,
      "cm"
    ),
    which_north = "true",
    style = north_arrow_nautical(
      text_family = font,
      text_size = fontsize,
      fill = c(
        "black",
        "aliceblue"
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  tema_mapa +
  
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      title.position = "left",
      override.aes = list(
        size = 2.5,
        alpha = 1
      )
    )
  )

ggsave(
  file.path(
    DIR_OUT_PIC,
    "fig_brasil_fonte_hibrido.png"
  ),
  g_fonte_br,
  width = FIG_WIDTH,
  height = 10,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# =============================================================================
# MAPA DA FONTE — SUL
# =============================================================================

g_fonte_sul <- ggplot() +
  
  geom_sf(
    data = sul_ufs,
    fill = "grey98",
    color = "grey50",
    linewidth = 0.40
  ) +
  
  geom_point(
    data = sta %>%
      filter(
        is_sul
      ),
    aes(
      x = long,
      y = lat,
      color = fonte_hibrido_rbar
    ),
    size = 2.3,
    alpha = 0.90,
    stroke = 0
  ) +
  
  scale_color_manual(
    name = "Fonte",
    values = cores_fonte,
    labels = labels_fonte
  ) +
  
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit(
      1.5,
      "mm"
    ),
    text_family = font,
    text_cex = 0.6,
    bar_cols = c(
      "black",
      "aliceblue"
    ),
    line_width = 0.4
  ) +
  
  annotation_north_arrow(
    location = "tr",
    width = unit(
      1,
      "cm"
    ),
    height = unit(
      1,
      "cm"
    ),
    which_north = "true",
    style = north_arrow_nautical(
      text_family = font,
      text_size = fontsize,
      fill = c(
        "black",
        "aliceblue"
      )
    )
  ) +
  
  coord_sf(
    xlim = BB_SUL[
      c(
        "xmin",
        "xmax"
      )
    ],
    ylim = BB_SUL[
      c(
        "ymin",
        "ymax"
      )
    ],
    expand = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  tema_mapa +
  
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    DIR_OUT_PIC,
    "fig_sul_fonte_hibrido.png"
  ),
  g_fonte_sul,
  width = FIG_WIDTH,
  height = 10,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# =============================================================================
# MAPA DE RBAR
# =============================================================================

g_rbar <- ggplot() +
  
  geom_sf(
    data = br,
    fill = "grey98",
    color = "grey50",
    linewidth = 0.40
  ) +
  
  geom_sf(
    data = ufs,
    fill = NA,
    color = "grey65",
    linewidth = 0.18
  ) +
  
  geom_point(
    data = sta,
    aes(
      x = long,
      y = lat,
      color = rbar
    ),
    size = 1.15,
    alpha = 0.90,
    stroke = 0
  ) +
  
  scale_color_viridis_c(
    option = "C",
    limits = c(
      0,
      1
    ),
    name = expression(bar(r)),
    breaks = c(
      0,
      0.2,
      0.3,
      0.5,
      0.75,
      1
    )
  ) +
  
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit(
      1.5,
      "mm"
    ),
    text_family = font,
    text_cex = 0.6,
    bar_cols = c(
      "black",
      "aliceblue"
    ),
    line_width = 0.4
  ) +
  
  annotation_north_arrow(
    location = "tr",
    width = unit(
      1,
      "cm"
    ),
    height = unit(
      1,
      "cm"
    ),
    which_north = "true",
    style = north_arrow_nautical(
      text_family = font,
      text_size = fontsize,
      fill = c(
        "black",
        "aliceblue"
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  tema_mapa +
  
  theme(
    legend.position = "bottom",
    legend.key.width = unit(
      28,
      "pt"
    ),
    legend.key.height = unit(
      5,
      "pt"
    )
  )

ggsave(
  file.path(
    DIR_OUT_PIC,
    "fig_rbar_brasil.png"
  ),
  g_rbar,
  width = FIG_WIDTH,
  height = 8,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# =============================================================================
# HISTOGRAMA DE RBAR
# =============================================================================

g_hist <- sta %>%
  mutate(
    regiao = if_else(
      is_sul,
      "Sul (RS/SC/PR)",
      "Resto do Brasil"
    )
  ) %>%
  ggplot(
    aes(
      x = rbar
    )
  ) +
  
  geom_histogram(
    bins = 30,
    fill = "grey75",
    color = "grey30",
    linewidth = 0.25
  ) +
  
  geom_vline(
    xintercept = RBAR_HI,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  facet_wrap(
    ~regiao,
    ncol = 1,
    scales = "free_y"
  ) +
  
  labs(
    x = expression(bar(r)),
    y = "Número de estações"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      linewidth = 0.3
    ),
    
    text = element_text(
      family = font,
      size = fontsize
    ),
    
    axis.text = element_text(
      family = font,
      size = fontsize
    ),
    
    axis.title = element_text(
      family = font,
      size = fontsize
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = font,
      size = fontsize
    ),
    
    panel.border = element_rect(
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

ggsave(
  file.path(
    DIR_OUT_PIC,
    "fig_hist_rbar.png"
  ),
  g_hist,
  width = FIG_WIDTH,
  height = 10,
  units = "cm",
  dpi = 600,
  bg = "white"
)

message("Tabelas: ", DIR_DF)
message("Figuras: ", DIR_OUT_PIC)

message(
  "Concluído — híbrido rbar ≥ ", RBAR_HI,
  " | último mês seco recalculado (SECO_FRAC = ",
  SECO_FRAC, ")."
)