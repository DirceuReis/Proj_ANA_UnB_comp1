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
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"
DIR_DADOS <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_PIC      <- file.path(DIR_OUT, "pictures")
DIR_OUT_PIC  <- file.path(DIR_PIC, "ano_hidro_rbar_hibrido")
DIR_ANALYZED <- file.path(DIR_OUT, "stations_analyzed")
DIR_DF_FALLBACK <- file.path(DIR_DADOS, "resultados", OUT_TAG, "dataframes")
DIR_ANALYZED_FALLBACK <- file.path(DIR_DADOS, "resultados", OUT_TAG, "stations_analyzed")

setwd(DIR_FLUXO)
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
  p1 <- file.path(DIR_DF, nome)
  p2 <- file.path(DIR_DF_FALLBACK, nome)
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  stop("Arquivo não encontrado: ", nome)
}

dir_analyzed <- function() {
  if (dir.exists(DIR_ANALYZED) &&
      length(list.files(DIR_ANALYZED, pattern = "_analyzed\\.rds$")) > 0L) {
    return(DIR_ANALYZED)
  }
  if (dir.exists(DIR_ANALYZED_FALLBACK)) return(DIR_ANALYZED_FALLBACK)
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
      best$end <<- as.integer(b)
      best$len <<- as.integer(len)
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
  list(mes = as.integer(vao$end), n_secos = n_secos)
}

# Medianas mensais a partir da série diária
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
      p = sum(p, na.rm = TRUE),
      n = n(),
      n_na = sum(is.na(p)),
      .groups = "drop"
    ) %>%
    mutate(lim_na = pmax(1, round(frac_na * n))) %>%
    filter(n_na <= lim_na)

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
      out[ok_th] <- as.integer(month(
        as.Date("1999-01-01") +
          round(((th[ok_th] + pi) %% (2 * pi)) * 365 / (2 * pi)) - 1
      ))
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

hy <- readRDS(F_HY) %>%
  filter(is.finite(lat), is.finite(long), is.finite(rbar), estado != "BO") %>%
  mutate(codigo = as.character(codigo))

if (!is.null(TEST_UF) && nzchar(TEST_UF)) {
  hy <- hy %>% filter(estado == TEST_UF)
}

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
message("Climatologia OK: ", nrow(clima), " / ", n)

# ---------------------------------------------------------------------------
# Híbrido
# ---------------------------------------------------------------------------
sta0 <- hy %>%
  select(any_of(c(
    "codigo", "nome", "estado", "lat", "long", "network", "year.size",
    "rbar", "theta", "mes_pico", "n_picos", "mes_inicio_oposto"
  ))) %>%
  inner_join(clima, by = "codigo")

# θ+π: usa mes_inicio_oposto se existir; senão calcula de theta
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
      confia_theta, as.integer(mes_oposto), as.integer(mes_relatorio)
    ),
    fonte_hibrido_rbar = if_else(confia_theta, "theta_pi", "ultimo_mes_seco"),
    mes_hibrido_nome = meses_pt[mes_hibrido_rbar],
    mes_relatorio_nome = meses_pt[mes_relatorio],
    mes_oposto_nome = meses_pt[mes_oposto],
    seco_frac = SECO_FRAC,
    is_sul = estado %in% UF_SUL
  ) %>%
  filter(is.finite(mes_hibrido_rbar))

rm(sta0)

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
  " (", round(100 * mean(!sta$confia_theta), 1), "%) → último mês seco (limiar ",
  100 * SECO_FRAC, "%)"
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

message("Brasil:"); print(resumo)
message("Sul:"); print(resumo_sul)

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

write_excel_csv(planilha, file.path(DIR_DF, "df_hidro_rbar_hibrido.csv"))
saveRDS(sta, file.path(DIR_DF, "df_hidro_rbar_hibrido.rds"))

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
write_excel_csv(tab_uf, file.path(DIR_DF, "df_hidro_rbar_hibrido_uf.csv"))

# ---------------------------------------------------------------------------
# Mapas (legenda circular + pontos maiores / figura mais compacta)
# ---------------------------------------------------------------------------
meses_lab <- c("Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
               "Jul", "Ago", "Set", "Out", "Nov", "Dez")

ufs <- read_state(year = 2020, showProgress = FALSE) %>%
  st_transform(4326) %>%
  mutate(estado = abbrev_state)
br <- read_country(year = 2020, showProgress = FALSE) %>% st_transform(4326)
sul_ufs <- ufs %>% filter(estado %in% UF_SUL)

tema_mapa <- theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey35"),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(4, 4, 2, 4)
  )

legenda_circular <- function() {
  d <- tibble(mes = 1:12, lab = meses_lab, xmin = 0:11, xmax = 1:12)
  ggplot(d) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0.58, ymax = 1.0, fill = factor(mes)),
      color = "white", linewidth = 0.35
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = 1.22, label = lab),
      size = 3.2, fontface = "bold", color = "grey15"
    ) +
    scale_fill_manual(values = unname(cores_mes), guide = "none") +
    coord_polar(theta = "x", start = -pi / 2 - pi / 12, direction = 1) +
    ylim(0, 1.45) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 2, 0))
}

mapa_mes <- function(df, col_mes, titulo, sub = NULL, bbox = NULL,
                      pt_size = 1.2) {
  d <- df %>%
    filter(is.finite(.data[[col_mes]]), is.finite(long), is.finite(lat)) %>%
    mutate(mes_f = factor(as.character(as.integer(.data[[col_mes]])),
                          levels = as.character(1:12)))
  g <- ggplot() +
    geom_sf(data = if (is.null(bbox)) br else sul_ufs,
            fill = "grey96", color = "grey35", linewidth = 0.3) +
    geom_sf(data = if (is.null(bbox)) ufs else sul_ufs,
            fill = NA, color = "grey55", linewidth = 0.18) +
    geom_point(
      data = d, aes(x = long, y = lat, color = mes_f),
      size = pt_size, alpha = 0.92, stroke = 0
    ) +
    scale_color_manual(values = cores_mes, drop = FALSE) +
    labs(title = titulo, subtitle = sub) +
    tema_mapa
  if (!is.null(bbox)) {
    g <- g + coord_sf(
      xlim = bbox[c("xmin", "xmax")],
      ylim = bbox[c("ymin", "ymax")],
      expand = FALSE
    )
  }
  g
}

salvar_mapa_legenda <- function(g_mapa, arquivo, w = 1800, h = 2100) {
  g_out <- g_mapa / legenda_circular() + plot_layout(heights = c(8.2, 1.8))
  ggsave(arquivo, g_out, width = w, height = h, units = "px", dpi = 180, bg = "white")
}

sub_txt <- paste0(
  "rbar ≥ ", RBAR_HI, " → θ+π; senão → último mês seco (≤mín+",
  100 * SECO_FRAC, "% amp.) | n = ", nrow(sta)
)

g_leg <- legenda_circular()

salvar_mapa_legenda(
  mapa_mes(sta, "mes_hibrido_rbar",
           "Ano hidrológico — híbrido rbar (Uniplu)", sub_txt, pt_size = 1.25),
  file.path(DIR_OUT_PIC, "fig_brasil_hibrido_rbar.png"),
  w = 1800, h = 2100
)

salvar_mapa_legenda(
  mapa_mes(sta %>% filter(is_sul), "mes_hibrido_rbar",
           "Sul — ano hidrológico híbrido rbar", sub_txt, BB_SUL, pt_size = 2.3),
  file.path(DIR_OUT_PIC, "fig_sul_hibrido_rbar.png"),
  w = 1700, h = 2000
)

tema_fonte <- tema_mapa + theme(legend.position = "right")

g_fonte_br <- ggplot() +
  geom_sf(data = br, fill = "grey96", color = "grey35", linewidth = 0.3) +
  geom_sf(data = ufs, fill = NA, color = "grey55", linewidth = 0.15) +
  geom_point(
    data = sta, aes(x = long, y = lat, color = fonte_hibrido_rbar),
    size = 1.2, alpha = 0.92, stroke = 0
  ) +
  scale_color_manual(
    values = c(theta_pi = "#1b9e77", ultimo_mes_seco = "#d95f02"),
    labels = c(
      theta_pi = paste0("θ+π (rbar ≥ ", RBAR_HI, ")"),
      ultimo_mes_seco = paste0(
        "Último mês seco (rbar < ", RBAR_HI, "; limiar ", 100 * SECO_FRAC, "%)"
      )
    ),
    name = "Fonte"
  ) +
  labs(
    title = "Fonte do híbrido rbar — Brasil",
    subtitle = paste0(
      "θ+π: ", sum(sta$fonte_hibrido_rbar == "theta_pi"),
      "  |  último mês seco: ", sum(sta$fonte_hibrido_rbar == "ultimo_mes_seco")
    )
  ) +
  tema_fonte
ggsave(
  file.path(DIR_OUT_PIC, "fig_brasil_fonte_hibrido.png"),
  g_fonte_br, width = 1800, height = 1900, units = "px", dpi = 180, bg = "white"
)

g_fonte_sul <- ggplot() +
  geom_sf(data = sul_ufs, fill = "grey96", color = "grey40", linewidth = 0.3) +
  geom_point(
    data = sta %>% filter(is_sul),
    aes(x = long, y = lat, color = fonte_hibrido_rbar),
    size = 2.3, alpha = 0.95, stroke = 0
  ) +
  scale_color_manual(
    values = c(theta_pi = "#1b9e77", ultimo_mes_seco = "#d95f02"),
    labels = c(
      theta_pi = paste0("θ+π (rbar ≥ ", RBAR_HI, ")"),
      ultimo_mes_seco = paste0(
        "Último mês seco (rbar < ", RBAR_HI, "; limiar ", 100 * SECO_FRAC, "%)"
      )
    ),
    name = "Fonte"
  ) +
  coord_sf(
    xlim = BB_SUL[c("xmin", "xmax")],
    ylim = BB_SUL[c("ymin", "ymax")],
    expand = FALSE
  ) +
  labs(
    title = paste0("Sul — fonte do híbrido rbar (", RBAR_HI, ")"),
    subtitle = paste0(
      "θ+π: ", sum(sta$is_sul & sta$fonte_hibrido_rbar == "theta_pi"),
      "  |  último mês seco: ",
      sum(sta$is_sul & sta$fonte_hibrido_rbar == "ultimo_mes_seco")
    )
  ) +
  tema_fonte
ggsave(
  file.path(DIR_OUT_PIC, "fig_sul_fonte_hibrido.png"),
  g_fonte_sul, width = 1700, height = 1900, units = "px", dpi = 180, bg = "white"
)

g_rbar <- ggplot() +
  geom_sf(data = br, fill = "grey96", color = "grey35", linewidth = 0.3) +
  geom_sf(data = ufs, fill = NA, color = "grey55", linewidth = 0.15) +
  geom_point(
    data = sta, aes(x = long, y = lat, color = rbar),
    size = 1.15, alpha = 0.92, stroke = 0
  ) +
  scale_color_viridis_c(
    option = "C", limits = c(0, 1), name = "rbar",
    breaks = c(0, 0.2, 0.3, 0.5, 0.75, 1)
  ) +
  labs(
    title = "Concentração dos extremos (rbar)",
    subtitle = paste0("Limiar do híbrido = ", RBAR_HI)
  ) +
  tema_fonte
ggsave(
  file.path(DIR_OUT_PIC, "fig_rbar_brasil.png"),
  g_rbar, width = 1800, height = 1900, units = "px", dpi = 180, bg = "white"
)

g_hist <- sta %>%
  mutate(regiao = if_else(is_sul, "Sul (RS/SC/PR)", "Resto do Brasil")) %>%
  ggplot(aes(x = rbar, fill = regiao)) +
  geom_histogram(bins = 30, alpha = 0.75, position = "identity", color = "white") +
  geom_vline(xintercept = RBAR_HI, linetype = "dashed", color = "grey20") +
  facet_wrap(~regiao, ncol = 1, scales = "free_y") +
  labs(
    title = "Distribuição de rbar",
    subtitle = paste0("Tracejado = limiar híbrido ", RBAR_HI),
    x = "rbar", y = "N postos"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(
  file.path(DIR_OUT_PIC, "fig_hist_rbar.png"),
  g_hist, width = 1400, height = 1700, units = "px", dpi = 150, bg = "white"
)

ggsave(
  file.path(DIR_OUT_PIC, "fig_legenda_meses.png"),
  g_leg, width = 900, height = 900, units = "px", dpi = 150, bg = "white"
)

message("Tabelas: ", DIR_DF)
message("Figuras: ", DIR_OUT_PIC)
message(
  "Concluído — híbrido rbar ≥ ", RBAR_HI,
  " | último mês seco recalculado (SECO_FRAC = ", SECO_FRAC, ")."
)
