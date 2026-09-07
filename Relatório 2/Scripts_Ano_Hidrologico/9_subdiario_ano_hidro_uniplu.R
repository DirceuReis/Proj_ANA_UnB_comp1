# Ano hidrológico dos postos subdiários via vizinhança Uniplu
#
# Para cada posto subdiário:
#   1) k vizinhos Uniplu (híbrido rbar) mais próximos
#   2) classifica MANCHA vs TRANSIÇÃO (f* moda e R̄ circular dos meses)
#   3) define o mês de início:
#        mancha    → moda dos meses dos vizinhos
#        transição → média circular ponderada (1/dist × rbar)
#
# Entrada:
#   df_hidro_rbar_hibrido.rds  (passo 6)
#   subdiario_br/  OU  df_subdiario_inventario.rds
# Saída:
#   df_subdiario_ano_hidro_uniplu.*
#   pictures/subdiario_uniplu/

library(dplyr)
library(tidyr)
library(readr)
library(sf)
library(ggplot2)
library(geobr)
library(FNN)
library(arrow)
library(patchwork)

# --- caminhos ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"
DIR_DADOS <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro"

OUT_TAG <- "BR"
DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_PIC      <- file.path(DIR_OUT, "pictures", "subdiario_uniplu")
DIR_SUBDIARIO <- file.path(DIR_DADOS, "subdiario_br")
DIR_DF_FALLBACK <- file.path(DIR_DADOS, "resultados", OUT_TAG, "dataframes")

setwd(DIR_FLUXO)
dir.create(DIR_DF, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PIC, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parâmetros
# ---------------------------------------------------------------------------
K_NN          <- 10L
F_MODA_MIN    <- 0.70
RBAR_CIRC_MIN <- 0.75
PESO_RBAR     <- TRUE   # peso = (1/dist) * max(rbar, 0.05)

meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")
meses_lab <- c("Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
               "Jul", "Ago", "Set", "Out", "Nov", "Dez")
cores_mes <- hsv(
  h = ((seq_len(12) - 1L) / 12 + 0.02) %% 1,
  s = 0.82, v = 0.93
)
names(cores_mes) <- as.character(seq_len(12))

mes_para_ang <- function(m) ((as.integer(m) - 1L) / 12) * 2 * pi
ang_para_mes <- function(th) {
  th <- th %% (2 * pi)
  as.integer(((th / (2 * pi)) * 12) %% 12) + 1L
}

media_circ_mes <- function(meses, pesos) {
  ok <- is.finite(meses) & is.finite(pesos) & pesos > 0
  if (!any(ok)) return(NA_integer_)
  m <- as.integer(meses[ok])
  w <- as.numeric(pesos[ok])
  ang <- mes_para_ang(m)
  ang_para_mes(atan2(sum(w * sin(ang)), sum(w * cos(ang))))
}

rbar_meses <- function(meses) {
  m <- as.integer(meses[is.finite(meses)])
  if (length(m) == 0L) return(NA_real_)
  ang <- mes_para_ang(m)
  sqrt(mean(cos(ang))^2 + mean(sin(ang))^2)
}

achar_rds <- function(nome) {
  p1 <- file.path(DIR_DF, nome)
  p2 <- file.path(DIR_DF_FALLBACK, nome)
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  NULL
}

# ---------------------------------------------------------------------------
# Inventário subdiário (coords)
# ---------------------------------------------------------------------------
F_INV <- achar_rds("df_subdiario_inventario.rds")
if (!is.null(F_INV)) {
  message("Inventário: ", F_INV)
  sub <- readRDS(F_INV) %>%
    filter(is.finite(lat), is.finite(long)) %>%
    mutate(gauge_code = as.character(gauge_code))
} else {
  message("Montando inventário a partir de ", DIR_SUBDIARIO, " ...")
  info_files <- list.files(
    DIR_SUBDIARIO, pattern = "_info\\.parquet$",
    recursive = TRUE, full.names = TRUE
  )
  if (length(info_files) == 0L) {
    stop("Sem inventário e sem *_info.parquet em ", DIR_SUBDIARIO)
  }
  inv_list <- vector("list", length(info_files))
  for (i in seq_along(info_files)) {
    f <- info_files[[i]]
    yr <- suppressWarnings(
      as.integer(sub(".*_(\\d{4})_info\\.parquet$", "\\1", basename(f)))
    )
    inf <- tryCatch(read_parquet(f), error = function(e) NULL)
    if (is.null(inf) || nrow(inf) == 0L) next
    inv_list[[i]] <- inf %>%
      transmute(
        gauge_code = as.character(gauge_code),
        nome = as.character(city),
        estado = as.character(state),
        lat = as.numeric(lat),
        long = as.numeric(long),
        time_step = as.integer(as.numeric(time_step)),
        network = as.character(network),
        year_zip = yr
      )
  }
  moda_int <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_integer_)
    ux <- sort(unique(as.integer(x)))
    ux[which.max(tabulate(match(as.integer(x), ux)))]
  }
  sub <- bind_rows(inv_list) %>%
    filter(is.finite(lat), is.finite(long), is.finite(time_step), time_step < 1440L) %>%
    group_by(gauge_code) %>%
    summarise(
      nome = dplyr::last(na.omit(nome)),
      estado = dplyr::last(na.omit(estado)),
      lat = dplyr::last(lat[is.finite(lat)]),
      long = dplyr::last(long[is.finite(long)]),
      time_step = moda_int(time_step),
      network = dplyr::last(na.omit(network)),
      n_years_zip = n_distinct(year_zip),
      .groups = "drop"
    )
  saveRDS(sub, file.path(DIR_DF, "df_subdiario_inventario.rds"))
  write_excel_csv(sub, file.path(DIR_DF, "df_subdiario_inventario.csv"))
}

# ---------------------------------------------------------------------------
# Uniplu (híbrido)
# ---------------------------------------------------------------------------
F_UNI <- achar_rds("df_hidro_rbar_hibrido.rds")
if (is.null(F_UNI)) stop("Falta df_hidro_rbar_hibrido.rds (rode o passo 6)")
message("Uniplu: ", F_UNI)

uni <- readRDS(F_UNI) %>%
  filter(
    is.finite(lat), is.finite(long),
    is.finite(mes_hibrido_rbar),
    mes_hibrido_rbar >= 1, mes_hibrido_rbar <= 12
  ) %>%
  mutate(
    codigo = as.character(codigo),
    mes_hib = as.integer(mes_hibrido_rbar),
    rbar = as.numeric(rbar)
  )

message("Uniplu: ", nrow(uni), " | Subdiário: ", nrow(sub))

# ---------------------------------------------------------------------------
# k-NN + classificação + mês de início
# ---------------------------------------------------------------------------
message("Buscando ", K_NN, " vizinhos Uniplu...")
xy_uni <- as.matrix(uni[, c("long", "lat")])
xy_sub <- as.matrix(sub[, c("long", "lat")])
nn <- FNN::get.knnx(data = xy_uni, query = xy_sub, k = K_NN)
idx <- nn$nn.index
dist_km <- nn$nn.dist * 111

message("Calculando mancha/transição e mês de início...")
n <- nrow(sub)
out_list <- vector("list", n)

for (i in seq_len(n)) {
  ii <- idx[i, ]
  dkm <- pmax(dist_km[i, ], 0.1)
  viz <- uni[ii, ]
  meses_hib <- viz$mes_hib

  tab <- table(meses_hib)
  moda <- as.integer(names(tab)[which.max(tab)])
  f_moda <- as.numeric(max(tab) / length(meses_hib))
  r_loc <- rbar_meses(meses_hib)

  classe <- if (isTRUE(f_moda >= F_MODA_MIN) || isTRUE(r_loc >= RBAR_CIRC_MIN)) {
    "mancha"
  } else {
    "transicao"
  }

  if (classe == "mancha") {
    mes_ini <- moda
    fonte <- "moda_mancha"
  } else {
    w <- 1 / dkm
    if (PESO_RBAR) w <- w * pmax(viz$rbar, 0.05)
    mes_ini <- media_circ_mes(meses_hib, w)
    fonte <- "media_circ_ponderada"
  }

  out_list[[i]] <- tibble(
    gauge_code = sub$gauge_code[i],
    estado = sub$estado[i],
    lat = sub$lat[i],
    long = sub$long[i],
    time_step = if ("time_step" %in% names(sub)) sub$time_step[i] else NA_integer_,
    network = if ("network" %in% names(sub)) sub$network[i] else NA_character_,
    n_viz = length(meses_hib),
    dist_nn_km = round(min(dkm), 2),
    dist_k_med_km = round(median(dkm), 2),
    mes_moda = moda,
    f_moda = round(f_moda, 3),
    rbar_local = round(r_loc, 3),
    classe_espacial = classe,
    mes_inicio = as.integer(mes_ini),
    fonte_inicio = fonte
  )
}

res <- bind_rows(out_list) %>%
  mutate(mes_inicio_nome = meses_pt[mes_inicio])

# ---------------------------------------------------------------------------
# Resumos / saídas
# ---------------------------------------------------------------------------
cat("\n=== Classificação espacial (k=", K_NN, ") ===\n", sep = "")
print(res %>% count(classe_espacial) %>% mutate(pct = round(100 * n / sum(n), 1)))
cat("\n=== Fonte do mês de início ===\n")
print(res %>% count(fonte_inicio) %>% mutate(pct = round(100 * n / sum(n), 1)))

tab_uf <- res %>%
  group_by(estado) %>%
  summarise(
    n = n(),
    pct_transicao = round(100 * mean(classe_espacial == "transicao"), 1),
    f_moda_med = round(median(f_moda), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_transicao))

write_excel_csv(res, file.path(DIR_DF, "df_subdiario_ano_hidro_uniplu.csv"))
saveRDS(res, file.path(DIR_DF, "df_subdiario_ano_hidro_uniplu.rds"))
write_excel_csv(tab_uf, file.path(DIR_DF, "df_subdiario_ano_hidro_uniplu_uf.csv"))

# ---------------------------------------------------------------------------
# Mapas
# ---------------------------------------------------------------------------
ufs <- read_state(year = 2020, showProgress = FALSE) %>% st_transform(4326)
br  <- read_country(year = 2020, showProgress = FALSE) %>% st_transform(4326)

tema <- theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey35"),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
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

g_mes <- ggplot() +
  geom_sf(data = br, fill = "grey96", color = "grey40", linewidth = 0.25) +
  geom_sf(data = ufs, fill = NA, color = "grey60", linewidth = 0.12) +
  geom_point(
    data = res %>% filter(is.finite(mes_inicio)),
    aes(x = long, y = lat, color = factor(mes_inicio, levels = 1:12)),
    size = 1.1, alpha = 0.9, stroke = 0
  ) +
  scale_color_manual(values = cores_mes, drop = FALSE) +
  labs(
    title = "Subdiário — mês de início (vizinhança Uniplu)",
    subtitle = paste0(
      "mancha → moda; transição → média circ. ponderada | k=", K_NN,
      " | n=", nrow(res)
    )
  ) +
  tema

ggsave(
  file.path(DIR_PIC, "fig_subdiario_mes_inicio_uniplu.png"),
  g_mes / legenda_circular() + plot_layout(heights = c(8.2, 1.8)),
  width = 1800, height = 2100, units = "px", dpi = 180, bg = "white"
)

g_cls <- ggplot() +
  geom_sf(data = br, fill = "grey96", color = "grey40", linewidth = 0.25) +
  geom_sf(data = ufs, fill = NA, color = "grey60", linewidth = 0.12) +
  geom_point(
    data = res,
    aes(x = long, y = lat, color = classe_espacial),
    size = 1.0, alpha = 0.9, stroke = 0
  ) +
  scale_color_manual(
    values = c(mancha = "#2166AC", transicao = "#D6604D"),
    name = NULL,
    labels = c(
      mancha = paste0("Mancha (", sum(res$classe_espacial == "mancha"), ")"),
      transicao = paste0("Transição (", sum(res$classe_espacial == "transicao"), ")")
    )
  ) +
  labs(
    title = "Classificação espacial (vizinhança Uniplu)",
    subtitle = paste0(
      "Mancha se f_moda ≥ ", F_MODA_MIN, " ou R̄ ≥ ", RBAR_CIRC_MIN,
      " | k = ", K_NN
    )
  ) +
  tema +
  theme(legend.position = "right")

ggsave(
  file.path(DIR_PIC, "fig_subdiario_mancha_transicao.png"),
  g_cls, width = 1800, height = 1900, units = "px", dpi = 180, bg = "white"
)

message("Tabelas: ", DIR_DF)
message("Figuras: ", DIR_PIC)
message("Concluído — mês de início via Uniplu (mancha/transição).")
