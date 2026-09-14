# Passo 13 — Máximos anuais Uniplu (diário) no ano hidrológico híbrido
#
# Usa a mesma fun_imax_wateryear do passo 12.
# Séries diárias → resolução base = 24 h; durações em HORAS = dias × 24.
#
# Ano hidrológico = mes_hibrido_rbar (passo 6).
# Durações: 1, 2, 3, 4, 5, 6, 7 e 10 dias.
#
# Entrada:
#   df_hidro_rbar_hibrido.rds
#   stations_analyzed/<codigo>_analyzed.rds  (dt, p)
# Saída:
#   dataframes/df_imax_uniplu_hidro.*
#   uniplu_cache/imax_by_uf/<UF>_imax.rds
#
# Uso:
#   Rscript Scripts_Fluxo_AnoHidrologico/13_uniplu_imax_hidro.R

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

# --- caminhos ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"
DIR_DADOS <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro"

OUT_TAG <- "BR"
TEST_UF <- NULL   # NULL = Brasil; ex. "CE"

DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_CACHE    <- file.path(DIR_OUT, "uniplu_cache", "imax_by_uf")
DIR_ANALYZED <- file.path(DIR_OUT, "stations_analyzed")
DIR_DF_FALLBACK <- file.path(DIR_DADOS, "resultados", OUT_TAG, "dataframes")
DIR_ANALYZED_FALLBACK <- file.path(DIR_DADOS, "resultados", OUT_TAG, "stations_analyzed")
DIR_FUNS <- file.path(DIR_FLUXO, "funs")

setwd(DIR_FLUXO)
dir.create(DIR_DF, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_CACHE, recursive = TRUE, showWarnings = FALSE)

source(file.path(DIR_FUNS, "fun_imax_wateryear.R"))

# Dias pedidos → horas (fun_imax_wateryear espera HORAS)
DURATIONS_DAY <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 10L)
DURATIONS_HR  <- DURATIONS_DAY * 24

meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")

achar_rds <- function(nome) {
  p1 <- file.path(DIR_DF, nome)
  p2 <- file.path(DIR_DF_FALLBACK, nome)
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  NULL
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
# Meta: mês de início = híbrido rbar
# ---------------------------------------------------------------------------
F_HIB <- achar_rds("df_hidro_rbar_hibrido.rds")
if (is.null(F_HIB)) stop("Falta df_hidro_rbar_hibrido.rds — rode o passo 6")
message("Híbrido: ", F_HIB)

hib <- readRDS(F_HIB) %>%
  filter(
    is.finite(mes_hibrido_rbar),
    mes_hibrido_rbar >= 1, mes_hibrido_rbar <= 12,
    is.finite(lat), is.finite(long)
  ) %>%
  mutate(
    codigo = as.character(codigo),
    mes_inicio = as.integer(mes_hibrido_rbar),
    mes_inicio_nome = meses_pt[mes_inicio]
  )

ufs <- sort(unique(hib$estado))
if (!is.null(TEST_UF) && nzchar(TEST_UF)) {
  ufs <- intersect(ufs, TEST_UF)
}
message(
  "Postos: ", nrow(hib),
  " | UFs: ", paste(ufs, collapse = ", "),
  " | durações (dias): ", paste(DURATIONS_DAY, collapse = ", "),
  " → horas: ", paste(DURATIONS_HR, collapse = ", ")
)

DIR_STA <- dir_analyzed()
message("Séries: ", DIR_STA)

# ---------------------------------------------------------------------------
# Loop por UF
# ---------------------------------------------------------------------------
all_imax <- list()

for (uf in ufs) {
  f_ckpt <- file.path(DIR_CACHE, paste0(uf, "_imax.rds"))
  if (file.exists(f_ckpt)) {
    message("[", uf, "] checkpoint — carregando")
    all_imax[[uf]] <- readRDS(f_ckpt)
    next
  }

  hib_uf <- hib %>% filter(estado == uf)
  if (nrow(hib_uf) == 0L) next

  message("[", uf, "] lendo ", nrow(hib_uf), " séries ...")
  data_ls <- vector("list", nrow(hib_uf))
  names(data_ls) <- hib_uf$codigo

  for (i in seq_len(nrow(hib_uf))) {
    cod <- hib_uf$codigo[i]
    f_sta <- file.path(DIR_STA, paste0(cod, "_analyzed.rds"))
    if (!file.exists(f_sta)) next
    sta <- tryCatch(readRDS(f_sta), error = function(e) NULL)
    if (is.null(sta) || !"p" %in% names(sta) || !"dt" %in% names(sta)) next

    # fun_imax_wateryear: datetime + rain_mm (diário → POSIXct à meia-noite)
    data_ls[[cod]] <- data.frame(
      datetime = as.POSIXct(as.Date(sta$dt), tz = "UTC"),
      rain_mm = as.numeric(sta$p),
      stringsAsFactors = FALSE
    )
  }

  data_ls <- data_ls[!vapply(data_ls, is.null, logical(1))]
  if (!length(data_ls)) {
    message("[", uf, "] nenhuma série lida")
    next
  }
  message("[", uf, "] séries OK: ", length(data_ls), "/", nrow(hib_uf))

  start_map <- setNames(hib_uf$mes_inicio, hib_uf$codigo)
  start_map <- start_map[names(data_ls)]

  imax_uf <- tryCatch(
    fun_imax_wateryear(
      data = data_ls,
      durations = DURATIONS_HR,
      start_month = start_map,
      which.mon = 1:12,
      names = c("datetime", "rain_mm")
    ),
    error = function(e) {
      message("[", uf, "] ERRO imax: ", conditionMessage(e))
      NULL
    }
  )
  rm(data_ls)
  invisible(gc())

  if (is.null(imax_uf) || nrow(imax_uf) == 0L) {
    message("[", uf, "] nenhum imax gerado")
    next
  }

  imax_uf <- imax_uf %>%
    mutate(
      estado = uf,
      d_days = d / 24,
      time_step_min = 1440L
    ) %>%
    left_join(
      hib_uf %>%
        select(
          codigo, nome, lat, long, network, year.size, rbar,
          mes_inicio, mes_inicio_nome, fonte_hibrido_rbar, classe_rbar
        ),
      by = c("gauge_code" = "codigo")
    )

  saveRDS(imax_uf, f_ckpt)
  all_imax[[uf]] <- imax_uf
  message(
    "[", uf, "] linhas=", nrow(imax_uf),
    " | postos=", n_distinct(imax_uf$gauge_code)
  )
}

# ---------------------------------------------------------------------------
# Consolida
# ---------------------------------------------------------------------------
df_imax <- bind_rows(all_imax)
if (is.null(df_imax) || nrow(df_imax) == 0L) {
  stop("Nenhum imax gerado")
}

tab_d <- df_imax %>%
  filter(is.finite(imax)) %>%
  group_by(d_days, d) %>%
  summarise(
    n = n(),
    n_postos = n_distinct(gauge_code),
    med_mm_h = round(median(imax), 4),
    p95_mm_h = round(quantile(imax, 0.95), 4),
    .groups = "drop"
  )
print(tab_d)

F_OUT_RDS <- file.path(DIR_DF, "df_imax_uniplu_hidro.rds")
F_OUT_CSV <- file.path(DIR_DF, "df_imax_uniplu_hidro.csv")
F_RESUMO  <- file.path(DIR_DF, "df_imax_uniplu_hidro_resumo.csv")

saveRDS(df_imax, F_OUT_RDS)
write_excel_csv(df_imax %>% mutate(date = as.character(date)), F_OUT_CSV)
write_excel_csv(tab_d, F_RESUMO)

message("Total linhas: ", nrow(df_imax))
message("Postos: ", n_distinct(df_imax$gauge_code))
message("UFs: ", n_distinct(df_imax$estado))
message("d (h): ", paste(sort(unique(df_imax$d)), collapse = ", "))
message("d_days: ", paste(sort(unique(df_imax$d_days)), collapse = ", "))
message("start_month = mes_hibrido_rbar | cache: ", DIR_CACHE)
message("Concluído — passo 13 (imax Uniplu via fun_imax_wateryear).")
