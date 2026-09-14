# Passo 12 — Intensidades máximas anuais (imax) no ano hidrológico subdiário
#
# Ano hidrológico = mes_inicio do passo 9 (vizinhança Uniplu):
#   mancha    → moda dos vizinhos
#   transição → média circular ponderada
#
# Pipeline por UF:
#   parquet subdiario_br → fill (fun_filter_set) → agrupa por time_step
#   → fun_imax_wateryear(start_month = mes_inicio)
#
# Entrada:
#   df_subdiario_ano_hidro_uniplu.rds   (passo 9)
#   df_subdiario_inventario.rds 
#   subdiario_br/<UF>/*_data.parquet
# Saída:
#   dataframes/df_imax_subdiario.*
#   subdiario_cache/imax_by_uf/<UF>_imax.rds
#
# Uso:
#   Rscript Scripts_Fluxo_AnoHidrologico/12_subdiario_imax.R
#

library(dplyr)
library(tidyr)
library(readr)
library(arrow)
library(lubridate)

# --- caminhos ---
DIR_FLUXO <- "Relatório 2/Scripts_Maximas"
DIR_DADOS <- "Relatório 2/Scripts_Ano_Hidrologico"

TEST_UF <- NULL   # NULL = Brasil; ex.: "CE", "RS"

OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) {
  "BR"
} else {
  TEST_UF
}

DIR_ANO_HIDRO <- file.path("Relatório 2/Scripts_Ano_Hidrologico")
DIR_OUT_HIDRO <- file.path(DIR_ANO_HIDRO,"resultados",OUT_TAG)
DIR_DF_HIDRO <- file.path(DIR_OUT_HIDRO, "dataframes")

DIR_OUT       <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF        <- file.path(DIR_OUT, "dataframes")
DIR_SUB_CACHE <- file.path(DIR_OUT, "subdiario_cache")
DIR_SUBDIARIO <- file.path(file.path("base/fonte/consolidado/subdiario_br"))
DIR_FUNS <- file.path(DIR_FLUXO, "funs")

dir.create(DIR_DF, recursive = TRUE, showWarnings = FALSE)

source(file.path(DIR_FUNS, "fun_filter_set.R"))
source(file.path(DIR_FUNS, "fun_group_ts.R"))
source(file.path(DIR_FUNS, "fun_imax_wateryear.R"))

# Durações (horas): sub-hora + 1–23 h
DURATIONS_HR <- c(
  c(10, 15, 20, 30, 40, 45, 50) / 60,
  1:24
)

meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")

achar_rds <- function(nome) {
  p <- file.path(DIR_DF_HIDRO, nome)
  if (file.exists(p)) return(p)
  NULL
}

# ---------------------------------------------------------------------------
# mes_inicio (passo 9) + metadados do inventário
# ---------------------------------------------------------------------------
F_INI <- achar_rds("df_subdiario_ano_hidro_uniplu.rds")
if (is.null(F_INI)) {
  stop("Falta df_subdiario_ano_hidro_uniplu.rds — rode o passo 9")
}
message("Mês de início: ", F_INI)
ini <- readRDS(F_INI)

F_INV <- achar_rds("df_subdiario_inventario.rds")
if (is.null(F_INV)) {
  # inventário legado com nomes próximos
  F_INV <- achar_rds("df_subdiario_mes_oposto_nearest.rds")
}
if (is.null(F_INV)) stop("Falta inventário subdiário (responsible/time_step)")
message("Inventário: ", F_INV)
inv <- readRDS(F_INV)

# Colunas mínimas do inventário
inv_cols <- intersect(
  c("gauge_code", "nome", "estado", "lat", "long", "time_step",
    "elevation", "network", "responsible", "n_years_zip",
    "year_first", "year_last"),
  names(inv)
)
inv <- inv %>%
  select(all_of(inv_cols)) %>%
  mutate(gauge_code = as.character(gauge_code))

inv_join <- inv %>%
  select(-any_of(c("estado", "lat", "long", "network"))) %>%
  rename(time_step_inv = time_step)

rastreio <- ini %>%
  mutate(gauge_code = as.character(gauge_code)) %>%
  filter(is.finite(mes_inicio), mes_inicio >= 1, mes_inicio <= 12) %>%
  select(
    gauge_code, estado, lat, long, time_step, network,
    classe_espacial, mes_inicio, fonte_inicio,
    mes_moda, f_moda, rbar_local, dist_nn_km
  ) %>%
  left_join(inv_join, by = "gauge_code") %>%
  mutate(
    mes_inicio = as.integer(mes_inicio),
    mes_inicio_nome = meses_pt[mes_inicio],
    metodo_ano_hidro = "uniplu_vizinhanca",
    time_step = as.integer(dplyr::coalesce(time_step, time_step_inv)),
    responsible = if ("responsible" %in% names(.)) {
      as.character(responsible)
    } else {
      "NA"
    },
    responsible = ifelse(is.na(responsible) | !nzchar(responsible),
                         "NA", responsible)
  ) %>%
  select(-any_of("time_step_inv"))

message(
  "Postos com mes_inicio: ", nrow(rastreio),
  " | mancha=", sum(rastreio$classe_espacial == "mancha", na.rm = TRUE),
  " | transição=", sum(rastreio$classe_espacial == "transicao", na.rm = TRUE)
)

saveRDS(rastreio, file.path(DIR_DF, "df_subdiario_mes_inicio_imax.rds"))
write_excel_csv(rastreio, file.path(DIR_DF, "df_subdiario_mes_inicio_imax.csv"))

ufs <- sort(unique(rastreio$estado))
if (!is.null(TEST_UF) && nzchar(TEST_UF)) {
  ufs <- intersect(ufs, TEST_UF)
}
message("UFs: ", paste(ufs, collapse = ", "))

F_CKPT <- file.path(DIR_SUB_CACHE, "imax_by_uf")
dir.create(F_CKPT, recursive = TRUE, showWarnings = FALSE)

F_OUT_RDS <- file.path(DIR_DF, "df_imax_subdiario.rds")
F_OUT_PQT <- file.path(DIR_DF, "df_imax_subdiario.parquet")
F_OUT_CSV <- file.path(DIR_DF, "df_imax_subdiario.csv")

# ---------------------------------------------------------------------------
# Loop por UF
# ---------------------------------------------------------------------------
all_imax <- list()

for (uf in ufs) {
  f_ckpt <- file.path(F_CKPT, paste0(uf, "_imax.rds"))
  if (file.exists(f_ckpt)) {
    message("[", uf, "] checkpoint — carregando")
    all_imax[[uf]] <- readRDS(f_ckpt)
    next
  }

  rast_uf <- rastreio %>% filter(estado == uf)
  if (nrow(rast_uf) == 0L) next

  data_files <- list.files(
    file.path(DIR_SUBDIARIO, uf),
    pattern = "_data\\.parquet$",
    full.names = TRUE
  )
  if (length(data_files) == 0L) {
    message("[", uf, "] sem parquet em subdiario_br/", uf)
    next
  }

  message("[", uf, "] lendo ", length(data_files),
          " arquivos | postos=", nrow(rast_uf))
  chunks <- vector("list", length(data_files))
  for (j in seq_along(data_files)) {
    d <- tryCatch(read_parquet(data_files[[j]]), error = function(e) NULL)
    if (is.null(d)) next
    chunks[[j]] <- d %>%
      transmute(
        gauge_code = as.character(gauge_code),
        rain_mm = as.numeric(rain_mm),
        datetime = as.POSIXct(datetime, tz = "UTC")
      )
  }
  raw <- bind_rows(chunks)
  rm(chunks)
  invisible(gc())

  raw <- raw %>% filter(gauge_code %in% rast_uf$gauge_code)
  if (nrow(raw) == 0L) {
    message("[", uf, "] sem cruzamento posto×dados")
    next
  }

  meta <- rast_uf %>%
    select(gauge_code, time_step, responsible, mes_inicio)
  raw <- raw %>%
    left_join(meta, by = "gauge_code") %>%
    filter(is.finite(time_step), time_step >= 1L)

  message("[", uf, "] preenchendo grades temporais (fun_filter_set)...")
  filled <- tryCatch(
    fun_filter_set(
      data = raw,
      daily = FALSE,
      col_names = c("gauge_code", "rain_mm", "datetime", "time_step", "responsible"),
      filter = FALSE
    ),
    error = function(e) {
      message("[", uf, "] ERRO fill: ", conditionMessage(e))
      NULL
    }
  )
  rm(raw)
  invisible(gc())

  if (is.null(filled) || nrow(filled) == 0L) {
    message("[", uf, "] nenhuma série preenchida")
    next
  }

  data.ls <- split(filled, filled$gauge_code)
  rm(filled)
  invisible(gc())

  data.by.time <- fun_group_ts(data.ls, ts_name = "time_step")
  time.steps <- names(data.by.time)
  message("[", uf, "] resoluções: ", paste(time.steps, collapse = ", "), " min")

  start_map <- setNames(as.integer(rast_uf$mes_inicio), rast_uf$gauge_code)

  imax.ls <- lapply(time.steps, function(ts) {
    current.ts <- data.by.time[[ts]]
    ds <- as.numeric(ts) / 60
    valid.durations <- DURATIONS_HR[(DURATIONS_HR / ds) %% 1 < 1e-6]

    if (length(valid.durations) == 0L) {
      message("[", uf, "] ts=", ts, " min — sem durações compatíveis")
      return(NULL)
    }

    message(
      "[", uf, "] ts=", ts, " min | durações=", length(valid.durations),
      " | postos=", length(current.ts)
    )

    current.imax <- lapply(current.ts, function(df) df[, c("datetime", "rain_mm")])
    names(current.imax) <- names(current.ts)
    sm <- start_map[names(current.imax)]

    out <- tryCatch(
      fun_imax_wateryear(
        data = current.imax,
        durations = valid.durations,
        start_month = sm,
        which.mon = 1:12,
        names = c("datetime", "rain_mm")
      ),
      error = function(e) {
        message("[", uf, "] ts=", ts, " ERRO imax: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(out) || nrow(out) == 0L) return(NULL)
    out$time_step_min <- as.integer(ts)
    out
  })

  imax_uf <- bind_rows(imax.ls[!vapply(imax.ls, is.null, logical(1))])
  rm(data.ls, data.by.time, imax.ls)
  invisible(gc())

  if (nrow(imax_uf) > 0L) {
    imax_uf$estado <- uf
    saveRDS(imax_uf, f_ckpt)
    all_imax[[uf]] <- imax_uf
    message(
      "[", uf, "] linhas=", nrow(imax_uf),
      " | postos=", n_distinct(imax_uf$gauge_code)
    )
  } else {
    message("[", uf, "] nenhum imax gerado")
  }
}

# ---------------------------------------------------------------------------
# Consolida
# ---------------------------------------------------------------------------
df_imax <- bind_rows(all_imax)
if (is.null(df_imax) || nrow(df_imax) == 0L) {
  stop("Nenhum imax gerado")
}

df_imax <- df_imax %>%
  left_join(
    rastreio %>%
      select(
        gauge_code, lat, long, time_step, mes_inicio, mes_inicio_nome,
        network, classe_espacial, fonte_inicio, metodo_ano_hidro,
        f_moda, rbar_local, dist_nn_km
      ),
    by = "gauge_code",
    suffix = c("", "_meta")
  ) %>%
  # start_month da função deve coincidir com mes_inicio
  mutate(
    mes_inicio = dplyr::coalesce(
      as.integer(mes_inicio),
      as.integer(start_month)
    )
  )

saveRDS(df_imax, F_OUT_RDS)
write_parquet(df_imax, F_OUT_PQT)
write_excel_csv(df_imax %>% mutate(date = as.character(date)), F_OUT_CSV)

message("Total linhas: ", nrow(df_imax))
message("Postos: ", n_distinct(df_imax$gauge_code))
message("UFs: ", n_distinct(df_imax$estado))
message("Resoluções (min): ", paste(sort(unique(df_imax$time_step_min)), collapse = ", "))
message("start_month = mes_inicio (passo 9) | cache: ", F_CKPT)
message("Concluído — passo 12 (imax subdiário).")
