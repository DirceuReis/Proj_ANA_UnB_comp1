# Extração e qualidade das séries diárias Uniplu
# 
# somente dados diários (time_step == 1440)  -> Hidroweb diário e INMET diário
# ignora subdiário (INMET horário, CEMADEN, telemetria)
# Uniplu não traz nível de consistência; usa a série diária publicada
# só postos com mais de 30 anos nos zips, e depois year.size > 30
# (ano bom = no máximo 20% de falhas, 73 dias)

library(dplyr)
library(tidyr)
library(stringr)
library(slider)
library(lubridate)
library(arrow)
library(readr)

# --- caminhos (saídas nesta pasta; Uniplu fora) ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"
DIR_DADOS <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

# Entrada (fora de Scripts_Fluxo_AnoHidrologico)
DIR_ZIP <- file.path(DIR_DADOS, "Dados_Uniplu")

# Saídas (dentro de Scripts_Fluxo_AnoHidrologico)
DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_YEAR     <- file.path(DIR_OUT, "stations", "uniplu_daily_year")
DIR_ANALYZED <- file.path(DIR_OUT, "stations_analyzed")
DIR_INV      <- DIR_DF   # inventário cache junto das demais saídas

TIME_STEP_DAILY <- "1440"
MIN_YEARS_ZIP <- 30
MAX_NA_DAYS <- 73   # 20% de 365 dias
MIN_YEARS_GOOD <- 30
SKIP_EXISTING_YEAR <- TRUE

setwd(DIR_FLUXO)

analyse_station_uniplu <- function(sta) {

  sta <- sta %>%
    mutate(
      codigo = as.character(codigo),
      dt     = as.Date(dt),
      p      = as.numeric(p)
    ) %>%
    filter(!is.na(dt)) %>%
    group_by(codigo, dt) %>%
    summarise(p = mean(p, na.rm = TRUE), .groups = "drop") %>%
    mutate(p = ifelse(is.nan(p), NA_real_, p))

  if (nrow(sta) < 365) {
    return(list(is.null = TRUE, sta = sta))
  }

  dt.min <- min(sta$dt, na.rm = TRUE)
  dt.max <- max(sta$dt, na.rm = TRUE)

  sta <- sta %>%
    complete(dt = seq.Date(dt.min, dt.max, by = "day")) %>%
    fill(codigo, .direction = "downup")

  return_null <- FALSE

  ano.m <- year(dt.max)
  uni.p <- length(unique(sta$p))
  dry.d <- sum(sta$p < 0.5, na.rm = TRUE) / nrow(sta)

  if (uni.p < 15) {
    return_null <- TRUE
  }

  if (is.finite(dry.d) && dry.d >= 0.995) {
    return_null <- TRUE
  }

  sta <- sta %>%
    mutate(
      p3  = slide_dbl(p, sum, .before = 0, .after = 2,  .complete = TRUE, na.rm = TRUE),
      p5  = slide_dbl(p, sum, .before = 0, .after = 4,  .complete = TRUE, na.rm = TRUE),
      p7  = slide_dbl(p, sum, .before = 0, .after = 6,  .complete = TRUE, na.rm = TRUE),
      p15 = slide_dbl(p, sum, .before = 0, .after = 14, .complete = TRUE, na.rm = TRUE),
      p30 = slide_dbl(p, sum, .before = 0, .after = 29, .complete = TRUE, na.rm = TRUE)
    )

  return(list(
    is.null = return_null,
    err     = tibble(codigo = unique(sta$codigo),
                     ano.m  = ano.m,
                     uni.p  = uni.p,
                     dry.d  = dry.d),
    sta     = sta
  ))
}

dir.create(DIR_DF,       recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_YEAR,     recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_ANALYZED, recursive = TRUE, showWarnings = FALSE)

zips <- list.files(DIR_ZIP, pattern = "\\.zip$", full.names = TRUE)
if (!is.null(TEST_UF)) {
  zips <- zips[grepl(paste0("^", TEST_UF, "_"), basename(zips))]
}

n_zips <- length(zips)
message("Zips a ler: ", n_zips)

# ---------------------------------------------------------------------------
# 1) Inventário: table_info de cada zip, só diário
# ---------------------------------------------------------------------------
f_inv <- file.path(DIR_INV, "uniplu_daily_inventory.rds")

if (file.exists(f_inv) && is.null(TEST_UF)) {
  inv <- readRDS(f_inv)
  message("Inventário já existe: ", nrow(inv), " linhas")
} else {
  inv_list <- vector("list", n_zips)
  for (i in seq_along(zips)) {
    zp <- zips[i]
    bn <- tools::file_path_sans_ext(basename(zp))
    parts <- strsplit(bn, "_", fixed = TRUE)[[1]]
    uf_zip <- parts[1]
    year_zip <- as.integer(parts[2])

    tmp <- tempfile()
    dir.create(tmp)
    ok <- tryCatch({
      unzip(zp, files = "table_info.parquet", exdir = tmp, junkpaths = TRUE)
      TRUE
    }, error = function(e) FALSE)

    if (!ok || !file.exists(file.path(tmp, "table_info.parquet"))) {
      unlink(tmp, recursive = TRUE)
      next
    }

    info <- read_parquet(file.path(tmp, "table_info.parquet")) %>%
      mutate(time_step = as.character(time_step)) %>%
      filter(time_step == TIME_STEP_DAILY) %>%
      transmute(
        codigo       = as.character(gauge_code),
        city         = as.character(city),
        estado       = as.character(state),
        lat          = as.numeric(lat),
        long         = as.numeric(long),
        time_step    = time_step,
        network      = as.character(network),
        responsible  = as.character(responsible),
        elevation    = as.numeric(elevation),
        uf_zip       = uf_zip,
        year_zip     = year_zip
      )

    inv_list[[i]] <- info
    unlink(tmp, recursive = TRUE)

    if (i %% 200 == 0) {
      message("inventário ", i, "/", n_zips)
    }
  }

  inv <- bind_rows(inv_list)
  if (is.null(TEST_UF)) {
    saveRDS(inv, f_inv)
  }
  message("Inventário diário: ", nrow(inv), " posto-ano")
}

sta_years <- inv %>%
  group_by(codigo) %>%
  summarise(
    n_years_zip = n_distinct(year_zip),
    year.first_zip = min(year_zip, na.rm = TRUE),
    year.last_zip  = max(year_zip, na.rm = TRUE),
    estado = dplyr::last(estado),
    nome   = dplyr::last(city),
    lat    = dplyr::last(lat),
    long   = dplyr::last(long),
    network = dplyr::last(network),
    responsible = dplyr::last(responsible),
    .groups = "drop"
  )

# mais de 30 anos nos arquivos zip (ainda sem olhar falhas)
keep <- sta_years %>%
  filter(n_years_zip > MIN_YEARS_ZIP)

saveRDS(sta_years, file.path(DIR_DF, "uniplu_daily_sta_years.rds"))
saveRDS(keep,      file.path(DIR_DF, "uniplu_daily_keep30.rds"))
message("Postos diários: ", nrow(sta_years),
        " | com >", MIN_YEARS_ZIP, " anos: ", nrow(keep))

keep_codes <- keep$codigo

# ---------------------------------------------------------------------------
# 2) Extrai table_data dos zips
# ---------------------------------------------------------------------------
for (i in seq_along(zips)) {
  zp <- zips[i]
  bn <- tools::file_path_sans_ext(basename(zp))
  f_out <- file.path(DIR_YEAR, paste0(bn, ".parquet"))

  if (SKIP_EXISTING_YEAR && file.exists(f_out)) {
    next
  }

  parts <- strsplit(bn, "_", fixed = TRUE)[[1]]
  year_zip <- as.integer(parts[2])

  tmp <- tempfile()
  dir.create(tmp)
  ok <- tryCatch({
    unzip(zp, files = c("table_info.parquet", "table_data.parquet"),
          exdir = tmp, junkpaths = TRUE)
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    unlink(tmp, recursive = TRUE)
    next
  }

  f_info <- file.path(tmp, "table_info.parquet")
  f_data <- file.path(tmp, "table_data.parquet")
  if (!file.exists(f_info) || !file.exists(f_data)) {
    unlink(tmp, recursive = TRUE)
    next
  }

  info <- read_parquet(f_info) %>%
    mutate(
      codigo    = as.character(gauge_code),
      time_step = as.character(time_step)
    ) %>%
    filter(time_step == TIME_STEP_DAILY, codigo %in% keep_codes)

  if (nrow(info) == 0) {
    unlink(tmp, recursive = TRUE)
    next
  }

  daily_codes <- unique(info$codigo)

  dat <- read_parquet(f_data) %>%
    transmute(
      codigo = as.character(gauge_code),
      p      = as.numeric(rain_mm),
      dt     = as.Date(datetime)
    ) %>%
    filter(codigo %in% daily_codes, !is.na(dt))

  if (nrow(dat) > 0) {
    write_parquet(dat, f_out)
  }

  unlink(tmp, recursive = TRUE)
  if (i %% 100 == 0) {
    message("extração ", i, "/", n_zips)
  }
}

# ---------------------------------------------------------------------------
# 3) Qualidade por posto
# ---------------------------------------------------------------------------
year_files <- list.files(DIR_YEAR, pattern = "\\.parquet$", full.names = TRUE)
if (!is.null(TEST_UF)) {
  year_files <- year_files[grepl(paste0("^", TEST_UF, "_"), basename(year_files))]
}
if (length(year_files) == 0) {
  stop("Nenhum parquet anual em ", DIR_YEAR,
       if (!is.null(TEST_UF)) paste0(" para TEST_UF=", TEST_UF) else "")
}
message("Parquets anuais: ", length(year_files),
        if (!is.null(TEST_UF)) paste0(" (", TEST_UF, ")") else "")

ds <- open_dataset(year_files)

# processa em blocos
cods <- keep_codes
chunk_n <- 150
n_cod <- length(cods)
n_chunk <- ceiling(n_cod / chunk_n)

df_sta <- tibble(
  codigo     = character(),
  year.first = integer(),
  year.last  = integer(),
  year.size  = integer(),
  estado     = character(),
  nome       = character(),
  lat        = double(),
  long       = double(),
  network    = character()
)

count.all <- 0
count.val <- 0

for (k in seq_len(n_chunk)) {
  i1 <- (k - 1) * chunk_n + 1
  i2 <- min(k * chunk_n, n_cod)
  chunk <- cods[i1:i2]

  raw <- ds %>%
    filter(codigo %in% chunk) %>%
    collect()

  for (cod in chunk) {
    count.all <- count.all + 1
    sta <- raw %>% filter(codigo == cod)

    sta.analysed <- tryCatch(
      analyse_station_uniplu(sta),
      error = function(e) NULL
    )

    if (!is.null(sta.analysed) && !isTRUE(sta.analysed$is.null)) {
      sta.ok <- sta.analysed$sta

      df <- sta.ok %>%
        mutate(ano = year(dt)) %>%
        group_by(ano) %>%
        summarise(n = sum(is.na(p)), .groups = "drop") %>%
        filter(n <= MAX_NA_DAYS)

      year.size <- nrow(df)
      if (year.size > MIN_YEARS_GOOD) {
        count.val <- count.val + 1
        saveRDS(sta.ok, file.path(DIR_ANALYZED, paste0(cod, "_analyzed.rds")))

        meta <- keep %>% filter(codigo == cod)
        df_sta <- df_sta %>%
          add_row(
            codigo     = cod,
            year.first = as.integer(min(df$ano, na.rm = TRUE)),
            year.last  = as.integer(max(df$ano, na.rm = TRUE)),
            year.size  = as.integer(year.size),
            estado     = meta$estado[1],
            nome       = meta$nome[1],
            lat        = meta$lat[1],
            long       = meta$long[1],
            network    = meta$network[1]
          )
      }
    }

    if (count.all %% 50 == 0) {
      message("qualidade ", count.all, "/", n_cod, " válidas: ", count.val)
    }
  }

  rm(raw)
  gc(verbose = FALSE)
}

saveRDS(df_sta, file.path(DIR_DF, "out_analise_uniplu.rds"))
message("Postos válidos (> ", MIN_YEARS_GOOD, " anos bons): ", nrow(df_sta))

