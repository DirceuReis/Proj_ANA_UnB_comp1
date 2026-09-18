## Extração, Caracterização e Séries Diárias (convencionais e automáticas) - UNIPLU-BR
##

##### instalar pacotes (rodar só uma vez) #####
#install.packages(c("dplyr","readr","ggplot2","arrow","lubridate","tidyr","stringr"))
# o script de figuras usa, além destes: sf, geobr, ggspatial, scales

##### 1. Library imports -----------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(lubridate)
library(tidyr)
library(stringr)

##### 2. Caminhos --------------------------------------------------------------------
# Pasta onde estão os arquivos .zip do UNIPLU-BR
base_path <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/Dados/UNIPLU_BR"

# Pasta de saída
out_path <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/00. Series"
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

##### 3. Função de leitura dos ZIPs -----------------------------------------------------------
read_UNIPLU_BR <- function(zip_path, table = "table_data") {
  parquet_file <- paste0(table, ".parquet")

  if (!file.exists(zip_path)) {
    message("Erro: arquivo não encontrado: ", zip_path)
    return(data.frame())
  }

  zip_contents <- tryCatch(
    unzip(zip_path, list = TRUE),
    error = function(e) {
      message("Erro ao abrir o ZIP: ", zip_path)
      return(NULL)
    }
  )

  if (is.null(zip_contents)) return(data.frame())
  if (!(parquet_file %in% zip_contents$Name)) {
    message("Erro: ", parquet_file, " não encontrado dentro do ZIP.")
    return(data.frame())
  }

  temp_dir <- tempfile("uniplu_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  unzip(zip_path, files = parquet_file, exdir = temp_dir)
  parquet_path <- file.path(temp_dir, parquet_file)

  df <- tryCatch(
    as.data.frame(arrow::read_parquet(parquet_path)),
    error = function(e) {
      message("Erro ao ler ", parquet_file, ": ", e$message)
      return(data.frame())
    }
  )

  return(df)
}

##### 4. Detectar arquivos disponíveis na pasta -----------------------------------------------------------
# Detecta automaticamente todos os arquivos "UF_ANO.zip" presentes na pasta.
# Ajuste os filtros abaixo caso queira restringir o processamento.
zip_files <- list.files(base_path, pattern = "^[A-Za-z]{2}_[0-9]{4}\\.zip$")

file_index <- tibble(
  file  = zip_files,
  state = str_extract(zip_files, "^[A-Za-z]{2}"),
  year  = as.integer(str_extract(zip_files, "[0-9]{4}"))
)

states_filter <- NULL   # ex.: c("AC", "RS")  |  NULL = todos os estados disponíveis
years_filter  <- NULL   # ex.: 2000:2020      |  NULL = todos os anos disponíveis

if (!is.null(states_filter)) file_index <- filter(file_index, state %in% states_filter)
if (!is.null(years_filter))  file_index <- filter(file_index, year  %in% years_filter)

file_index <- arrange(file_index, state, year)

cat("Arquivos encontrados na pasta:", length(zip_files), "\n")
cat("Arquivos selecionados para processamento:", nrow(file_index), "\n")
cat("Estados:", paste(unique(file_index$state), collapse = ", "), "\n")

##### 5. Conversão para série diária -----------------------------------------------------------
# Regra de agregação diária:
#  - Estações com time_step >= 1440 min (já diárias): mantém a data do próprio registro.
#  - Estações sub-diárias (time_step < 1440 min): aplica a convenção do "dia
#    hidrológico" 9h-9h (corte às 12h UTC), replicando a lógica usada no
#    exemplo da rede CEMADEN em "00. Leitura Uniplu.R", generalizada para
#    qualquer rede sub-diária.
LIMIAR_DIA_CHUVOSO <- 1  # mm; limiar (WMO) para contar um "dia chuvoso"

to_daily_generic <- function(df_data, df_info) {

  df <- df_data %>%
    left_join(df_info %>% select(gauge_code, time_step), by = "gauge_code") %>%
    mutate(time_step = suppressWarnings(as.numeric(time_step)))

  # --- estações já diárias ---
  df_daily_native <- df %>%
    filter(is.na(time_step) | time_step >= 1440) %>%
    mutate(date = as.Date(datetime)) %>%
    group_by(gauge_code, date) %>%
    summarise(rain_mm = sum(rain_mm, na.rm = TRUE), .groups = "drop")

  # --- estações sub-diárias: dia hidrológico 9h-9h ---
  df_sub <- df %>% filter(!is.na(time_step) & time_step < 1440)

  if (nrow(df_sub) > 0) {
    offset_sec <- 12 * 3600  # 12h UTC ~ 9h local (UTC-3)
    df_sub <- df_sub %>%
      mutate(
        datetime = as.POSIXct(datetime, tz = "UTC"),
        day_cut = as.POSIXct(
          floor((as.numeric(datetime) - offset_sec - 1) / 86400) * 86400 +
            offset_sec + 86400,
          origin = "1970-01-01", tz = "UTC"
        )
      ) %>%
      group_by(gauge_code, day_cut) %>%
      summarise(rain_mm = sum(rain_mm, na.rm = TRUE), .groups = "drop") %>%
      arrange(gauge_code, day_cut) %>%
      group_by(gauge_code) %>%
      mutate(rain_mm = lag(rain_mm, 1)) %>%
      ungroup() %>%
      filter(!is.na(rain_mm)) %>%
      rename(date = day_cut) %>%
      mutate(date = as.Date(date))
  } else {
    df_sub <- tibble(gauge_code = character(), date = as.Date(character()), rain_mm = numeric())
  }

  bind_rows(df_daily_native, df_sub) %>% arrange(gauge_code, date)
}

##### 6. Loop de leitura + processamento de todos os arquivos -----------------------------------------------------------
daily_list <- list()
info_list  <- list()

for (i in seq_len(nrow(file_index))) {
  f      <- file_index$file[i]
  path_i <- file.path(base_path, f)

  cat(sprintf("[%d/%d] Lendo %s...\n", i, nrow(file_index), f))

  df_data_i <- read_UNIPLU_BR(path_i, "table_data")
  df_info_i <- read_UNIPLU_BR(path_i, "table_info")

  if (nrow(df_data_i) == 0 || nrow(df_info_i) == 0) next

  info_list[[length(info_list) + 1]] <- df_info_i

  df_data_i$datetime <- as.POSIXct(df_data_i$datetime, tz = "UTC")
  daily_i <- to_daily_generic(df_data_i, df_info_i)

  if (nrow(daily_i) > 0) daily_list[[length(daily_list) + 1]] <- daily_i
}

# séries diárias consolidadas
df_daily <- bind_rows(daily_list) %>%
  distinct(gauge_code, date, .keep_all = TRUE) %>%
  arrange(gauge_code, date)

# metadados consolidados
df_info_total <- bind_rows(info_list) %>%
  distinct(gauge_code, .keep_all = TRUE)

cat("\nTotal de registros diários:", nrow(df_daily), "\n")
cat("Total de estações mapeadas:", nrow(df_info_total), "\n")

##### 7. Inventário das estações -----------------------------------------------------------
inventario <- df_daily %>%
  group_by(gauge_code) %>%
  summarise(
    data_inicio           = min(date),
    data_fim               = max(date),
    n_dias_com_dado        = n(),
    n_dias_periodo         = as.integer(data_fim - data_inicio) + 1,
    pct_completude         = round(100 * n_dias_com_dado / n_dias_periodo, 1),
    n_anos                 = round(n_dias_periodo / 365.25, 1),
    n_dias_chuva           = sum(rain_mm > LIMIAR_DIA_CHUVOSO, na.rm = TRUE),
    chuva_total_mm         = round(sum(rain_mm, na.rm = TRUE), 1),
    chuva_media_diaria_mm  = round(mean(rain_mm, na.rm = TRUE), 2),
    chuva_max_diaria_mm    = round(max(rain_mm, na.rm = TRUE), 1),
    data_chuva_max         = date[which.max(rain_mm)],
    .groups = "drop"
  ) %>%
  left_join(
    df_info_total %>%
      select(gauge_code, city, state, lat, long, elevation, time_step, network, responsible),
    by = "gauge_code"
  ) %>%
  relocate(city, state, lat, long, elevation, network, responsible, .after = gauge_code) %>%
  arrange(state, city)

write_csv(inventario, file.path(out_path, "inventario_estacoes.csv"))
cat("Inventário salvo em:", file.path(out_path, "inventario_estacoes.csv"), "\n")

##### 8. Séries diárias (salvar em disco) -----------------------------------------------------------
write_parquet(df_daily, file.path(out_path, "series_diarias.parquet"))
# alternativa em csv (arquivo maior, mas legível em qualquer software):
# write_csv(df_daily, file.path(out_path, "series_diarias.csv"))
cat("Séries diárias salvas em:", file.path(out_path, "series_diarias.parquet"), "\n")
