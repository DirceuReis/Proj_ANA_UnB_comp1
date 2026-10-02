## Extração e Caracterização (convencionais e automáticas) - UNIPLU-BR

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
dados <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/Dados/UNIPLU_BR"

# Pasta de saída
saida <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/00. Series"
dir.create(saida, showWarnings = FALSE, recursive = TRUE)

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

zip_files <- list.files(dados, pattern = "^[A-Za-z]{2}_[0-9]{4}\\.zip$")

file_index <- tibble(
  file  = zip_files,
  state = str_extract(zip_files, "^[A-Za-z]{2}"),
  year  = as.integer(str_extract(zip_files, "[0-9]{4}"))
)

states_filter <- NULL
years_filter  <- NULL

if (!is.null(states_filter)) file_index <- filter(file_index, state %in% states_filter)
if (!is.null(years_filter))  file_index <- filter(file_index, year  %in% years_filter)

file_index <- arrange(file_index, state, year)

cat("Arquivos encontrados na pasta:", length(zip_files), "\n")
cat("Arquivos selecionados para processamento:", nrow(file_index), "\n")
cat("Estados:", paste(unique(file_index$state), collapse = ", "), "\n")

##### 5. Resumo -----------------------------------------------------------
resumir_zip <- function(df_data, df_info, ano) {
  
  df_data %>%
    select(gauge_code, datetime, rain_mm) %>%
    group_by(gauge_code) %>%
    summarise(
      n_registros    = n(),
      n_na           = sum(is.na(rain_mm)),
      primeiro       = if (all(is.na(datetime))) as.POSIXct(NA, tz = "UTC") else min(datetime, na.rm = TRUE),
      ultimo         = if (all(is.na(datetime))) as.POSIXct(NA, tz = "UTC") else max(datetime, na.rm = TRUE),
      n_dias         = n_distinct(as.Date(datetime[!is.na(datetime)])),
      chuva_total_mm = sum(rain_mm, na.rm = TRUE),
      chuva_max      = suppressWarnings(max(rain_mm, na.rm = TRUE)),
      data_chuva_max = if (all(is.na(rain_mm)) || all(is.na(datetime))) as.POSIXct(NA, tz = "UTC")
      else datetime[which.max(rain_mm)],
      n_reg_chuva    = sum(rain_mm > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      chuva_max       = ifelse(is.finite(chuva_max), chuva_max, NA_real_),
      sem_dado_valido = n_registros == n_na,
      ano             = ano
    ) %>%
    left_join(
      df_info %>% 
        distinct(gauge_code, .keep_all = TRUE) %>%
        select(gauge_code, city, state, lat, long, elevation,
               time_step, network, responsible, UTC),
      by = "gauge_code"
    )
}

##### 6. Loop de leitura + resumo de todos os arquivos -----------------------------------------------------------

resumos <- vector("list", nrow(file_index))

for (i in seq_len(nrow(file_index))) {
  
  path_i <- file.path(dados, file_index$file[i])
  
  df_info_i <- read_UNIPLU_BR(path_i, "table_info")
  df_data_i <- read_UNIPLU_BR(path_i, "table_data")
  
  if (nrow(df_info_i) == 0 || nrow(df_data_i) == 0) next
  
  # o dado bruto não traz fuso; o instante é rotulado UTC, sem deslocamento
  df_data_i$datetime <- as.POSIXct(df_data_i$datetime, tz = "UTC")
  
  resumos[[i]] <- resumir_zip(df_data_i, df_info_i, file_index$year[i])
  
  if (i %% 100 == 0 || i == nrow(file_index))
    cat(sprintf("  %d/%d arquivos\n", i, nrow(file_index)))
}

resumo_por_ano <- bind_rows(resumos)
write_csv(resumo_por_ano, file.path(saida, "resumo_estacao_ano.csv"))

cat("Resumo por estação e ano:", nrow(resumo_por_ano), "linhas\n")
cat("Estação-anos sem nenhuma medição:", sum(resumo_por_ano$sem_dado_valido), "\n")

##### 7. Inventário das estações ----------------------------------------------------------------

inventario <- resumo_por_ano %>%
  group_by(gauge_code) %>%
  summarise(
    n_registros    = sum(n_registros),
    n_na           = sum(n_na),
    n_dias         = sum(n_dias),
    n_anos         = n_distinct(ano),                  # anos com algum registro
    ano_inicio     = min(ano),
    ano_fim        = max(ano),
    span_anos      = max(ano) - min(ano) + 1,
    data_inicio    = suppressWarnings(min(primeiro, na.rm = TRUE)),
    data_fim       = suppressWarnings(max(ultimo,   na.rm = TRUE)),
    chuva_total_mm = round(sum(chuva_total_mm, na.rm = TRUE), 1),
    chuva_max      = suppressWarnings(max(chuva_max, na.rm = TRUE)),
    data_chuva_max = if (all(is.na(chuva_max))) as.POSIXct(NA, tz = "UTC")
    else data_chuva_max[which.max(chuva_max)],
    n_reg_chuva    = sum(n_reg_chuva),
    .groups = "drop"
  ) %>%
  mutate(
    chuva_max       = ifelse(is.finite(chuva_max), chuva_max, NA_real_),
    data_inicio     = as.POSIXct(ifelse(is.finite(data_inicio), data_inicio, NA), tz = "UTC"),
    data_fim        = as.POSIXct(ifelse(is.finite(data_fim),    data_fim,    NA), tz = "UTC"),
    sem_dado_valido = n_registros == n_na              # estação sem nenhuma medição
  ) %>%
  left_join(
    resumo_por_ano %>%
      arrange(gauge_code, desc(ano)) %>%               # metadados do ano mais recente
      distinct(gauge_code, .keep_all = TRUE) %>%
      select(gauge_code, city, state, lat, long, elevation,
             time_step, network, responsible, UTC),
    by = "gauge_code"
  ) %>%
  arrange(state, city)

write_csv(inventario, file.path(saida, "inventario_estacoes.csv"))

cat("Inventário salvo:", nrow(inventario), "estações\n")
cat("Estações sem nenhuma medição:", sum(inventario$sem_dado_valido), "\n")
print(inventario %>% count(network, name = "estacoes") %>% arrange(desc(estacoes)))

