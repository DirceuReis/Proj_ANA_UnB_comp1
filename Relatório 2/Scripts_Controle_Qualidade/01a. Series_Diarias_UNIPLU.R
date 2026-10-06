## Séries Diárias - UNIPLU-BR
##

##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","stringr"))

##### 1. Library imports ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(stringr)

##### 2. Caminhos ----------------------------------------------------------------------------------------------------------
entrada_zip <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/Dados/UNIPLU_BR"
entrada_inv <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/00. Series"
saida       <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/01a. Series_Diarias"

parcial <- file.path(saida, "_parciais")   # um parquet por .zip; apagada no fim
dir.create(saida,   showWarnings = FALSE, recursive = TRUE)
dir.create(parcial, showWarnings = FALSE, recursive = TRUE)

t_inicio <- Sys.time()

##### 3. Parâmetros ----------------------------------------------------------------------------------------------------------
REDES_DIARIAS <- c("Hidroweb diário", "INMET diário")
REPROCESSAR   <- FALSE   # TRUE refaz tudo; FALSE pula os .zip já gravados
APAGAR_PARCIAIS <- TRUE  # apaga a pasta _parciais depois de juntar

##### 4. Função de leitura dos ZIPs ----------------------------------------------------------------------------------------------------------
read_UNIPLU_BR <- function(entrada_zip, table = "table_data") {
  parquet_file <- paste0(table, ".parquet")

  if (!file.exists(entrada_zip)) {
    message("Erro: arquivo não encontrado: ", entrada_zip)
    return(data.frame())
  }

  zip_contents <- tryCatch(
    unzip(entrada_zip, list = TRUE),
    error = function(e) {
      message("Erro ao abrir o ZIP: ", entrada_zip)
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

  unzip(entrada_zip, files = parquet_file, exdir = temp_dir)
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

##### 5. Inventário das estações diárias ----------------------------------------------------------------------------------------------------------
inventario_completo <- read_csv(file.path(entrada_inv, "inventario_estacoes.csv"),
                                col_types = cols(gauge_code = col_character(),
                                                 time_step  = col_character(),
                                                 .default   = col_guess()))

inventario_diario <- inventario_completo %>% filter(network %in% REDES_DIARIAS)

write_csv(inventario_diario, file.path(saida, "inventario_diario.csv"))

cat("Estações no inventário completo:", nrow(inventario_completo), "\n")
cat("Estações das redes diárias     :", nrow(inventario_diario), "\n")
print(count(inventario_diario, network, name = "estacoes"))

##### 6. Detectar arquivos disponíveis na pasta ----------------------------------------------------------------------------------------------------------
zip_files <- list.files(entrada_zip, pattern = "^[A-Za-z]{2}_[0-9]{4}\\.zip$")

file_index <- tibble(
  file  = zip_files,
  state = str_extract(zip_files, "^[A-Za-z]{2}"),
  year  = as.integer(str_extract(zip_files, "[0-9]{4}"))
)

states_filter <- NULL
years_filter  <- NULL

if (!is.null(states_filter)) file_index <- filter(file_index, state %in% states_filter)
if (!is.null(years_filter))  file_index <- filter(file_index, year  %in% years_filter)

file_index <- arrange(file_index, year, state)

cat("Arquivos encontrados:", length(zip_files), "| selecionados:", nrow(file_index), "\n")

##### 7. Obtém as estações diárias de cada .zip ----------------------------------------------------------------------------------------------------------
log_list <- vector("list", nrow(file_index))

for (i in seq_len(nrow(file_index))) {

  uf <- file_index$state[i]; ano <- file_index$year[i]
  arq_parcial <- file.path(parcial, sprintf("%s_%d.parquet", uf, ano))

  if (!REPROCESSAR && file.exists(arq_parcial)) {
    log_list[[i]] <- tibble(arquivo = file_index$file[i], estacoes = NA_integer_,
                            linhas = NA_integer_, situacao = "já processado")
    next
  }

  path_i    <- file.path(entrada_zip, file_index$file[i])
  df_info_i <- read_UNIPLU_BR(path_i, "table_info")

  if (nrow(df_info_i) == 0) {
    log_list[[i]] <- tibble(arquivo = file_index$file[i], estacoes = 0L, linhas = 0L,
                            situacao = "sem table_info")
    next
  }

  info_dia <- df_info_i %>% filter(network %in% REDES_DIARIAS)

  if (nrow(info_dia) == 0) {
    log_list[[i]] <- tibble(arquivo = file_index$file[i], estacoes = 0L, linhas = 0L,
                            situacao = "sem estação diária")
    next
  }

  df_data_i <- read_UNIPLU_BR(path_i, "table_data")

  if (nrow(df_data_i) == 0) {
    log_list[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_dia), linhas = 0L,
                            situacao = "sem table_data")
    next
  }

  # colunas pelo NOME: a ordem muda entre arquivos
  dados <- df_data_i %>%
    select(gauge_code, datetime, rain_mm) %>%
    filter(gauge_code %in% info_dia$gauge_code) %>%
    mutate(date = as.Date(as.POSIXct(datetime, tz = "UTC"))) %>%   # instante em UTC -> dia
    select(gauge_code, date, rain_mm) %>%
    arrange(gauge_code, date)

  if (nrow(dados) == 0) {
    log_list[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_dia), linhas = 0L,
                            situacao = "estações sem registro")
    next
  }

  write_parquet(dados, arq_parcial)

  log_list[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_dia),
                          linhas = nrow(dados), situacao = "ok")

  if (i %% 100 == 0 || i == nrow(file_index)) {
    el <- as.numeric(difftime(Sys.time(), t_inicio, units = "mins"))
    cat(sprintf("  %d/%d  (%.1f min; faltam ~%.1f min)\n",
                i, nrow(file_index), el, el / i * (nrow(file_index) - i)))
  }
}

log_processamento <- bind_rows(log_list)
write_csv(log_processamento, file.path(saida, "log_processamento.csv"))

##### 8. Elabora um parquet único ----------------------------------------------------------------------------------------------------------

etapa_juncao <- file.path(saida, "_juncao")
unlink(etapa_juncao, recursive = TRUE)

open_dataset(parcial) %>%
  write_dataset(etapa_juncao, format = "parquet",
                basename_template = "series_diarias-{i}.parquet")

arq_gerados <- list.files(etapa_juncao, pattern = "[.]parquet$", full.names = TRUE)
arq_final   <- file.path(saida, "series_diarias.parquet")

if (length(arq_gerados) == 1) {
  if (file.exists(arq_final)) file.remove(arq_final)
  file.rename(arq_gerados, arq_final)
  unlink(etapa_juncao, recursive = TRUE)
} else {
  stop("A junção gerou ", length(arq_gerados), " arquivos; esperado 1. Veja ", etapa_juncao)
}

if (APAGAR_PARCIAIS) unlink(parcial, recursive = TRUE)

##### 9. Resumo ----------------------------------------------------------------------------------------------------------
ds <- open_dataset(arq_final)

cat("\n===== RESUMO =====\n")
cat(sprintf("Arquivos lidos: %d | com estação diária: %d | sem: %d\n",
            nrow(log_processamento),
            sum(log_processamento$situacao == "ok"),
            sum(log_processamento$situacao == "sem estação diária")))
cat(sprintf("Série diária: %s registros | %d estações\n",
            format(ds$num_rows, big.mark = ".", decimal.mark = ","),
            nrow(inventario_diario)))