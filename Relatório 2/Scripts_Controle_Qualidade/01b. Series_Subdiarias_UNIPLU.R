## Séries Sub-Diárias - UNIPLU-BR
##
## Extrai dos .zip do UNIPLU-BR a série das redes sub-diárias - CEMADEN,
## Telemetria, INMET sub-diário e ICEA - e a consolida em um parquet único.
##
## Insumos: os .zip do UNIPLU-BR e, vindo do 00a, "inventario_estacoes.csv".
## Saídas:  "series_subdiarias.parquet", "inventario_subdiario.csv" e
##          "verificacao_files.csv".
##
##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","stringr"))

##### 1. Library imports ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(stringr)

##### 2. Caminhos ----------------------------------------------------------------------------------------------------------
entrada_zip <- "C:/UNB/Dados/UNIPLU_BR"
entrada_inv <- "C:/UNB/UNIPLU/00. Series"
saida       <- "C:/UNB/UNIPLU/01b. Series_Subdiarias"

parcial <- file.path(saida, "_parciais")   # um parquet por .zip
dir.create(saida,   showWarnings = FALSE, recursive = TRUE)
dir.create(parcial, showWarnings = FALSE, recursive = TRUE)

t_inicio <- Sys.time()

##### 3. Parâmetros ----------------------------------------------------------------------------------------------------------
REDES_DIARIAS      <- c("Hidroweb diário", "INMET diário")
REPROCESSAR        <- FALSE   # TRUE refaz tudo; FALSE pula os .zip já gravados
APAGAR_PARCIAIS    <- TRUE   # mantenha FALSE até conferir o resultado
SAIDA_PARTICIONADA <- FALSE   # TRUE: pastas state=UF/ano=AAAA | FALSE: parquet único

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

##### 5. Inventário das estações sub-diárias ----------------------------------------------------------------------------------------------------------
inventario_completo <- read_csv(file.path(entrada_inv, "inventario_estacoes.csv"),
                                col_types = cols(gauge_code = col_character(),
                                                 time_step  = col_character(),
                                                 .default   = col_guess()))

inventario_subdiario <- inventario_completo %>% filter(!network %in% REDES_DIARIAS)

write_csv(inventario_subdiario, file.path(saida, "inventario_subdiario.csv"))

cat("Estações no inventário completo   :", nrow(inventario_completo), "\n")
cat("Estações das redes sub-diárias    :", nrow(inventario_subdiario), "\n")
print(count(inventario_subdiario, network, name = "estacoes"))

##### 6. Detectar arquivos disponíveis na pasta ----------------------------------------------------------------------------------------------------------
zip_files <- list.files(entrada_zip, pattern = "^[A-Za-z]{2}_[0-9]{4}\\.zip$")

file_index <- tibble(
  file  = zip_files,
  state = str_extract(zip_files, "^[A-Za-z]{2}"),
  year  = as.integer(str_extract(zip_files, "[0-9]{4}"))
)

states_filter <- NULL   # ex.: c("SP", "AC")  |  NULL = todos os estados
years_filter  <- NULL   # ex.: 2014:2025     |  NULL = todos os anos

if (!is.null(states_filter)) file_index <- filter(file_index, state %in% states_filter)
if (!is.null(years_filter))  file_index <- filter(file_index, year  %in% years_filter)

file_index <- arrange(file_index, year, state)

cat("Arquivos encontrados:", length(zip_files), "| selecionados:", nrow(file_index), "\n")

##### 7. Obtém as estações sub-diárias de cada .zip ----------------------------------------------------------------------------------------------------------
verificacao_zip <- vector("list", nrow(file_index))

for (i in seq_len(nrow(file_index))) {

  uf <- file_index$state[i]; ano <- file_index$year[i]
  arq_parcial <- file.path(parcial, sprintf("%s_%d.parquet", uf, ano))

  if (!REPROCESSAR && file.exists(arq_parcial)) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = NA_integer_,
                                   linhas = NA_integer_, situacao = "já processado")
    next
  }

  path_i    <- file.path(entrada_zip, file_index$file[i])
  df_info_i <- read_UNIPLU_BR(path_i, "table_info")

  if (nrow(df_info_i) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = 0L, linhas = 0L,
                                   situacao = "sem table_info")
    next
  }

  info_sub <- df_info_i %>% filter(!network %in% REDES_DIARIAS)

  if (nrow(info_sub) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = 0L, linhas = 0L,
                                   situacao = "sem estação sub-diária")
    next
  }

  df_data_i <- read_UNIPLU_BR(path_i, "table_data")

  if (nrow(df_data_i) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_sub), linhas = 0L,
                                   situacao = "sem table_data")
    next
  }

  dados <- df_data_i %>%
    select(gauge_code, datetime, rain_mm) %>%
    filter(gauge_code %in% info_sub$gauge_code) %>%
    arrange(gauge_code, datetime)

  if (nrow(dados) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_sub), linhas = 0L,
                                   situacao = "estações sem registro")
    next
  }

  attr(dados$datetime, "tzone") <- "UTC"

  dados <- dados %>% mutate(state = uf, ano = ano)

  write_parquet(dados, arq_parcial)

  verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_sub),
                                 linhas = nrow(dados), situacao = "ok")

  if (i %% 100 == 0 || i == nrow(file_index)) {
    el <- as.numeric(difftime(Sys.time(), t_inicio, units = "mins"))
    cat(sprintf("  %d/%d  (%.1f min; faltam ~%.1f min)\n",
                i, nrow(file_index), el, el / i * (nrow(file_index) - i)))
  }
}

verificacao_files <- bind_rows(verificacao_zip)
write_csv(verificacao_files, file.path(saida, "verificacao_files.csv"))

##### 8. Junta os parciais ----------------------------------------------------------------------------------------------------------
if (SAIDA_PARTICIONADA) {

  destino <- file.path(saida, "series_subdiarias")
  unlink(destino, recursive = TRUE)

  open_dataset(parcial) %>%
    group_by(state, ano) %>%
    write_dataset(destino, format = "parquet", basename_template = "dados-{i}.parquet")

} else {

  etapa_juncao <- file.path(saida, "_juncao")
  unlink(etapa_juncao, recursive = TRUE)

  open_dataset(parcial) %>%
    write_dataset(etapa_juncao, format = "parquet",
                  basename_template = "series_subdiarias-{i}.parquet")

  arq_gerados <- list.files(etapa_juncao, pattern = "[.]parquet$", full.names = TRUE)
  arq_final   <- file.path(saida, "series_subdiarias.parquet")

  if (length(arq_gerados) == 1) {
    if (file.exists(arq_final)) file.remove(arq_final)
    file.rename(arq_gerados, arq_final)
    unlink(etapa_juncao, recursive = TRUE)
  } else {
    stop("A junção gerou ", length(arq_gerados), " arquivos; esperado 1. Veja ", etapa_juncao)
  }
}

if (APAGAR_PARCIAIS) unlink(parcial, recursive = TRUE)