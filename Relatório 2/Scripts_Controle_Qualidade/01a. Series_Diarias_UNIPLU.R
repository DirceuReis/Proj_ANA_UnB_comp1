## Séries Diárias - UNIPLU-BR
##
## Filtra da série global a série de base diária - Hidroweb e Inmet diário - e faz a análise de consistência da base:
##   Seção 6  - consistência do cadastro (metadado estável ao longo dos anos?)
##   Seção 10 - disponibilidade (anos com >= 80% de dados; estações com >= 30 anos)
##   Seção 11 - duplicatas (par estação-dia repetido?)
##
## Insumos: os .zip do UNIPLU-BR e, vindos do 00a, "inventario_estacoes.csv"
##          e "resumo_estacao_ano.csv".
## Saídas:   "series_diarias.parquet", "series_diarias_filtrada.parquet",
##           "inventario_diario.csv" e os CSVs de consistência e disponibilidade.
##
##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","tidyr","arrow","stringr"))

##### 1. Library imports ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)
library(arrow)
library(stringr)

##### 2. Caminhos ----------------------------------------------------------------------------------------------------------
entrada_zip <- "C:/UNB/Dados/UNIPLU_BR"
entrada_inv <- "C:/UNB/UNIPLU/00. Series"
saida       <- "C:/UNB/UNIPLU/01a. Series_Diarias"

parcial <- file.path(saida, "_parciais")
dir.create(saida,   showWarnings = FALSE, recursive = TRUE)
dir.create(parcial, showWarnings = FALSE, recursive = TRUE)

t_inicio <- Sys.time()

##### 3. Parâmetros ----------------------------------------------------------------------------------------------------------
REDES_DIARIAS <- c("Hidroweb diário", "INMET diário")
REPROCESSAR   <- FALSE   # TRUE refaz tudo; FALSE pula os .zip já gravados
APAGAR_PARCIAIS <- TRUE  # apaga a pasta _parciais depois de juntar

## exclusão por inconsistência de cadastro, na Seção 6
DIST_MAX_KM          <- 10   # acima disso a estação mudou de sítio, não de precisão
MIN_ANOS_MESMA_COORD <- 30   # anos exigidos numa mesma coordenada para salvá-la

## filtro de disponibilidade da Seção 10
PCT_ANO_MINIMO <- 0.80   # fração mínima de dias com dado para o ano valer
MIN_ANOS_BONS  <- 30     # anos válidos mínimos para a estação ser mantida

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

cat("Estações no inventário completo:", nrow(inventario_completo), "\n")
cat("Estações das redes diárias     :", nrow(inventario_diario), "\n")
print(count(inventario_diario, network, name = "estacoes"))

##### 6. Consistência do cadastro ao longo dos anos ------------------------------
resumo_diario <- read_csv(file.path(entrada_inv, "resumo_estacao_ano.csv"),
                          col_types = cols(gauge_code = col_character(),
                                           time_step  = col_character(),
                                           .default   = col_guess())) %>%
  filter(network %in% REDES_DIARIAS)

cat("Resumo estação-ano das redes diárias:", nrow(resumo_diario), "linhas\n")

## distância entre dois pontos quaisquer, em km
dist_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

## 6.1 Uma linha por estação-ano, comparada com o ano anterior --------------------------------
serie_cadastro <- resumo_diario %>%
  mutate(lat = round(lat, 6), long = round(long, 6)) %>%
  arrange(gauge_code, ano) %>%
  group_by(gauge_code) %>%
  mutate(passo_km    = round(dist_km(lag(lat), lag(long), lat, long), 3),
         delta_alt_m = elevation - lag(elevation),
         passo_ant   = lag(time_step),
         city_ant    = lag(city),
         uf_ant      = lag(state),
         mudou_coord = !is.na(passo_km)    & passo_km > 0,
         mudou_alt   = !is.na(delta_alt_m) & delta_alt_m != 0,
         mudou_passo = !is.na(passo_ant)   & passo_ant != time_step,
         mudou_city  = !is.na(city_ant)    & city_ant  != city,
         mudou_uf    = !is.na(uf_ant)      & uf_ant    != state) %>%
  ungroup()

## 6.2 Os anos em que algo mudou --------------------------------------------------------------
eventos_diarios <- serie_cadastro %>%
  filter(mudou_coord | mudou_alt | mudou_passo | mudou_city | mudou_uf) %>%
  transmute(gauge_code, ano,
            passo_km, delta_alt_m,
            passo_de = passo_ant, passo_para = time_step,
            city_de  = city_ant,  city_para  = city,
            uf_de    = uf_ant,    uf_para    = state,
            o_que_mudou = paste0(if_else(mudou_coord, "coordenada ", ""),
                                 if_else(mudou_alt,   "altitude ",   ""),
                                 if_else(mudou_passo, "passo ",      ""),
                                 if_else(mudou_city,  "municipio ",  ""),
                                 if_else(mudou_uf,    "UF",          "")) %>% trimws(),
            lat, long, elevation) %>%
  left_join(inventario_diario %>% select(gauge_code, network), by = "gauge_code") %>%
  arrange(desc(passo_km), desc(abs(delta_alt_m)))

## 6.3 Resumo por estação--------------------------------------------------------------------
cadastro_por_estacao <- serie_cadastro %>%
  group_by(gauge_code) %>%
  summarise(
    ## coordenada
    n_coord        = n_distinct(paste(lat, long)),
    n_trocas_coord = sum(mudou_coord),
    voltou_coord   = any(duplicated(rle(paste(lat, long))$values)),
    ## altitude
    n_alt          = n_distinct(elevation),
    alt_ini        = first(elevation),
    alt_fim        = last(elevation),
    amplitude_m    = suppressWarnings(max(elevation, na.rm = TRUE) -
                                        min(elevation, na.rm = TRUE)),
    n_trocas_alt   = sum(mudou_alt),
    voltou_alt     = any(duplicated(rle(as.character(elevation))$values)),
    ## passo
    n_passo         = n_distinct(time_step),
    n_trocas_passo  = sum(mudou_passo),
    sequencia_passo = paste(rle(time_step)$values, collapse = " -> "),
    ## município e UF
    n_city         = n_distinct(city),
    city_ini       = first(city),
    city_fim       = last(city),
    n_trocas_city  = sum(mudou_city),
    voltou_city    = any(duplicated(rle(city)$values)),
    n_uf           = n_distinct(state),
    uf_ini         = first(state),
    uf_fim         = last(state),
    n_trocas_uf    = sum(mudou_uf),
    .groups = "drop") %>%
  mutate(amplitude_m = ifelse(is.finite(amplitude_m), amplitude_m, NA_real_))

## 6.4 Maior distância entre quaisquer duas coordenadas ---------------------------------------
coords <- serie_cadastro %>%
  filter(!is.na(lat), !is.na(long)) %>%
  distinct(gauge_code, lat, long)

dist_max <- coords %>%
  semi_join(filter(cadastro_por_estacao, n_coord > 1), by = "gauge_code") %>%
  inner_join(coords, by = "gauge_code", suffix = c("_a", "_b"),
             relationship = "many-to-many") %>%
  mutate(d = dist_km(lat_a, long_a, lat_b, long_b)) %>%
  group_by(gauge_code) %>%
  summarise(dist_max_km = max(d), .groups = "drop")

cadastro_por_estacao <- cadastro_por_estacao %>%
  left_join(dist_max, by = "gauge_code") %>%
  mutate(dist_max_km = round(coalesce(dist_max_km, 0), 3)) %>%
  left_join(inventario_diario %>% select(gauge_code, network), by = "gauge_code") %>%
  relocate(network, .after = gauge_code) %>%
  relocate(dist_max_km, .after = voltou_coord)

## 6.5 Matriz larga do deslocamento -----------------------------------------------------------
## vazio = sem dado no ano; 0 = não mudou; dois valores iguais seguidos = saiu e voltou
deslocamento_diario <- serie_cadastro %>%
  semi_join(filter(cadastro_por_estacao, n_trocas_coord > 0), by = "gauge_code") %>%
  select(gauge_code, ano, passo_km) %>%
  pivot_wider(names_from = ano, values_from = passo_km,
              values_fn = ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
              names_sort = TRUE) %>%
  left_join(cadastro_por_estacao %>%
              select(gauge_code, network, city_fim, uf_fim, n_coord, n_trocas_coord,
                     voltou_coord, dist_max_km), by = "gauge_code") %>%
  relocate(network, city_fim, uf_fim, n_coord, n_trocas_coord, voltou_coord, dist_max_km,
           .after = gauge_code)

## 6.6 Quais campos mudaram em cada estação ---------------------------------------------------
campos <- c("city", "state", "lat", "long", "elevation", "time_step",
            "network", "responsible", "UTC")

## n_distinct conta NA como valor
n_valores <- resumo_diario %>%
  mutate(lat = round(lat, 6), long = round(long, 6)) %>%
  group_by(gauge_code) %>%
  summarise(across(all_of(campos), ~ n_distinct(.x, na.rm = TRUE)), .groups = "drop")

campos_muda <- as.matrix(n_valores[campos]) > 1
n_valores$campos_que_mudam <- apply(campos_muda, 1, function(x) paste(campos[x], collapse = ", "))
n_valores$n_campos_mudam   <- rowSums(campos_muda)

primeiroano_info <- resumo_diario %>%
  arrange(gauge_code, ano) %>%
  distinct(gauge_code, .keep_all = TRUE) %>%
  select(gauge_code, ano_ini = ano, lat_ini = lat, long_ini = long)

ultimoano_info <- resumo_diario %>%
  arrange(gauge_code, desc(ano)) %>%
  distinct(gauge_code, .keep_all = TRUE) %>%
  select(gauge_code, ano_fim = ano, lat_fim = lat, long_fim = long)

mudancas_diarias <- n_valores %>%
  filter(n_campos_mudam > 0) %>%
  select(gauge_code, n_campos_mudam, campos_que_mudam) %>%
  left_join(primeiroano_info,     by = "gauge_code") %>%
  left_join(ultimoano_info,       by = "gauge_code") %>%
  left_join(cadastro_por_estacao, by = "gauge_code") %>%
  left_join(inventario_diario %>% select(gauge_code, n_anos_com_registro, n_anos_medido),
            by = "gauge_code") %>%
  relocate(network, city_ini, uf_ini, city_fim, uf_fim,
           n_anos_com_registro, n_anos_medido, n_campos_mudam, campos_que_mudam,
           .after = gauge_code) %>%
  arrange(desc(dist_max_km))

write_csv(eventos_diarios,       file.path(saida, "cadastro_eventos_diarios.csv"))
write_csv(cadastro_por_estacao,  file.path(saida, "cadastro_por_estacao_diarias.csv"))
write_csv(deslocamento_diario,   file.path(saida, "deslocamento_por_ano_diarias.csv"))
write_csv(mudancas_diarias,      file.path(saida, "mudancas_cadastro_diarias.csv"))

cat("Estações diárias com mudança de cadastro:", nrow(mudancas_diarias),
    "de", nrow(inventario_diario), "\n")
cat("Anos-estação em que algo mudou:", nrow(eventos_diarios), "\n")
print(count(eventos_diarios, o_que_mudou, name = "ocorrencias") %>%
        arrange(desc(ocorrencias)) %>% head(12))

cat("\nEstações afetadas, por atributo:\n")
print(tibble(
  atributo = c("coordenada", "altitude", "passo", "municipio", "UF"),
  estacoes = c(sum(cadastro_por_estacao$n_trocas_coord > 0),
               sum(cadastro_por_estacao$n_trocas_alt   > 0),
               sum(cadastro_por_estacao$n_trocas_passo > 0),
               sum(cadastro_por_estacao$n_trocas_city  > 0),
               sum(cadastro_por_estacao$n_trocas_uf    > 0))))

cat("\nDeslocamento da coordenada (maior distância entre dois anos):\n")
print(cadastro_por_estacao %>%
        filter(dist_max_km > 0) %>%
        count(faixa = cut(dist_max_km, c(0, 0.1, 1, 10, 50, Inf),
                          labels = c("até 100 m", "100 m a 1 km", "1 a 10 km",
                                     "10 a 50 km", "mais de 50 km")),
              name = "estacoes"))

cat("\nAs que se deslocam mais de 1 km:\n")
print(mudancas_diarias %>% filter(dist_max_km > 1) %>%
        select(gauge_code, city_ini, uf_ini, city_fim, uf_fim, dist_max_km))

## 6.7 Exclusão por mudança de sítio ----------------------------------------------------------
## Deslocamento acima de DIST_MAX_KM é troca de local
## A estação só se salva se reunir MIN_ANOS_MESMA_COORD anos em uma única coordenada,
anos_por_coord <- serie_cadastro %>%
  filter(!is.na(lat), !is.na(long)) %>%
  count(gauge_code, lat, long, name = "anos_na_coord")

max_anos_coord <- anos_por_coord %>%
  group_by(gauge_code) %>%
  summarise(max_anos_coord = max(anos_na_coord), .groups = "drop")

estacoes_excluidas_cadastro <- cadastro_por_estacao %>%
  left_join(max_anos_coord, by = "gauge_code") %>%
  filter(dist_max_km > DIST_MAX_KM,
         !is.na(max_anos_coord),
         max_anos_coord < MIN_ANOS_MESMA_COORD) %>%
  transmute(gauge_code, network, city_ini, uf_ini, city_fim, uf_fim,
            dist_max_km, n_coord, max_anos_coord,
            motivo = sprintf("%.0f km de deslocamento e no máximo %d anos numa mesma coordenada",
                             dist_max_km, max_anos_coord))

write_csv(anos_por_coord,              file.path(saida, "anos_por_coordenada_diarias.csv"))
write_csv(estacoes_excluidas_cadastro, file.path(saida, "estacoes_excluidas_cadastro.csv"))

inventario_diario <- inventario_diario %>%
  filter(!gauge_code %in% estacoes_excluidas_cadastro$gauge_code)

write_csv(inventario_diario, file.path(saida, "inventario_diario.csv"))

cat(sprintf("\nExcluídas por mudança de sítio (> %g km e < %d anos numa mesma coordenada): %d\n",
            DIST_MAX_KM, MIN_ANOS_MESMA_COORD, nrow(estacoes_excluidas_cadastro)))
if (nrow(estacoes_excluidas_cadastro) > 0)
  print(estacoes_excluidas_cadastro %>%
          select(gauge_code, network, city_fim, uf_fim, dist_max_km, max_anos_coord))
cat("Inventário diário após a exclusão:", nrow(inventario_diario), "estações\n")

##### 7. Detectar arquivos disponíveis na pasta ----------------------------------------------------------------------------------------------------------
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

##### 8. Obtém as estações diárias de cada .zip ----------------------------------------------------------------------------------------------------------
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

  info_dia <- df_info_i %>% filter(gauge_code %in% inventario_diario$gauge_code)

  if (nrow(info_dia) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = 0L, linhas = 0L,
                            situacao = "sem estação diária")
    next
  }

  df_data_i <- read_UNIPLU_BR(path_i, "table_data")

  if (nrow(df_data_i) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_dia), linhas = 0L,
                            situacao = "sem table_data")
    next
  }

  dados <- df_data_i %>%
    select(gauge_code, datetime, rain_mm) %>%
    filter(gauge_code %in% info_dia$gauge_code) %>%
    mutate(date = as.Date(as.POSIXct(datetime, tz = "UTC"))) %>%   # instante em UTC -> dia
    select(gauge_code, date, rain_mm) %>%
    arrange(gauge_code, date)

  if (nrow(dados) == 0) {
    verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_dia), linhas = 0L,
                            situacao = "estações sem registro")
    next
  }

  write_parquet(dados, arq_parcial)

  verificacao_zip[[i]] <- tibble(arquivo = file_index$file[i], estacoes = nrow(info_dia),
                          linhas = nrow(dados), situacao = "ok")

  if (i %% 100 == 0 || i == nrow(file_index)) {
    el <- as.numeric(difftime(Sys.time(), t_inicio, units = "mins"))
    cat(sprintf("  %d/%d  (%.1f min; faltam ~%.1f min)\n",
                i, nrow(file_index), el, el / i * (nrow(file_index) - i)))
  }
}

verificacao_files <- bind_rows(verificacao_zip)
write_csv(verificacao_files, file.path(saida, "verificacao_files.csv"))

##### 9. Junção em um parquet único ----------------------------------------------------------------------------------------------------------
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

##### 10. Filtro de disponibilidade -------------------------------------------------------------------
## (i)  mantém o ano cuja fração de dias com dado alcança PCT_ANO_MINIMO = 80;
## (ii) mantém a estação que reúne pelo menos MIN_ANOS_BONS = 30.
cobertura <- open_dataset(arq_final) %>%
  filter(!is.na(date)) %>%
  group_by(gauge_code, ano = as.integer(year(date))) %>%   # as.integer: year() devolve int64
  summarise(dias_com_dado = sum(!is.na(rain_mm)), .groups = "drop") %>%
  collect() %>%
  mutate(dias_no_ano = as.integer(as.Date(sprintf("%d-12-31", ano)) -
                                    as.Date(sprintf("%d-01-01", ano))) + 1L,
         frac_ano    = dias_com_dado / dias_no_ano,
         ano_valido  = frac_ano >= PCT_ANO_MINIMO)

anos_validos <- cobertura %>% filter(ano_valido)

estacoes_validas <- anos_validos %>%
  count(gauge_code, name = "n_anos_validos") %>%
  filter(n_anos_validos >= MIN_ANOS_BONS)

pares_validos <- anos_validos %>%
  semi_join(estacoes_validas, by = "gauge_code") %>%
  select(gauge_code, ano)

write_csv(cobertura,        file.path(saida, "cobertura_estacao_ano.csv"))
write_csv(estacoes_validas, file.path(saida, "estacoes_filtro_disponibilidade.csv"))

## grava a série filtrada
etapa_filtro <- file.path(saida, "_filtro")
unlink(etapa_filtro, recursive = TRUE)

open_dataset(arq_final) %>%
  mutate(ano = as.integer(year(date))) %>%
  inner_join(pares_validos, by = c("gauge_code", "ano")) %>%
  select(gauge_code, date, rain_mm) %>%
  write_dataset(etapa_filtro, format = "parquet",
                basename_template = "series_diarias_filtrada-{i}.parquet")

arq_filtro_gerados <- list.files(etapa_filtro, pattern = "[.]parquet$", full.names = TRUE)
arq_filtrado       <- file.path(saida, "series_diarias_filtrada.parquet")

if (length(arq_filtro_gerados) == 1) {
  if (file.exists(arq_filtrado)) file.remove(arq_filtrado)
  file.rename(arq_filtro_gerados, arq_filtrado)
  unlink(etapa_filtro, recursive = TRUE)
} else {
  stop("O filtro gerou ", length(arq_filtro_gerados), " arquivos; esperado 1. Veja ", etapa_filtro)
}

cat(sprintf("\nFiltro (i) anos com >= %.0f%% de dados: %d de %d estação-anos\n",
            100 * PCT_ANO_MINIMO, nrow(anos_validos), nrow(cobertura)))
cat(sprintf("Filtro (ii) estações com >= %d anos válidos: %d de %d\n",
            MIN_ANOS_BONS, nrow(estacoes_validas), nrow(inventario_diario)))
cat(sprintf("Série filtrada: %s registros | %d estação-anos\n",
            format(open_dataset(arq_filtrado)$num_rows, big.mark = ".", decimal.mark = ","),
            nrow(pares_validos)))

##### 11. Duplicatas de estação, data e precipitação -----------------------------------------
## Confere, ano a ano, se a série filtrada tem mais de um registro para o mesmo par estação/dia
## Duplicata Exata - mesma estação, mesmo dia, mesma chuva.
## Duplicata Diff  - mesma estação, mesmo dia, chuva diferente.

ds_diario <- open_dataset(arq_filtrado)
mil <- function(x) format(x, big.mark = ".", decimal.mark = ",")

faixa <- ds_diario %>%
  summarise(ini = min(date, na.rm = TRUE), fim = max(date, na.rm = TRUE),
            sem_data = sum(is.na(date))) %>%   # usado no relatório da Seção 11
  collect()

anos_serie <- as.integer(format(faixa$ini, "%Y")):as.integer(format(faixa$fim, "%Y"))

lista <- vector("list", length(anos_serie))

for (i in seq_along(anos_serie)) {
  
  d <- ds_diario %>%
    filter(date >= as.Date(sprintf("%d-01-01", anos_serie[i])),
           date <= as.Date(sprintf("%d-12-31", anos_serie[i]))) %>%
    collect()
  if (nrow(d) == 0) next
  
  chave    <- paste(d$gauge_code, d$date)
  repetida <- duplicated(chave) | duplicated(chave, fromLast = TRUE)
  
  if (any(repetida))
    lista[[i]] <- d[repetida, ] %>%
    group_by(gauge_code, date) %>%
    summarise(n_linhas = n(), n_valores = n_distinct(rain_mm), .groups = "drop")
  
  rm(d, chave, repetida); invisible(gc())
}

duplicatas <- bind_rows(lista)

cat(sprintf("Registros: %s | período: %d a %d | sem data: %s\n",
            mil(ds_diario$num_rows), min(anos_serie), max(anos_serie), mil(faixa$sem_data)))

if (nrow(duplicatas) == 0) {
  cat("Nenhum par (estação, dia) repetido.\n")
} else {
  # chuva igual = duplicata exata; diferente = conflito
  duplicatas <- duplicatas %>% mutate(exata = n_valores == 1)
  write_csv(duplicatas, file.path(saida, "duplicatas_diarias.csv"))
  cat(sprintf("Pares repetidos: %s | exatos: %s | conflitantes: %s\n",
              mil(nrow(duplicatas)), mil(sum(duplicatas$exata)), mil(sum(!duplicatas$exata))))
  cat(sprintf("Linhas excedentes: %s | estações afetadas: %s\n",
              mil(sum(duplicatas$n_linhas) - nrow(duplicatas)),
              mil(n_distinct(duplicatas$gauge_code))))
}

##### 12. Resumo ----------------------------------------------------------------------------------------------------------
ds <- open_dataset(arq_final)

cat("\nArquivos lidos:", nrow(verificacao_files), "\n")
print(count(verificacao_files, situacao, name = "zips") %>% arrange(desc(zips)))
cat(sprintf("Série completa: %s registros | %d estações\n",
            format(ds$num_rows, big.mark = ".", decimal.mark = ","),
            nrow(inventario_diario)))
cat(sprintf("Série filtrada: %s registros | %d estações\n",
            format(open_dataset(arq_filtrado)$num_rows, big.mark = ".", decimal.mark = ","),
            nrow(estacoes_validas)))
