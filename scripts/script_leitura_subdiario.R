
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Instalar e carregar pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,
               tidyverse, # manipulação de dados
               lubridate, # manipulação de datas
               arrow,     # ler dados no formato 'parquet'
               ggplot2)   # visualização

# LER DADOS ---------------------------------------------------------------

# Caminho do diretório
path.dir <- "base/fonte/consolidado/subdiario/"


# Listar arquivos c/ extensão ".pqt" no diretório
# Os arquivos contém divisões por ano de duas tabelas: info e data

## Listar 'data'
files.data <- list.files(path = path.dir,           
                         pattern = "*.data.pqt",
                         full.names = TRUE)


## Listar 'info'
files.info <- list.files(path = path.dir,
                         pattern = "*.info.pqt",
                         full.names = TRUE)

# Ler dados
# list_data <- lapply(X = files_data, FUN = arrow::read_parquet, as_data_frame = TRUE)
list.data <- lapply(X = files.data, FUN = arrow::read_parquet)
list.info <- lapply(X = files.info, FUN = arrow::read_parquet)

# Juntar dados
df.data <- bind_rows(list.data)[-4]  # remover coluna '__intex_level_0__'
df.info <- bind_rows(list.info)[-10] # remover coluna '__intex_level_0__'

# Rodar somente se houver tanto colunas 'monitoring_time_step' e 'time_step'
df.info <- df.info %>% 
  mutate(time_step = coalesce(monitoring_time_step, as.numeric(time_step))) %>% 
  select(-monitoring_time_step)

# Manter somente 'df_data' e 'df_info' no ambiente
rm(list = setdiff(ls(), c("df.data", "df.info"))); invisible(gc())


# ORGANIZAR ---------------------------------------------------------------

# O arquivo 'df_info' contém estações repetidas, portanto precisamos manter
# somente um único registro de cada 'gauge_code'
df.info <- slice_tail(.data = df.info, n = 1, by = "gauge_code")

# Fuso horário
# Todos os dados do conjunto deveriam estar na zona UTC -3
# Isso corrige o shift de -2 ou -3 horas nos dados
df.data <- df.data %>% 
  mutate(
    datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), # alterar fuso
    time_step = df.info$time_step[match(gauge_code, df.info$gauge_code)],      # resolução temporal
    responsible = df.info$responsible[match(gauge_code, df.info$gauge_code)]   # operador
  ); invisible(gc()) # resolução temporal


# EXPORTAR BRUTO ----------------------------------------------------------

# Parquet
write_parquet(df.data, sink = "base/gerados/df_subdaily_data.pqt")
write_parquet(df.info, sink = "base/gerados/df_subdaily_info.pqt")

invisible(gc())

# PREENCHER DATAS ---------------------------------------------------------

# # Ler dados gerados em `rs_subdaily_read_parquet.R`
# df_data <- arrow::read_parquet("base/gerados/df_subdaily_data.parquet")
# df_info <- arrow::read_parquet("base/gerados/df_subdaily_info.parquet")
# 
# # Função completa
# source("scripts/funcoes/fun_fill_dates.R")
# list_subdaily_filled <- fun_fill_dates(df = df_data, col_names = names(df_data)); gc()
# df_subdaily_filled <- bind_rows(list_subdaily_filled); rm(list_subdaily_filled); gc()
# 
# # Salvar arquivo (deixar comentado)
# list_by_resp <- split(df_subdaily_filled, f = df_subdaily_filled$responsible); rm(df_subdaily_filled); gc()
# compression_codecs <- c("snappy", "gzip", "brotli", "zstd", "lz4", "lzo", "bz2")
# write_parquet(x = df_subdaily_filled,
#               sink = "base/gerados/df_subadaily_filled.parquet",
#               compression = compression_codecs[[4]],
#               compression_level = 9)