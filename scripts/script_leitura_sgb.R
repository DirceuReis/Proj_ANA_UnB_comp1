# Este script realiza a leitura dos dados disponibilizados pela Luna das estações do SGB.
# Registros estão disponíveis até 2025, mas não cobrem o ano todo como os últimos dados
# enviados pelo Filipe (UFPB).

# O conjunto do SGB consiste em estações com três tipos de sensor:
## 1. armazena dados sempre que chove considerando um intervalo mínimo que varia de
##    20 segundos até 20 minutos;
## 2. armazena dados sempre que chove, mas armazena os dados de forma acumulada;
## 3. registra o volume acumulado a cada 15 minutos.

# Todas as séries já estão com as datas preenchidas e as falhas registradas
# e apresentam o mesmo 'time_step' de 15 minutos.

rm(list = ls()); invisible(gc())

# Pacotes
pacman::p_load(pacman, tidyverse, arrow)

# Ler dados SGB (Luna) e transformar em .PARQUET
df.data <- readRDS(file = "base/fonte/bruto/subdiario/sgb/dados_sub2.RDS")
df.info <- readRDS(file = "base/fonte/bruto/subdiario/sgb/dados_sub_info2.rds")

# Formatar coluna 'gauge_code' como caractere e adicionar 0 no início
df.info <- df.info %>% 
  mutate(gauge_code = paste0("0", as.character(gauge_code)))

# Salvar arquivos
arrow::write_parquet(x = df.info, "base/fonte/consolidado/subdiario/RS_SGB_info.pqt")
arrow::write_parquet(x = df.data, "base/fonte/consolidado/subdiario/RS_SGB_data.pqt")