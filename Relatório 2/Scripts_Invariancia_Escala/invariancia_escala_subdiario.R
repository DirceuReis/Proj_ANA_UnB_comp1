
# Script para análise empírica da hipótese de invariância de escala simples
# em sentido amplo, utilizando os momentos não centrais de diferentes ordens
# para a faixa de validade de 10 minutos (quando disponível) até 24 horas.

rm(list = ls()); invisible(gc())

# PACOTES ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(boot)
library(arrow)
library(pbapply)
library(purrr)
library(parallelly)


# CAMINHOS ---------------------------------------------------------------

path.root <- "Relatório 2"
path.imax <- file.path(path.root, "Scripts_Maximas", "resultados", "BR", "dataframes", "imax_sem_iqr_subdiario_Uniplu.parquet")
path.inv.scale <- file.path(path.root, "Scripts_Invariancia_Escala", "resultados", "subdiarios", "invariancia_escala")
path.funs <- file.path(path.root, "Scripts_Invariancia_Escala", "funs", "fun_scale_invariance.R")
path.out <- file.path(path.root, "Scripts_Invariancia_Escala", "resultados", "subdiarios")


# CARREGAR FUNCOES -------------------------------------------------------

if(!file.exists(path.funs)) stop("Arquivo 'fun_scale_invariance.R' não encontrado.")
source(path.funs)


# LER DADOS --------------------------------------------------------------

# Intensidades maximas anuais com dados (UNIPLU-BR)
imax.preqc.df <- arrow::read_parquet(path.imax)

# Dados ja foram filtrados p/ 20% de dados e anos "incompletos"
na.accept <- 0.2   # por garantia
min.years <- 8     # pelo menos 8 anos completos
min.durations <- 3 # nro. minimo de duracoes p/ ajustar d-GEV

# Filtrar somente estacoes c/ pelo menos 'min.years' anos
imax.df <- imax.preqc.df |> 
  filter(na_prct_yr <= na.accept) |> 
  group_by(gauge_code) |> 
  mutate(n_years = n_distinct(wateryear), .after = wateryear) |> 
  filter(
    n_distinct(d) >= min.durations,
    n_distinct(wateryear) >= min.years
  ) |> 
  ungroup() |> 
  rename(imax = imax_mm_h) |> 
  select(gauge_code, wateryear, d, imax)


# AJUSTE INVARIANCIA SIMPLES SENTIDO AMPLO -------------------------------

# Argumentos
which.moments <- 1:3              # quais momentos usar de base
min.duration <- 3                 # minimo de duracoes diferentes p/ ajustar o modelo linear
duration.interval <- c(10/60, 24) # 10 minutos até 24 horas

# Ajuste simpels p/ recolher H, R² e os momentos estimados
inv.ls <- fit.simple.scaling(
  df.imax = imax.df,
  which.moment = which.moments,
  which.duration = duration.interval,
  min.duration = min.duuration
)

# Extrair momentos e coeficientes
moments.df <- inv.ls$scale.moments
scale.df <- inv.ls$scale.coefficient

if(!dir.exists(path.out)) dir.create(path.out, recursive = TRUE)
arrow::write_parquet(moments.df, file.path(path.out, "df_momentos_subdiario.pqt"))
arrow::write_parquet(scale.df, file.path(path.out, "df_expoente_escala_subdiario.pqt"))

# Ajuste com intervalos de confiança usando bootstrap
inv.boot.ls <- simple.scaling.boot(
  df.imax = imax.df,
  which.moment = which.moments,
  which.duration = duration.interval,
  min.duration = min.duration,
  d.target = NULL, # não calcular coeficientes de desagregação
  ci.boot.method = "perc",
  R.boot = 1e4,
  return.boot = TRUE
)
