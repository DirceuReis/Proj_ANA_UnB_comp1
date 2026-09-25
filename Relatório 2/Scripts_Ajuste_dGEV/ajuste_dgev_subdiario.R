
# Os scripts atuais da d-GEV utilizam a parametrização da distribuição
# em que o parâmetro de forma negativo representa a distribuição de 
# cauda pesada (Fréchet), mas o nome dos parâmetros foi alterado para
# corresponder à nomenclatura que estamos utilizando

# PACOTES ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(pbapply)
library(arrow)


# CAMINHOS ---------------------------------------------------------------

path.root <- "Relatório 2"
path.funs <- file.path(path.root, "Scripts_Ajuste_dGEV", "funs", "fit_dgev.R")
path.imax <- file.path(path.root, "Scripts_Maximas", "resultados", "BR", "dataframes", "df_imax_subdiario.parquet")
# path.imax <- file.path(path.root, "Scripts_Ajuste_dGEV", "resultados", "subdiarios", "dataframes", "df_imax_subdiario.pqt")
path.fit.ls <- file.path(path.root, "Scripts_Maximas", "resultados", "subdiarios", "fit")


# CARREGAR FUNCOES -------------------------------------------------------

if(!file.exists(path.funs)) stop("Arquivo 'fit_dgev.R' não encontrado.")
source(path.funs)


# LER DADOS --------------------------------------------------------------

# Intensidades maximas anuais com dados subdiarios (BR-SDR)
imax.preqc.df <- arrow::read_parquet(path.imax)

# Dados ja foram filtrados p/ 20% de dados e anos "incompletos"
na.accept <- 0.2   # por garantia
min.years <- 8    # pelo menos 8 anos completos
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
  select(gauge_code, wateryear, d, imax)


# AJUSTE d-GEV -----------------------------------------------------------

# Argumentos
durations <- NULL       # usa as duracoes disponiveis
prior <- c(-0.1, 0.122) # distribuicao a priori de Martins e Stedinger (2000)
cols <- names(imax.df)
scale.inv <- "wide"     # estimativas iniciais dos parametros de invariancia simples
return.period <- c(2, 5, 10, 25, 50, 100, 200) # tempos de retorno
conf <- 0.95            # nivel de confianca p/ intervalos de confianca
maxit <- 1e8            # maximo de iteracoes p/ convergencia
step <- 0.005           # passo p/ encontrar as raizes da funcao de verossimilhanca perfilada
lltol <- 0.005          # tolerancia p/ diferenca entre maximas log-verossimilhancas

# Ajustar d-GEV, estimar quantis e intervalos de confiança com verossimilhança perfilada
imax.ls <- split(imax.df, imax.df$gauge_code); invisible(gc())

fit.ls <- fitprof.dgev(
  data = imax.ls,
  prior.info = prior,
  return.period = return.period,
  durations = durations,
  cols = cols # deve seguir a ordem 1. gauge_code, 2. year, 3. duration, 4. imax
)

saveRDS(fit.ls, file.path(path.fit.ls, "df_dgev_subdiario.rds"))
