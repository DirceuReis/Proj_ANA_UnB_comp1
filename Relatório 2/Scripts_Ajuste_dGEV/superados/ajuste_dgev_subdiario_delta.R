
# Os scripts atuais da d-GEV utilizam a parametrização da distribuição
# em que o parâmetro de forma negativo representa a distribuição de 
# cauda pesada (Fréchet), mas o nome dos parâmetros foi alterado para
# corresponder à nomenclatura que estamos utilizando
#
# Este script ajusta a d-GEV com informação a priori vinda da distribuição
# Beta de Martins e Stedinger para o parâmetro de forma com média -0,10 e 
# variância 0,122². O ajuste é feito com Máxima Verossimilhança Generalizada
# e os intervalos de confiança são estimados com o método Delta.

# PACOTES ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(pbapply)
library(purrr)
library(arrow)


# CAMINHOS ---------------------------------------------------------------

path.root <- "Relatório 2"
path.funs <- file.path(path.root, "Scripts_Ajuste_dGEV", "funs", "delta_method", "fitdelta_dgev.R")
path.imax <- file.path(path.root, "Scripts_Maximas", "resultados", "BR", "dataframes", "imax_sem_iqr_subdiario_Uniplu.parquet")
path.fit.ls <- file.path(path.root, "Scripts_Ajuste_dGEV", "resultados", "subdiarios", "fit")


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
  select(gauge_code, wateryear, d, imax) # essa ordem é importante p/ definir `cols` adiante

# AJUSTE d-GEV -----------------------------------------------------------

# Argumentos
durations <- NULL       # usa as durações disponíveis
prior <- c(-0.1, 0.122) # distribuicao a priori de Martins e Stedinger (2000)
cols <- names(imax.df)  # nome das colunas (gauge_code, year, duration, imax)
scale.inv <- "wide"     # estimativas iniciais dos parâmetros da invariância simples
method = "L-BFGS-B"
return.period <- c(2, 5, 10, 25, 50, 100, 200) # tempos de retorno
conf <- 0.95            # nível de confiança

# Separar data.frames em listas por estação
imax.ls <- split(imax.df, imax.df$gauge_code)

# Ajustar d-GEV
fit.ls <- fitdelta.dgev(
  data = imax.ls,
  prior.info = prior,
  return.period = return.period,
  durations = durations,
  cols = cols,
  method = method,
  conf = conf
)

par <- purrr::map_dfr(fit.ls, "par")
idf <- purrr::map_dfr(fit.ls, "idf")

out.path <- file.path(path.fit.ls, "delta_method")
if(!dir.exists(out.path)) dir.create(out.path)

arrow::write_parquet(par, file.path(path.fit.ls, "delta_method", "df_subdiario_dgev_par.pqt"))
arrow::write_parquet(idf, file.path(path.fit.ls, "delta_method", "df_subdiario_dgev_idf.pqt"))
