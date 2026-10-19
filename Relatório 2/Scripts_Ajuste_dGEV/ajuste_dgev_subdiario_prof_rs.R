
# Os scripts atuais da d-GEV utilizam a parametrização da distribuição
# em que o parâmetro de forma negativo representa a distribuição de 
# cauda pesada (Fréchet), mas o nome dos parâmetros foi alterado para
# corresponder à nomenclatura que estamos utilizando.
#
# Este script ajusta a d-GEV com informação a priori vinda da distribuição
# Beta de Martins e Stedinger para o parâmetro de forma com média -0,10 e 
# variância 0,122². O ajuste é feito com Máxima Verossimilhança Generalizada
# e os intervalos de confiança são estimados com o método da verossimilhança
# perfilada.

rm(list = ls()); invisible(gc())

# PACOTES ----------------------------------------------------------------

library(dplyr)
library(tidyr)
library(pbapply)
library(lmom)
library(arrow)


# CAMINHOS ---------------------------------------------------------------

path.root <- "Relatório 2/Scripts_Ajuste_dGEV"
path.funs <- file.path(path.root, "funs", "profile_likelihood", "fitprof_dgev.R")
path.imax <- file.path(path.root, "resultados", "subdiarios", "dataframes", "df_imax_subdiario_ano_civil_rs.pqt")
path.fit.ls <- file.path(path.root, "resultados", "subdiarios", "fit", "profile_likelihood")


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
which.prior <- list(
  mle = NULL,               # s/ distribuicao a priori
  gmle = c(-0.1, 0.122),    # distribuicao a priori de Martins e Stedinger (2000)
  rgmle = c(0.0178, 0.0627) # distribuicao a priori regional c/ modelo GLS
) 
cols <- names(imax.df)
scale.inv <- "wide"     # estimativas iniciais dos parametros de invariancia simples
return.period <- c(2, 5, 10, 25, 50, 100, 200) # tempos de retorno
conf <- 0.95            # nivel de confianca p/ intervalos de confianca
maxit <- 1e8            # maximo de iteracoes p/ convergencia
step <- 0.01            # passo p/ encontrar as raizes da funcao de verossimilhanca perfilada
lltol <- 0.005          # tolerancia p/ diferenca entre maximas log-verossimilhancas

# Ajustar d-GEV, estimar quantis e intervalos de confiança com verossimilhança perfilada
imax.ls <- split(imax.df, imax.df$gauge_code); invisible(gc())

prior.idx <- 3
fit.ls <- fitprof.dgev(
  data = imax.ls,
  prior.info = which.prior[[prior.idx]],
  return.period = return.period,
  durations = durations,
  cols = cols, # deve seguir a ordem 1. gauge_code, 2. year, 3. duration, 4. imax
  step = step,
  profile.maxit = 500L
)

# Extrair e salvar resultados
if(!dir.exists(path.fit.ls)) dir.create(path.fit.ls)

# saveRDS(fit.ls, file.path(path.fit.ls, "ls_dgev_fit_subdiario.rds"))

ls <- lapply(list.files(path.fit.ls, full.names = TRUE, pattern = "*.rds"), readRDS)
par.df <- bind_rows(lapply(ls, `[[`, "par"))
idf.df <- bind_rows(lapply(ls, `[[`, "idf"))
profile.par <- bind_rows(lapply(ls, `[[`, "prof.par"))
profile.idf <- bind_rows(lapply(ls, `[[`, "prof.q"))

prefix <- "df_subdiario_dgev_rs"
prior.name <- names(which.prior)[prior.idx]

arrow::write_parquet(par.df, file.path(path.fit.ls, sprintf("%s_par_%s.pqt", prefix, prior.name))) # parametros
arrow::write_parquet(idf.df, file.path(path.fit.ls, sprintf("%s_idf_%s.pqt", prefix, prior.name))) # quantis
arrow::write_parquet(profile.par, file.path(path.fit.ls, sprintf("%s_profpar_%s.pqt", prefix, prior.name))) # verossimilhanca parametros
arrow::write_parquet(profile.idf, file.path(path.fit.ls, sprintf("%s_profidf_%s.pqt", prefix, prior.name))) # verossimilhanca quantis
