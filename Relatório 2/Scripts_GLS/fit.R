
# PACKAGES ----------------------------------------------------------------

if(!require("pacman")) install.packages("pacman"); pacman::p_load(pacman, tidyverse)


# READ DATA ---------------------------------------------------------------

# Clear environment
rm(list = ls()); invisible(gc())

# Read data
data <- arrow::read_parquet("base/gerados/df_daily_pmax.pqt")  # daily rainfall data
data <- data[!data$gauge_code %in% c("SBPA", "SBUG", "SBSM"),]; message(n_distinct(data$gauge_code), " gauges")

# Functions
source("scripts/funcoes/aux_functions.R")              # auxilliary functions for fitting the GEV
source("scripts/funcoes/fun_gev_fit_regional.R")       # fit GEV with LMOM, MLE and GMLEs
source("scripts/funcoes/fun_gev_ci_parameters.R")      # parameter CIs
source("scripts/funcoes/fun_gev_quantiles_regional.R") # GEV quantiles

transform.path <- list.dirs("base/gerados/regression_results/", recursive = FALSE)
# t <- 1 # transformation chosen (center), now there is only one
# transform.path <- paste0(transform.path[t], "/")
transform.path <- "base/gerados/regression_results/center/"

save.path <- paste0(transform.path, "fit/")
if(!dir.exists(save.path)) dir.create(save.path)
message("Directory: ", save.path)

# DISTRIBUTION FITTING AND QUANTILE ESTIMATION ----------------------------

# Return levels
return.lvl <- c(10.0, 20.0, 50.0, 100.0, 200.0, 500.0)

# Read files from appropriate transformation directory
files.regional <- list.files(transform.path, pattern = "_df_regional_kappa.pqt", full.names = TRUE)
names <- gsub("_df_regional_kappa\\.pqt", "", basename(files.regional))
ls.regional.kappa <- lapply(files.regional, arrow::read_parquet); names(ls.regional.kappa) <- names

purrr::walk(.x = seq_along(ls.regional.kappa), .progress = TRUE, .f = function(i){
  
  message(paste0(rep("-", getOption("width")), collapse = ""))
  
  name <- names(ls.regional.kappa[i])
  df.regional.kappa <- ls.regional.kappa[[i]]
  which.regression <- df.regional.kappa$regression[1]
  cat("Distribution fitting and quantile estimation using regional ",
          toupper(which.regression), " model\n>> ", name)
  
  # FIT GEV -------------------------------------------------------------
  
  # Fit GEV via L-MOM, MLE and GMLE with geophysical prior and regional prior
  df.fit <- fit.gev(df = data,
                    regional.models = df.regional.kappa,
                    data.names = c("gauge_code", "peak_24h"),
                    model.names = c("gauge_code", "mu_kappa_r", "sd_kappa_r"))
  
  # There still needs to be a step here to compute kappa variances using the hessian
  df.kappa.ci <- conf.gev(df = df.fit,
                          gauge.names = "gauge.code",
                          models = c("l.gev", "gl.gev", "gl.gev.r"),
                          pars = "kappa")
  
  
  # QUANTILES -----------------------------------------------------------
  
  df.quantiles <- quantile.gev(df = df.fit,
                               gauge.names = "gauge.code",
                               tempo.retorno = return.lvl,
                               signf.level = 0.05,
                               model = c("l.gev", "gl.gev", "gl.gev.r"))
  
  
  # SAVE RESULTS --------------------------------------------------------
  
  saveRDS(df.fit, paste0(save.path, name, "_df_fit.rds"))
  saveRDS(df.quantiles, paste0(save.path, name, "_df_quantiles.rds"))
  saveRDS(df.kappa.ci, paste0(save.path, name,"_df_kappa_ci.rds"))
  
})

message("Saving 'fit' results to: ", save.path)
