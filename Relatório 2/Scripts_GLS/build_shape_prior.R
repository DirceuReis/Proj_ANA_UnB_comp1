
rm(list = ls()); invisible(gc())
pacman::p_load(arrow)
source("scripts/funcoes/aux_functions.R")

# TRANSFORM EXPLANATORY VARIABLES -----------------------------------------

# Trying covariates without transformations
# Versão com variáveis dummy (z1, z2)
df.covariates <- arrow::read_parquet("base/gerados/covariates/df_covariates.pqt") %>% 
  arrange(gauge_code)
# df.covariates <- df.covariates %>% select(!z1)
# df.covariates <- df.covariates %>% select(!c(z1, z2))
gauges <- sort(unique(df.covariates$gauge_code))


# RUN OLS/WLS/GLS SCRIPT --------------------------------------------------

transformations <- c("none", "center", "scale", "log")
t <- 2 # generate final model only for centered variables

# purrr::walk(.x = seq_along(transformations), .f = function(t){
  
message(paste0(rep("_", getOption("width")), collapse = ""), "\n")
message("Running OLS/WLS/GLS script...\n")

save.path <- paste0("base/gerados/regression_results/", transformations[t], "/")

# Transform explanatory variables
covariates <- as.matrix(df.covariates[-1])
covariates <- covariate.transform(covariates, transformations[t], except = "z1")
# covariates <- covariate.transform(covariates, transformations[t])
covariates <- cbind("ones" = 1, covariates)
rownames(covariates) <- gauges

if(!dir.exists(save.path[t])) dir.create(save.path[t])


## BUILD LEAST SQUARES MODELS -------------------------------------------

# Run script to build LS models
source("scripts/ols_wls_gls_models.R")
source("scripts/funcoes/fun_regional_estimate.R") # call regional.est()


## ASSESS RESULTS -------------------------------------------------------

# Choose model with smallest AVP2 and estimate kappa and its variance
ls.regional.kappa <- regional.est(ls.regression = regression,
                                  covariates = covariates,
                                  conf.level = 0.95,
                                  select = "avp2",
                                  save = save.path)

# Save regression results
purrr::walk(.x = seq_along(ls.regional.kappa), .f = function(i){
  
  name <- unique(ls.regional.kappa[[i]]$regression)
  arrow::write_parquet(as.data.frame(ls.regional.kappa[[i]]), paste0(save.path, name, "_df_regional_kappa.pqt"))
  
})
  
# Choose regional average model and estimate kappa and its variance
ls.regional.average <- regional.est(ls.regression = regression,
                                    covariates = covariates,
                                    conf.level = 0.95,
                                    select = "average",
                                    save = save.path)

# Save regression results
purrr::walk(.x = seq_along(ls.regional.average), .f = function(i){

  name <- unique(ls.regional.average[[i]]$regression)
  arrow::write_parquet(as.data.frame(ls.regional.average[[i]]), paste0(save.path, "average_", name, "_df_regional_kappa.pqt"))

})
