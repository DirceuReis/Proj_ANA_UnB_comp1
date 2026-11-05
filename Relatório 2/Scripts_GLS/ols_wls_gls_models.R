
# rm(list = ls()); invisible(gc())

# PACKAGES ----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow, psych, janitor)
# psych::tr() is used to calculate the trace of a matrix
# this can be replaced with R base's sum(diag())


# FUNCTIONS ---------------------------------------------------------------

source("scripts/funcoes/fun_scov_matrix.R") # build sampling covariance matrix for GLS
source("scripts/funcoes/aux_functions.R ")   # fun.ex.var and fun.gls.var.delta

# READ AND PREPARE DATA ---------------------------------------------------

# Daily rainfall annual maximum series
data <- arrow::read_parquet("base/gerados/df_daily_pmax.pqt") %>% 
  filter(!gauge_code %in% c("SBPA", "SBSM", "SBUG")) %>% 
  arrange(gauge_code)

col.names <- names(data)

# Build predictor variables matrix
gauges <- sort(unique(data[[col.names[1]]]))

# Apply transformations
# First column (gauge_code) becomes 'ones' column

# # Trying covariates without transformations
# df.covariates <- arrow::read_parquet("base/gerados/covariates/df_covariates.pqt") %>%  # covariates
#   arrange(gauge_code)
# covariates <- df.covariates[-1]
# covariates <- as.matrix(covariates) # for no transformation, keep only this row
# covariates <- covariate.transform(covariates, "center")
# covariates <- cbind("ones" = 1, covariates)
# rownames(covariates) <- gauges

# Build WLS sampling covariance matrix and get lmom kappa estimates
scov <- build.scovmatrix(df = data, col = c("gauge_code", "peak_24h"), which.corr = 2)
df.lmom.fit <- scov[[1]]   # extract lmom estimates for kappa other things
scov.matrix <- scov[[2]]   # extract wls sampling covariance matrix
kappa <- df.lmom.fit %>%   # extract lmom shape estimates
  arrange(gauge_code) %>% 
  pull(kappa)

names(kappa) <- gauges

stopifnot("Mismatch in covariate rownames" = identical(rownames(covariates), gauges))
stopifnot("Mismatch in covariate rownames" = identical(rownames(scov.matrix), gauges))
stopifnot("Mismatch in covariate rownames" = identical(colnames(scov.matrix), gauges))

# Build models
# model.names <- c("ones", "elevation", "long", "lat")                      # variable names
model.names <- colnames(covariates)
p <- length(model.names) - 1                                              # number of variables
models <- cbind(rep(1, 2^p), do.call(expand.grid, rep(list(c(0, 1)), p))) # build models matrix (1s and 0s)
colnames(models) <- model.names # add names
n.predictors <- ncol(models)    # number of parameters
n.models <- nrow(models)        # number of models (different combinations of covariates)

# Set column names for results
res.col.names <- c(paste0("b", 0:(n.predictors - 1)),   # model parameters
                   paste0("se", 0:(n.predictors - 1)),  # model parameters standard error
                   "var_delta",                         # model error variance
                   "asve", "avp1", "avp2", "pseudo_r2", "sst", "sse", "ssr", # evaluation metrics
                   "k")                                 # number of parameters

# Create empty result matrices
ols.result <- matrix(0, nrow = n.models, ncol = length(res.col.names), dimnames = list(NULL, res.col.names)) # OLS
wls.result <- matrix(0, nrow = n.models, ncol = length(res.col.names), dimnames = list(NULL, res.col.names)) # WLS
gls.result <- matrix(0, nrow = n.models, ncol = length(res.col.names), dimnames = list(NULL, res.col.names)) # WLS


# BUILD MODELS ------------------------------------------------------------

# Run through each combination of covariates and build OLS, WLS and GLS
message("Total number of models: ", n.models, "\n")
for(kmodel in 1:n.models){
  
  # Numeric vector indicating which model to use
  model.idx <- which(models[kmodel,] > 0)
  
  # Extract predictor variables matrix X for the current model
  X <- covariates[, model.idx, drop = FALSE]
  m <- nrow(X) # number of gauges (changed from n)
  k <- ncol(X) # number of n=model parameters
  p <- k - 1   # number of covariates
  
  message("[", kmodel, "]    Covariates: ", paste(colnames(X), collapse = ", "))

  # # Pre-define dimensions for output parameters (b) and standard errors
  ols.res.par <- matrix(0, nrow = ncol(models), ncol = 1)    # OLS model parameters (beta.ols)
  wls.res.par <- matrix(0, nrow = ncol(models), ncol = 1)    # WLS model parameters (beta.wls)
  gls.res.par <- matrix(0, nrow = ncol(models), ncol = 1)    # GLS model parameters (beta.gls)
  ols.res.se.par <- matrix(0, nrow = ncol(models), ncol = 1) # OLS parameters standard error (sigma.delta.ols)
  wls.res.se.par <- matrix(0, nrow = ncol(models), ncol = 1) # WLS parameters standard error (sigma.delta.wls)
  gls.res.se.par <- matrix(0, nrow = ncol(models), ncol = 1) # GLS parameters standard error (sigma.delta.gls)
  
  Xt <- t(X) # X transposed
  

  # OLS --------------------------------------------------------------------
  
  # Compute parameter, variance, covariance and standard error
  ols.par <- solve(Xt %*% X) %*% Xt %*% kappa                             # matriz com parametros do modelo
  ols.residual <- kappa - X %*% ols.par                                   # aux residuals matrix estimates - model
  ols.var.delta <- as.numeric((t(ols.residual) %*% ols.residual)/(m - k)) # model error variance (antes tava m - p (errado))
  ols.cov.par <- ols.var.delta * solve(Xt %*% X)                          # parameter covariance matrix
  ols.se.par <- sqrt(diag(ols.cov.par))                                   # standard error
  
  # Performance metrics
  ols.diag.var <- ols.var.delta*diag(m)                      # model error variance * identity matrix
  ols.inv.lambda <- solve(ols.diag.var)                      # model error variance matrix Λ-1
  ols.sve <- X %*% solve(Xt %*% ols.inv.lambda %*% X) %*% Xt # sampling variance
  ols.sve <- ols.var.delta * (X %*% solve(Xt %*% X) %*% Xt)  # SVE = sigma.var * X * (X'X)^-1 * X'
  ols.asve <- sum(diag(ols.sve))/m                           # average sampling variance for new sites
  # ols.asve <- (ols.var.delta * k)/m
  ols.avp1 <- ols.var.delta + ols.asve                       # average variance of prediction
  term <- 2 * ols.var.delta * ols.sve %*% ols.inv.lambda
  ols.avp2 <- ols.avp1 - sum(diag(term))/m                   # average variance of prediction for old sites
  # ols.avp2 <- ols.var.delta - ols.asve
  
  # Pseudo-R^2
  ols.sse <- t(ols.residual) %*% ols.residual             # residual sum of squared errors
  mean.kappa <- as.numeric(apply(X, 2, mean) %*% ols.par) # average predicted kappa
  if(kmodel == 1) ols.sst <- t(kappa - mean.kappa) %*% (kappa - mean.kappa)          # total sum of squared deviation
  ols.ssr <- as.numeric(t(mean.kappa - X %*% ols.par) %*% (mean.kappa - X %*% ols.par)) # sum of squared erros from regression
  ols.pr2 <- 1 - ols.sse/(ols.sse + ols.ssr)              # pseudo-R^2
  

  # WLS --------------------------------------------------------------------
  
  # Extract sampling covariance matrix for WLS
  wls.scov <- scov.matrix*diag(m) # only diagonal
  
  # Check if var.delta = 0 is enough, if not, optimize for check.delta > 0
  sd.delta <- 0
  wls.est <- build.gls(sd.delta, m, wls.scov, kappa, X, k)
  check.delta <- wls.est$check.delta
  sd.delta <- wls.est$std.delta
  wls.par <- wls.est$b
  wls.lambda <- wls.est$lambda
  wls.inv.lambda <- wls.est$inv.lambda
  
  if(check.delta < 0){
    
    sd.delta <- 0
    
    } else{
    
      # Define check_delta as a function of just sd_delta (g)
      f.root <- function(g){build.gls(g, m, wls.scov, kappa, X, k)$check.delta}
      wls.est.root <- uniroot(f.root, c(0, 20)) # find root of objective function f_root
      sd.delta <- wls.est.root$root             # extract root -> optimized sd_delta
      
      # Run parameter function again for new sd_delta
      wls.est <- build.gls(sd.delta, m, wls.scov, kappa, X, k)
      check.delta <- wls.est$check.delta
      sd.delta <- wls.est$std.delta
      wls.par <- wls.est$b
      wls.lambda <- wls.est$lambda
      wls.inv.lambda <- wls.est$inv.lambda
      
    }
  
  # Covariance and standard error for WLS parameters
  wls.cov.par <- solve(Xt %*% wls.inv.lambda %*% X)
  wls.se.par <- sqrt(diag(wls.cov.par))
  wls.var.delta <- sd.delta^2
  
  # Performance metrics
  wls.sve <- X %*% wls.cov.par %*% Xt
  wls.asve <- tr(wls.sve)/m
  wls.avp1 <- wls.var.delta + wls.asve
  term <- 2 * wls.var.delta * wls.sve %*% wls.inv.lambda
  wls.avp2 <- wls.avp1 - tr(term)/m
  
  # Pseudo-R^2
  wls.sse <- m*wls.var.delta
  mean.kappa <- as.numeric(apply(X, 2, mean) %*% wls.par)
  if(kmodel == 1) wls.sst <- t(kappa - mean.kappa) %*% (kappa - mean.kappa)
  wls.ssr <- as.numeric(t(mean.kappa - X %*% wls.par) %*% (mean.kappa - X %*% wls.par))
  wls.pr2 <- 1 - wls.sse/(wls.sse + wls.ssr)
  

  # GLS --------------------------------------------------------------------
  
  # Check if var.delta = 0 is enough, if not, optimize for check.delta > 0
  sd.delta <- 0
  gls.est <- build.gls(sd.delta, m, scov.matrix, kappa, X, k)
  check.delta <- gls.est$check.delta
  sd.delta <- gls.est$std.delta
  gls.par <- gls.est$b
  gls.lambda <- gls.est$lambda
  gls.inv.lambda <- gls.est$inv.lambda
  
  if(check.delta < 0){
    
    sd.delta <- 0
    
  } else{
    
    # Define check_delta as a function of just sd_delta (g)
    f.root <- function(g){build.gls(g, m, scov.matrix, kappa, X, k)$check.delta}
    gls.est.root <- uniroot(f.root, c(0, 20)) # find root of objective function f_root
    sd.delta <- gls.est.root$root             # extract root -> optimized sd_delta
    
    # Run parameter function again for new sd_delta
    gls.est <- build.gls(sd.delta, m, scov.matrix, kappa, X, k)
    check.delta <- gls.est$check.delta
    sd.delta <- gls.est$std.delta
    gls.par <- gls.est$b
    gls.lambda <- gls.est$lambda
    gls.inv.lambda <- gls.est$inv.lambda
    
  }
  
  # Covariance and standard error for gls parameters
  gls.cov.par <- solve(Xt %*% gls.inv.lambda %*% X)
  gls.se.par <- sqrt(diag(gls.cov.par))
  gls.var.delta <- sd.delta^2
  
  # Performance metrics
  gls.sve <- X %*% gls.cov.par %*% Xt
  gls.asve <- tr(gls.sve)/m
  gls.avp1 <- gls.var.delta + gls.asve
  term <- 2 * gls.var.delta * gls.sve %*% gls.inv.lambda
  gls.avp2 <- gls.avp1 - tr(term)/m
  
  # Pseudo-R^2
  gls.sse <- m*gls.var.delta
  mean.kappa <- as.numeric(apply(X, 2, mean) %*% gls.par)
  if(kmodel == 1) gls.sst <- t(kappa - mean.kappa) %*% (kappa - mean.kappa)
  gls.ssr <- as.numeric(t(mean.kappa - X %*% gls.par) %*% (mean.kappa - X %*% gls.par))
  gls.pr2 <- 1 - gls.sse/(gls.sse + gls.ssr)
  

  # RESULTS ----------------------------------------------------------------
  
  # OLS
  ols.res.par[model.idx] <- ols.par
  ols.res.se.par[model.idx] <- ols.se.par
  ols.result[kmodel,] <- c(t(ols.res.par), t(ols.res.se.par),
                           ols.var.delta, ols.asve, ols.avp1, ols.avp2,
                           ols.pr2, ols.sst, ols.sse, ols.ssr, k)
  
  # WLS
  wls.res.par[model.idx] <- wls.par
  wls.res.se.par[model.idx] <- wls.se.par
  wls.result[kmodel,] <- c(t(wls.res.par), t(wls.res.se.par),
                           wls.var.delta, wls.asve, wls.avp1, wls.avp2,
                           wls.pr2, wls.sst, wls.sse, wls.ssr, k)
  
  # GLS
  gls.res.par[model.idx] <- gls.par
  gls.res.se.par[model.idx] <- gls.se.par
  gls.result[kmodel,] <- c(t(gls.res.par), t(gls.res.se.par),
                           gls.var.delta, gls.asve, gls.avp1, gls.avp2,
                           gls.pr2, gls.sst, gls.sse, gls.ssr, k)
  
}

regression <- list(ols = ols.result, wls = wls.result, gls = gls.result)

# Save results
purrr::walk(.x = seq_along(regression), .f = function(i){

  name <- names(regression)[i]
  arrow::write_parquet(as_tibble(regression[[i]]), paste0(save.path, name, "_regression.pqt"))

})

# Clear environment
keep.objects <- c("models", "ols.result", "wls.result", "gls.result", "scov.matrix", "data", "transformations",
                  "df.covariates", "covariates", "regression", "p" , "m", "gauges", "kappa", "model.names", "save.path")
keep.functions <- c("model.sigf", "build.gls", "covariate.transform")
rm(list = setdiff(ls(), c(keep.objects, keep.functions)))
invisible(gc())

message("OLS, WLS and GLS models built")
