
# This function takes a list of OLS, WLS and GLS result matrices (with parameters, standard errors and performance metrics),
# assess the models' significance by applying a t-test with function model.sigf() and chooses between the regional average
# model or the best ranking (smallest) 'avp2' (AVPold) to generate regional interest variable (GEV's shape, or kappa here)
# means and variances, returning list of data.frames for each regression model (OLS, WLS and GLS).

regional.est <- function(ls.regression, covariates, conf.level = 0.95, select = c("average", "avp2"), save = NULL){
  
  # Check 'ls.regression' type
  if(!inherits(ls.regression, "list")) stop("\nObject 'ls.regression' must be a list with OLS, WLS and GLS regression matrices inside.")
  if(length(select) > 1) stop("\nChoose either 'average' or 'avp2'.")
  
  # Call in functions
  tryCatch(source("scripts/funcoes/aux_functions.R"),
           error = function(e){
             stop("\nSpecify the correct directory for auxiliary function 'model.sigf()'. Look for 'aux_functions.R' RScript.\n",
                  ">> Error: ", e$message)})
  
  # Run through regressions and choose either the best AVP2 or the regional average model
  ls.regional.kappa <- lapply(X = seq_along(ls.regression), FUN = function(i){
    
    message(paste0(rep("-", getOption("width")), collapse = ""))
    
    which.regression <- names(ls.regression)[[i]]
    df.models <- as.data.frame(ls.regression[[i]])
    df.models$id <- 1:nrow(df.models)
    
    # Test parameters for significance and save significant ones
    m <- nrow(covariates)
    kmax <- log(nrow(df.models), 2) + 1
    signf.models <- model.sigf(models = df.models, m = m, conf.level = conf.level)
    signf.models$model <- which.regression # include model name in the end
    
    if(!is.null(save)){
      if(!dir.exists(save)) dir.create(save)
      message("Saving significant models to: ", save)
      arrow::write_parquet(signf.models, paste0(save, which.regression, "_regression_sigf.pqt"))
    }
    
    # Select model
    message("Model chosen by: ", select)
    model <- switch(select,
                    "average" = signf.models[1,],
                    "avp2" = signf.models[which.min(signf.models[,select]),])
    
    k <- model[["k"]] # number of model parameters in the chosen model
    
    
    ### BUILD MODELS FOR KAPPA --------------------------------------------
    
    # Select covariates respective to the chosen model
    model.idx <- which(unlist(model)[1:kmax] != 0) # parameter columns
    X <- covariates[, model.idx, drop = FALSE]     # select appropriate c ovariates for chosen model
    var.delta <- model$var_delta                   # extract model variance
    
    message("Resulting model variance was: ", round(var.delta, 5))
    message("Estimating interest variable with ", toupper(which.regression),
            " regression model with the following covariates:\n>> ", paste(colnames(X), collapse = ", "))
    
    
    # OLS, WLS and GLS have different covariance matrices depending on the
    # Lambda matrix, which accounts for different sources of errors and correlation
    cov.par <- switch(which.regression,
                      "ols" = {
                        # Compute OLS covariance parameter matrix using only model error variance
                        # var.delta * solve(crossprod(X)) # using crossprod can be more effective
                        cov.par <- var.delta * solve(t(X) %*% X)
                      },
                      "wls" = {
                        # Compute WLS covariance parameter matrix using only the sampling covariance matrix's diagonal
                        wls.scov <- scov.matrix*diag(m)
                        inv.lambda <- build.gls(std.delta = var.delta^0.5, m, wls.scov, kappa, X, k)$inv.lambda
                        # solve(crossprod(X, inv.lambda %*% X))
                        cov.par <- solve(t(X) %*% inv.lambda %*% X) # parameter covariance matrix
                      },
                      "gls" = {
                        # Compute GLS covariance parameter matrix using the complete sampling covariance matrix
                        inv.lambda <- build.gls(std.delta = var.delta^0.5, m, scov.matrix, kappa, X, k)$inv.lambda
                        # solve(crossprod(X, inv.lambda %*% X)) 
                        cov.par <- solve(t(X) %*% inv.lambda %*% X) # parameter covariance matrix
                      })
    
    # Extract estimated parameters and estimate interest variable and variance
    beta <- as.numeric(unlist(model[model.idx]))
    mu.kappa.r <- as.vector(X %*% beta)
    sd.kappa.r <- (var.delta + rowSums((X %*% cov.par) * X))^0.5 # rowSums gets xi %*% cov.par %*% t(xi)
    
    # Results
    df.regional.kappa <- cbind(gauge_code = gauges,
                               regression = which.regression,
                               as.data.frame(X),
                               mu_kappa_r = mu.kappa.r,
                               sd_kappa_r = sd.kappa.r)
    
    return(df.regional.kappa)
    
  })
  
  return(ls.regional.kappa)
  
}
