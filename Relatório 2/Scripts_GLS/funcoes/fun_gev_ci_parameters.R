# This function is implemented after using fun.gev.fit.regional() or fun.process() together
# with fun.ajuste.pmax() to estimate parameter variance in the fitting process based on the
# inverse hessinan matrix from the optimization process of the maximum likelihood methods

# The function returns variance and lower and upper limits for 95% (or other specified)
# confidence interval for all three parameters of the GEV distribution

conf.gev <- function(df,                                                # df.fit from fun.gev.fit.regional()
                     signf.level = 0.05,                                # for CI estimation
                     gauge.names = "gauge.code",                        # name of the column from df.fit with station names
                     par.preffix = "par.",                              # preffix which identifies columns with parameters
                     hess.preffix = "hess.",                            # preffix which identifies columns with hessian matrix
                     models = c("l.gu", "l.gev", "gl.gev", "gl.gev.r"), # maximum likelihood "models" from which to choose from
                     pars = c("csi", "alpha", "kappa")){                # which paremeters to estimate for
  
  # Create result storage
  results <- tibble()
  
  # Checking if models exist in 'df'
  par.cols <- paste0(par.preffix, models)
  hess.cols <- paste0(hess.preffix, models)
  if(any(par.cols %in% names(df)) == FALSE | any(hess.cols %in% names(df) == FALSE)){
    stop("One or more parameter or hessian columns could not be found in 'df'.\n", "names(df):\n", names(df),
         "Looking for: ", par.cols, "\nand ", hess.cols)
  }
  
  # Set indexes for getting parameters
  par.idx <- 1:3
  names(par.idx) <- c("csi", "alpha", "kappa")
  
  tryCatch({
    
    gauges <- df[[gauge.names]]
    
    # Go through each df.fit, extract parameters, hessian, variances and estimate CI
    for(i in seq_along(gauges)){
      
      # Current station
      gg <- df[[gauge.names]][i]            # current station name
      df.gg <- df[df[[gauge.names]] == gg,] # current row
      
      # Extract parameters
      for(mod in models){
        
        # Parameter and names
        par.col <- paste0(par.preffix, mod) # current parameter column name
        par <- df.gg[[par.col]][[1]]     # unlist parameter data.frame
        
        # Extract hessian
        hess.col <- paste0(hess.preffix, mod)     # current hessian column name
        hess <- as.matrix(df.gg[[hess.col]][[1]]) # extract hessian from listed cell
        if(any(dim(hess)) == 0) next              # skip if its empty
        if(mod == "l.gu") hess <- hess[1:2, 1:2]  # hessianas de l.gu sairam 3x3 (conferir depois)
        cov <- try(solve(hess), silent = TRUE)    # invert hessian -> covariance or inverse Fisher information matrix
        if(inherits(cov, "try-error")) next       # skip if error when solving
        
        # Define current variance output vector based on n.par
        n.par <- length(names(par)) # number of parameters in the given model
        sd.par <- rep(NA, n.par)    # vector of standard deviations
        ci.l <- rep(NA, n.par)      # vector of lower CI limits
        ci.u <- ci.l                # vector of upper CI limits
        
        # Extract parameters
        for(j in seq_along(names(par))){
          
          p <- names(par)[j]            # parameter name
          idx <- par.idx[[p]]           # matrix position to look for based on the current par
          var.par <- cov[[idx,idx]]     # get variance from diagonal and calculate standard deviation
          sd.par[j] <- var.par^0.5
          z <- 1 - signf.level/2        # CI percentile
          aux.ci <- qnorm(z)            # standard normal quantile
          ci.l[j] <- par[[p]] - aux.ci*sd.par[j] # lower limit
          ci.u[j] <- par[[p]] + aux.ci*sd.par[j] # upper limit
          
        }
        
        res.aux <- tibble(gauge.code = gg,
                          model = mod,
                          par.name = names(par),
                          par = unlist(par),
                          sd = sd.par,
                          ci.l = ci.l,
                          ci.u = ci.u)
        
        results <- bind_rows(results, res.aux)
        results <- results[results[["par.name"]] %in% pars,]
        
      }
      
    }
    
  }, error = function(e){
    
    stop("Processing stopped because of:\n", e)
    
  })
  
  return(results)
  
}