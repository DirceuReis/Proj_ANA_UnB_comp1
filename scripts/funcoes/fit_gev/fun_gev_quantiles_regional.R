# Funciton for computing GEV quantiles by selecting different models and
# giving a vector of return periods to calculate
# This function also returns quantile variances for quantiles when the method
# selected involves maximum likelihood estimators
fun.gev.quantile.regional <- function(df,
                                      tempo.retorno,
                                      gauge.names = "gauge_code",
                                      signf.level = 0.05,
                                      model = c("lmom.gu", "lmom.gev", "l.gu", "l.gev", "gl.gev", "gl.gev.r")){
  
  # Check col name
  if(is.null(df[[gauge.names]]) == TRUE){
    stop("None of the column names in 'df' are ", gauge.names)
  }
  
  # Timer
  start.time <- Sys.time()
  
  # Initial information
  estacoes <- df[[gauge.names]]
  n.estacoes <- length(estacoes)
  p <- 1 - 1/tempo.retorno
  n.p <- length(p)
  
  # Start resulting table
  df.q <- tibble()
  
  # Loop through stations
  for(i in seq_along(estacoes)){
    
    est <- estacoes[i]
    
    # cat("\nEstação", est, "em processamento...")
    
    # Set current station
    df.current <- df[df[[gauge.names]] == est,]
    series.length <- df.current$n.serie
    
    # Loop through models to determine which should be calculated
    for(mod in model){
      
      # Extract parameter vector for the current model
      par.col <- paste0("par.", mod)
      
      # Check if the column exists
      if(!par.col %in% names(df.current)){
        message(sprintf("Station %s: column %s not found. Skipping model %s.", est, par.col, mod))
        next
      }
      
      # Unlist parameters if the column exists and extract them
      par <- unlist(df.current[[par.col]])
      n.par <- length(par) # number of parameters for the current distribution
      
      # Select appropriate quantile function and extract parameters
      if(grepl("gu", mod)){ # check for the distribution suffix in the model name
        
        # Extract parameters
        csi <- par[1]
        alpha <- par[2]
        
        # Quantiles
        q <- fun.q.gu(p, param = par)
        
      } else {               # if it's not Gumbel, then use the GEV quantile function
        
        # Extract parameters
        csi <- par[1]
        alpha <- par[2]
        kappa <- par[3]
        
        # Quantiles
        q <- fun.q.gev(p, param = par)
        
      }
      
      # Estimate uncertainty using the Delta Method -> only for Maximum Likelihood models
      # Set NA for lmom models
      if(mod %in% c("l.gu", "l.gev", "gl.gev", "gl.gev.r")){
        
        # Extract hessian
        hess.col <- paste0("hess.", mod)
        if(!hess.col %in% names(df.current)){
          
          message(sprintf("Station %s: hessian column %s not found. Skipping uncertainty estimation for model %s.", est, hess.col, mod))
          
          # Set columns to NA
          sd.q <- rep(NA, n.p)
          ic.l <- rep(NA, n.p)
          ic.u <- rep(NA, n.p)
          
        } else {
          
          # If the hessian column exists, then estimate quantile variance
          # Extract and invert hessian matrix
          hess <- df.current[[hess.col]] %>% unlist %>% matrix(ncol = n.par, nrow = n.par, byrow = TRUE) # extract hessian for current model
          cov.q <- tryCatch(
            solve(hess),
            error = function(e){
              message(sprintf("\nStation %s: hessian is null as the station was not optimized in model %s → %s", est, mod, e$message))
              matrix(NA, nrow = length(par), ncol = length(par))
            })
          
          # Calculate variance via the Delta Method
          var.q <- rep(NA, n.p) # create empty vector
          for(prob in seq_along(p)){
            
            # Structure gradient matrix depending on the distribution
            if(grepl("gev", mod)){
              
              d.csi <- 1                                                                # dq/dξ
              d.alpha <- 1/kappa*(1 - (-log(p[prob]))^kappa)                            # dq/dα
              yp <- -log(p[prob])                                                       # to simplify dq/dκ expression
              d.kappa <- -alpha/(kappa^2)*(1 - yp^kappa) - alpha/kappa*yp^kappa*log(yp) # dq/dκ
              
              grad.q <- c(d.csi, d.alpha, d.kappa)
              
            }
            
            if(grepl("gu", mod)){
              
              d.csi <- 1                     # dq/dξ
              d.alpha <- -log(-log(p[prob])) # dq/dα
              
              grad.q <- c(d.csi, d.alpha)
              
            }
            
            # Variance
            var.q[prob] <- t(grad.q) %*% cov.q %*% grad.q
            
          }
          
          # Estimate confidence intervals
          sd.q <- sqrt(var.q)                  # quantile standard deviation
          delta.ci <- qnorm(1 - signf.level/2) # standard normal quantile for given confidence level
          ic.l <- q - delta.ci*sd.q            # lower limit of CI
          ic.u <- q + delta.ci*sd.q            # upper limit of CI
          
        }
        
      } else {
        
        # For lmom models: uncertainty values are not estimated
        sd.q <- rep(NA, n.p)
        ic.l <- rep(NA, n.p)
        ic.u <- rep(NA, n.p)
        
      }
      
      # Organize results for current station and model
      df.q.aux <- tibble(gauge.code = rep(est, n.p),
                         n.serie = rep(series.length, n.p),
                         model = mod,
                         tempo.retorno = tempo.retorno,
                         q = q,
                         sd.q = sd.q,
                         ic.l = ic.l,
                         ic.u = ic.u)
      
      # Combine results
      df.q <- rbind(df.q, df.q.aux)
      
    } # end model loop
    
  } # end station loop
  
  # Processing time and message
  procss.time <- round(difftime(Sys.time(), start.time, units = "mins"), 1)
  message(sprintf("\nProcessing complete!\nStations processed: %d\nDuration: %s min", n.estacoes, procss.time))
  gc()             # clean memory
  
  return(df.q)
  
}
