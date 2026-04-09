# This function fits the GEV distribution to precipitation timeseries
# Arguments are data, signf_level, regional_models, data_names and model_names:
# - data: tbl_df with AM series of a set of rain gauges having columns for station code and annual maximum;
# - regional_models: tbl_df with with regional information on kappa having columns for station code and 
#   the mean and standard deviation of kappa for each station to build the prior;
# - data_names: column names of the data tbl_df;
# - model_names: column names of the models tbl_df;
# - signf_level: significance levele

# Função com pblapply
# Para o processamento quando encontra um erro
fun.gev.fit.regional <- function(df,
                                 regional.models,
                                 data.names = c("gauge_code", "rain_mm"),
                                 model.names = c("gauge_code", "mu_kappa", "sd_kappa"),
                                 signf.level = 0.05){
  
  # Packages
  if(!require("pacman")) install.packages("pacman"); pacman::p_load(pacman, tidyverse, lmom, pbapply, beepr)
  
  # Check if df is tbl_df
  if(!inherits(df, c("data.frame", "data.table", "tibble"))){
    stop("Argument 'df' must be a data.frame (data.table or tibble), not a ", class(df), ".")
  }
  
  # Check if regional.models is tbl_df
  if(!inherits(regional.models, c("data.frame", "data.table", "tibble"))){
    stop("Argument 'regional.models' must be a data.frame (data.table or tibble), not a ", class(regional.models), ".")
  }
  
  # Function for fitting the GEV distribution via L-MOM, MLE and GMLE with geophysical and regional prior on kappa
  fun.fit <- function(df){
    
    tryCatch({
      
      # Extract basic information from single station
      gg <- df[[data.names[1]]][1]
      mu.kappa <- regional.models[regional.models[[model.names[1]]] == gg, model.names[2]]
      sd.kappa <- regional.models[regional.models[[model.names[1]]] == gg, model.names[3]]
      pmax <- df[[data.names[2]]]
      n <- length(pmax)
      p0 <- sum(pmax == 0)/n
      pmax.pos <- pmax[pmax > 0]
      
      # Fit GEV with L-MOM
      par.lmom.gev <- fun.lmom(xp = pmax.pos, dist = "gev")
      par.lmom.gu <- fun.lmom(xp = pmax.pos, dist = "gu")
      
      # Fit GEV with ML using L-MOM as starting point
      par.l.gev <- fun.max.l(param = unlist(par.lmom.gev), xp = pmax.pos, dist = "gev")
      par.l.gu <- fun.max.l(param = unlist(par.lmom.gu), xp = pmax.pos, dist = "gumbel")
      
      # Fit GEV with GMLE and geophysical prior (see Martins and Stedinger, 2000) and regional prior
      par.gl.gev <- fun.max.gl(param = unlist(par.lmom.gev), xp = pmax.pos, media.kappa = -0.10, desvpad.kappa = 0.122)
      par.gl.gev.r <- fun.max.gl(param = unlist(par.lmom.gev), xp = pmax.pos, media.kappa = mu.kappa, desvpad.kappa = sd.kappa)
      
      # Extract function maxima
      max.l.gu <- par.l.gu$value
      max.l.gev <- par.l.gev$value
      max.gl.gev <- par.gl.gev$value
      max.gl.gev.r <- par.gl.gev.r$value
      
      # Hessian matrix
      hess.l.gu <- par.l.gu$hess
      hess.l.gev <- par.l.gev$hess
      hess.gl.gev <- par.gl.gev$hess
      hess.gl.gev.r <- par.gl.gev.r$hess
      
      # Results
      res <- tibble(gauge.code = gg,
                    rain.mm = list(tibble(pmax = pmax)),
                    p0 = p0,
                    n.serie = n,
                    par.lmom.gev = list(tibble(par.lmom.gev)),
                    par.lmom.gu = list(tibble(par.lmom.gu)),
                    par.l.gev = list(tibble(par.l.gev[,-c(4,5)])),       # removing cols 'value' e 'hess'
                    par.l.gu = list(tibble(par.l.gu[,-c(3,4)])),         # removing cols 'value' e 'hess'
                    par.gl.gev = list(tibble(par.gl.gev[,-c(4,5)])),     # removing cols 'value' e 'hess'
                    par.gl.gev.r = list(tibble(par.gl.gev.r[,-c(4,5)])), # removing cols 'value' e 'hess'
                    max.l.gu = max.l.gu,
                    max.l.gev = max.l.gev,
                    max.gl.gev = max.gl.gev,
                    max.gl.gev.r = max.gl.gev.r,
                    hess.l.gu = hess.l.gu,
                    hess.l.gev = hess.l.gev,
                    hess.gl.gev = hess.gl.gev,
                    hess.gl.gev.r = hess.gl.gev.r,
                    .rows = 1)
      
      return(res)
      
    }, error = function(e){
      
      gg <- df[[data.names[1]]][1]
      message("\nSkipping station ", gg, " due to error: ", e$message)
      NULL
      
    })
    
  }
  
  # Apply function
  list.df <- split(df, df[[data.names[1]]])
  list.models <- split(regional.models, regional.models[[model.names[1]]])
  list.fitted <- pbapply::pblapply(X = list.df, FUN = fun.fit)
  df.fit <- bind_rows(list.fitted)
  beepr::beep(sound = 10)
  gc()
  
  return(df.fit)
  
}