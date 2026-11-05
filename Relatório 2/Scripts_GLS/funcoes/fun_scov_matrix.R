
# BUILD SAMPLING COVARIANCE MATRIX FOR WLS REGRESSION ---------------------

# This function builds a sampling covariance matrix for GEVs shape parameter
# for the WLS regression model. This means the function below returns a m x m
# matrix filled only with Var(kappa) in its main diagonal.

# An estimate of kappa that is almost independent was obtained by the method
# described by Chowdhury et al. (1991), Fill [Cornell PhD Thesis] (1994) and
# Madesn and Rosjberg (1997)

# Arguments of the function are:
# df: a data.frame containing Annual Maxima timeseries for one or more
# stations, with at least two columns: gauge_code and rain_mm, different
# names may be passed through the col argument;
# col: vector with two characters, indicating the names of the columns 
# containing in col[1] station names and in col[2] the hydrological variable
# of interest
# which.corr: which correlation model from the 'rhoij' list to use

build.scovmatrix <- function(df, col = c("gauge_code", "rain_mm"), which.corr){
  
  # Corrrelation models to estimate corr(kappa_i,kappa_j)
  source("scripts/correlation_model.R")
  message("Building correlation models and choosing model: ", which.corr)
  message("Building GLS sampling covariate matrix...")
  
  # Packages
  if(!require(pacman)){install.packages("pacman")}
  pacman::p_load(pacman, tidyverse, lmom, gsl)
  
  # Preliminary steps
  gauges <- sort(unique(df[[col[1]]])) # get unique station names
  m <- length(gauges)
  kappas <- rep(0, m)
  tau3 <- rep(0, m)
  n <- rep(0, m)
  
  # Get L-moments for all stations and fit GEV to extract kappa
  for(gg in seq_along(gauges)){
    
    # Get AM series for current gauge
    cod <- gauges[[gg]]
    curr.gg <- df[df[col[1]] == cod, col[2], drop = TRUE]
    n[gg] <- length(curr.gg)
    
    # Compute sample L-moments and ratios
    lmom <- lmom::samlmu(x = curr.gg, sort.data = TRUE) # sample l-moments
    tau3[gg] <- lmom[[3]]                               # l-skewness (t3)
    kappas[gg] <- lmom::pelgev(lmom)[[3]]               # fit GEV and get kappa
    
    # ADICIONAR:tau4 (curtose), precipitações máxima média de 1 e 3 dias
    
  }
  
  # Clear memory
  gc()
  
  # Average tau3 and kappa
  tau3.a <- mean(tau3)                          # average tau3
  c = 2/(3 + tau3.a) - log(2)/log(3)
  k <- 7.859*c + 2.9554*c^2                     # average kappa
  n.a <- mean(n)                                # average number of years
  
  # Tests Table 1 – Chowdhury et al. (1991)
  # k <- -0.2                                   # average kappa
  # tau3_a <- (2*(1 - 3^(-k)))/(1 - 2^(-k)) - 3 # average tau3 for testing
  # c = 2/(3 + tau3_a) - log(2)/log(3)
  # n_a <- 500                                  # average number of years for testing
  
  
  # L-SKEWNESS AND KAPPA VARIANCES ------------------------------------------
  
  # Compute var(tau3) and var(kappa) – Chowdhury et al. (1991) and Fill [Cornell PhD Thesis] (1994)
  r <- gamma(1 + 2*k)/(gamma(1 + k)^2) # R(k) function
  aux <- (1 - 2^(-k))^2                # auxiliary denominator that shows up a lot
  
  # Using the gsl::hyperg_2F1() function as MatLabs 'mhygfx'
  hy05 <- gsl::hyperg_2F1(k, 2*k, 1 + k, -0.5)
  hy13 <- gsl::hyperg_2F1(k, 2*k, 1 + k, -1/3)
  hy23 <- gsl::hyperg_2F1(k, 2*k, 1 + k, -2/3)
  
  # Elements fij from the symetric variance matrix for l-moment ratios T
  f11 <- (r - 1)/aux
  f22 <- (r*hy05 - 1)/(2^(2*k)*aux)
  f33 <- (r*hy23 - 1)/(3^(2*k)*aux)
  f12 <- 0.5*(r - 2^(1 + k) + 2^(2*k))/(2^(2*k)*aux)
  f13 <- 0.5*((r - 2*3^k)/(3^(2*k)*aux) - (r*hy05 - 2^(1 + k))/(2^(2*k)*aux))
  f23 <- 0.5*((r*hy13 - 2*(3/2)^k)/(3^(2*k)*aux) + 1/(2^(2*k)*aux))
  
  # Variance of tau3 and kappa
  var.tau3.a <- 1/n.a*((f11 - 4*f12 + 4*f22)*tau3.a^2 + 2*(f11 - 8*f12 + 12*f22 + 6*f13 - 12*f23)*tau3.a + f11 - 12*f12 + 36*f22 + 12*f13 - 72*f23 + 36*f33)
  dk.dtau3.a <- -(2.872/(tau3.a + 3))^2*(2.872/(tau3.a + 3) + 1) # partial kappa in relation to tau3
  var.tau3 <- (n.a/n)*var.tau3.a                                 # vector of tau3 variances
  var.k.a <- dk.dtau3.a^2*var.tau3.a                             # average variance of kappa
  var.k <- (n.a/n)*var.k.a                                       # variance of kappa for each station (vector)
  na.var.tau3 <- n.a*var.tau3.a
  
  # Tests Table 1 – Chowdhury et al. (1991) ----
  # Run below
  # k <- -0.2                                   # average kappa
  # tau3.a <- (2*(1 - 3^(-k)))/(1 - 2^(-k)) - 3 # average tau3 for testing
  # c = 2/(3 + tau3.a) - log(2)/log(3)
  # n.a <- 500                                  # average number of years for testing
  # # Then run 'var.tau.a' and then:
  # n.a*var.tau3.a # should return "0.4607" like in Table 1 for k = -0.2
  
  
  # BUILD GLS SCOV MATRIX ---------------------------------------------------
  
  # Interpolate mean kappa and get Delta power exponent
  # Martins and Stedinger's (2001) power exponents for cross-correlation
  pred.d <- approxfun(x = c(-0.3, -0.2, -0.1, 0, 0.1), y = c(2.4, 2.6, 2.8, 2.9, 3.1))
  d <- pred.d(mean(kappas))
  
  # Calculate estimated correlation between kappas
  rhoij <- rhoij[[which.corr]]
  corr.kappa <- sign(rhoij)*cfij*abs(rhoij)^d
  
  # Sampling variance matrix of kappa
  gls.scov.matrix <- corr.kappa * sqrt(var.k %*% t(var.k))
  
  return(list(df_summary = tibble(gauge_code = gauges, n = n, kappa = kappas, tau3 = tau3, var_tau3 = var.tau3, var_kappa = var.k),
              scovmatrix = gls.scov.matrix))
  
}