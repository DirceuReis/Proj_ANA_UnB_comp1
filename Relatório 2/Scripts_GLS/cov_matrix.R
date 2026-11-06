
# OBJECTS -----------------------------------------------------------------

# List of objects and what they are
# data: data.frame w/ AM series for different gauges in long format
# m: number of gauges (AM series)
# kappas: vector w/ lmom kappa estimates, length 'm'
# tau3: vector w/ lmom tau3 (l-skew) estimates, length 'm'


# PACKAGES ----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, readxl, janitor, lmom, gsl, arrow)


# CORRELATION -------------------------------------------------------------

# Source correlation RScript
source("scripts/correlation_model.r")


# READ DATA ---------------------------------------------------------------

# Daily files RS
# path_rs_daily <- file.choose()
# df_rs_daily <- read_excel(path_rs_daily)
path <- "base/df_daily_pmax.pqt"
data <- arrow::read_parquet(path) # daily rainfall annual maximum series
col.names <- names(data)[c(1,4,6,7,12,13,17)]
data <- as_tibble(clean_names(data[, col.names]))

# Martins and Stedinger's (2001) power exponents for cross-correlation
d.exponents <- data.frame(kappa = c(-0.3, -0.2, -0.1, 0, 0.1), delta = c(2.4, 2.6, 2.8, 2.9, 3.1))


# GET L-MOMENTS AND KAPPA -------------------------------------------------

gauges <- unique(data$gauge_code)
m <- length(gauges)
kappas <- rep(0, m)
tau3 <- rep(0, m)
n <- rep(0, m)

# Get L-moments for all stations and fit kappa
for(gg in seq_along(gauges)){
  
  # Get am series for current gauge
  cod <- gauges[[gg]]
  curr.gg <- data[data[col.names[1]] == cod, col.names[5], drop = TRUE] # vector
  n[gg] <- length(curr.gg)
  
  # Calculate sample lmoments and ratios
  lmom <- lmom::samlmu(x = curr.gg, sort.data = TRUE) # sample lmoments
  tau3[gg] <- lmom[[3]]                               # vector, l-skewness (3rd lmoment ratio)
  kappas[gg] <- lmom::pelgev(lmom)[[3]]               # vector, fit GEV and extract kappa
  
}

# Clear memory
rm(gg, lmom, cod, curr.gg); invisible(gc())

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

# Compute var(tau3) and var(kappa); see Chowdhury et al. (1991) and Fill [PhD Thesis Cornell] (1994)
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


# BUILD SAMPLE COVARIANCE MATRIX ------------------------------------------

# WLS: scov matrix consists only of its diagonal filled with estimated kappa variances
# Create empty matrix and fill diagonals with kappa variance
wls.scov.matrix <- matrix(0, ncol = m, nrow = m, dimnames = list(gauges, gauges))
diag(wls.scov.matrix) <- var.k

# GLS: scov matrix consists of
# diagonal: squared sample kappa variance
# others: estimated kappa correlation and sqrt kappa variances in i and j
rhoij <- rhoij[["m2"]]

# Interpolate mean kappa
pred.d <- approxfun(x = d.exponents$kappa, y = d.exponents$d)
d <- pred.d(mean(kappas))

# Calculate estimated correlation between kappas
corr.kappa <- sign(rhoij)*cfij*abs(rhoij)^d

# Sampling variance matrix of kappa
gls.scov.matrix <- corr.kappa * sqrt(var.k %*% t(var.k))


# TEST FUNCTION -----------------------------------------------------------

# Get fun_wls_scovmatrix()
source("scripts/funcoes/fun_scov_matrix.r")

# Test function
results <- fun_wls_scovmatrix(df = df_rs_daily, col = c("gauge_code", "peak_24h"))
wls_scovmatrix <- results[["scovmatrix"]]