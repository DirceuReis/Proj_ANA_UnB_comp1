# Script dedicated to building a model for sample correlation
# and distance between gauges


# PACKAGES ----------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman, tidyverse, arrow, sf, data.table, ggtext)


# FUNCTIONS ---------------------------------------------------------------

# Fit models with optim
model.fit <- function(m, data){
  
  # Preparationn
  model <- m[["model"]]         # expression
  start <- unlist(m[["start"]]) # starting parameters
  model.parts <- trimws(unlist(strsplit(model, "~"))) # split expression parts
  target <- model.parts[1]      # zij
  expression <- model.parts[2]  # expression
  
  # Objective function (RMSE)
  rmse <- function(par){
    
    env <- c(as.list(par), as.list(data))
    sim <- eval(parse(text = expression), envir = env)
    obs <- data[[target]]
    
    rmse <- sqrt(mean((obs - sim)^2, na.rm = TRUE))
    
  }
  
  # Optim
  optim.res <- optim(par = start, fn = rmse)
  par <- optim.res$par
  
  # Final expression
  pattern <- c(paste0("\\b(", paste(names(start), collapse = "|"), ")\\b"))
  template <- gsub(pattern, "{\\1}", model)
  
  # Predicted  values
  data.sim <- data.frame(dij = data[["dij"]])
  env <-  c(as.list(par), as.list(data.sim))
  data.sim$zij <- eval(parse(text = expression), envir = env)
  
  out <- list(formula = glue::glue_data(round(par, 3), template),
              model = model,
              par = par,
              rmse = optim.res$value,
              sim = data.sim,
              convergence = optim.res$convergence)
  
  return(out)
  
}

# Compute predicted zij from dij
zij.predict <- function(formula, par, dij){
  
  expression <- trimws(unlist(strsplit(formula, "~")))[2]
  env <- c(as.list(par), as.list(dij))
  zij <- eval(parse(text = expression), envir = env)
  
  return(zij)
  
}


# READ DATA ---------------------------------------------------------------

# Read daily AMS
data <- arrow::read_parquet("base/gerados/df_daily_pmax.pqt")
sf.gauges <- slice_tail(st_read("base/sig/sf_daily.gpkg", quiet = TRUE), n = 1, by = gauge_code)

# Temporary:
remove.gauges <- c("SBPA", "SBSM", "SBUG") # ICEA are not daily
data <- data %>% 
  filter(!gauge_code %in% remove.gauges) %>% 
  arrange(gauge_code)

sf.gauges <- sf.gauges[!sf.gauges$gauge_code %in% c("SBPA", "SBSM", "SBUG"),]
sf.gauges <- sf.gauges %>% 
  filter(!gauge_code %in% remove.gauges) %>% 
  mutate(gauge_code = sub("^0", "", gauge_code)) %>% 
  arrange(gauge_code)

st_crs(sf.gauges) # Sirgas 2000, reprojetar
sf.gauges <- sf::st_transform(sf.gauges, crs = 31982)

gauge.order <- sort(unique(sf.gauges$gauge_code))


# CORRELATION -------------------------------------------------------------

# Correlation between concurrent gauge AMS
# Store AMS in columns
dt <- as.data.table(data)            # to use dcast function from data.table package
cc.ams <- data.table::dcast(dt, year~gauge_code, value.var = "peak_24h") # cuncurrent AMS
is.data <- !is.na(cc.ams[, -1])      # remove year column, create a "mask"
nij <- crossprod(as.matrix(is.data)) # count how many common years: t(data) %*% data
rij <- cor(cc.ams[,-1], use = "pairwise.complete.obs") # correlation between cc.ams
zij <- log((1 + rij)/(1 - rij))      # fisher transform (logit) (ou log(base = 10)?) 
dij <- matrix(sf::st_distance(sf.gauges), nrow = nrow(sf.gauges))/1000 # distance between gauges [km]
n <- colSums(is.data)                # count TRUEs in each colum
cfij <- nij/sqrt(n %*% t(n))         # record length correction factor
names <- dimnames(nij)

# Build correlation moodel
# Once nij, rij and dij are pairwise operations, every data appears twice
u.matrix <- upper.tri(nij, diag = FALSE) # create upper triangular matrix mask

# Prepare data
data.obs <- data.frame(rij = rij[u.matrix],
                       zij = zij[u.matrix],
                       dij = dij[u.matrix],
                       nij = nij[u.matrix])

# Proposed models to smooth out correlation
models <- list(m1 = list(model = "zij ~ a + exp(b - c*dij)", start = list(a = min(zij), b = log(max(zij[u.matrix]) - min(zij[u.matrix])), c = 0.01)),
               m2 = list(model = "zij ~ exp(a - b*((dij^c - 1)/c))", start = list(a = 1, b = -0.1, c = 0.1)))

# Fit models
corr.res <- lapply(models, model.fit, data.obs)

# Final cross-correlation for each model fitte
rhoij <- lapply(corr.res, function(m){
  
  model <- m[["model"]]
  par <- m[["par"]]
  
  zij <- zij.predict(model, par, dij)        # compute zij matrix
  rhoij <- (exp(2*zij) - 1)/(exp(2*zij) + 1) # matrix
  diag(rhoij) <- 1                           # guarantee same-site corr to be 1
  dimnames(rhoij) <- names                   # put gauge names back
  
  stopifnot("Mismatch in covariate rownames" = identical(rownames(rhoij), gauge.order))
  stopifnot("Mismatch in covariate colnames" = identical(rownames(rhoij), gauge.order))
  
  return(rhoij)
  
})

# rm(list = setdiff(ls(), c("rhoij", "cfij"))); invisible(gc())

# VISUALIZATION -----------------------------------------------------------

# # Correlation Z-transformed models
# colors <- setNames(c("blue", "red3"), c(corr.res[["m1"]][["formula"]], corr.res[["m2"]][["formula"]]))
# 
# plot.corr <- ggplot(mapping = aes(x = dij, y = zij)) +
#   geom_point(data = data.obs, alpha = 0.5, color = "grey70") +
#   geom_line(data = corr.res[["m1"]][["sim"]], aes(color = corr.res[["m1"]][["formula"]]), linewidth = 0.9) +
#   geom_line(data = corr.res[["m2"]][["sim"]], aes(color = corr.res[["m2"]][["formula"]]), linewidth = 0.9) +
#   scale_color_manual(values = colors) +
#   labs(x = "Distance [km]", y = "Fisher Z-transformed cross-correlations", color = "") +
#   theme_minimal() +
#   theme(legend.position = c(1, 0.97),
#         legend.justification = c(1, 0.5),
#         legend.text = element_text(family = "mono", face = "bold"),
#         panel.border = element_rect(color = "black"),
#         text = element_text(family = "serif"))
# 
# ggsave(plot = plot.corr, filename = "artigo/figuras/plot.corr.png",
#        height = 12, width = 16, units = "cm", dpi = 300)



