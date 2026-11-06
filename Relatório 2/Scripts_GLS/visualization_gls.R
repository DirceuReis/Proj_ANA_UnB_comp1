
# PACKAGES ----------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman, tidyverse, extrafont, patchwork, sf, xtable, ggrepel, ggtext)
# install.package("distributional")


# READ DATA ---------------------------------------------------------------

# Clear environment
rm(list = ls()); invisible(gc())

root <- "Relatório 2/Scripts_GLS/gerados/regression_results/"
transform.path <- paste0(list.dirs(root, recursive = FALSE), "/")
t <- 1 # center
save.path <- transform.path[t]
fit.path <- paste0(save.path, "fit/")
transformation <- str_match(save.path, paste0(root, "(.*?)/"))[2]

message("Transformation applied to results being analysed: ", transformation)

plot.path <- paste0("Relatório 2/Scripts_GLS/figuras/")
if(!dir.exists(plot.path)) dir.create(plot.path)

# Files
files.regional <- list.files(save.path, pattern = "_df_regional_kappa.pqt", full.names = TRUE)
files.regression <- list.files(save.path, pattern = "_regression_sigf.pqt", full.names = TRUE)
files.fit <- list.files(fit.path, pattern = "_fit.rds", full.names = TRUE)
files.kappa.ci <- list.files(fit.path, pattern = "_kappa_ci.rds", full.names = TRUE)
files.quantiles <- list.files(fit.path, pattern = "_quantiles.rds", full.names = TRUE)

# Read data into lists
res.names <- gsub("_df_regional_kappa\\.pqt", "", basename(files.regional))
ls.regional.kappa <- lapply(files.regional, arrow::read_parquet); names(ls.regional.kappa) <- res.names
ls.results.sigf <- lapply(files.regression, arrow::read_parquet); names(ls.results.sigf) <- res.names[!grepl("^average", res.names)]
ls.fit <- lapply(files.fit, readRDS); names(ls.fit) <- res.names
ls.kappa.ci <- lapply(files.kappa.ci, readRDS); names(ls.kappa.ci) <- res.names
ls.quantiles <- lapply(files.quantiles, readRDS); names(ls.quantiles) <- res.names

# Bind
df.regional.kappa <- bind_rows(ls.regional.kappa, .id = "regression") %>% select(-c(mu_kappa_r, sd_kappa_r), everything()) 
df.fit <- bind_rows(ls.fit, .id = "regression")             # GEV parameters
df.kappa.ci <- bind_rows(ls.kappa.ci, .id = "regression")   # confidence intervals on kappa
df.quantiles <- bind_rows(ls.quantiles, .id = "regression") # estimated quantiles
reg.gls.all <- arrow::read_parquet("Relatório 2/Scripts_GLS/gerados/regression_results/center/gls_regression.pqt")

sf.rs <- sf::st_read("Relatório 2/Scripts_GLS/sig/sf_rs.gpkg")
sf.gauges <- sf::st_read("Relatório 2/Scripts_GLS/sig/sf_daily.gpkg")
sf.geomorfo <- sf::st_read("Relatório 2/Scripts_GLS/sig/geomorfologia.gpkg") %>% sf::st_transform(crs = st_crs(sf.rs))

if(!file.exists("Relatório 2/Scripts_GLS/sig/sf_hidro_rs.gpkg")){
  dir.basin <- "Relatório 2/Scripts_GLS/sig/hidrografia/"
  unzip("Relatório 2/Scripts_GLS/Bacias_RS_25.000.kmz", exdir = dir.basin)
  
  sf.hidro <- sf::st_read(paste0(dir.basin, "doc.kml")) %>%
    select(Name) %>% 
    janitor::clean_names() %>%
    mutate(code = substr(name, 1, 4)) %>% 
    arrange(name) %>% 
    slice_tail(n = 1, by = code) %>% 
    mutate(macro = substr(name, 1, 1),
           macro_name = case_when(
             macro == "L" ~ "Coastal",
             macro == "U" ~ "Uruguai",
             macro == "G" ~ "Guaíba"),
           name = sub(".* - Bacia Hidrográfica (do|dos)", "", name))
  
  sf::st_write(sf::st_transform(sf.hidro, crs = st_crs(sf.rs)), "Relatório 2/Scripts_GLS/sig/sf_hidro_rs.gpkg")
  
} else{
  sf.hidro <- sf::st_read("Relatório 2/Scripts_GLS/sig/sf_hidro_rs.gpkg")
}

data <- arrow::read_parquet("Relatório 2/Scripts_GLS/gerados/df_daily_pmax.pqt")  # daily rainfall data
data <- data[!data$gauge_code %in% c("SBPA", "SBUG", "SBSM"),]
df.covariates <- arrow::read_parquet("Relatório 2/Scripts_GLS/gerados/covariates/df_covariates.pqt")

# Source R functions
source("Relatório 2/Scripts_GLS/funcoes/aux_functions.R")              # fun.priori.beta()
source("Relatório 2/Scripts_GLS/funcoes/fun_gev_quantiles_regional.R") # GEV quantiles

which.regression <- "gls"

# SETTING UP DATA FOR VISUALIZATION ---------------------------------------

font <- "sans"
fontsize <- 8
textsize = fontsize/.pt
set_theme(
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = c(0, 0),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key.height = unit(fontsize, "pt"),
      strip.background = element_blank(),
      text = element_text(size = 8),
      panel.spacing.y = unit(0.5, "pt"),
    )
)

# Creating names only for regional average and best GLS models
names.factor <- c("l.gev", "gl.gev", "gl.gev.agls", "gl.gev.gls")
names(names.factor) <- c("MLE", "GMLE", "RGMLE-RMGLS", "RGMLE-GLS")

# Keep only one observation of models l.gev and gl.gev and only one observation
# as there are triplicates (one for each run of the regression functions)]
df.kappa.ci <- df.kappa.ci %>% 
  filter(model != "gl.gev.r") %>% 
  group_by(gauge.code, model) %>% 
  summarise_all(~first(.)) %>% # keep only one observation as they are all the same
  mutate(regression = NA_character_) %>% # remove regression model name as its not needed 
  bind_rows(df.kappa.ci %>% filter(model == "gl.gev.r")) %>% # join regression models
  mutate(model = paste0(gsub("r", "", model), replace(regression, is.na(regression), "")),
         model = gsub("average_", "a", model)) %>% 
  arrange(gauge.code) %>% 
  filter(!grepl("ols|wls", regression)) %>% 
  rowwise() %>% mutate(include_zero = between(0, ci.l, ci.u)) %>% ungroup()

# Compute sd reduction on kappa and quantile estimates
# ML estimates without regional information is used as the base model 
# Positive sd.reduction indicates alternative GMLE model reduced standard deviation
# compared to the standard MLE approach
df.kappa.sd.redct <- df.kappa.ci %>% 
  pivot_wider(id_cols = c(gauge.code),
              names_from = model,
              values_from = sd,
              names_prefix = "sd.") %>% 
  mutate(redct.gl.gev = (sd.l.gev - sd.gl.gev)/sd.l.gev,         # get differences
         redct.gl.gev.gls = (sd.l.gev - sd.gl.gev.gls)/sd.l.gev,
         redct.gl.gev.agls = (sd.l.gev - sd.gl.gev.agls)/sd.l.gev) %>% 
  pivot_longer(cols = c(starts_with("redct.")),
               names_to = "gml.model", # model being compared to standard MLE
               values_to = "sd.reduction") %>% 
  select(gauge.code, gml.model, sd.reduction) %>% 
  mutate(gml.model = str_remove(gml.model, "redct."))

# Investigate how kappa estimates vary within the same gauge across models
df.kappa.rank <- df.kappa.ci %>% 
  pivot_wider(id_cols = c(gauge.code),
              names_from = model,
              values_from = par,
              names_prefix = "kappa.") %>% 
  rowwise() %>% # to compute operations per row
  mutate(sd = sd(c_across(starts_with("kappa."))), # calculate standard deviation across kappa estimates
         range = max(c_across(starts_with("kappa."))) - min(c_across(starts_with("kappa."))),
         n.serie = df.fit$n.serie[match(gauge.code, df.fit$gauge.code)]) %>% 
  ungroup() %>% # get out of rowwise()
  mutate(across(.cols = starts_with("kappa.gl.gev"),        # targets all alternative models
                .fns = list(kappa = ~ .x,                   # kappa value itself
                            diff = ~ .x - kappa.l.gev,      # difference between kappa MLE and kappa GMLE alternatives   
                            decreased = ~ .x < kappa.l.gev, # TRUE if the alternative became more negative or less positive
                            sig.change = ~ sign(.x) != sign(kappa.l.gev)), # TRUE if the sign flipped (e.g., positive to negative)
                .names = "{.col}_{.fn}")) %>%  # creates 'kappa.gl.gev.r.diff', 'kappa.gl.gev.r.decreased', etc
  pivot_longer(cols = matches("kappa\\.gl\\.gev.*_"),
               names_to = c("gml.model", ".value"),
               names_sep = "_") %>% 
  mutate(gml.model = str_remove(gml.model, "kappa.")) %>% 
  select(-starts_with("kappa."))

# Quantiles
df.quantiles <- df.quantiles %>% 
filter(model != "gl.gev.r") %>% 
  group_by(gauge.code, model, tempo.retorno) %>% 
  summarise_all(~first(.)) %>% # keep only one observation as they are all the same
  mutate(regression = NA_character_) %>% # remove regression model name as its not needed 
  bind_rows(df.quantiles %>% filter(model == "gl.gev.r")) %>% # join regression models
  mutate(model = paste0(gsub("r", "", model), replace(regression, is.na(regression), "")),
         model = gsub("average_", "a", model)) %>% 
  filter(!grepl("ols|wls", regression)) %>% 
  arrange(gauge.code, tempo.retorno, model)

df.q.sd.redct <- df.quantiles %>% 
  pivot_wider(id_cols = c(gauge.code, tempo.retorno),
              names_from = model,
              values_from = sd.q,
              names_prefix = "sd.") %>% 
  mutate(redct.gl.gev = (sd.l.gev - sd.gl.gev)/sd.l.gev,         # get sd reduction between RGMLE and MLE
         redct.gl.gev.gls = (sd.l.gev - sd.gl.gev.gls)/sd.l.gev, # get sd reduction between RGMLE and GMLE
         redct.gl.gev.agls = (sd.l.gev - sd.gl.gev.agls)/sd.l.gev) %>% # get sd reduction between RGMLE and GMLE
  pivot_longer(cols = starts_with("redct."),
               names_to = "gml.model",
               values_to = "sd.reduction") %>% 
  mutate(gml.model = str_remove(gml.model, "redct.")) %>% 
  select(gauge.code, tempo.retorno, gml.model, sd.reduction) %>%
  left_join(df.kappa.rank %>% select(gauge.code, diff, decreased, sig.change, gml.model), by = c("gauge.code", "gml.model"))

df.quantiles <- df.quantiles %>% 
  filter(model != "gl.gev.r") %>% 
  group_by(gauge.code, model, tempo.retorno) %>% 
  slice_head(n = 1) %>% 
  mutate(regression = NA_character_) %>% 
  bind_rows(df.quantiles %>% filter(model == "gl.gev.r")) %>% 
  mutate(model = paste0(gsub("r", "", model), replace(regression, is.na(regression), "")),
         model = gsub("average_", "a", model),
         model = factor(model, levels = names.factor, labels = names(names.factor))) %>% 
  arrange(gauge.code, tempo.retorno) %>% select(-regression)

# Select a gauge where shape estimates did not vary much
# gg <- df.kappa.rank$gauge.code[which.min(df.kappa.rank$range)] # smallest kappa range
# gg <- df.kappa.rank$gauge.code[which.min(df.kappa.rank$sd)]    # smallest kappa sd
# gg <- "02955006" # gauge where quantile standard deviation increased the most

# To compare just MLE and AGLS or GLS use this:
# Remember to run everything below when changing which one to compare
gg <- df.kappa.ci %>% # select gauge with smallest difference between AGLS and MLE
  filter(model %in% c("l.gev", "gl.gev.agls")) %>% # compare with averarge GLS
  # filter(model %in% c("l.gev", "gl.gev.gls")) %>%    # compare with local GLS
  group_by(gauge.code, par.name) %>% 
  reframe(par_dif = abs(max(par) - min(par))) %>% 
  slice_min(par_dif) %>%
  pull(gauge.code)

# Check gg kappa estimates
df.kappa.ci %>% filter(gauge.code == gg)


# Frequency plot setup
# Return levels to estimate quantiles
return.level <- c(1.001, 1.005, 1.01, 1.02, 1.03, 1.05, 1.1, 1.2, 1.3, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 10, 20, 50, 100, 200, 500)

# Table of estimated parameters only for selected gauge where kappa didn't vary much
df.fq.par <- df.fit %>% 
  filter(gauge.code == gg) %>% 
  select(gauge.code, regression, rain.mm, n.serie, starts_with(c("par.", "hess.")) & contains("gev"))

# Table of observed frequncy: compute Weibull plotting position
df.fq.observed <- data %>% 
  filter(gauge_code == gg) %>% 
  select(gauge_code, peak_24h) %>% 
  arrange(desc(peak_24h)) %>% 
  mutate(n = length(gauge_code), i = row_number(), p = i/(n + 1), tr = 1/p)

# Estimate quantiles again for new return level vector
df.fq.estimated <- lapply(split(df.fq.par, df.fq.par$regression), function(regression){
  quantile.gev(df = regression,
               tempo.retorno = return.level,
               gauge.names = "gauge.code",
               model = c("l.gev", "gl.gev", "gl.gev.r"))
}) %>% bind_rows(.id = "regression")

df.fq.estimated <- df.fq.estimated %>% 
  filter(model != "gl.gev.r") %>% 
  group_by(gauge.code, model, tempo.retorno) %>% 
  slice_head(n = 1) %>% 
  mutate(regression = NA_character_) %>% 
  bind_rows(df.fq.estimated %>% filter(model == "gl.gev.r")) %>% 
  mutate(model = paste0(gsub("r", "", model), replace(regression, is.na(regression), "")),
         model = gsub("average_", "a", model),
         model = factor(model, levels = names.factor, labels = names(names.factor))) %>% 
  filter(!grepl("ols|wls", regression)) %>% 
  arrange(gauge.code, tempo.retorno) %>% select(-regression)

df.q.sd.redct %>% 
  filter(gml.model == "gl.gev.agls", tempo.retorno == 50) %>% 
  arrange(sd.reduction) %>% 
  head(3) %>%
  pull(gauge.code) -> more.uncertain.q


# VISUALIZATION -----------------------------------------------------------

# suppressMessages(extrafont::loadfonts())
# font <- "Consolas"
# # font <- "Times New Roman"
# # font <- "CM Roman"
digits <- 4
# fontsize <- 8
# textsize <- (fontsize - 2)/.pt
# 
col <- c("#648FFF", "#DC267F", "#FE6100", "#FFB000")
# col <- c("#FFB000", "#FE6100", "#DC267F",  "#648FFF")
colors <- setNames(col, names(names.factor))
# 
# set_theme(theme_minimal() +
#             theme(panel.border = element_rect(color = "black"),
#                   axis.text = element_text(color = "black", family = font, size = fontsize),
#           text = element_text(color = "black", family = font, size = fontsize)))

# caption <- paste("Covariate transformation:", transformation)


## REGRESSION MODEL VARIANCES ---------------------------------------------

df.reg.model <- bind_rows(ls.results.sigf) %>% 
  group_by(model) %>% 
  # slice_min(avp2)  # choose model with best AVP2
  slice_head(n = 1) # choose regional mean model

pos <- position_dodge(0.85)

# Comparison between OLS, WLS and GLS model variances and AVP2
plot.model.variance <- df.reg.model %>% na.omit() %>% 
  mutate(model = toupper(model),
         model = factor(model, c("OLS", "WLS", "GLS"))) %>% 
  select(avp2, var_delta, model) %>% 
  pivot_longer(cols = c(var_delta, avp2),
               names_to = "which_metric",
               values_to = "value") %>%
  ggplot(aes(x = model, y = value, fill = which_metric)) +
  geom_col(position = pos, width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = round(value, digits)), size = textsize, family = font, hjust = 0.5, vjust = -0.6, position = pos) +
  scale_fill_manual(values = c("var_delta" = "grey20", "avp2" = "grey80"),
                    labels = c("var_delta" = expression(sigma[delta]^2), "avp2" = "AVP")) +
  scale_y_continuous(limits = c(0, df.reg.model$avp2[df.reg.model$model == "ols"] + 0.001)) +
  labs(fill = "") +
  theme(legend.position =  c(1,1.08),
        legend.justification = c(1,1),
        legend.key.size = unit(3, "mm"), 
        legend.key.spacing.y = unit(1, "mm"),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_blank(),); plot.model.variance

ggsave(paste0(plot.path, "plot_model_variance.png"), plot.model.variance,
       width = 8, height = 6, units = "cm", dpi = 300)

# Changes in var_delta and avp2
m <- "ols"
(df.reg.model$avp2[df.reg.model$model == m]- df.reg.model$avp2[df.reg.model$model == "wls"])/df.reg.model$avp2[df.reg.model$model == m]*100
(df.reg.model$var_delta[df.reg.model$model == m]- df.reg.model$var_delta[df.reg.model$model == "wls"])/df.reg.model$var_delta[df.reg.model$model == m]*100

(df.reg.model$avp2[df.reg.model$model == m]- df.reg.model$avp2[df.reg.model$model == "gls"])/df.reg.model$avp2[df.reg.model$model == m]*100
(df.reg.model$var_delta[df.reg.model$model == m]- df.reg.model$var_delta[df.reg.model$model == "gls"])/df.reg.model$var_delta[df.reg.model$model == m]*100


## TABLE SUMMARY ----------------------------------------------------------

# Common explanatory variables throughout all regressions
lapply(ls.results.sigf, function(reg){
  
  general.res <- bind_rows(reg)
  valid <- colSums(general.res[,1:(ncol(general.res)-1)]) > 0
  general.res <- general.res[,1:(ncol(general.res)-1)][,valid]
  k <- sum(startsWith(names(general.res), "b"))
  idx <- as.numeric(gsub("b" ,"" , names(general.res)[1:k])) + 1
  valid.exp <- names(df.covariates)[idx]
  
})

# Create a LaTeX table with parameters, standard errors, model variances and avp2
regression.results <- ls.results.sigf[[which.regression]]
which.regression <- regression.results$model[1]
regression.names <- names(regression.results)                         # column names
res <- regression.results[,-length(regression.names)]

# R² calculated previously has been changed in Reis et al. (2020) to the following
res <- res %>% 
  mutate(pseudo_r2 = 1 - var_delta/var_delta[id == 1]) %>% 
  select(where(~any(. != 0)) & !c("sst", "sse", "ssr", "id", "k")) %>% 
  arrange(desc(avp2))

par.table.relatorio2 <- res |> 
  select(starts_with("b"), var_delta, asve, avp2) |> 
  rename(
    Constante = b0, "Z1" = b9, Aspecto = b1, MAM1d = b6, MAM3d = b7,
    "$\\sigma^2_\\delta$" = var_delta, ASVE = asve, AVP = avp2
  ) |> 
  print(n = Inf)

arrow::write_parquet(par.table.relatorio2, paste0(save.path, "latex_table_pars.pqt"))

k <- sum(startsWith(colnames(res), "b")) # number of parameters
n.col <- k*2                             # number of output columns

if(k == 1){
  
  out <- res
  
} else{
  
  out.names <- names(res)
  par.cols <- out.names[grepl("^(b|se)\\d+$", out.names)]
  suffixes <- gsub("^(b|se)", "", par.cols)
  
  ids <- unique(suffixes)
  par.se <- lapply(ids, function(i){
    b.col <- res[, paste0("b", i)]
    se.col <- res[, paste0("se",i)]
    val <- paste0("\\makecell{", t(round(b.col, 4)), " \\\\ (", t(round(se.col, 4)), ")}")
    return(val)
  })
  
  combine <- as.data.frame(do.call(cbind, par.se))
  colnames(combine) <- paste0("b", ids)
  out <- cbind(combine, res[,c("var_delta", "asve", "avp2", "pseudo_r2")])
  empty.str <- "\\makecell{0 \\\\ (0)}"
  out <- out %>% mutate(across(1:all_of(k), ~ str_replace(., fixed(empty.str), "-")))
  
}

# Scale performance metrics by 1000
scale <- 1000

# Pretty names
xtable.results <- out %>% 
  mutate(var_delta = round(var_delta*scale, 2),
         asve = round(asve*scale, 2),
         avp2 = round(avp2*scale,2)) %>% 
  select(b0, b9, b1, b6, b7, var_delta, asve, avp2, pseudo_r2) %>%
  rename("Constant" = "b0", "$Z_1$" = "b9", "$\\text{Aspect}_c$" = "b1", "$\\text{MAM1d}_c$" = "b6", "$\\text{MAM3d}_c$" = "b7",
         "$\\sigma_\\delta^2\\times 10^{-3}$" = "var_delta", 
         "ASV$\\times 10^{-3}$" = "asve",
         "AVP$\\times 10^{-3}$" = "avp2", 
         "$R_\\delta^2$" = "pseudo_r2")

table <- xtable::xtable(x = xtable.results,
                        caption = "GLS regression model results significant at 5\\% level alongide performance metrics.",
                        label = "tab-gls-results",
                        align = rep("c", ncol(xtable.results) + 1))

cat(print(table, type = "latex", include.rownames = FALSE, caption.placement = "top",
          sanitize.colnames.function = identity,
          sanitize.text.function = identity),
    file = paste0(save.path, "latex_table_pars.txt"))

# Check changes in avp2
as.numeric(diff(range(res$avp2))/min(res$avp2) ) %>% scales::percent(accuracy = 0.1)

# 
res

# Quick peak at the range of shape estimates
df.regional.kappa %>% 
  group_by(regression) %>% 
  summarise(kappa_range_l = range(mu_kappa_r)[1], # lower range
            kappa_range_u = range(mu_kappa_r)[2]) # upper range

df.regional.kappa %>% 
  filter(regression == "average_gls")

# Number of gauges with positive and negative shape estimates
df.regional.kappa %>% 
  filter(!grepl("average", regression),
         regression == which.regression) %>% 
  mutate(light_tail = if_else(mu_kappa_r > 0, TRUE, FALSE)) %>% 
  group_by(regression) %>% 
  count(light_tail) %>% 
  arrange(regression)

# Understanding aspect parameter
mean.aspect <- mean(df.covariates$aspect)

## RESIDUALS ACROSS RS ----------------------------------------------------

gls.regional.kappa <- df.regional.kappa[df.regional.kappa$regression == "average_gls",]
# gls.regional.kappa <- df.regional.kappa[df.regional.kappa$regression == "gls",]

residuals <- df.fit %>% 
  rowwise() %>% 
  mutate(kappa.lmom = unlist(par.lmom.gev)[[3]]) %>% 
  ungroup() %>% 
  select(gauge.code, rain.mm, n.serie, kappa.lmom) %>% 
  mutate(kappa.gls = gls.regional.kappa$mu_kappa_r[match(gauge.code, df.regional.kappa$gauge_code)],
         residual = kappa.lmom - kappa.gls)

sf.residual <- sf.gauges %>% 
  filter(!gauge_code %in% c("SBPA", "SBSM", "SBUG")) %>% 
  slice_tail(n = 1, by = gauge_code) %>% 
  mutate(gauge_code = gsub("^0", "", gauge_code),
         residual = residuals$residual[match(gauge_code, residuals$gauge.code)],
         breaks = cut(residual, c(-Inf, -0.1, 0.1, Inf)),
         sign = sign(residual),
         geomorfo = as.numeric(sf::st_within(., sf.geomorfo)))

sf.residual %>% 
  group_by(breaks) %>% 
  count()

sf.residual %>% 
  as.data.frame() %>% 
  group_by(geomorfo,sign) %>% 
  count()

sf.residual %>% 
  as.data.frame() %>% 
  group_by(sign) %>% 
  count()

ggplot() +
  # geom_sf(data = sf.rs, fill = "grey90") +
  geom_sf(data = sf.hidro, aes(fill = factor(macro_name)), alpha = 0.4) +
  geom_sf(data = sf.residual, aes(color = residual), size = 2) +
  # scale_fill_manual(name = "Geo", colors = c(""))
  scale_color_binned(name = "κ - xβ", palette = c("red", "lightpink", "steelblue", "darkblue"),
                     breaks = c(min(sf.residual$residual), -0.1, 0,  0.1, max(sf.residual$residual)),
                     labels = c(round(min(sf.residual$residual),2), -0.1, 0, 0.1, round(max(sf.residual$residual),2))) +
  theme(panel.border = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

## KAPPA ESTIMATES AND UNCERTAINTY ----------------------------------------

# Check if confidence intervals include 0
include.zero <- prop.table(table(df.kappa.ci$include_zero))
message(round(include.zero[[2]]*100, 2), "% out of ", nrow(df.kappa.ci), " total models evaluated predicts κ = 0 in their 95% CI")
df.kappa.ci %>% 
  filter(grepl("gls", model)) %>% 
  group_by(model) %>% 
  count(include_zero) %>% 
  arrange(model)

# Check range in kappa RMGLS estimates
df.kappa.ci %>% filter(model == "gl.gev.agls") %>% arrange(desc(par))

df.kappa.ci %>% 
  filter(model == "gl.gev.agls") %>%
  group_by(model) %>% 
  reframe(min(par), max(par))

# Which gauges do not include zero in their shape CI
df.kappa.ci %>% 
  filter(model == "gl.gev.agls") %>% 
  mutate(light_tail = par > 0) %>% 
  count(light_tail)

# Which are the gauges with largest kappa SD for MLE
df.kappa.ci %>%
  group_by(model) %>% 
  reframe(range = diff(range(sd)),
          min = min(sd),
          max = max(sd),
          iqr = quantile(sd, 0.75) - quantile(sd, 0.25)) %>% 
  mutate(iqr_redc_mle = (iqr[model == "l.gev"] - iqr)/iqr[model == "l.gev"]*100,
         iqr_redc_gmle = (iqr[model == "gl.gev"] - iqr)/iqr[model == "gl.gev"]*100,
         how_smaller_mle = iqr/iqr[model == "l.gev"]*100)

# Range in SD reduction
df.kappa.sd.redct %>% 
  filter(gml.model == "gl.gev.agls") %>% 
  pull(sd.reduction) %>% 
  range()

min.sd.redct.gg <- df.kappa.sd.redct %>% 
  group_by(gml.model) %>% 
  slice_min(sd.reduction) %>% 
  filter(gml.model == "gl.gev.agls") %>% 
  pull(gauge.code)

df.kappa.ci %>% 
  filter(gauge.code == min.sd.redct.gg) %>% 
  mutate(ci.width = ci.u - ci.l,
         width.redct = (ci.width[model == "l.gev"] - ci.width)/ci.width[model == "l.gev"]*100)


gls.model <- c("gl.gev.gls", "gl.gev.agls")

# Boxplots of kappa estimates from different models
plot.kappa.est <- df.kappa.ci %>%  
  # filter(model %in% c("l.gev", "gl.gev", gls.model)) %>%
  filter(model %in% c("l.gev", "gl.gev", gls.model[2])) %>% # comparar somente RMGLS
  mutate(model = factor(model, levels = names.factor, labels = names(names.factor))) %>% # recode changes individual factors
  ggplot(aes(x = model, y = par, fill = model)) +
  # ggplot(aes(x = model, y = sd, color = model)) +
  # geom_jitter(alpha = 0.4, size = 1.2, width = 0.2) +
  geom_jitter(alpha = 0.4, width = 0.2, color = "grey90") +
  # stat_boxplot(geom = "errorbar", width = 0.3, color = "black") +
  # stat_summary(color = "black", fun = "median", geom = "point", size = 1.8) +
  geom_boxplot(linewidth = 0.3, staplewidth = 0.3, staple.linewidth = 0.3, outlier.size = 0.8,
               median.linewidth = 0.3, whisker.linetype = "dashed") +
  scale_color_manual(values = colors, aesthetics = c("fill", "color")) +
  labs(x = "", y = "κ estimates") +
  theme(legend.position =  "none",
        # aspect.ratio = 1,
        plot.background = element_rect(color = "white")); plot.kappa.est

ggsave(paste0(plot.path, "plot_kappa_est.png"), plot.kappa.est,
       width = 12, height = 8, units = "cm", dpi = 300)


# Jitter boxplots of standard deviation on kappa across different models
plot.kappa.sd <- df.kappa.ci %>%  
  # filter(model %in% c("l.gev", "gl.gev", gls.model)) %>%
  filter(model %in% c("l.gev", "gl.gev", gls.model[2])) %>% # comparar somente RMGLS
  mutate(model = factor(model, levels = names.factor, labels = names(names.factor))) %>% # recode changes individual factors
  ggplot(aes(x = model, y = sd, fill = model)) +
  # ggplot(aes(x = model, y = sd, color = model)) +
  # geom_jitter(alpha = 0.4, size = 1.2, width = 0.2) +
  geom_jitter(alpha = 0.4, width = 0.2, color = "grey90") +
  # stat_boxplot(geom = "errorbar", width = 0.3, color = "black") +
  # stat_summary(color = "black", fun = "median", geom = "point", size = 1.8) +
  geom_boxplot(linewidth = 0.3, staplewidth = 0.3, staple.linewidth = 0.3, outlier.size = 0.8,
               median.linewidth = 0.3, whisker.linetype = "dashed") +
  scale_color_manual(values = colors, aesthetics = c("fill", "color")) +
  labs(x = "", y = "SD of κ estimates") +
  theme(legend.position =  "none",
        # aspect.ratio = 1,
        plot.background = element_rect(color = "white")); plot.kappa.sd

ggsave(paste0(plot.path, "plot_kappa_sd_jitter.png"), plot.kappa.sd,
       width = 12, height = 8, units = "cm", dpi = 300)


# Plot empirical cumulative distribution of reduction in kappa standard error across different models
sd.range <- diff(range(df.kappa.sd.redct$sd.reduction))
tol <- 0.03*sd.range
min.sd.reduction <- df.kappa.sd.redct %>% 
  filter(gml.model != "gl.gev.gls") %>%
  # filter(gml.model != "gl.gev.agls") %>% 
  mutate(gml.model = factor(gml.model, levels = names.factor, labels = names(names.factor))) %>% 
  group_by(gml.model) %>%
  summarise(min.reduction = min(sd.reduction)) %>% 
  mutate(label = scales::percent(min.reduction, accuracy = 0.1)) %>% 
  arrange(min.reduction) %>% 
  group_by(close = cumsum(c(TRUE, diff(min.reduction) > tol))) %>% 
  mutate(n = n(), offset = ifelse(n > 1, (row_number() - (n + 1)/2)*tol*sd.range, 0), y = min.reduction + offset) %>% 
  ungroup()

plot.sd.ecdf <- df.kappa.sd.redct %>% 
  filter(gml.model != "gl.gev.gls") %>%
  # filter(gml.model != "gl.gev.agls") %>%
  mutate(gml.model = factor(gml.model, levels = names.factor, labels = names(names.factor))) %>% 
  ggplot(aes(y = sd.reduction, color = gml.model)) +
  geom_hline(aes(yintercept = min.reduction, color = gml.model), min.sd.reduction, linetype = "dashed", alpha = 0.7) +
  geom_text(aes(x = Inf, y = y, label = label), min.sd.reduction, color = "black", hjust = -0.2,
            family = font, size = textsize) +
  stat_ecdf(geom = "point", size = 1.3, show.legend = TRUE, alpha = 0.8) +
  stat_ecdf(geom = "line", linewidth = 0.6, show.legend = TRUE) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = colors) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Decrease in κ SD from MLE", x = "Percentage of gauges", color = "") +
  theme(legend.position = c(0,1.12),
  # theme(legend.position = c(0,1.10),
        legend.justification = c(0,1),
        legend.key.spacing.y = unit(-3, "mm"),
        # aspect.ratio = 1,
        plot.background = element_rect(color = "white"),
        plot.margin = unit(c(1, 15, 1, 1), "mm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
        axis.text = element_text(color = "black")); plot.sd.ecdf

ggsave(paste0(plot.path, "plot_sd_ecdf.png"), plot.sd.ecdf,
       width = 12, height = 8, units = "cm", dpi = 300)


# Show resulting prior distributions
priors <- df.regional.kappa %>%
  filter(!grepl("ols|wls", regression)) %>% 
  filter(gauge_code == gg) %>%
  select(gauge_code, regression, mu_kappa_r, sd_kappa_r) %>% 
  bind_rows(data.frame(gauge_code = gg, regression = NA_character_, mu_kappa_r = -0.1, sd_kappa_r = 0.122)) %>% 
  mutate(regression = if_else(is.na(regression), "gl.gev", paste0("gl.gev.", regression)),
         regression = gsub("average_", "a", regression),
         regression = factor(regression, levels = names.factor, labels = names(names.factor)),
         # label = sprintf("%s (%0.4f, %0.4f\U00B2)", regression, mu_kappa_r, sd_kappa_r)) %>% 
         label = sprintf("%s (%0.4f, %0.4f\U00B2)", regression, -mu_kappa_r, sd_kappa_r)) %>%
  arrange(regression) %>% 
  mutate(label = factor(label, levels = unique(label)))

prior.colors <- setNames(colors[priors$regression], priors$label)

prior.densities <- priors |> 
  crossing(x = seq(-0.5, 0.5, length.out = 1000)) |> 
  rowwise() |> 
  mutate(
    mu_kappa_r = mu_kappa_r*(-1),
    y = fun.priori.beta(x, a = 0.5, mu = mu_kappa_r, sd = sd_kappa_r)
  ) |> 
  ungroup()

plot.prior.densities <- prior.densities %>% 
  # filter(regression != "RGMLE-RMGLS") %>%
  filter(regression != "RGMLE-GLS") %>%
  ggplot(aes(x = x, y = y, color = label, fill = label)) +
  geom_area(alpha = 0.3, position = "identity") +
  # geom_line(linewidth = 0.7) +
  # scale_color_manual(values = prior.colors, aesthetics = c("color", "fill"), name = paste("Gauge:", gg,"\nRegression B(\U03BC, \U03C3²)")) +
  scale_color_manual(values = prior.colors, aesthetics = c("color", "fill"), name = NULL) +
  labs(x = "ξ", y = "π(ξ)") +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.key.size = unit(3, "mm"),
        # legend.key.spacing.x = unit(0, "mm"),
        legend.key.spacing.y = unit(0, "mm"),
        legend.title = element_text(size = fontsize),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(color = "black")); plot.prior.densities

ggsave(paste0(plot.path, "plot_prior_densities.png"), plot.prior.densities,
       width = 15, height = 10, units = "cm", dpi = 600)

 # Combine the three above into a single plot
# plot.kappa <-  
#   plot.kappa.sd + labs(subtitle = "(a)") +
#   plot.sd.ecdf + labs(subtitle = "(b)") +
#   plot_layout(ncol = 1, heights = c(0.45, 0.55)); plot.kappa
# 
# layout <- "
# AAAAABBBBB
# ##CCCCCC##
# "
# 
# plot.kappa <- 
#   (plot.kappa.est + labs(subtitle = "(a)")) + 
#   (plot.kappa.sd + labs(subtitle = "(b)")) +
#   (plot.sd.ecdf  + labs(subtitle = "(c)")) +
#   plot_layout(design = layout); plot.kappa

top.row <- (plot.kappa.est + labs(subtitle = "(a)")) | (plot.kappa.sd + labs(subtitle = "(b)"))
bottom.row <- plot_spacer() | (plot.sd.ecdf + labs(subtitle = "(c)")) | plot_spacer()
bottom.row <- bottom.row + plot_layout(widths = c(1, 4, 1))
plot.kappa <- top.row / bottom.row

# plot.kappa <-  (plot.kappa.sd + plot.sd.ecdf)/plot.prior.densities +
#   plot_layout(heights = c(-1, -1)) &
#   plot_annotation(tag_levels = "a", tag_suffix = ")") &
#   theme(plot.tag = element_text(size = font.size),
#         text = element_text(size = font.size)); plot.kappa

ggsave(paste0(plot.path, "plot_kappa.png"), plot.kappa,
# ggsave(paste0(plot.path, "plot_kappa_all.png"), plot.kappa,
       # width = 10, height = 12, units = "cm", dpi = 300)
       width = 14, height = 12, units = "cm", dpi = 300)

## QUANTILE ANALYSIS AND FREQUENCY PLOT -----------------------------------

# gls.model <- "RGMLE-GLS"
gls.model <- "RGMLE-RMGLS" # check for appropriate 'gg'

freq.colors <- colors[-c(2,4)] # RGMLE-RMGLS
# freq.colors <- colors[-c(2,3)] # RMGLE-GLS
# names(freq.colors) <- paste("Predicted 95% CI", names(freq.colors))
# names(freq.colors) <- paste(names(freq.colors), "predicted 95% CI")
names(freq.colors) <- paste(names(freq.colors), "95% CI")
# freq.colors <- c("Observed" = "midnightblue", freq.colors)
freq.colors <- c("Observed" = "firebrick4", "Predicted" = "midnightblue", freq.colors)
freq.alpha <- setNames(c(1, 1, 0.1, 0.25), names(freq.colors))
freq.linetype <- setNames(c("blank", "blank", "dashed", "solid"), names(freq.colors))

gg.kappa <- df.kappa.ci %>% 
  filter(gauge.code == gg,
         par.name == "kappa",
         model %in% c("l.gev", "gl.gev.agls")) %>%
         # model %in% c("l.gev", "gl.gev.gls")) %>%
  mutate(model = factor(model, levels = names.factor, labels = names(names.factor)))

# Frequency plot
gls.ci.name <- paste(gls.model, "95% CI")
plot.frequency <- df.fq.estimated %>% 
  filter(model %in% c("MLE", gls.model)) %>%
  inner_join(gg.kappa %>% select(gauge.code, model, par), by = c("gauge.code", "model")) %>% 
  mutate(model = paste(model, "95% CI")) %>% 
  ggplot() +
  geom_ribbon(aes(x = tempo.retorno, ymin = ic.l, ymax = ic.u, color = model, fill = model, alpha = model, linetype = model), linewidth = 0.4) +
  geom_line(aes(x = tempo.retorno, y = q, color = "Predicted"), linetype = "solid", linewidth = 1) +
  geom_point(data = df.fq.observed, aes(x = tr, y = peak_24h, color = "Observed"), size = 1.2, shape = 8) +
  scale_x_log10(breaks = c(1, 10, 20, 30, 50, 100, 200, 500)) +
  # scale_color_manual(values = freq.colors, aesthetics = c("color", "fill"), name = paste("Gauge:", gg), limits = names(freq.colors)) +
  geom_text(aes(x = 0.8, y = Inf, label = paste("Gauge:", gauge.code)), stat = "unique", vjust = 2, # GLS
  # geom_text(aes(x = 0.8, y = 265, label = paste("Gauge:", gauge.code)), stat = "unique", # AGLS
            hjust = -0.035, family = font, size = unit(textsize, "mm")) +
  geom_text(aes(x = 0.8, y = Inf, label = paste('Record length:', n.serie, 'years')), stat = "unique", vjust = 3.5, # GLS
  # geom_text(aes(x = 0.8, y = 254, label = paste('Record length:', n.serie, 'years')), stat = "unique", # AGLS
            hjust = -0.025, family = font, size = unit(textsize, "mm")) +
  geom_text(data = ~subset(., model == "MLE 95% CI"),
            aes(x = 0.8, y = Inf, label = paste0("κ[MLE]==", round(par, 4))), parse = TRUE, stat = "unique", vjust = 4.3, # GLS
            # aes(x = 0.8, y = 243, label = paste0("κ[MLE]==", round(par, 4))), parse = TRUE, stat = "unique", # AGLS
            hjust = -0.05, family = font, size = unit(textsize, "mm")) +
  geom_text(data = ~subset(., model == gls.ci.name), # GLS
            aes(x = 0.8, y = Inf, label = paste0("κ[RGMLE]==", round(par, 4))), parse = TRUE, stat = "unique", vjust = 5.4, # GLS
            # aes(x = 0.8, y = 232, label = paste0("κ[RGMLE]==", round(par, 4))), parse = TRUE, stat = "unique",
            hjust = -0.05, family = font, size = unit(textsize, "mm")) +
  scale_color_manual(values = freq.colors, aesthetics = c("color", "fill"), name = "", limits = names(freq.colors)) +
  scale_alpha_manual(values = freq.alpha, name = "", limits = names(freq.colors)) +
  scale_linetype_manual(values = freq.linetype, name = "", limits = names(freq.colors)) +
  labs(x = "Return level [years]", y = "Rainfall depths [mm]") +
  theme(plot.background = element_rect(color = "white"),
        legend.position = "bottom",
        legend.key.size = unit(4.5, "mm"),
        legend.key.spacing.y = unit(1, "mm"),
        legend.spacing.y = unit(-8, "mm")); plot.frequency

ggsave(paste0(plot.path, "plot_frequency.png"), plot.frequency,
       width = 16, height = 12, units = "cm", dpi = 300)

# Changes in quantile standard deviation from MLE to RGMLE
df.q.sd.redct %>% 
  filter(gml.model %in% c("gl.gev.gls", "gl.gev.agls"),
         tempo.retorno == 10) %>% 
  group_by(gml.model) %>% 
  count(decreased) %>% 
  ungroup()

ci.redct <- df.quantiles %>% 
  filter(gauge.code == gg,
         # tempo.retorno == 100,
         model %in% c("MLE", "RGMLE-RMGLS")) %>% 
  mutate(width = ic.u - ic.l) %>% print()

df.quantiles %>% 
  filter(gauge.code == gg,
         model %in% c("MLE", "RGMLE-RMGLS")) %>% 
  mutate(width = ic.u - ic.l) %>% 
  group_by(tempo.retorno) %>% 
  reframe(ci_redct = (width[model == "MLE"] - width[model == "RGMLE-RMGLS"])/width[model == "MLE"]*100)

rl <- 50

# Gauges which quandile SD increased > 100%
gauges.out <- df.q.sd.redct %>% 
  mutate(gml.model = factor(gml.model, levels = names.factor, labels = names(names.factor))) %>% 
  filter(tempo.retorno == rl,
         gml.model == gls.model) %>% 
  # arrange(desc(sd.reduction)) %>%
  arrange(sd.reduction) %>%
  # filter(abs(sd.reduction) >= 1) %>% 
  head(3) %>%
  mutate(label = paste0(gauge.code, " (", round(sd.reduction, 3)*100, "%)")) %>% 
  arrange(desc(sd.reduction))

count.pos <- df.q.sd.redct %>% 
  mutate(gml.model = factor(gml.model, levels = names.factor, labels = names(names.factor))) %>% 
  filter(tempo.retorno == rl,
         gml.model == gls.model) %>% 
  mutate(quadrant = case_when(
    # RGMLE kappa < MLE kappa and RGMLE sd increased
    # meaning that RGMLE fitted a lower kappa and a heavier tail might have increased quantile sd
    diff < 0 & sd.reduction < 0 ~ "bottom-left",
    # RGMLE kappa > MLE kappa and RGMLE sd increased
    # meaning that RGMLE fitted a higher kappa and lighter tail, but despite that quantile sd increased
    diff > 0 & sd.reduction < 0 ~ "bottom-right",
    # RGMLE kappa < MLE kappa and RGMLE sd decreased
    # meaning that despite a heavier tail, quantile sd decreased
    diff < 0 & sd.reduction > 0 ~ "top-left",
    # RGMLE kappa > MLE kappa and RGMLE sd decreased
    # meaning that a lighter tail and quantile sd decreased
    diff > 0 & sd.reduction > 0 ~ "top-right"
  )) %>% 
  group_by(quadrant) %>%
  count() %>% 
  mutate(x = if_else(grepl("left", quadrant), -Inf, 0),
         y = if_else(grepl("top", quadrant), Inf, 0)) %>% print()

plot.q.sd <- df.q.sd.redct %>% 
  mutate(gml.model = factor(gml.model, levels = names.factor, labels = names(names.factor))) %>% 
  filter(tempo.retorno == rl,
         gml.model == gls.model,
         !gauge.code %in% gauges.out$gauge.code) %>% 
  ggplot(aes(x = diff, y = sd.reduction)) +
  geom_point(alpha = 0.4, color = "dodgerblue3", size = 1, shape = 16) +
  geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed", alpha = 0.7, color = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed", alpha = 0.7, color = "grey30") +
  annotate("curve", xend = -0.3, x = -0.2, yend = -1, y = -0.81, linewidth = 0.35, curvature = 0.1, arrow = arrow(type = "open", length = unit(1.5, "mm"))) +
  annotate("text", x = -0.22, y = -0.80, label = paste(gauges.out$label, collapse = "\n"), hjust = -0.2, family = font, size = textsize) +
  geom_text(data = count.pos, aes(x = x, y = y, label = paste(n, "gauges")), hjust = -0.1, vjust = 2, family = font, size = textsize) +
  scale_y_continuous(labels = scales::percent, breaks = seq(-1, 1, 0.25)) +
  labs(x = "(κ<sub>RGMLE</sub> - κ<sub>MLE</sub>)", y = "Quantile SD reduction") +
  # labs(x = "(κ<sub>RGMLE</sub> - κ<sub>MLE</sub>)", y = "(SD<sub>MLE</sub> - SD<sub>RGMLE</sub>)/SD<sub>MLE</sub>") +
  theme(legend.position = "bottom",
        plot.background = element_rect(color = "white"),                    
        aspect.ratio = 1, # ensure square plot shape when plotting
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown()); plot.q.sd

ggsave(paste0(plot.path, "plot_q_sd.png"), plot.q.sd,
       width = 13, height = 12, units = "cm", dpi = 300)

# Combine plots from frequency analysis and quantilsd.range
plot.quantile <-
  plot.q.sd + labs(subtitle = "(a)") +
  plot.frequency + labs(subtitle = "(b)") +
  plot_layout(ncol = 2, guides = "collect", widths = c(0.8, 1)) &
  theme(legend.position = "bottom",
        legend.key.size = unit(fontsize, "pt"),
        legend.margin = margin(t = -5, unit = "pt")); plot.quantile

ggsave(paste0(plot.path, "plot_quantile.png"), plot.quantile,
  width = 14, height = 7, units = "cm", dpi = 300)

# Quanto o valor de kappa diminuiu nas estações que ficaram de fora do scatter plot
df.kappa.ci %>% 
  filter(model %in% c("l.gev", "gl.gev.agls"),
         gauge.code %in% more.uncertain.q) %>% 
  group_by(gauge.code) %>% 
  reframe(drop = (par[model == "l.gev"] - par[model == "gl.gev.agls"])/par[model == "l.gev"]*100)

# Confidence intervals
n.subset <- 50
# subset.gauges <- sample(unique(df.fit$gauge.code), n.subset)
# subset.gauges <- df.fit %>% 
#   group_by(gauge.code) %>% 
#   slice_head(n = 1) %>% 
#   arrange(desc(n.serie)) %>% 
#   head(n.subset) %>% 
#   pull(gauge.code)

# Select gauges with largest differences in kappa estimates
subset.gauges <- df.kappa.rank %>% 
  filter(gml.model == "gl.gev.gls") %>% 
  arrange(desc(abs(diff)), n.serie) %>% 
  head(n.subset)

# order <- df.kappa.ci %>% 
#   filter(model == gls.model, par.name == "kappa") %>% 
#   arrange(sd) %>% 
#   pull(gauge.code) %>% 
#   unique()

order <- df.quantiles %>% 
  filter(model == gls.model, tempo.retorno == rl) %>% 
  mutate(width = ic.u - ic.l) %>%
  # arrange(width) %>%
  arrange(q) %>%
  pull(gauge.code) %>% 
  unique()
  
dodge.width <- position_dodge(0.8)

## Visualize differences in MLE and RGMLE with kappa change
plot.q.conf <- df.quantiles %>% 
  filter(gauge.code %in% subset.gauges$gauge.code,
         model %in% c("MLE", gls.model),
         tempo.retorno %in% c(10, 50, 100)) %>%
  mutate(gauge.code = factor(gauge.code, levels = order)) %>%
  # inner_join(subset.gauges %>% select(gauge.code, gml.model, diff, decreased), by = c("gauge.code", "gml.model"))
  filter(tempo.retorno == 50) %>%
  ggplot(aes(y = gauge.code, color = model)) +
  # facet_wrap(~tempo.retorno, ncol =c 1, scale = "free_y", label = as_labeller(function(x) paste0(x, "-yr"))) +
  geom_errorbar(aes(xmin = ic.l, xmax = ic.u), position = dodge.width, width = 0, linewidth = 0.6, alpha = 0.6) +
  geom_point(aes(x = q), position = dodge.width, size = 0.8) +
  geom_text(data = subset.gauges, aes(x = -Inf, y = gauge.code, label = ifelse(decreased == TRUE, "\U25BC", "\U25B2")), color = "red", hjust = -0.5, size = 2.5) +
  scale_color_manual(values = colors) +
  labs(y = "", x = "Quantile estimates [mm]", color = "") +
  theme_minimal() +
  theme(plot.background = element_rect(color = "white"),
        legend.position = "bottom",
        legend.margin = margin(t = -10),
        panel.border = element_rect(color = "black", fill = NA),
        # axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(family = font, color = "black", size = font.size)); plot.q.conf

ggsave(paste0(plot.path, "plot_q_conf_arrows.png"), plot.q.conf,
       width = 16, height = 18, units = "cm", dpi = 300)

# Visualize differences in MLE, GMLE and RGMLE
plot.q.conf <- df.quantiles %>% 
  filter(gauge.code %in% gauges.long,
         model %in% c("MLE", gls.model),
         tempo.retorno %in% c(10, 50, 100)) %>%
  mutate(gauge.code = factor(gauge.code, levels = order)) %>%
  filter(tempo.retorno == rl) %>%
  ggplot(aes(x = gauge.code, color = model)) +
  facet_wrap(~tempo.retorno, ncol = 1, scale = "free_y", label = as_labeller(function(x) paste0(x, "-yr"))) +
  geom_errorbar(aes(ymin = ic.l, ymax = ic.u), position = dodge.width, width = 0, linewidth = 0.6, alpha = 0.6) +
  geom_point(aes(y = q), position = dodge.width, size = 0.8) +
  scale_color_manual(values = colors) +
  labs(x = "", y = "Quantile estimates [mm]", color = "") +
  theme_minimal() +
  theme(plot.background = element_rect(color = "white"),
        legend.position = "bottom",
        legend.margin = margin(t = -10),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(family = font, color = "black", size = font.size)); plot.q.conf

ggsave(paste0(plot.path, "plot_q_conf.png"), plot.q.conf,
       width = 16, height = 18, units = "cm", dpi = 300)
