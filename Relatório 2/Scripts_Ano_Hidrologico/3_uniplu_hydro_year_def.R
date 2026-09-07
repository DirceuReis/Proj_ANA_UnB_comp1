# Ano hidrológico (POT + média circular) para as séries Uniplu
# Entrada: stations/uniplu_analyzed e dataframes/out_analise_uniplu.rds
# Saída:   dataframes/df_hydroyear_uniplu.rds

library(dplyr)
library(tidyr)
library(lubridate)
library(extremefit)
library(evd)

# --- caminhos (saídas nesta pasta; Uniplu fora) ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_ANALYZED <- file.path(DIR_OUT, "stations_analyzed")
MIN_YEARS_GOOD <- 30

setwd(DIR_FLUXO)

pot_station <- function(sta, serie, r = 1) {

  x <- sta[[serie]]
  x_fit <- x[is.finite(x)]

  threshold <- hill.adapt(x_fit)$Xadapt

  if (is.na(threshold) | is.nan(threshold)) {
    threshold <- quantile(x_fit, probs = .99, na.rm = TRUE)
  }

  peaks <- evd::clusters(sta[[serie]], u = threshold, r = r,
                         cmax = FALSE, plot = FALSE)

  n <- length(peaks)

  index_peaks <- vector(mode = "character", length = n)

  for (i in 1:length(peaks)) {
    index_peaks[i] <- names(peaks[[i]])[which.max(peaks[[i]])]
  }

  index_peaks <- as.integer(index_peaks)

  peaks <- sta[index_peaks, ]

  peaks <- peaks %>%
    mutate(
      jday  = yday(dt),
      theta = jday * 2 * pi / 365
    )

  return(tibble(
    codigo   = unique(sta$codigo),
    xbar     = mean(cos(peaks$theta)),
    ybar     = mean(sin(peaks$theta)),
    thetabar = atan(ybar / xbar),
    thetavar = var(peaks$theta),
    rbar     = sqrt(xbar^2 + ybar^2)
  ))
}

sta.inf <- readRDS(file.path(DIR_DF, "out_analise_uniplu.rds"))

sta.inf <- sta.inf %>%
  filter(year.size > MIN_YEARS_GOOD)

if (!is.null(TEST_UF)) {
  sta.inf <- sta.inf %>% filter(estado == TEST_UF)
  message("Filtro TEST_UF=", TEST_UF)
}

cods <- as.character(sta.inf$codigo)
n <- length(cods)
message("Postos no POT: ", n)

out <- tibble(
  codigo   = character(),
  xbar     = double(),
  ybar     = double(),
  thetabar = double(),
  thetavar = double(),
  rbar     = double()
)

for (i in 1:n) {

  f_sta <- file.path(DIR_ANALYZED, paste0(cods[i], "_analyzed.rds"))
  if (!file.exists(f_sta)) {
    next
  }

  sta <- readRDS(f_sta)

  tmp <- tryCatch(
    pot_station(sta, serie = "p", r = 1),
    error = function(e) NULL
  )

  if (!is.null(tmp)) {
    out <- rbind(out, tmp)
  }

  if (i %% 50 == 0) {
    print(paste0(cods[i], " total: ", i, "/", n, " válidas: ", nrow(out)))
  }
}

out <- out %>%
  mutate(quadrante = case_when(
    xbar >= 0 & ybar >= 0 ~ 1,
    xbar <  0 & ybar >= 0 ~ 2,
    xbar <  0 & ybar <  0 ~ 3,
    xbar >= 0 & ybar <  0 ~ 4
  ))

out <- out %>%
  mutate(theta = case_when(
    xbar >= 0 & ybar >= 0 ~   thetabar,
    xbar <  0 & ybar >= 0 ~   pi - abs(thetabar),
    xbar <  0 & ybar <  0 ~   pi + abs(thetabar),
    xbar >= 0 & ybar <  0 ~ 2 * pi - abs(thetabar)
  ))

out <- out %>%
  left_join(
    sta.inf %>% select(codigo, estado, nome, lat, long, network, year.size,
                       year.first, year.last),
    by = "codigo"
  )

saveRDS(out, file.path(DIR_DF, "df_hydroyear_uniplu.rds"))
message("Salvo df_hydroyear_uniplu.rds: ", nrow(out), " postos")

