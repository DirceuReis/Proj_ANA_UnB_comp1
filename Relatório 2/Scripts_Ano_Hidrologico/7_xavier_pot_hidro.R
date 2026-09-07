# Xavier — POT + híbrido rbar (passo 7) em R
#
# Mesma metodologia do Uniplu (passo 3):
#   limiar   = hill.adapt()$Xadapt; se NA → quantil 0.99
#   clusters = evd::clusters(..., r = 1)   # 1 dia ≤ u já separa
#   evento   = máximo de cada cluster
#   θ        = DOY * 2π/365 → rbar; mes_oposto = mês(θ+π)
#   híbrido  : rbar ≥ 0,30 → θ+π; senão → último mês seco (15%)
#

# Saída: Dados_Xavier/ano_hidro_out/xavier_ano_hidro_pot.nc

library(ncdf4)
library(terra)
library(extremefit)
library(evd)
library(lubridate)

# --- caminhos ---
DIR_XAVIER <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2024/Dados_Xavier"
DIR_OUT    <- file.path(DIR_XAVIER, "ano_hidro_out")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

F_MENSAL <- file.path(DIR_OUT, "xavier_pr_mensal.nc")
F_OUT    <- file.path(DIR_OUT, "xavier_ano_hidro_pot.nc")
F_CKPT   <- file.path(DIR_OUT, "xavier_ano_hidro_pot_ckpt.rds")
F_REL    <- file.path(DIR_OUT, "xavier_mes_relatorio.rds")

# --- controles ---
SCOPE <- "br"          # "sul" | "br"
CHUNK_LAT <- 3L
RESUME <- TRUE
REBUILD_MENSAL <- FALSE
R_CLUSTER <- 1L
Q_POT <- 0.99
RBAR_HI <- 0.30
SECO_FRAC <- 0.15
MIN_DAYS <- 365L * 10L
MIN_PEAKS <- 5L

BB_SUL <- c(xmin = -58.5, xmax = -47.0, ymin = -34.2, ymax = -21.8)

# ---------------------------------------------------------------------------
list_pr_files <- function() {
  ff <- list.files(DIR_XAVIER, pattern = "^pr_.*_v_3\\.2\\.2\\.nc$", full.names = TRUE)
  ff <- ff[!grepl("_RS_|_bacias|nordeste", basename(ff))]
  if (!length(ff)) stop("Nenhum pr_*.nc em ", DIR_XAVIER)
  sort(ff)
}

maior_vao_circular <- function(quiet) {
  n <- length(quiet)
  if (!any(quiet)) return(NULL)
  if (all(quiet)) return(c(start = 1L, end = n, len = n))
  runs <- list()
  i <- 1L
  while (i <= n) {
    if (quiet[i]) {
      s <- i
      while (i <= n && quiet[i]) i <- i + 1L
      runs[[length(runs) + 1L]] <- c(s, i - 1L)
    } else i <- i + 1L
  }
  best <- NULL
  best_len <- -1L
  add <- function(s, e, len) {
    if (len > best_len) {
      best_len <<- len
      best <<- c(start = s, end = e, len = len)
    }
  }
  wrap <- isTRUE(quiet[1] && quiet[n] && length(runs) >= 2L)
  if (wrap) {
    last <- runs[[length(runs)]]
    first <- runs[[1]]
    add(last[1], first[2], (n - last[1] + 1L) + first[2])
    runs <- if (length(runs) > 2L) runs[2:(length(runs) - 1L)] else list()
  }
  for (rg in runs) add(rg[1], rg[2], rg[2] - rg[1] + 1L)
  best
}

mes_relatorio_vec <- function(p, seco_frac = SECO_FRAC) {
  if (!any(is.finite(p))) return(NA_real_)
  pmin <- min(p, na.rm = TRUE)
  pmax <- max(p, na.rm = TRUE)
  amp <- pmax - pmin
  if (!is.finite(amp) || amp <= 0) return(as.numeric(which.min(p)))
  lim <- pmin + seco_frac * amp
  seco <- is.finite(p) & (p <= lim)
  n_sec <- sum(seco)
  if (n_sec == 0L || n_sec == 12L) return(as.numeric(which.min(p)))
  vao <- maior_vao_circular(seco)
  if (is.null(vao)) return(as.numeric(which.min(p)))
  as.numeric(vao[["end"]])
}

theta_to_month <- function(th) {
  if (!is.finite(th)) return(NA_real_)
  th <- th %% (2 * pi)
  offset <- as.integer(round(th * 365 / (2 * pi))) - 1L
  as.numeric(month(as.Date("1999-01-01") + offset))
}

# POT de uma série diária (= pot_station Uniplu)
pot_celula <- function(x, doy, r = R_CLUSTER) {
  finite <- is.finite(x)
  if (sum(finite) < MIN_DAYS) {
    return(c(n_picos = NA_real_, rbar = NA_real_, mes_oposto = NA_real_))
  }
  x_fit <- x[finite]

  threshold <- tryCatch(hill.adapt(x_fit)$Xadapt, error = function(e) NA_real_)
  if (is.na(threshold) || is.nan(threshold)) {
    threshold <- as.numeric(stats::quantile(x_fit, probs = Q_POT, na.rm = TRUE))
  }
  if (!is.finite(threshold) || threshold <= 0) {
    return(c(n_picos = NA_real_, rbar = NA_real_, mes_oposto = NA_real_))
  }

  peaks <- tryCatch(
    evd::clusters(x, u = threshold, r = r, cmax = FALSE, plot = FALSE),
    error = function(e) NULL
  )
  if (is.null(peaks) || !length(peaks)) {
    return(c(n_picos = 0, rbar = NA_real_, mes_oposto = NA_real_))
  }

  idx <- vapply(peaks, function(cl) {
    nm <- names(cl)
    if (is.null(nm) || !length(nm)) return(NA_integer_)
    as.integer(nm[[which.max(cl)]])
  }, integer(1))
  idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(x)]
  n_pk <- length(idx)

  if (n_pk < MIN_PEAKS) {
    return(c(n_picos = n_pk, rbar = NA_real_, mes_oposto = NA_real_))
  }

  th <- doy[idx] * 2 * pi / 365
  xbar <- mean(cos(th))
  ybar <- mean(sin(th))
  rbar <- sqrt(xbar * xbar + ybar * ybar)
  thetabar <- atan2(ybar, xbar) %% (2 * pi)
  mes_ops <- theta_to_month((thetabar + pi) %% (2 * pi))

  c(n_picos = n_pk, rbar = rbar, mes_oposto = mes_ops)
}

# ---------------------------------------------------------------------------
# Metadados dos NC diários
# ---------------------------------------------------------------------------
files <- list_pr_files()
message("Arquivos diários: ", length(files))
print(basename(files))

nc0 <- nc_open(files[[1]])
lat <- as.numeric(ncvar_get(nc0, "latitude"))
lon <- as.numeric(ncvar_get(nc0, "longitude"))
# pr[longitude, latitude, time]
stopifnot(identical(nc0$var$pr$dim[[1]]$name, "longitude"))
stopifnot(identical(nc0$var$pr$dim[[2]]$name, "latitude"))
stopifnot(identical(nc0$var$pr$dim[[3]]$name, "time"))
nc_close(nc0)
nlat <- length(lat)
nlon <- length(lon)
message("Grade: lon=", nlon, " lat=", nlat,
        " | lat[1]=", lat[1], " lat[n]=", lat[nlat])

# Template SpatRaster (lon × lat) — terra usa extent a partir dos eixos
dx <- mean(diff(lon))
dy <- mean(diff(lat))
r_tmpl <- rast(
  nrows = nlat, ncols = nlon,
  xmin = min(lon) - dx / 2, xmax = max(lon) + dx / 2,
  ymin = min(lat) - abs(dy) / 2, ymax = max(lat) + abs(dy) / 2,
  crs = "EPSG:4326"
)
# terra rows: y decresce do topo; se lat cresce, precisa flip
lat_increasing <- lat[nlat] > lat[1]

# ---------------------------------------------------------------------------
# Cache mensal + último mês seco
# ---------------------------------------------------------------------------
build_mensal <- function(files) {
  if (file.exists(F_MENSAL) && !REBUILD_MENSAL) {
    message("Usando cache mensal: ", F_MENSAL)
    return(rast(F_MENSAL))
  }
  if (REBUILD_MENSAL && file.exists(F_MENSAL)) unlink(F_MENSAL)

  message("Construindo cache mensal (terra::tapp) ...")
  parts <- list()
  for (f in files) {
    message("  ", basename(f))
    t0 <- Sys.time()
    r <- rast(f)
    rm <- tapp(r, "yearmonths", fun = "sum", na.rm = TRUE)
    parts[[length(parts) + 1L]] <- rm
    message("    ", nlyr(rm), " meses em ",
            round(difftime(Sys.time(), t0, units = "secs"), 1), " s")
    rm(r, rm); gc()
  }
  pr_m <- do.call(c, parts)
  if (!is.null(time(pr_m))) pr_m <- pr_m[[order(time(pr_m))]]
  writeCDF(pr_m, F_MENSAL, overwrite = TRUE, varname = "pr",
           longname = "monthly precipitation sum", unit = "mm")
  message("Gravado: ", F_MENSAL)
  pr_m
}

calcular_mes_relatorio <- function(pr_m) {
  if (file.exists(F_REL) && !REBUILD_MENSAL) {
    message("Usando mes_relatorio cache: ", F_REL)
    return(readRDS(F_REL))
  }
  message("Mediana climatológica mensal + último mês seco ...")
  tt <- time(pr_m)
  if (is.null(tt)) stop("Cache mensal sem time().")
  mons <- as.integer(format(tt, "%m"))

  clim <- rast(lapply(1:12, function(m) {
    app(pr_m[[which(mons == m)]], fun = median, na.rm = TRUE)
  }))
  names(clim) <- paste0("m", 1:12)

  V <- values(clim)  # ncell × 12
  out <- rep(NA_real_, nrow(V))
  land <- which(rowSums(is.finite(V) & V > 0.1, na.rm = TRUE) > 0)
  message("  células terra: ", length(land))
  t0 <- Sys.time()
  for (k in seq_along(land)) {
    j <- land[[k]]
    out[j] <- mes_relatorio_vec(as.numeric(V[j, ]))
    if (k %% 20000L == 0L) message("    ", k, "/", length(land))
  }
  message("  feito em ", round(difftime(Sys.time(), t0, units = "secs"), 1), " s")

  r_rel <- rast(clim, nlyrs = 1)
  values(r_rel) <- out
  names(r_rel) <- "mes_relatorio"
  saveRDS(r_rel, F_REL)
  r_rel
}

pr_m <- build_mensal(files)
r_rel <- calcular_mes_relatorio(pr_m)

# Células alvo
df_rel <- as.data.frame(r_rel, xy = TRUE, na.rm = TRUE)
names(df_rel)[3] <- "mes_relatorio"
if (SCOPE == "sul") {
  df_rel <- df_rel[
    df_rel$x >= BB_SUL["xmin"] & df_rel$x <= BB_SUL["xmax"] &
      df_rel$y >= BB_SUL["ymin"] & df_rel$y <= BB_SUL["ymax"],
  ]
} else if (SCOPE != "br") {
  stop("SCOPE deve ser 'sul' ou 'br'")
}
cells_target <- cellFromXY(r_rel, as.matrix(df_rel[, c("x", "y")]))
cells_target <- cells_target[is.finite(cells_target)]
n_target <- length(cells_target)
message("Escopo=", SCOPE, " | células alvo=", n_target)

# lat index (NC) de cada célula
xy_t <- xyFromCell(r_rel, cells_target)
ilat_of <- vapply(xy_t[, 2], function(y) which.min(abs(lat - y)), integer(1))
ilon_of <- vapply(xy_t[, 1], function(x) which.min(abs(lon - x)), integer(1))
names(ilat_of) <- as.character(cells_target)
names(ilon_of) <- as.character(cells_target)

lat_keys <- sort(unique(ilat_of))
lat_groups <- split(cells_target, ilat_of)

# ---------------------------------------------------------------------------
# Leitura fatia latitudinal: pr[lon, lat, time]
# ---------------------------------------------------------------------------
read_time_doy <- function(nc) {
  tunits <- nc$dim$time$units
  tvar <- ncvar_get(nc, "time")
  if (grepl("hours since", tunits, ignore.case = TRUE)) {
    origin <- trimws(sub("(?i)hours since\\s*", "", tunits))
    origin <- sub(" .*$", "", origin)
    times <- as.POSIXct(origin, tz = "UTC") + tvar * 3600
    as.Date(times)
  } else if (grepl("days since", tunits, ignore.case = TRUE)) {
    origin <- trimws(sub("(?i)days since\\s*", "", tunits))
    origin <- sub(" .*$", "", origin)
    as.Date(origin) + as.integer(tvar)
  } else {
    as.Date(tvar, origin = "1970-01-01")
  }
}

# Retorna array [nlon, nlat_slice, ntime_total] e doy
read_lat_block <- function(files, ilat0, ilat1) {
  nlat_s <- ilat1 - ilat0 + 1L
  pr_list <- list()
  doy_list <- list()
  for (f in files) {
    nc <- nc_open(f)
    dates <- read_time_doy(nc)
    ntime <- length(dates)
    # start/count na ordem das dims da variável: lon, lat, time
    arr <- ncvar_get(
      nc, "pr",
      start = c(1L, ilat0, 1L),
      count = c(nlon, nlat_s, ntime)
    )
    nc_close(nc)
    # 1 lat → [lon, time]; várias → [lon, lat, time]
    if (nlat_s == 1L) {
      dim(arr) <- c(nlon, 1L, ntime)
    }
    pr_list[[length(pr_list) + 1L]] <- arr
    doy_list[[length(doy_list) + 1L]] <- yday(dates)
  }
  ntime_all <- sum(vapply(pr_list, function(a) dim(a)[3], integer(1)))
  out <- array(NA_real_, dim = c(nlon, nlat_s, ntime_all))
  i0 <- 1L
  for (a in pr_list) {
    nt <- dim(a)[3]
    out[, , i0:(i0 + nt - 1L)] <- a
    i0 <- i0 + nt
  }
  list(pr = out, doy = unlist(doy_list, use.names = FALSE))
}

# ---------------------------------------------------------------------------
# Saída + checkpoint
# ---------------------------------------------------------------------------
out_n    <- rep(NA_real_, ncell(r_rel))
out_rbar <- rep(NA_real_, ncell(r_rel))
out_ops  <- rep(NA_real_, ncell(r_rel))
out_hib  <- rep(NA_real_, ncell(r_rel))
done     <- rep(FALSE, ncell(r_rel))

if (RESUME && file.exists(F_CKPT)) {
  message("Retomando: ", F_CKPT)
  ck <- readRDS(F_CKPT)
  out_n <- ck$out_n; out_rbar <- ck$out_rbar
  out_ops <- ck$out_ops; out_hib <- ck$out_hib
  done <- ck$done
}

message("POT (hill.adapt + evd::clusters r=", R_CLUSTER, ") ...")
t0_all <- Sys.time()
n_proc <- 0L
n_done0 <- sum(done[cells_target])

for (b in seq(1L, length(lat_keys), by = CHUNK_LAT)) {
  keys <- lat_keys[b:min(b + CHUNK_LAT - 1L, length(lat_keys))]
  cells_b <- unlist(lat_groups[as.character(keys)], use.names = FALSE)
  if (all(done[cells_b])) next

  ilat0 <- min(keys)
  ilat1 <- max(keys)
  message(
    "\nlat NC [", ilat0, ":", ilat1, "]  (",
    which(lat_keys == keys[1]), "/", length(lat_keys), ") ..."
  )
  t0 <- Sys.time()

  sl <- tryCatch(
    read_lat_block(files, ilat0, ilat1),
    error = function(e) {
      message("  ERRO leitura: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(sl)) next

  pr <- sl$pr
  doy <- sl$doy
  mes_rel_vals <- values(r_rel)[, 1]

  for (cell in cells_b) {
    if (done[cell]) next
    ilat <- ilat_of[[as.character(cell)]]
    ilon <- ilon_of[[as.character(cell)]]
    i_loc <- ilat - ilat0 + 1L
    x <- as.numeric(pr[ilon, i_loc, ])

    res <- tryCatch(
      pot_celula(x, doy),
      error = function(e) {
        c(n_picos = NA_real_, rbar = NA_real_, mes_oposto = NA_real_)
      }
    )

    mes_rel <- mes_rel_vals[cell]
    if (is.finite(res[["rbar"]]) && res[["rbar"]] >= RBAR_HI &&
        is.finite(res[["mes_oposto"]])) {
      mes_hib <- res[["mes_oposto"]]
    } else if (is.finite(mes_rel) && mes_rel >= 1 && mes_rel <= 12) {
      mes_hib <- mes_rel
    } else {
      mes_hib <- NA_real_
    }

    out_n[cell] <- res[["n_picos"]]
    out_rbar[cell] <- res[["rbar"]]
    out_ops[cell] <- res[["mes_oposto"]]
    out_hib[cell] <- mes_hib
    done[cell] <- TRUE
    n_proc <- n_proc + 1L

    if (n_proc %% 50L == 0L) {
      total_done <- n_done0 + n_proc
      elapsed <- as.numeric(difftime(Sys.time(), t0_all, units = "secs"))
      rate <- total_done / max(elapsed, 1e-6)
      eta <- (n_target - total_done) / max(rate, 1e-9)
      message(sprintf(
        "  %d/%d (%.1f%%)  %.2f cel/s  ETA %.1f min",
        total_done, n_target, 100 * total_done / n_target, rate, eta / 60
      ))
    }
  }

  saveRDS(
    list(out_n = out_n, out_rbar = out_rbar, out_ops = out_ops,
         out_hib = out_hib, done = done, scope = SCOPE),
    F_CKPT
  )
  message(
    "  faixa em ", round(difftime(Sys.time(), t0, units = "secs"), 1),
    " s | checkpoint ok"
  )
  rm(sl, pr); gc()
}

# ---------------------------------------------------------------------------
# Gravar NC
# ---------------------------------------------------------------------------
message("Gravando ", F_OUT, " ...")
mk <- function(v, nm) {
  r <- rast(r_rel)
  values(r) <- v
  names(r) <- nm
  # zera fora do escopo
  keep <- rep(FALSE, ncell(r))
  keep[cells_target] <- TRUE
  vv <- values(r)
  vv[!keep, 1] <- NA
  values(r) <- vv
  r
}

stk <- c(
  mk(out_n, "n_picos"),
  mk(out_rbar, "rbar"),
  mk(out_ops, "mes_oposto"),
  mk(out_hib, "mes_hibrido_rbar"),
  {
    r <- r_rel
    keep <- rep(FALSE, ncell(r))
    keep[cells_target] <- TRUE
    vv <- values(r)
    vv[!keep, 1] <- NA
    values(r) <- vv
    names(r) <- "mes_relatorio"
    r
  }
)
names(stk) <- c("n_picos", "rbar", "mes_oposto", "mes_hibrido_rbar", "mes_relatorio")

# split=TRUE grava uma variável por camada (permite rast(..., subds="rbar"))
writeCDF(
  stk, F_OUT, overwrite = TRUE, split = TRUE,
  longname = "Xavier POT (hill.adapt + evd) + hibrido rbar"
)

message(
  "Concluído. Escopo=", SCOPE,
  " | feitas=", sum(done[cells_target]), "/", n_target
)
message("Arquivo: ", F_OUT)