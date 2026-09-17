# =============================================================================
# FUNÇÃO — IMAX POR ANO HIDROLÓGICO
# =============================================================================
#
# Baseada na lógica histórica da fun_imax_agg(), com adaptações para:
#
# 1. ano hidrológico definido por start_month;
# 2. QC anual aplicado ao ano hidrológico completo;
# 3. períodos estruturalmente ausentes nas bordas contam como falha;
# 4. ano aceito se na_prct_yr <= na_accept;
# 5. roll_sum mantém na.rm = TRUE;
# 6. janelas 100% NA são invalidadas;
# 7. janelas temporalmente descontínuas são invalidadas;
# 8. janelas que atravessam a fronteira do ano hidrológico são invalidadas.
#
# =============================================================================

imax.wateryear <- function(data, durations, start_month = 1, which.mon = 1:12,
                           names = c("datetime","rain_mm"),
                           na_accept = 0.20, apply_qc = TRUE){
  
  if(!require(pacman)){
    message("Instalando pacote 'pacman'...")
    install.packages("pacman")
  }
  pacman::p_load(pacman, RcppRoll, pbapply, lubridate, data.table)
  
  if(!inherits(data,"list")) stop("'data' deve ser uma lista.")
  if(!is.numeric(durations) || any(!is.finite(durations)) || any(durations <= 0))
    stop("'durations' deve conter valores positivos em horas.")
  if(!is.numeric(na_accept) || length(na_accept) != 1 || !is.finite(na_accept) ||
     na_accept < 0 || na_accept > 1) stop("'na_accept' deve estar entre 0 e 1.")
  
  which.mon <- unique(as.integer(which.mon))
  if(any(which.mon < 1 | which.mon > 12)) stop("'which.mon' deve conter meses entre 1 e 12.")
  
  gauges <- names(data)
  if(is.null(gauges)) stop("'data' deve ser uma lista nomeada.")
  
  get_start_month <- function(gauge){
    if(length(start_month) == 1) return(as.integer(start_month))
    if(!is.null(names(start_month)) && gauge %in% names(start_month))
      return(as.integer(start_month[[gauge]]))
    stop("'start_month' deve ter comprimento 1 ou ser vetor nomeado pelos gauges.")
  }
  
  message("Agregando e extraindo intensidades máximas anuais para durações: ",
          paste(durations, collapse=", "))
  message("Meses analisados: ", paste(sort(which.mon), collapse=", "))
  message("QC anual: ",
          ifelse(apply_qc,
                 paste0("remover anos com mais de ", na_accept*100, "% de falhas"),
                 "não aplicado"))
  
  mon.filter <- paste(which.mon, collapse=",")
  track.env <- new.env()
  track.env$skipped_count <- 0L
  track.env$skipped_details <- list()
  track.env$qc_removed_count <- 0L
  track.env$qc_removed_details <- list()
  
  imax.gauge <- pbapply::pblapply(gauges, function(gauge){
    
    data.gauge <- data[[gauge]]
    if(!is.data.frame(data.gauge))
      stop("Elemento ", gauge, " deveria ser data.frame.")
    
    if(sum(is.element(names[1:2], names(data.gauge))) != 2)
      stop("Estação ", gauge, " não contém ", names[1], " e ", names[2], ".")
    
    start.month <- get_start_month(gauge)
    if(!is.finite(start.month) || start.month < 1 || start.month > 12)
      stop("start_month inválido para ", gauge, ".")
    
    data.gauge <- data.gauge[
      order(data.gauge[[names[1]]]),
      , drop=FALSE
    ]
    data.gauge <- data.gauge[
      !is.na(data.gauge[[names[1]]]),
      , drop=FALSE
    ]
    
    if(nrow(data.gauge) < 2) return(NULL)
    
    dates <- data.gauge[[names[1]]]
    depths <- as.numeric(data.gauge[[names[2]]])
    na.depths <- is.na(depths)
    mon.all <- lubridate::month(dates)
    
    d.base.hr <- as.numeric(dates[2] - dates[1], units="hours")
    if(!is.finite(d.base.hr) || d.base.hr <= 0)
      stop("Resolução temporal inválida em ", gauge, ".")
    
    ratio <- durations/d.base.hr
    if(any(abs(ratio-round(ratio)) > 1e-8))
      stop("Há duração não múltipla da resolução base de ",
           d.base.hr*60, " min em ", gauge, ".")
    
    shift <- ifelse(start.month > 1 & mon.all >= start.month, 1L, 0L)
    wateryear <- lubridate::year(dates) + shift
    mon.idx <- mon.all %in% which.mon
    years <- sort(unique(wateryear[mon.idx]))
    if(!length(years)) return(NULL)
    
    tz <- tryCatch(lubridate::tz(dates[1]), error=function(e) "UTC")
    if(is.na(tz) || !nzchar(tz)) tz <- "UTC"
    
    year.qc <- lapply(years, function(yr){
      
      if(start.month == 1){
        wy.start <- lubridate::make_datetime(yr, 1, 1, tz=tz)
        wy.end <- lubridate::make_datetime(yr+1, 1, 1, tz=tz)
      } else {
        wy.start <- lubridate::make_datetime(yr-1, start.month, 1, tz=tz)
        wy.end <- lubridate::make_datetime(yr, start.month, 1, tz=tz)
      }
      
      n.expected <- round(
        as.numeric(difftime(wy.end, wy.start, units="secs")) /
          (d.base.hr*3600)
      )
      
      idx <- wateryear == yr
      n.present <- sum(idx)
      n.na.present <- sum(na.depths[idx])
      n.missing.structural <- max(0L, as.integer(n.expected - n.present))
      n.na.total <- min(n.expected, n.na.present + n.missing.structural)
      
      na.prct.yr <- if(n.expected > 0) n.na.total/n.expected else NA_real_
      complete.wateryear <- n.present == n.expected
      qc.year <- is.finite(na.prct.yr) && na.prct.yr <= na_accept
      
      data.frame(
        wateryear=as.integer(yr),
        wy_start=wy.start,
        wy_end=wy.end,
        n_expected_yr=as.integer(n.expected),
        n_present_yr=as.integer(n.present),
        n_na_present_yr=as.integer(n.na.present),
        n_missing_structural_yr=as.integer(n.missing.structural),
        na_prct_yr=as.numeric(na.prct.yr),
        complete_wateryear=complete.wateryear,
        qc_year=qc.year
      )
    })
    
    year.qc <- data.table::rbindlist(year.qc, use.names=TRUE)
    
    if(isTRUE(apply_qc)){
      rejected <- year.qc[qc_year == FALSE]
      if(nrow(rejected) > 0){
        track.env$qc_removed_count <- track.env$qc_removed_count + nrow(rejected)
        for(i in seq_len(nrow(rejected))){
          pos <- length(track.env$qc_removed_details)+1L
          track.env$qc_removed_details[[pos]] <- sprintf(
            "Estação: %s | Ano: %s | Falhas: %.1f%%",
            gauge, rejected$wateryear[i], 100*rejected$na_prct_yr[i]
          )
        }
      }
    }
    
    imax.duration <- lapply(durations, function(d){
      
      k <- as.integer(round(d/d.base.hr))
      
      depth.agg <- RcppRoll::roll_sum(
        x=depths, n=k, na.rm=TRUE,
        fill=NA_real_, align="left"
      )
      
      na.window <- RcppRoll::roll_sum(
        x=as.integer(is.na(depths)), n=k,
        fill=NA_real_, align="left"
      )
      
      idx.ini <- seq_along(dates)
      idx.fim <- idx.ini + k - 1L
      valid.window <- idx.fim <= length(dates)
      
      date.fim <- rep(as.POSIXct(NA, tz=tz), length(dates))
      date.fim[valid.window] <- dates[idx.fim[valid.window]]
      
      all.na.window <- !is.na(na.window) & na.window == k
      
      span.h <- rep(NA_real_, length(dates))
      span.h[valid.window] <- as.numeric(
        difftime(
          date.fim[valid.window],
          dates[valid.window],
          units="hours"
        )
      )
      
      continuous.window <- rep(FALSE, length(dates))
      continuous.window[valid.window] <-
        abs(span.h[valid.window] - d.base.hr*(k-1L)) < 1e-8
      
      same.wateryear.window <- rep(FALSE, length(dates))
      same.wateryear.window[valid.window] <-
        wateryear[idx.ini[valid.window]] ==
        wateryear[idx.fim[valid.window]]
      
      depth.agg[
        all.na.window |
          !continuous.window |
          !same.wateryear.window
      ] <- NA_real_
      
      intensity.agg <- depth.agg/d
      
      imax.year <- lapply(years, function(yr){
        
        q <- year.qc[wateryear == yr]
        
        if(isTRUE(apply_qc) && !isTRUE(q$qc_year))
          return(NULL)
        
        year.idx <- wateryear == yr & mon.idx
        if(!any(year.idx)) return(NULL)
        
        intensity.yr <- intensity.agg[year.idx]
        date.yr <- dates[year.idx]
        depths.yr <- depths[year.idx]
        
        all.na.raw.year <- length(depths.yr) == 0 || all(is.na(depths.yr))
        all.na.intensity.year <- length(intensity.yr) == 0 || all(is.na(intensity.yr))
        
        if(all.na.raw.year || all.na.intensity.year){
          track.env$skipped_count <- track.env$skipped_count + 1L
          pos <- length(track.env$skipped_details)+1L
          track.env$skipped_details[[pos]] <- sprintf(
            "Estação: %s | Duração: %s h | Ano: %s",
            gauge, d, yr
          )
          return(NULL)
        }
        
        imax.yr <- max(intensity.yr, na.rm=TRUE)
        max.idx <- which(intensity.yr == imax.yr)[1]
        date.max <- date.yr[max.idx]
        n.max <- sum(intensity.yr == imax.yr, na.rm=TRUE)
        
        mon.max.idx <- wateryear == yr &
          mon.all == lubridate::month(date.max)
        
        na.prct.mon <- if(sum(mon.max.idx) > 0){
          sum(na.depths[mon.max.idx]) / sum(mon.max.idx)
        } else NA_real_
        
        data.frame(
          gauge_code=gauge,
          d=d,
          imax=imax.yr,
          date=date.max,
          wateryear=as.integer(yr),
          na_prct_yr=q$na_prct_yr,
          na_prct_mon=na.prct.mon,
          mon_filter=mon.filter,
          n_max=n.max,
          start_month=start.month,
          complete_wateryear=q$complete_wateryear,
          qc_year=q$qc_year,
          n_expected_yr=q$n_expected_yr,
          n_present_yr=q$n_present_yr,
          n_na_present_yr=q$n_na_present_yr,
          n_missing_structural_yr=q$n_missing_structural_yr
        )
      })
      
      imax.year <- imax.year[
        !vapply(imax.year, is.null, logical(1))
      ]
      
      if(!length(imax.year)) return(NULL)
      
      data.table::rbindlist(
        imax.year,
        use.names=TRUE,
        fill=TRUE
      )
    })
    
    imax.duration <- imax.duration[
      !vapply(imax.duration, is.null, logical(1))
    ]
    
    if(!length(imax.duration)) return(NULL)
    
    data.table::rbindlist(
      imax.duration,
      use.names=TRUE,
      fill=TRUE
    )
  })
  
  message("\nConfira sempre 'na_prct_yr' e 'na_prct_mon'.")
  
  if(track.env$qc_removed_count > 0){
    message("Anos removidos pelo QC: ", track.env$qc_removed_count)
  }
  
  if(track.env$skipped_count > 0){
    message("Anos/durações sem máximo válido: ", track.env$skipped_count)
  }
  
  imax.gauge <- imax.gauge[
    !vapply(imax.gauge, is.null, logical(1))
  ]
  
  if(!length(imax.gauge))
    return(data.frame())
  
  result <- data.table::rbindlist(
    imax.gauge,
    use.names=TRUE,
    fill=TRUE
  )
  
  data.table::setDF(result)
  return(result)
}
