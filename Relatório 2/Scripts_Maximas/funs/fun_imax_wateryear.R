#' Calcula intensidades máximas anuais por duração e ano hidrológico
#'
#' Percorre uma lista nomeada de séries temporais sub-diárias (um data.frame
#' por posto), aplica uma soma móvel de comprimento variável (RcppRoll) e
#' extrai o máximo de intensidade (mm/h) em cada ano hidrológico.
#'
#' Convenção de ano hidrológico adotada:
#'   Se m0 = 10 (outubro), então outubro/2020 pertence ao ano hidrológico
#'   2021 (recebe year + 1). O rótulo do ano hidrológico corresponde ao ano
#'   civil em que o período TERMINA.
#'   Exemplo: out/2020 – set/2021 → wateryear = 2021.
#'
#' @param data        Lista nomeada de data.frames (um por posto).
#' @param durations   Vetor numérico de durações em HORAS (ex.: c(1, 6, 24)).
#' @param start_month Inteiro único (mesmo mês para todos os postos) OU lista/
#'                    vetor nomeado pelo código do posto (mês por posto, 1–12).
#' @param which.mon   Inteiro(s) indicando os meses a incluir na busca do
#'                    máximo (default: 1:12, ou seja, o ano inteiro).
#'                    Útil para restringir a estação chuvosa.
#' @param names       Character(2): nomes das colunas de data/hora e
#'                    precipitação no data.frame de cada posto.
#'
#' @return data.frame com colunas:
#'   gauge_code, d, imax, date, wateryear, start_month,
#'   na_prct_yr, na_prct_mon, mon_filter, n_max.
#'   Retorna data.frame() vazio se nenhum posto produzir resultado.
#'
#' @details
#'   - Pacotes requeridos (devem estar instalados): RcppRoll, pbapply, lubridate.
#'   - Janelas com descontinuidade temporal são marcadas como NA antes de
#'     buscar o máximo, evitando somas espúrias sobre lacunas.
#'   - A resolução base é estimada pela MEDIANA dos intervalos consecutivos,
#'     tornando o cálculo robusto a eventuais falhas no início da série.
#'   - O contador `skipped` usa <<- (closure) — funciona corretamente para
#'     uso serial; NÃO paralelizar com mclapply/future sem ajuste.
fun_imax_wateryear <- function(data,
                               durations,
                               start_month,
                               which.mon = 1:12,
                               names = c("datetime", "rain_mm")) {

  # --- dependências (devem ser carregadas pelo script chamador) -------------
  if (!requireNamespace("RcppRoll",  quietly = TRUE)) stop("Instale RcppRoll")
  if (!requireNamespace("pbapply",   quietly = TRUE)) stop("Instale pbapply")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Instale lubridate")

  # --- validações de entrada ------------------------------------------------
  if (!inherits(data, "list")) {
    stop("'data' deve ser uma lista, e não: ", class(data), ".")
  }

  gauges <- base::names(data)

  if (is.null(gauges) || any(!nzchar(gauges))) {
    stop("'data' deve ser uma lista nomeada pelos códigos das estações.")
  }

  # --- auxiliar: mês inicial por posto -------------------------------------
  get.start.month <- function(gauge) {
    if (length(start_month) == 1L && is.null(base::names(start_month))) {
      return(as.integer(start_month[[1L]]))
    }
    if (!is.null(base::names(start_month)) &&
        gauge %in% base::names(start_month)) {
      return(as.integer(start_month[[gauge]]))
    }
    return(NA_integer_)
  }

  # --- contador de anos 100 % NA (closure serial) --------------------------
  skipped <- 0L

  # --- loop principal por posto --------------------------------------------
  imax.gauge <- pbapply::pblapply(X = gauges, FUN = function(gauge) {

    data.gauge <- data[[gauge]]

    if (!is.data.frame(data.gauge)) {
      stop("Elemento ", gauge, " não é data.frame: ", class(data.gauge))
    }
    if (sum(is.element(names[1:2], base::names(data.gauge))) != 2L) {
      stop("Posto ", gauge, " não contém colunas '",
           names[1], "' ou '", names[2], "'.")
    }

    m0 <- get.start.month(gauge)
    if (!is.finite(m0) || m0 < 1L || m0 > 12L) {
      warning("Posto ", gauge, " sem start_month válido — pulado.")
      return(NULL)
    }

    # Ordenar cronologicamente
    data.gauge <- data.gauge[order(data.gauge[[names[1]]]), ]

    dates  <- data.gauge[[names[1]]]
    depths <- as.numeric(data.gauge[[names[2]]])

    if (length(dates) < 2L) {
      warning("Posto ", gauge, " possui menos de duas observações.")
      return(NULL)
    }

    # Resolução temporal base — mediana dos intervalos (robusto a lacunas
    # no início da série)
    dt.hr <- as.numeric(diff(dates), units = "hours")
    d.base.hr <- median(dt.hr, na.rm = TRUE)

    if (!is.finite(d.base.hr) || d.base.hr <= 0) {
      warning("Resolução temporal inválida no posto ", gauge, ".")
      return(NULL)
    }

    # Alertar sobre irregularidades (não interrompe o processamento)
    if (any(abs(dt.hr - d.base.hr) > 1e-6)) {
      warning("Posto ", gauge, " possui intervalos temporais irregulares.")
    }

    months    <- lubridate::month(dates)
    na.depths <- is.na(depths)

    # Índices dos meses incluídos na análise
    mon.idx <- months %in% which.mon

    # Ano hidrológico: rótulo = ano civil em que o período TERMINA
    # (outubro/2020 com m0=10 → wateryear = 2021)
    wateryear <- lubridate::year(dates) +
      ifelse(m0 > 1L & months >= m0, 1L, 0L)

    years <- sort(unique(wateryear[mon.idx]))

    # Durações compatíveis com a resolução base
    ok.duration <- abs((durations / d.base.hr) %% 1) <= 1e-6
    if (any(!ok.duration)) {
      warning("Posto ", gauge, ": durações ignoradas (não múltiplas de ",
              round(d.base.hr * 60), " min): ",
              paste(durations[!ok.duration], collapse = ", "), " h.")
    }
    durations.gauge <- durations[ok.duration]
    if (length(durations.gauge) == 0L) return(NULL)

    # --- loop por duração --------------------------------------------------
    imax.duration <- lapply(X = durations.gauge, FUN = function(d) {

      n.window <- round(d / d.base.hr)  # número de registros na janela

      # Soma móvel sobre a série COMPLETA (align = "left": date retornada
      # é o instante de INÍCIO da janela máxima)
      depth.agg <- RcppRoll::roll_sum(
        x     = depths,
        n     = n.window,
        na.rm = TRUE,
        fill  = NA_real_,
        align = "left"
      )

      # Invalida janelas com descontinuidade temporal (ex.: falhas longas)
      if (n.window > 1L) {
        window.span <- c(
          as.numeric(
            difftime(
              dates[n.window:length(dates)],
              dates[1L:(length(dates) - n.window + 1L)],
              units = "hours"
            )
          ),
          rep(NA_real_, n.window - 1L)
        )
        expected.span <- (n.window - 1L) * d.base.hr
        depth.agg[abs(window.span - expected.span) > 1e-6] <- NA_real_
      }

      # Intensidade (mm/h)
      intensity.agg <- depth.agg / d

      # --- loop por ano hidrológico ----------------------------------------
      imax.year <- lapply(X = years, FUN = function(yr) {

        year.idx     <- wateryear == yr & mon.idx
        intensity.yr <- intensity.agg[year.idx]
        date.yr      <- dates[year.idx]

        # Percentual de NA no período do ano hidrológico filtrado
        na.prct.yr <- sum(na.depths[year.idx]) / sum(year.idx)

        # Ano sem nenhum valor válido
        if (all(is.na(depths[year.idx])) || all(is.na(intensity.yr))) {
          skipped <<- skipped + 1L
          return(data.frame(
            gauge_code  = gauge,
            d           = d,
            imax        = NA_real_,
            date        = date.yr[NA_integer_],
            wateryear   = yr,
            start_month = m0,
            na_prct_yr  = na.prct.yr,
            na_prct_mon = NA_real_,
            mon_filter  = paste(which.mon, collapse = ","),
            n_max       = 0L
          ))
        }

        # Máximo — registra número de empates e usa primeira ocorrência
        max.idx  <- which(intensity.yr == max(intensity.yr, na.rm = TRUE))
        n.max    <- length(max.idx)
        max.idx  <- max.idx[1L]

        imax.yr  <- intensity.yr[[max.idx]]
        date.max <- date.yr[[max.idx]]

        # Percentual de NA no mês do máximo (diagnóstico de qualidade)
        mon.max     <- lubridate::month(date.max)
        mon.max.idx <- wateryear == yr & months == mon.max
        na.prct.mon <- if (any(mon.max.idx)) {
          sum(na.depths[mon.max.idx]) / sum(mon.max.idx)
        } else {
          NA_real_
        }

        data.frame(
          gauge_code  = gauge,
          d           = d,
          imax        = imax.yr,
          date        = date.max,
          wateryear   = yr,
          start_month = m0,
          na_prct_yr  = na.prct.yr,
          na_prct_mon = na.prct.mon,
          mon_filter  = paste(which.mon, collapse = ","),
          n_max       = n.max
        )

      }) # fim loop anos

      do.call(rbind, imax.year)

    }) # fim loop durações

    do.call(rbind, imax.duration)

  }) # fim loop postos

  if (skipped > 0L) {
    message("Anos sem máximos válidos (posto/duração): ", skipped)
  }

  # Remove postos que retornaram NULL
  imax.gauge <- imax.gauge[!vapply(imax.gauge, is.null, logical(1L))]

  if (length(imax.gauge) == 0L) return(data.frame())

  do.call(rbind, imax.gauge)
}





# Intensidades máximas anuais no ano hidrológico (por posto)
#
# Correções vs versão original:
#   - start_month por gauge (mapa nomeado), não which.mon global
#   - na_prct_mon usa mon.max.idx (bug corrigido)
#   - wateryear: se m > 1 e month >= m → year+1, senão year

imax_wateryear <- function(data,
                           durations,
                           start_month,
                           names = c("datetime", "rain_mm"),
                           which_mon = NULL) {
  stopifnot(inherits(data, "list"))
  if (!requireNamespace("RcppRoll", quietly = TRUE)) {
    stop("Instale RcppRoll")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Instale data.table")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Instale lubridate")
  }

  gauges <- names(data)
  if (is.null(gauges) || any(!nzchar(gauges))) {
    stop("'data' precisa ser lista nomeada por gauge_code")
  }

  # start_month: escalar, vetor nomeado, ou lista
  sm_lookup <- function(g) {
    if (length(start_month) == 1L && is.null(names(start_month))) {
      return(as.integer(start_month))
    }
    if (!is.null(names(start_month)) && g %in% names(start_month)) {
      return(as.integer(start_month[[g]]))
    }
    if (is.list(start_month) && g %in% names(start_month)) {
      return(as.integer(start_month[[g]]))
    }
    NA_integer_
  }

  out_list <- vector("list", length(gauges))
  skipped <- 0L

  for (i in seq_along(gauges)) {
    gauge <- gauges[[i]]
    df <- data[[gauge]]
    if (!is.data.frame(df)) {
      stop("Elemento ", gauge, " não é data.frame")
    }
    if (sum(is.element(names[1:2], names(df))) != 2) {
      stop("Estação ", gauge, " sem colunas ", names[1], "/", names[2])
    }

    m0 <- sm_lookup(gauge)
    if (!is.finite(m0) || m0 < 1L || m0 > 12L) {
      warning("Estação ", gauge, " sem start_month válido — pulada")
      next
    }

    # Meses do ano hidrológico: m0, m0+1, ..., m0-1 (circular)
    if (is.null(which_mon)) {
      which_mon_g <- ((m0 - 1L + 0:11) %% 12L) + 1L
    } else {
      which_mon_g <- as.integer(which_mon)
    }

    dates_all <- df[[names[1]]]
    depths_all <- as.numeric(df[[names[2]]])
    mon_all <- lubridate::month(dates_all)
    mon_idx <- mon_all %in% which_mon_g
    dates <- dates_all[mon_idx]
    depths <- depths_all[mon_idx]
    mon_dates <- mon_all[mon_idx]
    na_depths <- is.na(depths)

    if (length(dates) < 2L) {
      next
    }

    shift <- ifelse(m0 > 1L & mon_dates >= m0, 1L, 0L)
    wateryear <- lubridate::year(dates) + shift
    years <- sort(unique(wateryear))

    d_base_hr <- as.numeric(difftime(dates[2], dates[1], units = "hours"))
    if (!is.finite(d_base_hr) || d_base_hr <= 0) {
      warning("Resolução inválida em ", gauge)
      next
    }
    if (any(abs((durations / d_base_hr) %% 1) > 1e-8)) {
      # filtra só durações compatíveis em vez de parar
      ok_d <- abs((durations / d_base_hr) %% 1) <= 1e-8
      durations_g <- durations[ok_d]
    } else {
      durations_g <- durations
    }
    if (length(durations_g) == 0L) {
      next
    }

    imax_d <- vector("list", length(durations_g))
    for (j in seq_along(durations_g)) {
      d <- durations_g[[j]]
      n_win <- as.integer(round(d / d_base_hr))
      depth_agg <- RcppRoll::roll_sum(
        x = depths, n = n_win, na.rm = TRUE, fill = NA_real_, align = "left"
      ) #date retornada é a data/hora de início da janela máxima.
      intensity_agg <- depth_agg / d

      rows_yr <- vector("list", length(years))
      for (k in seq_along(years)) {
        yr <- years[[k]]
        year_idx <- wateryear == yr
        if (!any(year_idx)) {
          next
        }
        if (all(is.na(depths[year_idx]))) {
          skipped <- skipped + 1L
          next
        }
        intensity_yr <- intensity_agg[year_idx]
        if (all(is.na(intensity_yr))) {
          skipped <- skipped + 1L
          next
        }
        date_yr <- dates[year_idx]
        na_prct_yr <- sum(is.na(depths[year_idx])) / sum(year_idx)
        max_idx <- which.max(intensity_yr)
        if (length(max_idx) < 1L || !is.finite(max_idx[1])) {
          next
        }
        max_idx <- max_idx[1]
        imax_yr <- intensity_yr[[max_idx]]
        date_max <- date_yr[[max_idx]]
        mon_max_idx <- year_idx & (mon_dates == lubridate::month(date_max))
        na_prct_mon <- if (any(mon_max_idx)) {
          sum(na_depths[mon_max_idx]) / sum(mon_max_idx)
        } else {
          NA_real_
        }

        rows_yr[[k]] <- data.frame(
          gauge_code = gauge,
          d = d,
          imax = imax_yr,
          date = date_max,
          wateryear = yr,
          start_month = m0,
          na_prct_yr = na_prct_yr,
          na_prct_mon = na_prct_mon,
          stringsAsFactors = FALSE
        )
      }
      imax_d[[j]] <- data.table::rbindlist(rows_yr, use.names = TRUE)
    }
    out_list[[i]] <- data.table::rbindlist(imax_d, use.names = TRUE)
  }

  if (skipped > 0L) {
    message("Anos 100% NA ignorados (estação/duração): ", skipped)
  }
  result <- data.table::rbindlist(out_list, use.names = TRUE)
  data.table::setDF(result)
}
