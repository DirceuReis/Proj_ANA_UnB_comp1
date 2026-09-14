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
