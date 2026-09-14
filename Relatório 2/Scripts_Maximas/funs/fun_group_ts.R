#' Agrupa lista de séries temporais por resolução (time_step)
#'
#' @param ls      Lista nomeada de data.frames (um por posto).
#' @param ts_name Nome da coluna com o passo temporal (em minutos).
#' @return Lista cujos elementos são sub-listas agrupadas por time_step.
fun_group_ts <- function(ls, ts_name = "time_step") {

  time_steps <- vapply(ls, function(df) {
    ts <- unique(df[[ts_name]])
    if (length(ts) > 1L) {
      warning("Data.frame com múltiplos '", ts_name,
              "' — usando o primeiro valor.")
    }
    ts[1L]
  }, numeric(1))

  split(ls, time_steps)
}
