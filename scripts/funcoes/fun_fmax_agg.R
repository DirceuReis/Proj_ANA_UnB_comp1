fun_fmax_agg <- function(data,
                         durations,
                         which.mon = 1:12,
                         names = c("datetime", "rain_mm"),
                         tz = "UTC",
                         strict_complete_blocks = TRUE,
                         anchor = c("year_start", "month_start", "day_start")) {
  
  anchor <- match.arg(anchor)
  
  if(!requireNamespace("pbapply", quietly = TRUE)) install.packages("pbapply")
  if(!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
  
  gauges <- names(data)
  
  out_gauge <- pbapply::pblapply(gauges, function(gauge){
    
    df <- data[[gauge]]
    if(!is.data.frame(df)) stop("Elemento ", gauge, " não é data.frame.")
    if(sum(is.element(names[1:2], names(df))) != 2){
      stop("Data.frame da estação ", gauge, " não contém ", names[1], " ou ", names[2], ".")
    }
    
    # ordenar + tz
    df <- df[order(df[[names[1]]]), ]
    dt <- as.POSIXct(df[[names[1]]], tz = tz)
    x  <- df[[names[2]]]
    
    # filtrar meses
    mon <- lubridate::month(dt)
    keep <- mon %in% which.mon
    dt <- dt[keep]
    x  <- x[keep]
    
    n <- length(dt)
    if(n < 2) return(NULL)
    
    # base step em segundos (moda)
    dsec_all <- diff(as.numeric(dt))
    dsec_all <- dsec_all[is.finite(dsec_all) & dsec_all > 0]
    if(length(dsec_all) == 0) stop("Não foi possível inferir resolução base na estação ", gauge)
    
    base_sec <- as.integer(names(sort(table(as.integer(round(dsec_all))), decreasing = TRUE))[1])
    
    dur_sec <- as.integer(round(durations * 3600))
    bad <- (dur_sec %% base_sec) != 0
    if(any(bad)){
      stop("Durações não-múltiplas do passo base (~", base_sec/60, " min) na estação ", gauge, ": ",
           paste0(durations[bad], collapse = ", "))
    }
    
    # trabalhar com segundos
    tnum <- as.numeric(dt)
    
    # ano/mês/dia como inteiros
    yr <- lubridate::year(dt)
    mo <- mon
    
    # âncora em segundos, vetorizada e feita 1 vez
    if(anchor == "year_start"){
      # 01/01 00:00 de cada ano
      anchor_time <- as.POSIXct(paste0(yr, "-01-01 00:00:00"), tz = tz)
    } else if(anchor == "month_start"){
      anchor_time <- as.POSIXct(paste0(yr, "-", sprintf("%02d", mo), "-01 00:00:00"), tz = tz)
    } else { # day_start
      anchor_time <- as.POSIXct(format(dt, "%Y-%m-%d 00:00:00"), tz = tz)
    }
    anchor_num <- as.numeric(anchor_time)
    
    # anos únicos
    years <- sort(unique(yr))
    
    res_all_d <- vector("list", length(durations))
    
    for(i in seq_along(durations)){
      d_hr <- durations[i]
      bs   <- dur_sec[i]             # tamanho do bloco em segundos
      
      # id do bloco fixo dentro da âncora
      block_id <- floor((tnum - anchor_num) / bs)
      
      # identificador "ano|bloco" para agregar
      key <- paste0(yr, "|", block_id)
      
      # soma por bloco (mm)
      sum_block <- tapply(x, key, function(v){
        if(strict_complete_blocks && anyNA(v)) return(NA_real_)
        sum(v, na.rm = TRUE)
      })
      
      # data de início do bloco: pegar o menor timestamp do bloco
      start_block <- tapply(tnum, key, min)
      
      # agora pegar, para cada ano, o maior bloco (ignorando NA)
      # separar a parte do ano da key sem reconstruir tabelas
      key_year <- as.integer(sub("\\|.*$", "", names(sum_block)))
      
      # iterar por ano
      out_year <- vector("list", length(years))
      
      for(j in seq_along(years)){
        y <- years[j]
        idx_y <- which(key_year == y)
        if(length(idx_y) == 0){
          out_year[[j]] <- data.frame(
            gauge_code = gauge, d = d_hr, fmax = NA_real_,
            date = as.POSIXct(NA, tz = tz), year = y,
            stringsAsFactors = FALSE
          )
          next
        }
        
        sb <- sum_block[idx_y]
        if(all(is.na(sb))){
          out_year[[j]] <- data.frame(
            gauge_code = gauge, d = d_hr, fmax = NA_real_,
            date = as.POSIXct(NA, tz = tz), year = y,
            stringsAsFactors = FALSE
          )
          next
        }
        
        # máximo do acumulado (mm)
        sb_safe <- replace(sb, is.na(sb), -Inf)
        kmax <- which.max(sb_safe)
        
        fmax_depth <- sb[kmax]
        date_start <- as.POSIXct(start_block[idx_y][kmax], origin = "1970-01-01", tz = tz)
        
        out_year[[j]] <- data.frame(
          gauge_code = gauge,
          d = d_hr,
          fmax = fmax_depth / d_hr,   # mm/h
          date = date_start,
          year = y,
          stringsAsFactors = FALSE
        )
      }
      
      res_all_d[[i]] <- do.call(rbind, out_year)
    }
    
    do.call(rbind, res_all_d)
  })
  
  do.call(rbind, out_gauge)
}