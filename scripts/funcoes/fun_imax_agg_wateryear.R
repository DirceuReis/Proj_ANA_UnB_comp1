
# FUNÇÃO ------------------------------------------------------------------

# Função para extrair séries de intensidades máximas anuais p/
# durações de interesse baseada no ano hidrológico informado

imax.wateryear <- function(data,             # lista contendo as séries das estações
                           durations,        # vetor c/ durações desde a duração base (em horas)
                           which.mon = 1:12, # subconjunto de meses se desejado
                           names = c("datetime", "rain_mm") # nome das colunas, nessa ordem
                           ){
  
  # Pacotes
  if(!require(pacman)){
    message("Instalando pacote 'pacman'...")
    install.packages("pacman")
  }
  pacman::p_load(pacman, RcppRoll, pbapply, lubridate, data.table)
  
  # Confere se o argumento 'data' é uma lista
  if(!inherits(data, "list")){
    stop("'data' deve ser uma lista, e não: ", class(data), ".\n")
  }
  
  message("Agregando e extraindo intensidades máximas anuais para durações: ", paste(durations, collapse = ", "))
  message("Início do ano hidrológico: ", which.mon[1])
  message("Meses analisados: ", paste(sort(which.mon), collapse = ", "))
  
  # Track run info
  track.env <- new.env()
  track.env$skipped_count <- 0
  track.env$skipped_details <- list()
  
  # Vetor de estações e outras variáveis globais
  gauges <- names(data)
  start.month <- which.mon[1]
  mon.filter <- deparse(which.mon)
  
  # Iterar sobre cada estação
  imax.gauge <- pbapply::pblapply(X = gauges, FUN = function(gauge){
    
    # Extrair data.frame da estação 'gauge'
    data.gauge <- data[[gauge]]
    
    # Confere se a lista contém data.frames
    if(!is.data.frame(data.gauge)){
      stop("Elements of 'data' must be data.frames. But element ", gauge, " contains: ", class(data.gauge))
    }
    
    # Usa is.element() parar conferir se os nomes do argumento 'names' são os mesmos nomes das colunas dos data.frames
    if(sum(is.element(names[1:2], names(data.gauge))) != 2){
      stop("Data.frame of station ", gauge, " does not contain $", names[1], " or $", names[2], ".")
    }
    
    # Processar uma estação
    mon.all <- lubridate::month(data.gauge[[names[1]]]) # calcular lubridate::month() só uma vez por estação
    mon.idx <- mon.all %in% which.mon  # TRUE/FALSE se mês está contindo em 'which.mon'
    data.gauge <- data.gauge[mon.idx,] # filtrar somente meses em 'which.mon'
    dates <- data.gauge[[names[1]]]    # extrair vetor c/ datas
    depths <- data.gauge[[names[2]]]   # extrair vetor c/ lâminas
    na.depths <- is.na(depths)         # pré-calcular NAs
    mon.dates <- mon.all[mon.idx]      # vetor de meses válidos p/ evitar recalcular posteriormente
    # shift <- ifelse(lubridate::month(dates) < start.month, 0, 1) # calcular deslocamento do ano hidrológico
    shift <- ifelse(start.month > 1 & lubridate::month(dates) >= start.month, 1, 0) # define whether a month gets into the current or the next year
    wateryear <- lubridate::year(dates) + shift # vetor c/ anos hidrológicos ajustados
    years <- sort(unique(wateryear))   # vetor c/ valores únicos dos anos hidrológicos definidos de 'gauge'
    
    # Conferir se durações são múltiplas da resolução 'base' da estação
    d.base.hr <- as.numeric(dates[2] - dates[1], units = "hours") # diferença entre 2ª e 1ª 'datetime'
    if(any((durations/d.base.hr) %% 1 > 1e-8)){
      stop("Pelo menos uma das durações de agregação não é múltipla da resolução base de ", d.base.hr*60, " min na estação ", gauge, ".")
    }
    
    # Iterar sobre cada duração
    imax.duration <- lapply(X = durations, FUN = function(d){
      
      # Vetor de somas móveis de precipitação
      depth.agg <- RcppRoll::roll_sum(x = depths,             # vetor c/ dados de precipitação
                                      n = round(d/d.base.hr), # tamanho da janela móvel
                                      na.rm = TRUE,           # calcular a soma mesmo c/ NAs no vetor 'depths'
                                      fill = NA,              # preencher as pontas c/ NA
                                      align = "left")         # alinhas janela no índice i (NAs no final)
      
      # Auxiliares
      # years <- unique(lubridate::year(dates)) # vetor c/ valores únicos dos anos da estação 'gauge'
      
      # Converter soma acumulada em intensidade/hora
      intensity.agg <- depth.agg/d
      
      # Calcular máximos anuais
      imax.year <- lapply(X = years, FUN = function(yr){
        
        # year.idx <- lubridate::year(dates) + shift == yr # índices p/ filtrar anos
        year.idx <- wateryear == yr # índices p/ filtrar anos (já foi adicionado o 'shift' anteriormente)
        intensity.yr <- intensity.agg[year.idx]          # filtrar intensidades em 'yr'
        
        # Conferir se o ano inteiro é NA
        if(all(is.na(depths[year.idx]))){
          track.env$skipped_count <- track.env$skipped_count + 1
          skip.row <- length(track.env$skipped_details) + 1
          track.env$skipped_details[[skip.row]] <-  sprintf("Estação: %s | Duração: %s h | Ano: %s", gauge, d, yr)
          return(NULL) # pula o ano inteiro caso seja NA
        }
        
        date.yr <- dates[year.idx]               # filtrar datas em 'yr'
        na.yr <- sum(is.na(depths)[year.idx])    # NAs em 'yr'
        na.prct.yr <- na.yr/length(intensity.yr) # porcentagem de NAs no ano
        max.idx <- which.max(intensity.yr)       # posição do máximo
        # n.max <- length(max.idx)                 # número de máximos
        imax.yr <- intensity.yr[[max.idx]][1]    # extrair máximo anual
        date.max <- date.yr[[max.idx]][1]        # extrair data do máximo anual
        # mon.idx <- year.idx & (lubridate::month(dates) == lubridate::month(date.max)) # filtrar meses de date.max
        mon.max.idx <- year.idx & (mon.dates == lubridate::month(date.max)) # filtrar meses de date.max
        # na.prct.mon <- sum(is.na(depths[mon.idx]))/sum(mon.idx)
        na.prct.mon <- sum(na.depths[mon.max.idx])/sum(mon.idx)
        
        out <- data.frame(
          "gauge_code" = gauge,   # código da estação
          "d" = d,                # duração [h]
          "imax" = imax.yr,       # intensidade máxima anual [mm/h]
          "date" = date.max,      # data em que ocorre o máximo             
          "wateryear" = yr,       # ano hidrológico em que ocorre o máximo (alterar dps p/ ano_hidro)
          "na_prct_yr" = na.prct.yr,         # porcentagem de falhas por ano
          "na_prct_mon" = na.prct.mon,       # porcentagem de falhas no mês do máximo
          "mon_filter" = mon.filter  # meses usados p/ análise
          # "n_max" = n.max]
        )
         
      }) # fim 'imax.year'
      
      # Transformar em  tbl_df
      # imax.year <- do.call(rbind, imax.year)
      data.table::rbindlist(imax.year, use.names = TRUE)
      
    }) # fim 'imax.duration'
    
    # Transformar em tbl_df
    # imax.duration <- do.call(rbind, imax.duration)
    data.table::rbindlist(imax.duration, use.names = TRUE)
    
  }) # fim 'imax.gauge'
  
  # message("Lembrete:\nEsta função não realiza filtros de qualidade. Havendo anos inteiros com falhas, a função retornará 'imax == 0'.")
  # message("Caso 'na_prct_yr' retorne 1, o ano da série era originalmente 100% preeenchido por NAs.")
  message("Confira sempre as colunas 'na_prct_yr' e 'na_prct_mon'.")
  
  if(track.env$skipped_count > 0){
    message("---- Relatório anos removidos (100% NA) ----")
    message(sprintf("Total de ocorrências (estação/ano/duração) ignoradas: %d", track.env$skipped_count))
    
    details <- unlist(track.env$skipped_details)
    
    if(length(details) > 10){
      message(paste(details[1:10], collapse = "\n"))
      message(sprintf("... e mais %d casos ocultados.", length(details) - 10))
    } else{
      message(paste(details, collapse = "\n"))
    }
    
  }
  
  # Resultado
  result <- data.table::rbindlist(imax.gauge, use.names = TRUE)
  return(data.table::setDF(result))
  
}
