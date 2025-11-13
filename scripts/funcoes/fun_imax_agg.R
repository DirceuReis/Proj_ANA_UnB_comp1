
fun_imax_agg <- function(data,             # lista contendo as séries das estações
                         durations,        # vetor c/ durações desde a duração base (em horas)
                         which.mon = 1:12, # subconjunto de meses se desejado
                         names = c("datetime", "rain_mm") # nome das colunas, nessa ordem
                         ){
  
  # Pacotes
  if(!require(pacman)){
    message("Instalando pacote 'pacman'...")
    install.packages("pacman")
  }
  pacman::p_load(pacman, RcppRoll, pbapply, lubridate)
  
  # Confere se o argumento 'data' é uma lista
  if(!inherits(data, "list")){
    stop("'data' deve ser uma lista, e não: ", class(data), ".\n")
  }
  
  # Vetor de estações
  gauges <- names(data)
  
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
      stop("Data.frame of station ", station, " does not contain $", names[1], " or $", names[2], ".")
    }
    
    # Processar uma estação
    mon.idx <- lubridate::month(data.gauge[[names[1]]]) %in% which.mon # TRUE/FALSE se mês está contindo em 'which.mon'
    data.gauge <- data.gauge[mon.idx,] # filtrar somente meses em 'which.mon'
    dates <- data.gauge[[names[1]]]    # extrair vetor c/ datas
    depths <- data.gauge[[names[2]]]   # extrair vetor c/ lâminas
    
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
                                      fill = NA,              # preencher as pontas c/ NA
                                      align = "left")         # alinhas janela no índice i (NAs no final)
      
      # Auxiliares
      years <- unique(lubridate::year(dates)) # vetor c/ valores únidos dos anos da estação 'gauge'
      
      # Converter soma acumulada em intensidade/hora
      intensity.agg <- depth.agg/d
      
      # Calcular máximos anuais
      imax.year <- lapply(X = years, FUN = function(yr){
        
        year.idx <- lubridate::year(dates) == yr # índices p/ filtrar anos
        
        intensity.yr <- intensity.agg[year.idx]  # filtrar intensidades em 'yr'
        date.yr <- dates[year.idx]               # filtrar datas em 'yr'
        na.yr <- sum(is.na(depths)[year.idx])    # NAs em 'yr'
        na.prct <- na.yr/length(intensity.yr)    # porcentagem de NAs no ano
        max.idx <- which.max(intensity.yr)       # posição do máximo
        n.max <- length(max.idx)                 # comprimento do 
        
        # Conferir se o ano inteiro é NA, senão extrair os máximos
        if(n.max <= 0){
          imax.yr <- as.numeric(NA)
          date.max <- as.numeric(NA)
        } else{
          imax.yr <- intensity.yr[[max.idx]][1] # extrair máximo anual
          date.max <- date.yr[[max.idx]][1]     # extrair data do máximo anual
        }
        
        out <- data.frame("gauge_code" = gauge, # código da estação
                          "d" = d,              # duração [h]
                          "imax" = imax.yr,     # intensidade máxima anual [mm/h]
                          "date" = date.max,    # data em que ocorre o máximo             
                          "year" = yr,          # ano em que ocorre o máximo (alterar dps p/ ano_hidro)
                          "na_prct" = na.prct,  # porcentagem de falhas por ano
                          "mon_filter" = deparse(which.mon), # meses usados p/ análise
                          "n_max" = n.max)      # número de máximos iguais
         
      }) # fim 'imax.year'
      
      # Transformar em  tbl_df
      imax.year <- do.call(rbind, imax.year)
      
    }) # fim 'imax.duration'
    
    # Transformar em tbl_df
    imax.duration <- do.call(rbind, imax.duration)
    
  }) # fim 'imax.gauge'
  
  return(do.call(rbind, imax.gauge))
  
}
