# Função p/ preencher e filtrar as séries de precipitação subdiária a partir dos argumentos:

## na_accpet: limiar de tolerância de falhas; e
## min_years: um número mínimo de anos.

# Essa função usa a função fun_fill_dates() que realiza o preeenchimento das datas
# filtra o conjunto preenchido conforme as especificações desejadas

fun_filter_set <- function(data,
                           daily = FALSE,
                           col_names = c("gauge_code", "rain_mm", "datetime", "time_step", "responsible"),
                           filter = FALSE,
                           na_accept = 0,
                           min_years = 0,
                           save_path = NULL
                           ){
  
  # Pacotes
  if(!require("pacman")) install.packages("pacman")
  pacman::p_load(tidyverse, fastmatch, beepr, pbapply, arrow)
  
  # Checagens preliminares
  # Check if df is a data.frame
  if(!inherits(data, c("data.frame", "data.table", "tibble"))){
    stop("\nArgument 'data' must be a data.frame (data.table or tibble), not a ", class(data), ".")
  }
  
  # Check if names(df) is the same as 'col.names' argument
  # Essa forma não confere se a ordem é a mesma
  if(sum(is.element(col_names[1:5], names(data))) != 5){
    stop("\nColumns of 'data' are not the same as:\n", paste(col_names, collapse = ", "))
  }
  
  # Check if 'datetime' column is of POSIX class
  if(!is.POSIXt(data[[col_names[3]]])){
    warning("\nReminder: 'time_step' must be in minutes.")
    stop("Column three must contain date information in POSIXt (ct or lt) format, instead it's ", class(col_names[3]), ".")
  }
  
  
  fun_fill_dates <- function(df, daily = FALSE, col_names){
    
    # Function to fill gauges time series based on their recording time step
    fun_fill <- function(df){
      
      # Extract gauge code and responsible based on the first observation
      gauge_code <- df[[col_names[1]]][1]
      resp <- df[[col_names[5]]][1]
      time_step <- unique(df[[col_names[4]]])*60 # convert to seconds to build new date sequence
      
      # Build new date sequence
      if(isFALSE(daily)){
        
        # Extract start and end dates
        # Adicionei 'na.rm = TRUE' nas funções min e max pq algumas estações do SGB estão com datetime == NA
        date_min <- lubridate::floor_date(min(df[[col_names[3]]], na.rm = TRUE), unit = "year")
        date_max <- lubridate::ceiling_date(max(df[[col_names[3]]], na.rm = TRUE), unit = "year") - 1
        
        if(!any(is.finite(c(date_min, date_max, time_step)))){
          warning("time_step: ", time_step)
          warning("date_min: ", date_min)
          warning("date_max: ", date_max)
          stop("\nSome values of 'datetime' or 'time_step' in gauge ", gauge_code, " (", resp, ")", " dates are NA.\n")
        }
        
        # Extract time_steps vector
        seq_date <- seq(from = date_min, to = date_max, by = time_step)
        
        # Fill new rainfall sequence by matching date between new sequence and df$datetime
        seq_rain <- df[[col_names[2]]][fastmatch::fmatch(x = seq_date, table = df[[col_names[3]]])]
        
      } else{
        
        # Extract start and end dates
        date <- lubridate::date(df[[col_names[3]]]) # remover 'time'
        date_min <- lubridate::floor_date(min(date), unit = "year")
        date_max <- lubridate::ceiling_date(max(date), unit = "year") - 1
        seq_date <- seq(from = date_min, to = date_max, by = "day")
        
        # Fill new rainfall sequence by matching date between new sequence and date(df$datetime)
        seq_rain <- df[[col_names[2]]][fastmatch::fmatch(x = seq_date, table = date)]
        
      }
      
      # Build final table
      df <- tibble(gauge_code = gauge_code,
                   rain_mm = seq_rain,
                   datetime = seq_date,
                   time_step = time_step,
                   responsible = resp)
      
      return(df)
      
    }
    
    # Função adaptada pra preencher completar as datas das estações do CEMADEN
    fun_fill_cemaden <- function(df){
      
      # Conferir resolução da estação
      time_step <- unique(df[[col_names[4]]])*60 # resolução padrão em segundos
      if(length(time_step) > 1){                 # conferir se há somente 1 passo de tempo padrão
        stop("Estação deve ter somente 1 único 'time_step' padrão, atualmente possui: ", length(time_step))}
      
      # Extrair gauge_code e responsible baseado na primeira observação
      gauge_code <- df[[col_names[1]]][1]
      resp <- df[[col_names[5]]][1]
      
      # Extrair datas de início e fim da série
      date_min <- lubridate::floor_date(min(df[[col_names[3]]]), unit = "year")
      date_max <- lubridate::ceiling_date(max(df[[col_names[3]]]), unit = "year") - 1
      
      # Construir nova sequência de datas
      df$delta <- difftime(df[[col_names[3]]], dplyr::lag(df[[col_names[3]]]), units = "mins")     # calcular deltas
      seq_date <- seq(from = date_min, to = date_max, by = time_step)                              # datas
      seq_rain <- df[[col_names[2]]][fastmatch::fmatch(x = seq_date, table = df[[col_names[3]]])]  # buscar rain_mm
      seq_delta <- df$delta[fastmatch::fmatch(x = seq_date, table = df[[col_names[3]]])] # buscar deltas
      
      # A nova sequência de deltas é importante p/ determinar a posição dos dados que devem ser ajustados
      # Construir tabela final (ainda vai ser corrigida)
      df <- tibble(gauge_code,
                   rain_mm = seq_rain,
                   datetime = seq_date,
                   time_step = time_step,
                   responsible = resp,
                   delta = seq_delta)
      
      # Corrigir observações c/ 10 < delta <= 60 min (atualmente estão c/ NA)
      # começando da segunda observação pq a primeira "delta" == NA
      # de forma vetorizada
      delta <- df$delta
      condition <- !is.na(delta) & delta > 10 & delta <= 60 # vetor de TRUE/FALSE
      valid <- which(condition)                                      # indices de condition == TRUE
      valid <- valid[valid >= 2]                                     # começa a partir da segunda linha
      interval <- as.integer((delta[valid]/10)) - 1
      
      # Encontrar quais índices devem ser corrigidos (zerados)
      zeros <- unlist(mapply(function(i, n){
        if(n > 0) seq(i - 1, max(1, i - n), -1)
        else interger(0)
      },  valid, interval))
      
      df[[col_names[2]]][zeros] <- 0
      
      # Reestabelecer nome original das colunas
      df <- df[-6]           # remover última coluna (delta)
      names(df) <- col_names # reeestabelecer nomes
      return(df)
      
    }
    
    # Separar os conjuntos
    df_others <- df[df[[col_names[5]]] != "CEMADEN",]  # ICEA, INMET, ANA
    df_cemaden <- df[df[[col_names[5]]] == "CEMADEN",] # CEMADEN
    
    # Outras instituições
    other_resp <- unique(df_others[[col_names[5]]])
    message("\nEstações ", paste(other_resp, collapse = ", "), "...")
    list_others <- split(x = df_others, f = df_others[[col_names[[1]]]])
    list_filled_others <- pbapply::pblapply(list_others, FUN = fun_fill)
    
    # CEMADEN
    message("\nEstações CEMADEN...")
    list_cemaden <- split(df_cemaden, f = df_cemaden[[col_names[[1]]]])
    list_filled_cemaden <- pbapply::pblapply(list_cemaden, FUN = fun_fill_cemaden)
    
    # Combinar listas
    list_filled <- c(list_filled_cemaden, list_filled_others)
    
    # Sound alert
    beepr::beep(sound = 10)
    invisible(gc())
    
    return(list_filled)
    
  }
  
  # Preencher datas
  message("\nPreenchendo datas...")
  list_filled <- fun_fill_dates(df = data, daily = FALSE, col_names = col_names)
  data <- bind_rows(list_filled)
  
  # Filtrar estações
  if(isTRUE(filter)){
    
    # Filtro: percentual máximo de falhas anual
    if(!is.numeric(na_accept)){
      stop("Argumento 'na_accept' não foi definido")
    }
    
    # Filtro: nro. mínimo de anos
    if(!is.numeric(min_years)){
      stop("Argumento 'min_years' não foi definido")
    }
    
    # Resumir e filtrar dados
    message("\nAplicando filtros...")
    data <- data %>% 
      mutate(year = lubridate::year(datetime)) %>% 
      group_by(gauge_code, year) %>% 
      mutate(na_prct = sum(is.na(rain_mm))/n()) %>% 
      ungroup() %>% 
      # filter(if_else(year == 2024, TRUE, na_prct <= na_accept)) %>% 
      filter(na_prct <= na_accept) %>%
      group_by(gauge_code) %>% 
      filter(n_distinct(year) >= min_years) %>% 
      ungroup()
    
    
  }
  
  # Filtrar estacoes
  if(is.null(save_path)){
    return(data)
  } else{
    arrow::write_parquet(x = data, sink = save_path)
  }
  
}
