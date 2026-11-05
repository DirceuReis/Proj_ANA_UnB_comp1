fill.dates <- function(df, # data.frame with all stations
                           daily = FALSE,
                           col_names = c("gauge_code", "rain_mm", "datetime", "time_step")){
  
  # Check if df is a data.frame
  if(!inherits(df, c("data.frame", "data.table", "tibble"))){
    stop("Argument 'df' must be a data.frame (data.table or tibble), not a ", class(df), ".")
  }
  
  # Check if names(df) is the same as 'col.names' argument
  # Essa forma não confere se a ordem é a mesma
  if(sum(is.element(col_names[1:4], names(df))) != 4){
    stop("Columns of 'df' are not the same as:\n", paste(col_names[1], col_names[2], col_names[3], col_names[4]))
  }
  
  # Check if 'datetime' column is of POSIX class
  if(!is.POSIXt(df[[col_names[3]]])){
    warning("Reminder: 'time_steps' must be in minutes.")
    stop("Column three must contain date information in POSIXt (ct or lt) format, instead it's ", class(col_names[3]), ".")
  }
  
  # Packages
  if(!require("pacman")) install.packages("pacman")
  pacman::p_load(tidyverse, lubridate, fastmatch, beepr, pbapply, parallel)
  
  # Function to fill gauges time series based on their recording time step
  fill <- function(df){
    
    # Extract gauge code and responsible based on the first observation
    gauge_code <- df[[col_names[1]]][1]
    time_step <- unique(df[[col_names[4]]])*60 # convert to seconds to build new date sequence
    
    # Build new date sequence
    if(isFALSE(daily)){
      
      # Extract start and end dates
      date_min <- lubridate::floor_date(min(df[[col_names[3]]]), unit = "year")
      date_max <- lubridate::ceiling_date(max(df[[col_names[3]]]), unit = "year") - 1
      
      # Extract time_steps vector
      # seq_date <- seq(from = date_min, to = date_max, by = time_step)
      seq_date <- seq(from = date_min, to = date_max, by = "day")
      
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
                 time_step = time_step)
    
    return(df)
    
  }
  
  # Outras instituições
  message("\nEstações ICEA, ANA, INMET...")
  ls <- split(x = df, f = df[[col_names[[1]]]])
  ls_filled <- pbapply::pblapply(ls, FUN = fill)
  
  # Sound alert
  beepr::beep(sound = 10)
  gc()
  
  return(ls_filled)
  
}
