# Função p/ preencher e filtrar as séries de precipitação subdiária a partir dos argumentos:

## na_accpet: limiar de tolerância de falhas; e
## min_years: um número mínimo de anos.

# Essa função usa a função fun_fill_dates() que realiza o preeenchimento das datas
# filtra o conjunto preenchido conforme as especificações desejadas

fun_filter_set <- function(data,
                           na_accept,
                           min_years,
                           aux_path = "scripts/funcoes/fun_fill_dates.R",
                           save_path = NULL
                           ){
  
  # Checagens
  if(!is.numeric(na_accept)){
    stop("Argumento 'na_accept' não foi definido")
  }
  
  if(!is.numeric(min_years)){
    stop("Argumento 'min_years' não foi definido")
  }
  
  
  
  # Pacotes
  if(!require(pacman)){
    message("Instalando pacote 'pacman'...")
    install.packages("pacman")
  }
  
  pacman::p_load(pacman, dplyr, arrow)
  
  # Ler função p/ preencher datas
  source(aux_path)
  
  # Preencher datas
  message("\nPreenchendo datas...")
  list_filled <- fun_fill_dates(df = data, col_names = names(data))
  data <- bind_rows(list_filled); gc()
  
  # Resumir e filtrar dados
  message("\nAplicando filtros...")
  valid_gauges <- data %>%
    mutate(year = lubridate::year(datetime)) %>%
    group_by(gauge_code, year) %>%
    summarise(fail_prct = sum(is.na(rain_mm))/n() * 100, .groups = "drop") %>% # calcular % falhas
    group_by(gauge_code) %>%
    summarise(valid_years = sum(fail_prct <= na_accept * 100)) %>% # nro anos que atende
    filter(valid_years >= min_years) %>%                           # o criterio de falhas
    pull(gauge_code); gc()                                         # extrair so codigo estacoes            
  
  df_valid <- data[data$gauge_code %in% valid_gauges,]
  
  # Filtrar estacoes
  if(is.null(save_path)){
    return(df_valid)
  } else{
    arrow::write_parquet(x = df_valid, sink = save_path)
  }
  
}
