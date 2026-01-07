# Calcular distância entre estações diárias e subdiárias
fun_close_gauges <- function(sf.x, sf.y, max.distance = as.numeric(NULL), gauge.name = "gauge_code"){
  
  # Pacotes
  if(!require(sf)){
    message("Instalando gerenciador de pacotes 'sf'...")
    install.packages("sf")
  }
  
  # Conferir se objetos são do tipo 'sf'
  if(!inherits(sf.x, "sf")) stop("'sf.x' must be and 'sf' object.")
  if(!inherits(sf.y, "sf")) stop("'sf.y' must be and 'sf' object.")
  
  # Conferir se objetos têm estão na mesma coordenada de referência
  if(sf::st_crs(sf.x) != sf::st_crs(sf.y)) stop("'sf.x' and 'sf.y' must have the same CRS.")
  
  # Calcular distância entre estações
  distance <- as.data.frame(sf::st_distance(sf.x, sf.y)) # matriz de distâncias
  
  # Extrair nomes
  gauges.subdaily <- sf.x[[gauge.name]]
  gauges.daily <- sf.y[[gauge.name]]
  rownames(distance) <- gauges.subdaily # nome das colunas
  colnames(distance) <- gauges.daily    # nome das linhas
  
  # Organizar dados
  distance <- lapply(X = gauges.daily, FUN = function(daily){
    
    subdaily.distance <- as.numeric(distance[[daily]]) # extrair vetor de distâncias
    distances <-                                       # dataframe:
      data.frame(gauge_daily = daily,                  # código estação diária
                 gauge_subdaily = gauges.subdaily,     # código estação subdiária
                 distance_m = subdaily.distance)       # distância entre estações [m]
    
  }) # fim 'distance'
  
  # Transformar em 'tbl_df'
  distance <- do.call(rbind, distance)
  
  # Filtrar se necessário
  if(!is.null(max.distance)) distance <- distance[distance$distance_m <= max.distance,]
  
  return(distance)
  
} # fim 'fun_close_gauges'