scatter <- function(pol, pts, cell.size){
  
  if(!any(c(inherits(pol, "sf"), inherits(pts, "sf")))) stop("Objetos 'pol' e 'pts' devem ser do tipo 'sf' e ter o mesmo CRS.")
  
  message("Tamanho da grade: ", cell.size[1], "° x ", cell.size[2], "°")
  
  # Gerar grade sobre um polígono sf
  grid <- sf::st_make_grid(x = pol, cellsize = cell.size, what = "polygons")
  sf.grid <- sf::st_sf(geometry = grid)
  sf.grid <- sf::st_intersection(sf.grid, pol)
  
  # Gerar pontos nos centróides da grade
  centroids <- sf::st_centroid(sf.grid)
  
  # Calcular distância entre os pontos e os centróides
  mat.dist <- sf::st_distance(x = centroids, y = pts)
  colnames(mat.dist) <- pts$gauge_code
  
  gauges <- mat.dist %>% 
    as_tibble() %>% 
    mutate(centroid = rownames(.)) %>% 
    pivot_longer(cols = !centroid,
                 names_to = "closest_gauge",
                 values_to = "distance") %>% 
    group_by(centroid) %>% 
    arrange(distance) %>% 
    slice(1) %>% 
    pull(closest_gauge) %>% 
    unique()
  
  message(length(gauges), " pontos selecionados")
  
  # return(pts[pts$gauge_code %in% gauges,])
  return(gauges)
  
}