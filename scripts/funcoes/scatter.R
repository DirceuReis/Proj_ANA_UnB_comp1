# https://gis.stackexchange.com/questions/20804/how-to-find-the-most-spread-out-locations

# Função p/ encontrar um subconjunto de pontos que melhor cobre uma área usando k-means
# p/ criar clusters e selecionando a estação que se encontra mais ao centro do cluster

scatter <- function(points, n.cluster){
  
  # Funções
  distance <- function(x,y) sqrt(sum((x - y)^2))
  f <- function(k) distance(x = centers[groups[k],], y = points[k,])
  
  # Calcular clusters
  clusters <- kmeans(points, n.cluster)
  groups <- clusters$cluster
  centers <- clusters$centers
  
  # Selecionar ponto mais próximo do centro do cluster
  epsilon <- sqrt(min(clusters$withinss))/1000 # medida de "proximidade"
  n <- dim(points)[1]
  radius <- apply(matrix(1:n), 1, f) + runif(n, max = epsilon)
  min.dist <- tapply(radius, groups, min)
  
  return(points[radius == min.dist[groups],])
  
}