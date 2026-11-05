# Função p/ transformar um data.frame e variáveis preditoras p/ cada estação
# Essa função transforma a primeira coluna (coluna geralmente com o nome
# das estações) em uma coluna com 1s e permite escolher um conjunto de colunas
# seguintes para aplicar log e outro conjunto para centralizar
covariate.transform <- function(covariate.matrix, transform = c("none", "center", "scale", "log", "boxcox"),
                                except = NULL){
  
  transform <- match.arg(transform)
  
  message("Covariates: ", paste(colnames(covariate.matrix), collapse = ", "))
  
  if(any(colnames(covariate.matrix) == "ones")) stop("Avoid including 'ones' column before running the transformations")
  if(!is.numeric(covariate.matrix)) stop("All variables should be numeric")
  
  out <- switch (transform,
                 "none" = {
                   message("Transformation applied: none")
                   covariate.matrix
                 },
                 "center" = {
                   message("Transformation applied: center")
                   apply(covariate.matrix, 2, scale, center = TRUE, scale = FALSE)
                 },
                 "scale" = {
                   message("Transformation applied: scale")
                   apply(covariate.matrix, 2, scale, center = TRUE, scale = TRUE)
                 },
                 "log" = {
                   message("Transformation applied: natural logarithm")
                   apply(covariate.matrix, 2, log) # lat and long should be positive
                   # apply(covariate.matrix, 2, function(x) sign(x)*log(abs(x)))
                 },
                 "boxcox" = {
                   stop("Not yet implemented. Choose from: NULL, center, scale or log")
                 }
  )
  
  if(!is.null(except)){
    message("Following covariates removed from transformation: ", paste(except, collapse = ", "))
    out[,except] <- covariate.matrix[,except]
  }
  
  return(out) # transformed (or nor) covariates matrix
  
}


# Function to estimate GLS model parameters beta and model error variance (Stedinger and Tasker, 1985)
build.gls <- function(std.delta,  # initial model error variance
                      n.gauges,   # number of stations
                      scovmatrix, # sample covariance matrix of theta
                      theta,      # hydrologic variable of interest
                      X,          # explanatory variables matrix 
                      k           # number of GLS parameters (or explanatory variables) + 1
                      ){
  
  # Build diagonal matrix only with model error variances (sigma_delta^2*Im)
  m <- n.gauges
  diag.var.delta <- matrix(0, ncol = m, nrow = m)
  for(i in 1:m) diag.var.delta[i,i] <- std.delta^2
  # diag.var.delta <- std.delta^2*diag()
  
  # Matriz de covariância do erro total
  lambda <- diag.var.delta + scovmatrix
  
  # Calcular estimador de beta e variancia do erro do modelo sigma_delta^2
  inv.lambda <- solve(lambda)                                             # inversa da matriz de covariancia dos erros
  b <- solve(t(X) %*% inv.lambda %*% X) %*% t(X) %*% inv.lambda %*% theta # estimador de beta – Stedinger e Tasker (1985)
  fo <- t(theta - X %*% b) %*% inv.lambda %*% (theta - X %*% b) - (m - k) # função objetivo fo = 0
  fo <- as.numeric(fo)
  
  return(list(check.delta = fo,
              std.delta = std.delta,
              b = b,
              lambda = lambda,
              inv.lambda = inv.lambda))
  
}

# Function to assess parameter significance by comparing estimated values to their standard errors
model.sigf <- function(models, m, conf.level = 0.95){
  
  pmax <- log(nrow(models), 2) # total amount of covariates tested
  kmax <- pmax + 1            # total amount of possible parameters
  par.matrix <- as.matrix(models[, 1:kmax])
  se.matrix <- as.matrix(models[,  (kmax + 1):(2*kmax)])

  # Ignore first column (average model) for testing
  par.test <- par.matrix[, -1]
  se.test <- se.matrix[, -1]
  
  # Student's t-test
  k <- models[["k"]]         # number of parameters in each model (row)
  t.stat <- par.test/se.test # (par - 0)/se
  pvalue <- 2*pt(abs(t.stat), df = m - k, lower.tail = FALSE) # p-value = 2*P[T > |t|]
  
  # Identify non-significant parameters and keep good models
  fails <- !is.na(pvalue) & (pvalue >= (1 - conf.level)) # pvalue below sigf.level rejects H0 (beta = 0)
  good.rows <- rowSums(fails) == 0
  sigf.models <- models[good.rows,]
  
  message("Two-sided t-test yielded ", nrow(sigf.models), 
          " significant models at ", (1 - conf.level)*100, 
          "% level out of ", nrow(models))
  
  return(sigf.models)
  
}


# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, lmom, optimx, beepr, ggplot2, sf, parallel, moments, patchwork)

# QUANTIL -----------------------------------------------------------------

# GEV
fun.q.gev <- function(p,     # probabilidade de não excedência
                      param  # vetor c/ parâmetros c(csi, alpha, kappa)
){                 
  
  csi <- param[1]; alpha <- param[2]; kappa <- param[3]
  
  # Quantil
  xp <- csi + alpha/kappa*(1 - (-log(p))^kappa)
  
  return(xp)  # quantil explícito em função dos parâmetros
  # e da probabilidade de não excedência
  
}

# Gumbel
fun.q.gu <- function(p,     # probabilidade de não excedência
                     param  # vetor c/ parâmetros c(csi, alpha)
){                 
  
  csi <- param[1]; alpha <- param[2]
  
  # Quantil
  xp <- csi - alpha*log(-log(p))
  
  return(xp)  # quantil explícito em função dos parâmetros
  # e da probabilidade de não excedência
  
}


# MOM-L -------------------------------------------------------------------

# Calcula os parâmetros iniciais que serão otimizados pelo
# Estimador de Máxima Verossimilhança (MLE) e
# Máxima Verossimilhança Generalizada (GMLE)

fun.lmom <- function(xp,          # coluna do dataframe c/ máximos anuais
                     dist = "gev" # "gumbel" ou "gev"
){
  
  # Carrega o pacote lmom
  pacman::p_load(lmom, pacman)
  
  # Estimativas dos momentos-L não-tendenciosas da amostra
  est.lmom <- lmom::samlmu(x = xp, sort.data = TRUE)
  
  # testar κ
  if(dist == "gumbel"){
    
    param <- lmom::pelgum(est.lmom)
    param <- data.frame(csi = param[1], alpha = param[2], row.names = NULL)
    
  } else{
    
    param <- lmom::pelgev(est.lmom)
    param <- data.frame(csi = param[1], alpha = param[2], kappa = param[3], row.names = NULL)
    
  }
  
  return(param)
  
}


# ML ----------------------------------------------------------------------

# Estimador de Máxima Verossimilhança - GEV
fun.l.gev.old <- function(param0, # vetor c/ os parâmetros, 1. csi 2. alpha, 3. kappa
                      xp      # vetor c/ a série AM observada
){
  
  # Paramentros Iniciais
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  
  n <- length(xp)
  
  # Testar parâmetro de escala (deve ser maior que zero)
  if(alpha < 0){
    return(1e6)
  }
  
  # Testar se kappa é negativo (Fréchet)
  if(kappa < 0){
    
    # Testar limite inferior
    if(min(xp) < (csi + alpha/kappa)){
      
      ln.L <- 1e6
      
    } else{
      
      # Função de verossimilhança GEV
      ln.L <- -(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa)))
      
    }
    
  } else{ # kappa positivo (Weibull)
    
    # Testar limite superior
    if(max(xp) > (csi + alpha/kappa)){
      i
      ln.L <- 1e6
      
    } else{
      
      # Função de verossimilhança GEV
      ln.L <- -(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa)))
      
    }
  }
  
  return(ln.L)
  
}

# Forma alternativa alterando a checagem dos limites da distribuição
fun.l.gev <- function(param0, xp){
  
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  n <- length(xp)
  
  zi <- 1 - kappa/alpha*(xp - csi)
  
  # Fallbacks
  if(alpha <= 0) return(1e6)
  if(any(zi <= 1e-8)) return(1e6 + abs(min(zi))*1e4)
  
  # Negative log-likelihood
  term1 <- (1/kappa - 1)*log(zi)
  term2 <- zi^(1/kappa)
  nll <- -(-n*log(alpha) + sum(term1 - term2))
  
  return(nll)
  
}

# Estimador de Máxima Verossimilhança - Gumbel
fun.l.gu <- function(param0, # vetor c/ os parâmetros 1. csi 2. alpha
                     xp      # vetor c/ a série AM observada
){
  
  # Parâmetros iniciaiss
  csi <- param0[1]
  alpha <- param0[2]
  
  n <- length(xp)
  
  # Testar parâmetro de escala
  if (alpha < 0){
    
    ln.L <- 1e6
    
  } else{
    
    # Função de verossimilhança (Gumbel)
    y <- (xp - csi)/alpha
    ln.L <- -(-n*log(alpha) - sum(y + exp(-y)))
    
  }
  
  return(ln.L) # naturalmente, o sinal do ln é negativo, então no momento q inverte o sinal pra minimizar a função,
  # precisa tirar reverter o sinal do máximo depois
  
}

# GML ---------------------------------------------------------------------

# Estimador de Máxima Verossimilhança Generalizada - GEV
# GMLE - priori Beta
# Problemas: deve o valor de erro deve ser negativo pra que o sinal trocado o torne positivo durante a maximização caso haja erro
fun.gl.beta.gev.old <- function(param0,        # vetor c/ os parâmetros, 1. csi 2. alpha, 3. kappa
                            media.kappa,   # -0,10
                            desvpad.kappa, # 0,122
                            xp){
  
  # Paramentros Iniciais
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  
  n <- length(xp)
  var.kappa <- desvpad.kappa^2
  
  # Distribuicao Beta
  a <- 0.5
  b <- media.kappa+a
  p <- b^2*(1-b)/var.kappa-b
  q <- p*(1/b-1)
  
  # Função log-verossimilhança
  if(abs(kappa) >= a) return(1e6) # adicionei essa restrição p/ evitar NaNs em log(a + kappa)
  if(alpha < 0) return(1e6)
  
  if(kappa < 0){
    if(min(xp) < (csi + alpha/kappa)){
      ln.GL <- 1e6
    }else{
      ln.GL <- (-(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa)) +
                    (p - 1)*log(a + kappa) + (q - 1)*log(a - kappa)))
    }
  }else{
    if(max(xp) > (csi+alpha/kappa)){
      ln.GL <- 1e6
    }else{
      ln.GL <- (-(- n*log(alpha) + sum((1/kappa - 1)*log(1 - kappa/alpha*(xp - csi)) - (1 - kappa/alpha*(xp - csi))^(1/kappa)) +
                    (p - 1)*log(a + kappa) + (q - 1)*log(a - kappa)))
    }
  }
  
  return(ln.GL) # naturalmente, o sinal do ln é negativo, então no momento q inverte o sinal pra minimizar a função,
  # precisa tirar reverter o sinal do máximo depois
  
}

fun.gl.beta.gev <- function(param0,        # vetor c/ os parâmetros, 1. csi 2. alpha, 3. kappa
                            media.kappa,   # -0,10
                            desvpad.kappa, # 0,122
                            xp){
  
  # Paramentros Iniciais
  csi <- param0[1]
  alpha <- param0[2]
  kappa <- param0[3]
  
  n <- length(xp)
  var.kappa <- desvpad.kappa^2
  
  # Distribuicao Beta
  a <- 0.5
  b <- media.kappa+a
  p <- b^2*(1-b)/var.kappa-b
  q <- p*(1/b-1)
  
  zi <- 1 - kappa/alpha*(xp - csi)
  
  # Fallbacks
  if(alpha <= 0) return(1e6)
  if(a - abs(kappa) <= 1e-8) return(1e6)
  if(any(zi <= 1e-8)) return(1e6 + abs(min(zi))*1e4)
  
  # Generalized log-likelihood
  term1 <- (1/kappa - 1)*log(zi)
  term2 <- zi^(1/kappa)
  lprior <- (p - 1)*log(a + kappa) + (q - 1)*log(a - kappa)
  ngll <- -(-n*log(alpha) + sum(term1 - term2) + lprior)
  
  return(ngll)
  
}



# PRIORIS -----------------------------------------------------------------

# Distribuições a priori p/ serem plotadas c/ ggplot2po
# Priori Beta
fun.priori.beta <- function(x, a = 0.5, mu, sd){
  
  b <- mu + a
  p <- b^2*(1 - b)/(sd^2) - b
  q <- p*(1/b - 1)
  
  fx <- (a + x)^(p - 1)*(a - x)^(q - 1)/(beta(p,q))
  
  return(fx)
  
} # sinais normais


# MAXIMIZAÇÃO ------------------------------------------------------------------

# A maximização das funções de verossimilhança  (GEV e Gumbel) é feita utilizando o pacote 'optimx' que, por padrão,
# minimiza a função, por isso negativamos as funções fun.l.gev e fun.l.gu.
# Lembrando que o sinal da função de verossimilhança é naturalmente negativo, então queremos o menor valor negativo (maior valor)
# quando buscamos a máxima verossimihança.

# O optimx retorna um data.frame c/ colunas 'p1' a 'pk' parâmetros e o valor mínimo (máximo) 'value' encontrado por cada método de otimização (linhas)

# Máxima verossimilhança (ML) -> tentativa de juntar todas as maximizações em uma única função
# Ainda não está pronta -> adicionar argumentos adicionais
fun.max.l <- function(param = c(NULL, NULL, NULL), # vetor vazio c/ 3 valores na ordem 1. csi 2. alpha 3. kappa
                      xp,                          # vetor c/ série AM observada
                      dist = "gev",                # "gev" ou "gumbel"
                      generalizada = FALSE        # aplicar a verossimilhança generalizada ou não à GEV
){
  
  # Pacotes
  p_load(pacman, optimx, dplyr)
  
  # Check GEV bounds
  check.bound <- function(par, x){
    
    csi <- par[1]; alpha <- par[2]; kappa <- par[3]
    min.z  <- min(1 - kappa/alpha*(x - csi))
    if(min.z <= 0){
      
      delta <- (tol - min.z)*alpha/(abs(kappa))
      if(kappa < 0) csi <- csi - delta # fréchet
      if(kappa > 0) csi <- csi + delta # weibull
      
    }
   
    param[1] <- csi
    return(param)
     
  }
  
  # Testar qual distribuição maximizar
  if(dist == "gev"){
    
    check.par <- check.bound(param, xp)
    
    if(generalizada == FALSE){
      
      max.l <-
        optim(par = check.par,
              fn = fun.l.gev,
              hessian = TRUE,
              xp = xp) %>% {
                # Organizar resultados finais
                tibble(csi = .$par[1],
                       alpha = .$par[2],
                       kappa = .$par[3],
                       value = -(.$value),
                       hess = list(as_tibble(.$hessian)))}
      
    }
    
    if(generalizada == TRUE){
      
      # Verificar se |kappa| > 0.5, caso sim, substituir pela média da distribuição a priori de Martins e Stedinger
      # if(abs(param[3]) > 0.5) param[3] <- media.kappa
      if(abs(check.par[3]) > 0.5) check.par[3] <- media.kappa
      
      # Otimização
      max.gl <- optim(par = check.par,
                      fn = fun.gl.beta.gev,
                      hessian = TRUE,
                      media.kappa = media.kappa,
                      desvpad.kappa = desvpad.kappa,
                      xp = xp) %>% {
                        # Extração dos resultados
                        tibble(csi = .$par[1],
                               alpha = .$par[2],
                               kappa = .$par[3],
                               value = -(.$value),
                               hess = list(as_tibble(.$hessian)))}
      
    }
    
  } else{
    
    if(dist == "gumbel"){
      
      max.l <-
        optim(par = param,
              fn = fun.l.gu,
              hessian = TRUE,
              xp = xp) %>% {
                # Organizar resultados finais
                tibble(csi = .$par[1],
                       alpha = .$par[2],
                       value = -(.$value),
                       hess = list(as_tibble(.$hessian)))}
      
    } else{errorCondition(message = " Escolha 'gumbel' ou 'gev'.")}
    
  }
  
  max.l <- unique(max.l)
  
  return(max.l)
  
}

# Máxima verossimilhança (ML)
fun.max.l.old <- function(param = c(NULL, NULL, NULL), # vetor vazio c/ 3 valores na ordem 1. csi 2. alpha 3. kappa
                      xp,                          # vetor c/ série AM observada
                      dist = "gev"                 # "gev" ou "gumbel"
){
  
  # Pacotes
  p_load(pacman, dplyr)
  
  # Testar qual distribuição maximizar
  if(dist == "gev"){
    
    # Verificar se |kappa| > 0.5, caso sim, substituir pela média da distribuição a priori de Martins e Stedinger
    if(param[3] == 0) param[3] <- -0.1
    
    max.l <-
      optim(par = param,
            fn = fun.l.gev,
            hessian = TRUE,
            xp = xp) %>% {
              # Organizar resultados finais
              tibble(csi = .$par[1],
                     alpha = .$par[2],
                     kappa = .$par[3],
                     value = -(.$value),
                     hess = list(as_tibble(.$hessian)))}
    
  } else{
    
    if(dist == "gumbel"){
      
      max.l <-
        optim(par = param,
              fn = fun.l.gu,
              hessian = TRUE,
              xp = xp) %>% {
                # Organizar resultados finais
                tibble(csi = .$par[1],
                       alpha = .$par[2],
                       value = -(.$value),
                       hess = list(as_tibble(.$hessian)))}
      
    } else{errorCondition(message = " Escolha 'gumbel' ou 'gev'.")}
    
  }
  
  max.l <- unique(max.l)      # força a retornar somente uma linha, caso não haja otimização
  # max.l$value <- -max.l$value # corrige sinal do máximo
  
  return(max.l)
  
}

# Ver depois se não da pra fazer uma função de maximizaçao única
# Máxima verossimilhança generalizada
fun.max.gl <- function(param = c(NULL, NULL, NULL),
                       xp,
                       media.kappa,
                       desvpad.kappa){
  
  # Pacotes
  p_load(pacman, optimx, dplyr)
  
  # Verificar se |kappa| > 0.5, caso sim, substituir pela média da distribuição a priori de Martins e Stedinger
  if(abs(param[3]) > 0.5 | param[3] == 0) param[3] <- media.kappa
  
  # Check GEV bounds
  check.bound <- function(par, x){
    
    csi <- par[1]; alpha <- par[2]; kappa <- par[3]
    min.z  <- min(1 - kappa/alpha*(x - csi))
    if(min.z <= 0){
      
      delta <- (tol - min.z)*alpha/(abs(kappa))
      if(kappa < 0) csi <- csi - delta # fréchet
      if(kappa > 0) csi <- csi + delta # weibull
      
    }
    
    param[1] <- csi
    return(param)
    
  }
  
  check.par <- check.bound(param, xp)
  
  # Otimização
  max.gl <- optim(par = check.par,
                  fn = fun.gl.beta.gev,
                  hessian = TRUE,
                  media.kappa = media.kappa,
                  desvpad.kappa = desvpad.kappa,
                  xp = xp) %>% {
                    # Extração dos resultados
                    tibble(csi = .$par[1],
                           alpha = .$par[2],
                           kappa = .$par[3],
                           value = -(.$value),
                           hess = list(as_tibble(.$hessian)))}
  
  return(max.gl)
  
} # retorna um dataframe com 5 colunas: csi, alpha, kappa, value e hess

# L-RATIO -----------------------------------------------------------------

fun.tao3 <- function(kappa){
  
  tao3 <- 2*(1 - 3^(-kappa))/(1 - 2^(-kappa)) - 3
  
  return(tao3)
  
}

fun.tao4 <- function(kappa){
  
  tao4 <- (1 - 6*2^(-kappa) + 10*3^(-kappa) - 5*4^(-kappa))/(1 - 2^(-kappa))
  
  return(tao4)
  
}

# DERIVADAS PARCIAIS ------------------------------------------------------

# Derivada parcial de yp em relação a csi: ∂yp/∂ξ = 1

# Derivada parcial de yp em relação a alpha
fun.deriv.alpha <- function(p, kappa){
  
  deriv.alpha <- 1/kappa*(1 - (-log(p))^kappa)
  return(deriv.alpha)
  
}

# Derivada parcial do quantil em relação a kappa
fun.deriv.kappa <- function(p, alpha, kappa){
  
  yp <- -log(p)
  deriv.kappa <- - alpha/(kappa^2)*(1 - yp^kappa) - alpha/kappa*yp^kappa*log(yp)
  return(deriv.kappa)
  
}

# RESUMO DAS ESTATÍSTICAS BÁSICAS -----------------------------------------

# Função que calcula a tabela resumo igual Papalexiou e Koutsoyiannis (2013)
summary.stats <- function(x){
  gg.stats <- c(
    n_years = length(x),
    mean = mean(x, na.rm = TRUE),
    p0 <- sum(x == 0)/length(x),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    t3 = lmom::samlmu(x)[[3]],
    t4 = lmom::samlmu(x)[[4]]
  )
}

fun.estat.bas <- function(x){
  c(min = min(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    mediana = median(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    media = mean(x, na.rm = TRUE),
    assimetria = skewness(x, na.rm = TRUE))}

fun.estat.bas2 <- function(x){
  c(min = min(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    mediana = median(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE), # adiciona o desvio
    media = mean(x, na.rm = TRUE),
    assimetria = skewness(x, na.rm = TRUE))}

