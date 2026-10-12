

# FUNCOES PRINCIPAIS DO AJUSTE -------------------------------------------


## FUNCOES INTERNAS ------------------------------------------------------

#' @title Função de estimativas iniciais do modelo d-GEV
start.dgev <- function(i.jd, durations, scale.inv){
  
  # Estimar parâmetros da GEV p/ cada duração
  par.d <- apply(i.jd, 2, function(imax){
    lmom::pelgev(lmom::samlmu(imax))
  })
  
  # Expoente 'H' (scale) de invariância de escala
  H <- switch(
    scale.inv,
    wide = {
      mom <- apply(i.jd, 2, function(imax) mean(imax))
      H <- -lm(log(mom) ~ log(durations))$coefficient[[2]]
    },
    strict = {
      H.mu <- -lm(log(par.d["sigma",]) ~ log(durations))$coefficient[[2]]
      H.sigma <- -lm(log(par.d["sigma",]) ~ log(durations))$coefficient[[2]]
      H <- mean(H.mu, H.sigma)
    }
  )
  
  if(H < 1e-6 || H > 1 - 1e-6) H <- 0.5
  theta <- 0
  sigma <- par.d[2,1]                          # estimate for first duration
  mut <- par.d[1,1]/sigma                      # rescaled mu
  sigma0 <- sigma*(durations[1] + theta)^H     # rescaled sigma
  xi <- mean(par.d[3,])                     # average xi
  # xi <- ifelse(abs(xi) >= 0.5, prior.info[1], xi)
  xi <- min(0.5 - 1e-6, max(-0.5 + 1e-6, xi))
  
  names <- c("mut", "sigma0", "xi", "H", "theta")
  
  return(setNames(c(mut, sigma0, xi, H, theta), names))
  
} # fim 'start.dgev()'

#' @title Função de log-verossimilhança da d-GEV
#' ACHO QUE NÃO PRECISARIA DE UM ARGUMENTO `durations` AQUI PQ
#' A DURAÇÃO JÁ É O NOME DAS COLUNAS DE `i.jd`
nll.dgev <- function(param, durations, i.jd, prior.info){
  
  # Parâmetros iniciais
  mut <- param[1]    # posição: mu~
  sigma0 <- param[2] # escala: sigma0
  xi <- param[3]  # forma
  H <- param[4]      # invariância de escala (expoente de escala - H ou eta)
  theta <- param[5]  # invariância de escala (parâmetro de posição - theta)

  if(abs(xi) < 1e-5) xi <- sign(xi) * 1e-5
  if(xi == 0) xi <- 1e-5
  
  N <- nrow(i.jd)
  D <- ncol(i.jd)
  
  # Limites teóricos
  if(any(durations + theta <= 0)) return(1e6 + abs(min(durations + theta))*1e4) # impedir offset negativo
  if(theta < 0) return(1e6 + abs(theta)*1e4)    # impedir offset negativo
  if(sigma0 <= 0) return(1e6 + abs(sigma0)*1e4) # impedir escala GEV nulo
  if(H < 1e-6 || H > 1 - 1e-6) return(1e6 + abs(H)*1e4)  # restrição invariância de escala
  
  # Calcular termos z.id = 1 - xi*(i.jd*(d + theta)^H)/sigma - mu~)
  # Função sweep() realiza uma operação nas linhas ou colunas de euma matriz
  A <- (durations + theta)^H      # termo de sigma(d)
  aux <- sweep(x = i.jd/sigma0,   # matriz c/ imax
               MARGIN  = 2,       # operação por colunas
               STATS = A,         # argumentos (nesse caso oq vai ser multiplicado por i.jd)
               FUN = "*")         # multiplicação
  
  z.id <- 1 - xi*(aux - mut) # variável reduzida vetorizada
  
  # Conferir limites teóricos
  # if(any(z.id <= 0)) return(1e8)
  if(any(z.id <= 1e-6)) return(1e6 + abs(min(z.id))*1e4)
  
  # Definir priori log-beta caso tenha sido informada
  if(!is.null(prior.info)){
    
    mu.xi <- prior.info[1]
    var.xi <- prior.info[2]^2
    
    a <- 0.5
    b <- mu.xi+a
    p <- b^2*(1-b)/var.xi-b
    q <- p*(1/b-1)
    
    if(a - abs(xi) <= 1e-6) return(1e6 + abs(xi)*1e4) # impedir log(<0)
    # if(a - abs(xi) <= 1e-6) return(1e6) # impedir log(<0)
    l.prior <- (p - 1)*log(a + xi) + (q - 1)*log(a - xi)
    
  } else{
    
    l.prior <- 0
    
  }
  
  
  # Log-verossimilhança
  nll <- -(-N*D*log(sigma0) + H*N*sum(log(durations + theta)) + (1/xi - 1)*sum(log(z.id)) - sum(z.id^(1/xi)) + l.prior)
  
  return(nll) # retornando negativo p/ "minimização"
  
} # fim 'nll.dgev()'

#' @title Função gradiente log-verossimilhança
grad.ll.dgev <- function(param, durations, i.jd, prior.info){
  
  # Parameters
  mut <- param[1]
  sigma0 <- param[2]
  xi <- param[3]
  H <- param[4]
  theta <- param[5]
  
  ds <- durations
  N <- nrow(i.jd)
  D <- ncol(i.jd)
  
  if(any(ds + theta <= 0) || sigma0 <= 0 || H < 0 || H > 1 - 1e-6){
    return(rep(1e6, 5))
  }
  
  # Variáveis aumuliares
  A <- (ds + theta)^H
  yA <- sweep(x = i.jd, MARGIN = 2, STATS = A, FUN = "*") # produto: i.jd * Ad
  B <- yA/sigma0 - mut                                    # i.jd*A/sigma0 - mut
  z.id <- 1 - xi*B
  
  if(any(z.id <= 0)) return(rep(1e6, 5))
  
  # Pre-calcular termos com z.id
  z.inv <- 1/z.id
  z.pow <- z.id^(1/xi - 1)
  ln.z <- log(z.id)
  z.pow.k <- z.id^(1/xi)
  
  # Derivadas parciais
  # Posição (mut)
  dl.dmut <- (1 - xi)*sum(z.inv) - sum(z.pow)
  
  # Escala (sigma0)
  dl.dsigma0 <- -(N*D)/sigma0 + (1 - xi)*sum(yA*z.inv/(sigma0^2)) - sum(z.pow*yA/(sigma0^2))
  
  # Forma (xi)
  if(!is.null(prior.info)){
    
    mu.xi <- prior.info[1]
    var.xi <- prior.info[2]^2
    
    a <- 0.5
    b <- mu.xi+a
    p <- b^2*(1-b)/var.xi-b
    q <- p*(1/b-1)
    if(a - abs(xi) <= 1e-6) return(rep(1e6, 5))
    
    dl.prior <- (p - 1)/(a + xi) - (q - 1)/(a - xi)
    
  } else{
    
    dl.prior <- 0
    
  }
  
  term1 <- (1/xi - 1)*B*z.inv + ln.z/(xi^2)
  term2 <- (1/xi)*((ln.z*z.pow.k)/xi + B*z.pow)
  dl.dxi <- -sum(term1) + sum(term2) + dl.prior
  
  # Expoente de invariância de escala (H)
  ln.dt <- log(ds + theta)
  aux1 <- sweep(yA/sigma0, MARGIN = 2, ln.dt, FUN = "*")
  dl.dH <- N*sum(ln.dt) + (xi - 1)*sum(z.inv*aux1) + sum(z.pow*aux1)
  
  # Deslocamento (theta)
  dt.inv <- 1/(ds + theta)
  aux2 <- sweep(yA/sigma0, MARGIN = 2, STATS = H*dt.inv, FUN = "*")
  dl.dtheta <- N*H*sum(dt.inv) + (xi - 1)*sum(z.inv*aux2) + sum(z.pow*aux2)
  
  par.names <- c("mut", "sigma0", "xi", "H", "theta")
  grad <- setNames(c(dl.dmut, dl.dsigma0, dl.dxi, dl.dH, dl.dtheta), par.names)
  return(-grad)
  
} # fim 'grad.ll.dgev()'


## FUNCOES PRINCIPAIS ----------------------------------------------------

#' @title Fit a duration dependent GEV with mamumum likelihood and a prior distribution
#' @details This function takes a `data.frame` with annual intensity mamumas for different durations and fits it a d-GEV distribution
#' for a single gauge. For applying it to multiple gauge's time series, lapply() or equivalent can be used.
#' @param data.gg a data.frame specific to a given gauge (gg) containing 1. 'year', 2. 'duration' and 3. 'imax' columns
#' @param yd.gg a matrix where every columns contains 'imax' for a given 'duration', with rows representing different years.
#' @param cols names(data.gg), following the same order specified above
#' @param prior.info a vector with regional mean and standard deviation of GEV's shape parameter to be used as prior knowledge in a
#' Beta distribution inside generalized mamumum likelihood estimation
#' @param method either Nelder-Mead or L-BFGS-B, the latter is tried first, if there is an issue with convergence, fallsback to NM
#' @param scale.inv either 'wide' or 'strict', doesn't have much effect on the end result, only the initial scale invariance exponent (H)
#' @param maxit mamumum number of iterations for 'optim()'
fdgev <- function(i.jd, # antes 'data.gg'
                  prior.info = NULL,
                  cols = c("gauge_code", "year", "duration", "imax"),
                  method = c("Nelder-Mead", "L-BFGS-B"),
                  scale.inv = c("wide", "strict"),
                  maxit = 1e6){
  
  
  #### CHECAGENS ----
  
  # Conferir argumentos
  # if(!is.data.frame(data.gg)) stop("Argumento 'data.gg' deve ser 'tbl_df' e não ", class(data.gg))
  # if(sum(is.element(cols, names(data.gg))) != length(cols)) stop("Colunas indicadas em 'cols' devem estar contidas em 'data.gg'\ncols: ", cols, "\nnames(data.gg): ", names(data.gg))
  if(!is.matrix(i.jd)) stop("'i.jd' must be an n.years x d matrix containing 'imax' values")
  
  method <- match.arg(method, c("Nelder-Mead", "L-BFGS-B"))
  scale.inv <- match.arg(scale.inv, c("wide", "strict"))
  
  #### CÁLCULOS ----
  
  # Configurar dados
  durations <- as.numeric(gsub("d_", "", colnames(i.jd)))
  
  # Definir limites para L-BFGS-B
  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  
  # Parâmetros iniciais (chute)
  start <- start.dgev(i.jd, durations, scale.inv)
  
  # Ajustar d-GEV
  # Usando "L-BFGS-B"
  if(method == "L-BFGS-B"){
    fit <- tryCatch(
      optim(par = start, fn = nll.dgev, gr = grad.ll.dgev, method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE,
            durations = durations, i.jd = i.jd, prior.info = prior.info,
            control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)),
      error = function(e){
        # warning("Erro na estação'", gg, "': ", e$message, ".\nOtimizando com Nelder-Mead.")
        warning("Error: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if(!is.null(fit) && fit$convergence == 0) return(fit)
  }
  
  # Usando Nelder-Mead
  fit <- tryCatch(
    optim(par = start, fn = nll.dgev, gr = NULL, method = "Nelder-Mead", hessian = TRUE,
          durations = durations, i.jd = i.jd, prior.info = prior.info,
          control = list(maxit = maxit)),
    error = function(e){
      # warning("Erro na estação'", gg, "': ", e$message, ".")
      warning("Error: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  # if(fit$value == 1e6) warning("\nErro na estação '", gg, "': não houve otimmização.\nRetornando parâmetros iniciais:\n", paste(round(start, 4), collapse = " "))
  # 
  # return(fit)
  
  
}


#' @title Estimar quantis de intensidade de chuva
#' @details
#' Esta função constrói as curvas de intensidade-duração-frequência a partir dos quantis estimados para diferentes durações e frequências
#' @param fit.obj o resultado da função `fit.dgev()`, contendo uma lista com um resultados padrão da função `optim()` para cada estação
#' @param durations um vetor com durações [h] para as quais serão calculadas as intensidades
#' @param return.period um vetor com períodos de retorno (inverso da probabildiade de excedência) [anos]
#' @returns um `data.frame` com intervalos de confiança estimados para quantil de cada duração, tempo de retorno e estação.
q.dgev <- function(fit.obj, durations, return.period){
  
  #### CHECAGENS ----
  
  if(!inherits(fit.obj, "list")) stop("fit.obj' deve conter uma lista de componentes do 'optim()' com ao menos: par, value e hessian")
  if(length(durations) <= 3) stop("Informe ao menos 3 durações em horas.")
  
  gauges <- names(fit.obj)
  
  #### FUNÇÕES ----
  
  # Calcular quantis
  get.q <- function(param, duration, p.nexceed){
    
    param <- unname(param)
    mut <- param[1]
    sigma0 <- param[2]
    xi <- param[3]
    H <- param[4]
    theta <- param[5]
    
    p <- p.nexceed; d <- duration
    
    # Calcular parâmetros convencionais da GEV
    sigma <- sigma0/(d + theta)^H
    mu <- mut*sigma
    
    qp <- mu + sigma/xi*(1 - (-log(p))^xi)
    
    return(qp)
    
  }
  
  # Calcular quantis p/ todas as estações
  q.all <- function(fit, durations, p.nexceed){
    
    par <- fit.obj$par; p <- p.nexceed; rp <- round(1/(1 - p), 0)
    ls.q <- lapply(durations, get.q, p.nexceed = p, param = par)
    idf.aux <- matrix(unlist(ls.q), length(p), length(durations), dimnames = list(rp , durations))
    
    return(idf.aux)
    
  }
  
  #### CÁLCULOS ----
  
  # Converter return.period em probabilidade de não excedência
  p <- 1 - 1/return.period
  idf <- q.all(fit.obj, durations, p)
  
  return(idf)
  
}


# INCERTEZA COM METODO DELTA ---------------------------------------------

## FUNCOES INTERNAS ------------------------------------------------------

# Fit and estimate quantiles
# source("R/funcoes/fdgev.R")

# Vectorized quantile function
q.dgev.vec <- function(param, d, p){
  
  mut <- param[1]; sigma0 <- param[2]; xi <- param[3]; H <- param[4]; theta <- param[5]
  A <- sigma0/(d + theta)^H
  B <- mut + (1/xi)*(1 - (-log(p))^xi)
  
  yp <- A*B
  
  return(unname(yp))
  
}

# Função gradiente p/ quantil
grad.q.dgev <- function(param, d, p){
  
  mut <- param[1]; sigma0 <- param[2]; xi <- param[3]; H <- param[4]; theta <- param[5]
  logp <- -log(p)
  A <- sigma0/(d + theta)^H
  B <- mut + (1/xi)*(1 - logp^xi)
  yp <- A*B
  
  dy.dmut <- A
  dy.dsigma0 <- B/(d + theta)^H
  dy.dxi <- -A/xi^2*(1 + logp^xi*(xi*log(logp) - 1))
  dy.dH <- -yp*log(d + theta)
  dy.dtheta <- -yp*H/(d + theta)
  
  rp <- 1/(1 - p)
  grad.names <- paste0(names(param), "_d", rep(d, each = length(param)), "_rp", rp)
  return(setNames(c(dy.dmut, dy.dsigma0, dy.dxi, dy.dH, dy.dtheta), grad.names))
  
}


## FUNCAO METODO DELTA PRINCIPAL -----------------------------------------


fitdelta.dgev <- function(data, prior.info = NULL,
                          durations, return.period,
                          cols = c("gauge_code", "year", "duration", "imax"),
                          method = c("L-BFGS-B", "Nelder-Mead"),
                          scale.inv = c("wide", "strict"),
                          conf = 0.95,
                          maxit = 1e6){
  
  #### CHECKINGS ----
  
  if(!inherits(data, "list")) stop("'data' deve ser uma lista com data.frames.")
  if(!length(prior.info) %in% c(0, 2)) stop("'prior.info' dever ser 'NULL' ou um vetor de comprimento 2.")
  method <- match.arg(method, c("L-BFGS-B", "Nelder-Mead"))
  scale.inv <- match.arg(scale.inv, c("wide", "strict"))
  
  #### RUN ----
  
  # Messages
  message("Fitting and estimating quantiles. Confidence intervals method: Delta (Wald)")
  message(">> Number of gauges: ", length(data))
  message(">> Optimizer method: ", method)
  message(">> Durations [h]: ", paste(durations, collapse = ", "))
  message(">> Return period [years]: ", paste(return.period, collapse = ", "))
  message(">> Confidence level: ", paste0(conf*100, "%"))

  probs <- 1 - 1/return.period
  
  ls.res <- pbapply::pblapply(data, function(data.gg){
    
    gg <- unique(data.gg[[cols[1]]])

    i.jd <- data.gg |>                  # create a matrix year x d, filled with imax
      pivot_wider(names_from = cols[3], # transformar durações em colunas
                  values_from = cols[4],
                  names_prefix = "d_") |>
      arrange(cols[2]) |>
      select(-cols[1:2]) |>
      na.omit() |> as.matrix()

    fit <- fdgev(i.jd, prior.info, names(data.gg), method, scale.inv)
    par <- fit$par
    hes <- fit$hessian

    # Durações
    if(is.null(durations)) durations <- sort(unique(data.gg[[cols[3]]]))
    
    # Intervalos de confiança parâmetros
    # Extrair elementos das diagonais
    cov.par <- try(solve(hes), silent = TRUE)
    if(inherits(cov.par, "try-error")) return(matrix(NA, length(par), length(par)))
    var.par <- diag(cov.par)
    sd.par <- sqrt(var.par)
    
    ci.par <- cbind(ci_lower = par - qnorm(1 - (1 - conf)/2)*sd.par, 
                    ci_upper = par + qnorm(1 - (1 - conf)/2)*sd.par)
    
    # Intervalos de confiança quantis
    grid <- expand.grid(d = durations, p = probs)
    ls.ci.q <- lapply(1:nrow(grid), function(i){
      
      d <- grid[i, "d"]
      p <- grid[i, "p"]
      
      q <- q.dgev.vec(par, d, p)
      grad <- grad.q.dgev(par, d, p)
      var.q <- t(grad) %*% cov.par %*% grad
      sd.q <- sqrt(var.q)
      ci.q.l <- q - qnorm(1 - (1 - conf)/2)*sd.q
      ci.q.u <- q + qnorm(1 - (1 - conf)/2)*sd.q
      
      out <- c(d = d, return_period = round(1/(1 - p)), q = q, sd = sd.q, ci_lower = ci.q.l, ci_upper = ci.q.u)
      
    }) # fim 'ls.ci.q'
    
    parameters <- bind_cols(
      gauge_code = gg, parname = names(par), par = par, ci.par, ci_type = "delta",
      prior = ifelse(is.null(prior.info), FALSE, TRUE)
    )
    
    quantiles <- bind_cols(
      gauge_code = gg,
      bind_rows(ls.ci.q),
      ci_type = "delta",
      prior = ifelse(is.null(prior.info), FALSE, TRUE)
    ) |>
      arrange(d, return_period)
    
    out <- list(par = parameters,
                idf = quantiles)
    
  }) # fim 'ls.res'
  
  return(ls.res)
  
}

