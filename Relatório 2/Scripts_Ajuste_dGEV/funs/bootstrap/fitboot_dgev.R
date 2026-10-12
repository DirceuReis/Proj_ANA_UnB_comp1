
# FUNCOES PRINCIPAIS DE AJUSTE -------------------------------------------

## FUNÇÕES INTERNAS ------------------------------------------------------

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
grad.nll.dgev <- function(param, durations, i.jd, prior.info){
  
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
  
} # fim 'grad.nll.dgev()'


## FUNÇÕES PRINCIPAIS -----------------------------------------------------

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
      optim(par = start, fn = nll.dgev, gr = grad.nll.dgev, method = "L-BFGS-B", lower = lower, upper = upper, hessian = TRUE,
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



# FUNCOES INCERTEZA BOOTSTRAP --------------------------------------------

#' @title Ajustar d-GEV e estimar quantis com intervalos de confiança bootstrap
#' @details
#' Função para ajustar o modelo d-GEV para séries de intensidades máximas anuais de diferentes durações,
#' estimar quantis de chuva intensa e intervalos de confiança baseados em reamostragens bootstrap e percentis.
#' @param data lista com data.frames de intensidades máximas anuais de diferentes estações
#' @param prior.info vetor com média e desvio padrão regionais para o parâmetro de forma da GEV (ou NULL)
#' @param durations vetor com durações [h] para estimar os quantis de chuvas intensas
#' @param return.perid vetor com períodos de retorno [anos] para estimar quantis de chuvas intensas
#' @param cols nome das colunas dos data.frames contidos em `data`
#' @param method otimizador passado em `optim()` (ou 'Nelder-Mead' ou 'L-BFGS-B')
#' @param scale.inv método para gerar estimativa inicial do expoente da relação de invariância de escala (ou 'strict' ou 'wide)
#' @param B número de réplicas bootstrap (recomenda-se B > 2000)
#' @param conf nível de confiaça para intervalo de confiança com percentis
#' @param maxit número máximo de iterações, passado ao `optim()`
fitboot.dgev <- function(
  data, prior.info = NULL,
  durations, return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  method = c("Nelder-Mead", "L-BFGS-B"),
  scale.inv = c("wide", "strict"),
  B = 3000,
  conf = 0.95,
  maxit = 1e6
){

  #### CHECKINGS ----

  if(!inherits(data, "list")) stop("'data' deve ser uma lista com data.frames.")
  if(!length(prior.info) %in% c(0, 2)) stop("'prior.info' dever ser 'NULL' ou um vetor de comprimento 2.")

  #### INTERNAL FUNCTIONS ----

  # Estimate confidence intervals with Bootstrap percentile method
  boot.dgev <- function(data.gg){

    # Resampling framework
    i.jd <- data.gg |>                 # create a matrix year x d, filled with imax
      pivot_wider(names_from = cols[3], # transformar durações em colunas
                  values_from = cols[4],
                  names_prefix = "d_") |>
      arrange(cols[2]) |>
      select(-cols[1:2]) |>
      na.omit() |> as.matrix()

    if(is.null(durations)) durations <- as.double(gsub("d_", "", colnames(i.jd)))

    N <- nrow(i.jd)

    # Generate a matrix of B x n.years
    idx.samples <- matrix(sample(seq_len(N), size = N*B, replace = TRUE), nrow = B, ncol = N) # resampled indices
    idx.original <- seq_len(N)                   # original indices
    indices <- rbind(unname(idx.original), idx.samples)  # bind original indices

    # Fit and estimate quantiles
    boot.dgev <- lapply(seq_len(B + 1), function(i){
      idx <- indices[i,] # extract row with indices, where i = 1 is the original data
      fit <- tryCatch( 
        fdgev(i.jd[idx, ], prior.info, cols, method, scale.inv, maxit),
        error = function(e) NULL
      )
      if(is.null(fit)) return(NULL)
      idf <- q.dgev(fit, durations, return.period) # matrix with quantiles
      q.names <- paste0(rep(return.period, length(durations)), "_", rep(durations, each = length(return.period)))
      return(c(fit$par, setNames(as.vector(idf), q.names)))
    }) # end of 'bood.dgev'

    # Organize results
    original <- boot.dgev[[1]]    # results from original data
    boot.samples <- boot.dgev[-1] # resulta from sampled data
    boot.samples <- boot.samples[!sapply(boot.samples, is.null)] # remove NULL results
    B.success <- length(boot.samples)           # number of 'good' replicates
    boot.matrix <- do.call(rbind, boot.samples) # build matrix
    npar <- 5                                   # number of d-GEV parameters
    boot.par <- boot.matrix[,1:npar, drop = FALSE] # parameter matrix
    boot.q <- boot.matrix[,(npar + 1):(npar + length(durations)*length(return.period))] # quantile  matrix

    # Process original estimates
    gg <- unique(data.gg[[cols[1]]])
    par <- original[1:npar]
    idf <- original[(npar + 1):(npar + length(durations)*length(return.period))]

    # Determine percentile confidence invervals
    ci.par <- apply(boot.par, 2, quantile, probs = c((1-conf)/2, 1 - (1-conf)/2), na.rm = TRUE) # parameter CI
    ci.q <- apply(boot.q, 2, quantile, probs = c((1-conf)/2, 1 - (1-conf)/2), na.rm = TRUE)     # quantile CI

    df.ci.par <- t(ci.par) |> # organize as tbl_df
      as_tibble(.name_repair = "minimal") |> 
      setNames(c("ci_lower", "ci_upper")) |> 
      mutate(parname = colnames(ci.par))

    df.ci.q <- t(ci.q) |> # organize as tbl_df
      as_tibble(.name_repair = "minimal") |> 
      setNames(c("ci_lower", "ci_upper")) |> 
      mutate(id = colnames(ci.q)) |> 
      separate(id, into = c("return_period", "d"), convert = TRUE)

    # Setup outputs
    parameters <- bind_cols(
      gauge_code = gg,
      par = unname(par),
      df.ci.par,
      ci_type = "bootstrap",
      prior = ifelse(is.null(prior.info), FALSE, TRUE),
      which_prior = deparse(prior.info)
    ) |> 
      select(gauge_code, parname, par, ci_lower, ci_upper, ci_type, prior, which_prior)

    quantiles <- bind_cols(
      gauge_code = gg,
      q = idf,
      df.ci.q,
      ci_type = "bootstrap",
      prior = ifelse(is.null(prior.info), FALSE, TRUE),
      which_prior = deparse(prior.info)
    ) |> 
      select(gauge_code, d, return_period, q, ci_lower, ci_upper, ci_type, prior, which_prior)

    out <- list(
      par = parameters,
      idf = quantiles,
      replicates = boot.matrix,
      success = sprintf("%g/%g", B.success, B)
    )

    return(out)

  }

  #### RUN ----

  # Messages
  message("Fitting and estimating quantiles. Confidence intervals method: Boostrap.")
  message(">> Number of gauges: ", length(data))
  message(">> Optimizer method: ", method)
  message(">> Durations [h]: ", paste(durations, collapse = ", "))
  message(">> Return period [years]: ", paste(return.period, collapse = ", "))
  message(">> Bootstrap replicates: ", B)
  message(">> Confidence level: ", paste0(conf*100, "%"))

  ls.res <- pbapply::pblapply(data, boot.dgev)

  return(ls.res)

}