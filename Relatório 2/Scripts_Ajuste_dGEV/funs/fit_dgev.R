
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


# INCERTEZA COM VEROSSIMILHANCA PERFILADA --------------------------------


## FUNÇÕES INTERNAS ------------------------------------------------------

# d-GEV negative log-likelihood function for all but one parameter
prof.nll.dgev <- function(par.free, which.par, par.value, i.jd, ds, prior.info){

  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[which.par] <- par.value          # fixed value for current parameter
  full[-which.par] <- par.free          # remaining parameters to be optimized

  nll.dgev(full, ds, i.jd, prior.info)

}

# Gradient function for d-GEV negative log-likelihood for all but one parameter
prof.grad.dgev <- function(par.free, which.par, par.value, i.jd, durations, prior.info){

  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[which.par] <- par.value          # fixed value for current parameter
  full[-which.par] <- par.free          # remaining parameters to be optimized
  gr <- grad.ll.dgev(full, durations, i.jd, prior.info)
  
  gr[-which.par]

}

# Find the 'roots' of the profile log-likelihood function
roots.proflik <- function(i.jd, ds, step, conf, par, i, ll, lower, upper, prior.info){

  val.i <- par[[i]]
  last.par <- par
  npar <- length(par)
  lim <- ll - qchisq(conf, 1)/2
  parnames <- c("mut", "sigma0", "xi", "H", "theta")

  result <- list()

  while(ll > lim && !is.na(ll)){

    last.ll <- ll
    if(i == 2){
      if(val.i <= lower[i]){
        ll <- NA_real_
        next
      }
    }
    if(i == 5 && step < 0){
      if(val.i <= lower[i]){
        ll <- NA_real_
        next
      }
    }

    fit.prof <- tryCatch(
      optim(
        par = last.par[-i],
        fn = function(internal.par){
          prof.nll.dgev(
            par.free = internal.par, which.par = i, par.value = val.i,
            i.jd = i.jd, ds = ds, prior.info = prior.info
          )
        },
        gr = NULL, method = "L-BFGS-B", lower = lower[-i], upper = upper[-i],
        control = list(maxit = 1e8, factr = 1e9, pgtol = 1e-8)
      ), error = function(e) NULL
    )

    if(!is.null(fit.prof) && fit.prof$convergence == 0 && fit.prof$value < 1e6){
      ll <- -fit.prof$value # save new profile value 'll'
    } else{
      fit.prof <- tryCatch(
        optim(
          par = last.par[-i],
          fn = function(internal.par){
            prof.nll.dgev(
              par.free = internal.par, which.par = i, par.value = val.i,
              i.jd = i.jd, ds = ds, prior.info = prior.info
            )
          },
          method = "Nelder-Mead", control = list(maxit = 1e8)
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10) && is.finite(fit.prof$value)){
        ll <- -fit.prof$value
        if(ll < -1e6) ll <- NA_real_
      } else{
        ll <- NA_real_
      }
    }

    full <- numeric(npar)
    if(!is.null(fit.prof) && !is.na(ll)){
      full[-i] <- fit.prof$par
    } else{
      full[-i] <- last.par[-i]
    }
    full[i] <- val.i

    val.i <- val.i + step               # update 'val.i' for next iteration
    slope <- (ll - last.ll)/step        # profile surface slope (check for inconsistency)
    last.par <- full                    # update 'last.par' with current results
    current.iter <- length(result) + 1  # update counter

    result[[current.iter]] <- c(full, ll, slope) # save current results

    # message(sprintf(
    #   "ci %-3s | %s: %0.3f | ll: %0.3f | slope: %0.3f",
    #   ifelse(step < 0, "low", "up"), parnames[i], val.i, ll, slope
    # ))

  } # end of 'while' search

  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(parnames, "loglik", "slope")

  return(prof.par)

}

# Estimate profile likelihood confidence intervals for d-GEV parameters
ciprof.par <- function(i.jd, ds, par, max.ll, conf, step, lltol, gg, prior.info){

  # message("Estimating parameter CIs...")
  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  parnames <- c("mut", "sigma0", "xi", "H", "theta")
  
  
  ls.profile <- lapply(X = seq_along(par), FUN = function(i){
  
    step.par <- max(abs(step*par[[i]]), step) + 1e-7 # passo redefinido para evitar xi == 0

    # Search left (lower CI limit)
    if(par[i] > lower[i] + step.par){
      left <- roots.proflik(i.jd, ds, -step.par, conf, par, i, max.ll, lower, upper, prior.info)
    } else{
      # message(sprintf("Optimal `%s` is at the edge of its domain`, returning `ci_lower == par[i]`", parnames[i]))
      left <- matrix(c(par, max.ll, 0), nrow = 1)
      left[1, i] <- lower[i]
      colnames(left) <- c(parnames, "loglik", "slope")
    }
    
    # Search right (upper CI limit)
    if(par[i] < upper[i] - step.par){
      right <- roots.proflik(i.jd, ds, step.par, conf, par, i, max.ll, lower, upper, prior.info)
    } else{
      # message(sprintf("Optimal `%s` is at the edge of its domain`, returning `ci_upper == par[i]`", parnames[i]))
      right <- matrix(c(par, max.ll, 0), nrow = 1)
      right[1, i] <- upper[i]
      colnames(right) <- c(parnames, "loglik", "slope")
    }
    
    proflik <- rbind(left, right)
    proflik <- proflik[order(proflik[,i]),]
    ci.lower <- min(proflik[,i], na.rm = TRUE)
    ci.upper <- max(proflik[,i], na.rm = TRUE)
    ll.max.proflik <- max(proflik[,"loglik"], na.rm = TRUE)
    par.max.proflik <- proflik[which.max(proflik[,"loglik"]), seq_len(length(par))]
    if(par[i] != 0){
      if(abs(par.max.proflik[i] - par[i])/par[i] > lltol){
        message <- 1               # new max found
        out.ll <- ll.max.proflik   # update max ll
        out.par <- par.max.proflik # update parameters
      } else{
        message <- 0     # keeping same max
        out.ll <- max.ll # same max ll
        out.par <- par   # same parameters
      }
    } else{
      message <- 0     # keeping same max
      out.ll <- max.ll # same max ll
      out.par <- par   # same parameters
    }
    
    # plot(proflik[,c(parnames[i], "loglik")], type = "p")
    # segments(x0 = par[i], x1 = par[i], y0 = max.ll, y1 = max.ll - qchisq(conf, 1)/2, lty = "dashed", col = "grey50")
    # segments(min(proflik[,i], na.rm = TRUE),  
    #   max(proflik[,i], na.rm = TRUE), 
    #   max.ll - qchisq(conf, 1)/2,
    #   max.ll - qchisq(conf, 1)/2,
    #   col = "red"
    # )
  
    profile <- proflik[,c(parnames[i], "loglik", "slope")]
    colnames(profile) <- c("par", "loglik", "slope")
    out <- list(
      ci = data.frame(
        gauge_code = gg,
        parname = parnames[i],
        par = unname(par[i]),
        ci_lower = ci.lower,
        ci_upper = ci.upper,
        ci_type = "proflik",
        rel_width = (ci.upper - ci.lower)/abs(out.par[i]),
        max_ll = max.ll,                  # original max ll
        lim = max.ll - qchisq(conf, 1)/2, # original chi-squared limit
        found_max = out.ll, # max(max.ll, max(proflik[,"loglik"]))
        prior = deparse(prior.info),
        n_profile = nrow(proflik),
        message = message, # if message == 1, 'prof.pars' should be checked
        conf = conf,
        step = step.par
      ),
      prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(length(par) + 1)]), # + 1 pega o loglik
      # prof.par = out.par,
      profile = cbind(
        gauge_code = gg, # keep as comment to avoid having to call it inside every function
        parname = parnames[i], as.data.frame(profile)
      )
    )
  
    return(out)
  
  })
  
  results <- bind_rows(lapply(ls.profile, `[[`, "ci"))
  prof.par <- bind_rows(lapply(ls.profile, `[[`, "prof.par"))
  profile <- bind_rows(lapply(ls.profile, `[[`, "profile"))
  
  if(any(results$message == 1)){
    new.par <- unlist(prof.par[which.max(prof.par$loglik), seq_len(length(par))])
    new.max <- max(profile[,"loglik"], na.rm = TRUE)
    results$par <- new.par
    results$max <- new.max
  }
  
  return(list(ci = results, profile = profile))

}

# d-GEV negative log-likelihood function specific to quantile CI estimation
prof.nll.dgev.q <- function(
  par.free, q0, p.nexceed, i.jd, 
  ds, # vetor com durações (nome coluna i.jd)
  duration, # escalar (1 valor somente) para calcular o quantil
  prior.info
){

  p <- p.nexceed; d <- duration
  sigma0 <- par.free[[1]]
  xi <- par.free[[2]]
  H <- par.free[[3]]
  theta <- par.free[[4]]

  A <- (d + theta)^H

  if(abs(xi) < 1e-5){
    mu0 <- q0*A + sigma0*log(-log(p)) # calcular quantil Gumbel p/ evitar erro numérico
  } else{
    mu0 <- q0*A - sigma0/xi*(1 - (-log(p))^xi)
  }

  mut <- mu0/sigma0
  full <- c(mut, par.free)

  return(nll.dgev(full, ds, i.jd, prior.info))

}

# Vectorized quantile function
q.dgev.vec <- function(param, d, p){
  
  mut <- param[1]; sigma0 <- param[2]; xi <- param[3]; H <- param[4]; theta <- param[5]
  A <- sigma0/(d + theta)^H
  B <- mut + (1/xi)*(1 - (-log(p))^xi)
  
  qp <- A*B
  
  return(unname(qp))
  
}

roots.proflik.q <- function(step, q0, i.jd, par, ds, d, p, conf, ll, lower, upper, prior.info){

  npar <- length(par)
  last.par.free <- par[-1]
  lim <- ll - qchisq(conf, 1)/2
  parnames <- c("mut", "sigma0", "xi", "H", "theta")

  result <- list()
  
  while(ll > lim && !is.na(ll)){
  
    fit.prof <- tryCatch(
      optim(
        par = last.par.free, fn = function(internal.par){
          prof.nll.dgev.q(
            par.free = internal.par,
            q0 = q0, p.nexceed = p, i.jd = i.jd,
            ds = ds, duration = d,
            prior.info = prior.info
          )
        },
        gr = NULL, method = "L-BFGS-B", lower = lower[-1], upper = upper[-1],
        control = list(maxit = 1e8, factr = 1e9, pgtol = 1e-8)
      ),
      error = function(e) NULL
    )
    
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
      ll <- -fit.prof$value
    } else{
      fit.prof <- tryCatch(
        optim(
          par = last.par.free, fn = function(internal.par){
            prof.nll.dgev.q(
              par.free = internal.par,
              q0 = q0, p.nexceed = p, i.jd = i.jd,
              ds = ds, duration = d,
              prior.info = prior.info
            )
          },
          method = "Nelder-Mead", control = list(maxit = 1e8)
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10) && is.finite(fit.prof$value)){
        ll <- -fit.prof$value
        if(ll < -1e6){
          ll <- NA_real_
        }
      } else{
        ll <- NA_real_
      }
    }
    
    if(!is.null(fit.prof) && !is.na(ll)) last.par.free <- fit.prof$par
    
    # # Check GEV limits
    sigma0 <- last.par.free[[1]]
    xi <- last.par.free[[2]]
    H <- last.par.free[[3]]
    theta <- last.par.free[[4]]

    A <- (d + theta)^H
    if(abs(xi) < 1e-5){
      mu0 <- q0*A + sigma0 * log(-log(p)) 
    } else{
      mu0 <- q0*A - (sigma0 / xi) * (1 - (-log(p))^xi)
    }
    mut <- mu0 / sigma0
    
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(d, q0, mut, last.par.free, ll)
    q0 <- q0 + step
  
    # message(sprintf(
    #   "p: %0.3f | d: %2d | q: %0.3f mm/h | ll: %0.3f | conv: %g",
    #   p, d, q0, ll, fit.prof[["convergence"]]
    # ))
  
  }
  
  rp <- 1/(1 - p)
  if(length(result) == 0) return(NULL)
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c("d", sprintf("q_%g", rp), parnames, "loglik")
  
  return(prof.par)

}

# Estimate profile likelihood confidence intervals for quantiles
ciprof.q <- function(i.jd, ds, durations, return.period, par, max.ll, conf, step, lltol, gg, prior.info){

  # message("Estimating quantile CIs...")

  # Parameter bounds
  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  
  # All combinations of duration and return period
  if(is.null(durations)){
    grid <- expand.grid(d = ds, rp = return.period) # use durations from the data given
  } else{
    grid <- expand.grid(d = durations, rp = return.period) # use user defined durations
  }

  # Loop through duration and return period combinations
  ls.profile <- lapply(seq_len(nrow(grid)), function(i){
  
    out <- tryCatch({
      
      d <- grid[i, "d"]
      rp <- grid[i, "rp"]
      # message(sprintf("d: %g h | return_period: %g years", d, rp))
      p <- 1 - 1/rp # non-excceedence probability
      q0 <- q.dgev.vec(par, d, p) # original quantile estimate
      step.q <- if(abs(q0) > 1e-5) step*abs(q0)*0.1 else step
      # step.q <- max(abs(step*q0), step) + 1e-7
      left <- roots.proflik.q(-step.q, q0, i.jd, par, ds, d, p, conf, max.ll, lower, upper, prior.info) # search right
      right <- roots.proflik.q(step.q, q0, i.jd, par, ds, d, p, conf, max.ll, lower, upper, prior.info) # search left
  
      proflik <- rbind(left, right)
      q.name <- sprintf("q_%g", rp)
      ci.lower <- min(proflik[,q.name], na.rm = TRUE)
      ci.upper <- max(proflik[,q.name], na.rm = TRUE)
      censored <- if(any(is.na(proflik[,"loglik"]))) 1 else 0
      neval <- nrow(proflik)
      
      if(neval <= 20){
        message <- 1 # possible error, few profile values
      } else{
        message <- 0 # might be ok, but check anyway
      }
      
      profile <- proflik[order(proflik[,q.name]),][,c(q.name, "loglik")]
      colnames(profile) <- c("q", "loglik")
      out <- list(
        ci = data.frame(
          gauge_code = gg,
          d = d,
          return_period = rp,
          q = q0,
          ci_lower = ci.lower,
          ci_upper = ci.upper,
          max = max.ll,
          lim = max.ll - qchisq(conf, 1)/2,
          prior = deparse(prior.info),
          n_profile = nrow(proflik),
          flag_shape = if(abs(par[3]) > 0.5) 1 else 0, # flag '1' whenever |xi| is larger than 0.
          censored_ci = censored,
          conf = conf,
          step = step.q,
          error = 0L,
          error_message = NA_character_
        ),
        # prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
        profile = cbind(
          gauge_code = gg, d = d,
          return_period = rp, as.data.frame(profile),
          error = 0L
        )
      )
      
      return(out)

    }, error = function(e){

      list(
        ci = data.frame(
          gauge_code = gg, 
          d = d,
          return_period = rp, 
          q = q0,
          ci_lower = NA_real_, 
          ci_upper = NA_real_,
          max = max.ll, 
          lim = max.ll - qchisq(conf, 1)/2,
          prior = deparse(prior.info),
          n_profile = NA_integer_,
          flag_shape = if(abs(par[3]) > 0.5) 1 else 0, # flag '1' whenever |xi| is larger than 0.
          censored_ci = 1L,
          conf = conf,
          step = step,
          error = 1L,
          error_message = conditionMessage(e)
        ),
        profile = data.frame(
          gauge_code = gg,
          d = d,
          return_period = rp,
          q = NA_real_,
          loglik = NA_real_,
          error = 1L
        )
      )

    })

    return(out)
    
  })
  
  results <- bind_rows(lapply(ls.profile, `[[`, "ci"))
  profile <- bind_rows(lapply(ls.profile, `[[`, "profile"))
  
  return(list(ci = results, profile = profile))

}


## FUNÇÃO VEROSSIMILHANÇA PRINCIPAL --------------------------------------

fitprof.dgev <- function(
  data, prior.info = NULL,
  return.period, durations = NULL,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.005, step = 0.005,
  method = "L-BFGS-B", scale.inv = "wide"
){

  #### CHECAGENS ----

  if(!is.list(data) || length(data) == 0){
    stop("`data` must be a non-empty list with data.frames.")
    warning("`cols` must contain column names matching the following information in the same order: 'gauge_code', 'year', 'duration', 'imax")
  }
  if(!is.numeric(return.period) || length(return.period) == 0) stop("`return.period` must have one or more numeric values.")
  if(!length(prior.info) %in% c(0,2)) stop("`prior.info` must be either 'NULL' or a vector of length == 2 with regional mean and standard deviation on shape")
  if(step <= 0) stop("A positive 'step' must be defined.")
  if(!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) stop("Confidence level `conf` must be between (0,1).")
  if(!is.numeric(lltol) || length(lltol) != 1 || lltol <= 0) stop("Profile likelihood change tolerance `lltol` musc be a positive number.")
  data.cols <- colnames(data[[1]])
  missing.cols <- setdiff(cols, data.cols)
  if(length(missing.cols) > 0) stop(sprintf("Missing columns in `data`: %s", paste(missing.cols, collapse = ", ")))
  
  #### CÁLCULOS ----

  message("Fitting d-GEV distribution, estimating quantiles and confidence intervals with 'profile likelihood'.")
  message(">> No. gauges: ", length(data))
  message(">> Durations [h]: ", paste(if(is.null(durations)) "all available" else durations, collapse = ", "))
  message(">> Return periods [years]: ", paste(return.period, collapse = ", "))
  message(">> Prior information: ", deparse(prior.info))
  message(">> Search step: ", step)
  message(">> Confidence level for profile likelihood CIs: ", paste0(conf*100, "%"))

  ls.res <- pbapply::pblapply(data, function(data.gg){
  # ls.res <- lapply(data, function(data.gg){

    gg <- unique(data.gg[[cols[1]]])

    # message(rep("=", getOption("width")))
    # message(sprintf("Gauge: %s", gg))
    # message(rep("-", getOption("width")))

    res <- tryCatch({

      i.jd <- data.gg |> 
        pivot_wider(
          names_from = cols[3],  # turn durations into columns
          values_from = cols[4], # fill with 'imax' values
          names_prefix = "d_"
        ) |> 
        arrange(cols[2]) |>   # arrange by 'year' (or 'wateryear')
        select(-cols[1:2]) |> # remove 'gauge_code' and 'year' columns
        na.omit() |> 
        as.matrix()
      
      ds <- as.integer(gsub("d_", "", colnames(i.jd))) # duration vector
      
      # This case does not require a loop thorugh all durations as they
      # are alll fitted in the same model
      fit <- fdgev(
        i.jd = i.jd,
        prior.info = prior.info,
        cols = cols,
        method = method,
        scale.inv = scale.inv
      )
      
      par <- fit$par
      npar <- length(par)
      
      convergence <- fit$convergence
      max.ll <- -fit$value
  
      # Confidence intervals
      ls.ci.par <- ciprof.par(i.jd, ds, par, max.ll, conf, step, lltol, gg, prior.info)
      ls.ci.q <- ciprof.q(i.jd, ds, durations, return.period, par, max.ll, conf, step, lltol, gg, prior.info)
  
      return(list(
        gauge_code = gg, 
        par = ls.ci.par, 
        quantiles = ls.ci.q, 
        error = 0L,
        error_message = NA_character_
      ))

    }, error = function(e){

      message(sprintf("Error in %s: %s", gg, conditionMessage(e)))
      return(list(
        gauge_code = gg,
        par = NULL,
        quantiles = NULL,
        error = 1L,
        error_message = conditionMessage(e)
      ))

    })

    return(res)

  }) # fim 'ls.res'

  ls.par <- lapply(ls.res, `[[`, "par")
  ls.quantiles <- lapply(ls.res, `[[`, "quantiles")
  ci.par <- bind_rows(lapply(ls.par, `[[`, "ci")); rownames(ci.par) <- NULL
  ci.q <- bind_rows(lapply(ls.quantiles, `[[`, "ci"))
  prof.par <- bind_rows(lapply(ls.par, `[[`, "profile"))
  prof.q <- bind_rows(lapply(ls.quantiles, `[[`, "profile"))

  return(list(
    par = ci.par,
    idf = ci.q,
    prof.par = prof.par,
    prof.q = prof.q
  ))

}