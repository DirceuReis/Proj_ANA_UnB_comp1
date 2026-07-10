
# FUNÇÕES INTERNAS -------------------------------------------------------

#' @title Estimativas iniciais com L-momentos
start.gev <- function(ams){
  sample.lmom <- lmom::samlmu(ams)
  par <- lmom::pelgev(sample.lmom)
  return(setNames(par, c("xi", "alpha", "kappa")))
}

#' @title Função de log-verossimilhança da GEV
nll.gev <- function(param, y, prior.info) {
  xi <- param[1]
  alpha <- param[2]
  kappa <- param[3]
  n <- length(y)

  if(!is.null(prior.info)) {
    mu.kappa <- prior.info[1]
    var.kappa <- prior.info[2]^2

    a <- 0.5
    b <- mu.kappa + a
    p <- b^2 * (1 - b) / var.kappa - b
    q <- p * (1 / b - 1)

    if(a - abs(kappa) <= 1e-6) return(1e6 + abs(kappa) * 1e4) # impedir log(<=0)
    l.prior <- (p - 1) * log(a + kappa) + (q - 1) * log(a - kappa)
  } else{
    l.prior <- 0
  }

  zi <- 1 - kappa / alpha * (y - xi)

  # Fallbacks
  # if(abs(kappa) < 1e-5) return(1e6 + 1/abs(kappa)) # wall-up infinitesimal "shape"
  # if(xi <= 1e-8) return(abs(1e6 + abs(xi) * 1e4))
  if(alpha <= 1e-8) return(1e6 + abs(alpha) * 1e4)
  if(any(zi <= 1e-8)) return(1e6 + abs(min(zi)) * 1e4)

  # Negative log-likelihood
  term1 <- (1 / kappa - 1) * log(zi)
  term2 <- zi^(1 / kappa)
  nll <- -(-n * log(alpha) + sum(term1 - term2) + l.prior)

  return(unname(nll))
}

#' @title Função gradiente da log-verossimilhança
grad.ll.gev <- function(param, y, prior.info) {
  xi <- param[1]
  alpha <- param[2]
  kappa <- param[3]
  N <- length(y)

  if (alpha <= 1e-8) {
    return(c(0, 1e6, 0)) # corrigir alpha (param[3]) e não os demais
  }

  # Avaliar se há distribuição a priori
  if (!is.null(prior.info)) {
    mu.kappa <- prior.info[1]
    var.kappa <- prior.info[2]^2

    a <- 0.5
    b <- mu.kappa + a
    p <- b^2 * (1 - b) / var.kappa - b
    q <- p * (1 / b - 1)

    if (a - abs(kappa) <= 1e-6) {
      return(rep(1e6, 3))
    }

    dl.prior <- (p - 1) / (a + kappa) - (q - 1) / (a - kappa)
  } else {
    dl.prior <- 0
  }

  # Termos simplificadores
  zi <- 1 - kappa / alpha * (y - xi)
  if (any(zi <= 1e-8)) {
    return(rep(1e6, 3))
  }
  inv.z <- 1 / zi
  z.inv.k <- zi^(1 / kappa)
  z.inv.km1 <- zi^(1 / kappa - 1)
  ln.z <- log(zi)

  # Derivadas em relação a zi
  dz.dxi <- kappa / alpha
  dz.dalpha <- kappa / alpha^2 * (y - xi)
  dz.dkappa <- -(y - xi) / alpha

  # Gradientes
  grad.xi <- -sum(
    (1 / kappa - 1) * inv.z * dz.dxi - 1 / kappa * z.inv.km1 * dz.dxi
  )
  grad.alpha <- N / alpha - sum((1 / kappa - 1) * inv.z * dz.dalpha - 1 / kappa * z.inv.km1 * dz.dalpha)
  grad.kappa <- -sum((1 / kappa - 1) * inv.z * dz.dkappa - ln.z / kappa^2 - z.inv.k * (1 / kappa * inv.z * dz.dkappa - ln.z / kappa^2)) -
    dl.prior

  return(c(grad.xi, grad.alpha, grad.kappa))
} # fim 'grad.ll.gev()'


# FUNÇÕES PRINCIPAIS -----------------------------------------------------

#' @title Função para ajuste da GEV com máxima verossimilhança
fgev <- function(y, prior.info = NULL, maxit = 1e8) {

  # Here, 'y' is a column of the y.id matrix
  # Get initial parameters with L-moments
  start <- start.gev(y)

  # Definir limites para L-BFGS-B
  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), -Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)

  # Usando Nelder-Mead
  fit <- tryCatch(
    {
      fit.obj <- optim(
        par = start, fn = nll.gev, method = "L-BFGS-B", gr = grad.ll.gev, lower = lower, upper = upper, hessian = TRUE,
        y = y, prior.info = prior.info, control = list(
          maxit = maxit,
          parscale = abs(start), # escala dos parâmetros
          factr = 1e10, # relaxa a tolerância da função
          pgtol = 1e-8, # adiciona tolerância no gradiente
          lmm = 20     # aumenta a memória da hessiana aproximada
          # trace = 1     # acompanhar o progresso do ajuste
        )
      )
    },
    error = function(e) {
      warning("Erro ao minimizar no processamento por duração: ", e)
      return(NULL)
    }
  )

  if(is.null(fit) || fit$convergence != 0){
    fit <- tryCatch(
      optim(
        par = start, fn = nll.gev, method = "Nelder-Mead", hessian = TRUE, 
        y = y, prior.info = prior.info, control = list(
          maxit = maxit,
          parscale = abs(start)
          # reltol = 1e-10
        )
      ),
      error = function(e){
        warning("Erro no processamento por duração: ", e)
        warning(sprintf("Chute inicial com L-momentos: %0.3f", start[[3]]))
        return(NULL)
      }
    )
  }

  return(fit)
}

#' @title Função para estimar quantis da GEV
qgev <- function(param, p.nexceed){
  p <- p.nexceed ; xi <- param[1]; alpha <- param[2]; kappa <- param[3]
  yp <- xi + alpha/kappa*(1 - (-log(p))^kappa)
  return(yp)
}
