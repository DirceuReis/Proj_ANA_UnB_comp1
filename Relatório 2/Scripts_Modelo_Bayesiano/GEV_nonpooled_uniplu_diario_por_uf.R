library(dplyr)
library(tidyr)
library(tibble)
library(rstan)

# =============================================================================
# CAMINHOS
# =============================================================================
DIR_FLUXO <- # pasta para salvar resultados
DIR_DADOS <- "Relatório 2/Scripts_Maximas"

OUT_TAG <- "BR"
TEST_UF <- "MG"       # NULL = Brasil inteiro
TEST_GAUGE <- NULL    # NULL = todas as estações da UF

DIR_OUT <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF <- file.path(DIR_DADOS, "resultados", OUT_TAG, "dataframes")
DIR_GEV <- file.path(DIR_OUT, "gev_nonpooled_by_uf")

dir.create(DIR_GEV, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# MCMC
# =============================================================================
N_CHAINS <- 5L
N_ITER <- 20000L
N_WARMUP <- 10000L
N_THIN <- 50L
ADAPT_DELTA <- 0.999
MAX_TREEDEPTH <- 17L
N_CORES <- min(N_CHAINS, parallel::detectCores())

# =============================================================================
# FILTROS
# =============================================================================
durations_days_allowed <- c(1, 2, 3, 4, 5, 6, 7, 10)
durations_hours_allowed <- 24 * durations_days_allowed
min_years <- 30L
min_n_durations <- length(durations_days_allowed)
D0_val <- 24

achar_rds <- function(nome) {
  p <- file.path(DIR_DF, nome)
  if (file.exists(p)) return(p)
  NULL
}

# =============================================================================
# 1. DADOS
# =============================================================================
F_IMAX <- achar_rds("df_imax_uniplu_hidro.rds")

if(is.null(F_IMAX)){
  stop("Falta df_imax_uniplu_hidro.rds — rode o passo 13")
}

message("Imax Uniplu: ",F_IMAX)

dados_agregados_brutos <- readRDS(F_IMAX)

dados_agregados <- dados_agregados_brutos %>%
  transmute(
    gauge_code = as.character(gauge_code) |> trimws(),
    duration_h = as.numeric(d),
    year = as.integer(wateryear),
    imax = as.numeric(imax),
    estado = as.character(estado)
  ) %>%
  mutate(duration_days = duration_h / 24) %>%
  filter(
    duration_h %in% durations_hours_allowed,
    !is.na(gauge_code), !is.na(duration_h),
    !is.na(year), !is.na(estado)
  )

if (!is.null(TEST_UF) && nzchar(TEST_UF)) {
  dados_agregados <- dados_agregados %>% filter(estado == TEST_UF)
  message("Filtro TEST_UF=", TEST_UF)
}

duracoes_ausentes <- setdiff(
  durations_hours_allowed,
  sort(unique(dados_agregados$duration_h))
)

if (length(duracoes_ausentes) > 0) {
  stop("Durações ausentes (h): ", paste(duracoes_ausentes, collapse = ", "))
}

dados_agregados <- dados_agregados %>%
  filter(is.finite(imax), imax > 0)

duplicatas_ano <- dados_agregados %>%
  count(gauge_code, duration_h, year) %>%
  filter(n > 1)

cat("Duplicatas estação-duração-ano:", nrow(duplicatas_ano), "\n")

max_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

dados_agregados <- dados_agregados %>%
  group_by(gauge_code, duration_h, duration_days, year, estado) %>%
  summarise(imax = max_safe(imax), .groups = "drop") %>%
  filter(is.finite(imax), imax > 0)

years_stats <- dados_agregados %>%
  group_by(gauge_code, duration_h) %>%
  summarise(n_years_valid = n_distinct(year), .groups = "drop")

resumo_filtro_estacoes <- years_stats %>%
  group_by(gauge_code) %>%
  summarise(
    n_duracoes_validas = sum(n_years_valid >= min_years),
    passa_filtro = n_duracoes_validas >= min_n_durations,
    .groups = "drop"
  )

estacoes_aprovadas <- resumo_filtro_estacoes %>%
  filter(passa_filtro) %>%
  pull(gauge_code)

cat("Estações aprovadas:", length(estacoes_aprovadas), "\n")

dados_modelo <- dados_agregados %>%
  filter(gauge_code %in% estacoes_aprovadas) %>%
  group_by(gauge_code, duration_h) %>%
  filter(n_distinct(year) >= min_years) %>%
  ungroup() %>%
  arrange(estado, gauge_code, duration_h, year)

if (!is.null(TEST_GAUGE) && nzchar(TEST_GAUGE)) {
  dados_modelo <- dados_modelo %>% filter(gauge_code == TEST_GAUGE)
  message("Filtro TEST_GAUGE=", TEST_GAUGE)
}

duration_lookup <- dados_modelo %>%
  distinct(duration_days, duration_h) %>%
  arrange(duration_h) %>%
  mutate(duration_id = row_number())

D_vals <- duration_lookup %>% arrange(duration_id) %>% pull(duration_h)
K <- length(D_vals)

if (!isTRUE(all.equal(D_vals, durations_hours_allowed))) {
  stop("Durações finais ≠ ", paste(durations_hours_allowed, collapse = ", "))
}

if (!(D0_val %in% D_vals)) stop("D0=24 h ausente")

ref_idx <- which(D_vals == D0_val)
ufs <- sort(unique(dados_modelo$estado))

message("UFs a processar: ", paste(ufs, collapse = ", "))
message("K=", K, " | obs=", nrow(dados_modelo))
print(duration_lookup)

# =============================================================================
# 2. MODELO STAN
# =============================================================================
modelo_nonpooled_diario_stan <- "
functions {
  real gev_lpdf(real y, real mu, real sigma, real xi) {
    real z = (y - mu) / sigma;

    if (fabs(xi) < 1e-12) {
      return -log(sigma) - z - exp(-z);
    } else {
      real t = 1 + xi * z;
      if (t <= 0) return negative_infinity();

      return -log(sigma)
             - (1 + 1 / xi) * log(t)
             - pow(t, -1 / xi);
    }
  }
}

data {
  int<lower=1> N_obs;
  vector[N_obs] y;
  int<lower=1> K;
  int<lower=1, upper=K> duration_id[N_obs];
  vector<lower=0>[K] D;
  real<lower=0> D0;
}

parameters {
  real alpha;
  real beta;
  real<lower=0, upper=1> H;
  real<lower=0, upper=10> theta;

  vector[K] mu_log;
  vector[K] sigma_log;

  vector<lower=-0.5, upper=0.5>[K] xi_d;
  real<lower=-0.5, upper=0.5> xi_bar;

  real<lower=0> tau2_mu;
  real<lower=0> tau2_sigma;
  real<lower=0> tau2_xi;
}

transformed parameters {
  vector[K] mu;
  vector[K] sigma;

  for (k in 1:K) {
    mu[k] = exp(mu_log[k]);
    sigma[k] = exp(sigma_log[k]);
  }
}

model {
  alpha ~ normal(0, 1e3);
  beta ~ normal(0, 1e3);

  H ~ beta(1, 1);
  theta ~ uniform(0, 10);

  tau2_mu ~ inv_gamma(0.01, 0.01);
  tau2_sigma ~ inv_gamma(0.01, 0.01);
  tau2_xi ~ inv_gamma(0.01, 0.01);

  {
    real a = -0.5;
    real b = 0.5;
    real p = 9;
    real q = 6;
    real lbetaf = lgamma(p) + lgamma(q) - lgamma(p + q);

    target += -lbetaf
      + (p - 1) * log(xi_bar - a)
      + (q - 1) * log(b - xi_bar)
      - (p + q - 1) * log(b - a);
  }

  for (k in 1:K) {
    real duration_term = log((D[k] + theta) / (D0 + theta));

    mu_log[k] ~ normal(
      alpha - H * duration_term,
      sqrt(tau2_mu)
    );

    sigma_log[k] ~ normal(
      beta - H * duration_term,
      sqrt(tau2_sigma)
    );

    xi_d[k] ~ normal(
      xi_bar,
      sqrt(tau2_xi)
    );
  }

  for (i in 1:N_obs) {
    int k = duration_id[i];

    target += gev_lpdf(
      y[i] |
      mu[k],
      sigma[k],
      xi_d[k]
    );
  }
}

generated quantities {
  vector[N_obs] log_lik;

  for (i in 1:N_obs) {
    int k = duration_id[i];

    log_lik[i] = gev_lpdf(
      y[i] |
      mu[k],
      sigma[k],
      xi_d[k]
    );
  }
}
"

message("Compilando Stan ...")
rstan_options(auto_write = TRUE)
options(mc.cores = N_CORES)

sm_nonpooled_diario <- rstan::stan_model(
  model_code = modelo_nonpooled_diario_stan
)

message("Stan OK.")

# =============================================================================
# 3. FUNÇÕES AUXILIARES
# =============================================================================
prepare_uf <- function(uf) {
  dados_uf <- dados_modelo %>% filter(estado == uf)
  
  gauge_lookup <- dados_uf %>%
    distinct(gauge_code) %>%
    arrange(gauge_code) %>%
    mutate(gauge_id = row_number())
  
  dados_stan <- dados_uf %>%
    left_join(gauge_lookup, by = "gauge_code") %>%
    left_join(duration_lookup, by = c("duration_days", "duration_h")) %>%
    arrange(gauge_id, duration_id, year)
  
  global_duration_stats <- dados_stan %>%
    group_by(duration_id) %>%
    summarise(
      median_imax = median(imax, na.rm = TRUE),
      sd_imax = sd(imax, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(tibble(duration_id = seq_len(K)), by = "duration_id") %>%
    arrange(duration_id)
  
  global_median_imax <- pmax(global_duration_stats$median_imax, 1e-3)
  global_sd_imax <- pmax(global_duration_stats$sd_imax, 1e-3)
  
  if (any(!is.finite(global_median_imax)) ||
      any(!is.finite(global_sd_imax))) {
    stop("UF ", uf, ": estatísticas globais incompletas.")
  }
  
  list(
    dados_stan = dados_stan,
    gauge_lookup = gauge_lookup,
    global_duration_stats = global_duration_stats,
    global_median_imax = global_median_imax,
    global_sd_imax = global_sd_imax
  )
}

get_station_data <- function(gauge_code, dados_stan) {
  dados_estacao <- dados_stan %>%
    filter(.data$gauge_code == .env$gauge_code) %>%
    arrange(duration_id, year)
  
  if (nrow(dados_estacao) == 0) {
    stop(sprintf("Estação %s não encontrada.", gauge_code))
  }
  
  dados_estacao
}

make_stan_data_station <- function(gauge_code, dados_stan) {
  dados_estacao <- get_station_data(gauge_code, dados_stan)
  duracoes_presentes <- sort(unique(dados_estacao$duration_id))
  
  if (!identical(duracoes_presentes, seq_len(K))) {
    stop(sprintf(
      "A estação %s não possui todas as %d durações.",
      gauge_code, K
    ))
  }
  
  list(
    N_obs = nrow(dados_estacao),
    y = as.numeric(dados_estacao$imax),
    K = K,
    duration_id = as.integer(dados_estacao$duration_id),
    D = as.numeric(D_vals),
    D0 = as.numeric(D0_val)
  )
}

make_nonpooled_diario_init <- function(gauge_code, ctx) {
  dados_estacao <- get_station_data(gauge_code, ctx$dados_stan)
  
  duration_stats <- dados_estacao %>%
    group_by(duration_id) %>%
    summarise(
      median_imax = median(imax, na.rm = TRUE),
      sd_imax = sd(imax, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(tibble(duration_id = seq_len(K)), by = "duration_id") %>%
    arrange(duration_id)
  
  median_imax <- ctx$global_median_imax
  sd_imax <- ctx$global_sd_imax
  
  idx_valid <- which(
    is.finite(duration_stats$median_imax) &
      is.finite(duration_stats$sd_imax) &
      duration_stats$median_imax > 0 &
      duration_stats$sd_imax > 0
  )
  
  median_imax[idx_valid] <- duration_stats$median_imax[idx_valid]
  sd_imax[idx_valid] <- pmax(
    ctx$global_duration_stats$sd_imax[idx_valid],
    1e-3
  )
  
  H0 <- runif(1, 0.45, 0.70)
  theta0 <- runif(1, 0.8, 1.2)
  
  alpha0 <- rnorm(
    1,
    mean = log(median_imax[ref_idx]),
    sd = 0.03
  )
  
  beta0 <- rnorm(
    1,
    mean = log(sd_imax[ref_idx]),
    sd = 0.03
  )
  
  duration_term <- log(
    (D_vals + theta0) /
      (D0_val + theta0)
  )
  
  xi_bar0 <- runif(1, 0.01, 0.04)
  
  xi_d0 <- pmin(
    pmax(
      rnorm(K, mean = xi_bar0, sd = 0.02),
      -0.49
    ),
    0.49
  )
  
  list(
    alpha = alpha0,
    beta = beta0,
    H = H0,
    theta = theta0,
    mu_log = alpha0 - H0 * duration_term,
    sigma_log = beta0 - H0 * duration_term,
    xi_bar = xi_bar0,
    xi_d = xi_d0,
    tau2_mu = runif(1, 0.02, 0.05),
    tau2_sigma = runif(1, 0.02, 0.05),
    tau2_xi = runif(1, 0.0025, 0.01)
  )
}

fit_nonpooled_diario_station <- function(gauge_code, ctx, seed = 123) {
  stan_data_station <- make_stan_data_station(
    gauge_code,
    ctx$dados_stan
  )
  
  rstan::sampling(
    object = sm_nonpooled_diario,
    data = stan_data_station,
    chains = N_CHAINS,
    cores = N_CORES,
    iter = N_ITER,
    warmup = N_WARMUP,
    thin = N_THIN,
    control = list(
      adapt_delta = ADAPT_DELTA,
      max_treedepth = MAX_TREEDEPTH
    ),
    seed = seed,
    init = function() {
      make_nonpooled_diario_init(gauge_code, ctx)
    },
    refresh = 50
  )
}

# =============================================================================
# 4. LOOP POR UF — RETOMADA AUTOMÁTICA
# =============================================================================
for (uf in ufs) {
  f_out <- file.path(DIR_GEV, paste0(uf, "_fits.rds"))
  f_meta <- file.path(DIR_GEV, paste0(uf, "_meta.rds"))
  
  message("\n========== UF ", uf, " ==========")
  
  ctx <- prepare_uf(uf)
  gauge_list <- sort(unique(ctx$dados_stan$gauge_code))
  
  message("Estações previstas na UF: ", length(gauge_list))
  
  # -------------------------------------------------------------------------
  # Carrega checkpoint existente
  # -------------------------------------------------------------------------
  if (file.exists(f_out)) {
    message("[", uf, "] checkpoint encontrado — carregando: ", f_out)
    
    fits_old <- readRDS(f_out)
    
    fits <- setNames(
      vector("list", length(gauge_list)),
      gauge_list
    )
    
    gauges_salvos <- intersect(
      names(fits_old),
      gauge_list
    )
    
    fits[gauges_salvos] <- fits_old[gauges_salvos]
    rm(fits_old)
    
  } else {
    message("[", uf, "] nenhum checkpoint anterior.")
    
    fits <- setNames(
      vector("list", length(gauge_list)),
      gauge_list
    )
  }
  
  # -------------------------------------------------------------------------
  # Verifica o que já terminou
  # -------------------------------------------------------------------------
  ja_concluidas <- vapply(
    gauge_list,
    function(g) inherits(fits[[g]], "stanfit"),
    logical(1)
  )
  
  gauge_faltando <- gauge_list[!ja_concluidas]
  
  message(
    "[", uf, "] concluídas: ",
    sum(ja_concluidas), " / ", length(gauge_list)
  )
  
  message(
    "[", uf, "] faltando ou com erro: ",
    length(gauge_faltando)
  )
  
  if (length(gauge_faltando) > 0) {
    message(
      "[", uf, "] serão processadas: ",
      paste(gauge_faltando, collapse = ", ")
    )
  }
  
  # -------------------------------------------------------------------------
  # Roda somente o que falta
  # -------------------------------------------------------------------------
  if (length(gauge_faltando) == 0L) {
    message("[", uf, "] todas as estações já possuem stanfit.")
    
  } else {
    
    for (j in seq_along(gauge_faltando)) {
      g <- gauge_faltando[[j]]
      pos_original <- match(g, gauge_list)
      
      message(sprintf(
        "[%s] restante %d/%d | estação %s | posição original %d/%d",
        uf,
        j,
        length(gauge_faltando),
        g,
        pos_original,
        length(gauge_list)
      ))
      
      fits[[g]] <- tryCatch(
        fit_nonpooled_diario_station(
          g,
          ctx,
          seed = 9000L + pos_original
        ),
        error = function(e) {
          message(
            "  ERRO estação ",
            g,
            ": ",
            conditionMessage(e)
          )
          e
        }
      )
      
      # Salva imediatamente depois de cada estação
      saveRDS(fits, f_out)
      
      if (inherits(fits[[g]], "stanfit")) {
        message("  OK — ", g, " concluída e salva.")
      } else {
        message("  ERRO — ", g, " não gerou stanfit.")
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # Status final
  # -------------------------------------------------------------------------
  status <- vapply(
    gauge_list,
    function(g) {
      x <- fits[[g]]
      if (inherits(x, "stanfit")) return("OK")
      if (inherits(x, "error")) return("ERRO")
      "PENDENTE"
    },
    character(1)
  )
  
  n_ok <- sum(status == "OK")
  n_err <- sum(status == "ERRO")
  n_pendente <- sum(status == "PENDENTE")
  
  saveRDS(
    list(
      uf = uf,
      gauge_lookup = ctx$gauge_lookup,
      duration_lookup = duration_lookup,
      D_vals = D_vals,
      D0 = D0_val,
      n_total = length(gauge_list),
      n_ok = n_ok,
      n_err = n_err,
      n_pendente = n_pendente,
      gauge_ok = gauge_list[status == "OK"],
      gauge_erro = gauge_list[status == "ERRO"],
      gauge_pendente = gauge_list[status == "PENDENTE"]
    ),
    f_meta
  )
  
  message(
    "[", uf, "] status final: ok=", n_ok,
    " | erro=", n_err,
    " | pendente=", n_pendente
  )
  
  if (n_ok == length(gauge_list)) {
    message("[", uf, "] UF COMPLETA.")
  } else {
    message(
      "[", uf, "] ainda faltam ",
      n_err + n_pendente,
      " estação(ões)."
    )
  }
}

message("Concluído — GEV non-pooled Uniplu por UF.")
message("Saída: ", DIR_GEV)
