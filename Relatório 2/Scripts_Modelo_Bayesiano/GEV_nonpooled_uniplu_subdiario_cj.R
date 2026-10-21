# =============================================================================
# PASSO 16 — GEV NON-POOLED SUBDIÁRIO CONJUNTO
# 10 min–24 h | INVARIÂNCIA DE ESCALA POR PARTES
# =============================================================================
#
# Um único modelo por estação:
#
#   LONGO: d >= 1 h
#   CURTO: d < 1 h
#
# D0 = 24 h
# D_break = 1 h
#
# A função de escala é contínua em D_break.
#
# Seleção:
#   * 24 durações de 1–24 h com >=8 anos em TODAS
#   * pelo menos 2 durações <1 h com >=8 anos
#
# QC20% já aplicado no PASSO 12.
# Este passo apenas audita o QC.
#
# imax em mm/h | D e theta em horas
# xi_d ∈ [-0.5,0.5]
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(rstan)

# =============================================================================
# 1. CAMINHOS
# =============================================================================

DIR_DADOS_FINAIS <- #Caminho dos dados Utilizados

F_IMAX <- file.path(
  DIR_DADOS_FINAIS,
  "imax_sem_iqr_subdiario_Uniplu.rds"
)

if (!file.exists(F_IMAX)) {
  stop("Arquivo de máximas não encontrado: ", F_IMAX)
}

# Pasta exclusiva para os ajustes com a nova base
DIR_GEV <- #

DIR_FIT <- file.path(
  DIR_GEV,
  "10min_24h_conjunto"
)

# O passo 5 salva o resumo da seleção nesta pasta
DIR_DF <- file.path(DIR_GEV, "selecao")

dir.create(
  DIR_FIT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  DIR_DF,
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# 2. CONFIGURAÇÕES
# =============================================================================

TEST_UF <- NULL
TEST_GAUGE <- NULL

MIN_YEARS <- 8L
QC_MAX_NA <- 0.20

DUR_LONGO_MIN <- (1:24)*60
DUR_CURTO_MIN <- c(10,15,20,30,40,45,50)
MIN_DUR_CURTO <- 2L

D0 <- 24
D_BREAK <- 1

N_CHAINS <- 5L
N_ITER <- 20000L
N_WARMUP <- 10000L
N_THIN <- 50L
ADAPT_DELTA <- 0.999
MAX_TREEDEPTH <- 17L
N_CORES <- 5L

# =============================================================================
# 3. LEITURA E AUDITORIA DA BASE STRICT_SEM_IQR
# =============================================================================

bruto <- readRDS(F_IMAX)

req <- c(
  "gauge_code",
  "state",
  "duration_min",
  "wateryear",
  "na_prct_yr",
  "imax_mm_h",
  "n_valid_years",
  "series_eligible"
)

faltam <- setdiff(req, names(bruto))

if (length(faltam)) {
  stop(
    "Colunas ausentes: ",
    paste(faltam, collapse = ", ")
  )
}

if (
  any(!is.finite(bruto$na_prct_yr)) ||
  any(bruto$na_prct_yr < 0) ||
  any(bruto$na_prct_yr >= QC_MAX_NA)
) {
  stop("Há anos com falhas inválidas ou iguais/superiores a 20%.")
}

if (!all(bruto$series_eligible %in% TRUE)) {
  stop("Há séries marcadas como não elegíveis.")
}

if (
  any(!is.finite(bruto$n_valid_years)) ||
  any(bruto$n_valid_years < MIN_YEARS)
) {
  stop("Há séries com menos de oito anos válidos.")
}

dados <- bruto %>%
  transmute(
    gauge_code = trimws(as.character(gauge_code)),
    estado = as.character(state),
    duration_min = as.integer(duration_min),
    duration_h = as.numeric(duration_min) / 60,
    year = as.integer(wateryear),
    imax = as.numeric(imax_mm_h),
    na_prct_yr = as.numeric(na_prct_yr)
  )

# Liberar o objeto original após a adaptação
rm(bruto)
invisible(gc())

# Mantém o critério do seu modelo atual: máximas positivas
dados <- dados %>%
  filter(
    is.finite(imax),
    imax > 0,
    is.finite(duration_h),
    duration_h > 0,
    !is.na(year),
    !is.na(gauge_code),
    nzchar(gauge_code)
  )

if (!is.null(TEST_UF) && nzchar(TEST_UF)) {
  dados <- dados %>%
    filter(estado == TEST_UF)
  
  message("TEST_UF = ", TEST_UF)
}

dup <- dados %>%
  count(gauge_code, duration_min, year) %>%
  filter(n > 1L)

if (nrow(dup)) {
  stop(
    "Há duplicatas estação–duração–ano: ",
    nrow(dup)
  )
}

message("Auditoria strict_sem_iqr: OK.")

message(
  "Observações disponíveis: ", nrow(dados),
  " | estações: ", n_distinct(dados$gauge_code)
)

# =============================================================================
# 4. SELEÇÃO — BLOCO LONGO
# 24 durações de 1–24 h com >=8 anos em TODAS
# =============================================================================

years_longo <- dados %>%
  filter(duration_min %in% DUR_LONGO_MIN) %>%
  distinct(gauge_code,estado,duration_min,year) %>%
  count(gauge_code,estado,duration_min,name="n_years")

resumo_longo <- years_longo %>%
  group_by(gauge_code,estado) %>%
  summarise(
    n_duracoes_longas=n_distinct(duration_min),
    n_longas_8anos=sum(n_years>=MIN_YEARS),
    min_anos_longo=min(n_years),
    .groups="drop"
  ) %>%
  mutate(
    passa_longo=n_duracoes_longas==24 &
      n_longas_8anos==24
  )

estacoes_longo <- resumo_longo %>%
  filter(passa_longo) %>%
  pull(gauge_code)

message("\nEstações com 1–24 h completas: ",length(estacoes_longo))

# =============================================================================
# 5. SELEÇÃO — RAMO SUB-HORÁRIO
# Pelo menos 2 durações <1 h com >=8 anos
# =============================================================================

years_curto <- dados %>%
  filter(
    gauge_code %in% estacoes_longo,
    duration_min %in% DUR_CURTO_MIN
  ) %>%
  distinct(gauge_code,duration_min,year) %>%
  count(gauge_code,duration_min,name="n_years")

pares_curto_validos <- years_curto %>%
  filter(n_years>=MIN_YEARS)

resumo_curto <- pares_curto_validos %>%
  group_by(gauge_code) %>%
  summarise(
    n_duracoes_curto=n_distinct(duration_min),
    min_anos_curto=min(n_years),
    med_anos_curto=median(n_years),
    max_anos_curto=max(n_years),
    .groups="drop"
  )

resumo_selecao <- resumo_longo %>%
  filter(passa_longo) %>%
  left_join(resumo_curto,by="gauge_code") %>%
  mutate(
    n_duracoes_curto=coalesce(n_duracoes_curto,0L),
    passa_curto=n_duracoes_curto>=MIN_DUR_CURTO,
    passa_modelo=passa_longo & passa_curto
  ) %>%
  arrange(desc(passa_modelo),estado,desc(n_duracoes_curto))

estacoes_aprovadas <- resumo_selecao %>%
  filter(passa_modelo) %>%
  pull(gauge_code)

message("Estações com ramo curto suficiente: ",length(estacoes_aprovadas))

write_excel_csv(
  resumo_selecao,
  file.path(DIR_DF,"filtro_nonpooled_subdiario_10min_24h_conjunto.csv")
)

# =============================================================================
# 6. BASE FINAL DO MODELO
# =============================================================================

dados_longo <- dados %>%
  filter(
    gauge_code %in% estacoes_aprovadas,
    duration_min %in% DUR_LONGO_MIN
  )

dados_curto <- dados %>%
  filter(
    gauge_code %in% estacoes_aprovadas,
    duration_min %in% DUR_CURTO_MIN
  ) %>%
  semi_join(
    pares_curto_validos,
    by=c("gauge_code","duration_min")
  )

dados_modelo <- bind_rows(dados_curto,dados_longo) %>%
  arrange(gauge_code,duration_min,year)

if(!is.null(TEST_GAUGE) && nzchar(TEST_GAUGE)){
  dados_modelo <- dados_modelo %>% filter(gauge_code==TEST_GAUGE)
  message("TEST_GAUGE = ",TEST_GAUGE)
}

message("\n==================================================")
message("BASE DO MODELO CONJUNTO")
message("==================================================")
message("Estações: ",n_distinct(dados_modelo$gauge_code))
message("Observações: ",nrow(dados_modelo))
message("Durações (min): ",
        paste(sort(unique(dados_modelo$duration_min)),collapse=", "))

bloco <- list(
  nome="Non-pooled conjunto 10 min–24 h",
  dados=dados_modelo,
  D0=D0,
  D_break=D_BREAK
)

# =============================================================================
# 7. MODELO STAN
# =============================================================================

modelo_nonpooled_conjunto_stan <- "
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
  int<lower=1,upper=K> duration_id[N_obs];

  vector<lower=0>[K] D;

  real<lower=0> D0;
  real<lower=0> D_break;
}

parameters {
  real alpha;
  real beta;

  real<lower=0,upper=1> H_long;
  real<lower=0,upper=1> H_short;

  real<lower=0,upper=10> theta_long;
  real<lower=0,upper=10> theta_short;

  vector[K] mu_log;
  vector[K] sigma_log;

  vector<lower=-0.5,upper=0.5>[K] xi_d;
  real<lower=-0.5,upper=0.5> xi_bar;

  real<lower=0> tau2_mu;
  real<lower=0> tau2_sigma;
  real<lower=0> tau2_xi;
}

transformed parameters {
  vector[K] eta;
  vector[K] mu;
  vector[K] sigma;

  real eta_break;

  eta_break =
    -H_long *
    log(
      (D_break + theta_long) /
      (D0 + theta_long)
    );

  for (k in 1:K) {

    if (D[k] >= D_break) {

      eta[k] =
        -H_long *
        log(
          (D[k] + theta_long) /
          (D0 + theta_long)
        );

    } else {

      eta[k] =
        eta_break
        -H_short *
        log(
          (D[k] + theta_short) /
          (D_break + theta_short)
        );
    }

    mu[k] = exp(mu_log[k]);
    sigma[k] = exp(sigma_log[k]);
  }
}

model {
  alpha ~ normal(0,1e3);
  beta ~ normal(0,1e3);

  H_long ~ beta(1,1);
  H_short ~ beta(1,1);

  theta_long ~ uniform(0,10);
  theta_short ~ uniform(0,10);

  tau2_mu ~ inv_gamma(0.01,0.01);
  tau2_sigma ~ inv_gamma(0.01,0.01);
  tau2_xi ~ inv_gamma(0.01,0.01);

  {
    real a = -0.5;
    real b = 0.5;
    real p = 9;
    real q = 6;

    real lbetaf =
      lgamma(p) +
      lgamma(q) -
      lgamma(p + q);

    target +=
      -lbetaf
      + (p - 1) * log(xi_bar - a)
      + (q - 1) * log(b - xi_bar)
      - (p + q - 1) * log(b - a);
  }

  for (k in 1:K) {

    mu_log[k] ~ normal(
      alpha + eta[k],
      sqrt(tau2_mu)
    );

    sigma_log[k] ~ normal(
      beta + eta[k],
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

# =============================================================================
# 8. COMPILAÇÃO
# =============================================================================

rstan_options(auto_write=TRUE)
options(mc.cores=N_CORES)

message("\nCompilando modelo conjunto...")
sm_nonpooled_conjunto <- rstan::stan_model(
  model_code=modelo_nonpooled_conjunto_stan
)
message("Stan compilado.")

# =============================================================================
# 9. CONTEXTO DE CADA ESTAÇÃO
# =============================================================================

make_station_context <- function(gauge_code,bloco){
  
  dados_estacao <- bloco$dados %>%
    filter(.data$gauge_code==.env$gauge_code) %>%
    arrange(duration_min,year)
  
  if(!nrow(dados_estacao))
    stop("Estação ",gauge_code," não encontrada.")
  
  duration_lookup <- dados_estacao %>%
    distinct(duration_min,duration_h) %>%
    arrange(duration_min) %>%
    mutate(duration_id=row_number())
  
  dados_estacao <- dados_estacao %>%
    left_join(
      duration_lookup,
      by=c("duration_min","duration_h")
    ) %>%
    arrange(duration_id,year)
  
  D_vals <- duration_lookup$duration_h
  K <- length(D_vals)
  
  n_long <- sum(D_vals>=bloco$D_break)
  n_short <- sum(D_vals<bloco$D_break)
  
  if(n_long!=24)
    stop("Estação ",gauge_code,": possui ",n_long,
         " durações longas; esperadas 24.")
  
  if(n_short<MIN_DUR_CURTO)
    stop("Estação ",gauge_code,": apenas ",n_short,
         " durações sub-horárias.")
  
  if(!any(abs(D_vals-bloco$D0)<1e-10))
    stop("Estação ",gauge_code,": D0=24 h ausente.")
  
  if(!any(abs(D_vals-bloco$D_break)<1e-10))
    stop("Estação ",gauge_code,": D_break=1 h ausente.")
  
  duration_stats <- dados_estacao %>%
    group_by(duration_id) %>%
    summarise(
      median_imax=median(imax),
      sd_imax=sd(imax),
      .groups="drop"
    ) %>%
    arrange(duration_id)
  
  if(any(!is.finite(duration_stats$median_imax)) ||
     any(!is.finite(duration_stats$sd_imax)) ||
     any(duration_stats$median_imax<=0) ||
     any(duration_stats$sd_imax<=0))
    stop("Estação ",gauge_code,": estatísticas iniciais inválidas.")
  
  list(
    gauge_code=gauge_code,
    dados=dados_estacao,
    duration_lookup=duration_lookup,
    D_vals=D_vals,
    K=K,
    duration_stats=duration_stats,
    D0=bloco$D0,
    D_break=bloco$D_break
  )
}

make_stan_data <- function(ctx){
  
  list(
    N_obs=nrow(ctx$dados),
    y=as.numeric(ctx$dados$imax),
    K=ctx$K,
    duration_id=as.integer(ctx$dados$duration_id),
    D=as.numeric(ctx$D_vals),
    D0=as.numeric(ctx$D0),
    D_break=as.numeric(ctx$D_break)
  )
}

# =============================================================================
# 10. INICIALIZAÇÃO
# =============================================================================

make_nonpooled_conjunto_init <- function(ctx){
  
  H_long0 <- runif(1,0.45,0.70)
  H_short0 <- runif(1,0.45,0.70)
  
  theta_long0 <- runif(1,0.8,1.2)
  theta_short0 <- runif(1,0.8,1.2)
  
  eta_break0 <-
    -H_long0 *
    log(
      (ctx$D_break + theta_long0) /
        (ctx$D0 + theta_long0)
    )
  
  eta0 <- ifelse(
    ctx$D_vals >= ctx$D_break,
    
    -H_long0 *
      log(
        (ctx$D_vals + theta_long0) /
          (ctx$D0 + theta_long0)
      ),
    
    eta_break0 -
      H_short0 *
      log(
        (ctx$D_vals + theta_short0) /
          (ctx$D_break + theta_short0)
      )
  )
  
  ref_idx <- which.min(abs(ctx$D_vals-ctx$D0))
  
  alpha0 <- rnorm(
    1,
    log(ctx$duration_stats$median_imax[ref_idx]),
    0.03
  )
  
  beta0 <- rnorm(
    1,
    log(ctx$duration_stats$sd_imax[ref_idx]),
    0.03
  )
  
  xi_bar0 <- runif(1,0.01,0.04)
  
  xi_d0 <- pmin(
    pmax(
      rnorm(ctx$K,xi_bar0,0.02),
      -0.49
    ),
    0.49
  )
  
  list(
    alpha=alpha0,
    beta=beta0,
    
    H_long=H_long0,
    H_short=H_short0,
    
    theta_long=theta_long0,
    theta_short=theta_short0,
    
    mu_log=alpha0+eta0,
    sigma_log=beta0+eta0,
    
    xi_bar=xi_bar0,
    xi_d=xi_d0,
    
    tau2_mu=runif(1,0.02,0.05),
    tau2_sigma=runif(1,0.02,0.05),
    tau2_xi=runif(1,0.0025,0.01)
  )
}

# =============================================================================
# 11. AJUSTE POR ESTAÇÃO
# =============================================================================

fit_station <- function(gauge_code,bloco,seed=123){
  
  ctx <- make_station_context(
    gauge_code,
    bloco
  )
  
  rstan::sampling(
    object=sm_nonpooled_conjunto,
    data=make_stan_data(ctx),
    
    chains=N_CHAINS,
    cores=N_CORES,
    
    iter=N_ITER,
    warmup=N_WARMUP,
    thin=N_THIN,
    
    control=list(
      adapt_delta=ADAPT_DELTA,
      max_treedepth=MAX_TREEDEPTH
    ),
    
    seed=seed,
    
    init=function()
      make_nonpooled_conjunto_init(ctx),
    
    refresh=50
  )
}

# =============================================================================
# 12. SALVAMENTO ATÔMICO
# =============================================================================

save_fit_atomic <- function(obj,arquivo){
  
  tmp <- paste0(arquivo,".tmp")
  
  if(file.exists(tmp))
    file.remove(tmp)
  
  saveRDS(obj,tmp)
  
  if(!file.exists(tmp) ||
     file.info(tmp)$size<=0)
    stop("Falha ao salvar: ",tmp)
  
  if(file.exists(arquivo))
    file.remove(arquivo)
  
  if(!file.rename(tmp,arquivo))
    stop("Não foi possível mover fit para: ",arquivo)
  
  invisible(TRUE)
}

# =============================================================================
# 13. LOOP COM RETOMADA
# =============================================================================

run_model <- function(
    bloco,
    dir_fit,
    inicio = 1L,
    fim = Inf
) {
  
  # Lista completa, sempre na mesma ordem
  gauge_list_completa <- sort(
    unique(bloco$dados$gauge_code)
  )
  
  total <- length(gauge_list_completa)
  
  if (!total) {
    stop("Nenhuma estação disponível.")
  }
  
  if (
    length(inicio) != 1L ||
    !is.finite(inicio) ||
    inicio < 1 ||
    inicio != floor(inicio) ||
    inicio > total
  ) {
    stop("Posição inicial inválida. Total disponível: ", total)
  }
  
  if (
    length(fim) != 1L ||
    is.na(fim) ||
    fim < inicio ||
    (is.finite(fim) && fim != floor(fim))
  ) {
    stop("Posição final inválida.")
  }
  
  fim <- as.integer(min(fim, total))
  inicio <- as.integer(inicio)
  
  gauge_list <- gauge_list_completa[
    seq.int(inicio, fim)
  ]
  
  dir.create(
    dir_fit,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  message(
    "Lote: posições ", inicio, "–", fim,
    " | estações neste lote: ", length(gauge_list),
    " | total aprovado: ", total
  )
  
  # Salva a identificação das estações deste lote
  write.csv(
    data.frame(
      posicao_global = seq.int(inicio, fim),
      gauge_code = gauge_list
    ),
    file.path(dir_fit, "estacoes_do_lote.csv"),
    row.names = FALSE
  )
  
  F_STATUS <- file.path(
    dir_fit,
    "status.csv"
  )
  
  if(file.exists(F_STATUS)){
    
    status <- read.csv(
      F_STATUS,
      stringsAsFactors=FALSE
    )
    
    status$gauge_code <-
      as.character(status$gauge_code)
    
  } else {
    
    status <- data.frame(
      gauge_code=gauge_list,
      posicao=seq_along(gauge_list),
      status="PENDENTE",
      data_hora=NA_character_,
      erro=NA_character_
    )
  }
  
  novas <- setdiff(
    gauge_list,
    status$gauge_code
  )
  
  if(length(novas)){
    
    status <- bind_rows(
      status,
      data.frame(
        gauge_code=novas,
        posicao=match(novas,gauge_list),
        status="PENDENTE",
        data_hora=NA_character_,
        erro=NA_character_
      )
    )
  }
  
  fit_file <- function(g)
    file.path(
      dir_fit,
      paste0(g,"_fit.rds")
    )
  
  ok_file <- function(g)
    file.path(
      dir_fit,
      paste0(g,"_fit.ok")
    )
  
  ja_ok <- vapply(
    gauge_list,
    function(g){
      
      f <- fit_file(g)
      
      file.exists(ok_file(g)) &&
        file.exists(f) &&
        file.info(f)$size>0
    },
    logical(1)
  )
  
  faltando <- gauge_list[!ja_ok]
  
  message("\n==================================================")
  message(bloco$nome)
  message("Estações previstas: ",length(gauge_list))
  message("Já concluídas: ",sum(ja_ok))
  message("Faltando: ",length(faltando))
  message("==================================================")
  
  for(j in seq_along(faltando)){
    
    g <- faltando[[j]]
    # Preserva a seed da estação independentemente do lote
    pos <- match(g, gauge_list_completa)
    
    message(sprintf(
      "[%d/%d do lote] estação %s | posição global %d/%d",
      j,
      length(faltando),
      g,
      pos,
      length(gauge_list_completa)
    ))
    
    fit_g <- tryCatch(
      fit_station(
        gauge_code=g,
        bloco=bloco,
        seed=9000L+pos
      ),
      error=function(e) e
    )
    
    idx <- match(
      g,
      status$gauge_code
    )
    
    if(inherits(fit_g,"stanfit")){
      
      f_fit <- fit_file(g)
      
      save_fit_atomic(
        fit_g,
        f_fit
      )
      
      file.create(
        ok_file(g)
      )
      
      status$status[idx] <- "OK"
      status$data_hora[idx] <-
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      
      status$erro[idx] <-
        NA_character_
      
      message("  OK — ",basename(f_fit))
      
    } else {
      
      erro <- conditionMessage(fit_g)
      
      status$status[idx] <- "ERRO"
      status$data_hora[idx] <-
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      
      status$erro[idx] <- erro
      
      message("  ERRO — ",erro)
    }
    
    write.csv(
      status,
      F_STATUS,
      row.names=FALSE
    )
    
    rm(fit_g)
    invisible(gc())
  }
  
  n_ok <- sum(
    vapply(
      gauge_list,
      function(g)
        file.exists(ok_file(g)) &&
        file.exists(fit_file(g)),
      logical(1)
    )
  )
  
  message(
    "\nConcluídas: ",
    n_ok,
    "/",
    length(gauge_list)
  )
  
  invisible(status)
}

# =============================================================================
# 14. EXECUÇÃO
# =============================================================================

status_conjunto <- run_model(
  bloco = bloco,
  dir_fit = file.path(DIR_FIT, "lote_101_200"),
  inicio = 101L,
  fim = 200L
)

message("\n==================================================")
message("PASSO 16 CONCLUÍDO")
message("==================================================")
message("Modelo conjunto 10 min–24 h: ",DIR_FIT)

