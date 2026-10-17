# =============================================================================
# COMPARAÇÃO DOS MODELOS BAYESIANOS REGIONAIS — RIO GRANDE DO SUL
#
# Objetivos:
# 1. Comparar o desempenho preditivo dos cinco modelos Bayesianos;
# 2. Comparar a estrutura posterior dos parâmetros;
# 3. Avaliar o efeito da estrutura de xi sobre as estimativas IDF;
# 4. Selecionar estações representativas para análise detalhada.
# =============================================================================


# =============================================================================
# 0. PACOTES
# =============================================================================

suppressPackageStartupMessages({
  
  library(rstan)
  library(loo)
  library(posterior)
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  
  library(ggplot2)
  library(patchwork)
  
  library(sf)
  library(geobr)
  library(arrow)
  
})


# =============================================================================
# 1. CAMINHOS
# =============================================================================

DIR_BASE <- "C:/Users/Adm/OneDrive/3. UFC/UFC - 2026/Art. 01_BHM"
DIR_SAVE <- "C:/Larissa"

DIR_RESULTADOS <- file.path(DIR_SAVE, "comparacao_modelos_bayesianos_RS")
DIR_TABELAS <- file.path(DIR_RESULTADOS, "tabelas")
DIR_FIGURAS <- file.path(DIR_RESULTADOS, "figuras")
DIR_GPKG <- file.path(DIR_RESULTADOS, "geopackage")

dir.create(DIR_RESULTADOS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABELAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIGURAS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_GPKG, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. ARQUIVOS DOS CINCO MODELOS
# =============================================================================

ARQUIVOS_MODELOS <- c(nonpooled = file.path( DIR_BASE, "modelo_nonpooled_xi_j_stan_metodologia_comloo.rds"),
                      bhm_xi_d = file.path( DIR_BASE, "fit_bhm_xi_d_metodologia_novo.rds" ),
                      bhm_xi_j = file.path( DIR_BASE, "fit_bhm_xi_j_metodologia_com loo.rds" ),
                      bhm_xi_jd = file.path( DIR_BASE, "fit_bhm_xi_jd_metodologia_com loo.rds"),
                      bhm_common = file.path( DIR_BASE, "fit_bhm_xi_common_metodologia_comloo.rds" ) )


# Verificação imediata
arquivos_ausentes <- ARQUIVOS_MODELOS[ !file.exists(ARQUIVOS_MODELOS) ]

if (length(arquivos_ausentes) > 0L) {
  stop("Os seguintes arquivos não foram encontrados:\n", paste( arquivos_ausentes, collapse = "\n" ))
}

F_METADATA <- file.path(DIR_BASE, "df_subdaily_info.pqt")

# =============================================================================
# 3. IDENTIFICAÇÃO DOS MODELOS
# =============================================================================

ORDEM_MODELOS <- c("bhm_common", "bhm_xi_d", "bhm_xi_j", "bhm_xi_jd", "nonpooled")
NOMES_MODELOS <- c(
  "bhm_common" = "Common ξ",
  "bhm_xi_d"   = "Duration-specific ξ",
  "bhm_xi_j"   = "Station-specific ξ",
  "bhm_xi_jd"  = "Station-duration ξ",
  "nonpooled"  = "Non-pooled"
)

# =============================================================================
# 4. PALETA DOS MODELOS
# =============================================================================

CORES_MODELOS <- c(
  "Common ξ"            = "#D55E00",
  "Duration-specific ξ" = "#E69F00",
  "Station-specific ξ"  = "#009E73",
  "Station-duration ξ"  = "#0072B2",
  "Non-pooled"          = "#CC79A7"
)
TR_GRID <- c( 2, 5, 10, 25, 50, 100, 200 )

# =============================================================================
# 5. PADRÃO GRÁFICO
# =============================================================================

font <- "sans"
fontsize <- 8
textsize <- fontsize / .pt

tema_relatorio <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    
    axis.ticks = element_line(
      linewidth = 0.3
    ),
    
    text = element_text(
      family = font,
      size = fontsize
    ),
    
    axis.text = element_text(
      family = font,
      size = fontsize
    ),
    
    axis.title = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.background = element_blank(),
    
    legend.key = element_blank(),
    
    legend.key.height = unit(
      fontsize,
      "pt"
    ),
    
    legend.title = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.text = element_text(
      family = font,
      size = fontsize
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = font,
      size = fontsize
    ),
    
    panel.spacing.y = unit(
      0.5,
      "pt"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

# =============================================================================
# 6. FUNÇÕES GERAIS
# =============================================================================

limpar_codigo_estacao <- function(x) {
  
  x %>%
    as.character() %>%
    str_remove("^root\\$") %>%
    str_remove("^obj\\$") %>%
    str_trim()
}

rotular_modelo <- function(x) {
  factor( NOMES_MODELOS[ as.character(x) ], levels = unname( NOMES_MODELOS[ ORDEM_MODELOS ] ) )
}

salvar_tabela <- function( objeto, nome ) {
  readr::write_excel_csv(objeto, file.path( DIR_TABELAS, nome ) )
}

salvar_figura <- function( figura, nome, largura = 15, altura = 10 ) {
  
  ggsave( filename = file.path( DIR_FIGURAS, nome ), plot = figura,
    width = largura, height = altura, units = "cm", dpi = 600 )
}

# =============================================================================
# 7. LEITURA E ORGANIZAÇÃO DOS AJUSTES STAN
# =============================================================================

is_stanfit <- function(x) {
  methods::is(x, "stanfit")
}

find_all_stanfits <- function( x, path = "root" ) {
  # Caso o próprio objeto já seja um stanfit
  if (is_stanfit(x)) {
    return( setNames( list(x), path ) )
  }
  
  # Se não for lista
  if (!is.list(x)) {
    return(list())
  }
  
  # Nomes dos elementos da lista
  nms <- names(x)
  if (is.null(nms)) {
    nms <- as.character( seq_along(x) )
  }
  
  # Procurar recursivamente
  out <- list()
  
  for (i in seq_along(x)) {
    
    new_path <- paste0( path, "$", nms[i] )
    out <- c( out, find_all_stanfits( x[[i]], path = new_path ) )
  }
  out
}

objetos_modelos <- purrr::map( ARQUIVOS_MODELOS, readRDS )
fits_encontrados <- purrr::map( objetos_modelos, find_all_stanfits)

n_fits_modelo <- purrr::map_int( fits_encontrados, length )

tabela_estrutura_modelos <- tibble( modelo = names(n_fits_modelo),
                                    nome_modelo = unname( NOMES_MODELOS[ names(n_fits_modelo) ] ),
                                    n_stanfit = unname( n_fits_modelo ) )
cat( "\nNúmero de objetos stanfit encontrados:\n" )
print( tabela_estrutura_modelos )

MODELOS_REGIONAIS <- c(
  "bhm_common",
  "bhm_xi_d",
  "bhm_xi_j",
  "bhm_xi_jd"
)

# Cada modelo regional deve ter um único ajuste conjunto

for (m in MODELOS_REGIONAIS) {
  
  if (length(fits_encontrados[[m]]) != 1L) {
    
    stop("Esperava exatamente 1 stanfit para o modelo ", m, ", mas encontrei ", length(fits_encontrados[[m]]), "." )
  }
}
# O non-pooled deve possuir vários ajustes, um por estação
if (length(fits_encontrados$nonpooled) == 0L) {
  stop( "Nenhum stanfit foi encontrado no modelo non-pooled." )
}

fits_regionais <- purrr::map( MODELOS_REGIONAIS, function(m) {
    fits_encontrados[[m]][[1]]
  }
)
names(fits_regionais) <- MODELOS_REGIONAIS
fits_nonpooled <- fits_encontrados$nonpooled

nomes_nonpooled_originais <- names(fits_nonpooled)
nomes_nonpooled_limpos <- limpar_codigo_estacao( nomes_nonpooled_originais)

# Em alguns RDS pode haver caminhos como:
# root$A801
# obj$A801
# root$obj$A801
#
# Pegamos sempre a última parte depois de "$".

nomes_nonpooled_limpos <- stringr::str_replace( nomes_nonpooled_limpos, "^.*\\$", "" )
nomes_nonpooled_limpos <- stringr::str_trim( nomes_nonpooled_limpos )

# Verificar duplicatas após a limpeza

if (anyDuplicated(nomes_nonpooled_limpos)) {
  
  duplicados <- unique( nomes_nonpooled_limpos[ duplicated( nomes_nonpooled_limpos ) ] )
  stop(
    "Foram encontrados códigos duplicados após limpar os nomes do non-pooled: ",
    paste( duplicados, collapse = ", " ) )
}
names(fits_nonpooled) <- nomes_nonpooled_limpos

gauge_codes <- sort( names( fits_nonpooled ) )
cat( "\nEstações no modelo non-pooled: ", length(gauge_codes), "\n", sep = "" )
cat( "Primeiras estações: ", paste( head( gauge_codes, 10 ), collapse = ", " ), "\n", sep = "" )

metadata_raw <- arrow::read_parquet( F_METADATA )
meta_estacoes <- metadata_raw %>%
  transmute( gauge_code = limpar_codigo_estacao( gauge_code ),
             long = as.numeric( long ),
             lat = as.numeric( lat ),
             elevation = if ( "elevation" %in% names(metadata_raw) ) { 
               as.numeric(elevation)
               } else {
                 NA_real_
                 },
             network = if ( "network" %in% names(metadata_raw) ) {
               as.character(network)
               } else {
                 NA_character_
                 },
             responsible = if ( "responsible" %in% names(metadata_raw) ) {
               as.character(responsible)
               } else {
                 NA_character_
                 }
            ) %>%
  filter( !is.na(gauge_code), is.finite(long), is.finite(lat) ) %>%
  distinct( gauge_code, .keep_all = TRUE )

estacoes_sem_metadata <- setdiff( gauge_codes, meta_estacoes$gauge_code )

if (length(estacoes_sem_metadata) > 0L) {
  warning(
    length(estacoes_sem_metadata),
    " estações do non-pooled não possuem coordenadas nos metadados."
  )
}

meta_estacoes_modelo <- meta_estacoes %>%
  filter( gauge_code %in% gauge_codes ) %>%
  arrange( gauge_code )
cat( "\nEstações com metadados válidos: ", nrow(meta_estacoes_modelo), " / ", length(gauge_codes), "\n", sep = "" )

duracoes_h <- 1:24
duracoes_por_modelo <- purrr::map(
  fits_regionais,
  ~ duracoes_h
)
# Verificar se todos os modelos regionais usam a mesma grade

duracoes_iguais <- all( purrr::map_lgl( duracoes_por_modelo, ~ identical( as.numeric(.x),
      as.numeric( duracoes_por_modelo[[1]] ) ) ) )

if (!duracoes_iguais) {
  stop("Os modelos regionais não utilizam a mesma grade de durações.")
}

D_GRID <- as.numeric(duracoes_por_modelo[[1]])

cat("\nDurações encontradas nos modelos:\n", paste( D_GRID, collapse = ", " ), " h\n", sep = "" )

library(rstan)
library(loo)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

diagnostico_mcmc_fit <- function(fit, modelo) {
  
  sum_fit <- summary(fit)$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parametro")
  
  # retirar lp__
  sum_pars <- sum_fit |>
    filter(parametro != "lp__")
  
  # Rhat
  n_rhat_101 <- sum( sum_pars$Rhat > 1.01, na.rm = TRUE )
  
  n_rhat_105 <- sum( sum_pars$Rhat > 1.05, na.rm = TRUE )
  
  rhat_max <- max( sum_pars$Rhat, na.rm = TRUE )
  
  # ESS
  ess_min <- min( sum_pars$n_eff, na.rm = TRUE )
  
  ess_mediana <- median( sum_pars$n_eff, na.rm = TRUE )
  
  # divergências
  sampler_params <- rstan::get_sampler_params( fit, inc_warmup = FALSE )
  
  divergencias <- sum( vapply( sampler_params, function(x) sum(x[, "divergent__"]), numeric(1) ) )
  
  tibble(
    modelo = modelo,
    rhat_max = rhat_max,
    n_rhat_maior_1.01 = n_rhat_101,
    n_rhat_maior_1.05 = n_rhat_105,
    ess_min = ess_min,
    ess_mediana = ess_mediana,
    divergencias = divergencias
  )
}

diagnostico_mcmc_regionais <- purrr::imap_dfr(
  fits_regionais,
  diagnostico_mcmc_fit
)

diagnostico_mcmc_regionais

diagnostico_mcmc_nonpooled_estacoes <- purrr::imap_dfr(
  fits_nonpooled,
  function(fit, gauge_code) {
    
    diagnostico_mcmc_fit(
      fit = fit,
      modelo = "nonpooled"
    ) %>%
      mutate(
        gauge_code = gauge_code,
        .after = modelo
      )
  }
)

diagnostico_mcmc_nonpooled_estacoes

diagnostico_mcmc <- bind_rows(
  diagnostico_mcmc_regionais,
  diagnostico_mcmc_nonpooled_estacoes
) %>%
  mutate(
    modelo = factor(
      modelo,
      levels = ORDEM_MODELOS
    )
  ) %>%
  arrange(modelo)

diagnostico_mcmc

salvar_tabela( diagnostico_mcmc, "diagnostico_mcmc_modelos.csv" )
salvar_tabela( diagnostico_mcmc_nonpooled_estacoes, "diagnostico_mcmc_nonpooled_por_estacao.csv" )

flags_diagnostico_nonpooled <- diagnostico_mcmc_nonpooled_estacoes %>%
  transmute(
    gauge_code,
    
    flag_divergencia = divergencias > 0,
    
    flag_rhat_105 = rhat_max > 1.05,
    
    flag_diagnostico =
      divergencias > 0 |
      rhat_max > 1.05,
    
    n_divergencias = divergencias,
    
    rhat_max_estacao = rhat_max,
    
    ess_min_estacao = ess_min
  )

diagnostico_estacoes_regionais <- purrr::imap_dfr(fits_regionais,
  function(fit, modelo) {

    pars_estacao <- c("alpha", "beta", "H", "theta")
    
    pars_disponiveis <- fit@sim$pars_oi
    
    candidatos_xi <- grep( "^xi", pars_disponiveis, value = TRUE)
    
    if (
      modelo %in% c(
        "bhm_xi_j",
        "bhm_xi_jd"
      ) &&
      length(candidatos_xi) > 0
    ) {
      
      # usar o primeiro parâmetro xi encontrado
      pars_estacao <- c(
        pars_estacao,
        candidatos_xi[1]
      )
    }
    
    resultados <- purrr::map_dfr( pars_estacao, function(par) {
        
        sm <- summary( fit, pars = par )$summary |>
          as.data.frame() |>
          tibble::rownames_to_column("parametro")
        
        sm |>
          filter(parametro != "lp__") |>
          mutate( parametro_base = par, modelo = modelo )
      }
    )
    
    # extrair station_id dos índices dos parâmetros
    resultados <- resultados |>
      mutate( station_id = case_when(
          
          # alpha[1], beta[1], H[1], theta[1], xi[1]
          stringr::str_detect( parametro, "\\[[0-9]+\\]$" ) ~ as.integer(
            stringr::str_match( parametro, "\\[([0-9]+)\\]" )[, 2] ),
          
          # xi[estacao,duracao]
          stringr::str_detect( parametro, "\\[[0-9]+,[0-9]+\\]$" ) ~ as.integer(
            stringr::str_match( parametro, "\\[([0-9]+)," )[, 2] ),
          
          TRUE ~ NA_integer_
        )
      ) |>
      filter(!is.na(station_id))
    
    # resumo por estação
    resultados |>
      group_by( modelo, station_id ) |>
      summarise( rhat_max_estacao = max( Rhat, na.rm = TRUE ),
                 ess_min_estacao = min( n_eff, na.rm = TRUE ),
                 flag_rhat_105 = any( Rhat > 1.05, na.rm = TRUE ), .groups = "drop" )
  }
)

divergencias_regionais <- purrr::imap_dfr( fits_regionais,
  function(fit, modelo) {
    sampler_params <- rstan::get_sampler_params( fit, inc_warmup = FALSE )
    n_div <- sum( vapply( sampler_params, function(x) {
          sum(x[, "divergent__"]) }, numeric(1) )
    )
    
    tibble(
      modelo = modelo,
      n_divergencias_modelo = n_div,
      flag_divergencia_modelo = n_div > 0
    )
  }
)

flags_diagnostico_regionais <- diagnostico_estacoes_regionais |>
  left_join(
    divergencias_regionais,
    by = "modelo"
  ) |>
  mutate(
    flag_diagnostico =
      flag_rhat_105 |
      flag_divergencia_modelo
  )

# =============================================================================
# 9. PSIS-LOO
# =============================================================================

cat("\nMODELOS REGIONAIS\n")

purrr::imap( fits_regionais, function(fit, modelo) {
    
    cat( modelo, ": ", paste( grep( "log_lik", fit@sim$pars_oi, value = TRUE ), collapse = ", " ), "\n" )
  }
)

cat("\nNON-POOLED — primeiras estações\n")

purrr::iwalk(
  fits_nonpooled[1:min(5, length(fits_nonpooled))],
  function(fit, gauge_code) {
    
    cat( gauge_code, ": ", paste( grep( "log_lik", fit@sim$pars_oi, value = TRUE ), collapse = ", " ), "\n" )
  }
)

extrair_log_lik <- function(fit) {
  
  pars <- fit@sim$pars_oi
  
  nome_loglik <- grep( "^log_lik$", pars, value = TRUE )
  
  if (length(nome_loglik) != 1L) {
    
    stop( "Não encontrei exatamente um parâmetro chamado log_lik. ", "Encontrados: ",
      paste( grep("log_lik", pars, value = TRUE), collapse = ", " ) )
  }
  
  x <- rstan::extract( fit, pars = nome_loglik, permuted = TRUE )[[nome_loglik]]
  
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  x
}

dim_loglik_regionais <- purrr::imap_dfr(
  fits_regionais,
  function(fit, modelo) {
    
    x <- extrair_log_lik(fit)
    
    tibble(
      modelo = modelo,
      n_draws = nrow(x),
      n_obs = ncol(x)
    )
  }
)

dim_loglik_regionais

dim_loglik_nonpooled <- purrr::imap_dfr(
  fits_nonpooled,
  function(fit, gauge_code) {
    
    x <- extrair_log_lik(fit)
    
    tibble(
      gauge_code = gauge_code,
      n_draws = nrow(x),
      n_obs = ncol(x)
    )
  }
)

dim_loglik_nonpooled

dim_loglik_nonpooled %>%
  summarise(
    n_estacoes = n(),
    min_draws = min(n_draws),
    max_draws = max(n_draws),
    total_obs = sum(n_obs),
    min_obs_estacao = min(n_obs),
    max_obs_estacao = max(n_obs)
  )

n_obs_regionais <- unique(
  dim_loglik_regionais$n_obs
)

n_obs_nonpooled <- sum(
  dim_loglik_nonpooled$n_obs
)

cat(
  "\nObservações modelos regionais:",
  paste(n_obs_regionais, collapse = ", "),
  "\nObservações non-pooled:",
  n_obs_nonpooled,
  "\n"
)

loo_regionais <- purrr::imap(
  fits_regionais,
  function(fit, modelo) {
    
    message("Calculando LOO: ", modelo)
    
    log_lik <- extrair_log_lik(fit)
    
    loo::loo(
      log_lik
    )
  }
)

loo_nonpooled_estacoes <- purrr::imap(
  fits_nonpooled,
  function(fit, gauge_code) {
    
    message(
      "LOO non-pooled: ",
      gauge_code
    )
    
    log_lik <- extrair_log_lik(fit)
    
    loo::loo(
      log_lik
    )
  }
)

diagnostico_loo_nonpooled_estacoes <- purrr::imap_dfr(
  loo_nonpooled_estacoes,
  function(x, gauge_code) {
    
    k <- x$diagnostics$pareto_k
    
    tibble( modelo = "nonpooled", gauge_code = gauge_code, n_obs = length(k),
            pareto_k_max = max( k, na.rm = TRUE ),
            n_k_maior_05 = sum( k > 0.5, na.rm = TRUE ),
            n_k_maior_07 = sum( k > 0.7, na.rm = TRUE ),
            n_k_maior_10 = sum( k > 1, na.rm = TRUE ),
            flag_pareto_07 = any( k > 0.7, na.rm = TRUE ),
            flag_pareto_10 = any( k > 1, na.rm = TRUE )
    )
  }
)

diagnostico_loo_nonpooled_estacoes

diagnostico_loo_nonpooled_estacoes %>%
  filter(flag_pareto_07) %>%
  arrange(desc(pareto_k_max))

table(dim_loglik_nonpooled$n_draws)

log_lik_nonpooled_lista <- purrr::map(
  fits_nonpooled,
  extrair_log_lik
)

log_lik_nonpooled <- do.call(
  cbind,
  log_lik_nonpooled_lista
)

dim(log_lik_nonpooled)

loo_nonpooled <- loo::loo(
  log_lik_nonpooled
)

loo_nonpooled

loo_modelos <- c(loo_regionais, list( nonpooled = loo_nonpooled ) )
names(loo_modelos)

tabela_loo <- purrr::imap_dfr(
  loo_modelos,
  function(x, modelo) {
    
    tibble( modelo = modelo,
            elpd_loo = x$estimates[ "elpd_loo", "Estimate" ],
            se_elpd_loo = x$estimates[ "elpd_loo", "SE" ],
            p_loo = x$estimates[ "p_loo", "Estimate" ],
            se_p_loo = x$estimates[ "p_loo", "SE" ],
            looic = x$estimates[ "looic", "Estimate" ],
            se_looic = x$estimates[ "looic", "SE" ]
            )
  }) %>%
  mutate( nome_modelo = unname( NOMES_MODELOS[modelo] ), .after = modelo )

comparacao_loo <- loo::loo_compare( loo_modelos )
comparacao_loo

comparacao_loo_df <- comparacao_loo %>%
  as.data.frame() %>%
  tibble::rownames_to_column( "modelo" ) %>%
  as_tibble() %>%
  mutate( nome_modelo = unname( NOMES_MODELOS[modelo] ), .after = modelo )

comparacao_loo_df

diagnostico_pareto_modelos <- purrr::imap_dfr(
  loo_modelos,
  function(x, modelo) {
    
    k <- x$diagnostics$pareto_k
    
    tibble(modelo = modelo,
           nome_modelo = unname( NOMES_MODELOS[modelo] ),
           n_obs = length(k), pareto_k_max = max( k, na.rm = TRUE ),
           n_k_maior_05 = sum( k > 0.5, na.rm = TRUE ),
           n_k_maior_07 = sum( k > 0.7, na.rm = TRUE ),
           n_k_maior_10 = sum( k > 1, na.rm = TRUE ),
           pct_k_maior_07 = 100 * mean( k > 0.7, na.rm = TRUE ) )
  }
)

diagnostico_pareto_modelos

salvar_tabela( tabela_loo, "loo_modelos.csv" )
salvar_tabela( comparacao_loo_df, "comparacao_loo_modelos.csv" )
salvar_tabela( diagnostico_pareto_modelos, "diagnostico_pareto_modelos.csv" )
salvar_tabela( diagnostico_loo_nonpooled_estacoes, "diagnostico_loo_nonpooled_por_estacao.csv" )
saveRDS( loo_modelos, file.path( DIR_RESULTADOS, "loo_modelos.rds" ) )

# =============================================================================
# 10. Parametros dos modelos
# =============================================================================

parametros_disponiveis <- function(fit) {
  fit@sim$pars_oi
}
extrair_primeiro_parametro <- function(fit, candidatos) {
  
  disponiveis <- parametros_disponiveis(fit)
  encontrado <- candidatos[candidatos %in% disponiveis]
  
  if (length(encontrado) == 0) {
    return(NULL)
  }
  nome <- encontrado[1]
  
  rstan::extract( fit, pars = nome, permuted = TRUE )[[nome]]
}

normalizar_estacao <- function(x, J) {
  
  if (is.null(x)) return(NULL)
  
  dx <- dim(x)
  
  # vetor: um valor por draw
  if (is.null(dx)) {
    
    return( matrix( rep(x, J), nrow = length(x), ncol = J )
    )
  }
  
  # draw x estação
  if (length(dx) == 2 && dx[2] == J) {
    return(x)
  }
  
  stop("Não consegui interpretar parâmetro por estação. Dimensões: ",
    paste(dx, collapse = " x ")
  )
}

normalizar_jk <- function(x, J, K) {
  
  if (is.null(x)) return(NULL)
  
  dx <- dim(x)
  
  if (is.null(dx)) {
    stop("Esperava um parâmetro com dimensões estação × duração.")
  }
  
  # draw x estação x duração
  if ( length(dx) == 3 && dx[2] == J && dx[3] == K ) {
    return(x)
  }
  
  # draw x duração x estação
  if ( length(dx) == 3 && dx[2] == K && dx[3] == J ) {
    return(aperm(x, c(1, 3, 2)))
  }
  
  stop( "Dimensões inesperadas para parâmetro estação-duração: ",
    paste(dx, collapse = " x ")
  )
}

obter_J <- function(fit) {
  
  H <- extrair_primeiro_parametro( fit, c("H", "H_j") )
  
  if (is.null(H)) {
    stop("Não encontrei H ou H_j no modelo.")
  }
  
  if (is.null(dim(H))) {
    stop("H não possui dimensão por estação.")
  }
  
  ncol(H)
}

extrair_parametros_estruturais <- function(fit) {
  
  J <- obter_J(fit)
  
  alpha <- extrair_primeiro_parametro( fit, c("alpha"))
  beta <- extrair_primeiro_parametro( fit, c("beta"))
  H <- extrair_primeiro_parametro( fit, c("H"))
  theta <- extrair_primeiro_parametro( fit, c("theta"))
  
  list(
    J = J,
    alpha = normalizar_estacao(alpha, J),
    beta = normalizar_estacao(beta, J),
    H = normalizar_estacao(H, J),
    theta = normalizar_estacao(theta, J)
  )
}

purrr::imap(
  fits_regionais,
  function(fit, modelo) {
    
    p <- extrair_parametros_estruturais(fit)
    
    cat(
      modelo,
      "| alpha:", paste(dim(p$alpha), collapse = " x "),
      "| beta:", paste(dim(p$beta), collapse = " x "),
      "| H:", paste(dim(p$H), collapse = " x "),
      "| theta:", paste(dim(p$theta), collapse = " x "),
      "\n"
    )
  }
)

resumir_parametro_estacao <- function(x, nome_parametro) {
  
  tibble( station_id = seq_len(ncol(x)), parametro = nome_parametro,
          mediana = apply( x, 2, median, na.rm = TRUE ),
          q025 = apply( x, 2, quantile, probs = 0.025, na.rm = TRUE ),
          q975 = apply( x, 2, quantile, probs = 0.975, na.rm = TRUE ) )
}

resumir_estruturais_modelo <- function(fit) {
  
  p <- extrair_parametros_estruturais(fit)
  
  bind_rows( resumir_parametro_estacao( p$alpha, "alpha" ),
             resumir_parametro_estacao( p$beta, "beta" ),
             resumir_parametro_estacao( p$H, "H" ),
             resumir_parametro_estacao( p$theta, "theta" )
             )
}

parametros_estruturais_regionais <- purrr::imap_dfr(
  fits_regionais,
  function(fit, modelo) {
    
    resumir_estruturais_modelo(fit) |>
      mutate(
        modelo = modelo,
        .before = 1
      )
  }
)

nrow(parametros_estruturais_regionais)

mapa_station_id <- tibble(
  station_id = seq_along(gauge_codes),
  gauge_code = gauge_codes
)

parametros_estruturais_regionais <-
  parametros_estruturais_regionais |>
  left_join(
    mapa_station_id,
    by = "station_id"
  ) |>
  relocate(
    gauge_code,
    .after = station_id
  )

resumir_vetor_posterior <- function(x, parametro) {
  
  tibble(
    parametro = parametro,
    mediana = median(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE)
  )
}

extrair_estruturais_nonpooled <- function(fit) {
  
  alpha <- extrair_primeiro_parametro(fit, c("alpha"))
  beta <- extrair_primeiro_parametro( fit, c("beta"))
  H <- extrair_primeiro_parametro( fit, c("H"))
  theta <- extrair_primeiro_parametro( fit, c("theta"))
  
  bind_rows(
    resumir_vetor_posterior(alpha, "alpha"),
    resumir_vetor_posterior(beta, "beta"),
    resumir_vetor_posterior(H, "H"),
    resumir_vetor_posterior(theta, "theta")
  )
}

parametros_estruturais_nonpooled <- purrr::imap_dfr(
  fits_nonpooled,
  function(fit, gauge_code) {
    
    extrair_estruturais_nonpooled(fit) |>
      mutate(
        modelo = "nonpooled",
        gauge_code = gauge_code,
        .before = 1
      )
  }
)

mapa_station_id <- tibble(
  station_id = seq_along(gauge_codes),
  gauge_code = gauge_codes
)

parametros_estruturais_regionais <-
  parametros_estruturais_regionais |>
  select(
    -any_of("gauge_code")
  ) |>
  left_join(
    mapa_station_id,
    by = "station_id"
  ) |>
  relocate(
    gauge_code,
    .after = station_id
  )

parametros_estruturais_nonpooled <-
  parametros_estruturais_nonpooled |>
  mutate(
    station_id = match(
      gauge_code,
      gauge_codes
    )
  ) |>
  relocate(
    station_id,
    .after = modelo
  )

parametros_estruturais_regionais <-
  parametros_estruturais_regionais |>
  left_join(
    flags_diagnostico_regionais,
    by = c(
      "modelo",
      "station_id"
    )
  )

parametros_estruturais_nonpooled <-
  parametros_estruturais_nonpooled |>
  left_join(
    flags_diagnostico_nonpooled,
    by = "gauge_code"
  )

parametros_regionais_base <-
  parametros_estruturais_regionais |>
  select(
    modelo,
    station_id,
    gauge_code,
    parametro,
    mediana,
    q025,
    q975
  )

parametros_nonpooled_base <-
  parametros_estruturais_nonpooled |>
  select(
    modelo,
    station_id,
    gauge_code,
    parametro,
    mediana,
    q025,
    q975
  )

parametros_estruturais_todos <- bind_rows(
  parametros_regionais_base,
  parametros_nonpooled_base
) |>
  mutate(
    nome_modelo = rotular_modelo(modelo)
  )

salvar_tabela(parametros_estruturais_todos, "parametros_estruturais_cinco_modelos.csv")
saveRDS( parametros_estruturais_todos, file.path( DIR_RESULTADOS, "parametros_estruturais_cinco_modelos.rds" ) )

plot_distribuicao_parametro <- function( dados, parametro_alvo, ylab ) {
  
  dados_plot <- dados |>
    filter( parametro == parametro_alvo )
  
  ggplot( dados_plot, aes( x = nome_modelo, y = mediana, fill = nome_modelo ) ) +
    geom_boxplot(
      width = 0.58,
      linewidth = 0.3,
      outlier.shape = NA,
      alpha = 0.75
    ) +
    geom_jitter(
      width = 0.12,
      height = 0,
      size = 0.75,
      alpha = 0.45
    ) +
    scale_fill_manual(
      values = CORES_MODELOS,
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = ylab
    ) +
    tema_relatorio +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      )
    )
}

p_alpha <- plot_distribuicao_parametro(
  parametros_estruturais_todos,
  "alpha",
  expression(alpha)
)

p_beta <- plot_distribuicao_parametro(
  parametros_estruturais_todos,
  "beta",
  expression(beta)
)

p_H <- plot_distribuicao_parametro(
  parametros_estruturais_todos,
  "H",
  expression(H)
)

p_theta <- plot_distribuicao_parametro(
  parametros_estruturais_todos,
  "theta",
  expression(theta)
)

estacoes_sf <- meta_estacoes_modelo |>
  filter(
    gauge_code %in% gauge_codes
  ) |>
  sf::st_as_sf(
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )

rs_sf <- geobr::read_state(
  code_state = "RS",
  year = 2020,
  simplified = TRUE,
  showProgress = FALSE
)

rs_sf <- sf::st_transform(
  rs_sf,
  4326
)

parametros_mapa <- parametros_estruturais_todos |>
  left_join(
    meta_estacoes_modelo |>
      select(
        gauge_code,
        long,
        lat
      ),
    by = "gauge_code"
  )

NOMES_MODELOS <- c(
  "bhm_common" = expression(xi ~ "— Comum"),
  "bhm_xi_d"   = expression(xi[d] ~ "— Por duração"),
  "bhm_xi_j"   = expression(xi[j] ~ "— Por estação"),
  "bhm_xi_jd"  = expression(xi[j*d] ~ "— Estação–duração"),
  "nonpooled"  = "Local (Non-pooled)"
)

NOMES_MODELOS_TEXTO <- c(
  "bhm_common" = "ξ — Comum",
  "bhm_xi_d"   = "ξ_d — Por duração",
  "bhm_xi_j"   = "ξ_j — Por estação",
  "bhm_xi_jd"  = "ξ_j,d — Estação–duração",
  "nonpooled"  = "Local (Non-pooled)"
)

parametros_estruturais_todos <- parametros_estruturais_todos |>
  mutate(
    nome_modelo = factor(
      NOMES_MODELOS_TEXTO[modelo],
      levels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS]
    )
  )

plot_distribuicao_parametro_modelos <- function(
    dados,
    parametro_alvo,
    xlab
) {
  
  dados_plot <- dados |>
    filter(parametro == parametro_alvo)
  
  ggplot(
    dados_plot,
    aes(x = mediana)
  ) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 30,
      fill = "grey80",
      color = "grey30",
      linewidth = 0.25
    ) +
    geom_density(
      linewidth = 0.6
    ) +
    facet_wrap(
      ~ nome_modelo,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      x = xlab,
      y = "Densidade"
    ) +
    tema_relatorio +
    theme(
      panel.border = element_rect(
        linewidth = 0.4,
        fill = NA
      ),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
}

p_H <- plot_distribuicao_parametro_modelos(
  parametros_estruturais_todos,
  "H",
  "Mediana posterior de H"
)

salvar_figura(
  p_H,
  "distribuicao_H_modelos.png",
  largura = 15,
  altura = 11
)

p_alpha <- plot_distribuicao_parametro_modelos(
  parametros_estruturais_todos,
  "alpha",
  expression("Mediana posterior de " * alpha)
)

p_beta <- plot_distribuicao_parametro_modelos(
  parametros_estruturais_todos,
  "beta",
  expression("Mediana posterior de " * beta)
)

p_theta <- plot_distribuicao_parametro_modelos(
  parametros_estruturais_todos,
  "theta",
  expression("Mediana posterior de " * theta)
)

salvar_figura(
  p_theta,
  "distribuicao_theta_modelos.png",
  largura = 15,
  altura = 11
)

library(sf)
library(gstat)
library(ggspatial)
library(classInt)
library(grid)

criar_mapa_idw_parametro <- function(
    dados,
    modelo_alvo,
    parametro_alvo,
    breaks,
    nome_legenda
) {
  
  # ----------------------------------------------------------
  # 1. Dados do modelo/parâmetro
  # ----------------------------------------------------------
  
  dados_mapa <- dados |>
    filter(
      modelo == modelo_alvo,
      parametro == parametro_alvo
    ) |>
    left_join(
      meta_estacoes_modelo |>
        select(
          gauge_code,
          lat,
          long
        ),
      by = "gauge_code"
    ) |>
    filter(
      is.finite(lat),
      is.finite(long),
      is.finite(mediana)
    )
  
  # ----------------------------------------------------------
  # 2. Estações em SIRGAS 2000
  # ----------------------------------------------------------
  
  estacoes_sf <- st_as_sf(
    dados_mapa,
    coords = c("long", "lat"),
    crs = 4674,
    remove = FALSE
  )
  
  # ----------------------------------------------------------
  # 3. RS no mesmo CRS
  # ----------------------------------------------------------
  
  rs_mapa <- st_transform(
    rs_sf,
    4674
  )
  
  rs_union <- st_as_sf(
    st_union(
      st_make_valid(rs_mapa)
    )
  )
  
  st_crs(rs_union) <- st_crs(4674)
  
  # ----------------------------------------------------------
  # 4. Grade
  # ----------------------------------------------------------
  
  grade_pts <- st_make_grid(
    rs_union,
    cellsize = 0.01,
    what = "centers",
    square = TRUE
  ) |>
    st_as_sf()
  
  st_crs(grade_pts) <- st_crs(rs_union)
  
  grade_pts <- grade_pts[
    lengths(
      st_intersects(
        grade_pts,
        rs_union
      )
    ) > 0,
  ]
  
  stopifnot(
    st_crs(estacoes_sf) == st_crs(grade_pts)
  )
  
  # ----------------------------------------------------------
  # 5. IDW
  # ----------------------------------------------------------
  
  idw_res <- gstat::idw(
    mediana ~ 1,
    locations = estacoes_sf,
    newdata = grade_pts,
    idp = 2
  ) |>
    rename(
      valor_interp = var1.pred
    )
  
  # ----------------------------------------------------------
  # 6. Classes
  # ----------------------------------------------------------
  
  idw_res <- idw_res |>
    mutate(
      classe = cut(
        valor_interp,
        breaks = breaks,
        include.lowest = TRUE,
        dig.lab = 2
      )
    )
  
  coords <- st_coordinates(idw_res)
  
  idw_plot <- idw_res |>
    st_drop_geometry() |>
    mutate(
      x = coords[, 1],
      y = coords[, 2]
    )
  
  # ----------------------------------------------------------
  # 7. BBOX expandido
  # ----------------------------------------------------------
  
  bb <- st_bbox(rs_union)
  
  dx <- bb["xmax"] - bb["xmin"]
  dy <- bb["ymax"] - bb["ymin"]
  
  xlim_expand <- c(
    bb["xmin"] - 0.03 * dx,
    bb["xmax"] + 0.08 * dx
  )
  
  ylim_expand <- c(
    bb["ymin"] - 0.03 * dy,
    bb["ymax"] + 0.08 * dy
  )
  
  # ----------------------------------------------------------
  # 8. Mapa
  # ----------------------------------------------------------
  
  ggplot() +
    
    # -------------------------------------------------------------------------
  # Interpolação IDW
  # -------------------------------------------------------------------------
  
  geom_raster(
    data = idw_plot,
    aes(
      x = x,
      y = y,
      fill = classe
    )
  ) +
    
    # -------------------------------------------------------------------------
  # Limite do RS
  # -------------------------------------------------------------------------
  
  geom_sf(
    data = rs_mapa,
    fill = NA,
    color = "grey45",
    linewidth = 0.25
  ) +
    
    geom_sf(
      data = rs_union,
      fill = NA,
      color = "grey25",
      linewidth = 0.45
    ) +
    
    # -------------------------------------------------------------------------
  # Estações
  # -------------------------------------------------------------------------
  
  geom_point(
    data = dados_mapa,
    aes(
      x = long,
      y = lat
    ),
    color = "black",
    size = 0.5
  ) +
    
    # -------------------------------------------------------------------------
  # Classes
  # -------------------------------------------------------------------------
  
  scale_fill_brewer(
    palette = "GnBu",
    #palette = "YlOrRd",
    name = nome_legenda,
    drop = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Escala
  # -------------------------------------------------------------------------
  
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit(1.5, "mm"),
    text_family = font,
    text_cex = 0.6,
    bar_cols = c(
      "black",
      "aliceblue"
    ),
    line_width = 0.4
  ) +
    
    # -------------------------------------------------------------------------
  # Norte
  # -------------------------------------------------------------------------
  
  annotation_north_arrow(
    location = "tr",
    width = unit(1, "cm"),
    height = unit(1, "cm"),
    pad_x = unit(0.30, "cm"),
    pad_y = unit(0.30, "cm"),
    which_north = "true",
    style = north_arrow_nautical(
      text_family = font,
      text_size = fontsize,
      fill = c("black", "aliceblue")
    )
  ) +
    
    # -------------------------------------------------------------------------
  # Extensão espacial
  # -------------------------------------------------------------------------
  
  coord_sf(
    crs = st_crs(4674),
    xlim = xlim_expand,
    ylim = ylim_expand,
    expand = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Eixos
  # -------------------------------------------------------------------------
  
  labs(
    x = "Longitude",
    y = "Latitude",
    subtitle = NOMES_MODELOS_TEXTO[[modelo_alvo]]
  ) +
    
    tema_relatorio +
    
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      
      legend.key.width = unit(12, "pt"),
      legend.key.height = unit(8, "pt"),
      legend.key.spacing.x = unit(2, "pt"),
      legend.spacing.x = unit(1, "pt")
    )
}

criar_breaks_quantis <- function(
    dados,
    parametro_alvo,
    probs = c(0, 0.20, 0.40, 0.60, 0.80, 1)
) {
  
  valores <- dados |>
    filter(
      parametro == parametro_alvo
    ) |>
    pull(mediana)
  
  breaks <- quantile(
    valores,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  )
  
  breaks <- round(breaks, 2)
  breaks <- unique(breaks)
  
  # garantir que extremos incluam todos os valores
  breaks[1] <- floor(min(valores, na.rm = TRUE) * 100) / 100
  breaks[length(breaks)] <-
    ceiling(max(valores, na.rm = TRUE) * 100) / 100
  
  breaks
}

breaks_H <- criar_breaks_quantis(
  parametros_estruturais_todos,
  "H"
)

breaks_H

mapas_H <- purrr::map(
  ORDEM_MODELOS,
  ~ criar_mapa_idw_parametro(
    dados = parametros_estruturais_todos,
    modelo_alvo = .x,
    parametro_alvo = "H",
    breaks = breaks_H,
    nome_legenda = expression(H)
  )
)

fig_mapa_H <- wrap_plots(
  mapas_H,
  ncol = 3,
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

salvar_figura(
  fig_mapa_H,
  "mapa_H_modelos.png",
  largura = 21,
  altura = 16
)

breaks_alpha <- criar_breaks_quantis(
  parametros_estruturais_todos,
  "alpha"
)

breaks_beta <- criar_breaks_quantis(
  parametros_estruturais_todos,
  "beta"
)

breaks_theta <- criar_breaks_quantis(
  parametros_estruturais_todos,
  "theta"
)

mapas_theta <- purrr::map(
  ORDEM_MODELOS,
  ~ criar_mapa_idw_parametro(
    dados = parametros_estruturais_todos,
    modelo_alvo = .x,
    parametro_alvo = "theta",
    breaks = breaks_theta,
    nome_legenda = expression(theta)
  )
)

fig_mapa_theta <- wrap_plots(
  mapas_theta,
  ncol = 3,
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

salvar_figura(
  fig_mapa_theta,
  "mapa_theta_modelos.png",
  largura = 21,
  altura = 16
)

# =============================================================================
# 11. PARÂMETROS GEV
# =============================================================================

D_GRID <- 1:24
D_REF <- 24
K <- length(D_GRID)

localizar_parametro <- function(fit, candidatos) {
  
  disponiveis <- fit@sim$pars_oi
  
  encontrado <- candidatos[
    candidatos %in% disponiveis
  ]
  
  if (length(encontrado) == 0L) {
    return(NULL)
  }
  
  encontrado[1]
}

extrair_mu_sigma_regionais <- function(fit) {
  
  J <- obter_J(fit)
  K <- length(D_GRID)
  
  # ----------------------------------------------------------
  # Tentar encontrar mu e sigma diretamente no Stan
  # ----------------------------------------------------------
  
  nome_mu <- localizar_parametro(
    fit,
    c(
      "mu",
      "mu_jd",
      "mu_jk",
      "mu_est"
    )
  )
  
  nome_sigma <- localizar_parametro(
    fit,
    c(
      "sigma",
      "sigma_jd",
      "sigma_jk",
      "sigma_est"
    )
  )
  
  mu <- NULL
  sigma <- NULL
  
  if (!is.null(nome_mu)) {
    
    x <- rstan::extract(
      fit,
      pars = nome_mu,
      permuted = TRUE
    )[[nome_mu]]
    
    mu <- tryCatch(
      normalizar_jk(x, J, K),
      error = function(e) NULL
    )
  }
  
  if (!is.null(nome_sigma)) {
    
    x <- rstan::extract(
      fit,
      pars = nome_sigma,
      permuted = TRUE
    )[[nome_sigma]]
    
    sigma <- tryCatch(
      normalizar_jk(x, J, K),
      error = function(e) NULL
    )
  }
  
  # ----------------------------------------------------------
  # Se não estiverem salvos:
  # reconstruir pela estrutura duração–intensidade
  # ----------------------------------------------------------
  
  if (is.null(mu) || is.null(sigma)) {
    
    p <- extrair_parametros_estruturais(fit)
    
    alpha <- p$alpha
    beta <- p$beta
    H <- p$H
    theta <- p$theta
    
    S <- nrow(alpha)
    
    mu_rec <- array(
      NA_real_,
      dim = c(S, J, K)
    )
    
    sigma_rec <- array(
      NA_real_,
      dim = c(S, J, K)
    )
    
    for (k in seq_len(K)) {
      
      d <- D_GRID[k]
      
      fator <- (
        (d + theta) /
          (D_REF + theta)
      )^(-H)
      
      mu_rec[, , k] <-
        exp(alpha) * fator
      
      sigma_rec[, , k] <-
        exp(beta) * fator
    }
    
    if (is.null(mu)) {
      mu <- mu_rec
    }
    
    if (is.null(sigma)) {
      sigma <- sigma_rec
    }
  }
  
  list(
    mu = mu,
    sigma = sigma
  )
}

extrair_xi_regionais <- function(
    fit,
    modelo
) {
  
  J <- obter_J(fit)
  K <- length(D_GRID)
  
  candidatos <- switch(
    modelo,
    
    "bhm_common" = c(
      "xi",
      "xi_bar",
      "xi_common"
    ),
    
    "bhm_xi_d" = c(
      "xi_d",
      "xi_k",
      "xi"
    ),
    
    "bhm_xi_j" = c(
      "xi_j",
      "xi",
      "xi_station"
    ),
    
    "bhm_xi_jd" = c(
      "xi_jd",
      "xi_jk",
      "xi"
    )
  )
  
  nome <- localizar_parametro(
    fit,
    candidatos
  )
  
  if (is.null(nome)) {
    
    stop(
      "Não encontrei xi no modelo ",
      modelo,
      ". Parâmetros com 'xi': ",
      paste(
        grep(
          "xi",
          fit@sim$pars_oi,
          value = TRUE
        ),
        collapse = ", "
      )
    )
  }
  
  x <- rstan::extract(
    fit,
    pars = nome,
    permuted = TRUE
  )[[nome]]
  
  dx <- dim(x)
  
  # ==========================================================
  # 1. xi comum
  # draw
  # ==========================================================
  
  if (
    is.null(dx) ||
    length(dx) == 1
  ) {
    
    x <- as.numeric(x)
    
    S <- length(x)
    
    out <- array(
      NA_real_,
      dim = c(S, J, K)
    )
    
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        out[, j, k] <- x
      }
    }
    
    return(out)
  }
  
  # ==========================================================
  # 2. xi_j ou xi_d
  # draw × dimensão
  # ==========================================================
  
  if (length(dx) == 2) {
    
    S <- dx[1]
    
    # --------------------------------------------------------
    # xi por estação
    # --------------------------------------------------------
    
    if (
      modelo == "bhm_xi_j" &&
      dx[2] == J
    ) {
      
      out <- array(
        NA_real_,
        dim = c(S, J, K)
      )
      
      for (k in seq_len(K)) {
        out[, , k] <- x
      }
      
      return(out)
    }
    
    # --------------------------------------------------------
    # xi por duração
    # --------------------------------------------------------
    
    if (
      modelo == "bhm_xi_d" &&
      dx[2] == K
    ) {
      
      out <- array(
        NA_real_,
        dim = c(S, J, K)
      )
      
      for (j in seq_len(J)) {
        out[, j, ] <- x
      }
      
      return(out)
    }
  }
  
  # ==========================================================
  # 3. xi estação × duração
  # ==========================================================

  if (
    modelo == "bhm_xi_jd" &&
    length(dx) == 3
  ) {

    return(
      normalizar_jk(
        x,
        J,
        K
      )
    )
  }

  stop(
    "Não consegui interpretar ",
    nome,
    " no modelo ",
    modelo,
    ". Dimensões = ",
    paste(dx, collapse = " x ")
  )
}

resumir_gev_regional <- function(
    fit,
    modelo
) {
  
  ms <- extrair_mu_sigma_regionais(fit)
  
  mu <- ms$mu
  sigma <- ms$sigma
  
  xi <- extrair_xi_regionais(
    fit,
    modelo
  )
  
  J <- dim(mu)[2]
  K <- dim(mu)[3]
  
  resultados <- vector(
    "list",
    J * K
  )
  
  contador <- 1L
  
  for (j in seq_len(J)) {
    
    for (k in seq_len(K)) {
      
      mu_x <- mu[, j, k]
      sigma_x <- sigma[, j, k]
      xi_x <- xi[, j, k]
      
      resultados[[contador]] <- tibble(
        
        modelo = modelo,
        station_id = j,
        duracao_h = D_GRID[k],
        
        mu_mediana =
          median(mu_x, na.rm = TRUE),
        
        mu_q025 =
          quantile(mu_x, 0.025, na.rm = TRUE),
        
        mu_q975 =
          quantile(mu_x, 0.975, na.rm = TRUE),
        
        sigma_mediana =
          median(sigma_x, na.rm = TRUE),
        
        sigma_q025 =
          quantile(sigma_x, 0.025, na.rm = TRUE),
        
        sigma_q975 =
          quantile(sigma_x, 0.975, na.rm = TRUE),
        
        xi_mediana =
          median(xi_x, na.rm = TRUE),
        
        xi_q025 =
          quantile(xi_x, 0.025, na.rm = TRUE),
        
        xi_q975 =
          quantile(xi_x, 0.975, na.rm = TRUE)
      )
      
      contador <- contador + 1L
    }
  }
  
  bind_rows(resultados)
}

parametros_gev_regionais <- purrr::imap_dfr(
  fits_regionais,
  resumir_gev_regional
)

parametros_gev_regionais <-
  parametros_gev_regionais |>
  left_join(
    mapa_station_id,
    by = "station_id"
  ) |>
  relocate(
    gauge_code,
    .after = station_id
  )

extrair_gev_nonpooled_fit <- function(
    fit,
    gauge_code
) {
  
  alpha <- as.numeric(
    extrair_primeiro_parametro(fit, "alpha")
  )
  
  beta <- as.numeric(
    extrair_primeiro_parametro(fit, "beta")
  )
  
  H <- as.numeric(
    extrair_primeiro_parametro(fit, "H")
  )
  
  theta <- as.numeric(
    extrair_primeiro_parametro(fit, "theta")
  )
  
  nome_xi <- localizar_parametro(
    fit,
    c(
      "xi",
      "xi_d",
      "xi_k",
      "xi_bar",
      "xi_j"
    )
  )
  
  xi <- rstan::extract(
    fit,
    pars = nome_xi,
    permuted = TRUE
  )[[nome_xi]]
  
  # garantir draw × duração
  if (is.null(dim(xi))) {
    
    xi <- matrix(
      xi,
      nrow = length(xi),
      ncol = length(D_GRID)
    )
    
  } else if (
    nrow(xi) == length(D_GRID) &&
    ncol(xi) == length(alpha)
  ) {
    
    xi <- t(xi)
  }
  
  purrr::map_dfr(
    seq_along(D_GRID),
    function(k) {
      
      d <- D_GRID[k]
      
      fator <- (
        (d + theta) /
          (D_REF + theta)
      )^(-H)
      
      mu <- exp(alpha) * fator
      sigma <- exp(beta) * fator
      
      xi_d <- as.numeric(
        xi[, k]
      )
      
      tibble(
        modelo = "nonpooled",
        gauge_code = gauge_code,
        duracao_h = d,
        
        mu_mediana = median(mu, na.rm = TRUE),
        mu_q025 = quantile(mu, 0.025, na.rm = TRUE),
        mu_q975 = quantile(mu, 0.975, na.rm = TRUE),
        
        sigma_mediana = median(sigma, na.rm = TRUE),
        sigma_q025 = quantile(sigma, 0.025, na.rm = TRUE),
        sigma_q975 = quantile(sigma, 0.975, na.rm = TRUE),
        
        xi_mediana = median(xi_d, na.rm = TRUE),
        xi_q025 = quantile(xi_d, 0.025, na.rm = TRUE),
        xi_q975 = quantile(xi_d, 0.975, na.rm = TRUE)
      )
    }
  )
}

parametros_gev_nonpooled <- purrr::imap_dfr(
  fits_nonpooled,
  extrair_gev_nonpooled_fit
)

parametros_gev_nonpooled <- purrr::imap_dfr(
  fits_nonpooled,
  extrair_gev_nonpooled_fit
) |>
  left_join(
    mapa_station_id,
    by = "gauge_code"
  ) |>
  relocate(
    station_id,
    .after = modelo
  )

parametros_gev_todos <- bind_rows(
  parametros_gev_regionais,
  parametros_gev_nonpooled
) |>
  mutate(
    nome_modelo = factor(
      NOMES_MODELOS_TEXTO[modelo],
      levels =
        NOMES_MODELOS_TEXTO[
          ORDEM_MODELOS
        ]
    )
  )

parametros_gev_todos |>
  group_by(modelo) |>
  summarise(
    mu_min = min(mu_mediana, na.rm = TRUE),
    mu_max = max(mu_mediana, na.rm = TRUE),
    
    sigma_min = min(sigma_mediana, na.rm = TRUE),
    sigma_max = max(sigma_mediana, na.rm = TRUE),
    
    xi_min = min(xi_mediana, na.rm = TRUE),
    xi_max = max(xi_mediana, na.rm = TRUE),
    
    .groups = "drop"
  )

resumo_gev_duracao <- parametros_gev_todos |>
  group_by(
    modelo,
    nome_modelo,
    duracao_h
  ) |>
  summarise(
    
    mu_med = median(
      mu_mediana,
      na.rm = TRUE
    ),
    
    mu_q05 = quantile(
      mu_mediana,
      0.05,
      na.rm = TRUE
    ),
    
    mu_q95 = quantile(
      mu_mediana,
      0.95,
      na.rm = TRUE
    ),
    
    sigma_med = median(
      sigma_mediana,
      na.rm = TRUE
    ),
    
    sigma_q05 = quantile(
      sigma_mediana,
      0.05,
      na.rm = TRUE
    ),
    
    sigma_q95 = quantile(
      sigma_mediana,
      0.95,
      na.rm = TRUE
    ),
    
    xi_med = median(
      xi_mediana,
      na.rm = TRUE
    ),
    
    xi_q05 = quantile(
      xi_mediana,
      0.05,
      na.rm = TRUE
    ),
    
    xi_q95 = quantile(
      xi_mediana,
      0.95,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

CORES_MODELOS_ID <- c(
  "bhm_common" = "#D55E00",
  "bhm_xi_d"   = "#E69F00",
  "bhm_xi_j"   = "#009E73",
  "bhm_xi_jd"  = "#0072B2",
  "nonpooled"  = "#CC79A7"
)
dados_A801 <- parametros_gev_todos |>
  filter(gauge_code == "A801")

p_mu_A801 <- ggplot(
  dados_A801,
  aes(
    x = duracao_h,
    y = mu_mediana,
    color = modelo,
    fill = modelo,
    group = modelo
  )
) +
  
  geom_ribbon(
    aes(
      ymin = mu_q025,
      ymax = mu_q975
    ),
    alpha = 0.12,
    color = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  scale_color_manual(
    values = CORES_MODELOS_ID,
    breaks = ORDEM_MODELOS,
    labels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS],
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = CORES_MODELOS_ID,
    breaks = ORDEM_MODELOS,
    labels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS],
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(1, 3, 6, 12, 18, 24)
  ) +
  
  labs(
    x = "Duração (h)",
    y = expression(mu),
    color = NULL,
    fill = NULL
  ) +
  
  tema_relatorio +
  
  theme(
    legend.position = "bottom"
  )

p_sigma_A801 <- ggplot(
  dados_A801,
  aes(
    x = duracao_h,
    y = sigma_mediana,
    color = modelo,
    fill = modelo,
    group = modelo
  )
) +
  
  geom_ribbon(
    aes(
      ymin = sigma_q025,
      ymax = sigma_q975
    ),
    alpha = 0.12,
    color = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  scale_color_manual(
    values = CORES_MODELOS_ID,
    breaks = ORDEM_MODELOS,
    labels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS],
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = CORES_MODELOS_ID,
    breaks = ORDEM_MODELOS,
    labels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS],
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(1, 3, 6, 12, 18, 24)
  ) +
  
  labs(
    x = "Duração (h)",
    y = expression(sigma),
    color = NULL,
    fill = NULL
  ) +
  
  tema_relatorio +
  
  theme(
    legend.position = "bottom"
  )

p_xi_A801 <- ggplot(
  dados_A801,
  aes(
    x = duracao_h,
    y = xi_mediana,
    color = modelo,
    fill = modelo,
    group = modelo
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  
  geom_ribbon(
    aes(
      ymin = xi_q025,
      ymax = xi_q975
    ),
    alpha = 0.12,
    color = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  scale_color_manual(
    values = CORES_MODELOS_ID,
    breaks = ORDEM_MODELOS,
    labels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS],
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = CORES_MODELOS_ID,
    breaks = ORDEM_MODELOS,
    labels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS],
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(1, 3, 6, 12, 18, 24)
  ) +
  
  labs(
    x = "Duração (h)",
    y = expression(xi),
    color = NULL,
    fill = NULL
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE
    ),
    fill = "none"
  ) +
  
  tema_relatorio +
  
  theme(
    legend.position = "bottom"
  )

p_xi_A801

salvar_figura(
  p_xi_A801,
  "parametro_GEV_xi_A801.png",
  largura = 17,
  altura = 8
)

diag_regionais_gpkg <- flags_diagnostico_regionais |>
  transmute(
    modelo,
    station_id,
    flag_diagnostico =
      flag_rhat_105 |
      flag_divergencia_modelo
  )

diag_nonpooled_gpkg <- flags_diagnostico_nonpooled |>
  left_join(
    mapa_station_id,
    by = "gauge_code"
  ) |>
  transmute(
    modelo = "nonpooled",
    station_id,
    flag_diagnostico
  )

diagnosticos_gpkg <- bind_rows(
  diag_regionais_gpkg,
  diag_nonpooled_gpkg
)

parametros_gev_gpkg <- parametros_gev_todos |>
  left_join(
    diagnosticos_gpkg,
    by = c(
      "modelo",
      "station_id"
    )
  )

parametros_gev_gpkg <- parametros_gev_gpkg |>
  left_join(
    meta_estacoes_modelo |>
      select(
        gauge_code,
        lat,
        long,
        elevation,
        network,
        responsible
      ),
    by = "gauge_code"
  )


parametros_gev_sf <- parametros_gev_gpkg |>
  sf::st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4674,
    remove = FALSE
  )

estruturais_gpkg <- parametros_estruturais_todos |>
  left_join(
    diagnosticos_gpkg,
    by = c(
      "modelo",
      "station_id"
    )
  ) |>
  left_join(
    meta_estacoes_modelo |>
      select(
        gauge_code,
        lat,
        long,
        elevation,
        network,
        responsible
      ),
    by = "gauge_code"
  )

estruturais_sf <- estruturais_gpkg |>
  sf::st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4674,
    remove = FALSE
  )

arquivo_gpkg_parametros <- file.path(
  DIR_GPKG,
  "RS_parametros_modelos.gpkg"
)

if (file.exists(arquivo_gpkg_parametros)) {
  file.remove(arquivo_gpkg_parametros)
}

sf::st_write(
  estruturais_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "parametros_estruturais",
  quiet = TRUE
)

sf::st_write(
  parametros_gev_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "parametros_gev",
  append = FALSE,
  quiet = TRUE
)

salvar_tabela(
  parametros_gev_gpkg,
  "parametros_GEV_cinco_modelos.csv"
)

saveRDS(
  parametros_gev_gpkg,
  file.path(
    DIR_RESULTADOS,
    "parametros_GEV_cinco_modelos.rds"
  )
)

dados_xi_24h <- parametros_gev_todos |>
  filter(duracao_h == 24)

breaks_xi_24h <- c(
  -0.45,
  -0.20,
  -0.05,
  0.05,
  0.20,
  0.40
)

rotulos_breaks <- function(breaks) {
  n <- length(breaks) - 1
  
  labs <- character(n)
  
  labs[1] <- paste0(
    "[",
    format(round(breaks[1], 2), nsmall = 2),
    ",",
    format(round(breaks[2], 2), nsmall = 2),
    "]"
  )
  
  if (n >= 2) {
    for (i in 2:n) {
      labs[i] <- paste0(
        "(",
        format(round(breaks[i], 2), nsmall = 2),
        ",",
        format(round(breaks[i + 1], 2), nsmall = 2),
        "]"
      )
    }
  }
  
  labs
}

criar_mapa_idw_xi_24h <- function(
    dados,
    modelo_alvo,
    breaks = breaks_xi_24h,
    mostrar_legenda = TRUE
) {
  
  labels_breaks <- rotulos_breaks(breaks)
  
  dados_mapa <- dados |>
    filter(
      modelo == modelo_alvo,
      duracao_h == 24
    ) |>
    left_join(
      meta_estacoes_modelo |>
        select(
          gauge_code,
          lat,
          long
        ),
      by = "gauge_code"
    ) |>
    filter(
      is.finite(lat),
      is.finite(long),
      is.finite(xi_mediana)
    )
  
  estacoes_sf <- st_as_sf(
    dados_mapa,
    coords = c("long", "lat"),
    crs = 4674,
    remove = FALSE
  )
  
  rs_mapa <- st_transform(
    rs_sf,
    4674
  )
  
  rs_union <- st_as_sf(
    st_union(
      st_make_valid(rs_mapa)
    )
  )
  
  st_crs(rs_union) <- st_crs(4674)
  
  grade_pts <- st_make_grid(
    rs_union,
    cellsize = 0.01,
    what = "centers",
    square = TRUE
  ) |>
    st_as_sf()
  
  st_crs(grade_pts) <- st_crs(rs_union)
  
  grade_pts <- grade_pts[
    lengths(
      st_intersects(
        grade_pts,
        rs_union
      )
    ) > 0,
  ]
  
  idw_res <- gstat::idw(
    xi_mediana ~ 1,
    locations = estacoes_sf,
    newdata = grade_pts,
    idp = 2
  ) |>
    rename(
      valor_interp = var1.pred
    ) |>
    mutate(
      classe = cut(
        valor_interp,
        breaks = breaks,
        labels = labels_breaks,
        include.lowest = TRUE,
        right = TRUE
      ),
      classe = factor(
        classe,
        levels = labels_breaks
      )
    )
  
  coords <- st_coordinates(idw_res)
  
  idw_plot <- idw_res |>
    st_drop_geometry() |>
    mutate(
      x = coords[, 1],
      y = coords[, 2]
    )
  
  bb <- st_bbox(rs_union)
  
  dx <- bb["xmax"] - bb["xmin"]
  dy <- bb["ymax"] - bb["ymin"]
  
  xlim_expand <- c(
    bb["xmin"] - 0.03 * dx,
    bb["xmax"] + 0.12 * dx
  )
  
  ylim_expand <- c(
    bb["ymin"] - 0.03 * dy,
    bb["ymax"] + 0.12 * dy
  )
  
  p <- ggplot() +
    
    # -------------------------------------------------------------------------
  # Interpolação IDW
  # -------------------------------------------------------------------------
  
  geom_raster(
    data = idw_plot,
    aes(
      x = x,
      y = y,
      fill = classe
    )
  ) +
    
    # -------------------------------------------------------------------------
  # Limite do RS
  # -------------------------------------------------------------------------
  
  geom_sf(
    data = rs_mapa,
    fill = NA,
    color = "grey45",
    linewidth = 0.25
  ) +
    
    geom_sf(
      data = rs_union,
      fill = NA,
      color = "grey25",
      linewidth = 0.45
    ) +
    
    # -------------------------------------------------------------------------
  # Estações
  # -------------------------------------------------------------------------
  
  geom_point(
    data = dados_mapa,
    aes(
      x = long,
      y = lat
    ),
    color = "black",
    size = 0.7
  ) +
    
    # -------------------------------------------------------------------------
  # Classes de xi
  # -------------------------------------------------------------------------
  
  scale_fill_manual(
    values = c(
      "#2166AC",
      "#67A9CF",
      "#F7F7F7",
      "#FDDBC7",
      "#B2182B"
    ),
    breaks = labels_breaks,
    limits = labels_breaks,
    name = expression(xi ~ "(24 h)"),
    drop = FALSE
  ) +
    
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "left"
      )
    ) +
    
    # -------------------------------------------------------------------------
  # Escala
  # -------------------------------------------------------------------------
  
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit(1.5, "mm"),
    text_family = font,
    text_cex = 0.6,
    bar_cols = c(
      "black",
      "aliceblue"
    ),
    line_width = 0.4
  ) +
    
    # -------------------------------------------------------------------------
  # Norte
  # -------------------------------------------------------------------------
  
  annotation_north_arrow(
    location = "tr",
    width = unit(1, "cm"),
    height = unit(1, "cm"),
    pad_x = unit(0.30, "cm"),
    pad_y = unit(0.30, "cm"),
    which_north = "true",
    style = north_arrow_nautical(
      text_family = font,
      text_size = fontsize,
      fill = c(
        "black",
        "aliceblue"
      )
    )
  ) +
    
    # -------------------------------------------------------------------------
  # Extensão espacial
  # -------------------------------------------------------------------------
  
  coord_sf(
    crs = st_crs(4674),
    xlim = xlim_expand,
    ylim = ylim_expand,
    expand = FALSE
  ) +
    
    # -------------------------------------------------------------------------
  # Eixos
  # -------------------------------------------------------------------------
  
  labs(
    x = "Longitude",
    y = "Latitude",
    subtitle = NOMES_MODELOS_TEXTO[[modelo_alvo]]
  ) +
    
    tema_relatorio +
    
    theme(
      legend.position = if (mostrar_legenda) "bottom" else "none",
      legend.justification = "center",
      legend.direction = "horizontal",
      
      legend.key.width = unit(12, "pt"),
      legend.key.height = unit(8, "pt"),
      legend.key.spacing.x = unit(2, "pt"),
      legend.spacing.x = unit(1, "pt")
    )
  
  p
}

# recria os mapas normalmente
mapas_xi_24h <- purrr::map(
  ORDEM_MODELOS,
  ~ criar_mapa_idw_xi_24h(
    dados = dados_xi_24h,
    modelo_alvo = .x,
    breaks = breaks_xi_24h,
    mostrar_legenda = TRUE
  )
)

library(cowplot)

legenda_xi <- cowplot::get_legend(
  mapas_xi_24h[[3]] +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

mapas_xi_24h_sem_legenda <- purrr::map(
  mapas_xi_24h,
  ~ .x +
    theme(
      legend.position = "none"
    )
)

painel_xi <- patchwork::wrap_plots(
  mapas_xi_24h_sem_legenda,
  ncol = 3
)

fig_mapa_xi_24h <- cowplot::plot_grid(
  painel_xi,
  legenda_xi,
  ncol = 1,
  rel_heights = c(1, 0.07)
)

fig_mapa_xi_24h

salvar_figura(
  fig_mapa_xi_24h,
  "mapa_xi_24h_modelos.png",
  largura = 21,
  altura = 16
)

# =============================================================================
# 12. MAPAS DE MU E SIGMA PARA 24 H
# =============================================================================

dados_gev_24h <- parametros_gev_todos |>
  filter(duracao_h == 24)

criar_breaks_gev_24h <- function(
    dados,
    variavel,
    probs = c(0, 0.20, 0.40, 0.60, 0.80, 1)
) {
  
  valores <- dados[[variavel]]
  
  valores <- valores[
    is.finite(valores)
  ]
  
  breaks <- quantile(
    valores,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  )
  
  breaks <- round(
    breaks,
    2
  )
  
  breaks <- unique(
    breaks
  )
  
  # garantir inclusão dos extremos
  breaks[1] <-
    floor(
      min(valores, na.rm = TRUE) * 100
    ) / 100
  
  breaks[length(breaks)] <-
    ceiling(
      max(valores, na.rm = TRUE) * 100
    ) / 100
  
  breaks
}

breaks_mu_24h <- criar_breaks_gev_24h(
  dados_gev_24h,
  "mu_mediana"
)

breaks_sigma_24h <- criar_breaks_gev_24h(
  dados_gev_24h,
  "sigma_mediana"
)

breaks_mu_24h
breaks_sigma_24h

criar_mapa_idw_gev_24h <- function(
    dados,
    modelo_alvo,
    variavel,
    breaks,
    nome_legenda,
    mostrar_legenda = TRUE
) {
  
  # ----------------------------------------------------------
  # 1. Rótulos das classes
  # ----------------------------------------------------------
  
  labels_breaks <- rotulos_breaks(breaks)
  
  # ----------------------------------------------------------
  # 2. Dados
  # ----------------------------------------------------------
  
  dados_mapa <- dados |>
    filter(
      modelo == modelo_alvo,
      duracao_h == 24
    ) |>
    left_join(
      meta_estacoes_modelo |>
        select(
          gauge_code,
          lat,
          long
        ),
      by = "gauge_code"
    ) |>
    filter(
      is.finite(lat),
      is.finite(long),
      is.finite(.data[[variavel]])
    )
  
  # ----------------------------------------------------------
  # 3. Estações
  # ----------------------------------------------------------
  
  estacoes_sf <- st_as_sf(
    dados_mapa,
    coords = c("long", "lat"),
    crs = 4674,
    remove = FALSE
  )
  
  # ----------------------------------------------------------
  # 4. RS
  # ----------------------------------------------------------
  
  rs_mapa <- st_transform(
    rs_sf,
    4674
  )
  
  rs_union <- st_as_sf(
    st_union(
      st_make_valid(rs_mapa)
    )
  )
  
  st_crs(rs_union) <- st_crs(4674)
  
  # ----------------------------------------------------------
  # 5. Grade
  # ----------------------------------------------------------
  
  grade_pts <- st_make_grid(
    rs_union,
    cellsize = 0.01,
    what = "centers",
    square = TRUE
  ) |>
    st_as_sf()
  
  st_crs(grade_pts) <- st_crs(rs_union)
  
  grade_pts <- grade_pts[
    lengths(
      st_intersects(
        grade_pts,
        rs_union
      )
    ) > 0,
  ]
  
  stopifnot(
    st_crs(estacoes_sf) == st_crs(grade_pts)
  )
  
  # ----------------------------------------------------------
  # 6. IDW
  # ----------------------------------------------------------
  
  formula_idw <- as.formula(
    paste0(
      variavel,
      " ~ 1"
    )
  )
  
  idw_res <- gstat::idw(
    formula = formula_idw,
    locations = estacoes_sf,
    newdata = grade_pts,
    idp = 2
  ) |>
    rename(
      valor_interp = var1.pred
    )
  
  # ----------------------------------------------------------
  # 7. Classes fixas
  # ----------------------------------------------------------
  
  idw_res <- idw_res |>
    mutate(
      classe = cut(
        valor_interp,
        breaks = breaks,
        labels = labels_breaks,
        include.lowest = TRUE,
        right = TRUE
      ),
      classe = factor(
        classe,
        levels = labels_breaks
      )
    )
  
  coords <- st_coordinates(idw_res)
  
  idw_plot <- idw_res |>
    st_drop_geometry() |>
    mutate(
      x = coords[, 1],
      y = coords[, 2]
    )
  
  # ----------------------------------------------------------
  # 8. BBOX expandido
  # ----------------------------------------------------------
  
  bb <- st_bbox(rs_union)
  
  dx <- bb["xmax"] - bb["xmin"]
  dy <- bb["ymax"] - bb["ymin"]
  
  xlim_expand <- c(
    bb["xmin"] - 0.03 * dx,
    bb["xmax"] + 0.12 * dx
  )
  
  ylim_expand <- c(
    bb["ymin"] - 0.03 * dy,
    bb["ymax"] + 0.12 * dy
  )
  
  # ----------------------------------------------------------
  # 9. Mapa
  # ----------------------------------------------------------
  
  ggplot() +
    
    geom_raster(
      data = idw_plot,
      aes(
        x = x,
        y = y,
        fill = classe
      )
    ) +
    
    geom_point(
      data = dados_mapa,
      aes(
        x = long,
        y = lat
      ),
      color = "black",
      size = 0.7
    ) +
    
    geom_sf(
      data = rs_mapa,
      fill = NA,
      color = "grey45",
      linewidth = 0.20
    ) +
    
    geom_sf(
      data = rs_union,
      fill = NA,
      color = "grey30",
      linewidth = 0.40
    ) +
    
    scale_fill_brewer(
      palette = "YlOrRd",
      breaks = labels_breaks,
      limits = labels_breaks,
      name = nome_legenda,
      drop = FALSE
    ) +
    
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "left"
      )
    ) +
    
    annotation_scale(
      location = "br",
      width_hint = 0.20,
      height = unit(1.5, "mm"),
      text_family = font,
      text_cex = 0.6,
      bar_cols = c(
        "black",
        "aliceblue"
      ),
      line_width = 0.4
    ) +
    
    annotation_north_arrow(
      location = "tr",
      width = unit(1, "cm"),
      height = unit(1, "cm"),
      which_north = "true",
      style = north_arrow_nautical(
        text_family = font,
        text_size = fontsize,
        fill = c(
          "black",
          "aliceblue"
        )
      )
    ) +
    
    coord_sf(
      crs = st_crs(4674),
      xlim = xlim_expand,
      ylim = ylim_expand,
      expand = FALSE
    ) +
    
    labs(
      x = NULL,
      y = NULL,
      subtitle = NOMES_MODELOS_TEXTO[[modelo_alvo]]
    ) +
    
    tema_relatorio +
    
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position =
        if (mostrar_legenda) "bottom" else "none",
      legend.justification = "center"
    )
}

mapas_mu_24h <- purrr::map(
  ORDEM_MODELOS,
  ~ criar_mapa_idw_gev_24h(
    dados = dados_gev_24h,
    modelo_alvo = .x,
    variavel = "mu_mediana",
    breaks = breaks_mu_24h,
    nome_legenda = expression(mu*" (24 h)"),
    mostrar_legenda = TRUE
  )
)

legenda_mu <- cowplot::get_legend(
  mapas_mu_24h[[3]] +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

mapas_mu_24h_sem_legenda <- purrr::map(
  mapas_mu_24h,
  ~ .x +
    theme(
      legend.position = "none"
    )
)

painel_mu <- patchwork::wrap_plots(
  mapas_mu_24h_sem_legenda,
  ncol = 3
)

fig_mapa_mu_24h <- cowplot::plot_grid(
  painel_mu,
  legenda_mu,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

fig_mapa_mu_24h

salvar_figura(
  fig_mapa_mu_24h,
  "mapa_mu_24h_modelos.png",
  largura = 18,
  altura = 14
)

mapas_sigma_24h <- purrr::map(
  ORDEM_MODELOS,
  ~ criar_mapa_idw_gev_24h(
    dados = dados_gev_24h,
    modelo_alvo = .x,
    variavel = "sigma_mediana",
    breaks = breaks_sigma_24h,
    nome_legenda = expression(sigma*" (24 h)"),
    mostrar_legenda = TRUE
  )
)

legenda_sigma <- cowplot::get_legend(
  mapas_sigma_24h[[3]] +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

mapas_sigma_24h_sem_legenda <- purrr::map(
  mapas_sigma_24h,
  ~ .x +
    theme(
      legend.position = "none"
    )
)

painel_sigma <- patchwork::wrap_plots(
  mapas_sigma_24h_sem_legenda,
  ncol = 3
)

fig_mapa_sigma_24h <- cowplot::plot_grid(
  painel_sigma,
  legenda_sigma,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

fig_mapa_sigma_24h

salvar_figura(
  fig_mapa_sigma_24h,
  "mapa_sigma_24h_modelos.png",
  largura = 18,
  altura = 14
)

# =============================================================================
# 12. NÍVEIS DE RETORNO
# =============================================================================

TR_GRID <- c(2, 5, 10, 25, 50, 100, 200)
modelo_idf_figura <- "bhm_xi_j"

qgev_retorno <- function(mu, sigma, xi, T) {
  
  p <- 1 - 1 / T
  z <- -log(p)
  
  sigma <- pmax(sigma, .Machine$double.eps)
  
  out <- ifelse(
    abs(xi) < 1e-8,
    mu - sigma * log(z),
    mu + (sigma / xi) * (z^(-xi) - 1)
  )
  
  out
}

calcular_niveis_retorno_regional <- function(fit, modelo) {
  
  ms <- extrair_mu_sigma_regionais(fit)
  xi <- extrair_xi_regionais(fit, modelo)
  
  mu <- ms$mu
  sigma <- ms$sigma
  
  J <- dim(mu)[2]
  K <- dim(mu)[3]
  
  resultados <- vector("list", J * K * length(TR_GRID))
  contador <- 1L
  
  for (j in seq_len(J)) {
    for (k in seq_len(K)) {
      for (Tret in TR_GRID) {
        
        rl <- qgev_retorno(
          mu = mu[, j, k],
          sigma = sigma[, j, k],
          xi = xi[, j, k],
          T = Tret
        )
        
        resultados[[contador]] <- tibble(
          modelo = modelo,
          station_id = j,
          duracao_h = D_GRID[k],
          return_period_years = Tret,
          median = median(rl, na.rm = TRUE),
          q025 = quantile(rl, 0.025, na.rm = TRUE),
          q975 = quantile(rl, 0.975, na.rm = TRUE)
        )
        
        contador <- contador + 1L
      }
    }
  }
  
  bind_rows(resultados)
}

niveis_retorno_regionais <- purrr::imap_dfr(
  fits_regionais,
  calcular_niveis_retorno_regional
) |>
  left_join(
    mapa_station_id,
    by = "station_id"
  ) |>
  relocate(
    gauge_code,
    .after = station_id
  )

calcular_niveis_retorno_nonpooled_fit <- function(
    fit,
    gauge_code
) {
  
  alpha <- as.numeric(
    extrair_primeiro_parametro(fit, "alpha")
  )
  
  beta <- as.numeric(
    extrair_primeiro_parametro(fit, "beta")
  )
  
  H <- as.numeric(
    extrair_primeiro_parametro(fit, "H")
  )
  
  theta <- as.numeric(
    extrair_primeiro_parametro(fit, "theta")
  )
  
  # ----------------------------------------------------------
  # Xi por duração
  # ----------------------------------------------------------
  
  nome_xi <- localizar_parametro(
    fit,
    c(
      "xi",
      "xi_d",
      "xi_k",
      "xi_bar",
      "xi_j"
    )
  )
  
  if (is.null(nome_xi)) {
    stop(
      "Não encontrei xi na estação ",
      gauge_code
    )
  }
  
  xi <- rstan::extract(
    fit,
    pars = nome_xi,
    permuted = TRUE
  )[[nome_xi]]
  
  # ----------------------------------------------------------
  # Garantir draw × duração
  # ----------------------------------------------------------
  
  if (is.null(dim(xi))) {
    
    # xi único: repetir para 24 durações
    xi <- matrix(
      xi,
      nrow = length(xi),
      ncol = length(D_GRID)
    )
    
  } else if (
    length(dim(xi)) == 2 &&
    nrow(xi) == length(alpha) &&
    ncol(xi) == length(D_GRID)
  ) {
    
    # já está draw × duração
    xi <- xi
    
  } else if (
    length(dim(xi)) == 2 &&
    nrow(xi) == length(D_GRID) &&
    ncol(xi) == length(alpha)
  ) {
    
    # duração × draw -> transpor
    xi <- t(xi)
    
  } else {
    
    stop(
      "Dimensão inesperada de xi na estação ",
      gauge_code,
      ": ",
      paste(dim(xi), collapse = " x ")
    )
  }
  
  # ----------------------------------------------------------
  # Níveis de retorno
  # ----------------------------------------------------------
  
  resultados <- vector(
    "list",
    length(D_GRID) * length(TR_GRID)
  )
  
  contador <- 1L
  
  for (k in seq_along(D_GRID)) {
    
    d <- D_GRID[k]
    
    fator <- (
      (d + theta) /
        (D_REF + theta)
    )^(-H)
    
    mu <- exp(alpha) * fator
    sigma <- exp(beta) * fator
    
    # xi específico da duração k
    xi_d <- as.numeric(
      xi[, k]
    )
    
    for (Tret in TR_GRID) {
      
      rl <- qgev_retorno(
        mu = mu,
        sigma = sigma,
        xi = xi_d,
        T = Tret
      )
      
      resultados[[contador]] <- tibble(
        modelo = "nonpooled",
        gauge_code = gauge_code,
        duracao_h = d,
        return_period_years = Tret,
        
        median = median(
          rl,
          na.rm = TRUE
        ),
        
        q025 = quantile(
          rl,
          0.025,
          na.rm = TRUE
        ),
        
        q975 = quantile(
          rl,
          0.975,
          na.rm = TRUE
        )
      )
      
      contador <- contador + 1L
    }
  }
  
  bind_rows(resultados)
}

niveis_retorno_nonpooled <- purrr::imap_dfr(
  fits_nonpooled,
  calcular_niveis_retorno_nonpooled_fit
) |>
  left_join(
    mapa_station_id,
    by = "gauge_code"
  ) |>
  relocate(
    station_id,
    .after = modelo
  )

niveis_retorno_todos <- bind_rows(
  niveis_retorno_regionais,
  niveis_retorno_nonpooled
) |>
  mutate(
    nome_modelo = factor(
      NOMES_MODELOS_TEXTO[modelo],
      levels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS]
    )
  )

niveis_retorno_todos |>
  count(modelo)

salvar_tabela(
  niveis_retorno_todos,
  "niveis_retorno_cinco_modelos.csv"
)

saveRDS(
  niveis_retorno_todos,
  file.path(
    DIR_RESULTADOS,
    "niveis_retorno_cinco_modelos.rds"
  )
)

diagnosticos_gpkg <- diagnosticos_gpkg |>
  left_join(
    mapa_station_id,
    by = "station_id"
  ) |>
  select(
    modelo,
    station_id,
    gauge_code,
    flag_diagnostico
  )

diagnosticos_sf <- diagnosticos_gpkg |>
  left_join(
    meta_estacoes_modelo |>
      select(
        gauge_code,
        lat,
        long,
        elevation,
        network,
        responsible
      ),
    by = "gauge_code"
  ) |>
  sf::st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4674,
    remove = FALSE
  )

diagnosticos_sf |>
  sf::st_drop_geometry() |>
  summarise(
    n = n(),
    n_sem_lat = sum(is.na(lat)),
    n_sem_long = sum(is.na(long)),
    n_flag = sum(
      flag_diagnostico %in% TRUE,
      na.rm = TRUE
    )
  )

arquivo_gpkg_parametros <- file.path(
  DIR_GPKG,
  "RS_parametros_modelos.gpkg"
)

if (file.exists(arquivo_gpkg_parametros)) {
  file.remove(arquivo_gpkg_parametros)
}

sf::st_write(
  estruturais_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "parametros_estruturais",
  quiet = TRUE
)

sf::st_write(
  parametros_gev_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "parametros_gev",
  append = TRUE,
  quiet = TRUE
)

sf::st_write(
  diagnosticos_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "diagnosticos",
  append = TRUE,
  quiet = TRUE
)

niveis_retorno_gpkg <- niveis_retorno_todos |>
  
  left_join(
    diagnosticos_gpkg |>
      select(
        modelo,
        station_id,
        flag_diagnostico
      ),
    by = c(
      "modelo",
      "station_id"
    )
  ) |>
  
  left_join(
    meta_estacoes_modelo |>
      select(
        gauge_code,
        lat,
        long,
        elevation,
        network,
        responsible
      ),
    by = "gauge_code"
  )

niveis_retorno_gpkg |>
  summarise(
    n = n(),
    n_estacoes = n_distinct(gauge_code),
    n_modelos = n_distinct(modelo),
    n_duracoes = n_distinct(duracao_h),
    n_TR = n_distinct(return_period_years),
    
    n_flag = sum(
      flag_diagnostico %in% TRUE,
      na.rm = TRUE
    ),
    
    sem_lat = sum(is.na(lat)),
    sem_long = sum(is.na(long))
  )

niveis_retorno_sf <- niveis_retorno_gpkg |>
  sf::st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4674,
    remove = FALSE
  )

arquivo_gpkg_retorno <- file.path(
  DIR_GPKG,
  "RS_periodos_retorno.gpkg"
)

if (file.exists(arquivo_gpkg_retorno)) {
  file.remove(arquivo_gpkg_retorno)
}

sf::st_write(
  niveis_retorno_sf,
  dsn = arquivo_gpkg_retorno,
  layer = "niveis_retorno",
  quiet = TRUE
)

# =============================================================================
# 3. ESTAÇÕES PARA A FIGURA IDF
# =============================================================================

modelo_idf_figura <- "bhm_xi_j"

# Xi por estação no modelo escolhido
xi_estacoes <- parametros_gev_todos |>
  filter(
    modelo == modelo_idf_figura,
    duracao_h == 24
  ) |>
  select(
    gauge_code,
    xi_mediana,
    xi_q025,
    xi_q975
  ) |>
  distinct()

# A801
sta_A801 <- xi_estacoes |>
  filter(
    gauge_code == "A801"
  ) |>
  mutate(
    tipo = "A801"
  )


# Cauda mais leve
sta_cauda_leve <- xi_estacoes |>
  filter(
    gauge_code != "A801"
  ) |>
  arrange(
    xi_mediana
  ) |>
  slice(1) |>
  mutate(
    tipo = "Cauda leve"
  )


# Cauda mais pesada
sta_cauda_pesada <- xi_estacoes |>
  filter(
    gauge_code != "A801"
  ) |>
  arrange(
    desc(xi_mediana)
  ) |>
  slice(1) |>
  mutate(
    tipo = "Cauda pesada"
  )

estacoes_exemplo <- bind_rows(
  sta_A801,
  sta_cauda_leve,
  sta_cauda_pesada
)

print(
  estacoes_exemplo,
  n = Inf,
  width = Inf
)


TR_FIG <- c(
  2,
  10,
  25,
  50,
  100,
  200
)

idf_exemplos_plot <- niveis_retorno_todos |>
  filter(
    modelo == modelo_idf_figura,
    gauge_code %in% estacoes_exemplo$gauge_code,
    return_period_years %in% TR_FIG
  ) |>
  inner_join(
    estacoes_exemplo |>
      select(
        gauge_code,
        tipo,
        xi_mediana
      ),
    by = "gauge_code"
  ) |>
  mutate(
    
    painel = paste0(
      gauge_code,
      "\n",
      tipo
    ),
    
    TR = factor(
      return_period_years,
      levels = TR_FIG
    )
  )

param_texto_idf <- parametros_estruturais_todos |>
  filter(
    modelo == modelo_idf_figura,
    gauge_code %in% estacoes_exemplo$gauge_code,
    parametro %in% c(
      "H",
      "theta"
    )
  ) |>
  select(
    gauge_code,
    parametro,
    mediana
  ) |>
  tidyr::pivot_wider(
    names_from = parametro,
    values_from = mediana
  ) |>
  left_join(
    estacoes_exemplo |>
      select(
        gauge_code,
        tipo,
        xi_mediana
      ),
    by = "gauge_code"
  ) |>
  mutate(
    
    painel = paste0(
      gauge_code,
      "\n",
      tipo
    ),
    
    texto_parametros = paste0(
      "\u03be = ",
      sprintf("%.3f", xi_mediana),
      "\nH = ",
      sprintf("%.3f", H),
      "\n\u03b8 = ",
      sprintf("%.3f", theta)
    )
  )

arquivo_gpkg_retorno <- file.path(
  DIR_GPKG,
  "RS_periodos_retorno.gpkg"
)

if (file.exists(arquivo_gpkg_retorno)) {
  file.remove(arquivo_gpkg_retorno)
}

sf::st_write(
  niveis_retorno_sf,
  dsn = arquivo_gpkg_retorno,
  layer = "niveis_retorno",
  quiet = TRUE
)

estacoes_exemplo |>
  select(
    gauge_code,
    tipo,
    xi_mediana,
    xi_q025,
    xi_q975
  )

modelo_idf_figura <- "bhm_xi_j"

idf_exemplos_plot <- niveis_retorno_todos |>
  filter(
    modelo == modelo_idf_figura,
    gauge_code %in% estacoes_exemplo$gauge_code,
    return_period_years %in% TR_FIG
  ) |>
  inner_join(
    estacoes_exemplo |>
      select(
        gauge_code,
        tipo,
        xi_mediana
      ),
    by = "gauge_code"
  ) |>
  mutate(
    painel = paste0(
      gauge_code,
      "\n",
      tipo
    ),
    
    TR = factor(
      return_period_years,
      levels = TR_FIG
    )
  )

param_texto_idf <- parametros_estruturais_todos |>
  filter(
    modelo == modelo_idf_figura,
    gauge_code %in% estacoes_exemplo$gauge_code,
    parametro %in% c(
      "H",
      "theta"
    )
  ) |>
  select(
    gauge_code,
    parametro,
    mediana
  ) |>
  tidyr::pivot_wider(
    names_from = parametro,
    values_from = mediana
  ) |>
  left_join(
    estacoes_exemplo |>
      select(
        gauge_code,
        tipo,
        xi_mediana
      ),
    by = "gauge_code"
  ) |>
  mutate(
    painel = paste0(
      gauge_code,
      "\n",
      tipo
    ),
    
    texto_parametros = paste0(
      "\u03be = ",
      sprintf("%.3f", xi_mediana),
      "\nH = ",
      sprintf("%.3f", H),
      "\n\u03b8 = ",
      sprintf("%.3f", theta)
    )
  )

cores_TR <- c(
  "2"   = "#D55E00",
  "10"  = "#E69F00",
  "25"  = "#009E73",
  "50"  = "#0072B2",
  "100" = "#CC79A7",
  "200" = "#7B3294"
)
estacoes_plot_idf <- estacoes_exemplo |>
  select(gauge_code) |>
  distinct()

idf_exemplos_plot <- niveis_retorno_todos |>
  filter(
    gauge_code %in% estacoes_plot_idf$gauge_code,
    return_period_years %in% TR_FIG
  ) |>
  mutate(
    
    gauge_code = factor(
      gauge_code,
      levels = estacoes_plot_idf$gauge_code
    ),
    
    modelo = factor(
      modelo,
      levels = ORDEM_MODELOS
    ),
    
    nome_modelo = factor(
      NOMES_MODELOS_TEXTO[as.character(modelo)],
      levels = NOMES_MODELOS_TEXTO[ORDEM_MODELOS]
    ),
    
    TR = factor(
      return_period_years,
      levels = TR_FIG
    )
  )
library(ggh4x)
p_idf_modelos <- ggplot(
  idf_exemplos_plot,
  aes(
    x = duracao_h,
    y = median,
    group = TR,
    color = TR
  )
) +
  
  geom_line(
    linewidth = 0.65
  ) +
  
  facet_grid2(
    rows = vars(nome_modelo),
    cols = vars(gauge_code),
    scales = "free_y",
    independent = "y"
  ) +
  
  scale_x_log10(
    breaks = c(
      1, 2, 3, 6, 12, 24
    ),
    labels = c(
      "1", "2", "3", "6", "12", "24"
    )
  ) +
  
  scale_y_log10(
    labels = scales::label_number(
      accuracy = 0.1
    )
  ) +
  
  scale_color_manual(
    values = cores_TR,
    breaks = as.character(TR_FIG),
    name = "Tempo de retorno (anos)"
  ) +
  
  labs(
    x = "Duração (h)",
    y = "Intensidade (mm/h)"
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "left"
    )
  ) +
  
  tema_relatorio +
  
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text.x = element_text(
      family = font,
      size = fontsize
    ),
    strip.text.y = element_text(
      family = font,
      size = fontsize
    ),
    panel.grid = element_blank(),
    panel.spacing = unit(4, "pt")
  )

p_idf_modelos
salvar_figura(
  p_idf_modelos,
  "curvas_IDF_tres_estacoes_cinco_modelos.png",
  largura = 15,
  altura = 20
)

TR_FIG_1 <- c(2, 25, 200)
TR_FIG_2 <- c(10, 50, 100)

plot_idf_modelos_ic <- function(dados, trs) {
  
  dados_plot <- dados |>
    filter(return_period_years %in% trs) |>
    mutate(
      TR = factor(
        return_period_years,
        levels = trs
      )
    )
  
  ggplot(
    dados_plot,
    aes(
      x = duracao_h,
      y = median,
      group = TR,
      color = TR,
      fill = TR
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = q025,
        ymax = q975
      ),
      alpha = 0.10,
      color = NA
    ) +
    
    geom_line(
      linewidth = 0.65
    ) +
    
    facet_grid2(
      rows = vars(nome_modelo),
      cols = vars(gauge_code),
      scales = "free_y",
      independent = "y"
    ) +
    
    scale_x_log10(
      breaks = c(1, 2, 3, 6, 12, 24),
      labels = c("1", "2", "3", "6", "12", "24")
    ) +
    
    scale_y_log10(
      labels = scales::label_number(
        accuracy = 0.1
      )
    ) +
    
    scale_color_manual(
      values = cores_TR[as.character(trs)],
      breaks = as.character(trs),
      name = "Tempo de retorno (anos)"
    ) +
    
    scale_fill_manual(
      values = cores_TR[as.character(trs)],
      breaks = as.character(trs),
      name = "Tempo de retorno (anos)"
    ) +
    
    labs(
      x = "Duração (h)",
      y = "Intensidade (mm/h)"
    ) +
    
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "left"
      ),
      fill = "none"
    ) +
    
    tema_relatorio +
    
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text.x = element_text(
        family = font,
        size = fontsize
      ),
      strip.text.y = element_text(
        family = font,
        size = fontsize
      ),
      panel.grid = element_blank(),
      panel.spacing = unit(4, "pt")
    )
}

p_idf_modelos_ic_1 <- plot_idf_modelos_ic(
  dados = idf_exemplos_plot,
  trs = TR_FIG_1
)

p_idf_modelos_ic_1

p_idf_modelos_ic_2 <- plot_idf_modelos_ic(
  dados = idf_exemplos_plot,
  trs = TR_FIG_2
)

p_idf_modelos_ic_2

salvar_figura(
  p_idf_modelos_ic_1,
  "curvas_IDF_tres_estacoes_cinco_modelos_IC_TR_2_25_200.png",
  largura = 15,
  altura = 20
)

salvar_figura(
  p_idf_modelos_ic_2,
  "curvas_IDF_tres_estacoes_cinco_modelos_IC_TR_10_50_100.png",
  largura = 15,
  altura = 20
)

resumo_estruturais_modelos <- parametros_estruturais_todos |>
  group_by(
    modelo,
    parametro
  ) |>
  summarise(
    n_estacoes = n_distinct(gauge_code),
    
    mediana = median(
      mediana,
      na.rm = TRUE
    ),
    
    q05 = quantile(
      mediana,
      0.05,
      na.rm = TRUE
    ),
    
    q25 = quantile(
      mediana,
      0.25,
      na.rm = TRUE
    ),
    
    q75 = quantile(
      mediana,
      0.75,
      na.rm = TRUE
    ),
    
    q95 = quantile(
      mediana,
      0.95,
      na.rm = TRUE
    ),
    
    minimo = min(
      mediana,
      na.rm = TRUE
    ),
    
    maximo = max(
      mediana,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  mutate(
    nome_modelo =
      NOMES_MODELOS_TEXTO[modelo]
  )

resumo_estruturais_modelos

resumo_H_theta <- parametros_estruturais_todos |>
  filter(
    parametro %in% c(
      "H",
      "theta"
    )
  ) |>
  group_by(
    modelo,
    parametro
  ) |>
  summarise(
    mediana_espacial =
      median(
        mediana,
        na.rm = TRUE
      ),
    
    IQR =
      IQR(
        mediana,
        na.rm = TRUE
      ),
    
    minimo =
      min(
        mediana,
        na.rm = TRUE
      ),
    
    maximo =
      max(
        mediana,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  mutate(
    nome_modelo =
      NOMES_MODELOS_TEXTO[modelo]
  )

resumo_H_theta

resumo_gev_24h <- parametros_gev_todos |>
  filter(
    duracao_h == 24
  ) |>
  group_by(
    modelo
  ) |>
  summarise(
    
    mu_med =
      median(mu_mediana, na.rm = TRUE),
    
    mu_q05 =
      quantile(mu_mediana, 0.05, na.rm = TRUE),
    
    mu_q95 =
      quantile(mu_mediana, 0.95, na.rm = TRUE),
    
    
    sigma_med =
      median(sigma_mediana, na.rm = TRUE),
    
    sigma_q05 =
      quantile(sigma_mediana, 0.05, na.rm = TRUE),
    
    sigma_q95 =
      quantile(sigma_mediana, 0.95, na.rm = TRUE),
    
    
    xi_med =
      median(xi_mediana, na.rm = TRUE),
    
    xi_q05 =
      quantile(xi_mediana, 0.05, na.rm = TRUE),
    
    xi_q95 =
      quantile(xi_mediana, 0.95, na.rm = TRUE),
    
    xi_min =
      min(xi_mediana, na.rm = TRUE),
    
    xi_max =
      max(xi_mediana, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  mutate(
    nome_modelo =
      NOMES_MODELOS_TEXTO[modelo]
  )

resumo_gev_24h

resumo_sinal_xi_24h <- parametros_gev_todos |>
  filter(
    duracao_h == 24
  ) |>
  group_by(
    modelo
  ) |>
  summarise(
    n_estacoes = n(),
    
    n_xi_negativo =
      sum(
        xi_mediana < 0,
        na.rm = TRUE
      ),
    
    n_xi_positivo =
      sum(
        xi_mediana > 0,
        na.rm = TRUE
      ),
    
    n_xi_proximo_zero =
      sum(
        abs(xi_mediana) <= 0.05,
        na.rm = TRUE
      ),
    
    perc_negativo =
      100 * n_xi_negativo / n_estacoes,
    
    perc_positivo =
      100 * n_xi_positivo / n_estacoes,
    
    perc_proximo_zero =
      100 * n_xi_proximo_zero / n_estacoes,
    
    .groups = "drop"
  ) |>
  mutate(
    nome_modelo =
      NOMES_MODELOS_TEXTO[modelo]
  )

resumo_sinal_xi_24h

resumo_xi_duracao <- parametros_gev_todos |>
  group_by(
    modelo,
    duracao_h
  ) |>
  summarise(
    xi_mediana_espacial =
      median(
        xi_mediana,
        na.rm = TRUE
      ),
    
    xi_q05 =
      quantile(
        xi_mediana,
        0.05,
        na.rm = TRUE
      ),
    
    xi_q95 =
      quantile(
        xi_mediana,
        0.95,
        na.rm = TRUE
      ),
    
    xi_min =
      min(
        xi_mediana,
        na.rm = TRUE
      ),
    
    xi_max =
      max(
        xi_mediana,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  mutate(
    nome_modelo =
      NOMES_MODELOS_TEXTO[modelo]
  )

resumo_xi_duracao

resumo_xi_duracoes_chave <- resumo_xi_duracao |>
  filter(
    duracao_h %in% c(
      1,
      3,
      6,
      12,
      24
    )
  )

resumo_xi_duracoes_chave

resumo_A801 <- parametros_gev_todos |>
  filter(
    gauge_code == "A801",
    duracao_h %in% c(
      1,
      6,
      12,
      24
    )
  ) |>
  select(
    modelo,
    duracao_h,
    mu_mediana,
    sigma_mediana,
    xi_mediana
  ) |>
  mutate(
    nome_modelo =
      NOMES_MODELOS_TEXTO[modelo]
  ) |>
  arrange(
    duracao_h,
    factor(
      modelo,
      levels = ORDEM_MODELOS
    )
  )

resumo_A801
