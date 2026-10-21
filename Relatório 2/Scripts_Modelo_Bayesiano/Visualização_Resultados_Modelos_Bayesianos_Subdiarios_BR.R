# =============================================================================
# RESULTADOS DOS MODELOS BAYESIANOS SUBDIÁRIOS — BRASIL
#
# Objetivos:
#   1. localizar todos os fits
#   2. identificar modelo conjunto ou complementar
#   3. extrair diagnósticos
#   4. extrair parâmetros estruturais
#   5. extrair parâmetros GEV por duração
#   6. posteriormente gerar mapas e GeoPackages
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(tibble)
library(rstan)
library(sf)

# =============================================================================
# 1. CAMINHOS
# =============================================================================

DIR_MODELOS <- "C:/Larissa/modelos_bayesianos"

DIR_CONJUNTO <- file.path(DIR_MODELOS, "strict_sem_iqr", "10min_24h_conjunto")

DIR_COMPLEMENTAR <- file.path(DIR_MODELOS, "1h_24h_complementar")

DIR_RESULTADOS <- file.path(DIR_MODELOS, "resultados_subdiarios")

DIR_GPKG <- file.path(DIR_RESULTADOS, "geopackage")

dir.create(DIR_RESULTADOS, recursive = TRUE, showWarnings = FALSE)

dir.create(DIR_GPKG, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. INVENTÁRIO DOS FITS
# =============================================================================

arquivos_conjunto <- list.files(DIR_CONJUNTO, pattern = "_fit\\.rds$", recursive = TRUE, full.names = TRUE)

arquivos_complementar <- list.files(DIR_COMPLEMENTAR, pattern = "_fit\\.rds$", recursive = TRUE,
                                    full.names = TRUE)

inventario_fits <- bind_rows(tibble(arquivo = arquivos_conjunto, tipo_modelo = "conjunto_10min_24h"), tibble(
  arquivo = arquivos_complementar, tipo_modelo = "complementar_1h_24h")) |>
  mutate(gauge_code = basename(arquivo) |> str_remove("_fit\\.rds$") |> trimws()) |>
  arrange(gauge_code, tipo_modelo)

inventario_fits |>
  count(tipo_modelo)

cat("\nTotal de fits encontrados:", nrow(inventario_fits), "\n")

cat("Estações distintas:", n_distinct(inventario_fits$gauge_code), "\n")

# =============================================================================
# 3. FUNÇÃO DE RESUMO POSTERIOR
# =============================================================================

resumo_posterior <- function(x) {
  x <- as.numeric(x)
  tibble(mediana = median(x, na.rm = TRUE), q025 = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
         q975 = quantile(x, 0.975, na.rm = TRUE, names = FALSE))
}

parametros_fit <- function(fit) {
  fit@sim$pars_oi
}

# =============================================================================
# 4. EXTRAIR PARÂMETROS ESTRUTURAIS
# =============================================================================

extrair_estruturais_fit <- function(arquivo, gauge_code, tipo_modelo) {
  fit <- readRDS(arquivo)
  pars <- fit@sim$pars_oi
  if (tipo_modelo == "conjunto_10min_24h") {
    nomes <- c("alpha", "beta", "H_long", "H_short", "theta_long", "theta_short", "xi_bar", "tau2_mu",
               "tau2_sigma", "tau2_xi")
  } else {
    nomes <- c("alpha", "beta", "H", "theta", "xi_bar", "tau2_mu", "tau2_sigma", "tau2_xi")
  }
  nomes <- intersect(nomes, pars)
  map_dfr(nomes, function(par) {
    draws <- rstan::extract(fit, pars = par, permuted = TRUE)[[par]]
    resumo_posterior(draws) |>
      mutate(gauge_code = as.character(gauge_code), tipo_modelo = tipo_modelo, parametro = par, .before = 1)
  })
}

parametros_estruturais_subdiarios <- pmap_dfr(inventario_fits, function(arquivo, tipo_modelo, gauge_code) {
  message(gauge_code, " | ", tipo_modelo)
  tryCatch(extrair_estruturais_fit(arquivo = arquivo, gauge_code = gauge_code, tipo_modelo = tipo_modelo),
           error = function(e) {
             warning("Erro em ", gauge_code, ": ", conditionMessage(e))
             NULL
           })
})

parametros_estruturais_subdiarios |>
  count(tipo_modelo, parametro) |>
  print(n = Inf)

F_IMAX <- paste0("C:/Larissa/sensibilidade_qc/dados_finais/", "imax_sem_iqr_subdiario_Uniplu.rds")

MIN_YEARS <- 8L

DUR_LONGO_MIN <- (1:24) * 60

DUR_CURTO_MIN <- c(10, 15, 20, 30, 40, 45, 50)

bruto <- readRDS(F_IMAX)

dados_duracao <- bruto |>
  transmute(gauge_code = trimws(as.character(gauge_code)), duration_min = as.integer(duration_min),
            duration_h = as.numeric(duration_min) / 60, year = as.integer(wateryear), imax = as.numeric(imax_mm_h)) |>
  filter(is.finite(imax),
         #imax > 0,
         is.finite(duration_h), duration_h > 0, !is.na(year), !is.na(gauge_code), nzchar(gauge_code))

rm(bruto)
invisible(gc())

years_curto <- dados_duracao |>
  filter(duration_min %in% DUR_CURTO_MIN) |>
  distinct(gauge_code, duration_min, year) |>
  count(gauge_code, duration_min, name = "n_years")

pares_curto_validos <- years_curto |>
  filter(n_years >= MIN_YEARS)

duracoes_conjunto <- bind_rows(
  # sub-horárias válidas
  pares_curto_validos |> transmute(gauge_code, duracao_min = duration_min),
  # 1–24 h
  inventario_fits |> filter(tipo_modelo == "conjunto_10min_24h") |> distinct(gauge_code) |> tidyr::crossing(
    duracao_min = DUR_LONGO_MIN)) |>
  semi_join(inventario_fits |> filter(tipo_modelo == "conjunto_10min_24h") |> distinct(gauge_code),
            by = "gauge_code") |>
  distinct(gauge_code, duracao_min) |>
  mutate(duracao_h = duracao_min / 60) |>
  arrange(gauge_code, duracao_min)

duracoes_complementar <- inventario_fits |>
  filter(tipo_modelo == "complementar_1h_24h") |>
  distinct(gauge_code) |>
  tidyr::crossing(duracao_min = DUR_LONGO_MIN) |>
  mutate(duracao_h = duracao_min / 60) |>
  arrange(gauge_code, duracao_min)

duracoes_por_fit <- bind_rows(duracoes_conjunto |> mutate(tipo_modelo = "conjunto_10min_24h"),
                              duracoes_complementar |> mutate(tipo_modelo = "complementar_1h_24h")) |>
  arrange(tipo_modelo, gauge_code, duracao_min)

library(dplyr)
library(purrr)
library(tibble)
library(rstan)

# =============================================================================
# 5.1. FUNÇÃO AUXILIAR — RESUMO POSTERIOR
# =============================================================================

resumir_draws <- function(x) {
  x <- as.numeric(x)
  tibble(mediana = median(x, na.rm = TRUE), q025 = quantile(x, probs = 0.025, na.rm = TRUE, names = FALSE),
         q975 = quantile(x, probs = 0.975, na.rm = TRUE, names = FALSE))
}

extrair_duracoes_fit <- function(gauge_code, tipo_modelo) {
  d <- duracoes_por_fit |>
    filter(.data$gauge_code == .env$gauge_code, .data$tipo_modelo == .env$tipo_modelo) |>
    arrange(duracao_min)
  if (!nrow(d)) {
    stop("Nenhuma duração encontrada para ", gauge_code, " | ", tipo_modelo)
  }
  d |>
    mutate(duration_id = row_number()) |>
    select(duration_id, duracao_min, duracao_h)
}

extrair_gev_fit <- function(fit, gauge_code, tipo_modelo) {
  duracoes <- extrair_duracoes_fit(gauge_code = gauge_code, tipo_modelo = tipo_modelo)
  draws <- rstan::extract(fit, pars = c("mu", "sigma", "xi_d"), permuted = TRUE)
  mu <- draws$mu
  sigma <- draws$sigma
  xi <- draws$xi_d
  K <- nrow(duracoes)
  if (is.null(dim(mu))) {
    mu <- matrix(mu, ncol = K)
  }
  if (is.null(dim(sigma))) {
    sigma <- matrix(sigma, ncol = K)
  }
  if (is.null(dim(xi))) {
    xi <- matrix(xi, ncol = K)
  }
  if (ncol(mu) != K || ncol(sigma) != K || ncol(xi) != K) {
    stop("Incompatibilidade entre K e os parâmetros: ", gauge_code, " | ", tipo_modelo, " | K esperado = ", K,
         " | mu = ", ncol(mu), " | sigma = ", ncol(sigma), " | xi = ", ncol(xi))
  }
  map_dfr(seq_len(K), function(k) {
    mu_res <- resumir_draws(mu[, k])
    sigma_res <- resumir_draws(sigma[, k])
    xi_res <- resumir_draws(xi[, k])
    tibble(gauge_code = as.character(gauge_code), tipo_modelo = tipo_modelo, duration_id = k, duracao_min =
             duracoes$duracao_min[k], duracao_h = duracoes$duracao_h[k], mu_mediana = mu_res$mediana, mu_q025 =
             mu_res$q025, mu_q975 = mu_res$q975, sigma_mediana = sigma_res$mediana, sigma_q025 = sigma_res$q025,
           sigma_q975 = sigma_res$q975, xi_mediana = xi_res$mediana, xi_q025 = xi_res$q025, xi_q975 = xi_res$q975)
  })
}

extrair_diagnostico_fit <- function(fit, gauge_code, tipo_modelo, max_treedepth_ref = 17L) {
  
  # ---------------------------------------------------------------------------
  # resumo Stan
  # ---------------------------------------------------------------------------
  sm <- summary(fit)$summary
  rhat <- sm[, "Rhat"]
  ess <- sm[, "n_eff"]
  rhat_validos <- rhat[is.finite(rhat)]
  ess_validos <- ess[is.finite(ess)]
  rhat_max <- if (length(rhat_validos)) {
    max(rhat_validos, na.rm = TRUE)
  } else {
    NA_real_
  }
  ess_min <- if (length(ess_validos)) {
    min(ess_validos, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # parâmetros do sampler
  # ---------------------------------------------------------------------------
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_divergencias <- sum(vapply(sampler, function(x) {
    sum(x[, "divergent__"], na.rm = TRUE)
  }, numeric(1)))
  n_treedepth_max <- sum(vapply(sampler, function(x) {
    sum(x[, "treedepth__"] >= max_treedepth_ref, na.rm = TRUE)
  }, numeric(1)))
  
  # ---------------------------------------------------------------------------
  # quantidade real de cadeias
  # ---------------------------------------------------------------------------
  n_chains <- length(sampler)
  
  # ---------------------------------------------------------------------------
  # flags
  # ---------------------------------------------------------------------------
  flag_divergencia <-
    is.finite(n_divergencias) &&
    n_divergencias > 0
  flag_rhat_105 <-
    is.finite(rhat_max) &&
    rhat_max > 1.05
  flag_cadeias <-
    is.na(n_chains) ||
    n_chains < 5L
  flag_diagnostico <-
    flag_divergencia |
    flag_rhat_105 |
    flag_cadeias
  tibble(gauge_code = as.character(gauge_code), tipo_modelo = tipo_modelo, n_chains = n_chains, rhat_max =
           rhat_max, ess_min = ess_min, n_divergencias = n_divergencias, n_treedepth_max = n_treedepth_max,
         flag_rhat_105 = flag_rhat_105, flag_divergencia = flag_divergencia, flag_cadeias = flag_cadeias,
         flag_diagnostico = flag_diagnostico)
}

processar_fit_subdiario <- function(arquivo, tipo_modelo, gauge_code) {
  tryCatch({
    fit <- readRDS(arquivo)
    gev <- extrair_gev_fit(fit = fit, gauge_code = gauge_code, tipo_modelo = tipo_modelo)
    diagnostico <- extrair_diagnostico_fit(fit = fit, gauge_code = gauge_code, tipo_modelo = tipo_modelo,
                                           max_treedepth_ref = 17L)
    list(gev = gev, diagnostico = diagnostico, erro = NULL)
  }, error = function(e) {
    list(gev = NULL, diagnostico = NULL, erro = tibble(gauge_code = as.character(gauge_code), tipo_modelo =
                                                         tipo_modelo, arquivo = arquivo, erro = conditionMessage(e)))
  })
}

resultados_fits <- pmap(inventario_fits, function(arquivo, tipo_modelo, gauge_code) {
  message(gauge_code, " | ", tipo_modelo)
  processar_fit_subdiario(arquivo = arquivo, tipo_modelo = tipo_modelo, gauge_code = gauge_code)
})

parametros_gev_subdiarios <- map_dfr(resultados_fits, "gev")

diagnosticos_subdiarios <- map_dfr(resultados_fits, "diagnostico")

erros_extracao_subdiarios <- map_dfr(resultados_fits, "erro")

parametros_gev_subdiarios |>
  group_by(tipo_modelo) |>
  summarise(n_estacoes = n_distinct(gauge_code), n_linhas = n(), duracao_min_min = min(duracao_min, na.rm = TRUE
  ), duracao_min_max = max(duracao_min, na.rm = TRUE), n_duracoes_distintas = n_distinct(duracao_min),
  .groups = "drop")

resumo_diagnosticos_subdiarios <-
  diagnosticos_subdiarios |>
  summarise(n_estacoes = n(), n_divergencia = sum(flag_divergencia, na.rm = TRUE), pct_divergencia = 100 * mean(
    flag_divergencia, na.rm = TRUE), n_rhat_105 = sum(flag_rhat_105, na.rm = TRUE), n_cadeias_incompletas =
      sum(flag_cadeias, na.rm = TRUE), n_flag_diagnostico = sum(flag_diagnostico, na.rm = TRUE),
    pct_flag_diagnostico = 100 * mean(flag_diagnostico, na.rm = TRUE), rhat_max = max(rhat_max, na.rm = TRUE),
    ess_min = min(ess_min, na.rm = TRUE), divergencias_total = sum(n_divergencias, na.rm = TRUE))

print(resumo_diagnosticos_subdiarios, width = Inf)

diagnosticos_subdiarios |>
  group_by(tipo_modelo) |>
  summarise(n_estacoes = n(), n_divergencia = sum(flag_divergencia, na.rm = TRUE), n_rhat_105 = sum(
    flag_rhat_105, na.rm = TRUE), n_cadeias_incompletas = sum(flag_cadeias, na.rm = TRUE), n_flag = sum(
      flag_diagnostico, na.rm = TRUE), .groups = "drop")

saveRDS(parametros_gev_subdiarios, file.path(DIR_RESULTADOS, "parametros_GEV_subdiarios.rds"))

saveRDS(diagnosticos_subdiarios, file.path(DIR_RESULTADOS, "diagnosticos_subdiarios.rds"))

write_excel_csv(diagnosticos_subdiarios, file.path(DIR_RESULTADOS, "diagnosticos_subdiarios.csv"))

if (nrow(erros_extracao_subdiarios) > 0) {
  write_excel_csv(erros_extracao_subdiarios, file.path(DIR_RESULTADOS, "erros_extracao_subdiarios.csv"))
}

extrair_estruturais_subdiario <- function(fit, gauge_code, tipo_modelo) {
  if (tipo_modelo == "conjunto_10min_24h") {
    especificacao <- tibble::tribble(~parametro_modelo, ~parametro,   ~ramo,
                                     "alpha",           "alpha",      "global", "beta",            "beta",       "global",
                                     "H_long",          "H",          "longo", "H_short",         "H",          "curto",
                                     "theta_long",      "theta",      "longo", "theta_short",     "theta",      "curto",
                                     "xi_bar",          "xi_bar",     "global", "tau2_mu",         "tau2_mu",    "global",
                                     "tau2_sigma",      "tau2_sigma", "global", "tau2_xi",         "tau2_xi",    "global")
  } else if (tipo_modelo == "complementar_1h_24h") {
    especificacao <- tibble::tribble(~parametro_modelo, ~parametro,   ~ramo,
                                     "alpha",           "alpha",      "global", "beta",            "beta",       "global",
                                     "H",               "H",          "longo", "theta",           "theta",      "longo",
                                     "xi_bar",          "xi_bar",     "global", "tau2_mu",         "tau2_mu",    "global",
                                     "tau2_sigma",      "tau2_sigma", "global", "tau2_xi",         "tau2_xi",    "global")
  } else {
    stop("tipo_modelo desconhecido: ", tipo_modelo)
  }
  pars_fit <- fit@sim$pars_oi
  especificacao |>
    filter(parametro_modelo %in% pars_fit) |>
    pmap_dfr(function(parametro_modelo, parametro, ramo) {
      draws <- rstan::extract(fit, pars = parametro_modelo, permuted = TRUE)[[parametro_modelo]]
      resumir_draws(draws) |>
        mutate(gauge_code = as.character(gauge_code), tipo_modelo = tipo_modelo, parametro_modelo =
                 parametro_modelo, parametro = parametro, ramo = ramo, .before = 1)
    })
}

quantil_gev <- function(mu, sigma, xi, T) {
  p <- 1 - 1 / T
  resultado <- numeric(length(mu))
  gumbel <- abs(xi) < 1e-8
  resultado[gumbel] <-
    mu[gumbel] -
    sigma[gumbel] *
    log(-log(p))
  resultado[!gumbel] <-
    mu[!gumbel] +
    sigma[!gumbel] /
    xi[!gumbel] *
    ((-log(p))^(-xi[!gumbel]) - 1)
  resultado
}

TR_SUBDIARIO <- c(2, 5, 10, 25, 50, 100, 200)

extrair_niveis_retorno_fit <- function(fit, gauge_code, tipo_modelo, TR = TR_SUBDIARIO) {
  duracoes <- extrair_duracoes_fit(gauge_code = gauge_code, tipo_modelo = tipo_modelo)
  draws <- rstan::extract(fit, pars = c("mu", "sigma", "xi_d"), permuted = TRUE)
  mu <- draws$mu
  sigma <- draws$sigma
  xi <- draws$xi_d
  K <- nrow(duracoes)
  if (ncol(mu) != K || ncol(sigma) != K || ncol(xi) != K) {
    stop("K incompatível em ", gauge_code, " | ", tipo_modelo)
  }
  purrr::map_dfr(seq_len(K), function(k) {
    purrr::map_dfr(TR, function(T) {
      rl <- quantil_gev(mu = mu[, k], sigma = sigma[, k], xi = xi[, k], T = T)
      resumo <- resumir_draws(rl)
      tibble(gauge_code = as.character(gauge_code), tipo_modelo = tipo_modelo, duration_id = k,
             duracao_min = duracoes$duracao_min[k], duracao_h = duracoes$duracao_h[k], periodo_retorno_anos = T,
             mediana = resumo$mediana, q025 = resumo$q025, q975 = resumo$q975)
    })
  })
}

processar_produtos_subdiarios <- function(arquivo, tipo_modelo, gauge_code) {
  tryCatch({
    fit <- readRDS(arquivo)
    estruturais <-
      extrair_estruturais_subdiario(fit = fit, gauge_code = gauge_code, tipo_modelo = tipo_modelo)
    niveis_retorno <-
      extrair_niveis_retorno_fit(fit = fit, gauge_code = gauge_code, tipo_modelo = tipo_modelo)
    list(estruturais = estruturais, niveis_retorno = niveis_retorno, erro = NULL)
  }, error = function(e) {
    list(estruturais = NULL, niveis_retorno = NULL, erro = tibble(gauge_code = as.character(gauge_code),
                                                                  tipo_modelo = tipo_modelo, arquivo = arquivo, erro = conditionMessage(e)))
  })
}

produtos_subdiarios <- pmap(inventario_fits, function(arquivo, tipo_modelo, gauge_code) {
  message(gauge_code, " | ", tipo_modelo)
  processar_produtos_subdiarios(arquivo = arquivo, tipo_modelo = tipo_modelo, gauge_code = gauge_code)
})

parametros_estruturais_subdiarios <- map_dfr(produtos_subdiarios, "estruturais")

niveis_retorno_subdiarios <- map_dfr(produtos_subdiarios, "niveis_retorno")

erros_produtos_subdiarios <- map_dfr(produtos_subdiarios, "erro")

cat("\nEstruturais — estações:", n_distinct(parametros_estruturais_subdiarios$gauge_code), "\n")

cat("Níveis de retorno — estações:", n_distinct(niveis_retorno_subdiarios$gauge_code), "\n")

cat("Erros:", nrow(erros_produtos_subdiarios), "\n")

saveRDS(parametros_estruturais_subdiarios, file.path(DIR_RESULTADOS, "parametros_estruturais_subdiarios.rds"))

saveRDS(niveis_retorno_subdiarios, file.path(DIR_RESULTADOS, "niveis_retorno_subdiarios.rds"))

write_excel_csv(parametros_estruturais_subdiarios, file.path(DIR_RESULTADOS,
                                                             "parametros_estruturais_subdiarios.csv"))

if (nrow(erros_produtos_subdiarios) > 0) {
  write_excel_csv(erros_produtos_subdiarios, file.path(DIR_RESULTADOS, "erros_produtos_subdiarios.csv"))
}

library(ggplot2)
library(grid)

font <- "sans"
fontsize <- 8

# =============================================================================
# H — COMBINAR O RAMO >= 1 h DOS DOIS TIPOS DE AJUSTE
# =============================================================================

dados_H_sub <- parametros_estruturais_subdiarios |>
  filter(parametro == "H") |>
  mutate(grupo = case_when(
    # Modelo conjunto — ramo longo
    tipo_modelo == "conjunto_10min_24h" & ramo == "longo" ~ "d ≥ 1 h",
    # Modelo complementar — também representa o ramo >= 1 h
    tipo_modelo == "complementar_1h_24h" & ramo == "longo" ~ "d ≥ 1 h",
    # Modelo conjunto — ramo sub-horário
    tipo_modelo == "conjunto_10min_24h" & ramo == "curto" ~ "d < 1 h", TRUE ~ NA_character_)) |>
  filter(!is.na(grupo)) |>
  mutate(grupo = factor(grupo, levels = c("d ≥ 1 h", "d < 1 h")))

p_H_sub <- ggplot(dados_H_sub, aes(x = mediana)) +
  
  geom_histogram(aes(y = after_stat(density)), bins = 40, boundary = 0, fill = "grey80", color = "grey30",
                 linewidth = 0.25) +
  
  geom_density(linewidth = 0.7) +
  
  facet_wrap(~grupo, ncol = 1, scales = "free_y") +
  
  labs(x = expression("Mediana posterior de " * H), y = "Densidade") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3), text = element_text(
    family = font, size = fontsize), axis.text = element_text(family = font, size = fontsize),
    axis.title = element_text(family = font, size = fontsize), strip.background = element_blank(),
    strip.text = element_text(family = font, size = fontsize), panel.border = element_rect(linewidth = 0.4),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_H_sub)

ggsave(file.path(DIR_RESULTADOS, "distribuicao_H_subdiario.png"), p_H_sub, width = 10, height = 9, units = "cm",
       dpi = 600)

dados_theta_sub <- parametros_estruturais_subdiarios |>
  filter(parametro == "theta") |>
  mutate(grupo = case_when(tipo_modelo == "conjunto_10min_24h" & ramo == "longo" ~ "d ≥ 1 h",
                           tipo_modelo == "complementar_1h_24h" & ramo == "longo" ~ "d ≥ 1 h", tipo_modelo == "conjunto_10min_24h" &
                             ramo == "curto" ~ "d < 1 h", TRUE ~ NA_character_)) |>
  filter(!is.na(grupo)) |>
  mutate(grupo = factor(grupo, levels = c("d ≥ 1 h", "d < 1 h")))

p_theta_sub <- ggplot(dados_theta_sub, aes(x = mediana)) +
  
  geom_histogram(aes(y = after_stat(density)), bins = 40, boundary = 0, fill = "grey80", color = "grey30",
                 linewidth = 0.25) +
  
  geom_density(linewidth = 0.7) +
  
  facet_wrap(~grupo, ncol = 1, scales = "free_y") +
  
  labs(x = expression("Mediana posterior de " * theta), y = "Densidade") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3), text = element_text(
    family = font, size = fontsize), axis.text = element_text(family = font, size = fontsize),
    axis.title = element_text(family = font, size = fontsize), strip.background = element_blank(),
    strip.text = element_text(family = font, size = fontsize), panel.border = element_rect(linewidth = 0.4),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_theta_sub)

ggsave(file.path(DIR_RESULTADOS, "distribuicao_theta_subdiario.png"), p_theta_sub, width = 10, height = 9,
       units = "cm", dpi = 600)

dados_xibar_sub <- parametros_estruturais_subdiarios |>
  filter(parametro == "xi_bar") |>
  mutate(modelo = recode(tipo_modelo, "conjunto_10min_24h" = "10 min–24 h", "complementar_1h_24h" = "1–24 h"),
         modelo = factor(modelo, levels = c("10 min–24 h", "1–24 h")))

p_xibar_sub <- ggplot(dados_xibar_sub, aes(x = mediana)) +
  
  geom_histogram(aes(y = after_stat(density)), bins = 40, boundary = 0, fill = "grey80", color = "grey30",
                 linewidth = 0.25) +
  
  geom_density(linewidth = 0.7) +
  
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4) +
  
  facet_wrap(~modelo, ncol = 1, scales = "free_y") +
  
  labs(x = expression("Mediana posterior de " * bar(xi)), y = "Densidade") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3), text = element_text(
    family = font, size = fontsize), axis.text = element_text(family = font, size = fontsize),
    axis.title = element_text(family = font, size = fontsize), strip.background = element_blank(),
    strip.text = element_text(family = font, size = fontsize), panel.border = element_rect(linewidth = 0.4),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_xibar_sub)

ggsave(file.path(DIR_RESULTADOS, "distribuicao_xi_bar_subdiario.png"), p_xibar_sub, width = 10, height = 9,
       units = "cm", dpi = 600)

library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(grid)

font <- "sans"
fontsize <- 8

# =============================================================================
# 1. COMPRIMENTO DAS SÉRIES (usando 24 h como referência)
# =============================================================================

comprimento_series <- dados_duracao |>
  filter(duration_min == 1440) |>
  distinct(gauge_code, year) |>
  count(gauge_code, name = "n_years")

F_META_SUB <- paste0("Relatório 2/Scripts_Maximas/", "resultados/BR/dataframes/",
                     "metadados_estacoes_subdiario_Uniplu.rds")

meta_sub <- readRDS(F_META_SUB)

names(meta_sub)

metadados_subdiarios <- meta_sub |>
  transmute(gauge_code = trimws(as.character(gauge_code)), UF = if ("state" %in% names(meta_sub)) {
    as.character(state)
  } else if ("UF" %in% names(meta_sub)) {
    as.character(UF)
  } else {
    NA_character_
  }, lat = dplyr::coalesce(suppressWarnings(as.numeric(if ("lat" %in% names(meta_sub)) lat else NA)),
                           suppressWarnings(as.numeric(if ("latitude" %in% names(meta_sub)) latitude else NA))),
  long = dplyr::coalesce(suppressWarnings(as.numeric(if ("long" %in% names(meta_sub)) long else NA)),
                         suppressWarnings(as.numeric(if ("longitude" %in% names(meta_sub)) longitude else NA)))) |>
  distinct(gauge_code, .keep_all = TRUE)

# =============================================================================
# 2. CANDIDATAS PARA AS IDFs SUBDIÁRIAS
#    - conjunto: 10 min–24 h
#    - complementar: 1–24 h
# =============================================================================

# =============================================================================
# 2. CANDIDATAS PARA AS IDFs SUBDIÁRIAS
#    Dois tipos de ajuste:
#    - conjunto: 10 min–24 h
#    - complementar: 1–24 h
#
# Critério:
#    1. comprimento da série
#    2. diversidade de cauda
#    3. distância espacial
#
# Maior peso para estações com séries longas
# =============================================================================

candidatas_idf_sub <- parametros_gev_subdiarios |>
  filter(tipo_modelo %in% c("conjunto_10min_24h", "complementar_1h_24h"), duracao_h == 24) |>
  
  select(gauge_code, tipo_modelo, xi_mediana) |>
  
  distinct() |>
  
  left_join(comprimento_series, by = "gauge_code") |>
  
  left_join(metadados_subdiarios |> select(gauge_code, UF, lat, long), by = "gauge_code") |>
  
  filter(is.finite(lat), is.finite(long), is.finite(xi_mediana), is.finite(n_years)) |>
  
  mutate(grupo_duracao = case_when(tipo_modelo == "conjunto_10min_24h" ~ "10 min–24 h",
                                   tipo_modelo == "complementar_1h_24h" ~ "1–24 h", TRUE ~ NA_character_))

candidatas_idf_sub <- candidatas_idf_sub |>
  group_by(tipo_modelo) |>
  
  mutate(q25_xi = quantile(xi_mediana, 0.25, na.rm = TRUE), q75_xi = quantile(xi_mediana, 0.75, na.rm = TRUE),
         classe_cauda = case_when(xi_mediana <= q25_xi ~ "Cauda leve", xi_mediana >= q75_xi ~ "Cauda pesada", TRUE ~
                                    "Próxima de Gumbel")) |>
  
  ungroup() |>
  
  select(-q25_xi, -q75_xi)

candidatas_idf_sub |>
  count(grupo_duracao, classe_cauda) |>
  print(n = Inf, width = Inf)

selecionar_estacoes_distantes <- function(dados, n, min_dist_km = 700) {
  if (nrow(dados) == 0) {
    return(dados)
  }
  dados_sf <- dados |>
    sf::st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)
  selecionadas <- dados_sf[0, ]
  for (i in seq_len(nrow(dados_sf))) {
    candidata <- dados_sf[i, ]
    if (nrow(selecionadas) == 0) {
      selecionadas <- candidata
    } else {
      distancias <- as.numeric(sf::st_distance(candidata, selecionadas)) / 1000
      if (all(distancias >= min_dist_km)) {
        selecionadas <- bind_rows(selecionadas, candidata)
      }
    }
    if (nrow(selecionadas) >= n) {
      break
    }
  }
  
  # Se a distância mínima impedir completar n,
  # completar pelas próximas séries mais longas
  if (nrow(selecionadas) < n) {
    restantes <- dados_sf |>
      filter(!gauge_code %in% selecionadas$gauge_code)
    faltam <-
      n -
      nrow(selecionadas)
    selecionadas <- bind_rows(selecionadas, restantes |> slice_head(n = faltam))
  }
  sf::st_drop_geometry(selecionadas)
}

candidatas_conjunto <- candidatas_idf_sub |>
  filter(tipo_modelo == "conjunto_10min_24h")

candidatas_complementar <- candidatas_idf_sub |>
  filter(tipo_modelo == "complementar_1h_24h")

N_CONJUNTO <- 250L
N_COMPLEMENTAR <- 250L

selecionar_exemplos_fit <- function(dados, n_total, min_dist_km = 700) {
  if (nrow(dados) == 0) {
    return(dados)
  }
  
  # ---------------------------------------------------------------------------
  # Primeiro candidato de cada classe:
  # sempre a estação de maior série dentro daquela classe
  # ---------------------------------------------------------------------------
  candidato_leve <- dados |>
    filter(classe_cauda == "Cauda leve") |>
    arrange(desc(n_years), xi_mediana) |>
    slice_head(n = 1)
  candidato_gumbel <- dados |>
    filter(classe_cauda == "Próxima de Gumbel") |>
    mutate(abs_xi = abs(xi_mediana)) |>
    arrange(desc(n_years), abs_xi) |>
    slice_head(n = 1) |>
    select(-abs_xi)
  candidato_pesada <- dados |>
    filter(classe_cauda == "Cauda pesada") |>
    arrange(desc(n_years), desc(xi_mediana)) |>
    slice_head(n = 1)
  candidatos_cauda <- bind_rows(candidato_leve, candidato_gumbel, candidato_pesada) |>
    distinct(gauge_code, .keep_all = TRUE) |>
    arrange(desc(n_years))
  
  # ---------------------------------------------------------------------------
  # Todas as demais candidatas, priorizando séries longas
  # ---------------------------------------------------------------------------
  restantes <- dados |>
    filter(!gauge_code %in% candidatos_cauda$gauge_code) |>
    mutate(abs_xi = abs(xi_mediana)) |>
    arrange(desc(n_years), abs_xi) |>
    select(-abs_xi)
  
  # ---------------------------------------------------------------------------
  # Pool:
  # diversidade de cauda primeiro,
  # depois demais séries mais longas
  # ---------------------------------------------------------------------------
  pool <- bind_rows(candidatos_cauda, restantes) |>
    distinct(gauge_code, .keep_all = TRUE)
  
  # ---------------------------------------------------------------------------
  # Finalmente considerar distância espacial
  # ---------------------------------------------------------------------------
  selecionar_estacoes_distantes(dados = pool, n = n_total, min_dist_km = min_dist_km)
}

sel_conjunto <- selecionar_exemplos_fit(dados = candidatas_conjunto, n_total = N_CONJUNTO, min_dist_km = 700)

sel_conjunto <- sel_conjunto[c(3,10,48),]

sel_complementar <- selecionar_exemplos_fit(dados = candidatas_complementar, n_total = N_COMPLEMENTAR,
                                            min_dist_km = 700)

sel_complementar <- sel_complementar[c(1,2,13),]

estacoes_idf_exemplo <- bind_rows(sel_complementar, sel_conjunto) |>
  
  distinct(gauge_code, tipo_modelo, .keep_all = TRUE) |>
  
  arrange(desc(n_years)) |>
  
  mutate(painel = paste0(UF, " - ", gauge_code), ordem_painel = row_number())

param_texto_idf_sub <- parametros_estruturais_subdiarios |>
  
  semi_join(estacoes_idf_exemplo |> select(gauge_code, tipo_modelo), by = c("gauge_code", "tipo_modelo")) |>
  
  filter(parametro %in% c("H", "theta")) |>
  
  mutate(nome_param = case_when(parametro == "H" & ramo == "curto" ~ "H_curto", parametro == "H" &
                                  ramo == "longo" ~ "H_longo", parametro == "theta" & ramo == "curto" ~ "theta_curto",
                                parametro == "theta" & ramo == "longo" ~ "theta_longo", TRUE ~ NA_character_)) |>
  
  filter(!is.na(nome_param)) |>
  
  select(gauge_code, tipo_modelo, nome_param, mediana) |>
  
  distinct() |>
  
  pivot_wider(names_from = nome_param, values_from = mediana)

param_texto_idf_sub <- param_texto_idf_sub |>
  
  left_join(parametros_gev_subdiarios |> filter(duracao_h == 24) |> semi_join(estacoes_idf_exemplo |> select(
    gauge_code, tipo_modelo), by = c("gauge_code", "tipo_modelo")) |> select(gauge_code, tipo_modelo,
                                                                             xi_mediana) |> distinct(), by = c("gauge_code", "tipo_modelo")) |>
  
  left_join(estacoes_idf_exemplo |> select(gauge_code, tipo_modelo, painel, n_years, classe_cauda,
                                           grupo_duracao, ordem_painel), by = c("gauge_code", "tipo_modelo")) |>
  
  arrange(ordem_painel)

param_texto_idf_sub <- param_texto_idf_sub |>
  mutate(painel = factor(painel, levels = estacoes_idf_exemplo$painel), texto_parametros = case_when(
    # -----------------------------------------------------------------------
    # Modelo conjunto: tem ramo curto e longo
    # -----------------------------------------------------------------------
    tipo_modelo == "conjunto_10min_24h" ~ paste0("ξ(24 h) = ", sprintf("%.3f", xi_mediana), "\nH<1 h = ",
                                                 sprintf("%.3f", H_curto), " | H≥1 h = ", sprintf("%.3f", H_longo), "\nθ<1 h = ", sprintf("%.3f",
                                                                                                                                          theta_curto), " | θ≥1 h = ", sprintf("%.3f", theta_longo), "\nAnos = ", n_years),
    # -----------------------------------------------------------------------
    # Modelo complementar: somente >= 1 h
    # -----------------------------------------------------------------------
    tipo_modelo == "complementar_1h_24h" ~ paste0("ξ(24 h) = ", sprintf("%.3f", xi_mediana), "\nH≥1 h = ",
                                                  sprintf("%.3f", H_longo), " | θ≥1 h = ", sprintf("%.3f", theta_longo), "\nAnos = ", n_years)))

# =============================================================================
# 10. DADOS DAS CURVAS IDF DAS ESTAÇÕES SELECIONADAS
# =============================================================================

idf_exemplos_sub <- niveis_retorno_subdiarios |>
  
  semi_join(estacoes_idf_exemplo |> select(gauge_code, tipo_modelo), by = c("gauge_code", "tipo_modelo")) |>
  
  filter(periodo_retorno_anos %in% TR_IDF_SUB) |>
  
  left_join(estacoes_idf_exemplo |> select(gauge_code, tipo_modelo, UF, painel, grupo_duracao, ordem_painel),
            by = c("gauge_code", "tipo_modelo")) |>
  
  mutate(TR = factor(periodo_retorno_anos, levels = TR_IDF_SUB), painel = factor(painel, levels =
                                                                                   estacoes_idf_exemplo$painel)) |>
  
  arrange(ordem_painel, periodo_retorno_anos, duracao_h)

estacoes_idf_exemplo |>
  select(grupo_duracao, UF, gauge_code, n_years, xi_mediana, classe_cauda) |>
  print(n = Inf, width = Inf)

idf_exemplos_sub |>
  group_by(painel, grupo_duracao) |>
  summarise(duracao_min = min(duracao_h, na.rm = TRUE), duracao_max = max(duracao_h, na.rm = TRUE), n_duracoes =
              n_distinct(duracao_h), .groups = "drop") |>
  print(n = Inf, width = Inf)

cores_TR <- c("2"   = "#D55E00", "5"   = "#E69F00", "10"  = "#009E73", "25"  = "#1B9E77", "50"  = "#1F78B4",
              "100" = "#CC79A7", "200" = "#7B3294")

plot_idf_sem_ic <- function(dados, dados_texto, TR_sel) {
  dados_plot <- dados |>
    filter(periodo_retorno_anos %in% TR_sel) |>
    mutate(TR = factor(periodo_retorno_anos, levels = TR_sel))
  ggplot(dados_plot, aes(x = duracao_h, y = mediana, group = TR, color = TR)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey40") +
    geom_text(data = dados_texto, aes(x = Inf, y = Inf, label = texto_parametros), inherit.aes = FALSE,
              hjust = 1.03, vjust = 1.05, family = font, size = fontsize / .pt, lineheight = 1.05) +
    facet_wrap(~painel, ncol = 2, scales = "free_y") +
    scale_x_log10(breaks = c(1/6, 0.25, 0.5, 1, 3, 6, 12, 24), labels = c("0.16", "0.25", "0.5", "1",
                                                                          "3", "6", "12", "24")) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
    scale_color_manual(values = cores_TR, breaks = as.character(TR_sel), name = "Tempo de retorno (anos)") +
    labs(x = "Duração (h)", y = "Intensidade (mm/h)") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "left")) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3), text = element_text(
      family = font, size = fontsize), axis.text = element_text(family = font, size = fontsize),
      axis.title = element_text(family = font, size = fontsize), strip.background = element_blank(),
      strip.text = element_text(family = font, size = fontsize), panel.border = element_rect(linewidth = 0.4),
      legend.position = "bottom", plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
}

plot_idf_com_ic <- function(dados, dados_texto, TR_sel) {
  dados_plot <- dados |>
    filter(periodo_retorno_anos %in% TR_sel) |>
    mutate(TR = factor(periodo_retorno_anos, levels = TR_sel))
  ggplot(dados_plot, aes(x = duracao_h, y = mediana, group = TR, color = TR, fill = TR)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey40") +
    geom_text(data = dados_texto, aes(x = Inf, y = Inf, label = texto_parametros), inherit.aes = FALSE,
              hjust = 1.03, vjust = 1.05, family = font, size = fontsize / .pt, lineheight = 1.05) +
    facet_wrap(~painel, ncol = 2, scales = "free_y") +
    scale_x_log10(breaks = c(1/6, 0.25, 0.5, 1, 3, 6, 12, 24), labels = c("0.16", "0.25", "0.5", "1",
                                                                          "3", "6", "12", "24")) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
    scale_color_manual(values = cores_TR, breaks = as.character(TR_sel), name = "Tempo de retorno (anos)") +
    scale_fill_manual(values = cores_TR, breaks = as.character(TR_sel), name = "Tempo de retorno (anos)") +
    labs(x = "Duração (h)", y = "Intensidade (mm/h)") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"), fill = "none") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3), text = element_text(
      family = font, size = fontsize), axis.text = element_text(family = font, size = fontsize),
      axis.title = element_text(family = font, size = fontsize), strip.background = element_blank(),
      strip.text = element_text(family = font, size = fontsize), panel.border = element_rect(linewidth = 0.4),
      legend.position = "bottom", plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
}

p_idf_sub_sem_ic <- plot_idf_sem_ic(dados = idf_exemplos_sub, dados_texto = param_texto_idf_sub,
                                    TR_sel = c(2, 5, 10, 25, 50, 100, 200))

print(p_idf_sub_sem_ic)

ggsave(file.path(DIR_RESULTADOS, "IDF_subdiarias_sem_IC.png"), p_idf_sub_sem_ic, width = 19, height = 19,
       units = "cm", dpi = 600)

p_idf_sub_ic_5_25_200 <- plot_idf_com_ic(dados = idf_exemplos_sub, dados_texto = param_texto_idf_sub,
                                         TR_sel = c(2, 10, 50, 200))

print(p_idf_sub_ic_5_25_200)

ggsave(file.path(DIR_RESULTADOS, "IDF_subdiarias_IC_5_25_200.png"), p_idf_sub_ic_5_25_200, width = 19,
       height = 19, units = "cm", dpi = 600)

p_idf_sub_ic_10_50_100 <- plot_idf_com_ic(dados = idf_exemplos_sub, dados_texto = param_texto_idf_sub,
                                          TR_sel = c(5, 25, 100))

print(p_idf_sub_ic_10_50_100)

ggsave(file.path(DIR_RESULTADOS, "IDF_subdiarias_IC_10_50_100.png"), p_idf_sub_ic_10_50_100, width = 19,
       height = 19, units = "cm", dpi = 600)

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

font <- "sans"
fontsize <- 8

# garantir um painel para cada estação
estacoes_idf_exemplo <- estacoes_idf_exemplo |>
  mutate(painel = if ("UF" %in% names(estacoes_idf_exemplo)) {
    paste0(UF, " - ", gauge_code)
  } else {
    gauge_code
  })

# dados dos parâmetros GEV apenas das 6 estações escolhidas
dados_parametros_idf <- parametros_gev_subdiarios |>
  semi_join(estacoes_idf_exemplo |> select(gauge_code), by = "gauge_code") |>
  left_join(estacoes_idf_exemplo |> select(gauge_code, painel), by = "gauge_code") |>
  select(gauge_code, painel, duracao_h, mu_mediana, mu_q025, mu_q975, sigma_mediana, sigma_q025, sigma_q975,
         xi_mediana, xi_q025, xi_q975) |>
  pivot_longer(cols = -c(gauge_code, painel, duracao_h), names_to = c("parametro", ".value"),
               names_pattern = "(mu|sigma|xi)_(mediana|q025|q975)") |>
  mutate(parametro = factor(parametro, levels = c("mu", "sigma", "xi")), painel = factor(painel,
                                                                                         levels = estacoes_idf_exemplo$painel))
library(ggh4x)
fig_parametros_6_idf <- ggplot(dados_parametros_idf, aes(x = duracao_h, y = mediana)) +
  
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "grey75", alpha = 0.45) +
  
  geom_line(color = "black", linewidth = 0.9) +
  
  geom_point(color = "black", size = 1.8) +
  
  ggh4x::facet_grid2(rows = vars(parametro), cols = vars(painel), scales = "free_y", independent = "y",
                     labeller = labeller(parametro = c(mu = "\u03BC", sigma = "\u03C3", xi = "\u03BE"))) +
  
  scale_x_log10(breaks = c(1/6, 0.25, 0.5, 1, 3, 6, 12, 24), labels = c("0.16", "0.25", "0.5", "1",
                                                                        "3", "6", "12", "24")) +
  
  labs(x = "Duração (h)", y = "Mediana posterior entre estações") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3), text = element_text(
    family = font, size = fontsize), axis.text = element_text(family = font, size = fontsize),
    axis.title = element_text(family = font, size = fontsize), strip.background = element_blank(),
    strip.text = element_text(family = font, size = fontsize), panel.border = element_rect(linewidth = 0.4),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(fig_parametros_6_idf)

ggsave(file.path(DIR_RESULTADOS, "parametros_GEV_6_estacoes_IDF.png"), fig_parametros_6_idf, width = 45,
       height = 15, units = "cm", dpi = 600)

library(dplyr)
library(tidyr)
library(sf)
library(gstat)
library(ggplot2)
library(ggspatial)
library(geobr)
library(scales)
library(grid)

font <- "sans"
fontsize <- 8

# =============================================================================
# METADADOS
# =============================================================================

meta_espacial_sub <- metadados_subdiarios |>
  transmute(gauge_code = trimws(as.character(gauge_code)), UF = as.character(UF), lat = as.numeric(lat),
            long = as.numeric(long)) |>
  filter(is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

brasil <- geobr::read_state(year = 2020) |>
  st_transform(4674)

brasil_union <- st_as_sf(st_union(st_make_valid(brasil)))

st_crs(brasil_union) <- st_crs(4674)

bbox_br <- st_bbox(brasil_union)

# =============================================================================
# ESTAÇÕES COM DIVERGÊNCIAS NO STAN
# =============================================================================

estacoes_divergencia <- diagnosticos_subdiarios |>
  transmute(gauge_code = as.character(gauge_code), tipo_modelo = as.character(tipo_modelo),
            n_divergencias = as.numeric(n_divergencias), flag_divergencia = as.logical(flag_divergencia)) |>
  
  left_join(meta_espacial_sub |> select(gauge_code, UF, lat, long), by = "gauge_code") |>
  
  filter(is.finite(lat), is.finite(long)) |>
  
  mutate(classe_divergencia = if_else(flag_divergencia, "Com divergência", "Sem divergência"),
         classe_divergencia = factor(classe_divergencia, levels = c("Sem divergência", "Com divergência")))

resumo_divergencias <- estacoes_divergencia |>
  count(classe_divergencia, name = "n_estacoes") |>
  mutate(percentual = 100 * n_estacoes / sum(n_estacoes))

print(resumo_divergencias, n = Inf, width = Inf)

estacoes_divergencia_sf <- estacoes_divergencia |>
  st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)

p_divergencias <- ggplot() +
  
  geom_sf(data = brasil, fill = "grey98", color = "grey65", linewidth = 0.20) +
  
  geom_sf(data = brasil_union, fill = NA, color = "grey30", linewidth = 0.40) +
  
  geom_sf(data = estacoes_divergencia_sf |> filter(classe_divergencia == "Sem divergência"), aes(
    color = classe_divergencia), size = 0.5, alpha = 0.60) +
  
  geom_sf(data = estacoes_divergencia_sf |> filter(classe_divergencia == "Com divergência"), aes(
    color = classe_divergencia), size = 0.5, alpha = 0.95) +
  
  scale_color_manual(name = NULL, values = c("Sem divergência" = "grey10", "Com divergência" = "#D55E00"),
                     drop = FALSE) +
  
  annotation_scale(location = "br", width_hint = 0.20, height = unit(1.5, "mm"), text_family = font,
                   text_cex = 0.6, bar_cols = c("black", "aliceblue"), line_width = 0.4) +
  
  annotation_north_arrow(location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),
                         pad_x = unit(0.30, "cm"), pad_y = unit(0.30, "cm"), which_north = "true", style = north_arrow_nautical(
                           text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))) +
  
  coord_sf(crs = st_crs(4674), xlim = c(bbox_br["xmin"], bbox_br["xmax"]), ylim = c(bbox_br["ymin"],
                                                                                    bbox_br["ymax"]), expand = FALSE) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), text = element_text(family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize), axis.title = element_text(family = font,
                                                                                            size = fontsize), legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(family = font, size = fontsize), legend.text = element_text(family = font,
                                                                                                size = fontsize), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

ggsave(filename = file.path(DIR_RESULTADOS, "mapa_estacoes_transicao_subdiario.png"), plot = p_divergencias,
       width = 15, height = 12, units = "cm", dpi = 600)

grade_pts <- st_make_grid(brasil_union, cellsize = 0.1, what = "centers", square = TRUE) |>
  st_as_sf()

grade_pts <- grade_pts[lengths(st_intersects(grade_pts, brasil_union)) > 0,]

grade_pts

interpolar_idw <- function(dados, coluna_valor, grade = grade_pts, idp = 2) {
  dados_sf <- dados |>
    filter(is.finite(.data[[coluna_valor]]), is.finite(lat), is.finite(long)) |>
    st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)
  formula_idw <- as.formula(paste0(coluna_valor, " ~ 1"))
  resultado <- gstat::idw(formula = formula_idw, locations = dados_sf, newdata = grade, idp = idp)
  resultado
}

dados_H_longo <- parametros_estruturais_subdiarios |>
  filter(parametro == "H", ramo == "longo") |>
  transmute(gauge_code, H = mediana) |>
  left_join(meta_espacial_sub |> select(gauge_code, lat, long), by = "gauge_code") |>
  filter(is.finite(H), is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

idw_H_longo <- interpolar_idw(dados = dados_H_longo, coluna_valor = "H")

idw_H_longo <- idw_H_longo |>
  rename(H_interp = var1.pred)

# ============================================================
# 8. H LONGO — CLASSES
# ============================================================

breaks_H_longo <- c(0.60, 0.75, 0.85, 0.90, 0.95, 1.00)

breaks_H_longo <- unique(breaks_H_longo)

idw_H_longo <- idw_H_longo |>
  mutate(classe_H = cut(H_interp, breaks = breaks_H_longo, include.lowest = TRUE, dig.lab = 3))

coords_H_longo <- st_coordinates(idw_H_longo)

idw_H_longo_plot <- idw_H_longo |>
  st_drop_geometry() |>
  mutate(x = coords_H_longo[, 1], y = coords_H_longo[, 2])

p_H_longo <- ggplot() +
  
  geom_raster(data = idw_H_longo_plot, aes(x = x, y = y, fill = classe_H)) +
  
  geom_sf(data = brasil, fill = NA, color = "grey45", linewidth = 0.20) +
  
  geom_sf(data = brasil_union, fill = NA, color = "grey30", linewidth = 0.40) +
  
  scale_fill_brewer(palette = "YlOrRd", name = expression(H), drop = FALSE) +
  
  annotation_scale(location = "br", width_hint = 0.20, height = unit(1.5, "mm"), text_family = font,
                   text_cex = 0.6, bar_cols = c("black", "aliceblue"), line_width = 0.4) +
  
  annotation_north_arrow(location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),
                         pad_x = unit(0.30, "cm"), pad_y = unit(0.30, "cm"), which_north = "true", style = north_arrow_nautical(
                           text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))) +
  
  coord_sf(crs = st_crs(4674), xlim = c(bbox_br["xmin"], bbox_br["xmax"]), ylim = c(bbox_br["ymin"],
                                                                                    bbox_br["ymax"]), expand = FALSE) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), text = element_text(family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize), axis.title = element_text(family = font,
                                                                                            size = fontsize), legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(family = font, size = fontsize), legend.text = element_text(family = font,
                                                                                                size = fontsize), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_H_longo)

ggsave(filename = file.path(DIR_RESULTADOS, "mapa_IDW_H_longo.png"), plot = p_H_longo, width = 15, height = 12,
       units = "cm", dpi = 600)

# ============================================================
# THETA LONGO — ESTAÇÕES
# ============================================================

dados_theta_longo <- parametros_estruturais_subdiarios |>
  filter(parametro == "theta", ramo == "longo") |>
  transmute(gauge_code, theta = mediana) |>
  left_join(meta_espacial_sub |> select(gauge_code, lat, long), by = "gauge_code") |>
  filter(is.finite(theta), is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

range(dados_theta_longo$theta, na.rm = TRUE)

quantile(dados_theta_longo$theta, probs = seq(0, 1, 0.2), na.rm = TRUE)

# ============================================================
# THETA LONGO — IDW
# ============================================================

idw_theta_longo <- interpolar_idw(dados = dados_theta_longo, coluna_valor = "theta")

idw_theta_longo <- idw_theta_longo |>
  rename(theta_interp = var1.pred)

range(idw_theta_longo$theta_interp, na.rm = TRUE)

quantile(idw_theta_longo$theta_interp, probs = seq(0, 1, 0.2), na.rm = TRUE)

breaks_theta_longo <- c(0.05, 0.30, 0.40, 0.50, 0.60, 2.00)

idw_theta_longo <- idw_theta_longo |>
  mutate(classe_theta = cut(theta_interp, breaks = breaks_theta_longo, include.lowest = TRUE, right = TRUE,
                            labels = c("0.05 – 0.30", "0.30 – 0.40", "0.40 – 0.50", "0.50 – 0.60", "0.60 – 2.00")))

coords_theta_longo <- st_coordinates(idw_theta_longo)

idw_theta_longo_plot <- idw_theta_longo |>
  st_drop_geometry() |>
  mutate(x = coords_theta_longo[, 1], y = coords_theta_longo[, 2])
p_theta_longo <- ggplot() +
  
  geom_raster(data = idw_theta_longo_plot, aes(x = x, y = y, fill = classe_theta)) +
  
  geom_sf(data = brasil, fill = NA, color = "grey45", linewidth = 0.20) +
  
  geom_sf(data = brasil_union, fill = NA, color = "grey30", linewidth = 0.40) +
  
  scale_fill_brewer(palette = "GnBu", name = expression(H), drop = FALSE) + 
  
  annotation_scale(location = "br", width_hint = 0.20, height = unit(1.5, "mm"), text_family = font,
                   text_cex = 0.6, bar_cols = c("black", "aliceblue"), line_width = 0.4) +
  
  annotation_north_arrow(location = "tr", width = unit(1, "cm"), height = unit(1, "cm"), pad_x = unit(0.30, "cm"
  ), pad_y = unit(0.30, "cm"), which_north = "true", style = north_arrow_nautical(text_family = font,
                                                                                  text_size = fontsize, fill = c("black", "aliceblue"))) +
  
  coord_sf(crs = st_crs(4674), xlim = c(bbox_br["xmin"], bbox_br["xmax"]), ylim = c(bbox_br["ymin"],
                                                                                    bbox_br["ymax"]), expand = FALSE) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), text = element_text(family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize), axis.title = element_text(family = font,
                                                                                            size = fontsize), legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(family = font, size = fontsize), legend.text = element_text(family = font,
                                                                                                size = fontsize), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_theta_longo)

ggsave(filename = file.path(DIR_RESULTADOS, "mapa_IDW_theta_longo.png"), plot = p_theta_longo, width = 14,
       height = 12, units = "cm", dpi = 600)

# ============================================================
# H CURTO — ESTAÇÕES
# ============================================================

dados_H_curto <- parametros_estruturais_subdiarios |>
  filter(parametro == "H", ramo == "curto") |>
  transmute(gauge_code, H = mediana) |>
  left_join(meta_espacial_sub |> select(gauge_code, lat, long), by = "gauge_code") |>
  filter(is.finite(H), is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

idw_H_curto <- interpolar_idw(dados = dados_H_curto, coluna_valor = "H")

idw_H_curto <- idw_H_curto |>
  rename(H_interp = var1.pred)

# ============================================================
# H CURTO — CLASSES
# ============================================================

breaks_H_curto <- c(0.63, 0.75, 0.80, 0.83, 0.86, 0.95)

idw_H_curto <- idw_H_curto |>
  mutate(classe_H = cut(H_interp, breaks = breaks_H_curto, include.lowest = TRUE, right = TRUE, labels = c(
    "0.63 – 0.75", "0.75 – 0.80", "0.80 – 0.83", "0.83 – 0.86", "0.86 – 0.95")))

coords_H_curto <- st_coordinates(idw_H_curto)

idw_H_curto_plot <- idw_H_curto |>
  st_drop_geometry() |>
  mutate(x = coords_H_curto[, 1], y = coords_H_curto[, 2])

p_H_curto <- ggplot() +
  
  # grade interpolada
  geom_raster(data = idw_H_curto_plot, aes(x = x, y = y, fill = classe_H)) +
  
  # limites estaduais
  geom_sf(data = brasil, fill = NA, color = "grey45", linewidth = 0.20) +
  
  # contorno do Brasil
  geom_sf(data = brasil_union, fill = NA, color = "grey30", linewidth = 0.40) +
  
  # legenda
  scale_fill_brewer(palette = "YlOrRd", name = expression(H), drop = FALSE) +
  
  # barra de escala
  annotation_scale(location = "br", width_hint = 0.20, height = unit(1.5, "mm"), text_family = font,
                   text_cex = 0.6, bar_cols = c("black", "aliceblue"), line_width = 0.4) +
  
  # seta norte
  annotation_north_arrow(location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),
                         pad_x = unit(0.30, "cm"), pad_y = unit(0.30, "cm"), which_north = "true", style = north_arrow_nautical(
                           text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))) +
  
  # enquadramento
  coord_sf(crs = st_crs(4674), xlim = c(bbox_br["xmin"], bbox_br["xmax"]), ylim = c(bbox_br["ymin"],
                                                                                    bbox_br["ymax"]), expand = FALSE) +
  
  # eixos
  labs(x = "Longitude", y = "Latitude") +
  
  # tema
  theme_bw() +
  
  theme(panel.grid = element_blank(), text = element_text(family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize), axis.title = element_text(family = font,
                                                                                            size = fontsize), legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(family = font, size = fontsize), legend.text = element_text(family = font,
                                                                                                size = fontsize), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_H_curto)

ggsave(filename = file.path(DIR_RESULTADOS, "mapa_IDW_H_curto.png"), plot = p_H_curto, width = 15, height = 12,
       units = "cm", dpi = 600)

dados_theta_curto <- parametros_estruturais_subdiarios |>
  filter(parametro == "theta", ramo == "curto") |>
  transmute(gauge_code, theta = mediana) |>
  left_join(meta_espacial_sub |> select(gauge_code, lat, long), by = "gauge_code") |>
  filter(is.finite(theta), is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

range(dados_theta_curto$theta, na.rm = TRUE)

quantile(dados_theta_curto$theta, probs = seq(0, 1, 0.2), na.rm = TRUE)

idw_theta_curto <- interpolar_idw(dados = dados_theta_curto, coluna_valor = "theta")

idw_theta_curto <- idw_theta_curto |>
  rename(theta_interp = var1.pred)

breaks_theta_curto <- c(0.04, 0.20, 0.27, 0.35, 0.45, 3.10)

idw_theta_curto <- idw_theta_curto |>
  mutate(classe_theta = cut(theta_interp, breaks = breaks_theta_curto, include.lowest = TRUE, right = TRUE,
                            labels = c("0.04 – 0.20", "0.20 – 0.27", "0.27 – 0.35", "0.35 – 0.45", "0.45 – 3.10")))

coords_theta_curto <- st_coordinates(idw_theta_curto)

idw_theta_curto_plot <- idw_theta_curto |>
  st_drop_geometry() |>
  mutate(x = coords_theta_curto[, 1], y = coords_theta_curto[, 2])

p_theta_curto <- ggplot() +
  
  geom_raster(data = idw_theta_curto_plot, aes(x = x, y = y, fill = classe_theta)) +
  
  geom_sf(data = brasil, fill = NA, color = "grey45", linewidth = 0.20) +
  
  geom_sf(data = brasil_union, fill = NA, color = "grey30", linewidth = 0.40) +
  
  scale_fill_brewer(palette = "GnBu", name = expression(H), drop = FALSE) +
  
  annotation_scale(location = "br", width_hint = 0.20, height = unit(1.5, "mm"), text_family = font,
                   text_cex = 0.6, bar_cols = c("black", "aliceblue"), line_width = 0.4) +
  
  annotation_north_arrow(location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),
                         pad_x = unit(0.30, "cm"), pad_y = unit(0.30, "cm"), which_north = "true", style = north_arrow_nautical(
                           text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))) +
  
  coord_sf(crs = st_crs(4674), xlim = c(bbox_br["xmin"], bbox_br["xmax"]), ylim = c(bbox_br["ymin"],
                                                                                    bbox_br["ymax"]), expand = FALSE) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), text = element_text(family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize), axis.title = element_text(family = font,
                                                                                            size = fontsize), legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(family = font, size = fontsize), legend.text = element_text(family = font,
                                                                                                size = fontsize), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_theta_curto)
ggsave(filename = file.path(DIR_RESULTADOS, "mapa_IDW_theta_curto.png"), plot = p_theta_curto, width = 15,
       height = 12, units = "cm", dpi = 600)

dados_xi_bar <- parametros_estruturais_subdiarios |>
  filter(parametro == "xi_bar") |>
  transmute(gauge_code, xi_bar = mediana) |>
  left_join(meta_espacial_sub |> select(gauge_code, lat, long), by = "gauge_code") |>
  filter(is.finite(xi_bar), is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

range(dados_xi_bar$xi_bar, na.rm = TRUE)

quantile(dados_xi_bar$xi_bar, probs = seq(0, 1, 0.2), na.rm = TRUE)

idw_xi_bar <- interpolar_idw(dados = dados_xi_bar, coluna_valor = "xi_bar")

idw_xi_bar <- idw_xi_bar |>
  rename(xi_interp = var1.pred)

# ============================================================
# XI_BAR — CLASSES
# ============================================================

idw_xi_bar <- idw_xi_bar |>
  mutate(classe_xi = cut(xi_interp, breaks = c(-Inf, -0.15, -0.05, 0.05, 0.15, Inf), labels = c("< -0.15",
                                                                                                "-0.15 – -0.05", "-0.05 – 0.05", "0.05 – 0.15", "> 0.15"), include.lowest = TRUE, right = TRUE))

coords_xi_bar <- st_coordinates(idw_xi_bar)

idw_xi_bar_plot <- idw_xi_bar |>
  st_drop_geometry() |>
  mutate(x = coords_xi_bar[, 1], y = coords_xi_bar[, 2])
p_xi_bar <- ggplot() +
  
  geom_raster(data = idw_xi_bar_plot, aes(x = x, y = y, fill = classe_xi)) +
  
  geom_sf(data = brasil, fill = NA, color = "grey45", linewidth = 0.20) +
  
  geom_sf(data = brasil_union, fill = NA, color = "grey30", linewidth = 0.40) +
  
  scale_fill_manual(name = expression(bar(xi)), values = c("< -0.15"       = "#2166AC",
                                                           "-0.15 – -0.05" = "#67A9CF", "-0.05 – 0.05"  = "#F7F7F7", "0.05 – 0.15"   = "#EF8A62",
                                                           "> 0.15"        = "#B2182B"), drop = FALSE) +
  
  annotation_scale(location = "br", width_hint = 0.20, height = unit(1.5, "mm"), text_family = font,
                   text_cex = 0.6, bar_cols = c("black", "aliceblue"), line_width = 0.4) +
  
  annotation_north_arrow(location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),
                         pad_x = unit(0.30, "cm"), pad_y = unit(0.30, "cm"), which_north = "true", style = north_arrow_nautical(
                           text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))) +
  
  coord_sf(crs = st_crs(4674), xlim = c(bbox_br["xmin"], bbox_br["xmax"]), ylim = c(bbox_br["ymin"],
                                                                                    bbox_br["ymax"]), expand = FALSE) +
  
  labs(x = "Longitude", y = "Latitude") +
  
  theme_bw() +
  
  theme(panel.grid = element_blank(), text = element_text(family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize), axis.title = element_text(family = font,
                                                                                            size = fontsize), legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(family = font, size = fontsize), legend.text = element_text(family = font,
                                                                                                size = fontsize), plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

print(p_xi_bar)

ggsave(filename = file.path(DIR_RESULTADOS, "mapa_IDW_xi_bar.png"), plot = p_xi_bar, width = 15, height = 12,
       units = "cm", dpi = 600)

# ============================================================
# 1. PASTA DOS GEOPACKAGES
# ============================================================

DIR_GPKG <- file.path(DIR_RESULTADOS, "geopackage")

dir.create(DIR_GPKG, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. FLAG FINAL DE DIAGNÓSTICO
# ============================================================

diagnosticos_subdiarios <- diagnosticos_subdiarios |>
  mutate(flag_diagnostico = flag_divergencia | flag_rhat_105)

# ============================================================
# 3. METADADOS ESPACIAIS
# ============================================================

meta_gpkg_sub <- meta_espacial_sub |>
  transmute(gauge_code = trimws(as.character(gauge_code)), UF = as.character(UF), lat = as.numeric(lat), long =
              as.numeric(long)) |>
  filter(is.finite(lat), is.finite(long)) |>
  distinct(gauge_code, .keep_all = TRUE)

# ============================================================
# 4. CAMADA — PARÂMETROS ESTRUTURAIS
# ============================================================

estruturais_gpkg <- parametros_estruturais_subdiarios |>
  transmute(gauge_code = as.character(gauge_code), tipo_modelo = as.character(tipo_modelo), parametro_modelo =
              as.character(parametro_modelo), parametro = as.character(parametro), ramo = as.character(ramo), mediana =
              as.numeric(mediana), q025 = as.numeric(q025), q975 = as.numeric(q975)) |>
  
  left_join(diagnosticos_subdiarios |> select(gauge_code, flag_diagnostico) |> distinct(gauge_code,
                                                                                        .keep_all = TRUE), by = "gauge_code") |>
  
  left_join(meta_gpkg_sub, by = "gauge_code") |>
  
  filter(is.finite(lat), is.finite(long)) |>
  
  st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)

# ============================================================
# 5. CAMADA — PARÂMETROS GEV
# ============================================================

gev_gpkg <- parametros_gev_subdiarios |>
  transmute(gauge_code = as.character(gauge_code), tipo_modelo = as.character(tipo_modelo), duracao_min =
              as.numeric(duracao_min), duracao_h = as.numeric(duracao_h), mu_mediana = as.numeric(mu_mediana), mu_q025 =
              as.numeric(mu_q025), mu_q975 = as.numeric(mu_q975), sigma_mediana = as.numeric(sigma_mediana),
            sigma_q025 = as.numeric(sigma_q025), sigma_q975 = as.numeric(sigma_q975), xi_mediana =
              as.numeric(xi_mediana), xi_q025 = as.numeric(xi_q025), xi_q975 = as.numeric(xi_q975)) |>
  
  left_join(diagnosticos_subdiarios |> select(gauge_code, flag_diagnostico) |> distinct(gauge_code,
                                                                                        .keep_all = TRUE), by = "gauge_code") |>
  
  left_join(meta_gpkg_sub, by = "gauge_code") |>
  
  filter(is.finite(lat), is.finite(long)) |>
  
  st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)

# ============================================================
# 6. CAMADA — DIAGNÓSTICOS
# ============================================================

diagnosticos_gpkg <- diagnosticos_subdiarios |>
  transmute(gauge_code = as.character(gauge_code), tipo_modelo = as.character(tipo_modelo), flag_diagnostico =
              as.logical(flag_diagnostico)) |>
  
  distinct(gauge_code, .keep_all = TRUE) |>
  
  left_join(meta_gpkg_sub, by = "gauge_code") |>
  
  filter(is.finite(lat), is.finite(long)) |>
  
  st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)

ARQ_PARAM_SUB <- file.path(DIR_GPKG, "BR_parametros_modelo_subdiario.gpkg")

if (file.exists(ARQ_PARAM_SUB)) {
  file.remove(ARQ_PARAM_SUB)
}

st_write(estruturais_gpkg, dsn = ARQ_PARAM_SUB, layer = "parametros_estruturais", quiet = TRUE)

st_write(gev_gpkg, dsn = ARQ_PARAM_SUB, layer = "parametros_gev", quiet = TRUE)

st_write(diagnosticos_gpkg, dsn = ARQ_PARAM_SUB, layer = "diagnosticos", quiet = TRUE)

# ============================================================
# 8. CAMADA — NÍVEIS DE RETORNO
# ============================================================

niveis_retorno_gpkg <- niveis_retorno_subdiarios |>
  transmute(gauge_code = as.character(gauge_code), tipo_modelo = as.character(tipo_modelo), duracao_min =
              as.numeric(duracao_min), duracao_h = as.numeric(duracao_h), periodo_retorno_anos = as.numeric(
                periodo_retorno_anos), mediana = as.numeric(mediana), q025 = as.numeric(q025), q975 = as.numeric(q975)
  ) |>
  
  left_join(diagnosticos_subdiarios |> select(gauge_code, flag_diagnostico) |> distinct(gauge_code,
                                                                                        .keep_all = TRUE), by = "gauge_code") |>
  
  left_join(meta_gpkg_sub, by = "gauge_code") |>
  
  filter(is.finite(lat), is.finite(long)) |>
  
  st_as_sf(coords = c("long", "lat"), crs = 4674, remove = FALSE)

# ============================================================
# 9. GEOPACKAGE — PERÍODOS DE RETORNO
# ============================================================

ARQ_RL_SUB <- file.path(DIR_GPKG, "BR_periodos_retorno_subdiario.gpkg")

if (file.exists(ARQ_RL_SUB)) {
  file.remove(ARQ_RL_SUB)
}

st_write(niveis_retorno_gpkg, dsn = ARQ_RL_SUB, layer = "niveis_retorno", quiet = TRUE)
