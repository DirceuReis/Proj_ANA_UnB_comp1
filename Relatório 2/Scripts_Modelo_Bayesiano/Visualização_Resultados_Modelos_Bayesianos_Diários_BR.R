library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)

# =============================================================================
# 1. CAMINHO DOS DIAGNÓSTICOS
# =============================================================================

DIR_DIAG <- "C:/Larissa/diagnosticos"


# =============================================================================
# 2. LOCALIZAR TODOS OS diagnostico_estacao.csv
#
# Estrutura esperada:
#
# C:/Larissa/diagnosticos/AM/estacoes/63000/diagnostico_estacao.csv
# C:/Larissa/diagnosticos/...
# =============================================================================

arquivos_diag <- list.files(
  DIR_DIAG,
  pattern = "^diagnostico_estacao\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "Arquivos de diagnóstico encontrados:",
  length(arquivos_diag),
  "\n"
)


# =============================================================================
# 3. LER TODOS OS DIAGNÓSTICOS
#
# gauge_code é forçado para character para não perder eventuais zeros iniciais.
# =============================================================================

ler_diagnostico <- function(arquivo) {
  
  x <- read_csv(
    arquivo,
    col_types = cols(
      .default = col_guess(),
      gauge_code = col_character(),
      UF = col_character()
    ),
    show_col_types = FALSE
  )
  
  x %>%
    mutate(
      arquivo_diagnostico = arquivo
    )
}


diagnosticos_brasil <- map_dfr(
  arquivos_diag,
  ler_diagnostico
)


# =============================================================================
# 4. GARANTIR UMA LINHA POR ESTAÇÃO
# =============================================================================

duplicadas <- diagnosticos_brasil %>%
  count(
    UF,
    gauge_code
  ) %>%
  filter(
    n > 1
  )

if (nrow(duplicadas) > 0) {
  
  print(
    duplicadas,
    n = Inf
  )
  
  stop(
    "Há estações com mais de um diagnostico_estacao.csv."
  )
}


# =============================================================================
# 5. CRIAR FLAGS
# =============================================================================

diagnosticos_brasil <- diagnosticos_brasil %>%
  mutate(
    
    # -------------------------------------------------------------------------
    # Divergências
    # -------------------------------------------------------------------------
    
    flag_divergencia =
      !is.na(n_divergences_retained) &
      n_divergences_retained > 0,
    
    # -------------------------------------------------------------------------
    # Rhat > 1.01
    #
    # É mantido como alerta para o diagnóstico e para o relatório,
    # mas NÃO será a flag espacial principal.
    # -------------------------------------------------------------------------
    
    alerta_rhat_101 =
      is.finite(rhat_max) &
      rhat_max > 1.01,
    
    # -------------------------------------------------------------------------
    # Rhat > 1.05
    #
    # Este será o critério usado como FLAG no produto espacial.
    # -------------------------------------------------------------------------
    
    flag_rhat_105 =
      is.finite(rhat_max) &
      rhat_max > 1.05,
    
    # -------------------------------------------------------------------------
    # Flag geral do diagnóstico
    #
    # 1 = possui pelo menos:
    #     divergência OU Rhat > 1.05
    # -------------------------------------------------------------------------
    
    flag_diagnostico =
      flag_divergencia |
      flag_rhat_105
  )


# =============================================================================
# 6. RESUMO GERAL — BRASIL
# =============================================================================

resumo_diagnostico_brasil <- diagnosticos_brasil %>%
  summarise(
    
    # -------------------------------------------------------------------------
    # Total
    # -------------------------------------------------------------------------
    
    n_estacoes =
      n(),
    
    # -------------------------------------------------------------------------
    # Divergências
    # -------------------------------------------------------------------------
    
    n_estacoes_com_divergencia =
      sum(
        flag_divergencia,
        na.rm = TRUE
      ),
    
    pct_estacoes_com_divergencia =
      100 *
      n_estacoes_com_divergencia /
      n_estacoes,
    
    n_divergencias_retidas_total =
      sum(
        n_divergences_retained,
        na.rm = TRUE
      ),
    
    max_divergencias_em_estacao =
      max(
        n_divergences_retained,
        na.rm = TRUE
      ),
    
    # -------------------------------------------------------------------------
    # Rhat > 1.01
    # -------------------------------------------------------------------------
    
    n_estacoes_rhat_gt_101 =
      sum(
        alerta_rhat_101,
        na.rm = TRUE
      ),
    
    pct_estacoes_rhat_gt_101 =
      100 *
      n_estacoes_rhat_gt_101 /
      n_estacoes,
    
    # -------------------------------------------------------------------------
    # Rhat > 1.05
    # -------------------------------------------------------------------------
    
    n_estacoes_rhat_gt_105 =
      sum(
        flag_rhat_105,
        na.rm = TRUE
      ),
    
    pct_estacoes_rhat_gt_105 =
      100 *
      n_estacoes_rhat_gt_105 /
      n_estacoes,
    
    # -------------------------------------------------------------------------
    # Divergência + Rhat > 1.05 simultaneamente
    # -------------------------------------------------------------------------
    
    n_estacoes_div_e_rhat_gt_105 =
      sum(
        flag_divergencia &
          flag_rhat_105,
        na.rm = TRUE
      ),
    
    # -------------------------------------------------------------------------
    # Qualquer flag considerada relevante para o produto
    # -------------------------------------------------------------------------
    
    n_estacoes_flag_diagnostico =
      sum(
        flag_diagnostico,
        na.rm = TRUE
      ),
    
    pct_estacoes_flag_diagnostico =
      100 *
      n_estacoes_flag_diagnostico /
      n_estacoes,
    
    # -------------------------------------------------------------------------
    # Maior Rhat observado
    # -------------------------------------------------------------------------
    
    rhat_max_brasil =
      max(
        rhat_max,
        na.rm = TRUE
      )
  )


cat(
  "\n============================================================\n",
  "RESUMO DO DIAGNÓSTICO — BRASIL\n",
  "============================================================\n"
)

print(
  resumo_diagnostico_brasil,
  width = Inf
)


# =============================================================================
# 7. RESUMO POR UF
# =============================================================================

resumo_diagnostico_uf <- diagnosticos_brasil %>%
  group_by(
    UF
  ) %>%
  summarise(
    
    n_estacoes =
      n(),
    
    n_com_divergencia =
      sum(
        flag_divergencia,
        na.rm = TRUE
      ),
    
    pct_com_divergencia =
      100 *
      n_com_divergencia /
      n_estacoes,
    
    divergencias_retidas_total =
      sum(
        n_divergences_retained,
        na.rm = TRUE
      ),
    
    n_rhat_gt_101 =
      sum(
        alerta_rhat_101,
        na.rm = TRUE
      ),
    
    pct_rhat_gt_101 =
      100 *
      n_rhat_gt_101 /
      n_estacoes,
    
    n_rhat_gt_105 =
      sum(
        flag_rhat_105,
        na.rm = TRUE
      ),
    
    pct_rhat_gt_105 =
      100 *
      n_rhat_gt_105 /
      n_estacoes,
    
    n_flag_diagnostico =
      sum(
        flag_diagnostico,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(pct_com_divergencia)
  )


print(
  resumo_diagnostico_uf,
  n = Inf,
  width = Inf
)


# =============================================================================
# 8. ESTAÇÕES COM DIVERGÊNCIA
# =============================================================================

estacoes_com_divergencia <- diagnosticos_brasil %>%
  filter(
    flag_divergencia
  ) %>%
  arrange(
    desc(n_divergences_retained),
    desc(rhat_max)
  ) %>%
  select(
    UF,
    gauge_code,
    thin,
    n_draws_retained,
    n_divergences_retained,
    pct_divergences_retained,
    chains_with_divergence,
    rhat_max,
    parameter_rhat_max,
    flag_rhat_105
  )


cat(
  "\n============================================================\n",
  "ESTAÇÕES COM DIVERGÊNCIA\n",
  "============================================================\n"
)

print(
  estacoes_com_divergencia,
  n = Inf,
  width = Inf
)


# =============================================================================
# 9. ESTAÇÕES COM RHAT > 1.05
# =============================================================================

estacoes_rhat_105 <- diagnosticos_brasil %>%
  filter(
    flag_rhat_105
  ) %>%
  arrange(
    desc(rhat_max)
  ) %>%
  select(
    UF,
    gauge_code,
    rhat_max,
    parameter_rhat_max,
    n_parameters_rhat_gt_101,
    n_parameters_rhat_gt_105,
    flag_divergencia
  )


cat(
  "\n============================================================\n",
  "ESTAÇÕES COM RHAT > 1.05\n",
  "============================================================\n"
)

print(
  estacoes_rhat_105,
  n = Inf,
  width = Inf
)


# =============================================================================
# 10. TABELA DE FLAGS PARA O GEOPACKAGE
#
# Esta é a tabela que você fará left_join() com seu objeto espacial.
# =============================================================================

flags_diagnostico <- diagnosticos_brasil %>%
  transmute(
    
    gauge_code =
      as.character(
        gauge_code
      ),
    
    UF,
    
    n_divergencias =
      n_divergences_retained,
    
    pct_divergencias =
      pct_divergences_retained,
    
    rhat_max,
    
    flag_divergencia =
      as.integer(
        flag_divergencia
      ),
    
    flag_rhat_105 =
      as.integer(
        flag_rhat_105
      ),
    
    flag_diagnostico =
      as.integer(
        flag_diagnostico
      )
  )


print(
  flags_diagnostico,
  n = 20,
  width = Inf
)


# =============================================================================
# 11. SALVAR PRODUTOS DO DIAGNÓSTICO
# =============================================================================

write_excel_csv(
  diagnosticos_brasil,
  file.path(
    DIR_DIAG,
    "diagnostico_modelo_Brasil.csv"
  )
)

saveRDS(
  diagnosticos_brasil,
  file.path(
    DIR_DIAG,
    "diagnostico_modelo_Brasil.rds"
  )
)

write_excel_csv(
  resumo_diagnostico_uf,
  file.path(
    DIR_DIAG,
    "resumo_diagnostico_por_UF.csv"
  )
)

write_excel_csv(
  flags_diagnostico,
  file.path(
    DIR_DIAG,
    "flags_diagnostico_modelo.csv"
  )
)

View(flags_diagnostico)

resumo_flags <- flags_diagnostico %>%
  summarise(
    n_estacoes = n(),
    
    n_flag_diagnostico =
      sum(flag_diagnostico == 1, na.rm = TRUE),
    
    pct_flag_diagnostico =
      100 * n_flag_diagnostico / n_estacoes,
    
    n_flag_divergencia =
      sum(flag_divergencia == 1, na.rm = TRUE),
    
    pct_flag_divergencia =
      100 * n_flag_divergencia / n_estacoes,
    
    n_flag_rhat_105 =
      sum(flag_rhat_105 == 1, na.rm = TRUE),
    
    pct_flag_rhat_105 =
      100 * n_flag_rhat_105 / n_estacoes,
    
    n_divergencia_e_rhat =
      sum(
        flag_divergencia == 1 &
          flag_rhat_105 == 1,
        na.rm = TRUE
      )
  )

print(
  resumo_flags,
  width = Inf
)

library(dplyr)
library(readr)
library(purrr)
library(sf)
library(ggplot2)
library(geobr)

# =============================================================================
# 1. LER OS METADADOS DE TODAS AS ESTAÇÕES
# =============================================================================

DIR_DIAG <- "C:/Larissa/diagnosticos"

arquivos_meta <- list.files(
  path = DIR_DIAG,
  pattern = "^metadados_estacao\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Arquivos de metadados encontrados:", length(arquivos_meta), "\n")

metadados_todos <- purrr::map_dfr(
  arquivos_meta,
  function(f) {
    tryCatch(
      {
        x <- readr::read_csv(f, show_col_types = FALSE)
        
        x %>%
          transmute(
            UF = as.character(UF),
            gauge_code = as.character(gauge_code),
            estado = dplyr::coalesce(as.character(estado), as.character(UF)),
            nome = dplyr::if_else("nome" %in% names(.), as.character(nome), NA_character_),
            lat = dplyr::coalesce(
              suppressWarnings(as.numeric(if ("lat" %in% names(.)) lat else NA)),
              suppressWarnings(as.numeric(if ("latitude" %in% names(.)) latitude else NA))
            ),
            long = dplyr::coalesce(
              suppressWarnings(as.numeric(if ("long" %in% names(.)) long else NA)),
              suppressWarnings(as.numeric(if ("longitude" %in% names(.)) longitude else NA))
            )
          ) %>%
          slice(1)
      },
      error = function(e) NULL
    )
  }
) %>%
  distinct(gauge_code, .keep_all = TRUE)

cat("Estações com metadados lidos:", nrow(metadados_todos), "\n")

# =============================================================================
# 2. JUNTAR COM AS FLAGS DO DIAGNÓSTICO
# =============================================================================
# Aqui assumo que você já tem um objeto flags_diagnostico
# com pelo menos:
# gauge_code, n_divergencias, flag_divergencia
#
# Se ainda não tiver flag_divergencia pronta:
# mutate(flag_divergencia = n_divergencias > 0)

estacoes_mapa <- metadados_todos %>%
  left_join(
    flags_diagnostico %>%
      transmute(
        gauge_code = as.character(gauge_code),
        n_divergencias = as.numeric(n_divergencias),
        flag_divergencia = as.logical(flag_divergencia)
      ),
    by = "gauge_code"
  ) %>%
  mutate(
    n_divergencias = coalesce(n_divergencias, 0),
    flag_divergencia = coalesce(flag_divergencia, FALSE)
  )

# conferir
estacoes_mapa %>%
  summarise(
    n_total = n(),
    n_div = sum(flag_divergencia, na.rm = TRUE),
    pct_div = 100 * mean(flag_divergencia, na.rm = TRUE)
  ) %>%
  print()

# =============================================================================
# 3. FILTRAR COORDENADAS VÁLIDAS
# =============================================================================

estacoes_mapa_validas <- estacoes_mapa %>%
  filter(
    is.finite(lat),
    is.finite(long)
  )

cat("Estações com coordenadas válidas:", nrow(estacoes_mapa_validas), "\n")

# =============================================================================
# 4. TRANSFORMAR EM SF
# =============================================================================

estacoes_sf <- estacoes_mapa_validas %>%
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4326,
    remove = FALSE
  )

estacoes_div_sf <- estacoes_sf %>%
  filter(flag_divergencia)

# =============================================================================
# 5. MAPA DO BRASIL
# =============================================================================

brasil <- geobr::read_state(
  year = 2020
) |>
  st_transform(4674)

brasil_union <- st_as_sf(
  st_union(
    st_make_valid(brasil)
  )
)

st_crs(brasil_union) <- st_crs(4674)

bbox_br <- st_bbox(brasil_union)

# =============================================================================
# 2. MAPA
# =============================================================================

library(ggplot2)
library(grid)

# ============================================================
# Tema
# ============================================================

library(ggplot2)
library(ggspatial)
library(sf)
library(grid)

font <- "sans"
fontsize <- 8
textsize <- fontsize / .pt

theme_set(
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      legend.background = element_blank(),
      legend.key.height = unit(fontsize, "pt"),
      strip.background = element_blank(),
      text = element_text(
        family = font,
        size = fontsize
      ),
      panel.spacing.y = unit(0.5, "pt")
    )
)

p_mapa <- ggplot() +
  
  # Estados
  geom_sf(
    data = brasil,
    fill = "grey98",
    color = "grey50",
    linewidth = 0.25
  ) +
  
  # Contorno externo
  geom_sf(
    data = brasil_union,
    fill = NA,
    color = "grey25",
    linewidth = 0.45
  ) +
  
  # Todas as estações
  geom_sf(
    data = estacoes_sf,
    aes(color = "Todas as estações"),
    size = 0.10,
    alpha = 0.70
  ) +
  
  # Estações com divergência
  geom_sf(
    data = estacoes_div_sf,
    aes(color = "Com divergência"),
    size = 0.35,
    alpha = 0.90
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Todas as estações" = "grey20",
      "Com divergência" = "darkorange2"
    ),
    breaks = c(
      "Todas as estações",
      "Com divergência"
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = c(1.8, 2.2),
        alpha = 1
      )
    )
  ) +
  
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit(1.5, "mm"),
    text_family = font,
    text_cex = 0.6,
    bar_cols = c("black", "aliceblue"),
    line_width = 0.4
  ) +
  
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
  
  coord_sf(
    crs = st_crs(4674),
    xlim = c(
      bbox_br["xmin"],
      bbox_br["xmax"]
    ),
    ylim = c(
      bbox_br["ymin"],
      bbox_br["ymax"]
    ),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text = element_text(
      family = font,
      size = fontsize
    ),
    
    axis.ticks = element_line(
      linewidth = 0.3
    ),
    
    axis.title = element_text(
      family = font,
      size = fontsize
    ),
    
    text = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.position = "bottom",
    legend.justification = "center",
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )
print(p_mapa)

ggsave(
  "C:/Larissa/mapa_estacoes_diagnostico_BR.png",
  p_mapa,
  width = 15,
  height = 15,
  units = "cm",
  dpi = 600
)


# =============================================================================
# SALVAR DADOS PARA O MAPA EM OUTRO COMPUTADOR
# =============================================================================

mapa_estacoes_csv <- estacoes_mapa %>%
  transmute(
    UF = as.character(UF),
    gauge_code = as.character(gauge_code),
    lat = as.numeric(lat),
    long = as.numeric(long),
    flag_divergencia = as.integer(flag_divergencia),
    n_divergencias = as.integer(n_divergencias)
  ) %>%
  filter(
    is.finite(lat),
    is.finite(long)
  )

readr::write_excel_csv(
  mapa_estacoes_csv,
  "C:/Larissa/mapa_estacoes_diagnostico_BR.csv"
)

cat(
  "Estações salvas:", nrow(mapa_estacoes_csv), "\n",
  "Com divergência:", sum(mapa_estacoes_csv$flag_divergencia == 1, na.rm = TRUE), "\n"
)


library(dplyr)
library(readr)
library(purrr)

DIR_DIAG <- "C:/Larissa/diagnosticos"

arquivos_parametros <- list.files(
  DIR_DIAG,
  pattern = "^parametros_posteriores\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "Arquivos encontrados:",
  length(arquivos_parametros),
  "\n"
)

parametros_brasil <- map_dfr(
  arquivos_parametros,
  function(f) {
    
    tryCatch(
      read_csv(
        f,
        show_col_types = FALSE
      ),
      error = function(e) NULL
    )
  }
)

cat(
  "Estações:",
  n_distinct(parametros_brasil$gauge_code),
  "\n"
)

cat(
  "Estações salvas:", nrow(mapa_estacoes_csv), "\n",
  "Com divergência:", sum(mapa_estacoes_csv$flag_divergencia == 1), "\n"
)

parametros_estruturais_br <- parametros_brasil %>%
  filter(
    parameter %in% c(
      "alpha",
      "beta",
      "H",
      "xi_bar",
      "tau2_mu",
      "tau2_sigma",
      "tau2_xi"
    )
  )

resumo_estruturais_br <- parametros_estruturais_br %>%
  group_by(parameter) %>%
  summarise(
    n_estacoes = n(),
    
    minimo =
      min(median, na.rm = TRUE),
    
    q05 =
      quantile(
        median,
        0.05,
        na.rm = TRUE
      ),
    
    q25 =
      quantile(
        median,
        0.25,
        na.rm = TRUE
      ),
    
    mediana =
      median(
        median,
        na.rm = TRUE
      ),
    
    q75 =
      quantile(
        median,
        0.75,
        na.rm = TRUE
      ),
    
    q95 =
      quantile(
        median,
        0.95,
        na.rm = TRUE
      ),
    
    maximo =
      max(
        median,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  resumo_estruturais_br,
  n = Inf,
  width = Inf
)

library(ggplot2)

dados_H <- parametros_brasil %>%
  filter(
    parameter == "H"
  )

font <- "sans"
fontsize <- 8

dados_H <- parametros_brasil %>%
  filter(
    parameter == "H"
  )

p_H <- ggplot(
  dados_H,
  aes(x = median)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    boundary = 0,
    fill = "grey80",
    color = "grey30",
    linewidth = 0.25
  ) +
  geom_density(
    linewidth = 0.6
  ) +
  labs(
    x = "Mediana posterior de H",
    y = "Densidade"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
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
    panel.border = element_rect(
      linewidth = 0.4
    ),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

print(p_H)

ggsave(
  "C:/Larissa/distribuicao_H_BR.png",
  p_H,
  width = 8,
  height = 6,
  units = "cm",
  dpi = 600
)

library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)
library(gstat)
library(grid)
library(classInt)

# ============================================================
# TEMA
# ============================================================

font <- "sans"
fontsize <- 8
textsize <- fontsize / .pt

theme_set(
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      legend.background = element_blank(),
      legend.key.height = unit(fontsize, "pt"),
      strip.background = element_blank(),
      text = element_text(
        family = font,
        size = fontsize
      ),
      panel.spacing.y = unit(0.5, "pt")
    )
)

# ============================================================
# 1. DADOS DAS ESTAÇÕES COM H
# ============================================================

dados_H <- parametros_brasil %>%
  filter(parameter == "H") %>%
  transmute(
    gauge_code = as.character(gauge_code),
    H = as.numeric(median)
  ) %>%
  left_join(
    metadados_todos %>%
      transmute(
        gauge_code = as.character(gauge_code),
        lat = as.numeric(lat),
        long = as.numeric(long)
      ) %>%
      distinct(gauge_code, .keep_all = TRUE),
    by = "gauge_code"
  ) %>%
  filter(
    is.finite(lat),
    is.finite(long),
    is.finite(H)
  )

estacoes_H_sf <- st_as_sf(
  dados_H,
  coords = c("long", "lat"),
  crs = 4674,
  remove = FALSE
)

# ============================================================
# 2. DISSOLVER O BRASIL EM UM ÚNICO POLÍGONO
# ============================================================

brasil_union <- st_as_sf(
  st_union(
    st_make_valid(brasil)
  )
)

grade_pts <- st_make_grid(
  brasil_union,
  cellsize = 0.1,
  what = "centers",
  square = TRUE
) %>%
  st_as_sf()

# manter apenas pontos dentro do Brasil
grade_pts <- grade_pts[
  lengths(
    st_intersects(
      grade_pts,
      brasil_union
    )
  ) > 0,
]

# ============================================================
# Conferências
# ============================================================

print(grade_pts)
print(st_geometry(grade_pts))
print(attr(grade_pts, "sf_column"))

# ============================================================
# 4. INTERPOLAÇÃO IDW
# ============================================================
# Observação:
# - idp = 2 é uma escolha comum
# - se quiser uma superfície mais suave, depois podemos testar outros valores

idw_res <- gstat::idw(
  formula = H ~ 1,
  locations = estacoes_H_sf,
  newdata = grade_pts,
  idp = 2.0
)

idw_res <- idw_res %>%
  rename(H_interp = var1.pred)

# ============================================================
# 5. CLASSIFICAR EM QUANTIS
# ============================================================

n_classes <- 4

quebras <- classInt::classIntervals(
  idw_res$H_interp,
  n = n_classes,
  style = "quantile"
)$brks

# garantir unicidade, caso haja empates em quantis
quebras <- unique(quebras)

quebras <- c(
  0,
  0.60,
  0.64,
  0.68,
  0.72,
  1
)

idw_res <- idw_res %>%
  mutate(
    classe_H = cut(
      H_interp,
      breaks = quebras,
      include.lowest = TRUE,
      dig.lab = 4
    )
  )

coords_idw <- st_coordinates(idw_res)

idw_plot <- idw_res %>%
  st_drop_geometry() %>%
  mutate(
    x = coords_idw[, 1],
    y = coords_idw[, 2]
  )

# conferir
summary(idw_plot$H_interp)
table(idw_plot$classe_H, useNA = "ifany")
font <- "sans"
fontsize <- 8

p_H_grid <- ggplot() +
  
  # ----------------------------------------------------------
# Superfície interpolada em células de 0.1°
# ----------------------------------------------------------

geom_raster(
  data = idw_plot,
  aes(
    x = x,
    y = y,
    fill = classe_H
  )
) +
  
  # ----------------------------------------------------------
# Limites estaduais
# ----------------------------------------------------------

geom_sf(
  data = brasil,
  fill = NA,
  color = "grey45",
  linewidth = 0.20
) +
  
  # ----------------------------------------------------------
# Contorno externo do Brasil
# ----------------------------------------------------------

geom_sf(
  data = brasil_union,
  fill = NA,
  color = "grey30",
  linewidth = 0.40
) +
  
  # ----------------------------------------------------------
# Classes de H
# ----------------------------------------------------------

scale_fill_brewer(
  palette = "YlOrRd",
  name = expression(H),
  drop = FALSE
) +
  
  # ----------------------------------------------------------
# Escala
# ----------------------------------------------------------

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
  
  # ----------------------------------------------------------
# Norte
# ----------------------------------------------------------

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
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    # 
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    
    text = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.position = "bottom",
    legend.justification = "center",
    
    legend.title = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.text = element_text(
      family = font,
      size = fontsize
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

print(p_H_grid)

# ============================================================
# 7. SALVAR
# ============================================================

ggsave(
  "C:/Larissa/mapa_H_interpolado_quantis.png",
  plot = p_H_grid,
  width = 15,
  height = 15,
  units = "cm",
  dpi = 600
)

library(dplyr)
library(tidyr)
library(readr)

# =============================================================================
# 1. PARÂMETROS ESTRUTURAIS QUE QUEREMOS EXPORTAR
#
# theta foi deliberadamente excluído da apresentação principal.
# =============================================================================

parametros_estruturais_exportar <- c(
  "H",
  "alpha",
  "beta",
  "xi_bar",
  "tau2_mu",
  "tau2_sigma",
  "tau2_xi"
)


# =============================================================================
# 2. UMA LINHA POR ESTAÇÃO
#
# Mantemos:
# - mediana posterior
# - limite inferior de 95%
# - limite superior de 95%
# =============================================================================

estruturais_mapa <- parametros_brasil %>%
  filter(
    parameter %in% parametros_estruturais_exportar
  ) %>%
  select(
    UF,
    gauge_code,
    parameter,
    median,
    q025,
    q975
  ) %>%
  mutate(
    gauge_code = as.character(gauge_code)
  ) %>%
  pivot_wider(
    names_from = parameter,
    values_from = c(
      median,
      q025,
      q975
    ),
    names_glue = "{parameter}_{.value}"
  )


# =============================================================================
# 3. ADICIONAR COORDENADAS
#
# metadados_todos foi o objeto que construímos a partir dos
# metadados_estacao.csv
# =============================================================================

estruturais_mapa <- estruturais_mapa %>%
  left_join(
    metadados_todos %>%
      transmute(
        gauge_code = as.character(gauge_code),
        lat = as.numeric(lat),
        long = as.numeric(long)
      ) %>%
      distinct(
        gauge_code,
        .keep_all = TRUE
      ),
    by = "gauge_code"
  ) %>%
  relocate(
    UF,
    gauge_code,
    lat,
    long
  )


# =============================================================================
# 4. CONFERÊNCIAS
# =============================================================================

cat(
  "Estações exportadas:",
  nrow(estruturais_mapa),
  "\n"
)

cat(
  "Estações sem coordenadas:",
  sum(
    !is.finite(estruturais_mapa$lat) |
      !is.finite(estruturais_mapa$long)
  ),
  "\n"
)

print(
  names(estruturais_mapa)
)


# =============================================================================
# 5. SALVAR CSV
# =============================================================================

readr::write_excel_csv(
  estruturais_mapa,
  "C:/Larissa/parametros_estruturais_BR_mapa.csv"
)

library(dplyr)
library(ggplot2)

dados_H <- parametros_brasil %>%
  filter(parameter == "H")

H_med <- median(dados_H$median, na.rm = TRUE)

H_q <- quantile(
  dados_H$median,
  probs = c(0.05, 0.95),
  na.rm = TRUE
)

p_H <- ggplot(
  dados_H,
  aes(x = median)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 45,
    fill = "grey80",
    color = "white"
  ) +
  geom_density(
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = H_med,
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = H_q,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = H_med,
    y = Inf,
    label = paste0("Mediana = ", round(H_med, 3)),
    vjust = 1.6,
    hjust = -0.1,
    size = 3.4
  ) +
  labs(
    x = "Mediana posterior de H",
    y = "Densidade"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_H)

font <- "sans"
fontsize <- 8

dados_xi <- parametros_brasil %>%
  filter(parameter == "xi_bar")

xi_med <- median(
  dados_xi$median,
  na.rm = TRUE
)

xi_q <- quantile(
  dados_xi$median,
  probs = c(0.05, 0.95),
  na.rm = TRUE
)

p_xi <- ggplot(
  dados_xi,
  aes(x = median)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 45,
    fill = "grey80",
    color = "grey30",
    linewidth = 0.25
  ) +
  geom_density(
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    x = expression("Mediana posterior de " * xi[bar]),
    y = "Densidade"
  ) +
  theme_bw() +
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
    
    panel.border = element_rect(
      linewidth = 0.4
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

print(p_xi)

ggsave(
  "C:/Larissa/distribuicao_xi_bar_BR.png",
  plot = p_xi,
  width = 8,
  height = 6,
  units = "cm",
  dpi = 600
)

library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)
library(gstat)
library(grid)

font <- "sans"
fontsize <- 8

# ============================================================
# 1. XI MÉDIO POR ESTAÇÃO
# ============================================================

dados_xi_mapa <- parametros_brasil %>%
  filter(parameter == "xi_bar") %>%
  transmute(
    gauge_code = as.character(gauge_code),
    xi_bar = as.numeric(median)
  ) %>%
  left_join(
    metadados_todos %>%
      transmute(
        gauge_code = as.character(gauge_code),
        lat = as.numeric(lat),
        long = as.numeric(long)
      ) %>%
      distinct(
        gauge_code,
        .keep_all = TRUE
      ),
    by = "gauge_code"
  ) %>%
  filter(
    is.finite(lat),
    is.finite(long),
    is.finite(xi_bar)
  )


# ============================================================
# 2. CONVERTER ESTAÇÕES PARA SF
# ============================================================

estacoes_xi_sf <- st_as_sf(
  dados_xi_mapa,
  coords = c("long", "lat"),
  crs = 4674,
  remove = FALSE
)


# ============================================================
# 3. POLÍGONO ÚNICO DO BRASIL
# ============================================================

brasil_union <- st_as_sf(
  st_union(
    st_make_valid(brasil)
  )
)


# ============================================================
# 4. GRADE REGULAR DE 0.1°
# ============================================================

grade_pts <- st_make_grid(
  brasil_union,
  cellsize = 0.1,
  what = "centers",
  square = TRUE
) %>%
  st_as_sf()

grade_pts <- grade_pts[
  lengths(
    st_intersects(
      grade_pts,
      brasil_union
    )
  ) > 0,
]


# ============================================================
# 5. INTERPOLAÇÃO IDW
# ============================================================

idw_xi <- gstat::idw(
  formula = xi_bar ~ 1,
  locations = estacoes_xi_sf,
  newdata = grade_pts,
  idp = 2
)

idw_xi <- idw_xi %>%
  rename(
    xi_interp = var1.pred
  )


# ============================================================
# 6. CLASSES
# ============================================================

idw_xi <- idw_xi %>%
  mutate(
    classe_xi = cut(
      xi_interp,
      breaks = c(
        -Inf,
        -0.15,
        -0.05,
        0.05,
        0.15,
        Inf
      ),
      labels = c(
        "< -0.15",
        "-0.15 – -0.05",
        "-0.05 – 0.05",
        "0.05 – 0.15",
        "> 0.15"
      ),
      include.lowest = TRUE,
      right = TRUE
    )
  )


# ============================================================
# 7. EXTRAIR COORDENADAS PARA O RASTER
# ============================================================

coords_xi <- st_coordinates(idw_xi)

idw_xi_plot <- idw_xi %>%
  st_drop_geometry() %>%
  mutate(
    x = coords_xi[, 1],
    y = coords_xi[, 2]
  )


# ============================================================
# 8. MAPA
# ============================================================

p_xi_mapa <- ggplot() +
  
  # grade interpolada
  geom_raster(
    data = idw_xi_plot,
    aes(
      x = x,
      y = y,
      fill = classe_xi
    )
  ) +
  
  # limites estaduais
  geom_sf(
    data = brasil,
    fill = NA,
    color = "grey45",
    linewidth = 0.20
  ) +
  
  # contorno do Brasil
  geom_sf(
    data = brasil_union,
    fill = NA,
    color = "grey30",
    linewidth = 0.40
  ) +
  
  # escala divergente
  scale_fill_manual(
    name = expression(bar(xi)),
    values = c(
      "< -0.15"       = "#2166AC",
      "-0.15 – -0.05" = "#67A9CF",
      "-0.05 – 0.05"  = "#F7F7F7",
      "0.05 – 0.15"   = "#EF8A62",
      "> 0.15"        = "#B2182B"
    ),
    drop = FALSE
  ) +
  
  # barra de escala
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
  
  # seta norte
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
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    
    text = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.position = "bottom",
    legend.justification = "center",
    
    legend.title = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.text = element_text(
      family = font,
      size = fontsize
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

print(p_xi_mapa)


# ============================================================
# 9. SALVAR
# ============================================================

ggsave(
  "C:/Larissa/mapa_xi_bar_interpolado_BR.png",
  plot = p_xi_mapa,
  width = 15,
  height = 15,
  units = "cm",
  dpi = 600
)

dados_outros <- parametros_brasil %>%
  filter(
    parameter %in% c(
      "alpha",
      "beta",
      "tau2_mu",
      "tau2_sigma",
      "tau2_xi"
    )
  ) %>%
  mutate(
    parameter = recode(
      parameter,
      alpha = "alpha",
      beta = "beta",
      tau2_mu = "tau² mu",
      tau2_sigma = "tau² sigma",
      tau2_xi = "tau² xi"
    )
  )

p_outros <- ggplot(
  dados_outros,
  aes(
    x = median,
    y = ""
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.alpha = 0.15
  ) +
  facet_wrap(
    ~parameter,
    scales = "free_x",
    ncol = 2
  ) +
  labs(
    x = "Mediana posterior",
    y = NULL,
    title = "Distribuição dos demais parâmetros estruturais"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_outros)

resumo_gev_duracao <- parametros_brasil %>%
  filter(
    parameter %in% c(
      "mu",
      "sigma",
      "xi_d"
    )
  ) %>%
  group_by(
    parameter,
    duration_h,
    duration_days
  ) %>%
  summarise(
    n_estacoes = n(),
    
    minimo = min(
      median,
      na.rm = TRUE
    ),
    
    q05 = quantile(
      median,
      0.05,
      na.rm = TRUE
    ),
    
    q25 = quantile(
      median,
      0.25,
      na.rm = TRUE
    ),
    
    mediana = median(
      median,
      na.rm = TRUE
    ),
    
    q75 = quantile(
      median,
      0.75,
      na.rm = TRUE
    ),
    
    q95 = quantile(
      median,
      0.95,
      na.rm = TRUE
    ),
    
    maximo = max(
      median,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    parameter,
    duration_h
  )

print(
  resumo_gev_duracao,
  n = Inf,
  width = Inf
)

library(dplyr)
library(ggplot2)

font <- "sans"
fontsize <- 8

dados_gev_resumo <- resumo_gev_duracao %>%
  mutate(
    parameter_plot = recode(
      parameter,
      mu = "mu",
      sigma = "sigma",
      xi_d = "xi"
    )
  )

p_gev <- ggplot(
  dados_gev_resumo,
  aes(
    x = duration_days,
    y = mediana
  )
) +
  
  # P05-P95
  geom_ribbon(
    aes(
      ymin = q05,
      ymax = q95
    ),
    alpha = 0.15
  ) +
  
  # Q25-Q75
  geom_ribbon(
    aes(
      ymin = q25,
      ymax = q75
    ),
    alpha = 0.30
  ) +
  
  # Mediana
  geom_line(
    linewidth = 0.6
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  facet_wrap(
    ~parameter_plot,
    scales = "free_y",
    nrow = 1,
    labeller = label_parsed
  ) +
  
  scale_x_continuous(
    breaks = c(
      1, 2, 3, 4, 5, 6, 7, 10
    )
  ) +
  
  labs(
    x = "Duração (dias)",
    y = "Mediana posterior entre estações"
  ) +
  
  theme_bw() +
  
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
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = font,
      size = fontsize
    ),
    
    panel.border = element_rect(
      linewidth = 0.4
    ),
    
    panel.spacing = unit(
      4,
      "pt"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

print(p_gev)

ggsave(
  "C:/Larissa/parametros_GEV_por_duracao.png",
  plot = p_gev,
  width = 15,
  height = 6,
  units = "cm",
  dpi = 600
)

dados_gev_plot <- parametros_brasil %>%
  filter(
    parameter %in% c(
      "mu",
      "sigma",
      "xi_d"
    )
  )

p_gev_box <- ggplot(
  dados_gev_plot,
  aes(
    x = factor(duration_days),
    y = median
  )
) +
  geom_boxplot(
    outlier.alpha = 0.08,
    width = 0.65
  ) +
  facet_wrap(
    ~parameter,
    scales = "free_y",
    nrow = 1,
    labeller = as_labeller(
      c(
        mu = "mu",
        sigma = "sigma",
        xi_d = "xi"
      ),
      label_parsed
    )
  ) +
  labs(
    x = "Duração (dias)",
    y = "Mediana posterior",
    title = "Distribuição dos parâmetros GEV por duração"
  ) +
  theme_minimal(base_size = 11)

print(p_gev_box)

library(dplyr)
library(tidyr)
library(readr)

# =============================================================================
# 1. PARÂMETROS ESTRUTURAIS
# =============================================================================

diagnosticos_gpkg <- flags_diagnostico |>
  transmute(
    UF = as.character(UF),
    gauge_code = as.character(gauge_code),
    flag_diagnostico = as.logical(flag_diagnostico)
  ) |>
  distinct(
    gauge_code,
    .keep_all = TRUE
  )

estruturais_gpkg <- parametros_brasil |>
  filter(
    parameter %in% c(
      "alpha",
      "beta",
      "H",
      "theta",
      "xi_bar",
      "tau2_mu",
      "tau2_sigma",
      "tau2_xi"
    )
  ) |>
  transmute(
    UF = as.character(UF),
    gauge_code = as.character(gauge_code),
    
    parametro = as.character(parameter),
    
    mediana = as.numeric(median),
    q025 = as.numeric(q025),
    q975 = as.numeric(q975)
  ) |>
  
  left_join(
    diagnosticos_gpkg |>
      select(
        gauge_code,
        flag_diagnostico
      ),
    by = "gauge_code"
  ) |>
  
  left_join(
    metadados_todos |>
      transmute(
        gauge_code = as.character(gauge_code),
        lat = as.numeric(lat),
        long = as.numeric(long)
      ) |>
      distinct(
        gauge_code,
        .keep_all = TRUE
      ),
    by = "gauge_code"
  )

estruturais_sf <- estruturais_gpkg |>
  filter(
    is.finite(lat),
    is.finite(long)
  ) |>
  st_as_sf(
    coords = c("long", "lat"),
    crs = 4674,
    remove = FALSE
  )

# =============================================================================
# 2. PARÂMETROS GEV POR DURAÇÃO
# =============================================================================

parametros_gev_gpkg <- parametros_brasil |>
  filter(
    parameter %in% c(
      "mu",
      "sigma",
      "xi_d"
    )
  ) |>
  
  mutate(
    UF = as.character(UF),
    gauge_code = as.character(gauge_code),
    
    parametro = recode(
      parameter,
      mu = "mu",
      sigma = "sigma",
      xi_d = "xi"
    )
  ) |>
  
  transmute(
    UF,
    gauge_code,
    duracao_h = as.numeric(duration_h),
    duracao_dias = as.numeric(duration_days),
    parametro,
    
    mediana = as.numeric(median),
    q025 = as.numeric(q025),
    q975 = as.numeric(q975)
  ) |>
  
  tidyr::pivot_wider(
    names_from = parametro,
    values_from = c(
      mediana,
      q025,
      q975
    ),
    names_glue = "{parametro}_{.value}"
  ) |>
  
  left_join(
    diagnosticos_gpkg |>
      select(
        gauge_code,
        flag_diagnostico
      ),
    by = "gauge_code"
  ) |>
  
  left_join(
    metadados_todos |>
      transmute(
        gauge_code = as.character(gauge_code),
        lat = as.numeric(lat),
        long = as.numeric(long)
      ) |>
      distinct(
        gauge_code,
        .keep_all = TRUE
      ),
    by = "gauge_code"
  )

parametros_gev_sf <- parametros_gev_gpkg |>
  filter(
    is.finite(lat),
    is.finite(long)
  ) |>
  st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4674,
    remove = FALSE
  )


# =============================================================================
# 7. SALVAR CSV
# =============================================================================

readr::write_excel_csv(
  parametros_gpkg,
  "C:/Larissa/parametros_modelo_BR_geopackage.csv"
)

arquivo_gpkg_parametros <- file.path(
  "C:/Larissa/",
  "BR_parametros_modelo.gpkg"
)

# 1. parâmetros estruturais
sf::st_write(
  estruturais_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "parametros_estruturais",
  quiet = TRUE
)

# 2. parâmetros GEV
sf::st_write(
  parametros_gev_sf,
  dsn = arquivo_gpkg_parametros,
  layer = "parametros_gev",
  append = TRUE,
  quiet = TRUE
)

library(dplyr)
library(readr)
library(purrr)
library(tidyr)

DIR_DIAG <- "C:/Larissa/diagnosticos"

# =============================================================================
# 1. JUNTAR TODAS AS TABELAS IDF
# =============================================================================

arquivos_idf <- list.files(
  DIR_DIAG,
  pattern = "^tabela_idf\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

idf_brasil <- map_dfr(
  arquivos_idf,
  function(f) {
    tryCatch(
      read_csv(
        f,
        show_col_types = FALSE
      ),
      error = function(e) NULL
    )
  }
)

cat(
  "Estações com IDF:",
  n_distinct(idf_brasil$gauge_code),
  "\n"
)


# =============================================================================
# 2. METADADOS COM COMPRIMENTO DA SÉRIE
# =============================================================================

arquivos_meta <- list.files(
  DIR_DIAG,
  pattern = "^metadados_estacao\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

metadados_idf <- map_dfr(
  arquivos_meta,
  function(f) {
    
    tryCatch(
      {
        x <- read_csv(
          f,
          show_col_types = FALSE
        )
        
        x %>%
          transmute(
            UF = as.character(UF),
            gauge_code = as.character(gauge_code),
            n_years = as.integer(n_years),
            lat = if ("lat" %in% names(.)) as.numeric(lat) else NA_real_,
            long = if ("long" %in% names(.)) as.numeric(long) else NA_real_
          ) %>%
          slice(1)
      },
      error = function(e) NULL
    )
  }
) %>%
  distinct(
    gauge_code,
    .keep_all = TRUE
  )


# =============================================================================
# 3. INCERTEZA RELATIVA DAS CURVAS DE TR = 100 ANOS
# =============================================================================

incerteza_idf <- idf_brasil %>%
  filter(
    return_period_years == 100
  ) %>%
  mutate(
    largura_relativa =
      (q975 - q025) / median
  ) %>%
  group_by(
    UF,
    gauge_code
  ) %>%
  summarise(
    inc_rel_100 =
      median(
        largura_relativa,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# =============================================================================
# 4. TABELA DE CANDIDATAS
# =============================================================================

metadados_idf <- metadados_idf %>%
  mutate(
    gauge_code = as.character(gauge_code),
    UF = as.character(UF)
  )

incerteza_idf <- incerteza_idf %>%
  mutate(
    gauge_code = as.character(gauge_code),
    UF = as.character(UF)
  )

candidatas_idf <- metadados_idf %>%
  left_join(
    incerteza_idf,
    by = c(
      "UF",
      "gauge_code"
    )
  ) %>%
  filter(
    is.finite(n_years),
    is.finite(inc_rel_100)
  ) %>%
  arrange(
    desc(n_years)
  )

print(
  candidatas_idf %>%
    slice_head(n = 30),
  n = 30,
  width = Inf
)

candidatas_idf <- metadados_idf %>%
  left_join(
    incerteza_idf,
    by = c(
      "UF",
      "gauge_code"
    )
  ) %>%
  filter(
    is.finite(n_years),
    is.finite(inc_rel_100)
  ) %>%
  arrange(
    desc(n_years)
  )

print(
  candidatas_idf %>%
    slice_head(n = 30),
  n = 30,
  width = Inf
)

limite_anos <- quantile(
  candidatas_idf$n_years,
  0.90,
  na.rm = TRUE
)

candidatas_longas <- candidatas_idf %>%
  filter(
    n_years >= limite_anos
  )

# Mais longa
sta_mais_longa <- candidatas_longas %>%
  arrange(
    desc(n_years)
  ) %>%
  slice(1)


# Banda mais estreita entre as longas
sta_banda_estreita <- candidatas_longas %>%
  arrange(
    inc_rel_100
  ) %>%
  filter(
    gauge_code != sta_mais_longa$gauge_code
  ) %>%
  slice(1)


# Banda mais larga entre as longas
sta_banda_larga <- candidatas_longas %>%
  arrange(
    desc(inc_rel_100)
  ) %>%
  filter(
    !gauge_code %in% c(
      sta_mais_longa$gauge_code,
      sta_banda_estreita$gauge_code
    )
  ) %>%
  slice(1)


# Caso próximo da incerteza mediana
mediana_incerteza <- median(
  candidatas_longas$inc_rel_100,
  na.rm = TRUE
)

sta_intermediaria <- candidatas_longas %>%
  mutate(
    distancia =
      abs(
        inc_rel_100 -
          mediana_incerteza
      )
  ) %>%
  filter(
    !gauge_code %in% c(
      sta_mais_longa$gauge_code,
      sta_banda_estreita$gauge_code,
      sta_banda_larga$gauge_code
    )
  ) %>%
  arrange(
    distancia
  ) %>%
  slice(1)


estacoes_exemplo <- bind_rows(
  sta_mais_longa %>%
    mutate(tipo = "Série mais longa"),
  
  sta_banda_estreita %>%
    mutate(tipo = "Baixa incerteza"),
  
  sta_intermediaria %>%
    mutate(tipo = "Incerteza intermediária"),
  
  sta_banda_larga %>%
    mutate(tipo = "Maior incerteza entre séries longas")
)

print(
  estacoes_exemplo,
  n = Inf,
  width = Inf
)

xi_estacao <- parametros_brasil %>%
  filter(
    parameter == "xi_d"
  ) %>%
  group_by(
    UF,
    gauge_code
  ) %>%
  summarise(
    xi_mediano_duracoes =
      median(
        median,
        na.rm = TRUE
      ),
    
    xi_min =
      min(
        median,
        na.rm = TRUE
      ),
    
    xi_max =
      max(
        median,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    gauge_code = as.character(gauge_code)
  )

candidatas_idf2 <- candidatas_idf %>%
  mutate(
    gauge_code = as.character(gauge_code)
  ) %>%
  left_join(
    xi_estacao,
    by = c(
      "UF",
      "gauge_code"
    )
  )

limite_anos <- quantile(
  candidatas_idf2$n_years,
  0.90,
  na.rm = TRUE
)

candidatas_longas <- candidatas_idf2 %>%
  filter(
    n_years >= limite_anos
  )

# Xi mais positivo entre séries longas
xi_positivo <- candidatas_longas %>%
  arrange(
    desc(xi_mediano_duracoes)
  ) %>%
  slice_head(n = 15)

print(
  xi_positivo,
  n = 15,
  width = Inf
)


# Xi mais negativo entre séries longas
xi_negativo <- candidatas_longas %>%
  arrange(
    xi_mediano_duracoes
  ) %>%
  slice_head(n = 15)

print(
  xi_negativo,
  n = 15,
  width = Inf
)

candidatas_longas <- candidatas_longas %>%
  mutate(
    regiao = case_when(
      UF %in% c("AC","AP","AM","PA","RO","RR","TO") ~ "Norte",
      UF %in% c("AL","BA","CE","MA","PB","PE","PI","RN","SE") ~ "Nordeste",
      UF %in% c("DF","GO","MT","MS") ~ "Centro-Oeste",
      UF %in% c("ES","MG","RJ","SP") ~ "Sudeste",
      UF %in% c("PR","RS","SC") ~ "Sul",
      TRUE ~ NA_character_
    )
  )

print(xi_positivo, n = 15, width = Inf)
print(xi_negativo, n = 15, width = Inf)

estacoes_exemplo <- tibble::tribble(
  ~UF,  ~gauge_code, ~tipo,
  "PR", "2549006",   "Série mais longa",
  "PE", "735041",    "Xi positivo elevado",
  "MG", "1842004",   "Xi negativo elevado",
  "SC", "2849002",   "Xi positivo - Sul",
  "CE", "739006",    "Xi negativo - série longa",
  "RN", "535011",    "Maior incerteza"
)

estacoes_exemplo <- estacoes_exemplo %>%
  left_join(
    candidatas_idf2,
    by = c("UF", "gauge_code")
  )

print(
  estacoes_exemplo,
  n = Inf,
  width = Inf
)

library(dplyr)
library(ggplot2)

# ============================================================
# 1. Selecionar as estações de exemplo
# ============================================================

estacoes_exemplo <- tibble::tribble(
  ~UF,  ~gauge_code, ~tipo,
  "PR", "2549006",   "Série mais longa",
  "PE", "735041",    "Xi positivo elevado",
  "SP", "2346057",   "Xi negativo elevado",
  "SC", "2849002",   "Xi positivo - Sul",
  "CE", "438051",    "Xi negativo - série longa",
  "RN", "535011",    "Maior incerteza"
)


# ============================================================
# 2. Preparar dados das IDFs
# ============================================================

idf_exemplos <- idf_brasil %>%
  mutate(
    gauge_code = as.character(gauge_code)
  ) %>%
  inner_join(
    estacoes_exemplo,
    by = c(
      "UF",
      "gauge_code"
    )
  ) %>%
  mutate(
    painel = paste0(
      UF,
      " - ",
      gauge_code,
      "\n",
      tipo
    ),
    TR = factor(
      return_period_years,
      levels = sort(unique(return_period_years))
    )
  )


# ============================================================
# 3. Conferir tempos de retorno disponíveis
# ============================================================

sort(
  unique(idf_exemplos$return_period_years)
)


# ============================================================
# 4. Gráfico
# ============================================================

p_idf <- ggplot(
  idf_exemplos,
  aes(
    x = duration_days,
    y = median,
    group = TR,
    color = TR,
    fill = TR
  )
) +
  
  # Intervalo de credibilidade de 95%
  geom_ribbon(
    aes(
      ymin = q025,
      ymax = q975
    ),
    alpha = 0.08,
    color = NA
  ) +
  
  # Mediana posterior
  geom_line(
    linewidth = 0.85
  ) +
  
  geom_point(
    size = 1.4
  ) +
  
  facet_wrap(
    ~painel,
    ncol = 2,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = c(
      1, 2, 3, 4, 5, 6, 7, 10
    )
  ) +
  
  labs(
    x = "Duração (dias)",
    y = "Intensidade (mm/h)",
    color = "Tempo de retorno\n(anos)",
    fill = "Tempo de retorno\n(anos)"  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_idf)

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

font <- "sans"
fontsize <- 8


# ============================================================
# 1. ESTAÇÕES DE EXEMPLO
# ============================================================

estacoes_exemplo <- tibble::tribble(
  ~UF,  ~gauge_code,
  "PR", "2549006",
  "PE", "735041",
  "MG", "1842004",
  "SC", "2849002",
  "CE", "739006",
  "RN", "535011"
)


# ============================================================
# 2. PARÂMETROS DAS ESTAÇÕES SELECIONADAS
# ============================================================

parametros_idf_exemplos <- parametros_brasil %>%
  mutate(
    gauge_code = as.character(gauge_code)
  ) %>%
  filter(
    gauge_code %in% estacoes_exemplo$gauge_code,
    parameter %in% c(
      "H",
      "xi_bar",
      "theta"
    )
  ) %>%
  select(
    UF,
    gauge_code,
    parameter,
    median
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = parameter,
    values_from = median
  )


# ============================================================
# 3. PREPARAR DADOS DAS IDFs
# ============================================================

idf_exemplos_plot <- idf_brasil %>%
  mutate(
    gauge_code = as.character(gauge_code)
  ) %>%
  inner_join(
    estacoes_exemplo,
    by = c(
      "UF",
      "gauge_code"
    )
  ) %>%
  filter(
    return_period_years %in% c(
      2,
      10,
      25,
      50,
      100,
      200
    )
  ) %>%
  mutate(
    
    # título simples do painel
    painel = paste0(
      UF,
      " - ",
      gauge_code
    ),
    
    # ordem da legenda
    TR = factor(
      return_period_years,
      levels = c(
        2,
        10,
        25,
        50,
        100,
        200
      )
    )
  )


# ============================================================
# 4. TEXTO DOS PARÂMETROS
# ============================================================

parametros_idf_exemplos <- parametros_idf_exemplos %>%
  mutate(
    painel = paste0(
      UF,
      " - ",
      gauge_code
    ),
    
    texto_parametros = paste0(
      "\u03be\u0304 = ",
      sprintf("%.3f", xi_bar),
      "\nH = ",
      sprintf("%.3f", H),
      "\n\u03b8 = ",
      sprintf("%.3f", theta)
    )
  )


# ============================================================
# 5. PALETA DOS TEMPOS DE RETORNO
# ============================================================

cores_TR <- c(
  "2"   = "#D55E00",
  "10"  = "#E69F00",
  "25"  = "#009E73",
  "50"  = "#0072B2",
  "100" = "#CC79A7",
  "200" = "#7B3294"
)


# ============================================================
# 6. GRÁFICO
# ============================================================

p_idf <- ggplot(
  idf_exemplos_plot,
  aes(
    x = duration_days,
    y = median,
    group = TR,
    color = TR
  )
) +
  
  # curvas IDF
  geom_line(
    linewidth = 0.65
  ) +
  
  # ----------------------------------------------------------
# eixo X logarítmico
# ----------------------------------------------------------

scale_x_log10(
  breaks = c(
    1,
    2,
    3,
    5,
    10
  ),
  labels = c(
    "1",
    "2",
    "3",
    "5",
    "10"
  )
) +
  
  # ----------------------------------------------------------
# eixo Y logarítmico
# ----------------------------------------------------------

scale_y_log10(
  labels = scales::label_number(
    accuracy = 0.1
  )
) +
  
  # ----------------------------------------------------------
# cores
# ----------------------------------------------------------

scale_color_manual(
  name = "Tempo de retorno (anos)",
  values = cores_TR,
  breaks = c(
    "2",
    "10",
    "25",
    "50",
    "100"
  )
) +
  
  # ----------------------------------------------------------
# parâmetros no canto superior direito
# ----------------------------------------------------------

geom_text(
  data = parametros_idf_exemplos,
  aes(
    x = Inf,
    y = Inf,
    label = texto_parametros
  ),
  inherit.aes = FALSE,
  hjust = 1.08,
  vjust = 1.20,
  family = font,
  size = fontsize / .pt,
  lineheight = 1.05
) +
  
  # ----------------------------------------------------------
# painéis
# ----------------------------------------------------------

facet_wrap(
  ~painel,
  ncol = 2,
  scales = "free_y"
) +
  
  labs(
    x = "Duração (dias)",
    y = "Intensidade (mm/h)"
  ) +
  
  # ----------------------------------------------------------
# legenda em uma única linha
# ----------------------------------------------------------

guides(
  color = guide_legend(
    nrow = 1,
    byrow = TRUE,
    title.position = "left",
    title.hjust = 0.5
  )
) +
  
  # ----------------------------------------------------------
# tema
# ----------------------------------------------------------

theme_bw() +
  
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
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      family = font,
      size = fontsize,
      face = "plain"
    ),
    
    panel.border = element_rect(
      linewidth = 0.4
    ),
    
    panel.spacing = unit(
      5,
      "pt"
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.title = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.text = element_text(
      family = font,
      size = fontsize
    ),
    
    legend.key.width = unit(
      14,
      "pt"
    ),
    
    legend.spacing.x = unit(
      4,
      "pt"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )


# ============================================================
# 7. VISUALIZAR
# ============================================================

print(p_idf)


# ============================================================
# 8. SALVAR
# ============================================================

ggsave(
  "C:/Larissa/IDF_estacoes_exemplo_sem_incerteza.png",
  plot = p_idf,
  width = 15,
  height = 15,
  units = "cm",
  dpi = 600
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(scales)

font <- "sans"
fontsize <- 8

# ============================================================
# 1. ESTAÇÕES DE EXEMPLO
# ============================================================

estacoes_exemplo <- tibble::tribble(
  ~UF,  ~gauge_code,
  "PR", "2549006",
  "PE", "735041",
  "MG", "1842004",
  "SC", "2849002",
  "CE", "739006",
  "RN", "535011"
)

# ============================================================
# 2. PARÂMETROS DAS ESTAÇÕES SELECIONADAS
# ============================================================

parametros_idf_exemplos <- parametros_brasil %>%
  mutate(
    gauge_code = as.character(gauge_code)
  ) %>%
  filter(
    gauge_code %in% estacoes_exemplo$gauge_code,
    parameter %in% c("H", "xi_bar", "theta")
  ) %>%
  select(
    UF,
    gauge_code,
    parameter,
    median
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = parameter,
    values_from = median
  ) %>%
  mutate(
    painel = paste0(UF, " - ", gauge_code),
    texto_parametros = paste0(
      "\u03be\u0304 = ", sprintf("%.3f", xi_bar),
      "\nH = ", sprintf("%.3f", H),
      "\n\u03b8 = ", sprintf("%.3f", theta)
    )
  )

# ============================================================
# 3. FUNÇÃO PARA GERAR O PLOT
# ============================================================

plot_idf_exemplos_ic <- function(
    idf_brasil,
    estacoes_exemplo,
    parametros_idf_exemplos,
    tr_escolhidos,
    cores_TR,
    arquivo_saida = NULL
) {
  
  idf_exemplos_plot <- idf_brasil %>%
    mutate(
      gauge_code = as.character(gauge_code)
    ) %>%
    inner_join(
      estacoes_exemplo,
      by = c("UF", "gauge_code")
    ) %>%
    filter(
      return_period_years %in% tr_escolhidos
    ) %>%
    mutate(
      painel = paste0(UF, " - ", gauge_code),
      TR = factor(
        return_period_years,
        levels = tr_escolhidos
      )
    )
  
  p_idf_ic <- ggplot(
    idf_exemplos_plot,
    aes(
      x = duration_days,
      y = median,
      group = TR,
      color = TR,
      fill = TR
    )
  ) +
    
    # faixa de incerteza
    geom_ribbon(
      aes(
        ymin = q025,
        ymax = q975
      ),
      alpha = 0.14,
      color = NA
    ) +
    
    # linha da mediana
    geom_line(
      linewidth = 0.70
    ) +
    
    geom_point(
      size = 1.2
    ) +
    
    # texto dos parâmetros
    geom_text(
      data = parametros_idf_exemplos,
      aes(
        x = Inf,
        y = Inf,
        label = texto_parametros
      ),
      inherit.aes = FALSE,
      hjust = 1.08,
      vjust = 1.20,
      family = font,
      size = fontsize / .pt,
      lineheight = 1.05
    ) +
    
    facet_wrap(
      ~painel,
      ncol = 2,
      scales = "free_y"
    ) +
    
    scale_x_log10(
      breaks = c(1, 2, 3, 5, 10),
      labels = c("1", "2", "3", "5", "10")
    ) +
    
    scale_y_log10(
      labels = label_number(accuracy = 0.1)
    ) +
    
    scale_color_manual(
      name = "Tempo de retorno (anos)",
      values = cores_TR
    ) +
    
    scale_fill_manual(
      name = "Tempo de retorno (anos)",
      values = cores_TR
    ) +
    
    labs(
      x = "Duração (dias)",
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
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      text = element_text(family = font, size = fontsize),
      axis.text = element_text(family = font, size = fontsize),
      axis.title = element_text(family = font, size = fontsize),
      strip.background = element_blank(),
      strip.text = element_text(
        family = font,
        size = fontsize,
        face = "plain"
      ),
      panel.border = element_rect(linewidth = 0.4),
      panel.spacing = unit(5, "pt"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(
        family = font,
        size = fontsize
      ),
      legend.text = element_text(
        family = font,
        size = fontsize
      ),
      legend.key.width = unit(14, "pt"),
      legend.spacing.x = unit(4, "pt"),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  print(p_idf_ic)
  
  if (!is.null(arquivo_saida)) {
    ggsave(
      arquivo_saida,
      plot = p_idf_ic,
      width = 15,
      height = 15,
      units = "cm",
      dpi = 600
    )
  }
  
  invisible(p_idf_ic)
}

cores_TR_1 <- c(
  "5"   = "#E69F00",
  "25"  = "#009E73",
  "200" = "#CC79A7"
)

p_idf_ic_1 <- plot_idf_exemplos_ic(
  idf_brasil = idf_brasil,
  estacoes_exemplo = estacoes_exemplo,
  parametros_idf_exemplos = parametros_idf_exemplos,
  tr_escolhidos = c(5, 25, 200),
  cores_TR = cores_TR_1,
  arquivo_saida = "C:/Larissa/IDF_estacoes_exemplo_IC_TR_5_25_200.png"
)

cores_TR_2 <- c(
  "10"  = "#009E73",
  "50"  = "#0072B2",
  "100" = "#CC79A7"
)

p_idf_ic_2 <- plot_idf_exemplos_ic(
  idf_brasil = idf_brasil,
  estacoes_exemplo = estacoes_exemplo,
  parametros_idf_exemplos = parametros_idf_exemplos,
  tr_escolhidos = c(10, 50, 100),
  cores_TR = cores_TR_2,
  arquivo_saida = "C:/Larissa/IDF_estacoes_exemplo_IC_TR_10_50_100.png"
)

library(dplyr)
library(sf)

# =============================================================================
# 1. PREPARAR METADADOS DAS ESTAÇÕES
# =============================================================================

meta_idf_gpkg <- metadados_idf |>
  transmute(
    UF = as.character(UF),
    gauge_code = as.character(gauge_code),
    n_years = as.integer(n_years),
    lat = as.numeric(lat),
    long = as.numeric(long)
  ) |>
  distinct(
    gauge_code,
    .keep_all = TRUE
  )

niveis_retorno_gpkg <- idf_brasil |>
  transmute(
    UF = as.character(UF),
    gauge_code = as.character(gauge_code),
    
    duracao_h = as.numeric(duration_h),
    duracao_dias = as.numeric(duration_days),
    
    periodo_retorno_anos = as.numeric(return_period_years),
    
    mediana = as.numeric(median),
    q025 = as.numeric(q025),
    q975 = as.numeric(q975)
  ) |>
  
  # ---------------------------------------------------------------------------
# diagnóstico: somente flag sintética
# ---------------------------------------------------------------------------

left_join(
  diagnosticos_gpkg |>
    select(
      gauge_code,
      flag_diagnostico
    ),
  by = "gauge_code"
) |>
  
  # ---------------------------------------------------------------------------
# metadados espaciais + comprimento da série
# ---------------------------------------------------------------------------

left_join(
  meta_idf_gpkg |>
    select(
      gauge_code,
      n_years,
      lat,
      long
    ),
  by = "gauge_code"
)

niveis_retorno_gpkg |>
  summarise(
    n = n(),
    
    n_estacoes =
      n_distinct(gauge_code),
    
    n_duracoes =
      n_distinct(duracao_dias),
    
    n_periodos_retorno =
      n_distinct(periodo_retorno_anos),
    
    sem_lat =
      sum(!is.finite(lat)),
    
    sem_long =
      sum(!is.finite(long)),
    
    sem_flag =
      sum(is.na(flag_diagnostico))
  )

niveis_retorno_sf <- niveis_retorno_gpkg |>
  filter(
    is.finite(lat),
    is.finite(long)
  ) |>
  sf::st_as_sf(
    coords = c(
      "long",
      "lat"
    ),
    crs = 4674,
    remove = FALSE
  )

arquivo_gpkg_retorno <- file.path(
  "C:/Larissa/",
  "BR_periodos_retorno.gpkg"
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

