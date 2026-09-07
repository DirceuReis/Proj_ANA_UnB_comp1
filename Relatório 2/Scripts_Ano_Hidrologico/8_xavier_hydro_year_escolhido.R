# Ano hidrológico Xavier (BR-DWGD) — método HÍBRIDO rbar (critério final)
#
# Regra (igual Uniplu):
#   se rbar >= 0,30 → θ+π (mes_oposto)
#   senão           → último mês seco (mes_relatorio; limiar 15% no NC)
#
# Pré-requisito:
#   Rscript Scripts_Fluxo_AnoHidrologico/7_xavier_pot_hidro.R
#   (hill.adapt + evd::clusters — mesma metodologia Uniplu)
#
# Entrada: Dados_Xavier/ano_hidro_out/xavier_ano_hidro_pot.nc
# Saída:   resultados/BR/dataframes/df_xavier_rbar_hibrido.*
#          resultados/BR/pictures/ano_hidro_xavier_hibrido/
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(geobr)
library(readr)
library(patchwork)
library(terra)

# --- caminhos ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"
DIR_XAVIER_NC <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2024/Dados_Xavier/ano_hidro_out"

OUT_TAG <- "BR"
DIR_OUT     <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF      <- file.path(DIR_OUT, "dataframes")
DIR_OUT_PIC <- file.path(DIR_OUT, "pictures", "ano_hidro_xavier_hibrido")

setwd(DIR_FLUXO)
dir.create(DIR_DF, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT_PIC, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Parâmetros
# ---------------------------------------------------------------------------
RBAR_HI <- 0.30
F_NC <- file.path(DIR_XAVIER_NC, "xavier_ano_hidro_pot.nc")
if (!file.exists(F_NC)) {
  stop(
    "Falta ", F_NC,
    "\nGere antes com: Rscript Scripts_Fluxo_AnoHidrologico/7_xavier_pot_hidro.R"
  )
}

meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")
meses_lab <- c("Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
               "Jul", "Ago", "Set", "Out", "Nov", "Dez")
cores_mes <- hsv(
  h = ((seq_len(12) - 1L) / 12 + 0.02) %% 1,
  s = 0.82, v = 0.93
)
names(cores_mes) <- as.character(seq_len(12))

UF_SUL <- c("RS", "SC", "PR")
BB_SUL <- c(xmin = -58.5, xmax = -47.0, ymin = -34.2, ymax = -21.8)

# ---------------------------------------------------------------------------
# Ler grade + recalcular híbrido (garante limiar adotado)
# ---------------------------------------------------------------------------
# O writeCDF do passo 7 pode gravar 5 bandas Z1=1..5 (sem subds nomeados).
# Ordem gravada em 7_xavier_pot_hidro.R:
#   1 n_picos | 2 rbar | 3 mes_oposto | 4 mes_hibrido_rbar | 5 mes_relatorio
message("Lendo: ", F_NC)
r_all <- rast(F_NC)
if (nlyr(r_all) < 5L) {
  stop("NC inesperado (nlyr=", nlyr(r_all), "). Rode de novo o passo 7.")
}

pick_layer <- function(r, nome, idx) {
  nm <- names(r)
  hit <- which(nm == nome | grepl(paste0("(^|_)", nome, "($|_|=)"), nm))
  if (length(hit) == 1L) return(r[[hit]])
  r[[idx]]
}

r_ops <- pick_layer(r_all, "mes_oposto", 3L)
r_rel <- pick_layer(r_all, "mes_relatorio", 5L)
r_bar <- pick_layer(r_all, "rbar", 2L)
r_hib_nc <- tryCatch(pick_layer(r_all, "mes_hibrido_rbar", 4L), error = function(e) NULL)

names(r_ops) <- "mes_oposto"
names(r_rel) <- "mes_relatorio"
names(r_bar) <- "rbar"

stk <- c(r_ops, r_rel, r_bar)
names(stk) <- c("mes_oposto", "mes_relatorio", "rbar")

df0 <- as.data.frame(stk, xy = TRUE, na.rm = FALSE)
# terra às vezes repete nomes; força colunas únicas
names(df0) <- make.unique(names(df0), sep = "_")
# se make.unique alterou, reconstrói pelos 3 layers + xy
if (!all(c("mes_oposto", "mes_relatorio", "rbar") %in% names(df0))) {
  names(df0) <- c("x", "y", "mes_oposto", "mes_relatorio", "rbar")
}

df <- df0 %>%
  dplyr::filter(
    is.finite(mes_oposto), is.finite(mes_relatorio), is.finite(rbar),
    mes_oposto >= 1, mes_oposto <= 12,
    mes_relatorio >= 1, mes_relatorio <= 12
  ) %>%
  dplyr::mutate(
    mes_oposto = as.integer(round(mes_oposto)),
    mes_relatorio = as.integer(round(mes_relatorio)),
    confia_theta = rbar >= RBAR_HI,
    mes_hibrido_rbar = dplyr::if_else(confia_theta, mes_oposto, mes_relatorio),
    fonte_hibrido_rbar = dplyr::if_else(confia_theta, "theta_pi", "ultimo_mes_seco"),
    mes_hibrido_nome = meses_pt[mes_hibrido_rbar]
  )

# Raster do híbrido recalculado (para mapas)
r_hib <- r_ops
v_ops <- values(r_ops)
v_rel <- values(r_rel)
v_rbar <- values(r_bar)
v_hib <- ifelse(
  is.finite(v_rbar) & v_rbar >= RBAR_HI & is.finite(v_ops),
  v_ops,
  ifelse(is.finite(v_rel), v_rel, NA_real_)
)
# Só onde há rbar finito (domínio POT)
v_hib[!is.finite(v_rbar)] <- NA_real_
values(r_hib) <- v_hib
names(r_hib) <- "mes_hibrido_rbar"

message(
  "Células: ", nrow(df),
  " | rbar≥", RBAR_HI, ": ",
  sum(df$confia_theta),
  " (", round(100 * mean(df$confia_theta), 1), "%) → θ+π"
)
message(
  "  rbar<", RBAR_HI, ": ",
  sum(!df$confia_theta),
  " (", round(100 * mean(!df$confia_theta), 1), "%) → último mês seco"
)

# ---------------------------------------------------------------------------
# Tabelas
# ---------------------------------------------------------------------------
resumo <- df %>%
  summarise(
    n_celulas = n(),
    rbar_med = median(rbar, na.rm = TRUE),
    pct_rbar_ge_limiar = round(100 * mean(rbar >= RBAR_HI), 1),
    pct_fonte_theta_pi = round(100 * mean(fonte_hibrido_rbar == "theta_pi"), 1),
    pct_fonte_ultimo_seco = round(100 * mean(fonte_hibrido_rbar == "ultimo_mes_seco"), 1)
  )

tab_mes <- df %>%
  count(mes_hibrido_rbar, name = "n") %>%
  mutate(
    mes_nome = meses_pt[mes_hibrido_rbar],
    pct = round(100 * n / sum(n), 1)
  ) %>%
  arrange(mes_hibrido_rbar)

write_excel_csv(resumo, file.path(DIR_DF, "df_xavier_rbar_hibrido_resumo.csv"))
write_excel_csv(tab_mes, file.path(DIR_DF, "df_xavier_rbar_hibrido_meses.csv"))
write_excel_csv(df, file.path(DIR_DF, "df_xavier_rbar_hibrido_grade.csv"))
saveRDS(df, file.path(DIR_DF, "df_xavier_rbar_hibrido_grade.rds"))

# Raster oficial do método escolhido (cópia leve no fluxo)
writeRaster(
  r_hib,
  filename = file.path(DIR_DF, "xavier_mes_hibrido_rbar.tif"),
  overwrite = TRUE
)

message("Brasil:"); print(as.data.frame(resumo))
message("Meses (híbrido):")
print(as.data.frame(tab_mes), row.names = FALSE)

# ---------------------------------------------------------------------------
# Mapas
# ---------------------------------------------------------------------------
ufs <- read_state(year = 2020, showProgress = FALSE) %>%
  st_transform(4326) %>%
  mutate(estado = abbrev_state)
br <- read_country(year = 2020, showProgress = FALSE) %>% st_transform(4326)
sul_ufs <- ufs %>% filter(estado %in% UF_SUL)

tema_mapa <- theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey35"),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(4, 4, 2, 4)
  )

legenda_circular <- function() {
  d <- tibble(mes = 1:12, lab = meses_lab, xmin = 0:11, xmax = 1:12)
  ggplot(d) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0.58, ymax = 1.0, fill = factor(mes)),
      color = "white", linewidth = 0.35
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = 1.22, label = lab),
      size = 3.2, fontface = "bold", color = "grey15"
    ) +
    scale_fill_manual(values = unname(cores_mes), guide = "none") +
    coord_polar(theta = "x", start = -pi / 2 - pi / 12, direction = 1) +
    ylim(0, 1.45) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 2, 0))
}

bbox_pad <- function(pol, pad_x = 0.6, pad_y = 0.4) {
  bb <- st_bbox(pol)
  list(
    xlim = c(as.numeric(bb["xmin"]) - pad_x, as.numeric(bb["xmax"]) + pad_x),
    ylim = c(as.numeric(bb["ymin"]) - pad_y, as.numeric(bb["ymax"]) + pad_y)
  )
}

raster_to_df_mes <- function(r) {
  as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
    rename(mes = 3) %>%
    mutate(mes_f = factor(as.character(as.integer(round(mes))),
                          levels = as.character(1:12)))
}

mapa_raster <- function(df_r, titulo, sub = NULL, pol_base, pol_contorno,
                        lim = NULL) {
  if (is.null(lim)) lim <- bbox_pad(pol_base)
  ggplot() +
    geom_raster(data = df_r, aes(x = x, y = y, fill = mes_f)) +
    geom_sf(data = pol_contorno, fill = NA, color = "grey35", linewidth = 0.22) +
    scale_fill_manual(values = cores_mes, drop = FALSE, na.value = "grey90") +
    coord_sf(xlim = lim$xlim, ylim = lim$ylim, expand = FALSE, clip = "off") +
    labs(title = titulo, subtitle = sub) +
    tema_mapa
}

salvar_mapa_legenda <- function(g_mapa, arquivo, w = 1800, h = 2100) {
  g_out <- g_mapa / legenda_circular() + plot_layout(heights = c(8.2, 1.8))
  ggsave(arquivo, g_out, width = w, height = h, units = "px", dpi = 180, bg = "white")
}

sub_txt <- paste0(
  "rbar ≥ ", RBAR_HI, " → θ+π; senão → último mês seco | células = ",
  format(nrow(df), big.mark = ".")
)

df_hib <- raster_to_df_mes(r_hib)
lim_br <- bbox_pad(br)
lim_sul <- list(
  xlim = BB_SUL[c("xmin", "xmax")],
  ylim = BB_SUL[c("ymin", "ymax")]
)

salvar_mapa_legenda(
  mapa_raster(
    df_hib,
    "Xavier — ano hidrológico híbrido rbar",
    sub_txt, br, ufs, lim_br
  ),
  file.path(DIR_OUT_PIC, "fig_xavier_brasil_hibrido_rbar.png"),
  w = 1800, h = 2100
)

df_hib_sul <- df_hib %>%
  filter(
    x >= BB_SUL["xmin"], x <= BB_SUL["xmax"],
    y >= BB_SUL["ymin"], y <= BB_SUL["ymax"]
  )

salvar_mapa_legenda(
  mapa_raster(
    df_hib_sul,
    "Xavier Sul — híbrido rbar",
    sub_txt, sul_ufs, sul_ufs, lim_sul
  ),
  file.path(DIR_OUT_PIC, "fig_xavier_sul_hibrido_rbar.png"),
  w = 1700, h = 2000
)

# Fonte do híbrido
df_fonte <- df %>%
  mutate(fonte_f = factor(
    fonte_hibrido_rbar,
    levels = c("theta_pi", "ultimo_mes_seco")
  ))

tema_fonte <- tema_mapa + theme(legend.position = "right")

g_fonte_br <- ggplot() +
  geom_raster(data = df_fonte, aes(x = x, y = y, fill = fonte_f)) +
  geom_sf(data = ufs, fill = NA, color = "grey35", linewidth = 0.2) +
  scale_fill_manual(
    values = c(theta_pi = "#1b9e77", ultimo_mes_seco = "#d95f02"),
    labels = c(
      theta_pi = paste0("θ+π (rbar ≥ ", RBAR_HI, ")"),
      ultimo_mes_seco = paste0("Último mês seco (rbar < ", RBAR_HI, ")")
    ),
    name = "Fonte"
  ) +
  coord_sf(xlim = lim_br$xlim, ylim = lim_br$ylim, expand = FALSE) +
  labs(
    title = "Xavier — fonte do híbrido rbar",
    subtitle = paste0(
      "θ+π: ", sum(df$fonte_hibrido_rbar == "theta_pi"),
      "  |  último mês seco: ", sum(df$fonte_hibrido_rbar == "ultimo_mes_seco")
    )
  ) +
  tema_fonte
ggsave(
  file.path(DIR_OUT_PIC, "fig_xavier_brasil_fonte_hibrido.png"),
  g_fonte_br, width = 1800, height = 1900, units = "px", dpi = 180, bg = "white"
)

g_fonte_sul <- ggplot() +
  geom_raster(
    data = df_fonte %>%
      filter(
        x >= BB_SUL["xmin"], x <= BB_SUL["xmax"],
        y >= BB_SUL["ymin"], y <= BB_SUL["ymax"]
      ),
    aes(x = x, y = y, fill = fonte_f)
  ) +
  geom_sf(data = sul_ufs, fill = NA, color = "grey35", linewidth = 0.3) +
  scale_fill_manual(
    values = c(theta_pi = "#1b9e77", ultimo_mes_seco = "#d95f02"),
    labels = c(
      theta_pi = paste0("θ+π (rbar ≥ ", RBAR_HI, ")"),
      ultimo_mes_seco = paste0("Último mês seco (rbar < ", RBAR_HI, ")")
    ),
    name = "Fonte"
  ) +
  coord_sf(xlim = lim_sul$xlim, ylim = lim_sul$ylim, expand = FALSE) +
  labs(title = paste0("Xavier Sul — fonte do híbrido (", RBAR_HI, ")")) +
  tema_fonte
ggsave(
  file.path(DIR_OUT_PIC, "fig_xavier_sul_fonte_hibrido.png"),
  g_fonte_sul, width = 1700, height = 1900, units = "px", dpi = 180, bg = "white"
)

# rbar
df_rbar <- as.data.frame(r_bar, xy = TRUE, na.rm = TRUE) %>% rename(rbar = 3)
g_rbar <- ggplot() +
  geom_raster(data = df_rbar, aes(x = x, y = y, fill = rbar)) +
  geom_sf(data = ufs, fill = NA, color = "grey35", linewidth = 0.2) +
  scale_fill_viridis_c(
    option = "C", limits = c(0, 1), name = "rbar",
    breaks = c(0, 0.2, 0.3, 0.5, 0.75, 1)
  ) +
  coord_sf(xlim = lim_br$xlim, ylim = lim_br$ylim, expand = FALSE) +
  labs(
    title = "Xavier — concentração dos extremos (rbar)",
    subtitle = paste0("Limiar do híbrido = ", RBAR_HI)
  ) +
  tema_fonte
ggsave(
  file.path(DIR_OUT_PIC, "fig_xavier_rbar_brasil.png"),
  g_rbar, width = 1800, height = 1900, units = "px", dpi = 180, bg = "white"
)

ggsave(
  file.path(DIR_OUT_PIC, "fig_legenda_meses.png"),
  legenda_circular(),
  width = 900, height = 900, units = "px", dpi = 150, bg = "white"
)

message("Tabelas: ", DIR_DF)
message("Figuras: ", DIR_OUT_PIC)
message("Concluído — Xavier híbrido rbar ≥ ", RBAR_HI, ".")
