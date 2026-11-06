## FIGURAS DO INVENTÁRIO - UNIPLU-BR
##
## Padrão dos mapas: fonte 8 pt, traços finos, legenda dentro do painel,
## enquadramento fixo pela bbox do país e CRS único em `epsg`.
## Insumos: 00a. inventario_estacoes.csv e resumo_estacao_ano.csv

##### 1. Library imports --------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(arrow)
library(sf)
library(geobr)
library(ggspatial)
library(scales)

##### 2. Caminhos  ----------------------------------------------------------------------------------------------------------
dados <- "C:/UNB/UNIPLU/00. Series"

##### 3. Padrão gráfico  ----------------------------------------------------------------------------------------------------
font     <- "sans"
fontsize <- 8
textsize <- fontsize / .pt     # pt -> mm, para geom_text/geom_sf_text/annotate
epsg     <- 5880               # SIRGAS 2000 / Policônica do Brasil (métrico)

## bbox do país projetada em 5880 tem razão largura/altura de 1,115
LARG_MAPA <- 16                # cm
ALT_MAPA  <- 14.5
LARG_GRAF <- 16
ALT_GRAF  <- 10

set_theme(
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      legend.position = "inside",
      legend.position.inside = c(0, 0),
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key.height = unit(fontsize, "pt"),
      strip.background = element_blank(),
      text = element_text(size = fontsize, family = font),
      panel.spacing.y = unit(0.5, "pt")
    )
)

escala_mapa <- annotation_scale(
  location = "br", width_hint = 0.17, height = unit(1.5, "mm"),
  text_family = font, text_cex = 0.6,
  bar_cols = c("black", "aliceblue"), line_width = 0.4
)

norte_mapa <- annotation_north_arrow(
  location = "tr", width = unit(1, "cm"), height = unit(1, "cm"),
  which_north = "true",
  style = north_arrow_nautical(text_family = font, text_size = fontsize,
                               fill = c("black", "white"))
)

tema_mapa <- theme(
  plot.margin            = margin(t = 10, r = 1, l = 1),
  panel.background       = element_rect(fill = "white"),
  legend.key.spacing.y   = unit(0, "pt"),
  legend.key.size        = unit(fontsize * 0.9, "pt"),
  legend.margin          = margin(t = -5, l = -2),
  legend.position        = "inside",
  legend.position.inside = c(0.03, 0.02),
  legend.justification   = c(0, 0)
)

## barra de cor vertical, usada nas Figuras 2 e 3
barra_cor <- function() {
  guide_colourbar(
    theme = theme(
      legend.key.height = unit(25, "mm"),
      legend.key.width  = unit(2.5, "mm"),
      legend.ticks      = element_line(color = "white", linewidth = 0.3),
      legend.frame      = element_rect(color = "grey30", linewidth = 0.2)),
    title.position = "top")
}

##### 4. Inventário e contorno dos estados  ------------------------------------------------------------------------
inventario <- read_csv(file.path(dados, "inventario_estacoes.csv"),
                       col_types = cols(gauge_code = col_character(),
                                        time_step  = col_character(),
                                        .default   = col_guess()))

resumo_por_ano <- read_csv(file.path(dados, "resumo_estacao_ano.csv"),
                           col_types = cols(gauge_code = col_character(),
                                            time_step  = col_character(),
                                            .default   = col_guess()))

mapa_dados    <- inventario %>% filter(!is.na(lat), !is.na(long))
mapa_dados_sf <- st_as_sf(mapa_dados, coords = c("long", "lat"), crs = 4674)

uf.sf   <- geobr::read_state(year = 2025, showProgress = FALSE)
br.sf   <- geobr::read_country(year = 2025, showProgress = FALSE)
bbox.br <- sf::st_bbox(br.sf)

uf.sf$area_km2 <- as.numeric(sf::st_area(uf.sf)) / 1e6   # área geodésica (m² -> km²)
mil <- label_number(big.mark = ".", decimal.mark = ",", accuracy = 1)

enquadra <- coord_sf(crs = st_crs(epsg), default_crs = st_crs(4674),
                     xlim = bbox.br[c(1, 3)], ylim = bbox.br[c(2, 4)])

cat("inventário:", nrow(inventario), "estações |", nrow(mapa_dados), "com coordenada\n")
cat("resumo por ano:", nrow(resumo_por_ano), "linhas estação-ano\n")

##### 5. Figura 1. Mapa das estações por rede  ------------------------------------------------------------------------

rede_cores <- c(
  "CEMADEN"         = "#1565c0",
  "Hidroweb diário" = "#2e7d32",
  "ICEA"            = "#f9a825",
  "INMET diário"    = "#6a1b9a",
  "INMET subdiário" = "#00838f",
  "Telemetria"      = "#b71c1c"
)
rede_formas <- c(
  "CEMADEN"         = 16,  # círculo
  "Hidroweb diário" = 17,  # triângulo
  "ICEA"            = 15,  # quadrado
  "INMET diário"    = 18,  # losango
  "INMET subdiário" = 8,   # asterisco
  "Telemetria"      = 25   # triângulo invertido
)

Figura1 <- ggplot() +
  geom_sf(data = mapa_dados_sf, aes(color = network, shape = network),
          size = 0.35, stroke = 0.08) +
  geom_sf(data = uf.sf, color = "grey30", fill = NA, linewidth = 0.1) +
  geom_sf(data = br.sf, color = "black",  fill = NA, linewidth = 0.4) +
  scale_color_manual(values = rede_cores,  name = "Rede") +
  scale_shape_manual(values = rede_formas, name = "Rede") +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  escala_mapa + norte_mapa + enquadra +
  labs(x = "Longitude", y = "Latitude") +
  tema_mapa

ggsave(file.path(dados, "Figura 1. Inventário.png"), Figura1,
       width = LARG_MAPA, height = ALT_MAPA, units = "cm", dpi = 300, bg = "white")
cat("Figura 1. Inventário.png salvo.\n")

##### 6. Figura 2. Densidade de estações por UF ------------------------------------------------------
uf_numero <- inventario %>% filter(!is.na(state)) %>% count(state, name = "n_estacoes")

uf_sf <- uf.sf %>%
  left_join(uf_numero, by = c("abbrev_state" = "state")) %>%
  mutate(km2_por_estacao = round(area_km2 / n_estacoes))

Figura2 <- ggplot(uf_sf) +
  geom_sf(aes(fill = km2_por_estacao), color = "grey30", linewidth = 0.1) +
  geom_sf(data = br.sf, color = "black", fill = NA, linewidth = 0.4) +
  geom_sf_text(aes(label = mil(km2_por_estacao)),
               size = fontsize * 0.85, size.unit = "pt", family = font) +
  scale_fill_gradient2(
    low = "blue", mid = "grey90", high = "red",
    midpoint = mean(range(uf_sf$km2_por_estacao, na.rm = TRUE, finite = TRUE)),
    labels = mil,
    name = "Área (km²)\npor estação",
    guide = barra_cor()) +
  escala_mapa + norte_mapa + enquadra +
  labs(x = "Longitude", y = "Latitude") +
  tema_mapa

ggsave(file.path(dados, "Figura 2. Densidade de estações.png"), Figura2,
       width = LARG_MAPA, height = ALT_MAPA, units = "cm", dpi = 300, bg = "white")
cat("Figura 2. Densidade de estações.png salvo.\n")

##### 7. Figura 3. Mapa de pontos: nº de anos em atividade  ---------------------------------------------------

TETO_ANOS <- 70

Figura3 <- ggplot() +
  geom_sf(data = mapa_dados_sf, aes(color = n_anos_com_registro),
          shape = 3, size = 0.3, stroke = 0.08) +
  geom_sf(data = uf.sf, color = "grey30", fill = NA, linewidth = 0.1) +
  geom_sf(data = br.sf, color = "black",  fill = NA, linewidth = 0.4) +
  scale_color_gradient2(
    low = "red", mid = "grey85", high = "blue",
    midpoint = TETO_ANOS / 2, limits = c(0, TETO_ANOS), oob = scales::squish,
    breaks = seq(0, TETO_ANOS, by = 10),
    name = "Número de\nanos em atividade",
    guide = barra_cor()) +
  escala_mapa + norte_mapa + enquadra +
  labs(x = "Longitude", y = "Latitude") +
  tema_mapa

ggsave(file.path(dados, "Figura 3. Mapa anos_atividade.png"), Figura3,
       width = LARG_MAPA, height = ALT_MAPA, units = "cm", dpi = 300, bg = "white")
cat("Figura 3. Mapa anos_atividade.png salvo \n")

##### 8. Figura 4. Estações por instituição responsável  ------------------------------------------------------------------------
LIMIAR_AGENCIA <- 200

responsavel_numero <- inventario %>%
  filter(!is.na(responsible), responsible != "") %>%
  count(responsible, name = "n") %>%
  filter(n > LIMIAR_AGENCIA) %>%
  arrange(desc(n))

Figura4 <- ggplot(responsavel_numero, aes(x = reorder(responsible, -n), y = n)) +
  geom_col(fill = "grey50") +
  geom_text(aes(label = mil(n)), vjust = -0.4,
            size = fontsize * 0.9, size.unit = "pt", family = font) +
  labs(x = "Operadora", y = "Número de estações no inventário") +
  scale_y_continuous(labels = mil, expand = expansion(mult = c(0, 0.15))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(dados, "Figura 4. Estações por instituição.png"), Figura4,
       width = LARG_GRAF, height = ALT_GRAF, units = "cm", dpi = 300, bg = "white")
cat("Figura 4. Estações por instituição salva. Agências acima do limiar:", nrow(responsavel_numero), "\n")

##### 9. Figura 5. Estações por ano, por resolução  ------------------------------------------------------------------------

ANO_MIN  <- 1855
ANO_MAX  <- 2025
PASSO_X  <- 5
COR_DIARIA  <- "#4C78A8"
COR_SUBDIARIA  <- "#F58518"

SUB <- c("CEMADEN", "Telemetria", "INMET subdiário", "ICEA")
resol_temporal <- inventario %>%
  mutate(passo_num = suppressWarnings(as.numeric(time_step)),
         resolucao = case_when(!is.na(passo_num) & passo_num <  1440 ~ "Sub-diária",
                               !is.na(passo_num) & passo_num >= 1440 ~ "Diária",
                               network %in% SUB                      ~ "Sub-diária",
                               TRUE                                  ~ "Diária")) %>%
  select(gauge_code, resolucao)
print(count(resol_temporal, resolucao))

soma_estacoes_ano <- resumo_por_ano %>%
  filter(!sem_dado_valido) %>%
  left_join(resol_temporal, by = "gauge_code") %>%
  distinct(gauge_code, ano, resolucao) %>%
  count(ano, resolucao, name = "estacoes")

write_csv(soma_estacoes_ano, file.path(dados, "estacoes_por_ano.csv"))

soma_estacoes_ano <- soma_estacoes_ano %>%
  filter(ano >= ANO_MIN, ano <= ANO_MAX) %>%
  complete(ano = ANO_MIN:ANO_MAX, resolucao = c("Diária", "Sub-diária"),
           fill = list(estacoes = 0)) %>%
  mutate(resolucao = factor(resolucao, c("Diária", "Sub-diária")))

cat("máximo de estações simultâneas:\n")
print(soma_estacoes_ano %>% group_by(resolucao) %>% slice_max(estacoes, n = 1) %>% ungroup())

Figura5 <- ggplot(soma_estacoes_ano, aes(ano, estacoes, fill = resolucao)) +
  geom_col(width = 0.75) +
  facet_wrap(~resolucao, ncol = 1, scales = "free_y",
             labeller = as_labeller(c("Diária" = "Resolução diária",
                                      "Sub-diária" = "Resolução sub-diária"))) +
  scale_fill_manual(values = c("Diária" = COR_DIARIA, "Sub-diária" = COR_SUBDIARIA), guide = "none") +
  scale_x_continuous(breaks = seq(ANO_MIN, ANO_MAX, by = PASSO_X), expand = expansion(mult = 0.01)) +
  scale_y_continuous(labels = mil, expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Ano", y = "Número de estações") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text  = element_text(face = 2, hjust = 0))

ggsave(file.path(dados, "Figura 5. Estações por ano por resolução.png"), Figura5,
       width = LARG_GRAF, height = 12, units = "cm", dpi = 300, bg = "white")
cat("Figura 5. Estações por ano por resolução salva.\n")
