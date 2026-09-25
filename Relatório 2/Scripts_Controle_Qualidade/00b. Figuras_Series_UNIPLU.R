## FIGURAS DO INVENTÁRIO - UNIPLU-BR
##
## Insumos:  00. Series/inventario_estacoes.csv
##           00. Series/series_diarias.parquet
##

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
path <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/00. Series"

##### 3. Insumos: inventário e contorno dos estados  ------------------------------------------------------------------------
inventario <- read_csv(file.path(path, "inventario_estacoes.csv"),
                       col_types = cols(gauge_code = col_character(), .default = col_guess()))
mapa_dados <- inventario %>% filter(!is.na(lat), !is.na(long))

brasil_sf <- geobr::read_state(year = 2020, showProgress = FALSE)
brasil_sf$area_km2 <- as.numeric(sf::st_area(brasil_sf)) / 1e6   # área geodésica (m² -> km²)
cat("inventário:", nrow(inventario), "estações |", nrow(mapa_dados), "com coordenada\n")

##### 4. Figura 1. Mapa das estações por rede  ------------------------------------------------------------------------

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

mapa_dados_sf <- st_as_sf(mapa_dados, coords = c("long", "lat"), crs = 4674)
estados_sf    <- brasil_sf %>% filter(abbrev_state %in% unique(mapa_dados$state))

Figura1 <- ggplot() +
  geom_sf(data = estados_sf, fill = "grey95", color = "grey50", linewidth = 0.4) +
  geom_sf(data = mapa_dados_sf, aes(color = network, shape = network),
          size = 2.2, alpha = 0.75) +
  scale_color_manual(values = rede_cores, name = "Rede") +
  scale_shape_manual(values = rede_formas, name = "Rede") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 22),
    legend.title    = element_text(size = 22, face = "bold"),
    axis.text.x     = element_text(size = 22, color = "black"),
    axis.text.y     = element_text(size = 22, color = "black"),
    axis.title      = element_text(size = 22, color = "black")
  ) +
  annotation_scale(location = "bl", width_hint = 0.3,
                   text_cex   = 1.7,                  
                   height     = unit(0.35, "cm"),  
                   line_width = 1.2,                
                   text_pad   = unit(0.2, "cm"))+
  annotation_north_arrow(
    location = "bl",
    pad_y  = unit(0.8, "cm"),
    height = unit(2.0, "cm"),    
    width  = unit(2.0, "cm"),
    style  = north_arrow_fancy_orienteering(
      text_size  = 22,           
      line_width = 1.2,         
      text_face  = "bold"))+
  coord_sf(xlim = c(min(mapa_dados$long) - 0.5, max(mapa_dados$long) + 0.5),
           ylim = c(min(mapa_dados$lat)  - 0.5, max(mapa_dados$lat)  + 0.5))

ggsave(file.path(path, "Figura 1. Inventário.png"), Figura1, width = 16, height = 12, dpi = 300)
cat("Figura 1. Inventário.png salvo.\n")

##### 5. Figura 2. Densidade de estações por UF ------------------------------------------------------
uf_numero <- inventario %>% filter(!is.na(state)) %>% count(state, name = "n_estacoes")

uf_sf <- brasil_sf %>%
  left_join(uf_numero, by = c("abbrev_state" = "state")) %>%
  mutate(n_estacoes      = ifelse(is.na(n_estacoes), 0, n_estacoes),
         km2_por_estacao = round(area_km2 / n_estacoes))

Figura2 <- ggplot(uf_sf) +
  geom_sf(aes(fill = km2_por_estacao), color = "white", linewidth = 0.2) +
  geom_sf_text(aes(label = mil(km2_por_estacao)), size = 20, size.unit = "pt") +
  scale_fill_gradient2(
    low = "blue", mid = "grey90", high = "red",
    midpoint = mean(range(uf_sf$km2_por_estacao, na.rm = TRUE, finite = TRUE)),
    labels = mil,
    name = "Área (km²)\npor estação",
    guide = guide_colourbar(
      theme = theme(
        legend.key.height = unit(20, "cm"),
        legend.key.width = unit(0.7,"cm"),
        legend.ticks = element_line(color = "white",linewidth = 0.4),
        legend.frame = element_rect(color = "grey30", linewidth = 0.3)),
      title.position = "top"))+
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 22),
    legend.title    = element_text(size = 22, face = "bold"),
    axis.text.x     = element_text(size = 22, color = "black"),
    axis.text.y     = element_text(size = 22, color = "black"),
    axis.title      = element_text(size = 22, color = "black")
  ) +
  annotation_scale(location = "bl", width_hint = 0.3,
                   text_cex   = 1.7,                  
                   height     = unit(0.35, "cm"),  
                   line_width = 1.2,                
                   text_pad   = unit(0.2, "cm"))+
  annotation_north_arrow(
    location = "bl",
    pad_y  = unit(0.8, "cm"),
    height = unit(2.0, "cm"),    
    width  = unit(2.0, "cm"),
    style  = north_arrow_fancy_orienteering(
      text_size  = 22,           
      line_width = 1.2,         
      text_face  = "bold"))+
  coord_sf(xlim = c(min(mapa_dados$long) - 0.5, max(mapa_dados$long) + 0.5),
           ylim = c(min(mapa_dados$lat)  - 0.5, max(mapa_dados$lat)  + 0.5))

ggsave(file.path(path, "Figura 2. Densidade de estações.png"), Figura2, width = 16, height = 12, dpi = 300)
cat("Figura 2. Densidade de estações.png.\n")

##### 6. Figura 3. Mapa de pontos: nº de anos em atividade  ---------------------------------------------------

TETO_ANOS <- 100

Figura3 <- ggplot() +
  geom_sf(data = brasil_sf, fill = "grey97", color = "grey40", linewidth = 0.3) +
  geom_point(data = mapa_dados, aes(x = long, y = lat, color = n_anos),
             size = 1, alpha = 0.7) +
  scale_color_gradient2(
    low = "red", mid = "grey85", high = "blue",
    midpoint = TETO_ANOS / 2, limits = c(0, TETO_ANOS), oob = scales::squish,
    breaks = seq(0, TETO_ANOS, by = 10),
    name = "Número de\nanos em atividade",
    guide = guide_colourbar(
      theme = theme(
        legend.key.height = unit(20, "cm"),
        legend.key.width = unit(0.7,"cm"),
        legend.ticks = element_line(color = "white",linewidth = 0.4),
        legend.frame = element_rect(color = "grey30", linewidth = 0.3)),
      title.position = "top"))+
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 22),
    legend.title    = element_text(size = 22, face = "bold"),
    axis.text.x     = element_text(size = 22, color = "black"),
    axis.text.y     = element_text(size = 22, color = "black"),
    axis.title      = element_text(size = 22, color = "black")
  ) +
  annotation_scale(location = "bl", width_hint = 0.3,
                   text_cex   = 1.7,                  
                   height     = unit(0.35, "cm"),  
                   line_width = 1.2,                
                   text_pad   = unit(0.2, "cm"))+
  annotation_north_arrow(
    location = "bl",
    pad_y  = unit(0.8, "cm"),
    height = unit(2.0, "cm"),    
    width  = unit(2.0, "cm"),
    style  = north_arrow_fancy_orienteering(
      text_size  = 22,           
      line_width = 1.2,         
      text_face  = "bold"))+
  coord_sf(xlim = c(min(mapa_dados$long) - 0.5, max(mapa_dados$long) + 0.5),
           ylim = c(min(mapa_dados$lat)  - 0.5, max(mapa_dados$lat)  + 0.5))

ggsave(file.path(path, "Figura 3. Mapa anos_atividade.png"), Figura3, width = 16, height = 12, dpi = 300)
cat("Figura 3. Mapa anos_atividade.png salvo \n")

##### 7. Figura 4. Estações por instituição responsável  ------------------------------------------------------------------------
LIMIAR_AGENCIA <- 200

responsavel_numero <- inventario %>%
  filter(!is.na(responsible), responsible != "") %>%
  count(responsible, name = "n") %>%
  filter(n > LIMIAR_AGENCIA) %>%
  arrange(desc(n))


Figura4 <- ggplot(responsavel_numero, aes(x = reorder(responsible, -n), y = n)) +
  geom_col(fill = "grey50") +
  geom_text(aes(label = mil(n)), vjust = -0.4, size = 22, size.unit = "pt") +
  labs(x = "Operadora", y = "Número de estações com dado") +
  scale_y_continuous(labels = mil, expand = expansion(mult = c(0, 0.15))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 22, color = "black"),
        axis.title  = element_text(size = 22, color = "black"))
    
    
ggsave(file.path(path, "Figura 4. Estações por instituição.png"), Figura4, width = 16, height = 12, dpi = 300)
cat("Figura 4. Estações por instituição salva. Agências acima do limiar:", nrow(responsavel_numero), "\n")

##### 8. Figura 5. Estações por ano, por resolução  ------------------------------------------------------------------------

ANO_MIN  <- 1855
ANO_MAX  <- 2025
PASSO_X  <- 5
RECONTAR <- TRUE
COR_DIA  <- "#4C78A8"
COR_SUB  <- "#F58518"

SUB <- c("CEMADEN", "Telemetria", "INMET subdiário", "ICEA")
resol_temporal <- inventario %>%
  mutate(passo_num = suppressWarnings(as.numeric(time_step)),
         resolucao = case_when(!is.na(passo_num) & passo_num <  1440 ~ "Sub-diária",
                               !is.na(passo_num) & passo_num >= 1440 ~ "Diária",
                               network %in% SUB                      ~ "Sub-diária",
                               TRUE                                  ~ "Diária")) %>%
  select(gauge_code, resolucao)
print(count(resol_temporal, resolucao))

arq_numero <- file.path(path, "estacoes_por_ano.csv")
if (RECONTAR) {
  t0 <- Sys.time()
  cat("buscando na série diária (163 milhões de linhas), ~15 s...\n")
  soma_estacoes_ano <- open_dataset(file.path(path, "series_diarias.parquet")) %>%
    filter(!is.na(rain_mm)) %>%
    mutate(ano = year(date)) %>%
    count(ano, gauge_code) %>%
    collect() %>%
    left_join(resol_temporal, by = "gauge_code") %>%
    count(ano, resolucao, name = "estacoes")
  write_csv(estacoes, arq_numero)
  cat(sprintf("contagem pronta em %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else {
  soma_estacoes_ano <- read_csv(arq_numero, show_col_types = FALSE)
}

# grade completa de anos: anos sem estação aparecem como zero
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
  scale_fill_manual(values = c("Diária" = COR_DIA, "Sub-diária" = COR_SUB), guide = "none") +
  scale_x_continuous(breaks = seq(ANO_MIN, ANO_MAX, by = PASSO_X), expand = expansion(mult = 0.01)) +
  scale_y_continuous(labels = label_number(big.mark = ".", decimal.mark = ",", accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Ano", y = "Número de estações") +
  theme_bw(base_size = 22) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 22, color = "black"),
        axis.text.y = element_text(color = "black",size = 22),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = 2, hjust = 0, size = 22))

ggsave(file.path(path, "Figura 5. Estações por ano por resolução.png"), Figura5, width = 16, height = 12, dpi = 300, bg = "white")
cat("Figura 5. Estações por ano por resolução salva.\n")
