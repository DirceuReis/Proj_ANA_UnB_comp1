## Figuras do QC Sub-Diário - UNIPLU-BR
##
## Mapeia as estações que sobreviveram ao controle de qualidade de cada rede
## sub-diária, em painéis.
##
## Insumos: de 02b1 a 02b4, "series_<rede>_final.parquet", pelas estações que
##          restaram; de 01b1 a 01b4, "inventario_<rede>.csv", pelas coordenadas.
## Saídas:  "Figura 6. Estações finais por rede.png".

##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","sf","geobr","ggspatial","ggplot2","scales"))

##### 1. Library imports ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(sf)
library(geobr)
library(ggspatial)
library(ggplot2)
library(scales)

##### 2. Caminhos ----------------------------------------------------------------------------------------------------------
saida <- "C:/UNB/UNIPLU/02b. Figuras_QC_Sub"

dir.create(saida, showWarnings = FALSE, recursive = TRUE)

##### 3. Padrão gráfico ----------------------------------------------------------------------------------------------------
font     <- "sans"
fontsize <- 8
textsize <- fontsize / .pt     # pt -> mm, para geom_text e annotate
epsg     <- 5880               # SIRGAS 2000 / Policônica do Brasil (métrico)

LARG_MAPA <- 16                # cm
ALT_MAPA  <- 15.5

theme_set(
  theme_bw() +
    theme(
      panel.grid       = element_blank(),
      axis.ticks       = element_line(linewidth = 0.3),
      strip.background = element_blank(),
      strip.text       = element_text(size = fontsize, family = font, hjust = 0),
      text             = element_text(size = fontsize, family = font),
      panel.spacing    = unit(2, "pt")
    )
)

escala_mapa <- annotation_scale(
  location = "br", width_hint = 0.25, height = unit(1, "mm"),
  text_family = font, text_cex = 0.5,
  bar_cols = c("black", "aliceblue"), line_width = 0.3
)

norte_mapa <- annotation_north_arrow(
  location = "tr", width = unit(0.7, "cm"), height = unit(0.7, "cm"),
  which_north = "true",
  style = north_arrow_nautical(text_family = font, text_size = fontsize * 0.7,
                               fill = c("black", "white"), line_width = 0.5)
)

##### 4. As redes, suas cores e seus arquivos --------------------------------------------------------------------
redes <- tibble(
  rede      = c("CEMADEN", "Telemetria", "INMET subdiário", "ICEA"),
  prefixo   = c("CEMADEN", "Telemetria", "INMETsub", "ICEA"),
  dir_serie = c("02b1. QC_Sub_CEMADEN", "02b2. QC_Sub_Telemetria",
                "02b4. QC_Sub_INMETsubdiario", "02b3. QC_Sub_ICEA"),
  dir_inv   = c("01b1. Series_Sub_CEMADEN", "01b2. Series_Sub_Telemetria",
                "01b4. Series_Sub_INMETsubdiario", "01b3. Series_Sub_ICEA")
)

## as mesmas cores do 00b, para que as figuras do relatório se correspondam
rede_cores <- c("CEMADEN"         = "#1565c0",
                "Telemetria"      = "#b71c1c",
                "INMET subdiário" = "#00838f",
                "ICEA"            = "#f9a825")

##### 5. Estações que restaram, com coordenada -------------------------------------------------------------------
estacoes <- bind_rows(lapply(seq_len(nrow(redes)), function(i) {

  arq_serie <- file.path("C:/UNB/UNIPLU", redes$dir_serie[i],
                         paste0("series_", redes$prefixo[i], "_final.parquet"))
  arq_inv   <- file.path("C:/UNB/UNIPLU", redes$dir_inv[i],
                         paste0("inventario_", redes$prefixo[i], ".csv"))

  if (!file.exists(arq_serie)) stop("Não achei ", arq_serie)
  if (!file.exists(arq_inv))   stop("Não achei ", arq_inv)

  ## as estações que sobraram ao fim do QC
  finais <- open_dataset(arq_serie) %>%
    distinct(gauge_code) %>%
    collect()

  inventario <- read_csv(arq_inv, col_types = cols(gauge_code = col_character(),
                                                   .default   = col_guess()))

  saiu <- inventario %>%
    select(gauge_code, lat, long) %>%
    inner_join(finais, by = "gauge_code") %>%
    filter(!is.na(lat), !is.na(long))

  if (nrow(saiu) < nrow(finais)) {
    cat("  ", redes$rede[i], ": ", nrow(finais) - nrow(saiu),
        " estações sem coordenada, fora do mapa\n", sep = "")
  }

  cat(sprintf("%-16s %4d estações finais | %4d no mapa\n",
              redes$rede[i], nrow(finais), nrow(saiu)))

  saiu$rede <- redes$rede[i]
  saiu
}))

contagem <- estacoes %>% count(rede, name = "estacoes")

rotulos <- setNames(
  sprintf("%s (%s)", contagem$rede,
          format(contagem$estacoes, big.mark = ".", decimal.mark = ",",
                 trim = TRUE)),
  contagem$rede)

estacoes$rede <- factor(estacoes$rede, levels = redes$rede)

estacoes_sf <- st_as_sf(estacoes, coords = c("long", "lat"), crs = 4674)

cat("\nTotal no mapa:", nrow(estacoes_sf), "estações\n")

##### 6. Contorno do país e dos estados ------------------------------------------------------------------------
uf.sf   <- geobr::read_state(year = 2025, showProgress = FALSE)
br.sf   <- geobr::read_country(year = 2025, showProgress = FALSE)
bbox.br <- sf::st_bbox(br.sf)

enquadra <- coord_sf(crs = st_crs(epsg), default_crs = st_crs(4674),
                     xlim = bbox.br[c(1, 3)], ylim = bbox.br[c(2, 4)])

##### 7. Figura 8. Estações finais por rede --------------------------------------------------------------------
Figura6 <- ggplot() +
  geom_sf(data = uf.sf, color = "grey70", fill = "grey97", linewidth = 0.08) +
  geom_sf(data = br.sf, color = "grey20", fill = NA, linewidth = 0.3) +
  geom_sf(data = estacoes_sf, aes(color = rede), size = 0.3, shape = 16) +
  scale_color_manual(values = rede_cores, guide = "none") +
  facet_wrap(~ rede, nrow = 2,
             labeller = labeller(rede = rotulos)) +
  escala_mapa + norte_mapa + enquadra +
  labs(x = NULL, y = NULL)

ggsave(file.path(saida, "Figura 8. Estações finais por rede.png"), Figura6,
       width = LARG_MAPA, height = ALT_MAPA, units = "cm", dpi = 300, bg = "white")

cat("Figura 8. Estações finais por rede.png salvo em", saida, "\n")