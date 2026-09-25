## Caracterização do Inventário Nacional - Estações Diárias - UNIPLU-BR
##
## Pré-requisito: rodar "00. Leitura_Diarias_UNIPLU.R" antes (usa
## "series_diarias.parquet" e "inventario_estacoes.csv", em "00. Series").

##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","ggplot2","stringr"))
#install.packages(c("sf","geobr"))

##### 1. Library imports ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(ggplot2)
library(stringr)
library(sf)
library(geobr)
library(scales)

##### 2. Caminhos ----------------------------------------------------------------------------------------------------------
entrada  <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/00. Series"
saida <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/01. Series_Diarias"
dir.create(saida, showWarnings = FALSE, recursive = TRUE)

##### 3. Filtragem séries diárias ---------------------------------------------------------------------------------------------------
inventario_completo <- read_csv(file.path(entrada, "inventario_estacoes.csv"), show_col_types = FALSE)

# Mantém apenas as redes diárias Hidroweb e INMET (12.640 estações)
inventario_diario <- inventario_completo %>% filter(network %in% c("Hidroweb diário", "INMET diário")) 

cat("Estações no inventário completo :", nrow(inventario_completo), "\n")
cat("Estações com time_step = 1440   :", nrow(inventario_diario), "\n")

write_csv(inventario_diario, file.path(saida, "inventario_diario.csv")) # salva inventário diário em excel

df_completo <- read_parquet(file.path(entrada, "series_diarias.parquet")) # lê o parquet criado em 00
df_diario <- df_completo %>% filter(gauge_code %in% inventario_completo$gauge_code) # filtra o parquet pra base diária

write_parquet(df_diario, file.path(saida, "series_diarias_1440.parquet")) #salva o parquet
cat("Registros diários (time_step = 1440):", nrow(df_diario), "de", nrow(df_completo), "no total\n")
cat("Série diária salva em:", file.path(saida, "series_diarias_1440.parquet"), "\n")

##### 4. Mapa estações de séries diárias -----------------------------------------------------------------------------------------
brasil_sf <- geobr::read_state(year = 2020, showProgress = FALSE)
brasil_sf$area_km2 <- as.numeric(sf::st_area(brasil_sf)) / 1e6
mapa_dados <- inventario_completo %>% filter(!is.na(lat), !is.na(long))
mapa_dados_sf <- st_as_sf(mapa_dados, coords = c("long", "lat"), crs = 4674)
estados_sf    <- brasil_sf %>% filter(abbrev_state %in% unique(mapa_dados$state))

rede_cores <- c(
  "Hidroweb diário" = "#2e7d32",
  "INMET diário"    = "#6a1b9a")

rede_formas <- c(
  "Hidroweb diário" = 17,  # triângulo
  "INMET diário"    = 18)

Figura6 <- ggplot() +
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

ggsave(file.path(saida, "Figura 6. Inventário Diário.png"), Figura6, width = 16, height = 12, dpi = 300)
cat("Figura 6. Inventário Diário.png salvo.\n")