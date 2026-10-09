# Mapa de setas do ano hidrológico (Uniplu)
# Entrada: dataframes/df_hydroyear_uniplu.rds
# Saída:   pictures/fig3_uniplu.png

library(dplyr)
library(ggplot2)
library(ggspatial)
library(lubridate)
library(sf)
library(RColorBrewer)
library(scales)
library(geobr)
library(patchwork)
library(grid)

# --- caminhos (saídas nesta pasta; shapefiles fora) ---
DIR_FLUXO <- "Relatório 2/Scripts_Ano_Hidrologico"
DIR_DADOS <- "Relatório 2/Scripts_Ano_Hidrologico"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

DIR_SHP      <- file.path(DIR_DADOS, "shp")
DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_PIC      <- file.path(DIR_OUT, "pictures")
DIR_SETAS_UF <- file.path(DIR_PIC, "fig3_setas_uf")

dir.create(DIR_PIC, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_SETAS_UF, recursive = TRUE, showWarnings = FALSE)

font <- "sans"
fontsize <- 8
textsize <- fontsize / .pt

FIG_WIDTH <- 15

tema_base <- theme_bw() +
  theme(panel.grid = element_blank(), axis.ticks = element_line(linewidth = 0.3),
        legend.background = element_blank(), legend.key = element_blank(),
        legend.key.height = unit(fontsize, "pt" ),
        strip.background = element_blank(),
        text = element_text( family = font, size = fontsize),
        axis.text = element_text(family = font, size = fontsize),
        axis.title = element_text(family = font, size = fontsize),
        legend.title = element_text(family = font, size = fontsize),
        legend.text = element_text(family = font, size = fontsize),
        strip.text = element_text(family = font, size = fontsize),
        panel.spacing.y = unit( 0.5, "pt"),
        plot.margin = margin( t = 5, r = 5, b = 5, l = 5)
        )

america <- st_read(file.path(DIR_SHP, "america_do_sul.gpkg"), quiet = TRUE)
brasil  <- america %>% filter(nome == "Brasil")

hidrografia2 <- st_read(file.path(DIR_SHP, "hidro_lvl2_EPSG5880.shp"), quiet = TRUE)
hidrografia3 <- st_read(file.path(DIR_SHP, "hidro_lvl3_EPSG5880.shp"), quiet = TRUE)
hidrografia4 <- st_read(file.path(DIR_SHP, "hidro_lvl4_EPSG5880.shp"), quiet = TRUE)

stations <- readRDS(file.path(DIR_DF, "df_hydroyear_uniplu.rds"))

stations <- stations %>%
  filter(!is.na(lat), !is.na(long), !is.na(theta), !is.na(rbar)) %>%
  filter(estado != "BO")

if (!is.null(TEST_UF)) {
  stations <- stations %>% filter(estado == TEST_UF)
}

if (nrow(stations) == 0) {
  stop("Nenhum posto em df_hydroyear_uniplu.rds",
       if (!is.null(TEST_UF)) paste0(" com estado=", TEST_UF,
         ". Rode scripts 2 e 3 depois de mudar TEST_UF.") else ".")
}

message("Mapa: ", nrow(stations), " postos",
        if (!is.null(TEST_UF)) paste0(" (", TEST_UF, ")") else "")

limite <- st_bbox(brasil)

stations <- stations %>%
  mutate(s0 = log(rbar))

stations <- stations %>%
  mutate(s1 = scales::rescale(rbar, to = c(0.05, 1.5)))

meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")

stations <- stations %>%
  mutate(
    s1 = scales::rescale(rbar, to = c(0.05, 1.5)),
    mes_color = meses_pt[
      as.integer(
        month(dmy("1-1-1999") + round(theta * 180 / pi) - 1)
      )
    ],
    mes_color = factor(mes_color, levels = meses_pt)
  )

#col_month <- brewer.pal(n = 6, name = "Spectral")
col_month <- c(
  "jan" = "#3B4CC0",
  "fev" = "#6F58C9",
  "mar" = "#9C4DC4",
  "abr" = "#C43A9A",
  "mai" = "#E34A6F",
  "jun" = "#F26B38",
  "jul" = "#F4A62A",
  "ago" = "#D8C63A",
  "set" = "#8FCB3C",
  "out" = "#35B779",
  "nov" = "#1FA4A9",
  "dez" = "#2A78B8"
)
col_month <- c(col_month, rev(col_month))

g_brasil <-
  ggplot() +
  geom_sf(data = america, fill = "grey95", color = "grey75", linewidth = 0.15) +
  geom_sf(data = brasil,  fill = "grey98", color = "grey50",, alpha = 0.4, linewidth = 0.4) +
  geom_sf(data = hidrografia2, color = "darkblue", alpha = 0.6, linewidth = 0.2) +
  geom_sf(data = hidrografia3, color = "darkblue", alpha = 0.3, linewidth = 0.15) +
  geom_sf(data = hidrografia4, color = "darkblue", alpha = 0.1, linewidth = 0.1) +
  geom_spoke(
    data = stations,
    aes(x = long, y = lat, angle = theta, color = mes_color, radius = s1),
    linewidth = 0.35,
    arrow = arrow(length = unit(0.08, "cm"))
  ) +
  scale_colour_manual(
    name = "Mês típico do extremo",
    values = col_month
  ) +
  scale_radius(
    labels = NULL,
    trans = "identity",
    range = c(0.5, 3),
    guide = "none"
  ) +
  coord_sf(
    xlim = limite[c(1, 3)],
    ylim = limite[c(2, 4)],
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude") +
  annotation_scale(
    location = "br",
    width_hint = 0.20,
    height = unit( 1.5, "mm"), text_family = font,
    text_cex = 0.6,
    bar_cols = c(
      "black",
      "aliceblue"
    ),
    line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr",
    width = unit(
      1,
      "cm"
    ),
    height = unit(
      1,
      "cm"
    ),
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
  
  labs( x = NULL, y = NULL ) +
  
  tema_base +
  
  theme(
    # sem coordenadas nos eixos
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    # legenda horizontal
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    
    legend.key.spacing.x = unit(
      fontsize,
      "pt"
    ),
    
    legend.spacing.x = unit(
      1,
      "pt"
    ),
    
    legend.margin = margin(
      t = -3
    )
  ) +
  
  guides(
  color = guide_legend(
    nrow = 2,
    byrow = TRUE,
    title.position = "top",
    override.aes = list(
      linewidth = 1.0
    )
  )
)


ggsave(
  filename = file.path(DIR_PIC, "fig3_uniplu.png"),
  plot = g_brasil,
  width = FIG_WIDTH,
  height = 12,
  units = "cm",
  dpi = 300
)

message("Salvo: ", file.path(DIR_PIC, "fig3_uniplu.png"))

# Segunda figura: recorte por estado com malha do geobr
meses_pt <- c("jan", "fev", "mar", "abr", "mai", "jun",
              "jul", "ago", "set", "out", "nov", "dez")
col_month_nome <- brewer.pal(n = 6, name = "Spectral")
col_month_nome <- c(col_month_nome, rev(col_month_nome))
names(col_month_nome) <- meses_pt

stations <- stations %>%
  mutate(
    mes_pico_lab = meses_pt[as.integer(month(
      as.Date("1999-01-01") + round(theta * 365 / (2 * pi)) - 1
    ))],
    mes_pico_lab = factor(mes_pico_lab, levels = meses_pt)
  )

ufs <- read_state(year = 2020, showProgress = FALSE) %>%
  st_transform(4326) %>%
  filter(abbrev_state %in% unique(stations$estado)) %>%
  mutate(estado = abbrev_state)

tema_mapa_uf <- tema_base +
  theme(
    axis.ticks = element_blank(),
    legend.position = "bottom"
  )

n_uf <- n_distinct(stations$estado)

plots_uf <- list()

for (uf in sort(unique(stations$estado))) {
  pol <- ufs %>% filter(abbrev_state == uf)
  sta_uf <- stations %>% filter(estado == uf)
  
  if (nrow(pol) == 0 || nrow(sta_uf) == 0) next
  
  bb  <- st_bbox(pol)
  pad <- max(1.6, max(sta_uf$s1, na.rm = TRUE) * 1.25)
  
  g_one <- ggplot() +
    geom_sf(data = pol, fill = "white", color = "black", linewidth = 0.35) +
    geom_spoke(
      data = sta_uf,
      aes(x = long, y = lat, angle = theta, color = mes_pico_lab, radius = s1),
      linewidth = 0.35,
      arrow = arrow(length = unit(0.08, "cm"))
    ) +
    scale_colour_manual(
      name = "Mês típico do extremo",
      values = col_month_nome,
      drop = FALSE
    ) +
    coord_sf(
      xlim = c(as.numeric(bb["xmin"]) - pad, as.numeric(bb["xmax"]) + pad),
      ylim = c(as.numeric(bb["ymin"]) - pad, as.numeric(bb["ymax"]) + pad),
      expand = FALSE,
      clip = "off"
    ) +
    labs(title = uf, x = "Longitude", y = "Latitude") +
    tema_mapa_uf +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  plots_uf[[uf]] <- g_one
  
  ggsave(
    filename = file.path(DIR_SETAS_UF, paste0(uf, ".png")),
    plot = g_one,
    width = FIG_WIDTH,
    height = 11,
    units = "cm",
    dpi = 300
  )
}

# Painel com todos os estados
ncol_panel <- min(4, max(1, n_uf))
nrow_panel <- ceiling(n_uf / ncol_panel)

g_painel <- wrap_plots(plots_uf, ncol = ncol_panel)

ggsave(
  filename = file.path(DIR_PIC, "fig3_uniplu_estados.png"),
  plot = g_painel,
  width = FIG_WIDTH,
  height = max(8, 3.8 * nrow_panel),
  units = "cm",
  dpi = 300
)

message("Salvo: ", file.path(DIR_PIC, "fig3_uniplu_estados.png"))
message("Salvo mapas por UF em: ", DIR_SETAS_UF)
