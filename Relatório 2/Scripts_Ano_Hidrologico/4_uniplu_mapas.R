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

# --- caminhos (saídas nesta pasta; shapefiles fora) ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"
DIR_DADOS <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

DIR_SHP      <- file.path(DIR_DADOS, "shp")
DIR_OUT      <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF       <- file.path(DIR_OUT, "dataframes")
DIR_PIC      <- file.path(DIR_OUT, "pictures")
DIR_SETAS_UF <- file.path(DIR_PIC, "fig3_setas_uf")

setwd(DIR_FLUXO)
dir.create(DIR_PIC, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_SETAS_UF, recursive = TRUE, showWarnings = FALSE)

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

stations <- stations %>%
  mutate(mes_color = letters[as.integer(month(dmy("1-1-1999") + round(theta * 180 / pi) - 1))])

col_month <- brewer.pal(n = 6, name = "Spectral")
col_month <- c(col_month, rev(col_month))

g <-
  ggplot() +
  geom_sf(data = america, fill = "grey88", linewidth = 0.1) +
  geom_sf(data = brasil,  fill = "grey78", alpha = 0.4, linewidth = 0.6) +
  geom_sf(data = hidrografia2, col = "darkblue", alpha = 0.6) +
  geom_sf(data = hidrografia3, col = "darkblue", alpha = 0.3) +
  geom_sf(data = hidrografia4, col = "darkblue", alpha = 0.1) +
  geom_spoke(
    data = stations,
    aes(x = long, y = lat, angle = theta, color = mes_color, radius = s1),
    linewidth = 0.4,
    arrow = arrow(length = unit(.1, "cm"))
  ) +
  scale_colour_manual(name = "site", values = col_month) +
  scale_radius(labels = NULL, trans = "identity", range = c(.5, 3), guide = "none") +
  coord_sf(
    xlim = limite[c(1, 3)],
    ylim = limite[c(2, 4)],
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude") +
  annotation_scale(location = "br", bar_cols = c("black", "white")) +
  annotation_north_arrow(
    location = "br", which_north = "true",
    pad_x = unit(0.7, "in"), pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering(
      fill = c("black", "white"),
      line_col = "grey20"
    )
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    title = element_text(size = 12, color = "black", face = "bold"),
    legend.key = element_rect(fill = "white"),
    legend.justification = "top",
    legend.text = element_text(size = 12, color = "black"),
    legend.spacing.y = unit(0.3, "cm")
  )

ggsave(
  filename = file.path(DIR_PIC, "fig3_uniplu.png"),
  plot = g,
  scale = 1,
  width = 3200,
  height = 2650,
  units = "px",
  dpi = 300
)

message("Salvo pictures/fig3_uniplu.png (", nrow(stations), " postos)")

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

tema_mapa <- theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    strip.background = element_rect(fill = "white"),
    legend.position = "bottom"
  )

n_uf <- n_distinct(stations$estado)
dir_uf <- DIR_SETAS_UF
dir.create(dir_uf, recursive = TRUE, showWarnings = FALSE)

plots_uf <- list()
for (uf in sort(unique(stations$estado))) {
  pol <- ufs %>% filter(abbrev_state == uf)
  sta_uf <- stations %>% filter(estado == uf)
  if (nrow(pol) == 0 || nrow(sta_uf) == 0) {
    next
  }
  bb <- st_bbox(pol)
  pad <- max(1.6, max(sta_uf$s1, na.rm = TRUE) * 1.25)
  g_one <- ggplot() +
    geom_sf(data = pol, fill = "white", color = "black", linewidth = 0.5) +
    geom_spoke(
      data = sta_uf,
      aes(x = long, y = lat, angle = theta, color = mes_pico_lab, radius = s1),
      linewidth = 0.4,
      arrow = arrow(length = unit(0.1, "cm"))
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
    tema_mapa +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  plots_uf[[uf]] <- g_one
  ggsave(
    filename = file.path(dir_uf, paste0(uf, ".png")),
    plot = g_one,
    width = 2000,
    height = 1800,
    units = "px",
    dpi = 150
  )
}

g_painel <- wrap_plots(plots_uf, ncol = min(4, max(1, n_uf)))

ggsave(
  filename = file.path(DIR_PIC, "fig3_uniplu_estados.png"),
  plot = g_painel,
  width = min(4800, 1400 * min(4, max(1, n_uf))),
  height = min(4000, 1200 * ceiling(n_uf / min(4, max(1, n_uf)))),
  units = "px",
  dpi = 150
)
message("Salvo pictures/fig3_uniplu_estados.png")
message("Salvo mapas por UF em pictures/fig3_setas_uf/")

