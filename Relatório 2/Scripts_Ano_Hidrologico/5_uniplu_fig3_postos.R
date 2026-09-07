# Figura 3 do relatório (um ponto por posto)
# Cada ponto = vetor médio (xbar, ybar) da estação; raio = rbar.
# Gera o círculo do recorte (Brasil ou UF) e um círculo por estado.
# Entrada: df_hydroyear_uniplu.rds
# Saída:   pictures/fig3_postos_uniplu.png
#          pictures/fig3_postos_uf/<UF>.png

library(dplyr)
library(ggplot2)
library(lubridate)

# --- caminhos (saídas nesta pasta) ---
DIR_FLUXO <- "C:/Users/laris/OneDrive/3. UFC/UFC - 2026/Analises_Relatório_ANA_Outubro/Ano_Hidro/Scripts_Fluxo_AnoHidrologico"

TEST_UF <- NULL   # NULL = Brasil; ex. "CE", "RS"
OUT_TAG <- if (is.null(TEST_UF) || !nzchar(TEST_UF)) "BR" else TEST_UF

DIR_OUT       <- file.path(DIR_FLUXO, "resultados", OUT_TAG)
DIR_DF        <- file.path(DIR_OUT, "dataframes")
DIR_PIC       <- file.path(DIR_OUT, "pictures")
DIR_POSTOS_UF <- file.path(DIR_PIC, "fig3_postos_uf")

setwd(DIR_FLUXO)
dir.create(DIR_PIC, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_POSTOS_UF, recursive = TRUE, showWarnings = FALSE)

get_scale_month <- function() {
  theta <- yday(dmy(paste0("1/", seq(1, 12, 1), "/1999"))) * 2 * pi / 365 - 0.5 * 2 * pi / 365
  theta <- theta - theta[1]
  tibble(xscale = cos(theta), yscale = sin(theta))
}

get_scale_text_month <- function(r) {
  theta <- yday(dmy(paste0("15/", seq(1, 12, 1), "/1999"))) * 2 * pi / 365 - 0.5 * 2 * pi / 365
  tibble(xscale = r * cos(theta), yscale = r * sin(theta))
}

plot_fig3_postos <- function(df, titulo) {

  xy_scale <- get_scale_month()
  x_circle <- cos(seq(0, 2 * pi, length.out = 200))
  y_circle <- sin(seq(0, 2 * pi, length.out = 200))
  txt_mes <- c("Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
               "Jul", "Ago", "Set", "Out", "Nov", "Dez")
  xy_text <- get_scale_text_month(r = 1.05)
  df_txt <- tibble(label = txt_mes, pos_x = xy_text$xscale, pos_y = xy_text$yscale)

  g <- ggplot() +
    geom_point(
      data = df, aes(x = xbar, y = ybar),
      color = "black", size = 1.6, alpha = 0.75
    )

  if (nrow(df) >= 8) {
    g <- g +
      geom_density_2d(
        data = df, aes(x = xbar, y = ybar),
        color = "blue", linewidth = 0.45, bins = 8
      )
  }

  g +
    geom_text(
      data = df_txt,
      aes(x = pos_x, y = pos_y, label = label),
      fontface = "bold", size = 12 / .pt
    ) +
    annotate(
      "text",
      x = 0.025,
      y = c(0.25, 0.5, 0.75) + 0.02,
      label = c("0.25", "0.50", "0.75"),
      size = 9 / .pt, hjust = 0, vjust = 0, fontface = "bold"
    ) +
    annotate("path", x = 1.00 * x_circle, y = 1.00 * y_circle,
             color = "black", linewidth = 0.5) +
    annotate("path", x = 0.75 * x_circle, y = 0.75 * y_circle,
             color = "gray", linewidth = 0.3) +
    annotate("path", x = 0.50 * x_circle, y = 0.50 * y_circle,
             color = "gray", linewidth = 0.3) +
    annotate("path", x = 0.25 * x_circle, y = 0.25 * y_circle,
             color = "gray", linewidth = 0.3) +
    geom_point(data = xy_scale, aes(x = xscale, y = yscale),
               color = "black", size = 1) +
    geom_segment(
      data = xy_scale,
      aes(x = 0, y = 0, xend = xscale, yend = yscale),
      color = "gray", linewidth = 0.3
    ) +
    coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
    labs(title = titulo) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.line  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
}

salvar_fig3 <- function(g, caminho) {
  ggsave(
    filename = caminho,
    plot = g,
    width = 3200,
    height = 2800,
    units = "px",
    dpi = 300
  )
}

stations <- readRDS(file.path(DIR_DF, "df_hydroyear_uniplu.rds")) %>%
  filter(is.finite(xbar), is.finite(ybar), is.finite(rbar), !is.na(estado))

if (!is.null(TEST_UF)) {
  stations <- stations %>% filter(estado == TEST_UF)
}

if (nrow(stations) == 0) {
  stop("Nenhum posto em df_hydroyear_uniplu.rds.")
}

uf_lab <- if (!is.null(TEST_UF)) TEST_UF else "Brasil"

g_all <- plot_fig3_postos(
  stations,
  paste0("Ocorrência dos eventos extremos de precipitação — ", uf_lab)
)
salvar_fig3(g_all, file.path(DIR_PIC, "fig3_postos_uniplu.png"))
message("Salvo fig3_postos_uniplu.png (", nrow(stations), " postos)")

ufs <- sort(unique(stations$estado))
for (uf in ufs) {
  df_uf <- stations %>% filter(estado == uf)
  g_uf <- plot_fig3_postos(
    df_uf,
    paste0("Ocorrência dos eventos extremos de precipitação — ", uf)
  )
  salvar_fig3(g_uf, file.path(DIR_POSTOS_UF, paste0(uf, ".png")))
  message("Salvo fig3_postos_uf/", uf, ".png (", nrow(df_uf), " postos)")
}
