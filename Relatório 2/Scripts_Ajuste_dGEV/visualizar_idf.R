library(ggplot2)
library(patchwork)

font <- "sans"
fontsize <- 8
textsize = fontsize/.pt
set_theme(
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key.height = unit(fontsize, "pt"),
    strip.background = element_blank(),
    text = element_text(size = 8),
    panel.spacing.y = unit(0.5, "pt"),
  )
)

gg <- "A801"
idf <- idf.df |> filter(gauge_code == gg)
par.gg <- par.df |> filter(gauge_code == gg)

idf.plot <- idf |> 
  filter(gauge_code == gg, return_period %in% c(5, 25, 200)) |> 
  ggplot(aes(x = d, color = factor(return_period), fill = factor(return_period))) +
  facet_wrap(~gauge_code) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linewidth = 0, show.legend = FALSE) +
  geom_line(aes(y = q), linewidth = 0.3) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("ξ = %7.3f", -par)), filter(par.gg, parname == "xi"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", par)), filter(par.gg, parname == "H"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + textsize*0.7, family = font, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ = %7.3f", par)), filter(par.gg, parname == "theta"), inherit.aes = FALSE,
    hjust = 1.1, vjust = 2.3 + 2*textsize*0.7, family = font, size = textsize
  ) +
  scale_x_continuous(
    name = "Tempo [h]",
    trans = "log10"
  ) +
  scale_y_continuous(
    name = "Intensidade [mm/h]",
    trans = "log10"
  ) +
  scale_color_grey(
    name = "Tempo de retorno", start = 0.9, end = 0.4, 
    labels = as_labeller(function(x) paste(x, "anos")),
    aesthetics = list("color", "fill")
  ); idf.plot


df <- imax.df |> filter(gauge_code == gg)

imax.plot <- df |> 
  filter(d %in% c(1, 6, 12, 24)) |> 
  ggplot(aes(x = wateryear, y = imax)) +
  facet_wrap(~factor(d), ncol = 1, scales = "free_y", labeller = as_labeller(function(x) paste(x, "h"))) +
  geom_line(linewidth = 0.3) +
  geom_point(shape = 21, color = "black", fill = "white", size = 0.55) +
  labs(x = NULL, y = NULL); imax.plot

idf.plot + imax.plot + plot_layout(nrow = 1, widths = c(1, 0.7))

ggsave(sprintf("C:/Users/Adm/Downloads/idf_series_%s.png", gg), width = 16, height = 9, units = "cm", dpi = 600)
ggsave("path", plot, width = 15, height = 8, units = "cm", dpi = 600)

# Ajuste individual `270210801A`
# Fit
i.jd <- imax.df |> 
  filter(gauge_code == gg) |> 
  pivot_wider(
    names_from = d,
    values_from = imax,
    names_prefix = "d_"
  ) |> 
  arrange(wateryear) |> 
  select(-c(wateryear, gauge_code)) |> 
  na.omit() |> 
  as.matrix()

ds <- as.double(gsub("d_", "", colnames(i.jd)))

fit <- fdgev(
  i.jd = i.jd,
  prior.info = prior,
  cols = c("gauge_code", "wateryear", "d", "imax"),
  method = "L-BFGS-B"
)

par <- fit$par
npar <- length(par)
convergence <- fit$convergence
max.ll <- -fit$value

cipar <- ciprof.par(i.jd, ds, par, max.ll, 0.95, 0.005, 0.005, gg, prior)
ciq <- ciprof.q(i.jd, ds, durations, return.period, par, max.ll, 0.95, 0.005, 0.005, gg, prior)

profile.par <- cipar$profile

profile.par |> 
  filter(parname == "xi") |> 
  ggplot(aes(x = par, y = loglik)) +
  geom_line()

profile.q <- ciq$profile


profile.q |> 
  filter(return_period == 200) |> 
  ggplot(aes(x = q, y = loglik)) +
  facet_wrap(~factor(d), scales = "free_x") +
  geom_line()

library(ggplot2)
library(ggspatial)
library(geobr)
library(rnaturalearth)
library(sf)

font <- "sans"
fontsize <- 8
textsize = fontsize/.pt
set_theme(
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key.height = unit(fontsize, "pt"),
    strip.background = element_blank(),
    text = element_text(size = 8),
    panel.spacing.y = unit(0.5, "pt"),
  )
)

br.sf <- geobr::read_state(year = 2025)
world.sf <- rnaturalearth::ne_countries(continent = "south america", scale = "large")

bbox.br <- st_bbox(br.sf)
bbox.br[3] <- -35

ggplot() +
  geom_sf(data = world.sf, color = "grey30", fill = "grey90") +
  geom_sf(data = br.sf, color = "grey10", fill = "lightgreen", linewidth = 0.2) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
    bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
    style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  coord_sf(crs = st_crs(4674), xlim = bbox.br[c(1,3)], ylim = bbox.br[c(2,4)]) +
  theme(plot.margin = margin(t = 10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.key.spacing.x = unit(fontsize*2, "pt"),
        legend.margin = margin(t = -5, l = -2),
        legend.spacing.x = unit(-4, "pt"),
        legend.position = "bottom")