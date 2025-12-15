# Script preliminar, preguiça de fazer 'leitura' e 'descricao' agora

# PACOTES -----------------------------------------------------------------

rm(list = ls()); invisible(gc())

if(!require(pacman)) install.packates("pacman")
pacman::p_load(pacman,
               arrow,
               tidyverse,
               lubridate,
               patchwork,
               ggplot2,
               ggtext,
               gghighlight,
               readxl,
               ggforce)


# FUNÇÕES -----------------------------------------------------------------

# Regressão linear e análise invariância de escala
source("scripts/funcoes/fun_scale_invariance.R")


# LER DADOS ---------------------------------------------------------------

# Dados diários

# Escolher colunas