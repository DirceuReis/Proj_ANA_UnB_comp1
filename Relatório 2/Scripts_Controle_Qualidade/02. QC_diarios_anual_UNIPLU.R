## Controle de Qualidade (QC) - ANUAL - UNIPLU-BR
##
## Insumos: "series_diarias_1440.parquet" e "inventario_1440.csv"
##
## Pré-requisito: rodar "00a. Leitura_Series_UNIPLU.R" e, em seguida,
## "01. Series_Diarias_1440_UNIPLU.R".

##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","lubridate","tidyr","data.table","ggplot2","patchwork","writexl"))

##### 1. Library imports -------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(lubridate)
library(tidyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(writexl)
library(sf)
library(geobr)
library(scales)
library(stringr)
library(ggspatial)

##### 2. Caminhos --------------------------------------------------------------
entrada_diarias  <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/01. Series_Diarias"
saida         <- "C:/Users/daniele.silva/OneDrive - RHAMA CONSULTORIA AMBIENTAL LTDA EPP/Área de Trabalho/UNB/R/UNIPLU/02a. QC_diarios_anual"
dir.create(saida, showWarnings = FALSE, recursive = TRUE)

# série diária bruta (sem QC) e inventário (filtrado por Hidroweb diário e INMET diario)
df_diario <- read_parquet(file.path(entrada_diarias, "series_diarias_1440.parquet"))
inventario_diario <- read_csv(
  file.path(entrada_diarias, "inventario_diario.csv"),
  col_types = cols(gauge_code = col_character(), .default = col_guess())
)

cat("Estações com series diarias   :", nrow(inventario_diario), "\n")
cat("Registros diários", nrow(df_diario), "\n")

##### 3. Filtro de série mínima (>= 30 anos)  ----------------------------------
LIMIAR_ANOS <- 30

estacoes_excluidas_anos <- inventario_diario %>%
  filter(is.na(n_anos) | n_anos < LIMIAR_ANOS) %>%
  arrange(desc(n_anos))

inventario_diario_final <- inventario_diario %>%
  filter(!is.na(n_anos) & n_anos >= LIMIAR_ANOS)

df_diario <- df_diario %>% filter(gauge_code %in% inventario_diario_final$gauge_code)

write_xlsx(
  estacoes_excluidas_anos %>% mutate(across(where(is.Date), as.character)),
  file.path(saida, "estacoes_excluidas_menos_30anos.xlsx")
)

cat(sprintf(
  "Estações excluídas (< %d anos de span ou n_anos ausente): %d de %d (%.1f%%)\n",
  LIMIAR_ANOS, nrow(estacoes_excluidas_anos), nrow(inventario_diario),
  100 * nrow(estacoes_excluidas_anos) / nrow(inventario_diario)
))

cat("Estações remanescentes (span >=", LIMIAR_ANOS, "anos):", nrow(inventario_diario_final), "\n")


## Mapa estações de séries diárias - para QC -----------------------------------
brasil_sf <- geobr::read_state(year = 2020, showProgress = FALSE)
brasil_sf$area_km2 <- as.numeric(sf::st_area(brasil_sf)) / 1e6
mapa_dados <- inventario_diario_final %>% filter(!is.na(lat), !is.na(long))
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
                   text_pad   = unit(0.2, "cm")) +
  annotation_north_arrow(
    location = "bl",
    pad_y  = unit(0.8, "cm"),
    height = unit(2.0, "cm"),
    width  = unit(2.0, "cm"),
    style  = north_arrow_fancy_orienteering(
      text_size  = 22,
      line_width = 1.2,
      text_face  = "bold")) +
  coord_sf(xlim = c(min(mapa_dados$long) - 0.5, max(mapa_dados$long) + 0.5),
           ylim = c(min(mapa_dados$lat)  - 0.5, max(mapa_dados$lat)  + 0.5))

ggsave(file.path(saida, "Figura 6. Inventário Diário (30 anos ou mais).png"), Figura6,
       width = 16, height = 12, units = "in", dpi = 300)

cat("Figura 6. Inventário Diário (30 anos ou mais).png salvo.\n")

##### 4. QC Básico -------------------------------------------------------------

LIMITE_SUPERIOR_MM <- 700

n_antes <- nrow(df_diario)

registros_excluidos_basico <- df_diario %>%
  filter(rain_mm < 0 | rain_mm > LIMITE_SUPERIOR_MM) %>%
  mutate(motivo = ifelse(rain_mm < 0, "negativo", paste0("> ", LIMITE_SUPERIOR_MM, " mm/dia"))) %>%
  left_join(inventario_diario_final %>% select(gauge_code, city, state), by = "gauge_code") %>%
  relocate(city, state, .after = gauge_code)

df_diario <- df_diario %>%
  mutate(rain_mm = ifelse(rain_mm < 0 | rain_mm > LIMITE_SUPERIOR_MM, NA_real_, rain_mm)) %>%
  filter(!is.na(rain_mm))

n_removido_basico <- nrow(registros_excluidos_basico)
n_restante_basico  <- nrow(df_diario)

write_xlsx(registros_excluidos_basico %>% mutate(date = format(date, "%Y-%m-%d")),
           file.path(saida, "qc_basico_registros_excluidos.xlsx"))

cat(sprintf(
  "\nQC Básico: %d de %d registros diários descartados (%.4f%% do total) por serem negativos ou > %d mm/dia.\n",
  n_removido_basico, n_antes, 100 * n_removido_basico / n_antes, LIMITE_SUPERIOR_MM
))
cat(sprintf("Registros diários restantes após o QC Básico: %d\n", n_restante_basico))

resumo_basico <- tibble(
  categoria = c("Válido", "Negativo", paste0("> ", LIMITE_SUPERIOR_MM, " mm/dia")),
  n = c(
    n_restante_basico,
    sum(registros_excluidos_basico$motivo == "negativo"),
    sum(registros_excluidos_basico$motivo != "negativo")
  )
) %>%
  mutate(percentual = round(100 * n / n_antes, 4))

write_xlsx(resumo_basico, file.path(saida, "qc_basico_resumo.xlsx"))

p_basico <- ggplot(resumo_basico, aes(x = categoria, y = pmax(n, 1))) +
  geom_col(fill = "grey50") +
  geom_text(aes(label = scales::comma(n, big.mark = ".", decimal.mark = ",")), vjust = -0.4, size = 3.5) +
  scale_y_log10(labels = function(x) scales::comma(x, big.mark = ".", decimal.mark = ",")) +
  labs(
    title = "Controle de Qualidade Básico - Estações Diárias (QC Anual)",
    subtitle = paste0("Limite físico: ", LIMITE_SUPERIOR_MM, " mm/dia"),
    x = NULL, y = "Número de registros diários (escala log)"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(saida, "fig_qc_basico_resumo.png"), p_basico, width = 7, height = 5, dpi = 600)
cat("Resumo do QC Básico (dados + gráfico) salvo em:", saida, "\n")

########################################################################

##### 5. QC Absoluto --------------------------------------------------------------------

## P - Disponibilidade de dados ------------------------------------------------------

dt <- as.data.table(df_diario)[, .(gauge_code, date, rain_mm)] #organiza os dados
dt[, ano := year(date)]

n_dup <- anyDuplicated(dt, by = c("gauge_code", "date")) #confere se há dias repetidos
if (n_dup > 0L) stop("Há estação-dia repetido na série diária", n_dup, ").")

anos_esperados <- inventario_diario_final %>% #lista os anos esperados com dado
  filter(!is.na(data_inicio), !is.na(data_fim)) %>%
  rowwise() %>%
  mutate(periodo = list(seq(year(as.Date(data_inicio)), year(as.Date(data_fim))))) %>%
  ungroup() %>%
  select(gauge_code, periodo) %>%
  unnest(periodo) %>%
  transmute(gauge_code, ano = periodo)

dias_com_dado <- dt[, .(n_dias_com_dado = .N), by = .(gauge_code, ano)] #dias com dado em cada estação-ano

qc_P <- merge(as.data.table(anos_esperados), dias_com_dado, by = c("gauge_code", "ano"), all.x = TRUE)
qc_P[is.na(n_dias_com_dado), n_dias_com_dado := 0L] # ano sem dado = 0 dias
qc_P[, n_dias_ano := fifelse((ano %% 4 == 0 & ano %% 100 != 0) | ano %% 400 == 0, 366L, 365L)]  # bissexto?
qc_P[, P := 100 * n_dias_com_dado / n_dias_ano]
if (any(qc_P$P < 0 | qc_P$P > 100)) stop("P fora do intervalo 0-100%.")

qc_P[, faixa_P := fcase(P == 0,  "0% (ano sem dado)",
                        P < 50,  "0-50%",
                        P < 80,  "50-80%",
                        P < 90,  "80-90%",
                        P < 95,  "90-95%",
                        P < 99,  "95-99%",
                        default = "99-100%")]

resumo_P <- qc_P[, .(estacao_anos = .N), by = faixa_P][
  order(match(faixa_P, c("0% (ano sem dado)", "0-50%", "50-80%", "80-90%", "90-95%", "95-99%", "99-100%")))][
    , pct := round(100 * estacao_anos / sum(estacao_anos), 2)]

write_xlsx(resumo_P, file.path(saida, "qc_P_resumo.xlsx"))
print(resumo_P)

## Q1 - Falhas --------------------------------------------------------------------------------

dt[, doy := yday(date)]   # dia do ano (1 a 365/366)

# maior sequência de dias sem dado, em cada estação-ano com algum dado
maior_falha_ano <- dt[, {
  n <- if ((ano[1] %% 4 == 0 & ano[1] %% 100 != 0) | ano[1] %% 400 == 0) 366L else 365L
  .(maior_falha = max(diff(c(0L, sort(doy), n + 1L)) - 1L))
}, by = .(gauge_code, ano)]

qc_Q1 <- merge(qc_P[, .(gauge_code, ano, n_dias_ano, n_dias_com_dado)], maior_falha_ano,
               by = c("gauge_code", "ano"), all.x = TRUE)
qc_Q1[is.na(maior_falha), maior_falha := n_dias_ano]   # ano sem dado: a falha é o ano inteiro
qc_Q1[, n_dias_sem_dado := n_dias_ano - n_dias_com_dado]
qc_Q1[, Q1 := 100 - 100 * (2 * n_dias_sem_dado + maior_falha) / n_dias_ano]

if (any(qc_Q1$Q1 < -200 | qc_Q1$Q1 > 100)) stop("Q1 fora do intervalo -200 a 100.")
if (any(qc_Q1$maior_falha > qc_Q1$n_dias_sem_dado)) stop("Maior falha maior que o total de dias sem dado.")

qc_Q1[, faixa_Q1 := fcase(Q1 == 100, "100 (sem falha)",
                          Q1 >= 90,  "90-100",
                          Q1 >= 80,  "80-90",
                          Q1 >= 50,  "50-80",
                          Q1 >= 0,   "0-50",
                          Q1 > -200, "negativo",
                          default = "-200 (ano sem dado)")]
resumo_Q1 <- qc_Q1[, .(estacao_anos = .N), by = faixa_Q1][
  order(match(faixa_Q1, c("-200 (ano sem dado)", "negativo", "0-50", "50-80", "80-90", "90-100", "100 (sem falha)")))][
    , pct := round(100 * estacao_anos / sum(estacao_anos), 2)]

write_xlsx(resumo_Q1, file.path(saida, "qc_Q1_resumo.xlsx"))
print(resumo_Q1)

## Q2 - Dia da semana --------------------------------------------------------------------------

LIMIAR_DIA_CHUVOSO <- 0.2

dt[, dow := wday(date)]

q2_ano <- dt[, {
  cont <- tabulate(dow[rain_mm >= LIMIAR_DIA_CHUVOSO], nbins = 7)   # dias chuvosos por dia da semana
  .(n_dias_chuvosos = sum(cont),
    cv_dia_semana   = if (mean(cont) > 0) sd(cont) / mean(cont) else 0)
}, by = .(gauge_code, ano)]

qc_Q2 <- merge(qc_P[, .(gauge_code, ano)], q2_ano, by = c("gauge_code", "ano"), all.x = TRUE)
qc_Q2[is.na(cv_dia_semana), `:=`(n_dias_chuvosos = 0L, cv_dia_semana = 0)]   # ano sem dado
qc_Q2[, Q2 := 100 - 100 * cv_dia_semana]

qc_Q2[, faixa_Q2 := fcase(Q2 >= 90, "90-100", Q2 >= 80, "80-90", Q2 >= 70, "70-80",
                          Q2 >= 50, "50-70", Q2 >= 0,  "0-50",  default = "negativo")]

resumo_Q2 <- qc_Q2[, .(estacao_anos = .N), by = faixa_Q2][
  order(match(faixa_Q2, c("negativo", "0-50", "50-70", "70-80", "80-90", "90-100")))][
    , pct := round(100 * estacao_anos / sum(estacao_anos), 2)]

write_xlsx(resumo_Q2, file.path(saida, "qc_Q2_resumo.xlsx"))
print(resumo_Q2)

## Q3 - Outliers -------------------------------------------------------------------------------

dt[, mes := month(date)]

limiar_outlier <- dt[rain_mm > 0,
                     .(q1 = quantile(rain_mm, 0.25, names = FALSE),
                       q3 = quantile(rain_mm, 0.75, names = FALSE)),
                     by = .(gauge_code, mes)]    

limiar_outlier[, limiar := q3 + 1.5 * (q3 - q1)]

dt[limiar_outlier, on = .(gauge_code, mes), limiar := i.limiar]
dt[, outlier := !is.na(limiar) & rain_mm > limiar]

n_outlier_ano <- dt[, .(n_outlier = sum(outlier)), by = .(gauge_code, ano)]

qc_Q3 <- merge(qc_P[, .(gauge_code, ano, n_dias_ano)], n_outlier_ano,
               by = c("gauge_code", "ano"), all.x = TRUE)
qc_Q3[is.na(n_outlier), n_outlier := 0L]                 # ano sem dado
qc_Q3[, Q3 := 100 * (n_dias_ano - n_outlier) / n_dias_ano]
if (any(qc_Q3$Q3 < 0 | qc_Q3$Q3 > 100)) stop("Q3 fora do intervalo 0-100.")

resumo_Q3 <- qc_Q3[, .(estacao_anos = .N, outliers_mediana = median(n_outlier),
                       outliers_p95 = quantile(n_outlier, 0.95), Q3_mediano = median(Q3))]

write_xlsx(resumo_Q3, file.path(saida, "qc_Q3_resumo.xlsx"))
print(resumo_Q3)

## Q - Índice de qualidade e classificação -----------------------------------------------------

# dias com chuva > 0 em cada estação-ano (para identificar anos falhos)
dias_chuva_positiva <- dt[, .(n_dias_chuva_positiva = sum(rain_mm > 0)), by = .(gauge_code, ano)]

# junta as quatro notas numa tabela só (uma linha por estação-ano)
qc_anual <- Reduce(function(a, b) merge(a, b, by = c("gauge_code", "ano"), all.x = TRUE),
                   list(qc_P[,  .(gauge_code, ano, n_dias_ano, n_dias_com_dado, P)],
                        qc_Q1[, .(gauge_code, ano, n_dias_sem_dado, maior_falha, Q1)],
                        qc_Q2[, .(gauge_code, ano, n_dias_chuvosos, Q2)],
                        qc_Q3[, .(gauge_code, ano, n_outlier, Q3)],
                        dias_chuva_positiva))

qc_anual[is.na(n_dias_chuva_positiva), n_dias_chuva_positiva := 0L]   # ano sem dado

qc_anual[, Q := (P + Q1 + Q2 + Q3) / 4]

qc_anual[, ano_falho := n_dias_com_dado >= 0.9 * n_dias_ano & n_dias_chuva_positiva == 0]

qc_anual[, rotulo := fcase(
  ano_falho,    "Muito Baixa",
  P >= 99 & Q >= 90, "Excelente",
  P >= 95 & Q >= 85, "Boa",
  P >= 90 & Q >= 80, "Aceitável",
  Q >= 50,           "Baixa",
  default =          "Muito Baixa")]

qc_anual[, qualidade_absoluta := fifelse(rotulo %in% c("Excelente", "Boa", "Aceitável"),
                                         "Alta Qualidade", "Baixa Qualidade")]
setorder(qc_anual, gauge_code, ano)

if (anyDuplicated(qc_anual, by = c("gauge_code", "ano")) > 0L) stop("Estação-ano repetido em qc_anual.")
if (anyNA(qc_anual[, .(P, Q1, Q2, Q3, Q)])) stop("Há nota faltando em qc_anual.")

resumo_rotulo <- qc_anual[, .(estacao_anos = .N), by = .(rotulo, qualidade_absoluta)][
  order(match(rotulo, c("Excelente", "Boa", "Aceitável", "Baixa", "Muito Baixa")))][
    , pct := round(100 * estacao_anos / sum(estacao_anos), 2)]

cat(sprintf("\nAnos degenerados (>= 90%% de cobertura, 0 dias de chuva): %d\n", sum(qc_anual$ano_degenerado)))
print(resumo_rotulo)

fwrite(qc_anual, file.path(saida, "controle_qualidade_anual.csv"))
write_xlsx(resumo_rotulo, file.path(saida, "classificacao_resumo.xlsx"))

##### 6. Série diária (só anos de Alta Qualidade) --------------------------------------

anos_HQ <- qc_anual[qualidade_absoluta == "Alta Qualidade", .(gauge_code, ano)]

serie_QC <- dt[anos_HQ, on = .(gauge_code, ano), nomatch = NULL, .(gauge_code, date, rain_mm)]
setorder(serie_QC, gauge_code, date)

write_parquet(serie_QC, file.path(saida, "series_diarias_QC.parquet"))

n_qc_removido <- nrow(dt) - nrow(serie_QC)
cat(sprintf(
  "QC Absoluto: %d de %d dias (%.2f%%) removidos por pertencerem a anos de Baixa Qualidade.\n",
  n_qc_removido, nrow(dt), 100 * n_qc_removido / nrow(dt)))
cat(sprintf("Série filtrada: %d dias | %d estação-anos | %d estações com ao menos 1 ano aprovado\n",
            nrow(serie_QC), nrow(anos_HQ), uniqueN(serie_QC$gauge_code)))

##### 7. Filtro final: estações com >= 30 anos de HQ -------------------------------

LIMIAR_ANOS_BONS <- 30

anos_bons_por_estacao <- qc_anual[qualidade_absoluta == "Alta Qualidade",
                                  .(n_anos_bons      = .N,
                                    primeiro_ano_bom = min(ano),
                                    ultimo_ano_bom   = max(ano)),
                                  by = gauge_code]

estacoes_avaliadas <- as.data.table(inventario_diario_final)[
  , .(gauge_code, city, state, lat, long, elevation, network, responsible,
      data_inicio, data_fim, n_anos)]

estacoes_avaliadas <- merge(estacoes_avaliadas, anos_bons_por_estacao, by = "gauge_code", all.x = TRUE)
estacoes_avaliadas[is.na(n_anos_bons), n_anos_bons := 0L]   # estação sem nenhum ano bom

estacoes_finais    <- estacoes_avaliadas[n_anos_bons >= LIMIAR_ANOS_BONS][order(-n_anos_bons)]
estacoes_excluidas <- estacoes_avaliadas[n_anos_bons <  LIMIAR_ANOS_BONS][order(-n_anos_bons)]

write_xlsx(estacoes_finais %>% mutate(across(where(is.Date), as.character)),
           file.path(saida, "estacoes_finais_30anos_altaqualidade.xlsx"))
write_xlsx(estacoes_excluidas %>% mutate(across(where(is.Date), as.character)),
           file.path(saida, "estacoes_excluidas_menos_30anos_bons.xlsx"))

cat(sprintf("\nEstações avaliadas no QC (span >= %d anos): %d\n", LIMIAR_ANOS, nrow(estacoes_avaliadas)))
cat(sprintf("Excluídas (< %d anos de Alta Qualidade)    : %d\n", LIMIAR_ANOS_BONS, nrow(estacoes_excluidas)))
cat(sprintf("ESTAÇÕES FINAIS (>= %d anos bons)           : %d\n", LIMIAR_ANOS_BONS, nrow(estacoes_finais)))
cat("\nDistribuição por UF:\n")
print(estacoes_finais[, .N, by = state][order(-N)])

##### 8. Série diária final: estações finais e anos de Alta Qualidade --------------------------

serie_diaria_final <- serie_QC[gauge_code %in% estacoes_finais$gauge_code]

serie_diaria_final <- merge(serie_diaria_final,
                     estacoes_finais[, .(gauge_code, city, state, lat, long, elevation, network, responsible)],
                     by = "gauge_code")

setcolorder(serie_diaria_final, c("gauge_code", "city", "state", "lat", "long", "elevation",
                           "network", "responsible", "date", "rain_mm"))

setorder(serie_diaria_final, gauge_code, date)

write_parquet(serie_diaria_final, file.path(saida, "series_diarias_final.parquet"))

cat(sprintf("Série final: %d dias | %d estações | %d estação-anos\n",
            nrow(serie_diaria_final), uniqueN(serie_diaria_final$gauge_code),
            uniqueN(serie_diaria_final[, .(gauge_code, year(date))])))

##### 9. Resultados por variável: série anual nacional ------------------------------------------

##### 6. Figuras - P, Q1, Q2 e Q3 ao longo do tempo --------------------------------------

painel_serie <- function(v, rotulo) {
  ggplot(qc_anual_medio, aes(x = ano, y = .data[[v]])) +
    geom_line(color = "black", linewidth = 0.3) +
    geom_point(color = "black", shape = 4, size = 0.9) +
    scale_x_continuous(breaks = seq(1850, 2030, by = 10)) +
    coord_cartesian(ylim = c(NA, 100)) +            # teto em 100, piso livre
    labs(x = "Ano", y = paste0(v, " (%)"), title = rotulo) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title  = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 14, color = "black", angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title  = element_text(size = 14, color = "black"))
}

Figura7 <- (painel_serie("P",  "(a)") |
                    painel_serie("Q1", "(b)")) /
  (painel_serie("Q2", "(c)") |
     painel_serie("Q3", "(d)"))

ggsave(file.path(saida, "Figura 7. Series anuais P Q1 Q2 Q3.png"), Figura7,
       width = 14, height = 8, units = "in", dpi = 300)

