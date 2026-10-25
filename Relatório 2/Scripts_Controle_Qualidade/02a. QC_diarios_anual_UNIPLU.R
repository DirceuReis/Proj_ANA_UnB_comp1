## Controle de Qualidade Diárias (QC) - ANUAL - UNIPLU-BR
##
## Verifica, por estação-ano, a série que o 01a já selecionou por disponibilidade e número de anos.
## Quatro índices - P, Q1, Q2 e Q3. Filtra por qualidade.
##
## Insumos: do 01a, "series_diarias_filtrada.parquet" e "inventario_diario.csv";
##          do geobr, os contornos de UF e país de 2025, baixados da internet.
## Saídas:  "series_diarias_QC.parquet" (anos de Alta Qualidade),
##          "series_diarias_final.parquet" (estações com >= 30 anos bons),
##          as Figuras 6 e 7 e as planilhas de resumo do QC.
##
##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","lubridate","data.table","ggplot2",
#                   "patchwork","writexl","sf","geobr","scales","ggspatial"))

##### 1. Library imports -------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(lubridate)
library(data.table)
library(ggplot2)
library(patchwork)
library(writexl)
library(sf)
library(geobr)
library(scales)
library(ggspatial)

##### 2. Caminhos --------------------------------------------------------------
entrada_diarias <- "C:/UNB/UNIPLU/01a. Series_Diarias"
saida           <- "C:/UNB/UNIPLU/02a. QC_diarios_anual"
dir.create(saida, showWarnings = FALSE, recursive = TRUE)

inventario_diario <- read_csv(
  file.path(entrada_diarias, "inventario_diario.csv"),
  col_types = cols(gauge_code = col_character(), .default = col_guess())
)

cat("Estações no inventário diário :", nrow(inventario_diario), "\n")

datas_para_texto <- function(d)
  dplyr::mutate(d, dplyr::across(
    dplyr::where(~ inherits(.x, "Date") | inherits(.x, "POSIXct")), as.character))

##### 3. Padrão gráfico --------------------------------------------------------
font     <- "sans"
fontsize <- 8
textsize <- fontsize / .pt     
epsg     <- 5880               

LARG_MAPA <- 16                
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

## contorno territorial e enquadramento, comuns às Figuras 6 e 8
uf.sf   <- geobr::read_state(year = 2025, showProgress = FALSE)
br.sf   <- geobr::read_country(year = 2025, showProgress = FALSE)
bbox.br <- sf::st_bbox(br.sf)

enquadra <- coord_sf(crs = st_crs(epsg), default_crs = st_crs(4674),
                     xlim = bbox.br[c(1, 3)], ylim = bbox.br[c(2, 4)])

rede_cores <- c(
  "Hidroweb diário" = "#2e7d32",
  "INMET diário"    = "#6a1b9a")

rede_formas <- c(
  "Hidroweb diário" = 17,  # triângulo
  "INMET diário"    = 18)  # losango

##### 4. Filtro de série mínima (>= 30 anos)  ----------------------------------
LIMIAR_ANOS <- 30

estacoes_excluidas_anos <- inventario_diario %>%
  filter(is.na(n_anos_medido) | n_anos_medido < LIMIAR_ANOS) %>%
  arrange(desc(n_anos_medido))

## A série já selecionada pelo 01a
df_diario <- open_dataset(file.path(entrada_diarias, "series_diarias_filtrada.parquet")) %>%
  collect()
setDT(df_diario)

inventario_diario_final <- inventario_diario %>%
  filter(!is.na(n_anos_medido) & n_anos_medido >= LIMIAR_ANOS,
         gauge_code %in% unique(df_diario$gauge_code))

cat("Registros diários lidos       :", format(nrow(df_diario), big.mark = "."), "\n")

write_xlsx(
  estacoes_excluidas_anos %>% datas_para_texto(),
  file.path(saida, "estacoes_excluidas_menos_30anos.xlsx")
)

cat(sprintf(
  "Estações excluídas (< %d anos de span ou n_anos ausente): %d de %d (%.1f%%)\n",
  LIMIAR_ANOS, nrow(estacoes_excluidas_anos), nrow(inventario_diario),
  100 * nrow(estacoes_excluidas_anos) / nrow(inventario_diario)
))
cat(sprintf("Estações na série selecionada pelo 01a                   : %d\n",
            uniqueN(df_diario$gauge_code)))
cat(sprintf("Estações avaliadas no QC (span >= %d anos e na série)     : %d\n",
            LIMIAR_ANOS, nrow(inventario_diario_final)))

##### 5. QC Básico -------------------------------------------------------------
LIMITE_SUPERIOR_MM <- 700

n_antes <- nrow(df_diario)

registros_excluidos_basico <- df_diario[rain_mm < 0 | rain_mm > LIMITE_SUPERIOR_MM] %>%
  mutate(motivo = ifelse(rain_mm < 0, "negativo", paste0("> ", LIMITE_SUPERIOR_MM, " mm/dia"))) %>%
  left_join(inventario_diario_final %>% select(gauge_code, city, state), by = "gauge_code") %>%
  relocate(city, state, .after = gauge_code)

df_diario <- df_diario[!is.na(rain_mm) & rain_mm >= 0 & rain_mm <= LIMITE_SUPERIOR_MM]

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
  geom_text(aes(label = scales::comma(n, big.mark = ".", decimal.mark = ",")),
            vjust = -0.4, size = fontsize * 0.9, size.unit = "pt", family = font) +
  scale_y_log10(labels = function(x) scales::comma(x, big.mark = ".", decimal.mark = ",")) +
  labs(
    title = "Controle de Qualidade Básico - Estações Diárias (QC Anual)",
    subtitle = paste0("Limite físico: ", LIMITE_SUPERIOR_MM, " mm/dia"),
    x = NULL, y = "Número de registros diários (escala log)"
  )

ggsave(file.path(saida, "fig_qc_basico_resumo.png"), p_basico,
       width = LARG_GRAF, height = ALT_GRAF, units = "cm", dpi = 300, bg = "white")
cat("Resumo do QC Básico (dados + gráfico) salvo em:", saida, "\n")

##### 6. QC Absoluto --------------------------------------------------------------------

## P - Disponibilidade de dados ------------------------------------------------------
dt <- df_diario
rm(df_diario); invisible(gc())

dt[, ano := year(date)]

n_dup <- anyDuplicated(dt, by = c("gauge_code", "date")) #confere se há dias repetidos
if (n_dup > 0L) stop("Há estação-dia repetido na série diária", n_dup, ").")

## Grade de anos medidas pela série
anos_esperados <- unique(dt[, .(gauge_code, ano)])

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
  P >= 95 & Q >= 90, "Excelente",
  P >= 90 & Q >= 85, "Boa",
  P >= 80 & Q >= 80, "Aceitável",
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

cat(sprintf("\nAnos sem chuva (>= 90%% de cobertura, 0 dias de chuva): %d\n", sum(qc_anual$ano_falho)))

fwrite(qc_anual, file.path(saida, "controle_qualidade_anual.csv"))
write_xlsx(resumo_rotulo, file.path(saida, "classificacao_resumo.xlsx"))

# média nacional por ano
qc_anual_medio <- qc_anual[, .(n_estacoes = uniqueN(gauge_code),
                               P = mean(P), Q1 = mean(Q1), Q2 = mean(Q2),
                               Q3 = mean(Q3), Q = mean(Q)), by = ano][order(ano)]

write_xlsx(qc_anual_medio, file.path(saida, "qc_serie_anual_nacional.xlsx"))

##### 7. Série diária (só anos de Alta Qualidade) --------------------------------------
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

##### 8. Filtro final: estações com >= 30 anos de HQ -------------------------------
LIMIAR_ANOS_BONS <- 30

anos_bons_por_estacao <- qc_anual[qualidade_absoluta == "Alta Qualidade",
                                  .(n_anos_bons      = .N,
                                    primeiro_ano_bom = min(ano),
                                    ultimo_ano_bom   = max(ano)),
                                  by = gauge_code]

estacoes_avaliadas <- as.data.table(inventario_diario_final)[
  , .(gauge_code, city, state, lat, long, elevation, network, responsible,
      data_inicio, data_fim, n_anos_medido)]

estacoes_avaliadas <- merge(estacoes_avaliadas, anos_bons_por_estacao, by = "gauge_code", all.x = TRUE)
estacoes_avaliadas[is.na(n_anos_bons), n_anos_bons := 0L]   # estação sem nenhum ano bom

estacoes_finais    <- estacoes_avaliadas[n_anos_bons >= LIMIAR_ANOS_BONS][order(-n_anos_bons)]
estacoes_excluidas <- estacoes_avaliadas[n_anos_bons <  LIMIAR_ANOS_BONS][order(-n_anos_bons)]

write_xlsx(estacoes_finais %>% datas_para_texto(),
           file.path(saida, "estacoes_finais_30anos_altaqualidade.xlsx"))
write_xlsx(estacoes_excluidas %>% datas_para_texto(),
           file.path(saida, "estacoes_excluidas_menos_30anos_bons.xlsx"))

cat(sprintf("\nEstações avaliadas no QC (span >= %d anos): %d\n", LIMIAR_ANOS, nrow(estacoes_avaliadas)))
cat(sprintf("Excluídas (< %d anos de Alta Qualidade)    : %d\n", LIMIAR_ANOS_BONS, nrow(estacoes_excluidas)))
cat(sprintf("ESTAÇÕES FINAIS (>= %d anos bons)           : %d\n", LIMIAR_ANOS_BONS, nrow(estacoes_finais)))
cat("\nDistribuição por UF:\n")
print(estacoes_finais[, .N, by = state][order(-N)])

##### 9. Série diária final: estações finais e anos de Alta Qualidade --------------------------
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

##### 10. Figura 6. P, Q1, Q2 e Q3 ao longo do tempo --------------------------------------
painel_indice <- function(v, rotulo) {
  ggplot(qc_anual_medio, aes(x = ano, y = .data[[v]])) +
    geom_line(color = "black", linewidth = 0.3) +
    geom_point(color = "black", shape = 4, size = 0.9) +
    scale_x_continuous(breaks = seq(1850, 2030, by = 10)) +
    coord_cartesian(ylim = c(NA, 100)) +            # teto em 100, piso livre
    labs(x = "Ano", y = paste0(v, " (%)"), title = rotulo) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

Figura6 <- (painel_indice("P",  "(a)") |
              painel_indice("Q1", "(b)")) /
  (painel_indice("Q2", "(c)") |
     painel_indice("Q3", "(d)"))

ggsave(file.path(saida, "Figura 6. Series anuais P Q1 Q2 Q3.png"), Figura6,
       width = LARG_GRAF, height = 12, units = "cm", dpi = 300, bg = "white")

##### 11. Figura 7. Estações da série diária final --------------------------------------
mapa_finais    <- estacoes_finais %>% filter(!is.na(lat), !is.na(long))
mapa_finais_sf <- st_as_sf(mapa_finais, coords = c("long", "lat"), crs = 4674)

Figura7 <- ggplot() +
  geom_sf(data = mapa_finais_sf, aes(color = network, shape = network),
          size = 0.35, stroke = 0.08) +
  geom_sf(data = uf.sf, color = "grey30", fill = NA, linewidth = 0.1) +
  geom_sf(data = br.sf, color = "black",  fill = NA, linewidth = 0.4) +
  scale_color_manual(values = rede_cores,  name = "Rede") +
  scale_shape_manual(values = rede_formas, name = "Rede") +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  escala_mapa + norte_mapa + enquadra +
  labs(x = "Longitude", y = "Latitude") +
  tema_mapa

ggsave(file.path(saida, "Figura 7. Estacoes da serie diaria final.png"), Figura7,
       width = LARG_MAPA, height = ALT_MAPA, units = "cm", dpi = 300, bg = "white")

cat("Figura 7. Estacoes da serie diaria final.png salvo.\n")