## Controle de Qualidade Sub-Diárias - INMET subdiário - UNIPLU-BR
##
## Verifica a série que o 01b4 já selecionou por cadastro, duplicatas e
## disponibilidade.
##
## Insumos: do 01b4, "series_INMETsub_selecionada.parquet" e "inventario_INMETsub.csv";
##          do 00a, "resumo_estacao_ano.csv", pelo time_step de cada ano.
## Saídas:  "series_INMETsub_qcbasico.parquet" (série sem os FAIL),
##          "qc_basico_resumo.xlsx", "qc_basico_registros_excluidos.xlsx",
##          "qc_basico_por_estacao.csv".

##### instalar pacotes (rode só uma vez) #####
#install.packages(c("dplyr","readr","arrow","data.table","writexl"))

##### 1. Library imports ----------------------------------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(arrow)
library(data.table)
library(writexl)

##### 2. Caminhos ----------------------------------------------------------------------------------------------------------
entrada_inv  <- "C:/UNB/UNIPLU/00. Series"
entrada_rede <- "C:/UNB/UNIPLU/01b4. Series_Sub_INMETsubdiario"
saida        <- "C:/UNB/UNIPLU/02b4. QC_Sub_INMETsubdiario"

dir.create(saida, showWarnings = FALSE, recursive = TRUE)

##### 3. Parâmetros ----------------------------------------------------------------------------------------------------------
PREFIXO <- "INMETsub"

## Limiar de SUSPEITA
ENVOLTORIA_10MIN <- 140

##### 4. QC Básico ----------------------------------------------------------------------------------------------------------
## U1  rain_mm negativo e finito == VALOR NEGATIVO
## U2  rain_mm não finito: NA, NaN, Inf == VALOR NÃO FINITO
## U3  rain_mm acima do teto da duração == VALOR ACIMA DO TETO
## U4  datetime posterior ao ano do arquivo == VALOR EM ANO DIFERENTE
## U5  gauge_code fora do catálogo da rede == SEM CADASTRO DE ESTAÇÃO

arq_serie <- file.path(entrada_rede, paste0("series_", PREFIXO, "_selecionada.parquet"))
ds <- open_dataset(arq_serie)

inventario <- read_csv(file.path(entrada_rede, paste0("inventario_", PREFIXO, ".csv")),
                       col_types = cols(gauge_code = col_character(),
                                        .default   = col_guess()))

codigos_rede <- inventario$gauge_code

cat("Registros na série selecionada:",
    format(ds$num_rows, big.mark = ".", decimal.mark = ","), "\n")
cat("Estações no inventário da rede:", nrow(inventario), "\n")

## 4.1 O teto de cada estação-ano ---------------------------------------------
## Teto absoluto por duração.
FATOR_LIMITE  <- 1.5
LIMITE_PADRAO <- 800   # mm, quando o timestep do ano não é conhecido

P_envoltoria <- function(timestep_min) 41.879 * (timestep_min)^0.5245

if (abs(ENVOLTORIA_10MIN - P_envoltoria(10)) > 0.5)
  stop("ENVOLTORIA_10MIN é ", ENVOLTORIA_10MIN, " mas P_envoltoria(10) dá ",
       round(P_envoltoria(10), 1), ". Os dois têm de concordar.")

cadastro <- read_csv(file.path(entrada_inv, "resumo_estacao_ano.csv"),
                     col_select = c(gauge_code, ano, time_step),
                     col_types  = cols(gauge_code = col_character(),
                                       time_step  = col_character(),
                                       .default   = col_guess()))

limites <- cadastro %>%
  filter(gauge_code %in% codigos_rede) %>%
  transmute(gauge_code,
            ano      = as.integer(ano),
            timestep = suppressWarnings(as.numeric(time_step)),
            teto_mm  = ifelse(is.na(timestep) | timestep <= 0,
                              LIMITE_PADRAO,
                              FATOR_LIMITE * P_envoltoria(timestep)),
            limite_dt = as.POSIXct(sprintf("%d-01-03", ano + 1L), tz = "UTC"))

cat("\nTeto aplicado, por timestep declarado:\n")
print(as.data.frame(limites %>%
  count(timestep, teto_mm = round(teto_mm, 1), name = "anos_estacao") %>%
  arrange(timestep)))

## 4.2 Linhas marcadas --------------------------------------------------------------
linhas_marcadas <- ds %>%
  left_join(limites, by = c("gauge_code", "ano")) %>%
  mutate(
    U1 = !is.na(rain_mm) & is.finite(rain_mm) & rain_mm < 0,
    U2 = is.na(rain_mm)  | !is.finite(rain_mm),
    U3 = !is.na(rain_mm) & is.finite(rain_mm) &
           !is.na(teto_mm) & rain_mm > teto_mm,
    U4 = !is.na(datetime) & !is.na(limite_dt) & datetime >= limite_dt,
    U5 = is.na(gauge_code) | !(gauge_code %in% codigos_rede))

contagem <- linhas_marcadas %>%
  summarise(U1 = sum(U1), U2 = sum(U2), U3 = sum(U3),
            U4 = sum(U4), U5 = sum(U5),
            sem_datetime = sum(is.na(datetime)),
            sem_teto     = sum(is.na(teto_mm))) %>%
  collect()

cat(sprintf("\nRegistros marcados: U1 %s | U2 %s | U3 %s | U4 %s | U5 %s\n",
            format(contagem$U1, big.mark = ".", decimal.mark = ","),
            format(contagem$U2, big.mark = ".", decimal.mark = ","),
            format(contagem$U3, big.mark = ".", decimal.mark = ","),
            format(contagem$U4, big.mark = ".", decimal.mark = ","),
            format(contagem$U5, big.mark = ".", decimal.mark = ",")))
cat("Registros sem datetime, que nenhum indicador barra:",
    format(contagem$sem_datetime, big.mark = ".", decimal.mark = ","), "\n")
cat("Registros da série sem teto definido no cadastro:",
    format(contagem$sem_teto, big.mark = ".", decimal.mark = ","), "\n")

## 4.3 Registros reprovados ----------------------------------------------------
reprovados <- linhas_marcadas %>%
  filter(U1 | U2 | U3 | U4 | U5) %>%
  select(gauge_code, datetime, ano, rain_mm, state, timestep, teto_mm,
         U1, U2, U3, U4, U5) %>%
  collect()

if (nrow(reprovados) > 0) {
  reprovados <- reprovados %>%
    mutate(qc_flag = case_when(
      U1 ~ "VALOR NEGATIVO",
      U2 ~ "VALOR NÃO FINITO",
      U3 ~ "VALOR ACIMA DO TETO",
      U4 ~ "VALOR EM ANO DIFERENTE",
      U5 ~ "SEM CADASTRO DE ESTAÇÃO")) %>%
    left_join(select(inventario, gauge_code, city), by = "gauge_code") %>%
    arrange(desc(rain_mm))

  cat("\nReprovados por flag:\n")
  print(as.data.frame(reprovados %>% count(qc_flag, name = "registros") %>%
                        arrange(desc(registros))))
  cat("\nOs 10 maiores valores reprovados:\n")
  print(as.data.frame(reprovados %>%
    select(gauge_code, city, state, datetime, rain_mm, timestep, teto_mm, qc_flag) %>%
    head(10)))

  write_xlsx(reprovados %>% mutate(datetime = format(datetime, "%Y-%m-%d %H:%M:%S")),
             file.path(saida, "qc_basico_registros_excluidos.xlsx"))

  por_estacao <- reprovados %>%
    count(gauge_code, city, state, qc_flag, name = "registros") %>%
    arrange(desc(registros))
  write_csv(por_estacao, file.path(saida, "qc_basico_por_estacao.csv"))
} else {
  cat("\nNenhum registro reprovado no QC Básico.\n")
}

## 4.4 Corte dos removidos da série ----------------------------------------------------
etapa <- file.path(saida, "_qcbasico")
unlink(etapa, recursive = TRUE)

linhas_marcadas %>%
  filter(!U1, !U2, !U3, !U4, !U5) %>%
  select(gauge_code, datetime, rain_mm, state, ano) %>%
  write_dataset(etapa, format = "parquet",
                basename_template = paste0("series_", PREFIXO, "_qcbasico-{i}.parquet"))

arq_gerados <- list.files(etapa, pattern = "[.]parquet$", full.names = TRUE)
arq_qc      <- file.path(saida, paste0("series_", PREFIXO, "_qcbasico.parquet"))

if (length(arq_gerados) == 1) {
  if (file.exists(arq_qc)) file.remove(arq_qc)
  file.rename(arq_gerados, arq_qc)
  unlink(etapa, recursive = TRUE)
} else {
  stop("O QC Básico gerou ", length(arq_gerados), " arquivos; esperado 1. Veja ", etapa)
}

serie_qc <- open_dataset(arq_qc)

n_removido <- ds$num_rows - serie_qc$num_rows
esperado   <- nrow(reprovados)
if (n_removido != esperado)
  stop("Seção 4: removidos ", n_removido, " registros, esperados ", esperado, ".")

## 4.5 Resumo ----------------------------------------------------------------
flags <- c("VALOR NEGATIVO", "VALOR NÃO FINITO", "VALOR ACIMA DO TETO",
           "VALOR EM ANO DIFERENTE", "SEM CADASTRO DE ESTAÇÃO")

resumo_basico <- tibble(
  categoria = c("VÁLIDO", flags),
  n = c(serie_qc$num_rows,
        vapply(flags, function(f) sum(reprovados$qc_flag == f), integer(1)))) %>%
  mutate(percentual = round(100 * n / ds$num_rows, 6))

write_xlsx(resumo_basico, file.path(saida, "qc_basico_resumo.xlsx"))

cat(sprintf("\nQC Básico: %s de %s registros descartados (%.6f%%)\n",
            format(n_removido, big.mark = ".", decimal.mark = ","),
            format(ds$num_rows, big.mark = ".", decimal.mark = ","),
            100 * n_removido / ds$num_rows))
cat(sprintf("Registros restantes: %s\n",
            format(serie_qc$num_rows, big.mark = ".", decimal.mark = ",")))

##### 5. Verificações por estação ----------------------------------------------------------------------------------------------------------
## S1 - mesmo valor não nulo repetido em 12 registros == CHUVA CONSTANTE
## S2 - mais de 99% de zeros no ano-estação == CHUVA SEM VARIAÇÃO
## S3 - pico isolado entre zeros, acima do limiar da estação == PICO ISOLADO
## S4 - acima da envoltória, fator 1 == ACIMA DA ENVOLTÓRIA
## S5 - acima do limiar da estação == OUTLIER ESTATÍSTICO
## S6 - período seco longo == PERÍODO SECO
S1_MIN_RUN     <- 12     # registros
S2_FRAC_ZERO   <- 0.99   # acima disso o ano-estação é marcado
S3_K           <- 3      # multiplicador do IQR
ESTACOES_BLOCO <- 100    # estações lidas por vez
S6_MIN_ANOS   <- 5     # anos mínimos para calcular o limiar de dias secos
S6_PADRAO_DIAS <- 270  # dias; limiar fixo de quem tem histórico curto
S6_K          <- 3     # multiplicador do IQR dos máximos anuais

serie <- open_dataset(arq_qc)
cat("Registros na série pós-QC Básico:",
    format(serie$num_rows, big.mark = ".", decimal.mark = ","), "\n")

## 5.1 S2 - ano-estação com excesso de zeros -----------------------------------------
por_ano <- serie %>%
  group_by(gauge_code, ano) %>%
  summarise(registros = n(), zeros = sum(rain_mm == 0)) %>%
  collect()

por_ano$fracao_zeros <- por_ano$zeros / por_ano$registros
s2 <- por_ano[por_ano$fracao_zeros > S2_FRAC_ZERO, ]

cat(sprintf("\nS2: %s de %s anos-estação com mais de %.0f%% de zeros, em %s estações\n",
            format(nrow(s2), big.mark = ".", decimal.mark = ","),
            format(nrow(por_ano), big.mark = ".", decimal.mark = ","),
            100 * S2_FRAC_ZERO,
            format(length(unique(s2$gauge_code)), big.mark = ".", decimal.mark = ",")))

## 5.2 S1 e S3, estação por estação -------------------------------------------
codigos <- sort(collect(distinct(select(serie, gauge_code)))$gauge_code)
blocos  <- split(codigos, ceiling(seq_along(codigos) / ESTACOES_BLOCO))

cat(sprintf("\n%s estações em %d blocos de até %d\n",
            format(length(codigos), big.mark = ".", decimal.mark = ","),
            length(blocos), ESTACOES_BLOCO))

lista_limiar <- list()
lista_s1     <- list()
lista_s3     <- list()

for (b in seq_along(blocos)) {
  cods <- blocos[[b]]
  x <- serie %>%
    filter(gauge_code %in% cods) %>%
    select(gauge_code, datetime, rain_mm) %>%
    collect()
  setDT(x)
  setorder(x, gauge_code, datetime)

  ## o limiar de cada estação, sobre os valores não nulos
  limiar <- x[rain_mm > 0,
              .(q3  = quantile(rain_mm, 0.75, names = FALSE),
                iqr = IQR(rain_mm),
                n_nao_nulos = .N),
              by = gauge_code]
  limiar[, limiar_mm := q3 + S3_K * iqr]
  lista_limiar[[b]] <- limiar

  ## S1: mesmo valor não nulo em registros seguidos
  x[, troca := is.na(shift(rain_mm)) | rain_mm != shift(rain_mm), by = gauge_code]
  x[, corrida := cumsum(troca), by = gauge_code]
  x[, tamanho := .N, by = .(gauge_code, corrida)]

  s1_bloco <- x[rain_mm > 0 & tamanho >= S1_MIN_RUN,
                .(gauge_code, datetime, rain_mm, corrida, tamanho)]
  lista_s1[[b]] <- s1_bloco

  x[, c("troca", "corrida", "tamanho") := NULL]

  ##S3: pico isolado entre zeros
  x[, anterior  := shift(rain_mm, 1L, type = "lag"),  by = gauge_code]
  x[, posterior := shift(rain_mm, 1L, type = "lead"), by = gauge_code]
  x <- limiar[, .(gauge_code, limiar_mm)][x, on = "gauge_code"]

  s3_bloco <- x[!is.na(anterior) & !is.na(posterior) &
                  anterior == 0 & posterior == 0 &
                  !is.na(limiar_mm) & rain_mm > limiar_mm,
                .(gauge_code, datetime, rain_mm, limiar_mm)]
  lista_s3[[b]] <- s3_bloco

  cat(sprintf("  bloco %2d: %s registros | S1 %s | S3 %s\n", b,
              format(nrow(x), big.mark = ".", decimal.mark = ","),
              format(nrow(s1_bloco), big.mark = ".", decimal.mark = ","),
              format(nrow(s3_bloco), big.mark = ".", decimal.mark = ",")))

  rm(x, limiar, s1_bloco, s3_bloco); invisible(gc(verbose = FALSE))
}

limiar_estacao <- rbindlist(lista_limiar)
s1 <- rbindlist(lista_s1)
s3 <- rbindlist(lista_s3)

## 5.3 S4 - valor acima da envoltória ---------------------------------------
s4 <- serie %>%
  left_join(limites[, c("gauge_code", "ano", "teto_mm")],
            by = c("gauge_code", "ano")) %>%
  filter(!is.na(teto_mm), rain_mm > teto_mm / FATOR_LIMITE) %>%
  select(gauge_code, datetime, ano, rain_mm, teto_mm) %>%
  collect()

s4$envoltoria_mm <- round(s4$teto_mm / FATOR_LIMITE, 1)

cat(sprintf("\nS4: %s registros acima da envoltória, em %s estações\n",
            format(nrow(s4), big.mark = ".", decimal.mark = ","),
            format(n_distinct(s4$gauge_code), big.mark = ".", decimal.mark = ",")))

## 5.4 S5 - outlier acima do limiar da estação -------------------------------
s5 <- serie %>%
  left_join(limiar_estacao[, .(gauge_code, limiar_mm)], by = "gauge_code") %>%
  filter(!is.na(limiar_mm), rain_mm > limiar_mm) %>%
  select(gauge_code, datetime, ano, rain_mm, limiar_mm) %>%
  collect()

cat(sprintf("\nS5: %s registros acima do limiar, em %s estações\n",
            format(nrow(s5), big.mark = ".", decimal.mark = ","),
            format(n_distinct(s5$gauge_code), big.mark = ".", decimal.mark = ",")))

## 5.5 S6 - período seco longo -----------------------------------------------
dia <- serie %>%
  mutate(data = as.Date(datetime)) %>%
  group_by(gauge_code, data) %>%
  summarise(total = sum(rain_mm)) %>%
  collect()
setDT(dia)
setorder(dia, gauge_code, data)

dia[, seco := total == 0]

dia[, novo := is.na(shift(seco)) | !seco | !shift(seco) |
              as.integer(data - shift(data)) != 1L, by = gauge_code]
dia[, periodo := cumsum(novo), by = gauge_code]

periodos <- dia[seco == TRUE,
                .(ini = min(data), fim = max(data), dias = .N),
                by = .(gauge_code, periodo)]

anos_do_periodo <- periodos[, .(ano = seq(year(ini), year(fim))),
                            by = .(gauge_code, periodo, dias)]
max_por_ano <- anos_do_periodo[, .(max_seco = max(dias)), by = .(gauge_code, ano)]

## anos sem nenhum período seco entram com zero
anos_estacao <- unique(dia[, .(gauge_code, ano = year(data))])
max_por_ano  <- merge(anos_estacao, max_por_ano,
                      by = c("gauge_code", "ano"), all.x = TRUE)
max_por_ano[is.na(max_seco), max_seco := 0L]

limiar_seco <- max_por_ano[, .(n_anos = .N,
                               q3  = quantile(max_seco, 0.75, names = FALSE),
                               iqr = IQR(max_seco)), by = gauge_code]
limiar_seco[, limiar_dias := q3 + S6_K * iqr]
limiar_seco[n_anos < S6_MIN_ANOS, limiar_dias := S6_PADRAO_DIAS]

periodos <- limiar_seco[, .(gauge_code, limiar_dias)][periodos, on = "gauge_code"]
periodos_longos <- periodos[dias > limiar_dias]

s6_dias <- dia[seco == TRUE][
  periodos_longos[, .(gauge_code, periodo)], on = .(gauge_code, periodo),
  nomatch = 0][, .(gauge_code, data)]

cat(sprintf("\nS6: %s períodos secos acima do limiar, em %s estações\n",
            format(nrow(periodos_longos), big.mark = ".", decimal.mark = ","),
            format(uniqueN(periodos_longos$gauge_code), big.mark = ".", decimal.mark = ",")))
cat(sprintf("     %s dias a marcar | estações no limiar fixo de %d dias: %s\n",
            format(nrow(s6_dias), big.mark = ".", decimal.mark = ","),
            S6_PADRAO_DIAS,
            format(sum(limiar_seco$n_anos < S6_MIN_ANOS), big.mark = ".", decimal.mark = ",")))

## 5.6 Resumo e gravação ------------------------------------------------------
cat(sprintf("\nS1: %s registros marcados em %s sequências, em %s estações\n",
            format(nrow(s1), big.mark = ".", decimal.mark = ","),
            format(nrow(unique(s1[, .(gauge_code, corrida)])), big.mark = ".", decimal.mark = ","),
            format(uniqueN(s1$gauge_code), big.mark = ".", decimal.mark = ",")))
cat(sprintf("S3: %s registros marcados em %s estações\n",
            format(nrow(s3), big.mark = ".", decimal.mark = ","),
            format(uniqueN(s3$gauge_code), big.mark = ".", decimal.mark = ",")))

if (nrow(s1) > 0) {
  cat("\nTamanho das sequências do S1:\n")
  print(as.data.frame(unique(s1[, .(gauge_code, corrida, tamanho)])[
    , .N, by = .(faixa = cut(tamanho, c(11, 12, 24, 48, 144, Inf),
      labels = c("12", "13-24", "25-48", "49-144", "acima de 144")))][order(faixa)]))
}

write_csv(as_tibble(limiar_estacao),
          file.path(saida, paste0("limiar_estacao_", PREFIXO, ".csv")))
write_csv(as_tibble(s2), file.path(saida, paste0("flag_s2_", PREFIXO, ".csv")))
write_csv(as_tibble(s1), file.path(saida, paste0("flag_s1_", PREFIXO, ".csv")))
write_csv(as_tibble(s3), file.path(saida, paste0("flag_s3_", PREFIXO, ".csv")))
write_csv(as_tibble(s4), file.path(saida, paste0("flag_s4_", PREFIXO, ".csv")))
write_csv(as_tibble(s5), file.path(saida, paste0("flag_s5_", PREFIXO, ".csv")))
write_csv(as_tibble(s6_dias), file.path(saida, paste0("flag_s6_dias_", PREFIXO, ".csv")))
write_csv(as_tibble(periodos_longos), file.path(saida, paste0("s6_periodos_", PREFIXO, ".csv")))
write_csv(as_tibble(limiar_seco), file.path(saida, paste0("s6_limiar_", PREFIXO, ".csv")))
write_csv(as_tibble(max_por_ano), file.path(saida, paste0("s6_maximo_anual_", PREFIXO, ".csv")))

##### 6. Consolidação dos rótulos ----------------------------------------------------------------------------------------------------------
final_s1  <- unique(as_tibble(s1[,  .(gauge_code, datetime)]))
final_s2  <- unique(as_tibble(as.data.table(s2)[, .(gauge_code, ano = as.integer(ano))]))
final_s3  <- unique(as_tibble(s3[,  .(gauge_code, datetime)]))
final_s4 <- unique(s4[, c("gauge_code", "datetime")])
final_s5 <- unique(s5[, c("gauge_code", "datetime")])
final_s6 <- unique(as_tibble(s6_dias))

final_s1$f_s1   <- TRUE
final_s2$f_s2   <- TRUE
final_s3$f_s3   <- TRUE
final_s4$f_s4 <- TRUE
final_s5$f_s5 <- TRUE
final_s6$f_s6 <- TRUE

etapa_cons <- file.path(saida, "_consolidada")
unlink(etapa_cons, recursive = TRUE)

serie %>%
  mutate(data = as.Date(datetime)) %>%
  left_join(final_s1,  by = c("gauge_code", "datetime")) %>%
  left_join(final_s2,  by = c("gauge_code", "ano")) %>%
  left_join(final_s3,  by = c("gauge_code", "datetime")) %>%
  left_join(final_s4, by = c("gauge_code", "datetime")) %>%
  left_join(final_s5, by = c("gauge_code", "datetime")) %>%
  left_join(final_s6, by = c("gauge_code", "data")) %>%
  mutate(qc_flag = case_when(
           !is.na(f_s1)  ~ "CHUVA CONSTANTE",
           !is.na(f_s2)  ~ "CHUVA SEM VARIAÇÃO",
           !is.na(f_s3)  ~ "PICO ISOLADO",
           !is.na(f_s4) ~ "ACIMA DA ENVOLTÓRIA",
           !is.na(f_s5) ~ "OUTLIER ESTATÍSTICO",
           !is.na(f_s6) ~ "PERÍODO SECO",
           TRUE          ~ "PASS")) %>%
  mutate(flagged = qc_flag != "PASS") %>%
  select(gauge_code, datetime, rain_mm, state, ano, qc_flag, flagged) %>%
  write_dataset(etapa_cons, format = "parquet",
                basename_template = paste0("series_", PREFIXO, "_qc-{i}.parquet"))

arq_cons_gerados <- list.files(etapa_cons, pattern = "[.]parquet$", full.names = TRUE)
arq_cons         <- file.path(saida, paste0("series_", PREFIXO, "_qc.parquet"))

if (length(arq_cons_gerados) == 1) {
  if (file.exists(arq_cons)) file.remove(arq_cons)
  file.rename(arq_cons_gerados, arq_cons)
  unlink(etapa_cons, recursive = TRUE)
} else {
  stop("A consolidação gerou ", length(arq_cons_gerados), " arquivos; esperado 1. Veja ", etapa_cons)
}

serie_qc_final <- open_dataset(arq_cons)

if (serie_qc_final$num_rows != serie$num_rows)
  stop("Seção 6: a consolidação mudou o número de registros, de ",
       serie$num_rows, " para ", serie_qc_final$num_rows,
       ". Nenhuma junção deveria duplicar linha.")

resumo_flags <- serie_qc_final %>%
  group_by(qc_flag) %>%
  summarise(registros = n(), estacoes = n_distinct(gauge_code)) %>%
  collect() %>%
  mutate(percentual = round(100 * registros / serie_qc_final$num_rows, 4)) %>%
  arrange(desc(registros))

cat("\nRótulo final, por registro:\n")
print(as.data.frame(resumo_flags))

write_csv(resumo_flags, file.path(saida, paste0("resumo_flags_", PREFIXO, ".csv")))

##### 7. Disponibilidade ----------------------------------------------------------------------------------------------------------
MIN_DISP <- 80        # % mínima de tempo coberto no ano
MIN_ANOS <- 8         # anos válidos mínimos por estação

## A duração de cada registro é o timestep declarado da própria estação-ano, e
## não uma constante de rede: esta rede transmite em intervalo fixo, que varia
## entre estações e ao longo dos anos. O limites da Seção 4.1 já traz esse
## timestep em minutos, lido do cadastro.
passo_ano <- as.data.table(limites)[, .(gauge_code, ano, passo_s = timestep * 60)]

sem_passo <- passo_ano[is.na(passo_s) | passo_s <= 0]
cat("Estação-anos sem timestep utilizável, fora da medição:", nrow(sem_passo), "\n")

passo_ano <- passo_ano[!is.na(passo_s) & passo_s > 0]

cat("\nTimestep declarado, em pares estação-ano:\n")
print(as.data.frame(passo_ano[, .(pares = .N),
                              by = .(minutos = passo_s / 60)][order(minutos)]))

## O sort descarta NA em silêncio, e o ano afetado nunca seria medido: confere
## antes de ordenar.
anos_disp <- collect(distinct(select(serie_qc_final, ano)))$ano
if (any(is.na(anos_disp))) {
  stop("Seção 7: há registros sem ano.")
}
anos_disp <- sort(anos_disp)
res_disp  <- vector("list", length(anos_disp))

for (k in seq_along(anos_disp)) {

  a <- anos_disp[k]

  x <- serie_qc_final %>%
    filter(ano == a, !is.na(datetime), !is.na(rain_mm), rain_mm >= 0) %>%
    select(gauge_code, datetime, rain_mm) %>%
    collect()
  setDT(x)
  if (nrow(x) == 0L) next

  ini_ano <- as.numeric(as.POSIXct(sprintf("%d-01-01", a),     tz = "UTC"))
  fim_ano <- as.numeric(as.POSIXct(sprintf("%d-01-01", a + 1), tz = "UTC"))

  x[, t   := as.numeric(datetime)]
  x[passo_ano[ano == a, .(gauge_code, passo_s)], on = "gauge_code", dur := i.passo_s]
  x <- x[!is.na(dur)]
  x[, `:=`(ini = pmax(t - dur, ini_ano), fim = pmin(t, fim_ano))]
  x <- x[ini < fim]
  x[, c("datetime", "rain_mm", "t", "dur") := NULL]

  setorder(x, gauge_code, ini)
  x[, fim_corrido := cummax(fim), by = gauge_code]
  x[, novo := ini > shift(fim_corrido, fill = -Inf), by = gauge_code]
  x[, grp  := cumsum(novo), by = gauge_code]

  res_disp[[k]] <- x[, .(ini = min(ini), fim = max(fim)), by = .(gauge_code, grp)][
    , .(coberto_s = sum(fim - ini),
        obs_ini   = min(ini),
        obs_fim   = max(fim)), by = gauge_code][, `:=`(ano = a, ano_s = fim_ano - ini_ano)]

  cat(sprintf("  %d: %s estações\n", a,
              format(nrow(res_disp[[k]]), big.mark = ".", decimal.mark = ",")))

  rm(x); invisible(gc(verbose = FALSE))
}

disp <- rbindlist(res_disp)
disp[, janela_s := obs_fim - obs_ini]   # período observado, só diagnóstico
disp[, disp     := round(100 * coberto_s / ano_s, 1)]
disp[, passa    := !is.na(disp) & disp >= MIN_DISP]
disp[, `:=`(horas_cobertas = round(coberto_s / 3600, 1),
            horas_ano      = round(ano_s     / 3600, 1),
            horas_janela   = round(janela_s  / 3600, 1))]

anos_bons <- disp[passa == TRUE]
contagem  <- anos_bons[, .(anos_bons = .N), by = gauge_code]
estacoes  <- contagem[anos_bons >= MIN_ANOS]
anos_fica <- anos_bons[gauge_code %in% estacoes$gauge_code,
                       .(gauge_code, ano = as.integer(ano))]

write_csv(as_tibble(disp[order(gauge_code, ano)][
            , .(gauge_code, ano, disp, passa,
                horas_cobertas, horas_ano, horas_janela)]),
          file.path(saida, paste0("disponibilidade_", PREFIXO, ".csv")))

cat("\nAnos bons por estação:\n")
print(as.data.frame(contagem[, .N, by = anos_bons][order(anos_bons)]))

etapa_fim <- file.path(saida, "_final")
unlink(etapa_fim, recursive = TRUE)

serie_qc_final %>%
  inner_join(anos_fica, by = c("gauge_code", "ano")) %>%
  write_dataset(etapa_fim, format = "parquet",
                basename_template = paste0("series_", PREFIXO, "_final-{i}.parquet"))

arq_fim_gerados <- list.files(etapa_fim, pattern = "[.]parquet$", full.names = TRUE)
arq_fim         <- file.path(saida, paste0("series_", PREFIXO, "_final.parquet"))

if (length(arq_fim_gerados) == 1) {
  if (file.exists(arq_fim)) file.remove(arq_fim)
  file.rename(arq_fim_gerados, arq_fim)
  unlink(etapa_fim, recursive = TRUE)
} else {
  stop("A seleção final gerou ", length(arq_fim_gerados), " arquivos. Veja ", etapa_fim)
}

cat(sprintf("Estações: %d -> %d (>= %d anos com >= %d%% de cobertura)\n",
            length(codigos), nrow(estacoes), MIN_ANOS, MIN_DISP))
