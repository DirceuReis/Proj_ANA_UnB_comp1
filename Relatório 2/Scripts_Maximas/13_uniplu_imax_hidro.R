# =============================================================================
# PASSO 13 — IMAX DIÁRIO POR ANO HIDROLÓGICO + QC20% — BRASIL
# =============================================================================
#
# Fluxo:
# stations_analyzed/<codigo>_analyzed.rds
#   -> regulariza grade diária
#   -> imax.wateryear()
#      * ano hidrológico = mes_hibrido_rbar
#      * QC <=20% no ano hidrológico completo
#      * falhas estruturais das bordas entram no QC
#      * janelas 100% NA, descontínuas ou cruzando ano hidrológico são inválidas
#
# Durações: 1, 2, 3, 4, 5, 6, 7 e 10 dias
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(arrow)
library(lubridate)

# =============================================================================
# 1. CONFIGURAÇÃO
# =============================================================================

DIR_FLUXO <- "Relatório 2/Scripts_Maximas"
DIR_ANO_HIDRO <- "Relatório 2/Scripts_Ano_Hidrologico"
DIR_DF_HIDRO <- file.path(DIR_ANO_HIDRO,"resultados","BR","dataframes")
DIR_OUT <- file.path(DIR_FLUXO,"resultados","BR")
DIR_DF <- file.path(DIR_OUT,"dataframes")
DIR_CACHE <- file.path(DIR_OUT,"uniplu_cache","imax_by_uf_hidro_qc20")
DIR_ANALYZED <- file.path(DIR_ANO_HIDRO,"resultados","BR","stations_analyzed")
DIR_FUNS <- file.path(DIR_FLUXO,"funs")

NA_ACCEPT_HYDRO <- 0.20
APPLY_QC_HYDRO <- TRUE
DURATIONS_DAY <- c(1L,2L,3L,4L,5L,6L,7L,10L)
DURATIONS_HR <- DURATIONS_DAY*24
meses_pt <- c("jan","fev","mar","abr","mai","jun","jul","ago","set","out","nov","dez")

dir.create(DIR_DF,recursive=TRUE,showWarnings=FALSE)
dir.create(DIR_CACHE,recursive=TRUE,showWarnings=FALSE)

source(file.path(DIR_FUNS,"fun_imax_wateryear.R"))
if(!exists("imax.wateryear")) stop("Função 'imax.wateryear()' não encontrada.")

achar_rds <- function(nome){
  p <- file.path(DIR_DF_HIDRO,nome)
  if(file.exists(p)) p else NULL
}

if(!dir.exists(DIR_ANALYZED)) stop("Pasta stations_analyzed não encontrada: ",DIR_ANALYZED)

# =============================================================================
# 2. MÊS INICIAL DO ANO HIDROLÓGICO
# =============================================================================

F_HIB <- achar_rds("df_hidro_rbar_hibrido.rds")
if(is.null(F_HIB)) stop("Falta df_hidro_rbar_hibrido.rds — rode o passo 6.")
message("Híbrido: ",F_HIB)

hib <- readRDS(F_HIB) %>%
  filter(
    is.finite(mes_hibrido_rbar),mes_hibrido_rbar>=1,mes_hibrido_rbar<=12,
    is.finite(lat),is.finite(long)
  ) %>%
  mutate(
    codigo=trimws(as.character(codigo)),
    mes_inicio=as.integer(mes_hibrido_rbar),
    mes_inicio_nome=meses_pt[mes_inicio]
  )

ufs <- sort(unique(na.omit(hib$estado)))

message("Postos: ",nrow(hib),
        " | UFs: ",paste(ufs,collapse=", "),
        " | durações: ",paste(DURATIONS_DAY,collapse=", ")," dias")

# =============================================================================
# 3. PROCESSAMENTO POR UF
# =============================================================================

all_imax <- list()

for(uf in ufs){
  
  f_ckpt <- file.path(DIR_CACHE,paste0(uf,"_imax.rds"))
  
  if(file.exists(f_ckpt)){
    message("\n[",uf,"] checkpoint encontrado — carregando")
    all_imax[[uf]] <- readRDS(f_ckpt)
    next
  }
  
  hib_uf <- hib %>% filter(estado==uf)
  if(!nrow(hib_uf)) next
  
  message("\n==================================================")
  message("[",uf,"] postos=",nrow(hib_uf))
  message("==================================================")
  
  data_ls <- vector("list",nrow(hib_uf))
  names(data_ls) <- hib_uf$codigo
  
  # -------------------------------------------------------------------------
  # Ler e regularizar séries diárias
  # -------------------------------------------------------------------------
  
  for(i in seq_len(nrow(hib_uf))){
    
    cod <- hib_uf$codigo[i]
    f_sta <- file.path(DIR_ANALYZED,paste0(cod,"_analyzed.rds"))
    if(!file.exists(f_sta)) next
    
    sta <- tryCatch(readRDS(f_sta),error=function(e){
      message("[",uf,"] erro lendo ",cod,": ",conditionMessage(e)); NULL
    })
    if(is.null(sta) || !"dt" %in% names(sta) || !"p" %in% names(sta)) next
    
    x <- data.frame(
      date=as.Date(sta$dt),
      rain_mm=as.numeric(sta$p)
    ) %>%
      filter(!is.na(date)) %>%
      arrange(date) %>%
      distinct(date,.keep_all=TRUE)
    
    if(nrow(x)<2) next
    
    # Grade diária regular entre primeira e última observação
    grade <- data.frame(date=seq(min(x$date),max(x$date),by="day")) %>%
      left_join(x,by="date")
    
    data_ls[[cod]] <- data.frame(
      datetime=as.POSIXct(grade$date,tz="UTC"),
      rain_mm=grade$rain_mm
    )
  }
  
  data_ls <- data_ls[!vapply(data_ls,is.null,logical(1))]
  
  if(!length(data_ls)){
    message("[",uf,"] nenhuma série válida.")
    next
  }
  
  message("[",uf,"] séries OK: ",length(data_ls),"/",nrow(hib_uf))
  
  # -------------------------------------------------------------------------
  # Mês inicial específico de cada estação
  # -------------------------------------------------------------------------
  
  start_map <- setNames(as.integer(hib_uf$mes_inicio),hib_uf$codigo)
  start_map <- start_map[names(data_ls)]
  keep <- is.finite(start_map)
  data_ls <- data_ls[keep]
  start_map <- start_map[keep]
  names(start_map) <- names(data_ls)
  
  if(!length(data_ls)) next
  
  # -------------------------------------------------------------------------
  # IMAX + QC20 no ano hidrológico
  # -------------------------------------------------------------------------
  
  imax_uf <- tryCatch(
    imax.wateryear(
      data=data_ls,
      durations=DURATIONS_HR,
      start_month=start_map,
      which.mon=1:12,
      names=c("datetime","rain_mm"),
      na_accept=NA_ACCEPT_HYDRO,
      apply_qc=APPLY_QC_HYDRO
    ),
    error=function(e){
      message("[",uf,"] ERRO imax: ",conditionMessage(e))
      NULL
    }
  )
  
  rm(data_ls); invisible(gc())
  
  if(is.null(imax_uf) || !nrow(imax_uf)){
    message("[",uf,"] nenhum IMAX gerado.")
    next
  }
  
  # -------------------------------------------------------------------------
  # Metadados
  # -------------------------------------------------------------------------
  
  imax_uf <- imax_uf %>%
    mutate(
      estado=uf,
      d_days=d/24,
      time_step_min=1440L,
      qc_na_accept=NA_ACCEPT_HYDRO,
      qc_applied=APPLY_QC_HYDRO
    ) %>%
    left_join(
      hib_uf %>%
        select(codigo,nome,lat,long,network,year.size,rbar,mes_inicio,
               mes_inicio_nome,fonte_hibrido_rbar,classe_rbar),
      by=c("gauge_code"="codigo")
    )
  
  qc_uf <- imax_uf %>%
    distinct(gauge_code,wateryear,start_month,na_prct_yr,
             complete_wateryear,qc_year)
  
  message("[",uf,"] linhas=",nrow(imax_uf),
          " | postos=",n_distinct(imax_uf$gauge_code),
          " | posto-anos aceitos=",nrow(qc_uf))
  
  saveRDS(imax_uf,f_ckpt)
  all_imax[[uf]] <- imax_uf
}

# =============================================================================
# 4. CONSOLIDAÇÃO BRASIL
# =============================================================================

df_imax <- bind_rows(all_imax)
if(!nrow(df_imax)) stop("Nenhum IMAX diário foi gerado.")

n_div <- df_imax %>%
  filter(is.finite(start_month),is.finite(mes_inicio),start_month!=mes_inicio) %>%
  nrow()
if(n_div>0) warning(n_div," linhas com start_month != mes_inicio.")

# =============================================================================
# 5. SALVAR BASE IMAX
# =============================================================================

F_OUT_RDS <- file.path(DIR_DF,"df_imax_uniplu_hidro.rds")
F_OUT_PQT <- file.path(DIR_DF,"df_imax_uniplu_hidro.parquet")
F_OUT_CSV <- file.path(DIR_DF,"df_imax_uniplu_hidro.csv")

saveRDS(df_imax,F_OUT_RDS)
write_parquet(df_imax,F_OUT_PQT)

df_imax %>%
  mutate(across(any_of(c("date","wy_start","wy_end")),as.character)) %>%
  write_excel_csv(F_OUT_CSV)

# =============================================================================
# 6. RESUMO DAS DURAÇÕES
# =============================================================================

tab_d <- df_imax %>%
  filter(is.finite(imax)) %>%
  group_by(d_days,d) %>%
  summarise(
    n=n(),
    n_postos=n_distinct(gauge_code),
    med_mm_h=round(median(imax),4),
    p95_mm_h=round(quantile(imax,0.95),4),
    .groups="drop"
  )

print(tab_d,n=Inf)

write_excel_csv(
  tab_d,
  file.path(DIR_DF,"df_imax_uniplu_hidro_resumo.csv")
)

# =============================================================================
# 7. RESUMO DO QC
# =============================================================================

qc_final <- df_imax %>%
  distinct(gauge_code,estado,wateryear,start_month,na_prct_yr,
           complete_wateryear,qc_year,n_expected_yr,n_present_yr,
           n_na_present_yr,n_missing_structural_yr)

message("\n==================================================")
message("PASSO 13 — DIÁRIO BRASIL CONCLUÍDO")
message("==================================================")
message("Linhas IMAX: ",format(nrow(df_imax),big.mark="."))
message("Postos: ",n_distinct(df_imax$gauge_code))
message("UFs: ",n_distinct(df_imax$estado))
message("Posto-anos aceitos QC20: ",nrow(qc_final))
message("Anos aceitos incompletos: ",
        sum(!qc_final$complete_wateryear,na.rm=TRUE))
message("NA anual: min=",
        round(100*min(qc_final$na_prct_yr,na.rm=TRUE),2),
        "% | mediana=",
        round(100*median(qc_final$na_prct_yr,na.rm=TRUE),2),
        "% | máx=",
        round(100*max(qc_final$na_prct_yr,na.rm=TRUE),2),"%")

# =============================================================================
# 8. SELEÇÃO PARA O MODELO
# 8 durações (1,2,3,4,5,6,7,10 dias) e > 30 anos em TODAS
# =============================================================================

sel_diario_qc20 <- df_imax %>%
  filter(d_days %in% DURATIONS_DAY, is.finite(imax)) %>%
  distinct(gauge_code, estado, d_days, wateryear) %>%
  count(gauge_code, estado, d_days, name="n_anos_validos") %>%
  group_by(gauge_code, estado) %>%
  summarise(
    n_duracoes=n_distinct(d_days),
    n_duracoes_30anos=sum(n_anos_validos > 30),
    min_anos=min(n_anos_validos),
    med_anos=median(n_anos_validos),
    max_anos=max(n_anos_validos),
    .groups="drop"
  ) %>%
  mutate(
    passa=n_duracoes==length(DURATIONS_DAY) &
      n_duracoes_30anos==length(DURATIONS_DAY)
  ) %>%
  arrange(desc(passa), estado, desc(min_anos))

estacoes_aprovadas <- sel_diario_qc20 %>% filter(passa)

message("\n==================================================")
message("SELEÇÃO 8 DURAÇÕES × > 30 ANOS")
message("==================================================")
message("Estações aprovadas: ", nrow(estacoes_aprovadas))

print(estacoes_aprovadas, n=Inf, width=Inf)

write_excel_csv(
  sel_diario_qc20,
  file.path(DIR_DF, "selecao_estacoes_diario_hidro_qc20_mais30anos.csv")
)

# =============================================================================
# 9. RESUMO POR UF
# =============================================================================

resumo_uf <- sel_diario_qc20 %>%
  group_by(estado) %>%
  summarise(
    n_avaliadas=n(),
    n_aprovadas=sum(passa),
    .groups="drop"
  ) %>%
  arrange(desc(n_aprovadas))

print(resumo_uf, n=Inf)

write_excel_csv(
  resumo_uf,
  file.path(DIR_DF, "resumo_selecao_diario_por_uf_mais30anos.csv")
)

message("\nArquivos principais:")
message(F_OUT_RDS)
message(F_OUT_PQT)
message(F_OUT_CSV)
message("Cache por UF: ",DIR_CACHE)
