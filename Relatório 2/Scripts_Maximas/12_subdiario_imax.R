# =============================================================================
# PASSO 12 — IMAX SUBDIÁRIO POR ANO HIDROLÓGICO + QC20% — BRASIL
# =============================================================================
#
# Fluxo:
# dados brutos
#   -> fun_filter_set(filter = FALSE)
#      * CEMADEN: gaps 10 < delta <= 60 min são preenchidos com zero
#   -> time_step: segundos -> minutos
#   -> fun_group_ts()
#   -> imax.wateryear()
#      * ano hidrológico específico de cada estação
#      * QC <= 20% no ano hidrológico completo
#      * falhas estruturais das bordas entram no QC
#      * janelas 100% NA, descontínuas ou cruzando ano hidrológico são inválidas
#
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
DIR_DF_HIDRO <- file.path(DIR_ANO_HIDRO, "resultados", "BR", "dataframes")
DIR_OUT <- file.path(DIR_FLUXO, "resultados", "BR")
DIR_DF <- file.path(DIR_OUT, "dataframes")
DIR_CACHE <- file.path(DIR_OUT, "subdiario_cache", "imax_by_uf_hidro_qc20")
DIR_SUBDIARIO <- file.path("base", "fonte", "consolidado", "subdiario_br")
DIR_FUNS <- file.path(DIR_FLUXO, "funs")

NA_ACCEPT_HYDRO <- 0.20
APPLY_QC_HYDRO <- TRUE
DURATIONS_HR <- c(c(10,15,20,30,40,45,50)/60, 1:24)
meses_pt <- c("jan","fev","mar","abr","mai","jun","jul","ago","set","out","nov","dez")

dir.create(DIR_DF, recursive=TRUE, showWarnings=FALSE)
dir.create(DIR_CACHE, recursive=TRUE, showWarnings=FALSE)

source(file.path(DIR_FUNS, "fun_filter_set.R"))
source(file.path(DIR_FUNS, "fun_group_ts.R"))
source(file.path(DIR_FUNS, "fun_imax_wateryear.R"))
if(!exists("imax.wateryear")) stop("Função 'imax.wateryear()' não encontrada.")

achar_rds <- function(nome){
  p <- file.path(DIR_DF_HIDRO, nome)
  if(file.exists(p)) p else NULL
}

# =============================================================================
# 2. MÊS INICIAL + INVENTÁRIO
# =============================================================================

F_INI <- achar_rds("df_subdiario_ano_hidro_uniplu.rds")
if(is.null(F_INI)) stop("Falta df_subdiario_ano_hidro_uniplu.rds — rode o passo 9.")
ini <- readRDS(F_INI) %>% mutate(gauge_code=trimws(as.character(gauge_code)))

F_INV <- achar_rds("df_subdiario_inventario.rds")
if(is.null(F_INV)) F_INV <- achar_rds("df_subdiario_mes_oposto_nearest.rds")
if(is.null(F_INV)) stop("Falta inventário subdiário.")
inv <- readRDS(F_INV) %>% mutate(gauge_code=trimws(as.character(gauge_code)))

inv_cols <- intersect(
  c("gauge_code","nome","estado","lat","long","time_step","elevation",
    "network","responsible","n_years_zip","year_first","year_last"),
  names(inv)
)
inv <- inv %>% select(all_of(inv_cols))
if(!"responsible" %in% names(inv)) inv$responsible <- NA_character_
if(anyDuplicated(inv$gauge_code)){
  warning("Gauge duplicado no inventário; mantendo primeira ocorrência.")
  inv <- inv %>% distinct(gauge_code, .keep_all=TRUE)
}

inv_join <- inv %>%
  select(-any_of(c("estado","lat","long","network"))) %>%
  rename(time_step_inv=time_step)

# =============================================================================
# 3. RASTREIO DAS ESTAÇÕES
# =============================================================================

rastreio <- ini %>%
  filter(is.finite(mes_inicio), mes_inicio >= 1, mes_inicio <= 12) %>%
  select(gauge_code,estado,lat,long,time_step,network,classe_espacial,
         mes_inicio,fonte_inicio,mes_moda,f_moda,rbar_local,dist_nn_km) %>%
  left_join(inv_join, by="gauge_code") %>%
  mutate(
    gauge_code=trimws(as.character(gauge_code)),
    mes_inicio=as.integer(mes_inicio),
    mes_inicio_nome=meses_pt[mes_inicio],
    metodo_ano_hidro="uniplu_vizinhanca",
    time_step=dplyr::coalesce(as.numeric(time_step), as.numeric(time_step_inv)),
    network=toupper(trimws(as.character(network))),
    responsible=toupper(trimws(as.character(responsible))),
    # CRÍTICO: recuperar CEMADEN/INMET/ANA quando responsible estiver ausente
    responsible=ifelse(is.na(responsible) | !nzchar(responsible) | responsible=="NA",
                       network, responsible),
    responsible=ifelse(is.na(responsible) | !nzchar(responsible), "NA", responsible)
  ) %>%
  select(-any_of("time_step_inv"))

ufs <- sort(unique(na.omit(rastreio$estado)))
message("UFs: ", paste(ufs, collapse=", "))
message("Postos com mês inicial: ", nrow(rastreio))
print(rastreio %>% count(network,responsible,time_step,sort=TRUE), n=Inf)

saveRDS(rastreio, file.path(DIR_DF, "df_subdiario_mes_inicio_imax.rds"))
write_excel_csv(rastreio, file.path(DIR_DF, "df_subdiario_mes_inicio_imax.csv"))

# =============================================================================
# 4. PROCESSAMENTO POR UF
# =============================================================================

all_imax <- list()

for(uf in ufs){
  
  f_ckpt <- file.path(DIR_CACHE, paste0(uf, "_imax.rds"))
  
  if(file.exists(f_ckpt)){
    message("\n[",uf,"] checkpoint encontrado — carregando")
    all_imax[[uf]] <- readRDS(f_ckpt)
    next
  }
  
  rast_uf <- rastreio %>% filter(estado == uf)
  data_files <- list.files(file.path(DIR_SUBDIARIO,uf),
                           pattern="_data\\.parquet$", full.names=TRUE)
  
  if(!length(data_files)){
    message("[",uf,"] sem arquivos subdiários.")
    next
  }
  
  message("\n==================================================")
  message("[",uf,"] arquivos=",length(data_files)," | postos=",nrow(rast_uf))
  message("==================================================")
  
  chunks <- lapply(data_files, function(f){
    d <- tryCatch(read_parquet(f), error=function(e){
      message("[",uf,"] erro lendo ",basename(f),": ",conditionMessage(e)); NULL
    })
    if(is.null(d)) return(NULL)
    d %>% transmute(
      gauge_code=trimws(as.character(gauge_code)),
      rain_mm=as.numeric(rain_mm),
      datetime=as.POSIXct(datetime,tz="UTC")
    )
  })
  
  raw <- bind_rows(chunks)
  rm(chunks); invisible(gc())
  
  meta <- rast_uf %>%
    select(gauge_code,time_step,responsible,mes_inicio) %>%
    distinct(gauge_code,.keep_all=TRUE)
  
  raw <- raw %>%
    filter(gauge_code %in% meta$gauge_code) %>%
    left_join(meta,by="gauge_code") %>%
    filter(is.finite(time_step),time_step>0,!is.na(datetime)) %>%
    arrange(gauge_code,datetime)
  
  if(!nrow(raw)){
    message("[",uf,"] sem dados válidos.")
    next
  }
  
  message("[",uf,"] linhas=",format(nrow(raw),big.mark="."),
          " | postos=",n_distinct(raw$gauge_code))
  print(raw %>% distinct(gauge_code,time_step,responsible) %>%
          count(responsible,time_step,sort=TRUE), n=Inf)
  
  # -------------------------------------------------------------------------
  # Preenchimento: NÃO aplicar QC civil aqui
  # -------------------------------------------------------------------------
  
  message("[",uf,"] preenchendo grades temporais...")
  filled <- fun_filter_set(
    data=raw, daily=FALSE,
    col_names=c("gauge_code","rain_mm","datetime","time_step","responsible"),
    filter=FALSE
  )
  rm(raw); invisible(gc())
  
  if(is.null(filled) || !nrow(filled)){
    message("[",uf,"] nenhuma série preenchida.")
    next
  }
  
  # fun_filter_set antiga devolve time_step em segundos
  filled <- filled %>% mutate(time_step=as.numeric(time_step)/60)
  message("[",uf,"] resoluções: ",
          paste(sort(unique(filled$time_step)),collapse=", ")," min")
  
  # -------------------------------------------------------------------------
  # Agrupar por estação e resolução
  # -------------------------------------------------------------------------
  
  data.ls <- split(filled,filled$gauge_code)
  rm(filled); invisible(gc())
  
  data.by.time <- fun_group_ts(data.ls,ts_name="time_step")
  time.steps <- names(data.by.time)
  start_map <- setNames(as.integer(rast_uf$mes_inicio),rast_uf$gauge_code)
  
  # -------------------------------------------------------------------------
  # IMAX por resolução
  # -------------------------------------------------------------------------
  
  imax.ls <- lapply(time.steps,function(ts){
    
    current.ts <- data.by.time[[ts]]
    ds <- as.numeric(ts)/60
    ratio <- DURATIONS_HR/ds
    valid.durations <- DURATIONS_HR[abs(ratio-round(ratio)) <= 1e-8]
    if(!length(valid.durations)) return(NULL)
    
    message("[",uf,"] ts=",ts," min | postos=",length(current.ts),
            " | durações=",length(valid.durations))
    
    current.imax <- lapply(current.ts,function(df)
      df[,c("datetime","rain_mm"),drop=FALSE])
    names(current.imax) <- names(current.ts)
    
    sm <- start_map[names(current.imax)]
    keep <- is.finite(sm)
    current.imax <- current.imax[keep]
    sm <- sm[keep]
    names(sm) <- names(current.imax)
    if(!length(current.imax)) return(NULL)
    
    out <- imax.wateryear(
      data=current.imax,
      durations=valid.durations,
      start_month=sm,
      which.mon=1:12,
      names=c("datetime","rain_mm"),
      na_accept=NA_ACCEPT_HYDRO,
      apply_qc=APPLY_QC_HYDRO
    )
    
    if(is.null(out) || !nrow(out)) return(NULL)
    out$time_step_min <- as.numeric(ts)
    out$qc_na_accept <- NA_ACCEPT_HYDRO
    out$qc_applied <- APPLY_QC_HYDRO
    out
  })
  
  imax.ls <- imax.ls[!vapply(imax.ls,is.null,logical(1))]
  if(!length(imax.ls)){
    message("[",uf,"] nenhum IMAX gerado.")
    next
  }
  
  imax_uf <- bind_rows(imax.ls)
  imax_uf$estado <- uf
  
  qc_uf <- imax_uf %>%
    distinct(gauge_code,wateryear,start_month,na_prct_yr,
             complete_wateryear,qc_year)
  
  message("[",uf,"] IMAX=",format(nrow(imax_uf),big.mark="."),
          " | postos=",n_distinct(imax_uf$gauge_code),
          " | posto-anos aceitos=",nrow(qc_uf))
  
  saveRDS(imax_uf,f_ckpt)
  all_imax[[uf]] <- imax_uf
  
  rm(data.ls,data.by.time,imax.ls,imax_uf,qc_uf)
  invisible(gc())
}

# =============================================================================
# 5. CONSOLIDAÇÃO BRASIL
# =============================================================================

df_imax <- bind_rows(all_imax)
if(!nrow(df_imax)) stop("Nenhum IMAX foi gerado.")

df_imax <- df_imax %>%
  left_join(
    rastreio %>%
      select(gauge_code,lat,long,time_step,mes_inicio,mes_inicio_nome,
             network,responsible,classe_espacial,fonte_inicio,
             metodo_ano_hidro,f_moda,rbar_local,dist_nn_km),
    by="gauge_code", suffix=c("","_meta")
  ) %>%
  mutate(
    mes_inicio=dplyr::coalesce(as.integer(mes_inicio),as.integer(start_month)),
    qc_na_accept=as.numeric(qc_na_accept),
    qc_applied=as.logical(qc_applied)
  )

n_div <- df_imax %>%
  filter(is.finite(start_month),is.finite(mes_inicio),start_month!=mes_inicio) %>%
  nrow()
if(n_div>0) warning(n_div," linhas com start_month != mes_inicio.")

# =============================================================================
# 6. SALVAR BASE IMAX
# =============================================================================

F_OUT_RDS <- file.path(DIR_DF,"df_imax_subdiario.rds")
F_OUT_PQT <- file.path(DIR_DF,"df_imax_subdiario.parquet")
F_OUT_CSV <- file.path(DIR_DF,"df_imax_subdiario.csv")

saveRDS(df_imax,F_OUT_RDS)
write_parquet(df_imax,F_OUT_PQT)

df_imax %>%
  mutate(across(any_of(c("date","wy_start","wy_end")),as.character)) %>%
  write_excel_csv(F_OUT_CSV)

# =============================================================================
# 7. RESUMO DO QC
# =============================================================================

qc_final <- df_imax %>%
  distinct(gauge_code,estado,wateryear,start_month,na_prct_yr,
           complete_wateryear,qc_year,n_expected_yr,n_present_yr,
           n_na_present_yr,n_missing_structural_yr)

message("\n==================================================")
message("PASSO 12 — BRASIL CONCLUÍDO")
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
# 24 durações de 1–24 h e >= 8 anos em TODAS
# =============================================================================

sel_hidro_qc20 <- df_imax %>%
  filter(d %in% 1:24,is.finite(imax)) %>%
  distinct(gauge_code,estado,d,wateryear) %>%
  count(gauge_code,estado,d,name="n_anos_validos") %>%
  group_by(gauge_code,estado) %>%
  summarise(
    n_duracoes=n_distinct(d),
    n_duracoes_8anos=sum(n_anos_validos>=8),
    min_anos=min(n_anos_validos),
    med_anos=median(n_anos_validos),
    max_anos=max(n_anos_validos),
    .groups="drop"
  ) %>%
  mutate(passa=n_duracoes==24 & n_duracoes_8anos==24) %>%
  arrange(desc(passa),estado,desc(min_anos))

estacoes_aprovadas <- sel_hidro_qc20 %>% filter(passa)

message("\n==================================================")
message("SELEÇÃO 24 DURAÇÕES × >= 8 ANOS")
message("==================================================")
message("Estações aprovadas: ",nrow(estacoes_aprovadas))
print(estacoes_aprovadas,n=Inf,width=Inf)

write_excel_csv(sel_hidro_qc20,
                file.path(DIR_DF,"selecao_estacoes_hidro_qc20.csv"))

# =============================================================================
# 9. RESUMO POR UF
# =============================================================================

resumo_uf <- sel_hidro_qc20 %>%
  group_by(estado) %>%
  summarise(
    n_avaliadas=n(),
    n_aprovadas=sum(passa),
    .groups="drop"
  ) %>%
  arrange(desc(n_aprovadas))

print(resumo_uf,n=Inf)
write_excel_csv(resumo_uf,file.path(DIR_DF,"resumo_selecao_por_uf.csv"))

# =============================================================================
# 10. AUDITORIA CEMADEN
# =============================================================================

diag_cemaden_final <- df_imax %>%
  filter(responsible=="CEMADEN") %>%
  distinct(gauge_code,estado,wateryear,na_prct_yr) %>%
  group_by(gauge_code,estado) %>%
  summarise(
    n_anos_aceitos=n(),
    na_mediano_pct=100*median(na_prct_yr,na.rm=TRUE),
    na_max_pct=100*max(na_prct_yr,na.rm=TRUE),
    .groups="drop"
  )

message("\nCEMADEN: ",n_distinct(diag_cemaden_final$gauge_code),
        " postos com algum ano aceito | ",
        sum(diag_cemaden_final$n_anos_aceitos>=8),
        " com >=8 anos.")

write_excel_csv(diag_cemaden_final,
                file.path(DIR_DF,"auditoria_cemaden_hidro_qc20.csv"))

message("\nArquivos principais:")
message(F_OUT_RDS)
message(F_OUT_PQT)
message(F_OUT_CSV)
message("Cache por UF: ",DIR_CACHE)