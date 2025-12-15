
# COMENTÁRIOS -------------------------------------------------------------

# Este script realiza a leitura e organização das séries de precipitação diária consolidadas.
# Pela UFPB nos foi fornecido um conjunto de séries de máximos anuais das estações com qualidade
# alta, armazenadas em um arquivo Excel de nome 'Máximos anuais para IDF.xlsx'.

# Como precisamos das séries completas, olharemos também para as séries brutas, e as filtraremos
# de acordo com os anos e as estações contidas em 'Máximos anuais para IDF.xlsx'.
# Esses dados brutos estão armazenados em 'df_daily_data.parquet', um arquivo gerado a partir
# dos arquivos .HDF disponibilizados pela UFPB.

# Esses dados também precisarão ter suas séries de datas completadas, para que seja possível visualizar
# buracos existentes onde não houve registro.

# PACOTES -----------------------------------------------------------------

rm(list = ls()); invisible(gc())

if(!require(pacman)){
  message("Instalando gerenciador de pacotes 'pacman'...")
  install.packages("pacman")
}

# Pacotes
pacman::p_load(pacman, tidyverse, arrow)


# FUNÇÕES -----------------------------------------------------------------

# Função p/ completar datas das séries
source("scripts/funcoes/fun_fill_dates.R")


# Caminhos
path.pdmax <- "base/fonte/consolidado/diario/Máximos anuais para IDF.xlsx" # precipitação diária máxima anual
path.bruto <- "base/fonte/bruto/diario/df_daily_data.parquet" # séries completas estações diárias
path.info <- "base/fonte/bruto/diario/df_daily_info.parquet"  # outras informações estações diárias

# Ler arquivos
df_daily_max <- readxl::read_excel(path.pdmax)   # precipitação diária máxima anual
df_daily_data <- arrow::read_parquet(path.bruto) # séries completas estações diárias
df_daily_info <- arrow::read_parquet(path.info)  # informações estações diárias
df_daily_data <- df_daily_data[-c(5,6)]          # remover últimas duas colunas
df_daily_info <- df_daily_info[-c(8,9)]          # remover últimas duas colunas


# ORGANIZAR DADOS ---------------------------------------------------------

# Código estações consolidadas
gauges <- unique(df_daily_max$gauge_code)

# Filtrar séries completas e reforçar classes
df_daily_data <- df_daily_data %>% 
  filter(gauge_code %in% gauges) %>% 
  mutate(time_steps = 24*60)

# Organizar colunas
df_daily_data <- df_daily_data[c("gauge_code", "rain_mm", "datetime", "time_steps", "responsible")]

ls.daily.data <- split(x = df_daily_data, f = df_daily_data$gauge_code)
ls.daily.max <- split(x = df_daily_max, f = df_daily_max$gauge_code)

ls.daily.data <- lapply(X = gauges, FUN = function(gauge){
  
  df.data <- ls.daily.data[[gauge]] # extrair séries
  df.max <- ls.daily.max[[gauge]]   # atuais
  
  years.max <- df.max$year                        # extrair anos série consolidado
  years.data <- lubridate::year(df.data$datetime) # extrair anos série bruta
  years.idx <- years.data %in% years.max          # extrair posição dos anos comuns
  df.data <- df.data[years.idx,]                  # filtrar série bruta                 
  
}) # fim 'ls.daily.data'
# Warning: tz(): Don't know how to compute timezone for object of class NULL; returning "UTC". 

names(ls.daily.data) <- gauges

# Transformar em 'tbl_df'
df.daily.data <- bind_rows(ls.daily.data)
rm(ls.daily.max, ls.daily.data); invisible(gc())

# O arquivo 'df_info' contém estações repetidas, portanto precisamos manter
# somente um único registro de cada 'gauge_code'
df.daily.info <- slice_tail(.data = df_daily_info, n = 1, by = "gauge_code")


# COMPLETAR DATAS ---------------------------------------------------------

ls.daily.filled <- fun_fill_dates(df = df.daily.data, col_names = names(df.daily.data), daily = TRUE)
df.daily.filled <- bind_rows(ls.daily.filled)


# VISUALIZAR --------------------------------------------------------------

# Gráfico com séries
n.gauges <- length(gauges)
ggplot(df.daily.filled, aes(x = datetime, y = rain_mm)) +
  facet_wrap(~gauge_code) +
  geom_line(linewidth = 0.05) +
  geom_vline(data = df.daily.data[is.na(df.daily.data$rain_mm),], aes(xintercept = datetime), color = "yellow", alpha = 0.5, linewidth = 0.01) +
  theme_minimal() +
  theme(plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing.y = unit(0.1, "lines"),
        strip.text = element_text(size = 6),
        text = element_text(size = 8, family = "serif", color = "black"))


# SALVAR ------------------------------------------------------------------

# Parquet
arrow::write_parquet(df.daily.info, sink = "base/gerados/df_daily_info.parquet")
arrow::write_parquet(df.daily.data, sink = "base/gerados/df_daily_data.parquet")
arrow::write_parquet(df.daily.filled, sink = "base/gerados/df_daily_filled.parquet")
