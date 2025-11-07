
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); gc()

# Instalar e carregar pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,
               tidyverse, # manipulação de dados
               lubridate, # manipulação de datas
               arrow,     # ler dados no formato 'parquet'
               ggplot2)   # visualização

# LER DADOS ---------------------------------------------------------------

# Caminho do diretório
path_dir <- "base/fonte/consolidado/subdiario/"


# Listar arquivos c/ extensão ".parquet" no diretório
# Os arquivos contém divisões por ano de duas tabelas: info e data

## Listar 'data'
files_data <- list.files(path = path_dir,           
                         pattern = "*.data.parquet",
                         full.names = TRUE)


## Listar 'info'
files_info <- list.files(path = path_dir,
                         pattern = "*.info.parquet",
                         full.names = TRUE)

# Ler dados
# list_data <- lapply(X = files_data, FUN = arrow::read_parquet, as_data_frame = TRUE)
list_data <- lapply(X = files_data, FUN = arrow::read_parquet)
list_info <- lapply(X = files_info, FUN = arrow::read_parquet)

# Juntar dados
df_data <- bind_rows(list_data)[-4]  # remover coluna '__intex_level_0__'
df_info <- bind_rows(list_info)[-10] # remover coluna '__intex_level_0__'

# Manter somente 'df_data' e 'df_info' no ambiente
rm(list = setdiff(ls(), c("df_data", "df_info"))); gc()


# ORGANIZAR ---------------------------------------------------------------

# O arquivo 'df_info' contém estações repetidas, portanto precisamos manter
# somente um único registro de cada 'gauge_code'
df_info <- slice_tail(.data = df_info, n = 1, by = "gauge_code")

# Fuso horário
# Todos os dados do conjunto deveriam estar na zona UTC -3
# Isso corrige o shift de -2 ou -3 horas nos dados
df_data <- df_data %>% 
  mutate(
    datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),        # alterar fuso
    time_steps = df_info$monitoring_time_step[match(gauge_code, df_info$gauge_code)], # resolução temporal
    responsible = df_info$responsible[match(gauge_code, df_info$gauge_code)]          # operador
  ); gc() # resolução temporal


# EXPORTAR BRUTO ----------------------------------------------------------

# Parquet
write_parquet(df_data, sink = "base/gerados/df_subdaily_data.parquet")
write_parquet(df_info, sink = "base/gerados/df_subdaily_info.parquet")

gc()

# PREENCHER DATAS ---------------------------------------------------------

# Ler dados gerados em `rs_subdaily_read_parquet.R`
df_data <- arrow::read_parquet("base/gerados/df_subdaily_data.parquet")
df_info <- arrow::read_parquet("base/gerados/df_subdaily_info.parquet")

# Função completa
source("scripts/funcoes/fun_fill_dates.R")
list_subdaily_filled <- fun_fill_dates(df = df_data, col_names = names(df_data)); gc()
df_subdaily_filled <- bind_rows(list_subdaily_filled); rm(list_subdaily_filled); gc()

# Salvar arquivo (deixar comentado)
list_by_resp <- split(df_subdaily_filled, f = df_subdaily_filled$responsible); rm(df_subdaily_filled); gc()
compression_codecs <- c("snappy", "gzip", "brotli", "zstd", "lz4", "lzo", "bz2")
write_parquet(x = df_subdaily_filled,
              sink = "base/gerados/df_subadaily_filled.parquet",
              compression = compression_codecs[[4]],
              compression_level = 9)


# ANÁLISE DAS SÉRIES ------------------------------------------------------

# Tabela resumo
df_summary <- 
  df_subdaily_filled %>% 
  mutate(year = lubridate::year(datetime),        # criar coluna c/ ano
         rain_na = is.na(rain_mm)) %>%            # TRUE/FALSE se rain_mm é NA
  group_by(gauge_code, year) %>%                  # agrupar por estacao e por ano
  summarise(rain_yr = sum(rain_mm, na.rm = TRUE), # acumulado anual
            n_fail = sum(rain_na),                # número de observações c/ NA por ano
            n_obs = n(),                          # número total de observações
            fail_prct = n_fail/n_obs*100)         # percentual de falhas

# Adicionar operador responsável
df_summary <- merge(df_summary, df_info[, c("gauge_code", "responsible")], by = "gauge_code", all.x = TRUE)


# VISUALIZAÇÃO ------------------------------------------------------------

# Analisar uma série em específico
gauge_cemaden <- "430920901A"
df_subdaily_filled %>% 
  filter(gauge_code == gauge_cemaden) %>% 
  ggplot(aes(x = datetime, y = rain_mm)) +
  geom_line(linewidth = 0.5) +
  labs(x = "", y = "Precipitação [mm]", col = "") +
  theme_minimal() +
  theme(legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        legend.key.width = unit(1, "cm")); gc()

# Qtas estações p/ uma dada combinação de critérios de qualidade
# Limites
thresholds <- seq(10, 40, 10)
min_years <- 1:10

# Avaliar quantas estacoes atendem os criterios de falha e nro mínimo de anos 
analise_series <-
  expand_grid(threshold = thresholds,                              # combinar possibilidades
              min_years = min_years) %>%                           # entre threshold e min_years
  mutate(n_gauges = purrr::map2_int(threshold, min_years, ~{  
    df_summary %>% 
      filter(fail_prct <= .x) %>% 
      count(gauge_code) %>% 
      filter(n >= .y) %>% 
      nrow()
  })) # extrair nro estacoe atendem criterios

# Visualização
cores <- c("green4", "lightgreen", "orange", "#F10000")
analise_series %>% 
  mutate(threshold = factor(threshold, levels = thresholds)) %>% 
  ggplot(aes(x = min_years, y = n_gauges, fill = threshold)) +
  geom_col() +
  geom_text(aes(label = n_gauges), vjust = -0.3, size = 3) +
  facet_wrap(~threshold, nrow = 4, labeller = labeller(threshold = function(x) paste0(x, "% falhas"))) +
  scale_x_continuous(breaks = seq(min(analise_series$min_years), max(analise_series$min_years), by = 1)) +
  scale_y_continuous(limits = c(0, 430)) +
  scale_fill_manual(values = cores) +
  labs(x = "Número de anos por estação",
       y = "Número de estações") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        legend.position =  "none",
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(size = 10))