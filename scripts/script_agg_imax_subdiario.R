
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); gc()

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,    # gerenciador de pacotes
               arrow,     # leitura de dados 'parquet'
               lubridate, # manipulação de datas
               parallel,
               beepr      # alerta sonoro
               )


# LER DADOS ---------------------------------------------------------------

# Importar
df_data <- arrow::read_parquet("base/gerados/df_subdaily_data.parquet")

# Critérios p/ filtrar séries
na_accept <- 0 # 0.2
min_years <- 0 # 8

# Funções
source("scripts/funcoes/fun_filter_set.R")        # preencher datas e filtrar
source("scripts/funcoes/fun_group_by_timestep.R") # agrupar por time_step
source("scripts/funcoes/fun_imax_agg.R")          # agregar por durações e calcular intensidade máxima anual

# Preencher datas e filtrar conforme 'na_accept' e 'min_years'
df_data <- fun_filter_set(data = df_data, filter = FALSE); invisible(gc())

# AGREGAR DURAÇÕES --------------------------------------------------------

# Agrupar séries por 'time_steps'
data.ls <- split(df_data, f = df_data$gauge_code)             # lista de data.frames de duração
data.by.time <- fun_group_ts(data.ls, ts_name = "time_steps") # agrupar em durações diferentes numa lista

# Durações
d.subdaily <- c(1/6, 1/4, 0.5, seq(1, 23)) # durações subdiárias
d.daily <- c(1, 2, 3, 4, 5, 7, 10)         # durações diárias
durations <- c(d.subdaily, d.daily*24)     # durações
time.steps <- names(data.by.time)          # lista de time_steps

# Aplicar função 'fun_imax_agg' p/ cada grupo de durações
imax.ls <- lapply(X = time.steps, FUN = function(ts){
  
  current.ts <- data.by.time[[ts]] # extrair listas de estações c/ resolução 'ts'
  ds <- as.numeric(ts)/3600        # segundos p/ horas
  valid.durations <- durations[(durations/ds) %% 1 < 1e-8] # selecionar somente durações compatíveis c/ a resolução 'ts'
  
  message(paste0("\nProcessando conjunto de estações com resolução de ", ds*60, " minutos..."))
  
  out <- fun_imax_agg(data = current.ts,
                      durations = valid.durations,
                      which.mon = 1:12,
                      names = c("datetime", "rain_mm"))
  
}) # fim 'imax.ls'

rm(df_data, data.ls, data.by.time); invisible(gc())

# Nomear lista
names(imax.ls) <- time.steps

# Salvar resultados
# saveRDS(object = imax.ls, file = "base/gerados/imax_ls.rds")
arrow::write_parquet(x = bind_rows(imax.ls), sink = "base/gerados/df_imax.parquet")

# Ler arquivo c/ intensidades máximas
df_imax <- arrow::read_parquet(file = "base/gerados/df_imax.parquet")
