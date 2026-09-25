
# Este script lê e salva os dados brutos por estado, separados em anos, em um único data.frame
# em arquivo .PARQUET/.PQT


# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Instalar e carregar pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  pacman,
  tidyverse, # manipulação de dados
  arrow     # ler dados no formato 'parquet'
)


# FUNÇÃO -----------------------------------------------------------------

# Ler, juntar e salvar tabela única por estado para 'info' e 'data' em .PARQUET
read_subdaily <- function(path, outdir, ext = c(".pqt", ".parquet")){

  files <- list.files(path, pattern = ext)
  ext <- match.arg(ext, c(".pqt", ".parquet"))
  if(length(files) == 0) stop("Nenhum arquivo encontrado com a extensão especificada. Defina 'ext' como '.pqt' ou '.parquet'.")
  
  # Listar arquivos
  files.data <- list.files(path = file.path(path), pattern = paste0("*.data", ext), full.names = TRUE) # listar 'data'
  files.info <- list.files(path = file.path(path), pattern = paste0("*.info", ext), full.names = TRUE) # listar 'info

  # Ler dados
  data.ls <- lapply(X = files.data, FUN = arrow::read_parquet)
  info.ls <- lapply(X = files.info, FUN = arrow::read_parquet)

  # Juntar dados
  data.df <- bind_rows(data.ls)
  info.df <- bind_rows(info.ls)

  # O arquivo 'df_info' contém estações repetidas, portanto precisamos manter
  # somente um único registro de cada 'gauge_code'
  info.df <- slice_tail(.data = info.df, n = 1, by = "gauge_code")
  uf <- unique(info.df$state)

  # Fuso horário
  # Todos os dados do conjunto deveriam estar na zona UTC -3
  # Isso corrige o shift de -2 ou -3 horas nos dados
  data.df <- data.df |> 
    mutate(
      datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"), # alterar fuso
      time_step = info.df$time_step[match(gauge_code, info.df$gauge_code)],      # resolução temporal
      responsible = info.df$responsible[match(gauge_code, info.df$gauge_code)],  # operador
      state = uf
    ); invisible(gc()) # resolução temporal
  
  # Verificar diretório de saída
  if(!dir.exists(outdir)){
    warning(sprintf("Diretório 'outdir' não existe. Criando '%s'", outdir))
    dir.create(outdir)
  }

  # Parquet
  write_parquet(data.df, sink = paste0(outdir, tolower(uf), "_subdaily_data.pqt"), compression = "gzip", compression_level = 5)
  write_parquet(info.df, sink = paste0(outdir, tolower(uf), "_subdaily_info.pqt"), compression = "gzip", compression_level = 5)
  invisible(gc())
  message(sprintf("[%s] Séries subdiárias salvas em '%s'", uf, outdir))

}


# LER DADOS ---------------------------------------------------------------

# Caminho do diretório
path.raw <- "base/fonte/consolidado/subdiario_br/"
ufs <- list.dirs(path.raw, full.names = FALSE, recursive = FALSE) |> print()
# path <- paste0(path.raw, ufs[23])
outdir <- "Relatório 2/Scripts_Ajuste_dGEV/resultados/subdiarios/binded_raw/"

purrr::walk(.x = ufs, .f = \(uf) read_subdaily(path = paste0(path.raw, uf), outdir = outdir, ext = ".parquet"))
