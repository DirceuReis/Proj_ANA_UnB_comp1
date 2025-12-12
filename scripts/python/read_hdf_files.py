# Pacotes
import pandas as pd
import nexusformat.nexus as nx
import pyarrow as pa
import pyarrow.parquet as pq
import glob
import os

# Diretório dos arquivos
daily_path_ana = "base/fonte/bruto/diario/ANA/"
daily_path_inmet = "base/fonte/bruto/diario/INMET/"

# Ler caminho completo arquivos ANA e INMET
files_daily_ana = glob.glob(os.path.join(daily_path_ana, '*.h5'))
files_daily_inmet = glob.glob(os.path.join(daily_path_inmet, '*.h5'))

# Visualizar estrutura do HDF
hdf_file = nx.nxload(files_daily_inmet[1])
print(hdf_file.tree)

# Importar e salvar arquivos dentro de uma lista
table_data_daily_ana = {}   # lista vazia p/ armazenar 'table_data'
table_info_daily_ana = {}   # lista vazia p/ armazenar 'table_info'
table_data_daily_inmet = {} # lista vazia p/ armazenar 'table_data
table_info_daily_inmet = {} # lista vazia p/ armazenar 'table_info'

# Ler séries ANA
for file in files_daily_ana:
  key = os.path.basename(file)
  table_data_daily_ana[key] = pd.read_hdf(file, "table_data")
  table_info_daily_ana[key] = pd.read_hdf(file, "table_info")

# Ler séries INMET
for file in files_daily_inmet:
  key = os.path.basename(file)
  table_data_daily_inmet[key] = pd.read_hdf(file, "table_data")
  table_info_daily_inmet[key] = pd.read_hdf(file, "table_info")

# Juntar dados
table_data_ana = pd.concat(table_data_daily_ana)     # 'table_data' ANA
table_info_ana = pd.concat(table_info_daily_ana)     # 'table_info' ANA
table_data_inmet = pd.concat(table_data_daily_inmet) # 'table_data' inmet
table_info_inmet = pd.concat(table_info_daily_inmet) # 'table_info' inmet

# Adicionar coluna 'responsible' aos 'table_data' e juntar
table_data_ana["responsible"] = "ANA"
table_data_inmet["responsible"] = "INMET"

# Juntar dados data e info
table_data = pd.concat([table_data_ana, table_data_inmet])
table_data["datetime"] = pd.to_datetime(table_data["datetime"])
table_info = pd.concat([table_info_ana, table_info_inmet])
table_data_arrow = pa.Table.from_pandas(table_data)
table_info_arrow = pa.Table.from_pandas(table_info)

# Salvar arquivos em PARQUET
pq.write_table(table_data_arrow, "base/gerados/df_daily_data.parquet")
pq.write_table(table_info_arrow, "base/gerados/df_daily_info.parquet")
