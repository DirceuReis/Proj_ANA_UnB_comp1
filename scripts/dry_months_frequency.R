
# MÁXIMOS ANUAIS ----------------------------------------------------------

# Pacotes
pacman::p_load(pacman, tidyverse, arrow)

# Funções
source("scripts/funcoes/fun_imax_agg_wateryear.R")
source("scripts/funcoes/fun_filter_set.R")
source("scripts/funcoes/fun_group_ts.R")

# Contar frequência de meses
parse.months <- function(data.d){
  
  df <- data.d %>% 
    select(-c(mon_filter, n_max)) %>% 
    mutate(month = factor(month(date))) %>% 
    group_by(gauge_code, wateryear, month) %>% 
    reframe(mon_freq = n())
  
  m <- as.matrix(table(df$gauge_code, df$month))
  
}

# Ler dados
# Rodar abaixo se não encontrar arquivo "df_imax_wyr.pqt"
if(file.exists("base/gerados/df_imax_.pqt")){
 
  df.imax <- arrow::read_parquet("base/gerados/df_imax_.pqt") # alterar p/ janeiro/agosto/outubro
  
}else{
   
  df.data <- arrow::read_parquet("base/gerados/df_subdaily_data.pqt")
  
  # Preencher datas
  df.data <- fun_filter_set(data = df.data,
                            daily = FALSE,
                            col_names = names(df.data),
                            filter = FALSE); invisible(gc())
  
  # Criar lista por "time_step"
  ls.data <- split(df.data, f = df.data$gauge_code)
  data.by.time <- fun_group_ts(ls.data, ts_name = "time_step"); invisible(gc())
  
  d.subhour <- c(10, 15, 20, 30, 40, 45, 50)
  d.subdaily <- seq(1, 23)
  d.daily <- 1:10
  durations <- c(d.subhour/60, d.subdaily, d.daily*24)
  time.steps <- names(data.by.time)
  wateryear <- c(8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7)
  
  # Extrair intensidades máximas anuais
  imax.ls <- lapply(X = time.steps, FUN = function(ts){
    
    current.ts <- data.by.time[[ts]]
    ds <- as.numeric(ts)/3600
    valid.durations <- durations[(durations/ds) %% 1 < 1e-8]
    
    message(paste0("Processando conjunto de estações com resolução de ", ds*60, " minutos..."))
    
    out <- imax.wateryear(data = current.ts,
                          durations = valid.durations,
                          which.mon = wateryear,
                          names = c("datetime", "rain_mm")) 
    
  })
  
  df.imax <- dplyr::bind_rows(imax.ls)
  filename <- paste0("base/gerados/df_imax_", lubridate::month(wateryear[1], label = TRUE), ".pqt")
  arrow::write_parquet(dplyr::bind_rows(imax.ls), filename)
} 


# FREQUÊNCIA DE NÃO OCORRÊNCIA --------------------------------------------

# Ler dados df.imax do ano civil
df.imax.janeiro <- arrow::read_parquet("base/gerados/analise_sazonalidade/df_imax_jan.pqt")

# Calcular a frequência de meses em que ocorre chuva
# Fazer uma lista por duração
imax.ds <- split(df.imax, df.imax$d)
imax.month <- lapply(imax.ds, function(data.d){
  
  m <- parse.months(data.d)                         # calcular freq meses por estação
  sum.months <- apply(X = m, MARGIN = 2, FUN = sum) # somar por colunas (meses)
  df <- rownames_to_column(as_tibble(sum.months), var = "month") %>% 
    mutate(month = factor(month, levels = 1:12),
           d = data.d$d[1],
           freq = value/sum(value)) %>%
    arrange(month) %>% 
    group_by(d) %>% 
    mutate(is_min = ifelse(value == min(value), TRUE, FALSE)) %>% 
    ungroup()
  
}) %>% bind_rows()

# Durações de interesse
durations <- as.character(c(1, 6, 12, 24, 72))

ggplot(imax.month[imax.month$d %in% durations,], aes(x = month, y = value)) +
  facet_wrap(~d, ncol = 1, scales = "free_y", labeller = labeller(d = function(x) paste0(x, " h"))) +
  geom_col(aes(fill = is_min), color = "black", linewidth = 0.3) +
  geom_text(aes(label = scales::label_percent(accuracy = 0.1)(freq)),
            size = 3, family = "Times New Roman", vjust = -0.5) +
  scale_fill_manual(values = c("TRUE" = "darkgoldenrod2", "FALSE" = "cadetblue3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Meses", y = "Número de meses com `imax`") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black"),
        legend.position = "non",
        text = element_text(family = "Times New Roman", color = "black"))

ggsave(filename = "figuras/sazonalidade/plot_freq_meses.png",
       height = 20, width = 16, units = "cm", dpi = 300)


# DIFERENTES ANOS HIDROLÓGICOS --------------------------------------------

df.imax.janeiro <- arrow::read_parquet("base/gerados/df_imax_janeiro.pqt")
df.imax.outubro <- arrow::read_parquet("base/gerados/df_imax_outubro.pqt")
df.imax.agosto <- arrow::read_parquet("base/gerados/df_imax_agosto.pqt")

parse.months <- function(data.d){
  
  df <- data.d %>% 
    select(-c(mon_filter, n_max)) %>% 
    mutate(month = factor(month(date))) %>% 
    group_by(gauge_code, wateryear, month) %>% 
    reframe(mon_freq = n())
  
  m <- as.matrix(table(df$gauge_code, df$month))
  
}

imax.wyr <- list(
  df.imax.civil %>% mutate(ano = "Janeiro", wateryear = year),
  df.imax.outubro %>% mutate(ano = "Outubro"),
  df.imax.agosto %>% mutate(ano = "Agosto")
)

imax.month <- lapply(imax.wyr, function(df.imax){
  
  imax.ds <- split(df.imax, df.imax$d)
  imax.month <- lapply(imax.ds, function(data.d){
    
    # Em todos os anos c/ imax == 0, a data é sempre no primeiro
    # dia do ano hidrológico definido, então é necessário filtrar
    # imax > 0 p/ evitar essa distorção
    data.d <- data.d[data.d$imax > 0,]
    
    m <- parse.months(data.d)                         # calcular freq meses por estação
    sum.months <- apply(X = m, MARGIN = 2, FUN = sum) # somar por colunas (meses)
    df <- rownames_to_column(as_tibble(sum.months), var = "month") %>% 
      mutate(month = factor(month, levels = 1:12),
             d = data.d$d[1],
             ano = data.d$ano[1],
             freq = value/sum(value)) %>%
      arrange(month) %>% 
      group_by(d) %>% 
      mutate(is_min = ifelse(value == min(value), TRUE, FALSE)) %>% 
      ungroup()
    
  }) %>% bind_rows()
  
}) %>% bind_rows() %>% 
  mutate(ano = factor(ano, levels = c("Janeiro", "Agosto", "Outubro")))

# Remover falhas 

unique(imax.month$ano)

durations <- as.character(c(1, 6, 12, 24, 72))

ggplot(imax.month[imax.month$d %in% durations,], aes(x = month, y = value, fill = ano)) +
  facet_wrap(~d, ncol = 1, scales = "free_y", labeller = labeller(d = function(x) paste0(x, " h"))) +
  geom_col(color = "black", position = position_dodge(width = 0.9), linewidth = 0.3) +
  # geom_text(aes(label = value), position = position_dodge(width = 1), size = 2.5, family = "Consolas", vjust = -0.5) +
  geom_text(data = subset(imax.month[imax.month$d %in% durations,], is_min == TRUE),
            aes(x = month, y = value, label = "▼", group = ano), color = "red3",
            position = position_dodge(width = 0.99), vjust = -0.6, hjust = 0.5, size = 2.5, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_manual(values = c("Janeiro" = "grey10", "Agosto" = "grey50", "Outubro" = "grey90")) +
  labs(x = "Meses", y = "Número de meses com `imax`", fill = "") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black"),
        legend.position = "bottom",
        text = element_text(family = "Times New Roman", color = "black"))

ggsave(filename = "figuras/sazonalidade/plot_freq_meses.png",
       height = 20, width = 16, units = "cm", dpi = 300)
