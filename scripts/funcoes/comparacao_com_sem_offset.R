rm(list = ls()); invisible(gc())

pacman::p_load(pacman, tidyverse, arrow)

df0 <- arrow::read_parquet("base/gerados/scale_invariance/df_scale_block_boot.pqt")
df1 <- arrow::read_parquet("base/gerados/scale_invariance/df_scale_block_boot_offset.pqt")

gauges <- unique(df0$gauge_code)
duration.intervals <- unique(df0$d)
gauge <- "430120601A" # estação p/ avaliação

df0 <- df0 %>%
  filter(gauge_code == gauge,
         d == duration.intervals[[1]])

df1 <- df1 %>%
  filter(gauge_code == gauge,
         which_moment == 1,
         d == duration.intervals[[1]])

df.mom <- arrow::read_parquet("base/gerados/df_imax.parquet") %>% 
  filter(gauge_code == gauge) %>% 
  group_by(d) %>% 
  summarise(mom = mean(imax))

scale_0 <- df0[df0$which_variable == "scale",][["statistic"]]
scale_1 <- df1[df1$which_variable == "scale",][["statistic"]]
offset_0 <- 0
offset_1 <- df1[df1$which_variable == "offset",][["statistic"]]

hat0 <- df.mom %>% 
  mutate(which_model = factor(0),
         log_id = log(mom),
         log_iD = log(max(mom)),
         scale = scale_0,
         offset = offset_0,
         hat_log_id = -scale*log(d + offset) + scale*log(max(d) + offset) + log_iD)

hat1 <- df.mom %>% 
  mutate(which_model = factor(1),
         log_id = log(mom),
         log_iD = log(max(mom)),
         scale = scale_1,
         offset = offset_1,
         hat_log_id = -scale*log(d + offset) + scale*log(max(d) + offset) + log_iD)

ggplot() +
  geom_line(aes(x = log(d + offset), y = hat_log_id, color = which_model), data = hat0) +
  geom_line(aes(x = log(d + offset), y = hat_log_id, color = which_model), data = hat1) +
  geom_point(aes(x = log(d), y = log_id), data = hat0)
