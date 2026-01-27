fn1.map <- function(R.boot, gauge){
  
  map.teste <- gauge %>%
    select(gauge_code, d, imax, year) %>%
    group_by(gauge_code) %>%
    nest(data = c(d, year, imax)) %>% # criar um tibble com informações p/ fun.ml.boot
    mutate(
      boot = purrr::map(.x = data,                                # criar coluna c/ resultados
                        .f = ~boot::boot(data = .x,               # usar dados em .x = data
                                         statistic = fun.lm.boot, # função a ser calculada
                                         R = R.boot,              # nro. replicas
                                         strata = .x$d) %>%       # reamostrar por durações
                          broom::tidy(conf.int = TRUE) %>%        # extrair infos e intervalos de confiança
                          mutate(d = first(.x$d),                 # adicionar coluna duração
                                 n_year = n_distinct(.x$year),    # adicionar coluna nro. anos
                                 which_moment = row_number()))) %>%  # adicionar coluna momento calculado
    unnest(boot) %>%                                              # tirar 'boot' do nest
    rename_with(~gsub(".", "_", .x, fixed = TRUE)) %>%            # renomear colunas
    rename(scale = statistic, ci_lower = conf_low, ci_upper = conf_high) %>%
    select(-data)
  
}

teste.fn1.map <- fn1.map(R.boot, gauge)

fn2.nomap <- function(R.boot, gauge){
  
  boot.obj <- boot::boot(data = gauge,
                         statistic = fun.lm.boot,
                         R = R.boot,
                         strata = gauge$d)
  
  boot.ci <- broom::tidy(boot.obj, conf.int = TRUE) %>% 
    mutate(d = first(gauge$d),
           n_year = n_distinct(gauge$year),
           which_moment = row_number())
  
}

teste.fn2.nomap <- fn2.nomap(R.boot, gauge)

pacman::p_load(rbenchmark)

boot.tests <- list(map = expression(fn1.map(R.boot, gauge)),
                   nomap = expression(fn2.nomap(R.boot, gauge)))

R.boot <- 10000

benchmark <- do.call(rbenchmark::benchmark,
                     c(boot.tests, list(replications = 1,
                                        columns = c("test", "elapsed", "replications"),
                                        order = "elapsed"))); benchmark
