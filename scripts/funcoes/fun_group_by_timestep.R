# Essa função agrupa uma lista com séries temeporais em diferentes resoluções
# temporais em listas 

fun_group_ts <- function(ls, ts_name){
  
  time_steps <- sapply(ls, function(df){
    ts <- unique(df[[ts_name]])
    if(length(ts) > 1) warning("Um ou mais data.frames possuem múltiplos 'time_steps'. Utilizando primeiro valor.")
    return(ts[1])
  })
  
  grouped_list <- split(ls, time_steps)
  return(grouped_list)
  
}