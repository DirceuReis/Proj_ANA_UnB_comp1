# Essa função estima os intervalos de confiança usando o método Bootstrap
# com o pacote 'boot' para estimativas de L-curtose e L-assimetria do método dos
# L-momentos usando o pacote 'boot'.

fun_boot_lratio <- function(data,                  # tbl_df: c/ máximos anuais p/ durações diferentes
                            rep = 2000,            # nro. réplicas bootstrap
                            which_lratio = c(3,4), # quais lmom ratio estimar: um (ou todos) de c(3,4)
                            na_accept = 0.2,       # percentual máximo de falhas anual
                            min_years = 8,         # nro. mínimo de anos
                            signf = 0.05,          # nível de significância
                            ci_type = "all",       # norm, stud, basic, perc, bca
                            # return_rep = FALSE,    # retornar as réplicas bootstrap
                            col_names = c("imax", "d", "gauge_code")){
  
  # Pacotes
  if(!require(pacman)) install.packages("pacman")
  pacman::p_load(pacman, lmom, boot, pbapply)
  
  # Função que calcula os intervalos bootstrap p/ 1 estação
  fun_boot_aux <- function(data, rep, which_lratio, signf, ci_type, col_names){
    
    # Função p/ realizar a reamostragem e calcular a estatística de interesse (lratio)
    fun_lratios <- function(data, idx, which_lratio){
      
      # Reamostragem
      d <- data[idx,]
      
      # Razões de momentos-L
      lratio <- lmom::samlmu(d[[col_names[1]]])[which_lratio] # extrair valores nas posicoes
      
    }
    
    # Repetir p/ cada duração
    ds <- unique(data[[col_names[2]]]) # extrair vetor de durações de cada estação
    ls_ci <- list()
    
    # Calcular ICs p/ cada duração
    for(d in ds){
      
      # Extrair imax p/ uma duração por vez
      df <- data[data[[col_names[2]]] == d,]
      
      # Bootstrap
      boot <- boot::boot(data = df,
                         statistic = fun_lratios,
                         R = rep,
                         which_lratio = which_lratio)
      
      # Estimar ICs
      for(i in seq_along(which_lratio)){
        
        # Nomes
        ratio_idx <- which_lratio[i]
        ratio_name <- paste0("tau", ratio_idx)
        
        # IC
        ci <- boot::boot.ci(boot.out = boot,
                            conf = 1 - signf,
                            type = ci_type,
                            index = i)
        
        # Valor da estimativa original
        t0 <- ci$t0
        
        # Extrair IC p/ tipo especificado em "ci_type"
        if(ci_type == "norm") {
          ci_bounds <- ci$normal[2:3]
        } else if(ci_type == "perc") {
          ci_bounds <- ci$percent[4:5]
        } else if(ci_type == "basic") {
          ci_bounds <- ci$basic[4:5]
        } else if(ci_type == "stud") {
          ci_bounds <- ci$student[4:5]
        } else if(ci_type == "bca") {
          ci_bounds <- ci$bca[4:5]
        }
        
        # Resultados linha
        aux <- data.frame(gauge_code = df[[col_names[3]]][1],
                          ds = d,
                          lratio = ratio_name,
                          value = t0,
                          ci_l = ci_bounds[1],
                          ci_u = ci_bounds[2])
        
        ls_ci[[paste0("d", d, "_", ratio_name)]] <- aux
        
      }
      
    }
    
    # Resultados
    res <- do.call(rbind, ls_ci)
    rownames(res) <- NULL
    
    return(res)
    
  }
  
  # Aplicar a função fun_boot_aux p/ todas as estações de 'data'
  ls_data <- split(x = data, f = data[[col_names[3]]]) # dividir 'data' por estação
  res <- pbapply::pblapply(ls_data, function(df){
    
    tryCatch({
      
      fun_boot_aux(data = df,
                   rep = rep,
                   which_lratio = which_lratio,
                   signf = 0.05,
                   ci_type = ci_type,
                   col_names = col_names)
      
    }, error = function(e){
      
      message("Pulando estação '", unique(df[[col_names[3]]]), "':\n", e)
      return(NULL)
      
    })
    
  })
  
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  return(res)
  
}
