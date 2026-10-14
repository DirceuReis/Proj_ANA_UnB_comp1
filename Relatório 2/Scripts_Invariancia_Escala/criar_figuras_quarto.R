library(stringr)
path.figuras <- "Relatório 2/Scripts_Invariancia_Escala/figuras/"
figuras <- list.files(path.figuras, pattern = "plot_invariancia_escala_momentos_*", full.names = TRUE)
ufs <- str_extract(figuras, "[^_]+(?=\\.[^.]*$)")

uf.sf <- geobr::read_state(year = 2025)

cap <- "Retas de regressão linear entre momentos não centrais de ordem $q=$ 1, 2 e 3 e histogramas de $R^2$ do ajuste representando a linearidade dos modelos ("
string <- paste0("![", cap, uf.sf$name_state[order(uf.sf$abbrev_state)], ")](", figuras, ")", "{#fig-inv-escala-", tolower(ufs), " fig-align=\"center\" width=\"15cm\"}")
string <- gsub("Relatório 2/Scripts_Invariancia_Escala/", "", string)

for(i in seq_along(string)){
  cat(string[i])
  cat("\n")
  cat("\n")
}

