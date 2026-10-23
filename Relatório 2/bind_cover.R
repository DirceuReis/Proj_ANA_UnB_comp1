if(!require(pak)) install.packages("pak")
pak::pkg_install("officer")
library(officer)

cover <- read_docx("Relatório 2/.capa.docx")
body <- read_docx("Relatório 2/output_docx/Relatório-2.docx")
final <- body_add_docx(cover, src = "Relatório 2/output_docx/Relatório-2.docx")
print(final, target = "Relatório 2/output_docx/Relatório 2 (com capa).docx")
