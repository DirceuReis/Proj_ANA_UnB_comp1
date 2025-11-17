# Organização dos scripts de execução:
Os códigos executados neste projeto contam com uma versão em R acompanhada de uma versão escrita em QMD -- utilizada para gerar as páginas no site -- que apresentam um nível de detalhe maior sobre o funcionamento do código.

Todos os dados/resultados produzidos neste projeto foram gerados a partir dos scripts em R (pasta scripts/). A ordem de execução é a seguinte:
-   script_leitura_subdiario.R;
-   script_agg_imax_subdiario.R;
-   script_descricao_subdiario.R;
-   ...

# Git

## Resolução de conflitos

### Problema: 
`git pull` -> conflito no `merge` de arquivos.

**Resolução**: paciência e calma, leia primeiro a mensagem de erro.
```bash
git status                     # mostra quais são os problemas
git add <arquivo1 c/ conflito> # adicionar arquivos p/ commit
git add <arquivo2 c/ conflito> #...
# Se forem vários arquivos
git add .
# os conflitos foram resolvidos, mas o merge não foi concluído ainda (precisa do commit)
git status                     # ver como está a situação
git commit -m"mensagem"        # fazer o git commit p/ terminar de fazer o merge e resolver os conflitos
git pull                       # "Already up to date" (nem sei se precisava)
git push                       # enviar mudanças

```
