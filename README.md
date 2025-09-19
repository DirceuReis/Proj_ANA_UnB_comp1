# Git

## Resolução de conflitos

### Problema: 
`git pull` -> conflito no `merge` de arquivos.

**Resolução**: paciência e calma, leia primeiro a mensagem de erro.
```{bash}
git status                     # mostra quais são os problemas
git add <arquivo1 c/ conflito> # adicionar arquivos p/ commit
git add <arquivo2 c/ conflito> #...
# os conflitos foram resolvidos, mas o merge não foi concluído ainda (precisa do commit)
git status                     # ver como está a situação
git commit -m"mensagem"        # fazer o git commit p/ terminar de fazer o merge e resolver os conflitos
git pull                       # "Already up to date" (nem sei se precisava)
git push                       # enviar mudanças

```
