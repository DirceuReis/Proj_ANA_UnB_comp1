# Organização dos scripts de execução

Os códigos executados neste projeto contam com uma versão em R acompanhada de uma versão escrita em QMD — utilizada para gerar as páginas no site — que apresentam um nível de detalhe maior sobre o funcionamento do código.

Todos os dados/resultados produzidos (pasta base/gerados/) neste projeto foram gerados a partir dos scripts em R (pasta scripts/) que, por sua vez, utilizam funções localizadas em scripts/funcoes/.

## Investigações preliminares e caracterizações

-   **script_leitura_subdiario.R**: combina os dados em formato ".pqt" consolidados, realiza correções de fuso-horário na coluna `datetime` e exporta os arquivos **df_subdaily_data.parquet** e **df_subdaily_info.parquet**;
-   **script_agg_imax_subdiario.R**: utiliza as funções `fun_group_by_ts()` e `fun_imax_agg()` ou `imax.wateryear()` que é a versão mais recente --- que internamente utiliza a função `fun_filter_set()` que completa das datas, incluindo `NA`s quando não houve registro --- para agregar as lâminas de precipitação no arquivo **df_subdaily_data.pqt** em diferentes durações e extrair os máximos anuais para cada uma, exportando o arquivo **df_imax.pqt** com os resultados do processamento;
-   **script_descricao_subdiario.R**: realiza a descrição dos dados subdiários, análise de falha. Atualmente falta incluir análise bootstrap das razões de L-momentos (L-assimetria e L-curtose) e a obtenção de estatísticas básicas do conjunto para diferentes durações;
-   **script_duracao_subdiario.R**: análise empírica da hipótese de invariância de escala. Estima o expoente de escala da hipótese em seu sentido amplo (*wide-sense*) a partir dos momentos não centrais e gera figuras para avaliação utilizando a função `fun_scale_invariance()`;
-   **script_scale_block_bootstrap.R**: expande a análise feita no script anterior calculando intervalos de confiança nas estimativas do expoente de escala ($H$) e de posição ($\theta$), além dos coeficientes de desagregação, utilizando bootstrap por blocos. O parâmetro $\theta$ pode ou não ser incluído na estimativa através do argumento `offset = FALSE/TRUE` gerando dois resultados: **scale_invariance/df_scale_block_boot.pqt** e **scale_invariance/df_scale_block_boot_offset.pqt** respectivamente;
-   **script_coef_desagregacao.R**: a partir dos resultados em **df_scale_block_boot.pqt** e **df_scale_block_boot_offset.pqt**, avalia os coeficientes de desagregação de lâminas de precipitação, comparando-os com coeficientes convencionais da CETESB.

### Ordem de execução

1.  **script_leitura_subdiario.R**: lê arquivos "RS\_\<ano\>\_info/data.pqt" na pasta "base/fonte/consolidado/subdiario" e gera um único arquivo compilando todas as estações e todos os anos;
2.  **script_agg_imax_subdiario.R**: cria a tabela de intensidades máximas anuais de diferentes durações "df_imax.pqt" que é utilizado nas análises subsequentes.

## Construção de modelos de Intensidade-Duração-Frequência

-   Bayesiano, dGEV, dGEV + GMLE...

# Git

## Resolução de conflitos

### Problema: `git pull` gera erro de `merge`

**Resolução**: paciência e calma, leia primeiro a mensagem de erro.

``` bash
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

Há diferenças entre os arquivos em `origin/main` e `main` (local)? Veja quais são as diferenças:

``` bash
git diff main origin/main # mostra quais são as diferenças no console (exceto arquivos binários)
```

### Problema: o branch atual está vários commits adiantados do `origin/main`

**Resolução**: reestabeleça os commits anteriores para a versão do commit que não tinha problemas. Uma das fontes principais desse problema é quando é adicionado ao `commit` um arquivo com mais de 100 MB por engano. Este arquivo permanece no `commit`, mas trava o `git push`. O problema é que o `commit` permanece e agora está adiantado com relação ao `origin/main`. A solução mais prática que encontrei até o momento é voltar atrás os `commits` até que esteja na mesma versão do repositório remoto. **Importante**: isso não significa que, ao voltar atrás, as atualizações feitas no repositório local serão desfeitas, elas são mantidas, o que é desfeito são os `commits`.

Voltamos atrás usando `git reset:`

``` bash
git reset ---soft HEAD~n # em que 'n' é quantas vezes queremos voltar atrás
```

Se as mensagens de aviso disserem "Repositório atual 3 commits adiantado do principal", então usamos `HEAD~3`.
