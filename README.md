# Organização dos scripts de execução:
Os códigos executados neste projeto contam com uma versão em R acompanhada de uma versão escrita em QMD -- utilizada para gerar as páginas no site -- que apresentam um nível de detalhe maior sobre o funcionamento do código.

Todos os dados/resultados produzidos (pasta base/gerados/) neste projeto foram gerados a partir dos scripts em R (pasta scripts/) que, por sua vez, utilizam funções localizadas em scripts/funcoes/.

## Investigações preliminares e caracterizações
-   **script_leitura_subdiario.R**: combina os dados em formato ".parquet" consolidados, realiza correções de fuso-horário na coluna `datetime` e exporta os arquivos **df_subdaily_data.parquet** e **df_subdaily_info.parquet**;
-   **script_agg_imax_subdiario.R**: utiliza as funções `fun_group_by_timestep()` e `fun_imax_agg()` --- que internamente utiliza a função `fun_filter_set()` que completa das datas, incluindo `NA`s quando não houve registro --- para agregar as lâminas de precipitação no arquivo **df_subdaily_data.parquet** em diferentes durações e extrair os máximos anuais para cada uma, exportando o arquivo **df_imax.parquet** com os resultados do processamento;
-   **script_descricao_subdiario.R**: realiza a descrição dos dados subdiários, análise de falha. Atualmente falta incluir análise bootstrap das razões de L-momentos (L-assimetria e L-curtose) e a obtenção de estatísticas básicas do conjunto para diferentes durações;
-   **script_duracao_subdiario.R**: análise empírica da hipótese de invariância de escala. Estima o expoente de escala da hipótese em seu sentido amplo (*wide-sense*) a partir dos momentos não centrais e gera figuras para avaliação utilizando a função `fun_scale_invariance()`;
-   **script_scale_block_bootstrap.R**: expande a análise feita no script anterior calculando intervalos de confiança nas estimativas do expoente de escala ($H$) e de posição ($\theta$), além dos coeficientes de desagregação, utilizando bootstrap por blocos. O parâmetro $\theta$ pode ou não ser incluído na estimativa através do argumento `offset = FALSE/TRUE` gerando dois resultados: **scale_invariance/df_scale_block_boot.pqt** e **scale_invariance/df_scale_block_boot_offset.pqt** respectivamente;
-   **script_coef_desagregacao.R**: a partir dos resultados em **df_scale_block_boot.pqt** e **df_scale_block_boot_offset.pqt**, avalia os coeficientes de desagregação de lâminas de precipitação, comparando-os com coeficientes convencionais da CETESB.

## Construção de modelos de Intensidade-Duração-Frequência
-   Bayesiano, dGEV, dGEV + GMLE...

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
