As imagens nesta pasta precisam ser geradas novamente, pois elas usam estações que não foram filtradas propriamente. O filtro utilizado anteriormente selecionava estações em um vetor `gauges` da seguinte maneira:

1.  Filtro em cada ano de cada estação para selecionar somente os anos em que `na_prct <= na.accept`;
2.  Filtro em cada estação, contando o número de anos após o filtro em (1) e mantendo somente as estações em que `n_years >= min.years`.

O problema era na hora de filtrar as estações, isso porque o filtro que estava sendo utilizado era:

``` {r}
df.imax[df.imax$gauge_code %in% gauges,]
```

Então mesmo se as estações naturalmente apresentassem anos com falhas elevadas acima de `na.accept`, estes anos não seriam removidos.

As estações na pasta `8 anos 20 prct` já apresentam o filtro corrigido, agora adicionando:

```{r}
df.imax[df.imax$na_prct <= na.accept & df.imax$gauge_code %in% gauges,]
```
