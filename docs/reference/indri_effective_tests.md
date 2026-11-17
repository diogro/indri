# Effective number of independent traits

Galwey's (2009) effective number of tests, from the eigenvalues `lambda`
of the Spearman correlation matrix of the traits:
`(sum(sqrt(lambda)))^2 / sum(lambda)`. It is 1 for perfectly correlated
traits and `ncol(Y)` for independent ones, and is the factor
`adjust = "effective"` multiplies single-trait p-values by in
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md).

## Usage

``` r
indri_effective_tests(Y, X = NULL)
```

## Arguments

- Y:

  Numeric matrix, individuals by traits.

- X:

  Optional fixed-effects design, as in
  [`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md).

## Value

A number between 1 and `ncol(Y)`.

## Details

With a design `X`, the correlation is that of the residuals of `Y` on
`X`: the single-trait tests see the traits after the fixed effects, and
covariates shared by the traits (sex, age, line) would otherwise make
them look more correlated than they are to the tests.

## References

Galwey, N. W. (2009). A new measure of the effective number of tests, a
practical tool for comparing families of non-independent significance
tests. Genetic Epidemiology 33, 559-568.
[doi:10.1002/gepi.20408](https://doi.org/10.1002/gepi.20408)
