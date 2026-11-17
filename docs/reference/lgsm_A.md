# LG/SM F3 intercross: relationship matrix

The additive relationship matrix of the 1611 F3 mice of
[lgsm_f3](https://diogro.github.io/indri/reference/lgsm_f3.md), from
[lgsm_pedigree](https://diogro.github.io/indri/reference/lgsm_pedigree.md)
with the F2 parents as unrelated, non-inbred founders: 1 on the
diagonal, 1/2 between full sibs, 1/4 between half sibs (sires with more
than one mate) and 0 otherwise. The four sires missing in the source
come from the markers (`lgsm_f3$sire_inferred`).

## Usage

``` r
lgsm_A
```

## Format

A 1611 x 1611 numeric matrix with the F3 IDs (as character) as row and
column names, in the order of `lgsm_f3$ID`. It is positive definite
(smallest eigenvalue 1/2). Take the rows and columns of the mice in an
analysis, as
[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md) and
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
match by name.

## Source

Built by `data-raw/lgsm.R` from
[lgsm_pedigree](https://diogro.github.io/indri/reference/lgsm_pedigree.md).

## Examples

``` r
table(lgsm_A[upper.tri(lgsm_A)])
#> 
#>       0    0.25     0.5 
#> 1289838     752    6265 
```
