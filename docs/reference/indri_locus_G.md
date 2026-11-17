# Genetic covariance matrices from locus effects

The additive (`G`) and dominance (`D`) covariance matrices of a set of
traits implied by the additive and dominance effects of individual loci,
their allele frequencies and the linkage disequilibrium between them
(Melo, Marroig and Wolf 2019, Eqs. 7-14). Each locus contributes through
its vector of pleiotropic effects, so `G` and `D` show how the loci
found by
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md),
or by any other mapping method, build up the genetic covariance between
the traits.

## Usage

``` r
indri_locus_G(
  a,
  d = NULL,
  genotypes = NULL,
  p = NULL,
  lambda = NULL,
  ld = TRUE,
  traits_in = c("columns", "rows")
)
```

## Arguments

- a:

  Additive effects, loci by traits (a matrix, a numeric data frame, or a
  data frame with one column of locus names; a vector for one trait), or
  an `indri_gwas` scan (whose `a`, and `d` with `dominance = TRUE`, are
  used). Loci without an estimate are dropped.

- d:

  Dominance effects, in the same layout as `a`; `NULL` for none.

- genotypes:

  Individuals by loci, 0/1/2 dosages of the allele the effects refer to,
  for the allele frequencies and the linkage disequilibrium. May have
  more loci than `a`.

- p:

  Allele frequencies, one per locus, instead of `genotypes`.

- lambda:

  Linkage disequilibrium between the loci (loci by loci; the diagonal is
  ignored), instead of the estimate from `genotypes`.

- ld:

  If `FALSE`, ignore linkage disequilibrium (only the per-locus sums).

- traits_in:

  Where the traits are in `a` and `d`: `"columns"` (the default: one row
  per locus, one column per trait) or `"rows"`.

## Value

A list with `G` and `D` (traits by traits), `alpha` (the average
effects, loci by traits), `p` and `loci` (the loci used).

## Details

For locus `j` with allele frequency `p_j` (of the allele the 0/1/2
dosage counts; `q_j = 1 - p_j`), additive effects `a_j` and dominance
effects `d_j` (vectors over the traits), the average effect is
`alpha_j = a_j + (q_j - p_j) d_j` (Eq. 7), and

- `G = sum_j 2 p_j q_j alpha_j alpha_j' + sum_{j != x} 2 lambda_jx alpha_j alpha_x'`
  (Eqs. 9-10),

- `D = sum_j (2 p_j q_j)^2 d_j d_j' + sum_{j != x} 4 lambda_jx^2 d_j d_x'`
  (Eqs. 11-12),

where `lambda_jx` is the linkage disequilibrium between loci `j` and
`x`. From an unphased genotype matrix it is the composite estimate
`lambda_jx = Cov(g_j, g_x) / 2` of the dosages.

## Conventions for effects from other software

The formulas assume one convention for the effects; estimates from other
packages have to be put on it, or `G` and `D` are wrong without any
error.

- **Effect allele.** `a` is the effect of one copy of the allele that
  the 0/1/2 genotypes count (half the difference between the two
  homozygotes), and `p` is that allele's frequency. GWAS software
  reports each beta for a named allele (for example `A1` in PLINK 2's
  `--glm` output, `ALLELE1` in regenie, `allele1` in GEMMA): the
  genotypes must count that allele, or that locus's `a` must be negated.
  Check this per locus; tools differ in which allele they report.

- **Dominance coding.** `d` is the deviation of the heterozygote from
  the midpoint of the homozygotes, which is what a regression on the
  dosage plus a heterozygote indicator (1 for genotype 1, 0 otherwise)
  estimates (as
  [`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
  does, and as the LL = +1, LS = 0, SS = -1 coding of
  [lgsm_genotypes](https://diogro.github.io/indri/reference/lgsm_genotypes.md)
  with its `D` column). Orthogonal or frequency-centred dominance
  codings (Cockerham's, for instance) estimate other quantities and must
  be converted first. Most GWAS software estimates no dominance: leave
  `d = NULL`.

- **Units.** Effects in the units of the traits give `G` and `D` in
  those units squared. Standardised betas give matrices on the
  standardised scale.

- **Shape.** Effects are loci by traits (one row per locus); a matrix
  with one row per trait needs `traits_in = "rows"`. A data frame may
  carry one non-numeric column of locus names.

Loci are matched to the columns of `genotypes`, and to `d`, `p` and
`lambda`, by name whenever both sides are named. Matching to `genotypes`
by position warns, as does matching `d` by position when only one side
is named.

## References

Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic perspective on
multivariate variation, pleiotropy, and evolution. Journal of Heredity
110, 479-493.
[doi:10.1093/jhered/esz011](https://doi.org/10.1093/jhered/esz011)

## See also

[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)

## Examples

``` r
# one locus with effects on three traits
indri_locus_G(matrix(c(0.5, -0.2, 0.1), 1), d = matrix(c(0.3, 0, 0), 1),
              p = 0.4)$G
#>           [,1]     [,2]     [,3]
#> [1,]  0.150528 -0.05376  0.02688
#> [2,] -0.053760  0.01920 -0.00960
#> [3,]  0.026880 -0.00960  0.00480
```
