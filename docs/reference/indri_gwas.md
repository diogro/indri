# Genome-wide association scan with the assembled covariances

Tests each SNP for association with the traits in `Y`, either jointly
(one multi-trait "any effect" test per SNP) or trait by trait, in a
linear mixed model with relationship matrix `A` or in a linear model.
The multi-trait mixed model uses the genetic and residual covariances
`G` and `E`, usually from
[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md), which is
what makes it feasible for many traits.

## Usage

``` r
indri_gwas(
  Y,
  genotypes,
  fit = NULL,
  G = NULL,
  E = NULL,
  A = NULL,
  X = NULL,
  snp_info = NULL,
  mode = c("multi_trait", "single_trait"),
  model = c("lmm", "lm"),
  single_trait_variance = c("reml", "fit"),
  delta = NULL,
  dominance = FALSE,
  gxe = NULL,
  adjust = c("none", "bonferroni", "effective"),
  n_permutations = 0L,
  fdr = NULL,
  seed = NULL,
  cores = 1L,
  chunk_size = 1000L
)
```

## Arguments

- Y:

  Numeric matrix, individuals by traits, with no missing values.

- genotypes:

  Numeric matrix, individuals by SNPs (typically 0/1/2). Missing values
  are replaced by the SNP's mean; monomorphic SNPs get `NA`.

- fit:

  An `indri_fit` from
  [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md)
  supplying `G` and `E`, or an `indri_fit_env` from
  [`fit_G_env()`](https://diogro.github.io/indri/reference/fit_G_env.md)
  for an environment-specific background; see that section.

- G, E:

  Genetic and residual covariance matrices of the traits, in the units
  of `Y`, instead of `fit`. Matched to the columns of `Y` by name when
  both are named. For the multi-trait LM, `E` replaces the residual
  covariance.

- A:

  Relationship matrix of the individuals; required for `model = "lmm"`.

- X:

  Fixed-effects design with full column rank that can represent a
  constant, as for
  [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md);
  defaults to an intercept. Principal components go here.

- snp_info:

  Optional data frame with one row per column of `genotypes` (for
  example `snp`, `chrom`, `pos`), carried into the results.

- mode:

  `"multi_trait"` for one joint test per SNP, `"single_trait"` for one
  test per SNP and trait.

- model:

  `"lmm"` (with `A`) or `"lm"`.

- single_trait_variance:

  For the single-trait LMM, `"reml"` fits each trait's variance ratio,
  `"fit"` takes it from `G` and `E`.

- delta:

  Optional fixed variance ratio for the LMMs, instead of the fitted one:
  a number for the multi-trait test, one value per trait for the
  single-trait test.

- dominance:

  If `TRUE`, also fit a dominance (heterozygote) effect for each SNP;
  see the Dominance section. Needs genotypes coded 0, 1 and 2.

- gxe:

  Names of columns of `X` (environments, such as sex, diet or a
  covariate) whose interaction with each SNP's additive effect is also
  fitted; see the Genotype-by-environment section. A factor enters as
  the indicator columns `X` has for it.

- adjust:

  Single-trait p-value adjustment for the number of traits:
  `"bonferroni"` multiplies by `P`, `"effective"` by
  [`indri_effective_tests()`](https://diogro.github.io/indri/reference/indri_effective_tests.md)
  of the residuals after `X`; both cap at 1.

- n_permutations:

  Number of genotype permutations for empirical p-values (`0`, the
  default, for none); see the Permutations section.

- fdr:

  Optional false discovery rate for a p-value threshold estimated from
  the permutations; without `n_permutations` it runs 10 permutations
  (LiMMBo's number).

- seed:

  Optional seed for the permutations; the caller's random stream is
  restored afterwards.

- cores:

  Number of processes for the permuted scans (forked; one on Windows).

- chunk_size:

  Number of SNPs rotated at a time.

## Value

An object of class `indri_gwas`: `pvalues` (a vector over SNPs, or SNPs
by traits), `a` (the additive effects, SNPs by traits, in the units of
`Y` per allele; not to be confused with the relationship matrix `A`),
`delta`, `pvalues_adjust` and `n_tests` when adjusted, `snp_info`,
`mode`, `model` and `call`. With `dominance = TRUE` also
`pvalues_dominance`, `d` (the dominance effects: the heterozygote's
deviation from the homozygotes' midpoint) and `dominance_estimable`;
`pvalues_adjust` then adjusts the joint test. With `gxe`, also
`pvalues_gxe`, `gxe` (SNPs by traits by environments), `gxe_estimable`
and `gxe_terms`; `pvalues_adjust` then adjusts the joint test. With
permutations, `pvalues_empirical`, `qvalues_empirical` and
`n_permutations`, and with `fdr`, `fdr` and `fdr_threshold`.

## Models

For `N` individuals, `P` traits and one SNP `g`,
`Y = X B + g beta' + U + Eps`. Every test profiles an overall scale
`sigma^2` and is a likelihood-ratio test:

- multi-trait LMM:
  `vec(U + Eps) ~ N(0, sigma^2 (G (x) A + delta E (x) I))`, `P` degrees
  of freedom. `G` and `E` set the shape; the ratio `delta` is fitted
  once by REML under the null (no SNP) and held fixed for every SNP.
  With `G` and `E` fitted to the same data it should be near 1, and is
  returned so that it can be checked.

- multi-trait LM: `sigma^2 E (x) I`, `P` degrees of freedom, with `E` by
  default the covariance of the residuals of `Y` on `X`.

- single-trait LMM: per trait, `sigma^2 (A + delta I)`, one degree of
  freedom, with `delta` fitted by REML under the null
  (`single_trait_variance = "reml"`) or taken as `E[k, k] / G[k, k]`
  (`"fit"`).

- single-trait LM: ordinary least squares.

## Dominance

With `dominance = TRUE` each SNP enters as its additive dosage (0/1/2)
and a heterozygote indicator (1 for genotype 1, 0 for both homozygotes).
`pvalues` is then the joint test of both terms (2 degrees of freedom per
trait), `pvalues_dominance` the test of the dominance term with the
additive one in the null (1 per trait), and `a` and `d` the additive and
dominance effects of the joint model. A SNP observed in only two
genotype classes has a dominance term collinear with its dosage: it
keeps the additive-only test in `pvalues`, gets `NA` dominance results,
and is `FALSE` in `dominance_estimable`. Missing genotypes are imputed
by the column mean in each coding. The variance ratio is fitted without
any SNP, as without dominance.

## Genotype-by-environment

With `gxe`, each SNP's dosage `g` also enters multiplied by each named
column `e` of `X` (`g * e`), next to the environments' main effects,
which are in `X` already. `pvalues` is then the joint test of the SNP's
main effect and its interactions (and dominance), with
`P * (1 + n_gxe [+ 1])` degrees of freedom for `P` traits, and
`pvalues_gxe` tests the interactions with the main effect (and
dominance) in the null. `a` is the main effect (the SNP's effect where
every named environment is 0) and `gxe` holds the interaction effects
(SNPs by traits by environments). With `dominance = TRUE` the dominance
term is a main effect only and `pvalues_dominance` tests it given the
rest. A SNP whose interactions cannot be estimated (for instance every
carrier in one environment) keeps the test without them, gets `NA`
interaction results, and is `FALSE` in `gxe_estimable`. The covariance
is the same as without `gxe`: GxE in the background (a polygenic G x E,
or a residual variance that differs between environments) is not
modelled.

## Environment-specific background

With `fit` from
[`fit_G_env()`](https://diogro.github.io/indri/reference/fit_G_env.md),
the polygenic background has its own genetic covariance in each group of
a discrete environment and between groups, and its own residual
covariance in each: the covariance of the traits of individuals `i` and
`j` is `A_ij G[(env_i), (env_j)]` plus, for `i = j`, `delta E[(env_i)]`.
This is the background to use with `gxe` on that environment: unmodelled
polygenic GxE inflates the interaction tests. The covariance is no
Kronecker product, so the scan rotates by a generalized
eigendecomposition of its genetic and residual parts, once per connected
component of `A`: fast for families (a block-diagonal `A`), cubic in the
number of individuals times traits for a dense `A`. `delta` is fitted by
REML as usual; `single_trait_variance` does not apply.

## Permutations

With `n_permutations > 0` the genotype rows are permuted across
individuals and the scan rerun, keeping the null variance ratio of the
observed scan. `pvalues_empirical` is then
`(1 + #{permuted p <= observed p}) / (n_permutations + 1)`, per SNP (and
trait); with `dominance = TRUE` it refers to the joint test.

The permuted scans also estimate the false discovery rate, per SNP (for
single-trait scans, from each SNP's smallest adjusted p-value, or raw
without `adjust`): at a threshold `t`,
`FDR(t) = (mean # permuted p <= t per permutation) / (# observed p <= t)`.
`qvalues_empirical` is, for each SNP, the smallest `FDR(t)` over
thresholds at or above its own p-value, and with `fdr`, `fdr_threshold`
is the largest observed p-value with `FDR(t) <= fdr` (`NA` if none): the
SNPs at or below it are the discoveries at that rate. (LiMMBo's
"empirical FDR" was the `fdr`-quantile of the permuted p-values, which
is close to `fdr` itself and so an uncorrected threshold.) Permutations
are drawn in the main process (`set.seed(seed)`, then `n_permutations`
calls to `sample.int(N)`), so results do not depend on `cores`.

**Relatedness.** Permuting genotypes across individuals treats them as
exchangeable, which related individuals are not: the permuted null
breaks the link between genotype and relatedness that the mixed model
accounts for. The empirical p-values follow LiMMBo and are a check on
the analytic ones, not a replacement for them in structured samples.

The covariance is diagonalised once (eigenvectors of `A`, then of `G`
relative to `E`), so the scan costs one `N x N` rotation of the
genotypes plus `O(N P)` per SNP.

## References

Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
a simple, scalable approach for linear mixed models in high-dimensional
genetic association studies. bioRxiv, preprint.
[doi:10.1101/255497](https://doi.org/10.1101/255497)

## See also

[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md),
[`indri_effective_tests()`](https://diogro.github.io/indri/reference/indri_effective_tests.md)
