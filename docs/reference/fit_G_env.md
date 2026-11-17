# Environment-specific genetic covariance (experimental)

Estimates the genetic covariance of a set of traits within and between
the groups of a discrete environment (sex, diet, site), treating each
trait in each environment as its own trait ("character states"). Each
individual is observed only in its own environment, so the genetic
covariance between environments is identified through relatives in
different environments, and the residual covariance between environments
is not identified at all: it is zero by construction.

## Usage

``` r
fit_G_env(Y, env, X = NULL, A, engine = "sommer", ...)
```

## Arguments

- Y:

  Numeric matrix, individuals by traits.

- env:

  The environment of each individual: a vector or factor with one value
  per row of `Y` and at least two levels, or the name of a column of
  `X`.

- X, A:

  As for [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md).

- engine:

  Block engine; `"sommer"` by default (the stanimal engine imputes every
  unobserved character state and is much slower here). `"wombat"` holds
  the unidentified residual covariances between environments at zero
  instead of estimating them.

- ...:

  Passed to
  [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md).

## Value

An object of class `indri_fit_env`: `G` and `E` (character states by
character states, in the units of `Y`), `corr_between` (traits by pairs
of environments: each trait's genetic correlation between two
environments), `levels`, `traits`, `states` (the character-state names),
`env` (each individual's environment, named by individual) and `fit`
(the [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md) fit
of the character states). Pass it to
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
as `fit` to scan with this covariance.

## Details

The response is reshaped to one column per trait and environment
(`"trait:level"`), each observed only in that environment's individuals,
and [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md) fits
it with its usual assembly and engines. Unless `subsets`, `s` or
`min_cooccurrence` is given, the subsets are trait-paired: one block per
pair of traits, holding both traits in every environment
(`P (P - 1) / 2` blocks of `2 L` states for `P` traits in `L`
environments;
[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md)'s usual
subsets for fewer than three traits). Each cross-environment covariance
is then estimated in every block that holds the trait, next to its own
within-environment variances. `G` is the assembled genetic covariance of
the character states: its diagonal blocks are the within-environment `G`
and its off-diagonal blocks the genetic covariance between environments
(a cross-environment genetic correlation below 1 is
genotype-by-environment interaction in the background). `E` is block
diagonal: within each environment, the residual covariance of the traits
after `X` minus that environment's block of `G` (scaled for the
relatedness within it), floored like
[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md)'s `E`.

The environment's main effect is carried by the character states' own
means, so a column of `X` named by `env` is dropped; every other column
of `X` must vary within each environment.

This is a trial: the cross-environment covariances rest on relatives in
different environments only. In simulations (900 sibs of both sexes, 10
and 30 traits; STATUS.md)
[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md)'s usual
subsets gave cross-environment genetic correlations **too low**
(0.67-0.77 where the truth was 1, 0.38-0.45 where it was 0.5): most
blocks held a trait in one environment only, and the projection onto a
positive semi-definite `G` pulled the few cross-environment entries
towards zero. Trait-paired subsets give 0.84-0.85 and 0.49-0.50. Scans
with this covariance
([`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
with this fit) remain somewhat conservative for interaction tests
(genomic inflation 0.91-0.99), not anti-conservative.

## See also

[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md),
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
