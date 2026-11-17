# Fit the assembled blocks to the closest positive semi-definite matrix

The average of overlapping subset estimates is not a covariance matrix:
it is a patchwork, and nothing constrains it to be positive
semi-definite. This fits the full matrix that is closest to all the
subset estimates at once,

## Usage

``` r
indri_project_psd(
  blocks,
  subsets,
  K = NULL,
  init = NULL,
  weights = c("none", "inverse_variance"),
  control = list()
)
```

## Arguments

- blocks:

  List of subset fits, each with a `G` matrix and, when
  `weights = "inverse_variance"`, a `se` matrix or `draws` array.

- subsets:

  Integer matrix of trait indices, one subset per row.

- K:

  Number of traits; inferred from `subsets` when `NULL`.

- init:

  Starting matrix. Defaults to the average of the blocks.

- weights:

  `"none"` for LiMMBo's equal weights, or `"inverse_variance"` to weight
  each block's entries by their inverse uncertainty.

- control:

  Optional list of controls for
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `matrix`, the fitted PSD matrix, plus `convergence`, `value`
and the raw `optim` result.

## Details

`C* = argmin over C >= 0 of sum_j || W_j * (C[I_j, I_j] - G_j) ||_F^2`

by writing `C = L L'` with `L` lower triangular, which makes every
candidate positive semi-definite by construction, and optimising the
`K (K + 1) / 2` free entries of `L` with L-BFGS-B. The average is used
as the starting point. The per-entry weights `W_j` are all one unless
`weights = "inverse_variance"`, in which case they come from each
block's `se` or `draws` (see
[`indri_assemble()`](https://diogro.github.io/indri/reference/indri_assemble.md)).

## Examples

``` r
# Three pairwise estimates can imply an indefinite full matrix.
subsets <- t(utils::combn(3, 2))
blocks <- list(
  list(G = matrix(c(1, 0.9, 0.9, 1), 2)),
  list(G = matrix(c(1, -0.9, -0.9, 1), 2)),
  list(G = matrix(c(1, 0.9, 0.9, 1), 2))
)
average <- indri_assemble(blocks, subsets)$G
eigen(average, symmetric = TRUE)$values
#> [1]  1.9  1.9 -0.8
fit <- indri_project_psd(blocks, subsets)
fit$matrix
#>            [,1]      [,2]       [,3]
#> [1,]  1.1599995 0.5799982 -0.5799998
#> [2,]  0.5799982 1.1599979  0.5800003
#> [3,] -0.5799998 0.5800003  1.1599990
fit$convergence
#> [1] 0
eigen(fit$matrix, symmetric = TRUE)$values
#> [1] 1.740000e+00 1.739997e+00 7.683854e-12
```
