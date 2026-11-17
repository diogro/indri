# Average the subset estimates into one matrix

Every entry of the assembled matrix is a weighted mean over the subsets
that contain both of its traits. With `weights = "none"` all subsets
count equally, which is what LiMMBo does; `"inverse_variance"` weights
each estimate by the inverse of its reported uncertainty. REML engines
report a per-entry standard error; the stanimal engine reports posterior
draws, whose per-entry variance supplies the weight.

## Usage

``` r
indri_assemble(
  blocks,
  subsets,
  K = NULL,
  weights = c("none", "inverse_variance")
)
```

## Arguments

- blocks:

  List of subset fits, each with a `G` matrix and optionally a `se`
  matrix of the same shape.

- subsets:

  Integer matrix of trait indices, one subset per row.

- K:

  Number of traits; inferred from `subsets` when `NULL`.

- weights:

  `"none"` or `"inverse_variance"`. The latter needs each block to carry
  a `se` matrix or a `draws` array (used for its per-entry variance when
  `se` is absent); blocks that carry neither are treated as weight 1.

## Value

A list with `G` (the averaged matrix) and `counts` (how many subsets
contributed to each entry).

## Details

This is the averaging step, not the estimate: the average of overlapping
blocks is generally not positive semi-definite.
[`indri_project_psd()`](https://diogro.github.io/indri/reference/indri_project_psd.md)
is what turns it into a covariance matrix.

## Examples

``` r
# Two estimates of traits 1 and 2, plus blocks covering trait 3.
subsets <- rbind(c(1, 2), c(1, 2), c(1, 3), c(2, 3))
blocks <- list(
  list(G = matrix(c(1, 0.8, 0.8, 1), 2), se = matrix(0.1, 2, 2)),
  list(G = matrix(c(1, 0.2, 0.2, 1), 2), se = matrix(0.4, 2, 2)),
  list(G = matrix(c(1, 0.3, 0.3, 1), 2), se = matrix(0.2, 2, 2)),
  list(G = matrix(c(1, 0.4, 0.4, 1), 2), se = matrix(0.2, 2, 2))
)
equal <- indri_assemble(blocks, subsets, K = 3)
weighted <- indri_assemble(blocks, subsets, K = 3,
                           weights = "inverse_variance")
equal$G
#>      [,1] [,2] [,3]
#> [1,]  1.0  0.5  0.3
#> [2,]  0.5  1.0  0.4
#> [3,]  0.3  0.4  1.0
equal$counts
#>      [,1] [,2] [,3]
#> [1,]    3    2    1
#> [2,]    2    3    1
#> [3,]    1    1    2
weighted$G[1, 2]  # Closer to the more precise first block.
#> [1] 0.7647059
```
