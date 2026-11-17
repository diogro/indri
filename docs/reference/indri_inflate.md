# Place a subset covariance matrix into the full trait space

Place a subset covariance matrix into the full trait space

## Usage

``` r
indri_inflate(small, index, K, fill = NA_real_)
```

## Arguments

- small:

  Square subset covariance matrix.

- index:

  Traits the subset covers, as indices into the full set.

- K:

  Number of traits in the full set.

- fill:

  Value for the entries the subset does not cover.

## Value

A `K` by `K` matrix.

## Examples

``` r
block <- matrix(c(1, 0.4, 0.4, 2), nrow = 2)
indri_inflate(block, index = c(1, 3), K = 4)
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.0   NA  0.4   NA
#> [2,]   NA   NA   NA   NA
#> [3,]  0.4   NA  2.0   NA
#> [4,]   NA   NA   NA   NA
indri_inflate(block, index = c(1, 3), K = 4, fill = 0)
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.0    0  0.4    0
#> [2,]  0.0    0  0.0    0
#> [3,]  0.4    0  2.0    0
#> [4,]  0.0    0  0.0    0
```
