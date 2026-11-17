# Trait subsets for the block estimator

Draws the trait subsets on which the small models are fitted. This is
the set-cover construction used by LiMMBo from version 1.0 on: each
subset is built from a window of trait pairs, and the number of subsets
is chosen so that every pair is covered at least `min_cooccurrence`
times.

## Usage

``` r
indri_subsets(K, s, min_cooccurrence = 3L, seed = NULL)
```

## Arguments

- K:

  Number of traits.

- s:

  Subset size parameter, at least 2. Controls how many trait pairs each
  subset covers.

- min_cooccurrence:

  Minimum number of distinct subsets every trait pair must appear in.

- seed:

  Optional seed for the subset construction. The caller's random state
  is left untouched.

## Value

An integer matrix with one distinct subset per row, holding 1-based
indices of the traits in that subset. Every pair is in at least
`min_cooccurrence` distinct subsets, capped at the number of distinct
subsets of this width that contain it (one, for two-trait subsets); the
coverage reached is stored in the `min_cooccurrence` attribute and the
pair counts in the `coverage` attribute.

## Details

Three things differ from LiMMBo's implementation. The permutation that
orders the pairs is seeded here, so the result is reproducible (LiMMBo
draws it from the global numpy generator and is not), and the caller's
random seed is restored afterwards. The coverage is guaranteed: LiMMBo's
passes cover each pair only `ceiling(min_cooccurrence / 2)` times by
construction and leave the rest to chance, which often falls short;
indri adds greedily built subsets of the same width until every pair
reaches `min_cooccurrence`. And coverage counts distinct subsets: a
subset is never repeated, since a pair estimated twice in the same trait
context adds nothing.

Note that the subset size is not exactly `s`. A subset covers `s - 1`
trait pairs, which is up to `2 * (s - 1)` distinct traits; the subsets
therefore come out larger than `s`, and the assembled blocks use
`ncol(subsets)` rather than `s`.

## Examples

``` r
subsets <- indri_subsets(K = 6, s = 3, min_cooccurrence = 3)
subsets
#>       [,1] [,2] [,3] [,4]
#>  [1,]    1    3    4    5
#>  [2,]    1    2    4    6
#>  [3,]    2    3    4    5
#>  [4,]    2    3    4    6
#>  [5,]    1    3    5    6
#>  [6,]    1    2    4    5
#>  [7,]    3    4    5    6
#>  [8,]    1    2    3    6
#>  [9,]    1    2    3    5
#> [10,]    1    2    5    6
#> attr(,"coverage")
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    7    5    4    3    5    4
#> [2,]    5    7    4    4    4    4
#> [3,]    4    4    7    4    5    4
#> [4,]    3    4    4    6    4    3
#> [5,]    5    4    5    4    7    3
#> [6,]    4    4    4    3    3    6
#> attr(,"min_cooccurrence")
#> [1] 3
ncol(subsets)  # Actual block width; it can be larger than s.
#> [1] 4
coverage <- attr(subsets, "coverage")
coverage
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    7    5    4    3    5    4
#> [2,]    5    7    4    4    4    4
#> [3,]    4    4    7    4    5    4
#> [4,]    3    4    4    6    4    3
#> [5,]    5    4    5    4    7    3
#> [6,]    4    4    4    3    3    6
min(coverage[upper.tri(coverage)])
#> [1] 3
```
