# Relationship matrices in the forms the engines want

The engines do not agree on how the relationship matrix arrives: a
Bayesian model takes it dense or as its eigendecomposition, sommer takes
its inverse marked as such, and ASReml takes a sparse three-column
inverse. Build the forms once per run so that a loop over many trait
subsets pays for each one only once.

## Usage

``` r
indri_a_forms(A, need = c("dense", "inverse", "sparse_inverse", "eigen"))
```

## Arguments

- A:

  Square relationship matrix.

- need:

  Which forms to build. `"dense"` is always returned.

## Value

A list with some of `dense`, `inverse`, `sparse_inverse` and `eigen`.
`eigen` is `list(Q, lambda)`, the eigenvectors and eigenvalues of `A` in
ascending order, as `stanimal::lmm_animal(eigen_A = )` takes them (and
in the order stanimal would compute them itself, so a block fitted with
it is identical to one fitted without).

## Examples

``` r
# A relationship matrix for two pairs of full siblings.
A <- kronecker(diag(2), matrix(c(1, 0.5, 0.5, 1), 2))
forms <- indri_a_forms(A, need = c("inverse", "eigen"))
forms$inverse
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  1.3333333 -0.6666667  0.0000000  0.0000000
#> [2,] -0.6666667  1.3333333  0.0000000  0.0000000
#> [3,]  0.0000000  0.0000000  1.3333333 -0.6666667
#> [4,]  0.0000000  0.0000000 -0.6666667  1.3333333
#> attr(,"inverse")
#> [1] TRUE
forms$eigen$lambda
#> [1] 0.5 0.5 1.5 1.5

# Reconstruct A from the shared eigendecomposition.
Q <- forms$eigen$Q
Q %*% diag(forms$eigen$lambda) %*% t(Q)
#>      [,1] [,2] [,3] [,4]
#> [1,]  1.0  0.5  0.0  0.0
#> [2,]  0.5  1.0  0.0  0.0
#> [3,]  0.0  0.0  1.0  0.5
#> [4,]  0.0  0.0  0.5  1.0
```
