# Block engines

The per-subset fit is a plugin. An engine is a list with a `name`, a
`capabilities` list and a `fit` function with the signature
`fit(Y, X, a_forms, control)` returning at least `G` and `E` (on the
standardised scale, in the column order of `Y`), plus either posterior
`draws` (an `s` by `s` by `D` array; stanimal) or `se` and `vcov`, the
joint sampling covariance of the lower-triangle G entries (the REML
engines), from which
[`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md) draws.

## Usage

``` r
indri_engine(name = "stanimal")

indri_engines()
```

## Arguments

- name:

  Engine name.

## Value

An engine list.

`indri_engines()` returns the names of the implemented engines.

## Details

LiMMBo delegates its subset fits to limix; indri delegates to stanimal,
sommer, ASReml-R or the WOMBAT program, because the right engine depends
on the number of individuals, on whether the data have gaps, on whether
a licence is available, and on whether the caller wants a posterior or a
likelihood. WOMBAT is a standalone REML program (Linux and Windows, free
for academic use); the `"wombat"` engine finds it through
`engine_args = list(wombat = <path>)`, `options(indri.wombat = <path>)`
or the PATH, `engine_args = list(wombat_dir = <dir>)` keeps a block's
WOMBAT files for inspection, and `engine_args = list(wombat_rank = q)`
fits G at reduced rank q (the residual stays full rank).

## Examples

``` r
# Inspect the registry without fitting a model or loading an engine.
indri_engines()
#> [1] "stanimal" "sommer"   "asreml"   "wombat"  
engine <- indri_engine("sommer")
engine$name
#> [1] "sommer"
engine$capabilities
#> $draws
#> [1] FALSE
#> 
#> $se
#> [1] TRUE
#> 
#> $vcov
#> [1] TRUE
#> 
#> $a_form
#> [1] "inverse"
#> 
#> $missing
#> [1] TRUE
#> 
#> $default_s
#> [1] 3
#> 
#> $default_cooccurrence
#> [1] 1
#> 
```
