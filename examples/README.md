# Examples

Runnable scripts for `indri`. Run them from the repository root, because the
figures are written to `examples/figures/`:

```sh
Rscript examples/six_traits.R
```

The scripts need `indri` installed, plus `stanimal`, `cmdstanr` and CmdStan:

```sh
R CMD INSTALL package
# or, during development
Rscript -e 'devtools::install_local("package")'
```

| Script | Shows |
|---|---|
| `six_traits.R` | The whole pipeline on six traits: a known 6 x 6 genetic covariance from `stanimal::random_covariance()`, simulated on a sibship relationship matrix, estimated from 16 four-trait subsets, compared entry by entry with the truth, and plotted. About three minutes with two chains per block |

## What to look at in `six_traits.R`

- **`pair coverage`** — every one of the 15 trait pairs is estimated at least
  three times. A pair the subsets never cover has no estimate at all, which is
  why `indri_subsets()` is built around coverage rather than around drawing
  subsets at random and hoping.
- **The correlation table and the right-hand scatter panel** — the estimate
  tracks both the sign and the magnitude of the truth. At 200 individuals the
  correlations sit slightly closer to zero than the truth, because the block
  fits carry stanimal's LKJ(2) prior on the trait correlations. That shrinkage is
  the prior's, not the subsetting's: stanimal's own full six-trait model shrinks
  the same way and, on this simulation, is no closer to the truth.
- **`min eigenvalue of G`** — positive. The plain average of overlapping subset
  blocks is not a covariance matrix, and `indri_project_psd()` is what makes it
  one.
- **The genetic variances** — `fit$G` is in the units of the traits, so its
  diagonal compares directly with the true genetic variances.

## Making it faster or slower

```r
# fewer subsets, quicker, noisier: each pair covered about once
fit_G(Y, A = A, s = 3L, min_cooccurrence = 1L, ...)

# weaker prior, less shrinkage towards zero correlations
fit_G(Y, A = A, s = 3L, engine_args = list(lkj_prior = 1), ...)
```

The block fits are the cost. `cores` at the `fit_G()` level parallelises across
subsets; `cores` inside `engine_args` parallelises the chains of a single subset.

| `validation_12_traits.R` | The M1 gate: 12 traits, a known covariance, indri against the full twelve-trait model, per-entry error against pair coverage. One seed per invocation (`INDRI_STUDY_SEEDS=2`), roughly eight minutes each |
