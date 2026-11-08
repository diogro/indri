# Examples

Runnable scripts for `indri`. Run them from the repository root (results go to
`results/`, figures to `examples/figures/`), with `indri` installed plus the
engines a script uses:

```sh
R CMD INSTALL package
Rscript examples/six_traits.R
```

Scripts with settings document them in their header as `INDRI_*` environment
variables; the main ones are below. `*_CORES` is a process budget for
`fit_G()`, not a worker count.

| script | what it does | settings |
|---|---|---|
| `six_traits.R` | the whole pipeline on six traits against a known G: coverage, per-pair correlations, variances, a plot. About three minutes | none |
| `ratones_G.R` | real data: the ratones mice shipped with stanimal (329 mice, 35 skull traits, pedigree), G against P | `INDRI_RATONES_ENGINE`, `_CORES`, `_ITER`, `_CHAINS`, `_SEED`, `_TAG` |
| `compare_engines.R` | one simulated dataset on stanimal, sommer and asreml, against the truth and each other | `INDRI_ENGINE_K`, `_CORES`, `_CORES_SOMMER`, `_CORES_ASREML`, `_S*`, `_SOMMER`, `_ASREML` |

The studies used to choose indri's defaults (validation against the full
model, uncertainty coverage, block sizes, REML convergence, eigendecomposition
timing) are no longer here; their results are in [STATUS.md](../STATUS.md)
and the scripts in git history at `7b67eab`.

## Reading `six_traits.R`

- **Pair coverage**: every trait pair is estimated in several distinct subsets.
  A pair no subset covers has no estimate, which is why `indri_subsets()` is
  built around coverage rather than random subsets.
- **The correlations** track the truth in sign and magnitude, sitting slightly
  closer to zero: the block fits carry a prior on the correlations, and at a
  few hundred individuals some shrinkage remains.
- **`min eigenvalue of G`** is positive: the average of overlapping blocks is
  not a covariance matrix, and `indri_project_psd()` makes it one.
- **The genetic variances**: `fit$G` is in the units of the traits, so its
  diagonal compares directly with the true variances.

## Making a fit faster or slower

```r
# fewer subsets, quicker, noisier: each pair in at least one subset
fit_G(Y, A = A, min_cooccurrence = 1L, ...)

# a REML engine instead of the sampler, much faster at large N
fit_G(Y, A = A, engine = "sommer", cores = 8)

# parallel blocks; for stanimal, one chain per block lets cores go to blocks
fit_G(Y, A = A, cores = 16, engine_args = list(chains = 1, refresh = 0))
```

The block fits are the cost. `cores` at the `fit_G()` level is shared between
blocks and, for stanimal, the chains inside each block.
