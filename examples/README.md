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
| `gwas.R` | a genome scan on simulated sibs: `fit_G()` (sommer) for G and E, then `indri_gwas()` multi-trait, single-trait and with dominance, on four planted SNPs; calibration, permutations, and Manhattan/QQ plots from `plot()`. Under a minute | `INDRI_GWAS_CORES`, `_SEED` |
| `gwas_calibration.R` | null scans of simulated sibs: genomic inflation and type-I error of every test, with indri's G and E and the truth | `INDRI_GWAS_STUDY_P`, `_SEEDS`, `_N_FAM`, `_M`, `_CORES` |
| `gwas_power.R` | power of the multi-trait test against single-trait tests, by the pattern of a SNP's effects across traits, and of the dominance test | `INDRI_GWAS_STUDY_P`, `_SEEDS`, `_N_FAM`, `_REPS`, `_H`, `_CORES` |
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

## Reading `gwas.R`

- **`delta` near 1** says the `G` and `E` from `fit_G()` describe the data
  being scanned. The scan re-fits this one ratio under the null, so a `G` and
  `E` off in scale would still give a calibrated test, but a value far from 1
  is worth a look.
- **Genomic inflation near 1** on the null SNPs: the tests are calibrated
  with relatives in the sample, because `A` is in the model.
- **Each test has its SNP.** The joint test is strongest for small effects on
  many traits, especially against their correlation (`rs_pleio`); the
  single-trait tests for a large effect on one trait (`rs_one`, where the
  joint test spends 10 degrees of freedom); the dominance test for an effect
  the additive coding cannot see (`rs_dom`).
- **Dominance needs all three genotypes in numbers.** At low allele frequency
  the minor homozygote is rare, the dominance and additive codings are nearly
  collinear, and both effects become noisy.

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
