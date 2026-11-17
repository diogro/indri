# Estimate the additive genetic covariance matrix

Fits a multivariate animal model on many small trait subsets and
reassembles the subset estimates into one positive semi-definite matrix.
The estimator is LiMMBo's: the speed comes from giving each subset its
own nuisance residual covariance, so the many small fits replace one
large one.

## Usage

``` r
fit_G(
  Y,
  X = NULL,
  A,
  s = NULL,
  min_cooccurrence = NULL,
  subsets = NULL,
  engine = "stanimal",
  engine_args = list(),
  assembly = c("ls", "average"),
  weights = c("none", "inverse_variance"),
  missing = c("model", "complete_block"),
  spread = c("mixture", "average"),
  cores = NULL,
  n_draws = 1000L,
  seed = NULL,
  progress = getOption("indri.progress", interactive())
)
```

## Arguments

- Y:

  Phenotype matrix, individuals in rows and traits in columns. Gaps
  (`NA`) are allowed with every engine; see `missing`.

- X:

  Fixed-effect design, one row per individual. Defaults to an intercept.
  It must have full column rank and be able to represent a constant (an
  intercept, or the indicators of every level of a factor).
  [`stanimal::animal_design()`](https://diogro.github.io/stanimal/reference/animal_design.html)
  builds `Y`, `X` and `A` from a formula, a data frame and a pedigree:
  `d <- animal_design(cbind(t1, t2) ~ sex + age, dat, pedigree = ped); fit_G(d$Y, d$X, d$A)`.

- A:

  Relationship matrix, individuals in the same order as `Y`. When both
  have row names they are matched by name.

- s:

  Subset size parameter. Controls how many trait pairs each subset
  covers, not the number of traits in a subset; see
  [`indri_subsets()`](https://diogro.github.io/indri/reference/indri_subsets.md).
  Defaults to `min(floor(sqrt(K)), 3)`, at least 2 and kept narrow
  enough that a block never covers every trait: blocks of up to four
  traits. Larger blocks shrank the stanimal correlations more (and
  undercovered) and gave the REML engines nothing for their extra cost.

- min_cooccurrence:

  Minimum number of distinct subsets every trait pair must appear in.
  The default depends on the engine: 3 for `"stanimal"`, 1 for the REML
  engines, where three distinct contexts per pair measured no better
  than one (same accuracy and interval coverage at half the cost).
  Identical subsets are never repeated, since a pair estimated twice in
  the same trait context adds nothing, so the value is capped at what
  the block width allows (1 for two-trait blocks).

- subsets:

  Optional precomputed subset matrix, from
  [`indri_subsets()`](https://diogro.github.io/indri/reference/indri_subsets.md).

- engine:

  Block engine: `"stanimal"` (Bayesian, posterior draws, models gaps in
  `Y`), `"sommer"`, `"asreml"` or `"wombat"` (REML, with uncertainty
  from each block's sampling covariance; much faster for large N).
  `"wombat"` runs the external WOMBAT program; see
  [`indri_engine()`](https://diogro.github.io/indri/reference/indri_engine.md).

- engine_args:

  List of arguments passed to the engine's fit function.

- assembly:

  `"ls"` fits the closest positive semi-definite matrix to all the
  subset estimates; `"average"` keeps the plain average, nudged onto the
  positive semi-definite cone.

- weights:

  `"none"` for LiMMBo's equal weights, or `"inverse_variance"` to weight
  each subset estimate by its reported uncertainty (per-entry standard
  error for REML engines, per-entry posterior variance for the
  draws-based stanimal engine). Weights apply to both the average and
  the projection.

- missing:

  What to do with gaps in `Y`. `"model"` (the default) keeps every
  observed value and lets the engine handle the gaps: stanimal
  integrates them out, and the REML engines use each individual's
  observed traits (a missing value is one missing record, not a missing
  individual). `"complete_block"` drops incomplete individuals once,
  before subsetting; it discards information (with 15% of values missing
  at random, a third of the individuals and about 20% wider standard
  errors).

- spread:

  How `G_draws` and `corrG_draws` combine the block posteriors.
  `"mixture"` (the default) takes each entry of each draw from one
  covering block at random; `"average"` averages draw `d` across the
  covering blocks. Blocks are fitted to the same individuals, so
  averaging their independent chains is overconfident by roughly the
  square root of a pair's coverage (40% coverage of the true
  correlations at nominal 95%, against 95-97% for the mixture at the
  default settings). See
  [`indri_assemble_draws()`](https://diogro.github.io/indri/reference/indri_assemble_draws.md).

- cores:

  Total number of concurrent processes for the whole fit, not the number
  of subsets to run at once. Defaults to 1, except for `"asreml"`, which
  defaults to 2: every worker checks out its own ASReml licence seat,
  and a standard licence allows two at once, on one machine at a time.
  When asreml is already loaded in the session (which usually holds a
  seat) the default drops to 1 and the blocks are fitted in-session. Use
  `cores = 1` with a single-seat licence; more workers than free seats
  stop with a licence error. The budget is divided between blocks and
  chains: a block fit supplies `chains` processes by default (four for
  `"stanimal"`, its own default), so `cores = 4` with the default engine
  arguments runs one block at a time with four chains in parallel, while
  `cores = 4` with `engine_args = list(chains = 1)` runs four blocks at
  once. `chains` is never changed here; only how many run concurrently.

- n_draws:

  Number of draws per block for the REML engines, which report a joint
  sampling covariance of each block's G entries (the inverse
  average-information matrix) instead of posterior draws: indri draws
  from that normal approximation and assembles the draws like
  stanimal's. Ignored for `"stanimal"`.

- seed:

  Optional seed for the subset construction. The engine's sampler seed
  is passed through `engine_args`.

- progress:

  Show a progress bar of the block fits (blocks done, elapsed time, an
  estimate of the time left) on the console. Defaults to
  `getOption("indri.progress", interactive())`: on in an interactive
  session, off under `Rscript`. With parallel blocks the bar updates
  after each load-balanced round of about four blocks per worker. Engine
  output (for example CmdStan's own progress) is printed as before;
  silence it with
  `engine_args = list(refresh = 0, show_messages = FALSE)` for stanimal.

## Value

An object of class `indri_fit`. The estimate is `G`, in the units of
`Y`, with `corrG` its correlation matrix; `G_bar` (the average before
projection) and `sigma_G` are in the units of `Y` too. `E` is the
residual covariance, for
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md):
the covariance of the residuals of `Y` on `X` minus the assembled `G`
(scaled for the relatedness in `A`), with eigenvalues floored at 1e-3 of
the largest. `E_bar` is the average of the block residual covariances,
for reference: at many traits it is not a covariance matrix, and a
projection of it miscalibrates the association tests. `blocks` holds the
raw standardised-scale block fits. `G_draws` and `corrG_draws` (traits x
traits x draws) carry the assembled block uncertainty, from the
posterior draws of `"stanimal"` blocks or from the normal sampling
approximation of the REML blocks (see `n_draws`); use `corrG_draws` for
intervals on correlations, since it mixes each block's own correlation
draws and stays within \[-1, 1\]. Mixed draws are pointwise: their
entries are drawn independently, so a single draw need not be positive
semi-definite, and they support per-entry intervals, not functions of
the whole matrix. Neither is a posterior over the joint matrix.

## Details

Traits are standardised once, using every observed value, before the
subsets are drawn. Every engine therefore sees the same numbers and
estimates the same matrix on the same scale, whatever it does
internally; the fixed effects in `X` are handled inside each block fit.
Assembly and the positive semi-definite projection also run on that
scale, and the assembled covariances (`G`, `G_bar`, `G_draws`,
`sigma_G`) are then returned in the units of `Y`. The raw block fits in
`blocks` stay on the standardised scale: do not combine them directly
with the assembled outputs (`y_sd` converts between the two).

## References

Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
a simple, scalable approach for linear mixed models in high-dimensional
genetic association studies. bioRxiv, preprint.
[doi:10.1101/255497](https://doi.org/10.1101/255497)

## Examples

``` r
# \donttest{
# This example fits six small Stan models and needs stanimal and CmdStan.
if (requireNamespace("stanimal", quietly = TRUE) &&
    requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.na(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  K <- 4L
  N <- 60L
  # Twelve families of five full siblings.
  A <- kronecker(diag(12), matrix(0.5, 5, 5))
  diag(A) <- 1
  G_true <- matrix(0.2, K, K)
  diag(G_true) <- 0.8
  E_true <- diag(0.6, K)
  X <- cbind(intercept = 1, group = rep(c(0, 1), length.out = N))
  genetic <- t(chol(A)) %*% matrix(rnorm(N * K), N, K) %*% chol(G_true)
  residual <- matrix(rnorm(N * K), N, K) %*% chol(E_true)
  Y <- genetic + residual + X[, "group"] %o% c(0.2, -0.2, 0.1, 0)
  colnames(Y) <- paste0("trait", seq_len(K))

  # These short chains are for a quick example, not reliable inference.
  # For real work, use more chains and iterations and check diagnostics.
  fit <- fit_G(Y, X = X, A = A, subsets = t(utils::combn(K, 2)),
               engine = "stanimal", cores = 1,
               engine_args = list(chains = 2, iter = 400,
                                  refresh = 0, show_messages = FALSE))
  fit$corrG
  fit$coverage
  fit$diagnostics

  # G is in the units of Y; the standardised-scale estimate, if needed:
  G_std <- diag(1 / fit$y_sd) %*% fit$G %*% diag(1 / fit$y_sd)
  fit$sigma_G

  # Pointwise correlation spread, not a joint posterior interval.
  quantile(fit$corrG_draws[1, 2, ], c(0.025, 0.975))
}
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 209, column 4 to column 45)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/RtmpDmKWl7/model-175bd6da677d7.stan', line 207, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#>       2.5%      97.5% 
#> -0.6590653  0.9108255 
# }
```
