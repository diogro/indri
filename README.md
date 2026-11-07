# indri

Estimates the **additive genetic covariance matrix G** (traits × traits) for
large trait sets — without ever fitting one giant multivariate animal model.

It does this LiMMBo's way: split the traits into many small, overlapping
subsets, fit a small animal model on each, then reassemble the subset estimates
into a single positive semi-definite matrix. The per-subset fit is pluggable, so
the same estimator runs on a Bayesian sampler or on a REML engine. See the
[original LiMMBo paper](#reference) for the method.

> **Status:** working package, M1/M2/M4 measured. `G` itself is validated
> (within ~2× the full model's error at 12 traits, and better than it on one
> 30-trait seed); the default averaged uncertainty is overconfident, while
> `spread = "mixture"` covers ~89% at nominal 95% — see the caveats below and
> [STATUS.md](STATUS.md).

## Install

```r
# from GitHub
remotes::install_github("diogro/indri")
```

You also need at least one block engine, plus the data it requires:

| Engine | Package | Need | Returns |
|---|---|---|---|
| `"stanimal"` (default) | [`stanimal`](https://github.com/meyer-lab-cshl/stanimal) | CmdStan | posterior draws + assembled uncertainty |
| `"sommer"` | [`sommer`](https://cran.r-project.org/package=sommer) | dense `A` | REML point estimate + SEs |
| `"asreml"` | ASReml-R (licensed) | sparse `A` | REML point estimate + SEs |

The default engine is `"stanimal"`; the REML engines are much faster for
large N but drop incomplete rows and (for asreml) need a licence.

## Quick start

```r
library(indri)
set.seed(1)

# A known 6-trait genetic covariance and a small sibship design.
G_true <- stanimal::random_covariance(
  n_traits = 6, method = "partial", min_var = 0.8, max_var = 1.4
)
E_true <- stanimal::random_covariance(
  n_traits = 6, method = "lkj", lkj_prior = 2, min_var = 0.5, max_var = 0.9
)
A <- kronecker(diag(20), matrix(0.5, 8, 8))
diag(A) <- 1
Y <- t(chol(A)) %*% matrix(stats::rnorm(nrow(A) * 6), nrow(A), 6) %*% chol(G_true) +
     matrix(stats::rnorm(nrow(A) * 6), nrow(A), 6) %*% chol(E_true)
colnames(Y) <- paste0("trait", 1:6)

# Estimate G from trait subsets. `s` controls how many trait pairs each subset
# covers (the subset itself has a few more traits). It defaults to at most 4.
fit <- fit_G(
  Y, A = A, engine = "stanimal",
  engine_args = list(chains = 1, iter = 600, refresh = 0, show_messages = FALSE)
)

fit             # summary print
fit$G           # assembled G, on the standardised scale
fit$corrG       # ... as correlations
fit$y_sd        # per-trait scales used for standardisation

# back to the original units
G_original <- diag(fit$y_sd) %*% fit$G %*% diag(fit$y_sd)
```

Add fixed effects with a design matrix: `fit_G(Y, X = X, A = A, ...)`, where
`X` has one row per individual (intercept and covariates are fine).

## What you get back

`fit_G()` returns an `indri_fit` object:

| slot | what it is |
|---|---|
| `G`, `corrG` | the assembled estimate, covariance and correlation scale |
| `G_draws` | block posterior draws assembled by `spread` (`"average"` or `"mixture"`; stanimal only) |
| `coverage` | how many subsets covered each trait pair |
| `subsets` | the trait subsets that were fitted |
| `diagnostics` | per-block convergence (Rhat for stanimal, flags for REML) |
| `y_mean`, `y_sd` | the global standardisation, for back-transforming `G` |

`print()`, `summary()` and `plot()` are provided.

## Caveats worth knowing before you use it

- **Blocks are small by default.** Every engine fits overlapping blocks of up
  to four traits (`s = 3`): stanimal with each pair in three distinct blocks
  and `lkj_prior = 1`, the REML engines with each pair in at least one. Larger
  blocks were measured to shrink stanimal's correlations (and undercover) and
  to give the REML engines nothing for their cost.
- **asreml needs a licence seat per worker.** It runs two workers by default
  (a standard licence allows two processes, on one machine at a time); pass
  `cores = 1` if yours has one seat.

- **Scale.** `G` is reported on the standardised scale; multiply back with
  `diag(y_sd) %*% G %*% diag(y_sd)` for the original units.
- **Uncertainty: use `spread = "mixture"`.** Blocks are fitted to the same
  data, so averaging their draws (the current default) narrows the band by
  about the square root of a pair's coverage: ~40% coverage at nominal 95% on
  a 12-trait simulation. The mixture assembly covers ~89% (the full model:
  ~80%). Neither is a joint posterior, and mixed draws are pointwise — a mixed
  matrix need not be positive semi-definite.
- **Point estimate is the headline.** The assembled `G` itself is the
  well-validated output; at 30 traits it beat the full model on both accuracy
  and wall clock in the current measurements.

## Repository map

- [`NOTES.md`](NOTES.md) — design, engine contract, verified call forms.
- [`STATUS.md`](STATUS.md) — where the work stands and what's measured.
- [`AGENTS.md`](AGENTS.md) — guidance for coding agents working on the package.
- [`examples/`](examples/README.md) — runnable studies: the M1 validation gate,
  the cross-engine comparison, the block-size sweep, the uncertainty coverage
  study, and a real-data example (the `ratones` mice).

## Scope

In scope: `G` (traits × traits additive genetic covariance) from a phenotype
matrix `Y`, a design `X` and a known relationship matrix `A`.

Out of scope for now: `E`/`P` reporting, association testing, breeding values,
missing-data imputation, model comparison.

## Reference

Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
a simple, scalable approach for linear mixed models in high-dimensional genetic
association studies. *bioRxiv*, preprint.
<https://www.biorxiv.org/content/10.1101/255497v1>
