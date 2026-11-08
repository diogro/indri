# indri

Estimates the **additive genetic covariance matrix G** (traits × traits) for
large trait sets, without fitting one giant multivariate animal model.

It does this LiMMBo's way: split the traits into many small, overlapping
subsets, fit a small animal model on each, and reassemble the estimates into a
single positive semi-definite matrix. The per-subset fit is pluggable, so the
same estimator runs on a Bayesian sampler or on a REML engine, and every engine
comes with per-entry uncertainty. See the [LiMMBo paper](#reference).

> **Status (0.2.0):** working and validated on simulations of up to 35 traits.
> The assembled `G` is as accurate as the full model where the full model can
> be fitted; the 95% bands cover the true correlations about 95-97% of the
> time. See [STATUS.md](STATUS.md) for the numbers.

## Install

```r
remotes::install_github("diogro/indri")
```

You also need at least one block engine:

| Engine | Package | Needs | Uncertainty |
|---|---|---|---|
| `"stanimal"` (default) | [`stanimal`](https://github.com/diogro/stanimal) (>= 0.8.0 recommended) | CmdStan | posterior draws |
| `"sommer"` | [`sommer`](https://cran.r-project.org/package=sommer) | | REML sampling distribution |
| `"asreml"` | [ASReml-R](https://vsni.co.uk/software/asreml-r/) | a licence (one seat per worker) | REML sampling distribution |

Every engine accepts gaps (`NA`) in `Y` and uses each individual's observed
traits. stanimal gives a posterior; the REML engines are much faster for large
N.

## Quick start

```r
library(indri)
set.seed(1)

# A known 6-trait genetic covariance on a sibship design.
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

fit <- fit_G(
  Y, A = A, engine = "stanimal",
  engine_args = list(chains = 1, iter = 600, refresh = 0, show_messages = FALSE)
)

fit                                            # summary print
fit$G                                          # G, in the units of Y
fit$corrG                                      # genetic correlations
apply(fit$corrG_draws, 1:2, quantile, 0.975)   # per-entry upper 95% bound
```

Fixed effects go in a design matrix with one row per individual, which must
include an intercept (or every level of a factor). stanimal's
`animal_design()` builds all three inputs from a formula, a data frame and a
pedigree:

```r
d <- stanimal::animal_design(cbind(t1, t2, t3) ~ sex + age, data = dat,
                             pedigree = ped)
fit <- fit_G(d$Y, X = d$X, A = d$A)
```

## What you get back

| slot | what it is |
|---|---|
| `G`, `corrG`, `sigma_G` | the assembled estimate (`G` and `sigma_G` in the units of `Y`) |
| `G_draws`, `corrG_draws` | assembled draws (traits × traits × draws), for per-entry intervals |
| `G_bar` | the average of the blocks before the PSD projection |
| `subsets`, `coverage` | the trait subsets fitted, and how many covered each pair |
| `diagnostics` | per-block convergence; uncertified blocks listed by trait |
| `blocks` | the raw block fits, on the internal standardised scale |
| `y_mean`, `y_sd` | that standardisation |

`print()`, `summary()` and `plot()` are provided. A progress bar of the block
fits shows in interactive sessions (`progress = TRUE` forces it).

## Caveats worth knowing

- **Intervals are per entry.** The draws are assembled from dependent blocks:
  use them for intervals on individual covariances and correlations, not for
  functions of the whole matrix (eigenvalues, integration). Measured coverage
  at nominal 95%: 0.96-0.98 for the REML engines, about 0.95 for stanimal.
- **Strong correlations are shrunk a little**, so integration is understated
  rather than overstated: the eigenvalue variance of `corrG` came out at about
  75-85% of the truth's at 35 traits.
- **Blocks are small by default**: up to four traits each, every pair in three
  distinct blocks for stanimal (fitted with `lkj_prior = 1`) and in at least one
  for the REML engines. Larger blocks shrank stanimal's correlations and gave
  the REML engines nothing for their cost.
- **asreml needs a licence seat per worker.** It runs two workers by default; a
  standard licence allows two processes on one machine at a time. Use
  `cores = 1` with a single seat.
- **Units.** `G` is in the units of `Y` since 0.2.0; fits saved by 0.1.0 hold a
  standardised `G`.

## Repository map

- [`NOTES.md`](NOTES.md): design, engine contract and the reasons for the
  defaults.
- [`STATUS.md`](STATUS.md): what has been measured and what comes next.
- [`AGENTS.md`](AGENTS.md): guidance for coding agents working on the package.
- [`examples/`](examples/README.md): worked examples and the validation studies.

## Scope

In scope: `G` from a phenotype matrix `Y`, a design `X` and a known
relationship matrix `A`. Out of scope for now: `E`/`P` reporting, association
testing, breeding values, imputation, model comparison.

## Reference

Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
a simple, scalable approach for linear mixed models in high-dimensional genetic
association studies. *bioRxiv*, preprint.
<https://www.biorxiv.org/content/10.1101/255497v1>
