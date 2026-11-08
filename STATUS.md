# Status and next steps

Last updated: 2026-11-07. The design is in [NOTES.md](NOTES.md); this file is
where the work stands, what has been measured, and what to do next.

## Where things stand

| repo | branch | remote | state |
|---|---|---|---|
| `indri` | `main` | `diogro/indri` (**private**) | **0.2.0**: three engines, uncertainty for all of them, `G` in the units of `Y` |
| `stanimal` | `master` | `diogro/stanimal` | **0.8.1** (`eigen_A`, CmdStan output fix); Laplace labels on `master` as development, unreleased |
| `limmbo` | `indri-session-fixes` | `diogro/limmbo` (fork of `meyer-lab-cshl/limmbo`) | pipeline fixes and a pixi environment; reference only |

```sh
git clone https://github.com/diogro/indri.git
git clone git@github.com:diogro/stanimal.git
git clone git@github.com:diogro/limmbo.git && cd limmbo && git switch indri-session-fixes
R CMD INSTALL stanimal/package            # >= 0.8.0 shares A's decomposition
R CMD INSTALL indri/package
Rscript -e 'cmdstanr::install_cmdstan()'  # once, if CmdStan is absent
cd indri && Rscript -e 'devtools::test("package")'
```

About 2,850 test assertions pass (asreml tests skip without a licence).

| milestone | state |
|---|---|
| **M0** subsets, assembly, projection, engine registry | done |
| **M1** stanimal engine, validated against the full model | done, gate passed |
| **M2** sommer and asreml engines | done; the two agree with each other within ~0.01 |
| **M3** shared eigendecomposition of `A` for stanimal blocks | done (stanimal 0.8.0) |
| **M4** uncertainty, weights, diagnostics | done: mixture draws for every engine, REML from each block's `vcov`, coverage measured; weights measured; convergence reported |
| **M5** performance | partly: parallel blocks, core budget, analytic projection gradient, progress bar; no caching |

## Current defaults

Blocks of up to four traits (`s = min(floor(sqrt(K)), 3)`); each pair in 3
distinct blocks for stanimal and at least 1 for the REML engines; stanimal
blocks with `lkj_prior = 1`; `spread = "mixture"`; equal weights; least-squares
projection; asreml on 2 workers. NOTES.md §9 lists the evidence for each.

## What has been measured

All studies simulate a known G with `stanimal::random_covariance(method =
"partial")` on a sibship design (families of five, N = 300 unless stated) and
score `corrG` against the truth. The study scripts named below were removed
from `examples/` once their decisions were made; they are in git history at
`7b67eab` (`git show 7b67eab:examples/uncertainty_coverage.R`), and their
tracked outputs are in `results/`.

### Accuracy against the full model (M1 gate)

12 traits, s = 4, three seeds (`validation_12_traits.R`): indri's median
per-entry error was 0.0895 against 0.0964 for the full 12-trait model (ratio
0.93), better on q90 (0.213 vs 0.248) and relative Frobenius error (0.236 vs
0.283), at 4.6x the wall clock. The full model's LKJ prior shrinks harder in 12
dimensions than in a block (REVIEW.md §3), so part of that edge is the prior.

### Engines compared (`compare_engines.R`, seed 1)

| K | stanimal | sommer | asreml | median abs error (stanimal / sommer / asreml) |
|---|---|---|---|---|
| 10 | 0.963 | 0.942 | 0.932 | 0.069 / 0.122 / 0.129 |
| 16 | 0.966 | 0.938 | 0.937 | 0.103 / 0.103 / 0.093 |
| 30 | 0.951 | 0.931 | 0.930 | 0.106 / 0.100 / 0.101 |

Correlation of correlations with the truth (REML at s = 3). stanimal has the
best overall agreement and the smallest q90 at every K; the REML engines match
or beat it on the median at 16 and 30 traits, and agree with each other.

### Uncertainty (`uncertainty_coverage.R`)

95% coverage of the true correlations, mean band width, slope of estimate on
truth (1 = no shrinkage) and median error:

| K, N | method | coverage | width | slope | median err |
|---|---|---|---|---|---|
| 35, 300 | sommer, s = 3, c = 1 | 0.964 | 0.719 | 0.933 | 0.102 |
| 35, 300 | **stanimal, s = 3, LKJ 1 (default)** | **0.959** | 0.574 | 0.848 | **0.096** |
| 35, 300 | stanimal, s = 3, LKJ 2 | 0.937 | 0.554 | 0.808 | 0.104 |
| 35, 300 | stanimal, s = 4, LKJ 1 | 0.935 | 0.549 | 0.809 | 0.101 |
| 35, 300 | stanimal, s = 4, LKJ 2 (old default) | 0.899 | 0.533 | 0.774 | 0.109 |
| 30, 1000 | sommer / asreml, s = 3, c = 1 | 0.969 / 0.971 | 0.34 | | 0.049 |
| 12, 300 | sommer / asreml, s = 3, c = 3 | 0.970 / 0.977 | 0.69 / 0.74 | | 0.103 / 0.106 |
| 12, 300 | sommer / asreml, pairs | 0.957 / 0.960 | 0.69 | | 0.107 |
| 12, 300 | stanimal, s = 3, LKJ 1 | 0.952 | 0.55 | 0.879 | 0.092 |
| 12, 300 | stanimal, s = 4, LKJ 2, mixture | 0.894 | 0.52 | 0.795 | 0.098 |
| 12, 300 | stanimal, s = 4, LKJ 2, *averaged draws* | 0.404 | 0.16 | | |
| 12, 300 | full stanimal model | 0.801 | 0.49 | | |

- **Averaging block draws is overconfident** by about `sqrt(coverage)`
  (blocks share the data, their chains do not); the mixture fixes it.
- **The REML bands are slightly conservative** (0.96-0.98): the estimate
  averages several contexts and is less noisy than one block, while the band
  carries one block's sampling variance.
- **stanimal's undercoverage was shrinkage** (every miss had the truth further
  from zero than the band); 4-trait blocks and LKJ 1 bring it to nominal.
- c = 3 and c = 1 give the same REML coverage, width and accuracy; a
  recentred mixture was no better (REML) or worse (stanimal).
- The full 35-trait stanimal model was dropped as a reference: two of its six
  chains were still in warmup after 71 minutes.

### Integration (K = 35)

Eigenvalue variance (`evolqg::CalcEigenVar()`) of `corrG` over that of the true
correlation matrix (0.23): sommer 0.83, stanimal (default) 0.74, stanimal
s = 4, LKJ 2 0.64. Every estimate understates integration, in step with the
shrinkage of the strongest correlations; noise alone would overstate it.

Where it is lost, on the saved fits (slope at blocks / average / projection):
sommer 0.972 / 0.973 / 0.933, stanimal 0.872 / 0.859 / 0.848. stanimal's loss is
in the blocks (the prior). For REML the projection moves a strongly indefinite
average (smallest eigenvalue -0.52) onto the PSD cone. Alternatives to the
projection, scored on the same fits:

| sommer, K = 35 | slope | integration | median err | err, strongest 10% |
|---|---|---|---|---|
| least squares (current) | 0.933 | 0.83 | 0.102 | 0.053 |
| least squares, inverse-variance | 0.926 | 0.83 | 0.102 | 0.056 |
| nearest correlation (Higham) | 0.944 | 0.80 | 0.102 | 0.048 |
| eigenvalue clipping | 0.890 | 0.90 | 0.103 | 0.079 |

None is better on every measure (clipping raises integration by distorting the
strongest correlations), and the projection inflates the variances by under
2%. **The current projection stays.** Understating integration is the
preferred direction of error.

### REML block size (`block_size_sweep.R`)

Paired within seed against pairs (c = 1): sommer s = 3 cut the median error by
0.007 at K = 16 (5 of 6 seeds) at the same cost, s = 4 and 5 a little more at
2-5x the time; at K = 30 the differences were 0.001-0.003. asreml s = 3:
-0.011 at K = 16 (2 seeds), level at K = 30, 0.6x the time. Larger blocks
certify convergence less often (asreml K = 30: 99% / 80% / 28% at s = 2 / 3 / 4)
without losing accuracy.

### Non-converged REML pairs (`reml_convergence.R`)

2.4% of sommer pairs and 3.8% of asreml pairs report no convergence. They are
as accurate as the rest (asreml median error 0.099 vs 0.107; sommer 0.115 vs
0.098, same q90), sommer reaches the same optimum from every start (the
log-likelihood moves by ~1e-6 when it stops), and no refit (more iterations,
restart, h2 = 0.5 or univariate starts, a 3-trait block) improves accuracy; the
larger block makes it worse. They are uncertified optima of hard pairs (strong
true correlations, boundaries), and are reported in `fit$diagnostics`.

### Weights

12 traits, seed 1: inverse-variance (median error 0.0894, q90 0.209) vs equal
weights (0.0866, 0.213), a wash; equal weights stay the default.

### Shared eigendecomposition of A (`eigen_A_timing.R`)

One six-trait block, four chains, reference BLAS: 107 s per block at N = 4000
with stanimal decomposing `A` in each chain, 50 s with indri's shared
decomposition; 7 s saved at N = 2000, nothing below 1000. The models and their
draws are unchanged.

## Changes and fixes

- **indri 0.2.0 (breaking):** `G`, `G_bar`, `G_draws` and `sigma_G` are returned
  in the units of `Y`; estimation stays on the standardised scale and nothing
  is re-projected. Scripts that multiplied `fit$G` by `diag(y_sd)` must stop.
  Fits saved by 0.1.0 (at or before `4e8a719`) hold a standardised `G`.
- **Subset coverage** (`0837352`): `indri_subsets()` did not guarantee
  `min_cooccurrence` (inherited from LiMMBo). It now tops up, and counts and
  keeps only distinct subsets.
- **sommer trait order** (`4b5748a`): blocks came back in alphabetical trait
  order, so variances (pairs) or covariances (larger blocks) landed on the
  wrong traits in any run with 10 or more traits. Every sommer number above
  postdates the fix; earlier ones are superseded.
- **CmdStan output collision**: fits from the same random state within a minute
  shared an output file; each stanimal block now has its own directory
  (indri), and stanimal 0.8.1 gives every run a unique basename.
- **asreml licence**: two processes on one machine at a time. indri defaults to
  two workers, drops to one when asreml is already loaded, and reports a seat
  failure clearly. A local asreml run locks out a remote one.
- **Convergence reporting**: `diagnostics$converged` is `NA`, not `TRUE`, when no
  block can tell; unconverged blocks are listed by trait.
- **asreml with covariates**: asreml's own budget (13 iterations, one update)
  stopped far from the optimum once `X` held covariates (G off by 0.12 where
  sommer was right; centring the covariate did not help). The engine now runs
  `maxit = 50` and up to 5 updates while unconverged, which matched sommer
  exactly at no cost when a fit converges early. The convergence study above
  had an intercept only, where non-convergence was benign.
- **Gaps in Y with the REML engines**: they were flagged as unable to handle
  gaps, so `fit_G()` refused them or dropped incomplete individuals. Both
  handle gaps natively (sommer's long format skips the missing record, asreml
  keeps partly observed rows) and agree exactly; `missing = "model"` is now
  the default for every engine. Keeping the gaps cost 7% in standard error at
  15% missing, against 29% for dropping incomplete individuals.
- **Design checks**: `fit_G()` refuses an `X` that cannot represent a constant,
  is rank deficient, has no columns or has missing values, as stanimal does;
  without a constant sommer and asreml fitted different models.

## Next

1. **Missing data across the assembly.** Every engine now handles gaps within a
   block (the REML engines agree exactly on gapped data); still unmeasured:
   coverage and accuracy of the assembled G when traits have unequal gaps, and
   stanimal's imputation against REML on the same gapped data.
2. **Real data: ratones (in progress, 2026-11-07).** 35 skull traits, 329
   mice, AGE + SEX + LIN (`ratones_G.R`). sommer and stanimal (LKJ 1, 1000
   iterations) agree on the structure (`corrG` correlation 0.96), with the
   correlations stronger than P's (mean |r| 0.385 REML, 0.300 stanimal, 0.238
   P) and G more integrated than P. But the data constrain G poorly: `corrG`
   has 17 (REML) and 9 (stanimal) zero eigenvalues of 35, 30% of sommer
   blocks are uncertified and 213 of 599 stanimal blocks had Rhat >= 1.05 at
   1000 iterations. **This dataset seems to benefit from LKJ 2's
   regularisation**; revisit the stanimal default (or a data-size rule for
   it) with this in mind. **Longer chains do not fix it**: at 3000 iterations
   215 of 599 blocks still have Rhat >= 1.05 (max 1.28), with the same 9 zero
   eigenvalues, so those block posteriors are poorly identified under LKJ 1,
   not under-sampled. Results in `results/ratones_*` (tags `_sommer`,
   `_stanimal_iter3000`; untagged is LKJ 1 at 1000 iterations).
3. **Large-N overhead.** With the decomposition shared, about 12 s per block
   remains at N = 4000, probably CmdStan reading the N x N `Q` as JSON per
   block. Measure before optimising.
4. **stanimal review**, set aside for later: its own `lkj_prior` default, and
   releasing the Laplace labels.

## Running the studies

The remaining examples (`six_traits.R`, `ratones_G.R`, `compare_engines.R`)
document their settings in their headers; see `examples/README.md`. To rerun a
removed study, restore it from `7b67eab` first:

```sh
git show 7b67eab:examples/uncertainty_coverage.R > examples/uncertainty_coverage.R
INDRI_COVERAGE_SEEDS=1 INDRI_COVERAGE_CORES=48 \
  nohup Rscript examples/uncertainty_coverage.R > results/cov_seed1.log 2>&1 &
```

`*_CORES` is a process budget for `fit_G()`, not a worker count; run one seed
per process. `results/validation_12_traits_*.csv` and
`results/uncertainty_coverage_*.csv` are tracked, everything else under
`results/` is ignored. AGENTS.md has the practical traps (connection limits,
ports, licence seats, process-name matching).
