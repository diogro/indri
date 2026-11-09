# Status and next steps

Last updated: 2026-11-08. The design is in [NOTES.md](NOTES.md); this file is
where the work stands, what has been measured, and what to do next.

## Where things stand

| repo | branch | remote | state |
|---|---|---|---|
| `indri` | `main` | `diogro/indri` (**private**) | **0.2.0**: three engines, uncertainty for all of them, `G` in the units of `Y` |
| `indri` | `gwas` | `diogro/indri` | development: `indri_gwas()` (multi- and single-trait, dominance, permutations, plots), `E` in `fit_G()`, the LG/SM F3 data; not merged |
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

About 3,300 test assertions pass on `gwas` (asreml tests skip without a
licence).

| milestone | state |
|---|---|
| **M0** subsets, assembly, projection, engine registry | done |
| **M1** stanimal engine, validated against the full model | done, gate passed |
| **M2** sommer and asreml engines | done; the two agree with each other within ~0.01 |
| **M3** shared eigendecomposition of `A` for stanimal blocks | done (stanimal 0.8.0) |
| **M4** uncertainty, weights, diagnostics | done: mixture draws for every engine, REML from each block's `vcov`, coverage measured; weights measured; convergence reported |
| **M5** performance | partly: parallel blocks, core budget, analytic projection gradient, progress bar; no caching |
| **GWAS** association scans (`gwas` branch) | done: the LiMMBo port (GWAS_PLAN.md A0-A5), checked against dense GLS and limmbo; calibration and power measured; `E` redefined after the first calibration run |

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

### Association scans (`gwas` branch, 2026-11-08)

`indri_gwas()` on simulated sibs (100 families of 5, N = 500), G from
`random_covariance(method = "partial")`, E LKJ(2), sex in `X`, `fit_G()` with
sommer; 20 replicates. Results in `results/gwas_calibration_P*.csv` and
`results/gwas_power_P*.csv` (`gwas_calibration.R`, `gwas_power.R`).

**Calibration** (5000 null SNPs per replicate):

| P | test | lambda_gc | type-I 1e-2 / 1e-3 / 1e-4 | delta |
|---|---|---|---|---|
| 10 | multi-trait LMM, indri G and E | 0.986 | 0.0076 / 0.00065 / 0.00007 | 1.04 |
| 10 | multi-trait LMM, true G and E | 1.005 | 0.0101 / 0.0010 / 0.00013 | 0.99 |
| 10 | single-trait LMM | 1.006 | 0.0104 / 0.0011 / 0.00012 | |
| 10 | multi-trait LM | 1.358 | 0.076 / 0.018 / 0.0044 | |
| 10 | single-trait LM | 1.501 | 0.035 / 0.0074 / 0.0015 | |
| 30 | multi-trait LMM, indri G and E | 0.976 | 0.0058 / 0.00048 / 0.00005 | 1.06 |
| 30 | multi-trait LMM, true G and E | 1.004 | 0.0110 / 0.0011 / 0.00010 | 1.00 |
| 30 | single-trait LMM | 1.008 | 0.0103 / 0.0010 / 0.00011 | |
| 30 | multi-trait LM | 1.233 | 0.080 / 0.017 / 0.0035 | |

- The mixed models are calibrated (indri's `G` and `E` slightly
  conservative); the linear models are inflated by the relatedness they
  ignore. The worst replicate of the multi-trait LMM had 2-3 of 5000 null
  SNPs at 1e-4 (0.5 expected).
- **The first run found `E` broken at 30 traits.** With `E` assembled from
  the blocks, the multi-trait LMM ran at 0.0066 at 1e-4 (66x nominal), worst
  replicate 0.051, with `delta` 0.30-0.94. The block average of `E` was
  strongly indefinite (smallest eigenvalue -2.9 on seed 4) and its projection
  had eigenvalues near 5e-6 against a truth of 0.006-0.08; true `G` with that
  `E` failed, fitted `G` with true `E` was calibrated. `fit_G()`'s `E` is now
  the residual covariance after `X` minus `c G` (floored at 1e-3 of the
  largest eigenvalue), which gives the table above; on a balanced sibship it
  equals the full-model REML `E` to 5e-4. 10 traits were fine either way.

**Power** at alpha = 1e-5, effect `h` = share of `trace(G + E)` explained by
the SNP, 200 causal SNPs per cell (multi-trait LMM / best single-trait LMM
with the effective-tests adjustment):

| pattern | P = 10, h = 0.25% | P = 10, h = 0.5% | P = 30, h = 0.1% | P = 30, h = 0.25% |
|---|---|---|---|---|
| one trait | 0.28 / 0.06 | 0.81 / 0.43 | 0.86 / 0.08 | 1.00 / 0.77 |
| three traits | 0.30 / 0 | 0.79 / 0.05 | 0.72 / 0 | 1.00 / 0.12 |
| all, same sign | 0.36 / 0 | 0.75 / 0.01 | 0.81 / 0 | 1.00 / 0 |
| all, random signs | 0.33 / 0 | 0.77 / 0.05 | 0.79 / 0 | 1.00 / 0.01 |
| against the phenotypic correlation | 0.98 / 0 | 1.00 / 0.08 | 1.00 / 0 | 1.00 / 0.01 |
| dominance (dominance test) | 0.27 | 0.68 | 0.61 | 0.97 |

- The joint test is the more powerful one in every pattern here, even for an
  effect on a single trait, and most against the correlation. **That edge is
  inflated by the simulation:** `random_covariance(method = "partial")` gives
  a nearly singular `G` at these sizes (smallest eigenvalue ~1e-8 at 30
  traits), so each trait's genetic value is almost determined by the others
  and conditioning on them removes most of its noise. The phenotypic
  correlations are moderate (median |r| 0.16 at 30 traits).
- Bonferroni and effective-tests adjustment gave nearly the same power.
- The dominance test finds overdominant SNPs that the additive tests miss
  (multi-trait additive 0.02-0.32 in the same cells).

**Second round** (LKJ(2) G, well conditioned; 20 replicates unless stated):

| P, N | test | lambda_gc | type-I 1e-2 / 1e-3 / 1e-4 | delta |
|---|---|---|---|---|
| 30, 500 | multi-trait LMM, indri G and E | 0.993 | 0.0067 / 0.00044 / 0.00004 | 1.13 |
| 30, 500 | same, E from a Schafer-Strimmer shrunk residual covariance | 0.993 | 0.0075 / 0.00051 / 0.00003 | 1.03 |
| 100, 300 | multi-trait LMM, indri G and E | 0.902 | 0.00011 / 0.00001 / 0 | 1.49 |
| 100, 300 | same, shrunk E | 0.883 | 0.00025 / 0.00001 / 0 | 0.96 |
| 100, 300 | multi-trait LMM, true G and E | 1.006 | 0.0116 / 0.0013 / 0.0001 | 0.99 |

- **Shrinking E buys nothing.** At 30 traits both are calibrated; at 100
  traits and 300 individuals both are very conservative (1% of the nominal
  rate at 1e-2) while the truth is calibrated, so the loss is in `G`, which
  is noise-dominated at that size (below).
- Power with the LKJ `G`: the joint test still leads in every pattern, by
  less for a single trait (30 traits, h = 0.25%: 0.87 against 0.675; 100
  traits, h = 0.1%: 0.66 against 0.38), and keeps its power at 100 traits
  despite being conservative there. `h` is a share of `trace(G + E)`, so the
  same `h` is a larger effect per trait at more traits: compare tests within
  a column, not across trait counts.

**Shrinking the assembled genetic correlations** (`shrinkage_G.R`, sommer,
N = 300; Ledoit-Wolf intensity from the per-entry variance of
`corrG_draws`; scores against the true genetic correlations):

| G, K | method | lambda | median err | strongest 10% err | slope | rel. Frobenius | integration |
|---|---|---|---|---|---|---|---|
| partial, 35 | none | | **0.120** | **0.090** | 0.94 | **0.41** | 1.05 |
| partial, 35 | towards zero | 0.27 | 0.149 | 0.267 | 0.69 | 0.44 | 0.57 |
| partial, 35 | towards phenotypic | 0.78 | 0.164 | 0.352 | 0.58 | 0.49 | 0.41 |
| LKJ, 35 | none | | 0.137 | **0.153** | 0.86 | 1.23 | 2.24 |
| LKJ, 35 | towards phenotypic | 1.00 | **0.093** | 0.178 | 0.47 | **0.84** | 0.65 |
| LKJ, 100 | none | | 0.121 | 0.140 | 0.84 | 1.79 | 3.90 |
| LKJ, 100 | towards phenotypic | 1.00 | **0.065** | **0.114** | 0.48 | **0.98** | 0.92 |

(Towards the mean correlation behaved like towards zero.)

- **With a strong genetic structure, shrinkage only hurts.** The
  intensity from the draws is too large there (the REML bands are
  conservative, 0.96-0.98 coverage), and the strongest correlations are
  pulled in by a factor of 3.
- **With a weak structure the estimate is mostly noise.** The relative
  error of the unshrunk `corrG` exceeds 1, and it *overstates* integration
  2.2x (35 traits) to 3.9x (100 traits): the opposite of the strong-`G`
  case above (STATUS "Integration", 0.75-0.85). The intensity goes to 1 for
  every target, so "shrinkage" is replacement, and the phenotypic
  correlations beat the assembled `G` on median and Frobenius error. That
  is a statement about these data sizes, not a reason to report `P`.
- **No shrinkage option for now.** A useful one would need an intensity
  that is not overconservative for strong structure, and a check of what it
  does to the association tests at K close to N.

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
- **Convergence warning** (stanimal): names the diagnostics (rank-normalised
  split-Rhat >= 1.05 on G, from `posterior::rhat()`, and divergent
  transitions, now recorded per block), says the likelihood may not identify
  the genetic correlations strongly enough, and suggests a larger
  `lkj_prior` instead of longer chains.
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
5. **The multi-trait test at K close to N** (100 traits, N = 300) is very
   conservative with indri's `G` (1% of the nominal rate at 1e-2): the
   assembled `G` is noise-dominated there. Try a regularised `G` in the scan
   (shrunk towards the phenotypic correlations, or a factor `G`) and score
   calibration and power; and record when `fit_G()`'s `G` should be trusted
   (per-entry draw variance against the spread of the estimate).
6. **Shrinkage intensity.** The Ledoit-Wolf intensity from the REML draws
   overshrinks a strong `G`; a cross-validated or bootstrap intensity, or
   stanimal's posterior draws, may do better. Bending (Hayes and Hill 1981;
   Meyer and Kirkpatrick 2010) remains the quantitative-genetics framing.

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
