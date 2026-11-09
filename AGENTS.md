# AGENTS.md — guidance for coding agents working on `indri`

`indri` is an R package that estimates the additive genetic covariance matrix
**G** (traits × traits) from a phenotype matrix `Y`, a design `X` and a
relationship matrix `A`. It fits a multivariate animal model on many small
trait subsets and reassembles the estimates into one positive semi-definite
matrix (LiMMBo's idea). The per-subset fit is a **pluggable engine**:
`stanimal` (Bayesian), `sommer` (REML) or `asreml` (REML, licensed).
The assembled `G`, with a residual covariance `E`, also drives genome-wide
association scans (`indri_gwas()`, since 0.3.0; design and decisions in
[`GWAS_PLAN.md`](GWAS_PLAN.md)).

## Read first

- [`STATUS.md`](STATUS.md): where the work stands, what has been measured, and
  the next steps. The source of truth for priorities and numbers.
- [`NOTES.md`](NOTES.md): the design, the engine contract and the reasons
  behind the defaults.

## Commands

```sh
R CMD INSTALL package                        # install
Rscript -e 'devtools::document("package")'   # regenerate man/ after roxygen edits
Rscript -e 'devtools::test("package")'       # the suite (asreml tests skip w/o licence)
```

Run from the repository root. The suite is the spec (~3,300 assertions):
subsets, assembly, projection, the three engines, parallelism, scale and units,
uncertainty, diagnostics, progress, association scans and the shipped data. **Reinstall before running tests that use
parallel blocks**: PSOCK workers load the installed package.

## Layout

| path | role |
|---|---|
| `package/R/fit_G.R` | `fit_G()`: standardisation, defaults, parallel blocks, REML draws, assembly, units, diagnostics, progress bar |
| `package/R/subsets.R` | `indri_subsets()`: set cover, greedy top-up, distinct subsets |
| `package/R/assemble.R` | `indri_assemble()`, `indri_assemble_draws()` (`"mixture"` / `"average"`), weights |
| `package/R/project_psd.R` | `indri_project_psd()`: least-squares PSD fit, analytic gradient |
| `package/R/a_forms.R` | forms of `A`: dense, inverse, sparse inverse, eigen |
| `package/R/engines.R` | registry, `init` checks, normal draws from `vcov` |
| `package/R/engine-stanimal.R`, `engine-sommer.R`, `engine-asreml.R` | the engines |
| `package/R/methods.R` | `print`, `summary`, `plot` for fits |
| `package/R/gwas.R` | `indri_gwas()`: rotation, REML `delta`, dominance, permutations and q-values; `indri_effective_tests()` |
| `package/R/gwas-methods.R` | `plot()` (ggplot2 Manhattan/QQ) and `summary()` for scans |
| `package/R/data.R`, `package/data-raw/lgsm.R` | the LG/SM F3 data sets and the script that builds them from `../mouseGrowthQTLs` |
| `package/tests/testthat/fixtures/limmbo/` | limmbo's association results on its test data, with the script that made them |
| `examples/` | studies and worked examples; see `examples/README.md` |

## The engine contract

An engine is a list with a `name`, `capabilities` and `fit(Y, X, a_forms,
control)`:

- `capabilities`: `draws`, `se`, `vcov` (logical), `a_form` (`"dense"`,
  `"inverse"`, `"sparse_inverse"`, `"eigen"`), `missing`, and optionally
  `default_s`, `default_cooccurrence`, `default_cores`, `default_chains`.
- `fit()` returns `G` and `E` (s × s, standardised scale, **in the column order
  of `Y`**), then `draws` (s × s × D) or `se` plus `vcov` (joint covariance of
  the lower-triangle G entries, row by row), and `converged`.

Adding an engine: one `engine-*.R`, a registry entry in `engines.R`, a test file
(with trait names that do not sort in column order), and a `Suggests:` line.
Never load an engine at package load; check with `requireNamespace()` in `fit()`.

## Hard-won rules

Each of these cost real debugging time; several are invisible in the code.

1. **Never `mclapply()` around cmdstanr.** Forked children conflict with its
   `processx` supervisor and the run dies silently. Blocks use PSOCK workers.
   (Forking is fine for pure-R post-processing.)
2. **`X` is a trap as an argument name** in `lapply`/`parLapply` workers (their
   first formal); the block worker uses `design`.
3. **Keep the projection's analytic gradient.** A numerical gradient costs
   `K (K + 1) / 2` evaluations per step and looked like an engine hang at K = 30.
4. **Engine output follows the column order of `Y`.** sommer returned blocks in
   alphabetical order (`factor()` sorts levels), which put variances and
   covariances on the wrong traits; trait levels follow the columns and the
   output is reordered by name.
5. **sommer**: `OMP_NUM_THREADS = 1` (set before R starts, or in the script
   before sommer loads); bare `vsm`/`usm`/`ism` bound in the formula
   environment; `Gu` dimnames matching the `id` levels; filter the injected
   `cores` out of `mmes_control()`. Single traits use formulas without `trait`.
6. **asreml**: `vm()` matches `A^-1` through a `rowNames` attribute; G and E come
   from `summary(fit)$varcomp` by name, `vcov` from `fit$ai`; two passes
   (`asreml()`, then `update()`), with a budget of `maxit = 50` and up to 5
   updates while unconverged: asreml's own 13 iterations stopped far from the
   optimum with covariates in `X` (`test-engine-asreml.R` checks it against
   sommer). Every process needs a licence seat: two on
   one machine at a time, and loading asreml takes one, so never load it in a
   main process that then starts two workers (check installation with
   `system.file(package = "asreml")`).
7. **stanimal**: the engine passes `lkj_prior = 1` unless `engine_args` sets it,
   and one shared eigendecomposition of `A` as `eigen_A`, in ascending order
   (blocks stay bit-identical to a direct `lmm_animal()` fit). Each block writes
   CmdStan output to its own directory, deleted after reading: cmdstanr names
   files from the minute, the PID and one `runif()`, so fits from the same
   random state collide.
8. **Decide capabilities in the main process**, not in `fit()`: workers load the
   installed package and may not have new helpers.
9. **`E` is not assembled from the blocks.** Their average is strongly
   indefinite at many traits (smallest eigenvalue -2.9 at 30), its projection
   has eigenvalues near zero, and the multi-trait scan ran at 66x its nominal
   rate at 1e-4. `fit_G()`'s `E` is the residual covariance after `X` minus
   `c G` (`c = tr(M A) / (n - J)`), floored at 1e-3 of the largest
   eigenvalue; `E_bar` keeps the block average for reference.
10. **The limmbo parity tests pass LIMIX's fitted `delta`** (`log_delta.csv`):
    LIMIX grid-searches it (refined for the multi-trait test, not for the
    single-trait one), indri fits it by REML. They check the scan, not the
    optimiser; `test-gwas-dense.R` checks indri's `delta` against dense REML.
    LIMIX stays out of user documentation.
11. **Cap BLAS threads before R starts** for anything forked or parallel:
    `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`. Setting them
    inside R is too late for BLAS; 80 forked study workers once ran 28,000
    threads. Count cores as processes x threads (`nlwp` in `ps`).

## `fit_G()` behaviour to preserve

- **`X` is checked like stanimal's**: full column rank, able to represent a
  constant, finite (`.indri_check_design()`); point users to
  `stanimal::animal_design()`.
- **One standardisation, one projection.** Traits are standardised once; fitting,
  weights, assembly, the projection and every draw run on that scale; only the
  assembled covariances (`G`, `G_bar`, `G_draws`, `sigma_G`) are converted to the
  units of `Y` at the end, never re-projected. `corrG` comes from the
  standardised `G`; `blocks` stay standardised.
- **Never average block draws to get a spread** (overconfident by ~`sqrt(c)`);
  `spread = "mixture"` is the default, and `corrG_draws` mixes each block's own
  correlation draws. Mixed draws are pointwise, not valid matrices.
- **REML uncertainty** is a joint normal approximation from each block's `vcov`,
  drawn in the main process (`n_draws`, seeded from `seed`).
- **Defaults come from the engine** (`.indri_engine_defaults()`,
  `.indri_default_s()`, `.indri_default_cores()`): s = 3 cap, c = 3 stanimal /
  1 REML, asreml 2 workers. Change them only with a new measurement (the
  studies behind them are in git history at `7b67eab`; STATUS.md).
- **`min_cooccurrence` counts distinct subsets**; `indri_subsets()` tops up to
  the guarantee, never repeats a subset, and caps c by block width. Keep the
  test that sweeps seeds: a single seed once passed while the guarantee was
  broken.
- **Non-converged REML blocks are reported, not refitted** (measured to be
  uncertified optima). stanimal blocks with rank-normalised split-Rhat >= 1.05
  on G or with divergent transitions warn, naming the diagnostic and
  suggesting a larger `lkj_prior` (longer chains did not help on ratones);
  `options(indri.warn_rhat = FALSE)` silences it, and the test helper sets
  it.
- **Equal weights and the least-squares projection** stay the defaults; their
  alternatives were measured and lost or tied.
- **Progress bar**: `progress = getOption("indri.progress", interactive())`; with
  parallel blocks it sends load-balanced rounds, without it one `parLapply()`.

## Running studies

The studies fit many models; run them on a machine with spare cores, synced
through git (commit and push, then fast-forward the remote checkout).

```sh
# one run per process, detached
for k in 10 16 30; do
  INDRI_ENGINE_K=$k INDRI_ENGINE_CORES=48 \
    nohup Rscript examples/compare_engines.R > results/compare_k$k.log 2>&1 &
done
```

Traps:

- **Thread caps** (rule 11): launch with
  `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`; the GWAS and
  shrinkage studies warn when the BLAS cap is missing.
- **Never wait on, or kill by, a process-name pattern from inside
  `ssh '...'`.** `pgrep -f`, `ps | grep "[x]yz"` and `pkill -f` match the remote
  shell running the command, whose command line contains the pattern: wait
  loops never end and `pkill` kills its own session. Wait for output files
  (newer than a marker written at launch) and kill by PID. `Rscript` re-execs as
  `R --file=...`, so match `file=examples/...`.
- **More than ~125 parallel blocks** needs `Rscript --max-connections=1024`
  (R allows 128 connections; each PSOCK worker holds one). If `makeCluster()`
  fails with "port ... cannot be opened", another service owns a port in R's
  random range (11000-11999): set `R_PARALLEL_PORT`.
- **Write stdout and data streams separately**: anything a remote command prints
  corrupts a `ssh host 'tar cf - ...' | tar xf -` transfer. Batch transfers into
  one connection; many connections can trip SSH rate limits.
- **Results**: `results/validation_12_traits_*.csv`,
  `results/uncertainty_coverage_*.csv`, `results/gwas_calibration_*.csv`,
  `results/gwas_power_*.csv` and `results/shrinkage_G_*.csv` are tracked. On a remote checkout, copy
  them back and commit locally, then discard the remote copies before the next
  fast-forward. A rename step must check its source exists, or a failed run
  passes off an old file as new.

## Conventions

- Keep changes surgical; regenerate `man/` when roxygen changes.
- Add a test for anything that failed or that encodes a gotcha.
- Update `NOTES.md` when a design decision changes and `STATUS.md` when a study
  is measured.
- Commit messages: a plain-English line, a body when needed.
