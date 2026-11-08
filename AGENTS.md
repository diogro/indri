# AGENTS.md — guidance for coding agents working on `indri`

`indri` is an R package that estimates the additive genetic covariance matrix
**G** (traits × traits) from a phenotype matrix `Y`, a design `X` and a
relationship matrix `A`. It fits a multivariate animal model on many small
trait subsets and reassembles the estimates into one positive semi-definite
matrix (LiMMBo's idea). The per-subset fit is a **pluggable engine**:
`stanimal` (Bayesian), `sommer` (REML) or `asreml` (REML, licensed).

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

Run from the repository root. The suite is the spec (~2,850 assertions):
subsets, assembly, projection, the three engines, parallelism, scale and units,
uncertainty, diagnostics, progress. **Reinstall before running tests that use
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
| `package/R/methods.R` | `print`, `summary`, `plot` |
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
  1 REML, asreml 2 workers. Change them only with a new measurement
  (`block_size_sweep.R`, `uncertainty_coverage.R`).
- **`min_cooccurrence` counts distinct subsets**; `indri_subsets()` tops up to
  the guarantee, never repeats a subset, and caps c by block width. Keep the
  test that sweeps seeds: a single seed once passed while the guarantee was
  broken.
- **Non-converged REML blocks are reported, not refitted** (measured to be
  uncertified optima). stanimal blocks with Rhat >= 1.05 warn
  (`options(indri.warn_rhat = FALSE)`; the test helper sets it).
- **Equal weights and the least-squares projection** stay the defaults; their
  alternatives were measured and lost or tied.
- **Progress bar**: `progress = getOption("indri.progress", interactive())`; with
  parallel blocks it sends load-balanced rounds, without it one `parLapply()`.

## Running studies

The studies fit many models; run them on a machine with spare cores, synced
through git (commit and push, then fast-forward the remote checkout).

```sh
# one seed per process, detached
for s in 1 2 3; do
  INDRI_COVERAGE_SEEDS=$s INDRI_COVERAGE_CORES=48 \
    nohup Rscript examples/uncertainty_coverage.R > results/cov_seed$s.log 2>&1 &
done
```

Traps:

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
- **Results**: `results/validation_12_traits_*.csv` and
  `results/uncertainty_coverage_*.csv` are tracked. On a remote checkout, copy
  them back and commit locally, then discard the remote copies before the next
  fast-forward. A rename step must check its source exists, or a failed run
  passes off an old file as new.

## Conventions

- Keep changes surgical; regenerate `man/` when roxygen changes.
- Add a test for anything that failed or that encodes a gotcha.
- Update `NOTES.md` when a design decision changes and `STATUS.md` when a study
  is measured.
- Commit messages: a plain-English line, a body when needed.
