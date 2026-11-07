# AGENTS.md — guidance for coding agents working on `indri`

`indri` is an R package that estimates the additive genetic covariance matrix
**G** (traits × traits) from a phenotype matrix `Y`, a design `X` and a
relationship matrix `A`. It fits a multivariate animal model on many small
trait subsets and reassembles the subset estimates into one positive
semi-definite matrix (LiMMBo's idea). The per-subset fit is a **pluggable
engine**: `stanimal` (Bayesian, posterior draws), `sommer` (REML) or `asreml`
(REML, sparse `A`).

## Read first

- [`STATUS.md`](STATUS.md) — where the work stands, what has been measured, and
  the next steps. This is the source of truth for priorities.
- [`NOTES.md`](NOTES.md) — the design, the engine contract, and the verified
  call forms for each engine. The engine sections record behaviour discovered
  by experiment; do not "simplify" it away.

## Commands

```sh
R CMD INSTALL package                     # install the package
Rscript -e 'devtools::document("package")'  # regenerate man/ after editing roxygen
Rscript -e 'devtools::test("package")'      # run the suite (asreml tests skip w/o licence)
```

Run the tests from the repository root. The suite is the spec: ~1,760
assertions cover the subset sampler, the assembly, the PSD projection, the three engines,
parallelism, scale consistency and diagnostics.

## Layout

| path | role |
|---|---|
| `package/R/fit_G.R` | public entry point `fit_G()`; global standardisation, core budget, diagnostics |
| `package/R/subsets.R` | `indri_subsets()` — deterministic trait-pair coverage sampler, with a greedy top-up that guarantees `min_cooccurrence` |
| `package/R/assemble.R` | `indri_assemble()`, `indri_assemble_draws()` (`"average"` / `"mixture"`), `.indri_weights()` |
| `package/R/project_psd.R` | `indri_project_psd()` — closest PSD matrix by least squares, analytic gradient |
| `package/R/a_forms.R` | relationship-matrix forms: dense, inverse, sparse inverse |
| `package/R/engines.R` | engine registry (`indri_engine()`, `indri_engines()`) |
| `package/R/engine-stanimal.R` / `engine-sommer.R` / `engine-asreml.R` | the three engines |
| `package/R/methods.R` | `print` / `summary` / `plot` methods for `indri_fit` |
| `package/tests/testthat/` | the spec |
| `examples/` | runnable studies; see `examples/README.md` |

## The engine contract

An engine is a list returned by `indri_engines()`-registry constructors with:

- `capabilities`: `draws` (bool), `se` (bool), `a_form` (`"dense"`, `"inverse"`
  or `"sparse_inverse"`), `missing` (bool), optionally `default_chains`.
- `fit(Y, X, a_forms, control)`: fits one subset and returns a list with `G`,
  `E`, and one of `draws` (an `s × s × D` array) or `se` (an `s × s` matrix),
  plus `converged` and any engine extras (`loglik`, `time`, `n_draws`).

Adding an engine = one `engine-*.R` file + a registry entry in `engines.R` + a
test file + a `Suggests:` line in `DESCRIPTION` if the dependency is optional.
Keep engines **not imported at package load**; check with `requireNamespace()`
at call time.

## Hard-won rules (read before touching the engines)

These cost real debugging time. The tests encode several of them; the rest are
documented because they are invisible from reading the code.

1. **Never `mclapply()` around cmdstanr.** Forked children conflict with the
   `processx` supervisor and kill the run silently. Block parallelism is PSOCK
   (`parallel::makeCluster` + `parLapply`).
2. **`X` is a trap as a variable name.** It is the first formal of
   `lapply`/`mclapply`/`parLapply`, so a worker parameter called `X` silently
   shifts the arguments. The block worker uses `design` for this reason.
3. **The PSD projection must keep its analytic gradient.** The objective
   optimises `K (K + 1) / 2` free entries of `L` with `C = L L'`; a numerical
   gradient evaluates the objective ~K²/2 times per iteration over every block
   and made K = 30 look like an engine "hang". See `project_psd.R`.
4. **sommer specifics.** `OMP_NUM_THREADS` must be set in the *shell before R
   starts* (`Sys.setenv` is too late); the formula must use **bare**
   `vsm`/`usm`/`ism` bound in a formula environment, never namespace-qualified;
   `Gu` (the inverse of `A`) needs dimnames matching the `id` levels; the
   injected `cores` budget must be filtered out of `mmes_control()`.
5. **asreml specifics.** `vm()` matches the inverse of `A` through a `rowNames`
   attribute; G/R components are extracted from `summary(fit)$varcomp` by name
   (`"^vm\("` and `"^units.*trait!trait"`), with `se` from the `std.error`
   column; the fit is two-pass (`asreml()` then `stats::update()`); tests skip
   when unlicensed.
6. **stanimal specifics.** It returns `draws`, not `se`; `fit_G()` derives
   inverse-variance weights from the per-entry posterior variance when `se` is
   absent. Block Rhat comes from `posterior::rhat`, and only for `CmdStanMCMC`.

## `fit_G()` behaviour to preserve

- Traits are **globally standardised once** before subsetting, so every engine
  sees identical numbers; the returned `G` is on that standardised scale
  (back-transform with `diag(y_sd) %*% G %*% diag(y_sd)`).
- Default subset size is capped:
  `.indri_default_s(K) = max(2, min(K - 1, floor(sqrt(K)), 4))`.
- `weights = "inverse_variance"` now feeds **both** the average and the
  projection (`indri_project_psd(..., weights = ...)`). It was measured to be a
  wash at N = 300, so equal weights stay the default — do not make weighting the
  default without a new measurement.
- **Never average block draws to get a spread.** Blocks are fitted to the same
  individuals, so their errors are nearly perfectly correlated while their
  chains are independent; averaging draw `d` over `c` blocks divides the spread
  by ~`sqrt(c)`. Measured: averaged draws cover 0.40 at nominal 0.95, the
  mixture assembly (`spread = "mixture"`: each entry from one random covering
  block) 0.89, the full model 0.80. Mixed draws are pointwise (a mixed matrix
  need not be PSD), and none of these is a joint posterior. The default is
  still `"average"` pending a `corrG_draws` output — see STATUS.md.
- `indri_subsets()` guarantees every pair is covered `min_cooccurrence` times;
  LiMMBo's construction alone does not (it only guarantees `ceiling(c / 2)`).
  Keep the top-up, and keep the test that sweeps seeds — a single seed passed
  while the guarantee was broken.

## Running studies (remote machine)

Long runs go to `meyer` via SSH; background Copilot commands on the remote
stall, so run remote commands over `ssh` and pull results back with `scp`.

The checkout is `~/projects/indri` (plural). It syncs through GitHub only:
commit and push locally first.

```sh
# sync + reinstall (meyer has no local changes, so a fast-forward suffices;
# check `git status --short` there first if in doubt)
ssh meyer 'cd ~/projects/indri && git pull --ff-only origin main && R CMD INSTALL package'

# run a study (async via the tool), e.g. the M1 gate or the coverage study
ssh meyer 'cd ~/projects/indri && INDRI_STUDY_SEEDS=1 INDRI_STUDY_CORES=32 Rscript examples/validation_12_traits.R'
ssh meyer 'cd ~/projects/indri && INDRI_COVERAGE_SEEDS=1,2,3 INDRI_COVERAGE_CORES=32 Rscript examples/uncertainty_coverage.R'

# long runs: one seed per process, detached so a dropped connection is harmless
ssh meyer 'cd ~/projects/indri && for s in 1 2 3; do INDRI_COVERAGE_SEEDS=$s INDRI_COVERAGE_CORES=48 nohup Rscript examples/uncertainty_coverage.R > results/uncertainty_coverage_seed$s.log 2>&1 & done'

# pull results back in ONE connection
ssh meyer 'cd ~/projects/indri/results && tar cf - <files>' | tar xf - -C results/
```

Remote gotchas:

- **meyer rate-limits SSH.** A burst of separate `scp`/`ssh` calls gets the
  host to refuse connections for a minute or so. Batch commands into one
  `ssh` and transfer files with one `tar` stream.
- **`pgrep -f <script>` matches its own wait loop.** A remote
  `while pgrep -f uncertainty_coverage.R; do sleep 30; done` never exits,
  because the loop's own command line contains the pattern. Wait on the
  summary CSVs appearing, or use `pgrep -f '^Rscript examples/...'`.
- The 12-trait coverage study takes ~3 minutes per seed at 48 cores.

Tracked results (`results/validation_12_traits_*.csv`,
`results/uncertainty_coverage_*.csv`) are whitelisted in `.gitignore`; the
benchmark overwrites them on meyer, so pull the CSVs back, commit them
locally, and discard meyer's modified copies (`git checkout -- results/` there)
before the next fast-forward. A *new* results file that you committed from
meyer's output is still untracked there and blocks the fast-forward ("would be
overwritten"); check it matches the committed copy (`md5sum`) and remove it.

## Conventions

- Keep changes surgical; update the roxygen docs and regenerate `man/` when the
  public API changes.
- Add tests for anything that failed or that encodes a discovered gotcha.
- Update `NOTES.md`/`STATUS.md` when a decision changes or a study is measured.
- Commit messages are plain English, one line, with a body when needed.
