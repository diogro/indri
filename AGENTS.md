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

## Running studies

The studies in `examples/` fit many Stan models; run anything beyond the test
suite on a machine with spare cores, not a laptop. If you use a remote
machine, sync it through git (commit and push locally, then fast-forward the
remote checkout) rather than copying files.

```sh
# one seed per process, detached so a dropped connection does not kill the run
for s in 1 2 3; do
  INDRI_COVERAGE_SEEDS=$s INDRI_COVERAGE_CORES=48 \
    nohup Rscript examples/uncertainty_coverage.R \
    > results/uncertainty_coverage_seed$s.log 2>&1 &
done
```

The `INDRI_*_CORES` variables are a total process budget for `fit_G()`, not a
worker count; see STATUS.md, "Running the studies". The 12-trait coverage study
takes about 3 minutes per seed at 48 cores.

Gotchas:

- **`pgrep -f <script>` matches its own wait loop.** A
  `while pgrep -f uncertainty_coverage.R; do sleep 30; done` run through
  `ssh` or `bash -c` never exits, because the loop's own command line contains
  the pattern. Wait for the summary CSVs to appear instead, or anchor the
  pattern (`pgrep -f '^Rscript examples/'`).
- When pulling results back from a remote machine, batch them into one
  transfer (`ssh host 'tar cf - <files>' | tar xf -`); many separate
  connections can trip SSH rate limiting.

Tracked results (`results/validation_12_traits_*.csv`,
`results/uncertainty_coverage_*.csv`) are whitelisted in `.gitignore`. A study
overwrites them where it runs, so on a remote checkout: copy the CSVs back,
commit them locally, then discard the remote's copies (`git checkout --
results/` for modified files; remove new ones after checking they match the
committed copy) before the next fast-forward, which otherwise refuses to
overwrite them.

## Conventions

- Keep changes surgical; update the roxygen docs and regenerate `man/` when the
  public API changes.
- Add tests for anything that failed or that encodes a discovered gotcha.
- Update `NOTES.md`/`STATUS.md` when a decision changes or a study is measured.
- Commit messages are plain English, one line, with a body when needed.
