# Status and next steps

Last updated: 2026-11-07. Design lives in [NOTES.md](NOTES.md); this file is the
short account of where the work stands, what has been measured, and what to do
next.

## Where things stand

| repo | branch | remote | what is on it |
|---|---|---|---|
| `indri` | `main` | `diogro/indri` (**private**) | the package, tests, examples, NOTES/STATUS |
| `stanimal` | `indri-support` | `diogro/stanimal` | `animalModelGQR.stan`, `eigen_A`, `test-eigen-A.R` (+480/-1 over master) |
| `limmbo` | `indri-session-fixes` | `diogro/limmbo` (fork of `meyer-lab-cshl/limmbo`) | pipeline fixes, pixi environment, AGENTS.md (+140/-34 over master) |

All three are pushed and the working trees are clean. Heads at the time of
writing: `indri` `51bac4d`, `stanimal/indri-support` `bdd30cf`, `limmbo`
`indri-session-fixes` `4edd228`; check with `git ls-remote --heads origin`
rather than trusting this line.

Cloning on another machine — note that two of the branches are not the default
branch of their repo:

```sh
git clone https://github.com/diogro/indri.git
git clone git@github.com:diogro/stanimal.git && cd stanimal && git switch indri-support
git clone git@github.com:diogro/limmbo.git  && cd limmbo  && git switch indri-session-fixes
```

`stanimal` also carries `master` (`39ece71`, the branch `indri-support` is cut
from) and older `feature/0.7.0`/`roadmap` branches. In the `limmbo` checkout,
`origin` is the fork; add the upstream if a PR is wanted:

```sh
git remote add upstream https://github.com/meyer-lab-cshl/limmbo
```

160+ test assertions pass in `indri` (the asreml tests skip where no licence is
available); 14 in the stanimal branch's new test file.

| Milestone | State |
|---|---|
| **M0** subset sampler, assembly, PSD projection, engine registry | done |
| **M1** `fit_G()` on the stanimal engine, validated against the full model | **done, gate passed** |
| **M2** sommer and ASReml engines | **done**: both implement the block contract; the two REML engines agree with each other and with stanimal |
| **M3** wire stanimal's `eigen_A` variant into an engine | prototype ready upstream, not wired |
| **M4** uncertainty, weights, diagnostics | partly: `G_draws` and inverse-variance weighting exist and are unit-tested, nothing validated |
| **M5** performance | partly: parallel blocks and the core budget work; the PSD projection now uses an analytical gradient (the former K=30 bottleneck), no caching |

## What has been measured

M1's gate: 12 traits, N = 300, s = 4, 46 six-trait subsets,
`min_cooccurrence = 3`, three seeds, `results/validation_12_traits_*.csv`.

| seed | indri median per-entry abs corr error | full 12-trait model | ratio |
|---|---|---|---|
| 1 | 0.0895 | 0.0912 | 0.98 |
| 2 | 0.0829 | 0.0964 | 0.86 |
| 3 | 0.1079 | 0.1557 | 0.69 |
| median | 0.0895 | 0.0964 | **0.93** |

Also better on q90 error (0.213 vs 0.248), relative Frobenius (0.236 vs 0.283)
and correlation-of-correlations (0.965 vs 0.960), at **4.6x the wall clock**
(5.9 vs 1.3 minutes).

The six-trait example ([examples/six_traits.R](examples/six_traits.R)) gives a
relative error of 0.237, against 0.369 for stanimal's full six-trait fit on the
same data.

### The cross-engine comparison

[examples/compare_engines.R](examples/compare_engines.R) fits one simulated
dataset with indri on the stanimal, sommer and asreml engines and reports the
correlation-of-correlations against the truth.

| K (seed 1, N = 300) | indri+stanimal | indri+sommer | indri+asreml | full stanimal |
|---|---|---|---|---|
| 10 | 0.963 | 0.943 | 0.942 | — |
| 30 | 0.951 | 0.929 | — | 0.911 |

At K = 10 the two REML engines agree with each other to 1.000 and with stanimal
to 0.98. At K = 30, indri+stanimal (0.951, 2.5 min) beats the full 30-trait
model (0.911, 24 min, reduced to 250 iterations) on both accuracy and wall
clock. The K = 30 full-model comparison ran without sommer/asreml because the
PSD projection was still the numerical-gradient version at the time; asreml at
K = 30 has not been measured.

### What those numbers prove, and what they do not

Proven: the pipeline is correct. A block from `fit_G()` is *bit-identical* to
`lmm_animal()` fitted on the same globally standardised columns
(`test-scale-consistency.R`), two blocks that share traits agree about them, the
two REML engines agree with each other and with stanimal, and at K = 30 the
assembled matrix beats the full model on both accuracy and time.

Not proven: how general that crossover is. It rests on a single K = 30 seed, so
the honest statement is "at K = 30 the method wins on one dataset", not yet "the
method wins". Repeat across seeds (and push to K = 40) before citing a
threshold.

## Next, in priority order

1. **Multi-seed crossover (K = 30 and 40).** The single-seed K = 30 result is
   promising but is one dataset. Run `compare_engines.R` across several seeds and
   at K = 40 to turn it into a statement. At K = 30 a stanimal run costs ~2.5
   minutes per seed at 64 cores.
2. **Uncertainty is unvalidated.** `G_draws` is documented as "an assembled
   spread, not a posterior". That is currently an assertion. M4's coverage study
   (fraction of true entries inside the 95% band, reported per coverage count)
   is what turns it into a claim, or into an honest caveat.
3. **Two code paths have never run end to end**: `weights =
   "inverse_variance"` (unit tests only; the study's weights run is gated behind
   `INDRI_STUDY_WEIGHTS=1` and has not been executed) and the missing-data policy
   (no fit with gaps in `Y` has gone through `fit_G()`).
4. **asreml at K = 30, and sommer's per-block cost.** asreml is validated at
   K = 10 (agrees with sommer to 1.000) but unmeasured at K = 30; sommer's
   `mmes` per-block cost degrades on near-boundary bivariate pairs at large K,
   which is worth characterising before recommending block sizes.
5. **M3, stanimal's `eigen_A`.** Prototype ready upstream, not wired. The shared
   eigen-decomposition is the path to making the stanimal engine cheap enough at
   large N.

## Running the studies

```sh
# one seed per process; each costs roughly 8 CPU-minutes
INDRI_STUDY_SEEDS=1 INDRI_STUDY_CORES=2 Rscript examples/validation_12_traits.R
INDRI_STUDY_SEEDS=2 INDRI_STUDY_CORES=2 Rscript examples/validation_12_traits.R
INDRI_STUDY_SEEDS=3 INDRI_STUDY_CORES=2 Rscript examples/validation_12_traits.R

# the optional inverse-variance comparison, on seed 1
INDRI_STUDY_SEEDS=1 INDRI_STUDY_WEIGHTS=1 Rscript examples/validation_12_traits.R
```

Each run writes `results/validation_12_traits_{summary,entries}_seed<N>.csv`;
those files are tracked so results survive a dead run. Everything else under
`results/` is ignored.

`INDRI_STUDY_CORES` is the total process budget passed to `fit_G()`, not a
worker count. Each block supplies `chains` processes (four by default for
stanimal), so the budget is split: `cores = 8` with the default arguments is one
block with four chains in parallel, while `cores = 8` with
`engine_args = list(chains = 2)` is four blocks at a time. `chains` itself is
never changed by `fit_G()`.

The cross-engine comparison has its own entry point:

```sh
INDRI_ENGINE_K=10 Rscript examples/compare_engines.R
INDRI_ENGINE_K=30 INDRI_ENGINE_CORES=64 Rscript examples/compare_engines.R
```

`INDRI_ENGINE_S*`, `INDRI_ENGINE_CORES*` and `INDRI_ENGINE_{SOMMER,ASREML}`
flags select per-engine block sizes, core budgets and which engines run; see the
script header. Results go to `results/compare_engines_*_seed<N>.csv`
(git-ignored).

## On a bigger machine

- **The 15 CPU-minute cap was a property of the previous sandbox.** Long
  CPU-bound runs there died silently with exit code 1 and nothing in the log
  between chain messages. Verify whether the new machine has anything similar
  before assuming chunking is needed: `Rscript -e 'while(TRUE) {}'` for ten
  minutes tells you. If the cap is gone, the study can run all three seeds in
  one invocation.
- **Prefer processes to threads.** With many cores, run several
  `Rscript examples/validation_12_traits.R` processes with different
  `INDRI_STUDY_SEEDS` rather than raising `cores` on one. The seed files are
  written separately and can be combined afterwards.
- **PSOCK workers are the only safe block parallelism.** `parallel::mclapply()`
  around cmdstanr kills runs; that is not a guess, it was reproduced.
- **Memory** grows with block size, not with K: with `s = 4` each block is a
  six-trait model. At K = 40 the blocks are larger, so watch resident memory per
  worker rather than total trait count.
- Installing on a fresh machine:

  ```sh
  # stanimal, from the branch with the eigen_A prototype (or master for stock)
  git -C stanimal checkout indri-support
  R CMD INSTALL stanimal/package

  R CMD INSTALL indri/package
  Rscript -e 'cmdstanr::install_cmdstan()'   # once, if CmdStan is absent

  # verify
  cd indri && Rscript -e 'devtools::test("package")'
  ```

## Open decisions

From NOTES section 11, still open: default engine; whether stanimal stays
`Suggests:` (currently yes); whether the factor/horseshoe prior is exposed for
block fits; the missing-data default (`"error"` today); GPL-3 licence; and
whether to create the GitHub remote (`gh` is authenticated as `diogro`).

## Two failure modes worth remembering

Both cost real debugging time and neither is visible from reading the code:

- **`parallel::mclapply()` around cmdstanr kills the run silently.** The forked
  children conflict with the `processx` supervisor. Blocks use PSOCK workers.
- **`X` as an argument name in `lapply`/`mclapply`/`parLapply`** silently shifts
  the arguments, because `X` is their first formal. The block worker's parameter
  is called `design` for this reason.
- **The PSD projection used a numerical gradient.** It optimises `K (K + 1) / 2`
  free entries, so each L-BFGS-B iteration evaluated the objective ~K²/2 times,
  each looping over every block. At K = 30 with 872 bivariate sommer blocks this
  looked like the sommer engine "hanging": the fits were 38 s and the projection
  was 30 min. The fix is an analytical gradient.
- **sommer multithreads internally.** `mmes()` spawns OpenMP threads even for
  bivariate blocks; combined with PSOCK block parallelism it oversubscribes and
  slows each fit ~2.6x. Pin `OMP_NUM_THREADS = 1` and let `cores` drive the
  parallelism.
