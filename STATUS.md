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

All three are pushed. Heads of the two side branches at the time of writing:
`stanimal/indri-support` `bdd30cf`, `limmbo` `indri-session-fixes` `4edd228`;
check with `git ls-remote --heads origin` rather than trusting this line.

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

About 1,760 test assertions pass in `indri` (the asreml tests skip where no licence is
available); 14 in the stanimal branch's new test file.

| Milestone | State |
|---|---|
| **M0** subset sampler, assembly, PSD projection, engine registry | done |
| **M1** `fit_G()` on the stanimal engine, validated against the full model | **done, gate passed** |
| **M2** sommer and ASReml engines | **done**: both implement the block contract; the two REML engines agree with each other and with stanimal |
| **M3** wire stanimal's `eigen_A` variant into an engine | prototype ready upstream, not wired |
| **M4** uncertainty, weights, diagnostics | mostly: inverse-variance weights measured (a wash, equal stays default); interval coverage measured — the averaged `G_draws` is overconfident (0.40 at nominal 0.95), the mixture assembly reaches 0.89 against 0.80 for the full model (see below). Left: a `corrG_draws` output and the `spread` default |
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

## Weights: measured (2026-11-07)

The M4 weights question is closed for the current sample size. On the 12-trait
seed-1 benchmark, inverse-variance weighting does not beat equal weights:

| method | median abs err | q90 abs err | corr of corr | rel Frobenius |
|---|---|---|---|---|
| indri, equal weights | 0.0866 | 0.213 | 0.965 | 0.237 |
| indri, inverse-variance | 0.0894 | 0.209 | 0.964 | 0.237 |
| full 12-trait model | 0.0943 | 0.248 | 0.959 | 0.285 |

The two indri variants agree within Monte Carlo noise, and both beat the full
model. Conclusion: equal weights stay the default (paper fidelity);
inverse-variance remains available but is not advertised as better — noisy
per-block uncertainties and dependent blocks mean the weighting adds nothing at
N = 300.

## Uncertainty coverage: measured (2026-11-07)

The M4 uncertainty question is now measured, not asserted. On the 12-trait
simulation (`examples/uncertainty_coverage.R`, three seeds), the 95% band of
indri's assembled `G_draws` covers the true correlation in **~43%** of entries
(per seed: 53%, 48%, 29%), against **~81%** for the full model (85%, 82%, 77%).
Neither is pointwise calibrated at the nominal 95%, but indri's assembled spread
is substantially more overconfident: its mean band width is ~0.17 versus ~0.46
for the full model, and coverage *falls* as a pair's coverage count rises (the
`10+` band runs 25-35%), because averaging more blocks narrows the spread.

Conclusion at the time: `G_draws` is a relative spread, not a credible
interval. *Superseded by the next subsection*: the narrowing is an artefact of
averaging, and the mixture assembly removes it.

### The fix: a mixture assembly (2026-11-07, six seeds)

The overconfidence is the averaging. Draw `d` averaged over `c` independent
block chains has ~1/c of a block's variance, but the blocks see the same data
and their errors are nearly perfectly correlated. Scoring three assemblies of
the *same* fits (`f3abee8`, seeds 1-6, 396 pairs, correlation scale):

| arm | 95% coverage | mean width | coverage 7-9 / 10+ |
|---|---|---|---|
| `average` (old `G_draws`) | 0.40 | 0.16 | 0.51 / 0.33 |
| `mixture` (entry from one random covering block) | **0.89** | 0.52 | 0.91 / 0.88 |
| `inflated` (average, deviations x sqrt(coverage)) | 0.90 | 0.51 | 0.92 / 0.89 |
| full 12-trait model | 0.80 | 0.49 | — |

- The mechanism is confirmed: scaling the average by sqrt(coverage) alone
  recovers the coverage, and coverage no longer falls with pair coverage.
- `mixture` and `inflated` agree, so the between-block spread is small; the
  gain is entirely from not dividing by `c`.
- Mixture bands are about as wide as the full model's but cover more (0.89 vs
  0.80), consistent with indri's point estimates being better centred.
- Still below nominal 95%. The full model undercovers too, so part of the gap
  is likely shared (prior shrinkage at N = 300, short chains: blocks run one
  chain of 400 iterations).
- Mixed draws are pointwise: entries are drawn independently, so a mixed
  matrix need not be PSD. Fine for per-entry intervals, not for functionals of
  the whole matrix (eigenvalues, response to selection).

`fit_G(spread = "mixture")` is available; the default is still `"average"`.
Per-seed results: `results/uncertainty_coverage_summary_seed{1..6}.csv`.

## Subset coverage bug: fixed (2026-11-07)

`indri_subsets()` did not guarantee `min_cooccurrence`: the LiMMBo
construction covers each pair only `ceiling(c / 2)` times by design, and about
half of all settings fell short (always for s = 2, the default below K = 9).
It now tops up with greedy subsets. The M1 runs (K = 12, s = 4) were not
affected. See [REVIEW.md](REVIEW.md) §1.

## Next, in priority order

0. **Finish the mixture assembly.** What was measured is the mixture of each
   block's *correlation* draws; `fit_G(spread = "mixture")` mixes
   *covariance* draws, whose `cov2cor()` can leave [-1, 1]. Add a
   `corrG_draws` output built per block, as in the study, then make
   `spread = "mixture"` the default (recommended) and document both outputs as
   pointwise intervals (a mixed draw need not be PSD).
1. **A REML full-model baseline.** The stanimal full model's LKJ prior shrinks
   each correlation harder as K grows (prior sd 0.38 at 6 traits, 0.28 at 12,
   0.18 at 30), so "indri beats the full model" partly compares priors. Refit
   the K = 12 and K = 30 references with sommer or asreml on all traits. See
   [REVIEW.md](REVIEW.md) §3.
2. **Missing data end to end.** The missing-data policy has still never gone
   through `fit_G()` with gaps in `Y`; `"complete_block"` is implemented but
   unexercised on real fits.
3. **asreml at K = 30, and sommer's per-block cost.** asreml is validated at
   K = 10 (agrees with sommer to 1.000) but unmeasured at K = 30; sommer's
   `mmes` per-block cost degrades on near-boundary bivariate pairs at large K,
   which is worth characterising before recommending block sizes.
4. **M3, stanimal's `eigen_A`.** Prototype ready upstream, not wired. The shared
   eigen-decomposition is the path to making the stanimal engine cheap enough at
   large N.

*Dropped*: the multi-seed K=30/40 crossover. That run was about how general
the accuracy/time crossover is, not about correctness, so without it the K = 30
result is a single-seed observation: do not cite "the method wins at K = 30".

## Running the studies

```sh
# one seed per process; each costs roughly 8 CPU-minutes
INDRI_STUDY_SEEDS=1 INDRI_STUDY_CORES=2 Rscript examples/validation_12_traits.R
INDRI_STUDY_SEEDS=2 INDRI_STUDY_CORES=2 Rscript examples/validation_12_traits.R
INDRI_STUDY_SEEDS=3 INDRI_STUDY_CORES=2 Rscript examples/validation_12_traits.R

# the optional inverse-variance comparison, on seed 1
INDRI_STUDY_SEEDS=1 INDRI_STUDY_WEIGHTS=1 Rscript examples/validation_12_traits.R

# the uncertainty coverage study (average / mixture / inflated indri bands vs
# the full-model posterior); one seed per process, ~3 minutes each on meyer
INDRI_COVERAGE_SEEDS=1 INDRI_COVERAGE_CORES=48 Rscript examples/uncertainty_coverage.R
```

The coverage study also saves its fits (block draws, subsets, truth, full-model
correlation draws) to `results/uncertainty_coverage_fits_seed<N>.rds`
(git-ignored), so a new assembly can be scored without refitting.

Each run writes `results/validation_12_traits_{summary,entries}_seed<N>.csv`
and `results/uncertainty_coverage_summary_seed<N>.csv`; those files are tracked
so results survive a dead run. Everything else under `results/` is ignored.

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
the `spread` default (`"average"` today, `"mixture"` recommended — see Next,
item 0).

## Failure modes worth remembering

Each cost real debugging time and none is visible from reading the code:

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
