# Status and next steps

Last updated: 2026-11-07. Design lives in [NOTES.md](NOTES.md); this file is the
short account of where the work stands, what has been measured, and what to do
next.

## Where things stand

| repo | branch | remote | what is on it |
|---|---|---|---|
| `indri` | `main` | `diogro/indri` (**private**) | the package, tests, examples, NOTES/STATUS |
| `stanimal` | `master` | `diogro/stanimal` | 0.8.0: every animal model takes the eigendecomposition of `A` as data; `lmm_animal(eigen_A =)` |
| `limmbo` | `indri-session-fixes` | `diogro/limmbo` (fork of `meyer-lab-cshl/limmbo`) | pipeline fixes, pixi environment, AGENTS.md (+140/-34 over master) |

All three are pushed. Heads at the time of writing: `stanimal/master`
`9e0d831` (0.8.0), `limmbo` `indri-session-fixes` `4edd228`; check with
`git ls-remote --heads origin` rather than trusting this line.

Cloning on another machine — note that the limmbo branch is not the default
branch of its repo:

```sh
git clone https://github.com/diogro/indri.git
git clone git@github.com:diogro/stanimal.git
git clone git@github.com:diogro/limmbo.git  && cd limmbo  && git switch indri-session-fixes
```

`stanimal`'s `indri-support` branch has been fast-forwarded into `master` and
is no longer needed; it also carries older `feature/0.7.0`/`roadmap` branches. In the `limmbo` checkout,
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
| **M3** share stanimal's eigendecomposition of `A` across blocks | **done**: stanimal takes `Q`/`lambda` as data in every animal model (stanimal 0.8.0); indri decomposes `A` once per `fit_G()` call and passes it as `eigen_A`. Saves one O(N^3) decomposition per block: ~60 s per block at N = 4000 |
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

*The indri+sommer column predates the sommer ordering fix (`4b5748a`) and is
affected by it at both K (trait10+ sorts before trait2); asreml and stanimal
are not. Rerun before quoting sommer here.*

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

## REML uncertainty, and the sommer ordering bug (2026-11-07)

**Bug, fixed (`4b5748a`).** The sommer engine built `factor(trait)`, whose
levels sort alphabetically, and sommer returns its matrices in level order,
so any block whose trait names do not sort in column order (`trait2`,
`trait10`) came back permuted: variances swapped for pairs, covariances on
the wrong pairs for larger blocks. Every sommer run with 10 or more traits
before that commit is affected: the K = 10/16/30 engine comparisons, the first
sommer block-size sweep, and the heritabilities (not the correlations) of the
sommer convergence study. asreml and stanimal were not affected.

**The corrected sweep overturns "pairs only".** sommer's apparent collapse at
s >= 3 was the bug. Paired against s = 2 within seed (c = 1):

| engine, K | s = 3 (4-trait blocks) | s = 4 | s = 5 |
|---|---|---|---|
| sommer, 16 (6 seeds) | median err -0.007 (5/6 better), same time | -0.009 (4/6), 2x | -0.014 (5/6), 5x |
| sommer, 30 (4 seeds) | -0.001 (1/4), same time | -0.003 (4/4), 1.4x | -0.003 (2/4), 2.6x |
| asreml, 16 (2 seeds) | -0.011 (2/2), 0.6x time | -0.004 (1/2) | |

**Defaults now** (`4d5959b`): REML engines s = 3, c = 3, where c counts
*distinct* trait combinations. `indri_subsets()` never repeats a subset (a
pair estimated twice in the same context adds nothing) and caps c at what the
block width allows (1 for pairs). Default s is clamped so blocks stay
narrower than K.

**REML uncertainty** (`49f0657`): each REML block reports the joint sampling
covariance of its G entries (sommer's natural-scale parameter vcov; asreml's
inverse AI matrix, identical numbers), `fit_G()` draws from that normal
approximation (`n_draws`, seeded in the main process) and assembles the draws
like stanimal's. `examples/uncertainty_coverage.R`, 12 traits, 6 seeds, 396
pairs, the same data as the stanimal study:

| method | 95% coverage | mean width | median abs err of corrG |
|---|---|---|---|
| sommer, s = 2, c = 1 | 0.957 | 0.687 | 0.107 |
| asreml, s = 2, c = 1 | 0.960 | 0.687 | 0.107 |
| **sommer, s = 3, c = 3 (default)** | **0.970** | 0.690 | **0.103** |
| **asreml, s = 3, c = 3 (default)** | **0.977** | 0.744 | **0.106** |
| stanimal mixture, s = 4, c = 3 | 0.894 | 0.519 | 0.098 |
| full stanimal model | 0.801 | 0.490 | |

- The REML bands are calibrated with pairs and slightly conservative with
  three contexts per pair: the assembled estimate averages three contexts, so
  it is less noisy than any one block, while the band carries one block's
  sampling variance (widths are unchanged, accuracy improves). Getting
  exactly 0.95 would need the correlation between contexts, which share all
  the data; conservative is the honest side to err on.
- A **recentred mixture** (each draw a block's deviation from its own mean,
  added to `corrG`) was scored on the saved draws and rejected: REML 0.975
  (no change), stanimal 0.881 (worse).
- stanimal's bands are narrower and undercover; its point estimate is the
  most accurate. The difference is the prior, not the assembly.

## Non-converged REML pairs: measured (2026-11-07)

`examples/reml_convergence.R`: every pair fitted once (s = 2, c = 1), then
each pair whose fit reports no convergence refitted five ways. sommer K = 16
(6 seeds) and K = 30 (4 seeds); asreml K = 16 (2 seeds) and K = 30 (1 seed).

| | asreml: 26 of 678 pairs (3.8%) | sommer: 59 of 2,470 (2.4%) |
|---|---|---|
| abs error, converged vs not (median) | 0.107 vs 0.099 | 0.098 vs 0.115 (same q90) |
| more iterations | 42% then certified; estimate moves 0.004 | stops at ~29 anyway; no change |
| restart / h2 = 0.5 / univariate starts | 31-73% certified; estimate moves 0.006-0.017; error no better | identical estimate from every start |
| pair inside a 3-trait block | never certified; no better | worse: 38 pairs worse, 19 better |

- **"Not converged" means "optimum reached, not certified".** sommer stalls
  ("AI-first optimization stalled ... convergence was not certified") with
  the log-likelihood changing by ~1e-6 per step, and three very different
  starts end at the same log-likelihood (to 5e-6) and the same estimate. The
  flagged pairs are the hard ones (median true |r| 0.60 vs 0.38; some at
  |r| ~ 1 or h2 near a boundary), not badly estimated ones.
- **No refit improves accuracy**, and the occasional larger block makes it
  worse. Dropping a pair is not an option either: with c = 1 its entry would
  then rest on nothing.
- **Decision: report, don't refit.** `fit$diagnostics$unconverged` lists the
  blocks and their traits, `print()` shows the count ("not certified" for
  REML), and only stanimal blocks with Rhat >= 1.05 raise a warning
  (`options(indri.warn_rhat = FALSE)` silences it). `diagnostics$converged` is
  now `NA`, not `TRUE`, when no block can tell (REVIEW §4).
- The REML engines now fit single-trait blocks (formulas without the `trait`
  factor) and take starting values (`engine_args = list(init = list(G =,
  E =))`); asreml also takes `maxit` and `updates`.

## REML block size: measured, defaults changed (2026-11-07)

`examples/block_size_sweep.R`, s in {2, 3, 4, 5} x c in {1, 3}, same data
within each seed. sommer: K = 16 (6 seeds), K = 30 (4 seeds). asreml: K = 16
(seeds 1-2), K = 30 (seed 1). Correlation scale; `conv` is the share of blocks
whose fit reports convergence.

| engine, K | s (block width) | conv | median abs err | corr of corr |
|---|---|---|---|---|
| asreml, 16 | 2 (2) / 3 (4) / 4 (6) | 0.91 / 0.55 / 0.09 | 0.111 / 0.100 / 0.107 | 0.950 / 0.953 / 0.945 |
| asreml, 30 | 2 / 3 / 4 | 0.99 / 0.80 / 0.28 | 0.101 / 0.103 / 0.111 | 0.930 / 0.930 / 0.923 |
| sommer, 16 | 2 / 3 / 4 / 5 | 0.90 / 0.57 / 0.12 / 0.02 | 0.106 / 0.355 / 0.371 / 0.383 | 0.950 / 0.20 / -0.04 / 0.02 |
| sommer, 30 | 2 / 3 / 4 / 5 | 0.98 / 0.80 / 0.36 / 0.08 | 0.105 / 0.256 / 0.282 / 0.313 | 0.929 / 0.33 / 0.25 / 0.13 |

(c = 1 rows; c = 3 is within 0.01 everywhere.)

*Superseded: the sommer rows were corrupted by the ordering bug, and the
conclusions below were revised. See "REML uncertainty, and the sommer
ordering bug" above for the corrected sweep and the current defaults.*

- **No gain from larger blocks for REML.** asreml is level across s (the s = 3
  edge at K = 16 is 2 seeds and within noise), and its convergence falls with
  block size. sommer's multi-trait `mmes` fits mostly fail to converge and the
  estimate collapses. Pairs are the only setting where both engines converge.
- **c = 1 with pairs.** At s = 2 every pair was fitted exactly c times on the
  same data, so the repeats are identical numbers: c = 1 matches c = 3 in every
  seed (|difference| < 0.0002) at a third of the blocks (sommer K = 30: 66.7 s
  to 24.5 s; asreml K = 30: 351 s to 118 s).
- **Defaults now depend on the engine** (`.indri_engine_defaults()`): sommer and
  asreml use s = 2, c = 1; stanimal keeps `min(floor(sqrt(K)), 4)` and c = 3
  (its repeats are independent chains, and its larger blocks are measured).
- **asreml runs two workers by default.** The licence allows two concurrent
  processes, on one machine at a time; a third, or any process on a second
  machine while the first holds a seat, fails with "All licenses in use (-22)".
  `fit_G()` drops to one in-session worker when asreml is already loaded (the
  session then holds a seat) and turns a worker licence failure into an
  actionable error.
- Earlier engine comparisons used s = 2 for the REML engines explicitly, so
  their numbers are unaffected; their c = 3 only tripled the cost.

## M3, shared eigendecomposition: measured (2026-11-07)

`examples/eigen_A_timing.R`: one six-trait block, four chains, 200 iterations,
median of two runs, reference BLAS/LAPACK.

| N | old: Stan, per chain | new: R, once per fit | shared `eigen_A` | slowest chain |
|---|---|---|---|---|
| 500 | 7.5 s | 7.6 s | 7.5 s | 5.7 s |
| 1000 | 12.0 s | 12.3 s | 11.4 s | 9.0 s |
| 2000 | 32.9 s | 34.8 s | 27.4 s | 23.0 s |
| 4000 | 107 s | 111 s | **50 s** | 38.3 s |

- The models are unchanged: the slowest chain takes the same time in every
  mode. All the difference is fixed overhead outside sampling.
- For a single stanimal fit, "new" is no faster than "old": stock stanimal ran
  its four per-chain decompositions in parallel, so the wall-clock cost was
  about one decomposition either way (the saving there is CPU, 4x).
- Sharing removes the decomposition from every block: ~60 s per block at
  N = 4000, ~7 s at N = 2000, nothing measurable below N = 1000. A 46-block
  fit at N = 4000 drops from roughly 82 to 39 block-minutes. The decomposition
  grows as N^3 and sampling roughly as N, so the ratio grows with N.
- `eigen()` here runs on reference LAPACK (62 s at N = 4000); an optimised
  BLAS would shrink both the cost and the saving.
- Blocks fitted with the shared decomposition are bit-identical to a direct
  `lmm_animal()` fit (`test-scale-consistency.R`), and stanimal fits are
  bit-identical to 0.7.4 (a custom model that decomposes `A` the old way
  reproduces the built-in one draw for draw).

## Subset coverage bug: fixed (2026-11-07)

`indri_subsets()` did not guarantee `min_cooccurrence`: the LiMMBo
construction covers each pair only `ceiling(c / 2)` times by design, and about
half of all settings fell short (always for s = 2, the default below K = 9).
It now tops up with greedy subsets. The M1 runs (K = 12, s = 4) were not
affected. See [REVIEW.md](REVIEW.md) §1.

## Next, in priority order

0. *Done (2026-11-07):* `spread = "mixture"` is the default and `fit_G()`
   returns `corrG_draws` (mixed per-block correlation draws, the measured
   0.89). The REML engines still return no uncertainty; a Gaussian ensemble
   from their `se` would need the same mixture treatment and its own
   coverage check.
1. *Dropped:* a REML full-model baseline. Full multi-trait REML fits are too
   slow and rarely converge without curated starting values, so the
   full-model comparisons stay Bayesian, with the LKJ caveat of REVIEW.md §3.
2. **Missing data end to end.** The missing-data policy has still never gone
   through `fit_G()` with gaps in `Y`; `"complete_block"` is implemented but
   unexercised on real fits.
3. **asreml at K = 30, and sommer's per-block cost.** asreml is validated at
   K = 10 (agrees with sommer to 1.000) but unmeasured at K = 30; sommer's
   `mmes` per-block cost degrades on near-boundary bivariate pairs at large K,
   which is worth characterising before recommending block sizes.
4. **The next N-side cost: data transfer.** With the decomposition shared, the
   remaining per-block overhead at N = 4000 is ~12 s, most of it CmdStan
   writing and reading the N x N `Q` as JSON for every block. Worth measuring
   before optimising.

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
# the full-model posterior); one seed per process, ~3 minutes each at 48 cores
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
  # stanimal >= 0.8.0 (master) shares the eigendecomposition across blocks;
  # older versions still work, without the saving
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
