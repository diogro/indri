# Code and plan review — 2026-11-07

A review of the package (`package/R`), NOTES.md, STATUS.md and the M4 coverage
study, at `943b366`. Findings are ranked by how much they change what the
results mean. Overall the code is in good shape: the engine contract is clean,
the global standardisation is careful, and the analytical PSD gradient is
correct (derivation checked against the objective).

Status of the findings (updated 2026-11-07):

| # | finding | state |
|---|---|---|
| 1 | subset sampler misses `min_cooccurrence` | **fixed** (`0837352`) |
| 2 | averaged `G_draws` overconfident | **diagnosed and fixed** in `f3abee8` / `96f0894`: mixture 0.89 vs average 0.40; left: `corrG_draws`, flip the default |
| 3 | LKJ prior confounds the full-model comparison | open — REML full-model baseline is STATUS item 1 |
| 4 | doc drift | docs fixed; asreml `engine_args` / formula and the all-NA diagnostics flag still open |

## 1. Bug: `indri_subsets()` does not guarantee `min_cooccurrence`

Swept 20 seeds over K in {6, 8, 10, 12, 20, 30, 40}, s in {2, 3, 4},
c in {1, ..., 5}. About half of the configurations leave some trait pair
covered fewer than `c` times. At `c = 3`:

| s | violations (of 20 seeds per K) | worst min coverage |
|---|---|---|
| 2 | every K, every seed | 2 |
| 3 | 1/20 at K = 8, 9/20 at K = 12, 20/20 at K >= 20 | 2 |
| 4 | 1-5/20 at K >= 20 | 2 |

This hits the defaults: `.indri_default_s()` gives s = 2 for K < 9 and s = 3
for K in 9-15. The 12-trait validation used `s = 4` explicitly and was not
affected.

Cause: `n_covers <- ceiling(min_cooccurrence / 2)`. Each pass of windows covers
each pair exactly once by construction, so the guarantee is only
`ceiling(c / 2)`; the remainder is incidental co-membership. The halving is
inherited from LiMMBo's `multiple_set_covers_all`
(`number_of_covers // 2 + 1`), so LiMMBo has the same defect.

The existing test (K = 8, s = 3, c = 3, seed = 1) passes on a lucky seed.

Fix: top up with greedily built subsets until every pair reaches `c`; have
`fit_G()` check the coverage of user-supplied subsets; test across a grid of
seeds rather than one.

**Fixed (2026-11-07).** `indri_subsets()` now appends greedily built subsets
of the same width until every pair reaches `c`. The base construction is
unchanged and the top-up only appends rows, so configurations that already met
`c` get identical subsets. That includes the M1 runs (K = 12, s = 4, c = 3:
no top-up on any of 20 seeds), so the published validation numbers stand.
Cost: 0-2 extra blocks at s >= 3 (10-15 at K = 30, s = 3, c = 5), but
+45% at s = 2, where a block is a single pair. Also fixed: the padding step
used `sample(available, need)`, which samples from `1:available` when only
one trait is available and can duplicate a trait. `fit_G()` now warns when
hand-built `subsets` leave a pair uncovered. The test loops over K, s, c and
five seeds.

## 2. `G_draws` is overconfident by roughly sqrt(coverage)

`indri_assemble_draws()` averages draw `d` across the blocks that contain an
entry. The draws come from independent chains, so the averaged spread has
variance about `var_block / c`. But every block is fitted to the same
individuals and the same trait columns, so their estimation errors are close to
perfectly correlated, and the real uncertainty of the average is close to
`var_block`.

Prediction for `examples/uncertainty_coverage.R`: indri bands narrower than
they should be by about `sqrt(c)` (2-3x at the coverage counts seen), with
undercoverage *worsening* as pair coverage rises. The full model will not show
this.

Alternative to test as a second arm: a mixture assembly — for each draw `d` and
entry, take the value from one randomly chosen covering block rather than the
average. That keeps the per-block spread and adds between-block spread.

**Confirmed (`4f70b8d`).** The coverage study measured ~43% coverage for the
indri band against ~81% for the full model, with the `10+` coverage band at
25-35% — the predicted pattern.

**Mixture measured (`f3abee8`, six seeds).** Mixture 0.89 coverage, width
0.52; sqrt(coverage)-inflated average 0.90; old average 0.40; full model 0.80.
Coverage is flat across pair-coverage bands. Details in STATUS.md.

The same dependence is a plausible explanation for the weights study coming out
as a wash: inverse-variance weighting pays off when errors are independent.

Smaller points on `G_draws`:

- it is centred on the unweighted `G_bar`, not on the reported projected `G`;
- it ignores `weights`.

Both should be fixed or documented.

## 3. The comparison against the full model is confounded by the LKJ prior

stanimal puts LKJ(eta) on the full correlation matrix. Under LKJ the marginal
sd of each correlation is `1 / sqrt(2 eta + d - 1)`, so the prior shrinks
harder as dimension grows. With eta = 1:

| model | dimension | prior sd of each correlation |
|---|---|---|
| indri block (s = 4) | 6 | 0.38 |
| full model | 12 | 0.28 |
| full model | 30 | 0.18 |

At N = 300, "indri beats the full 30-trait model" is partly a weaker-prior
effect, not only the method. The K = 30 full fit also ran 250 iterations; its
Rhat should be checked before it is used as a reference.

Fairer reference: a REML full model (sommer or asreml on all K traits), or an
explicit statement that the two Bayesian fits carry different effective
priors.

STATUS drops the multi-seed crossover because "the implementation is accepted
as correct", but that run was about generality of the claim, not correctness.
Dropping it is fine as long as nothing cites "the method wins at K = 30".

## 4. Plan and doc drift (minor)

- **M3** *(docs fixed: NOTES §5 now records route B)*: NOTES §5 records route A (indri owns `indriBlock.stan`) as chosen;
  what was built is route B (the `eigen_A` hook on stanimal's `indri-support`
  branch). Depending on it requires merging that branch into stanimal
  `master`. By the §5 cost table, M3 only matters above N of about 4000.
- **STATUS** *(fixed)*: "create the GitHub remote" is listed as open but the remote
  exists; the recorded head `51bac4d` is stale; "Two failure modes worth
  remembering" lists four.
- **asreml engine**: silently ignores `engine_args`; its no-columns fallback
  builds `cbind(...) ~ 0 +`, which is not a valid formula.
- **`.indri_diagnostics()`**: reports `converged = TRUE` when every block
  returned NA (the rstan path), so failures there go unnoticed.

## Suggested order

1. ~~Fix the sampler and its test.~~ Done.
2. ~~Run the coverage study with both draw assemblies.~~ Done: mixture wins.
   Follow-up: `corrG_draws`, then `spread = "mixture"` as the default.
3. Re-baseline the K = 30 claim against a REML full model.

Then the existing STATUS items (missing data end to end, asreml at K = 30).
