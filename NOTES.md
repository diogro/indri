# indri — design notes and implementation plan

Status: pre-implementation. Repo initialized (`main`, no commits yet); only this
file plus `README.md`/`.gitignore`. Everything below is a design decision taken
far enough to code from, not code.

Scope decision of record: **G only.**

---

## 1. What this package is

`indri` estimates the additive genetic covariance matrix **G** (K x K) from a
phenotype matrix **Y** (N x K), a fixed-effect design **X** (N x J) and a known
relationship matrix **A** (N x N), for trait counts where a single multivariate
animal model is either too slow or too poorly conditioned.

It is the child of two parents:

- **LiMMBo** (Python, `../limmbo`) — estimates G by fitting REML on many small
  trait subsets and reassembling the blocks. Cheap, assumption-light, no
  uncertainty.
- **stanimal** (R, `../stanimal`) — full Bayesian animal models in Stan, with an
  unstructured or sparse-factor G. Honest uncertainty, expensive in the number of
  traits.

`indri` takes the subsetting idea from the first and the Bayesian block fits from
the second.

Name: *indri* is a lemur — same genus as the animal-ish names the parent projects
already use, and short enough to be a valid R package name. Checked free on CRAN,
PyPI and GitHub.

---

## 2. Scope

In scope:

- `G` (K x K additive genetic covariance), on the original response scale.
- Prior choices that affect `G`.
- Diagnostics for the subsetting (pair coverage, per-block convergence).

Out of scope, for now:

- Reporting `E`, `P`, `h2` (E is estimated inside blocks as a nuisance, not
  summarised).
- Association testing, breeding values, `size_shape()`, imputation API,
  LOO/PPC/model comparison.
- Any frequentist API parity with `limmbo` (no CLI).

---

## 3. Why this must be two-stage

This is the load-bearing design argument. Worth re-reading before changing the
architecture.

LiMMBo's speedup does not come from fitting on fewer columns. It comes from
**decentralising E**: each subset gets its own small nuisance covariance and its
own `O(s^5)` optimisation, so the total is `b * O(N s^4 + s^5)` with `s` small
and fixed.

The tempting alternative — keep one global `E`, but evaluate the likelihood on
subsets (a composite likelihood) — buys nothing. Covering all trait pairs forces
`b ~ K^2 / s^2`, so

```
sum_j O(N s^2 + s^3)  ~  O(N K^2)
```

which is the cost of the *full* model, with a noisier objective and the same
number of parameters. Same price, worse estimator.

Therefore a faithful Bayesian translation must also decentralise `E`, which
forces a two-stage structure:

```
stage 1   b independent small Bayesian fits      ->  G_j draws (s x s)
stage 2   assemble those blocks into one K x K   ->  G, G_draws
```

Stage 2 is where LiMMBo's averaging and PSD projection live. Stage 1 is where
stanimal works.

### What the output is, and what it is not

`G` (the point estimate) is a least-squares projection of block means onto the
PSD cone — the same estimand LiMMBo targets.

`G_draws` is an **ensemble spread, not a joint posterior**. Blocks share traits,
so block posteriors are dependent; averaging draw `d` across overlapping blocks
is a heuristic that carries block-level uncertainty through, not a posterior over
the joint matrix. Document this in the man page and in anything published.
Measuring what it actually covers (section 9) is the honest way to justify it.

---

## 4. Architecture

```
Y, X, A
  |
  +-> indri_subsets(K, s, c)          b x s index matrix, every pair covered >= c
  |
  +-> block fits (parallel over blocks)
  |     stanimal::lmm_animal(Y[, I_j], X, A)   ->  G_j draws, Rhat, divergences
  |
  +-> assemble_G()
        inflate + average overlaps        -> G_bar      (K x K)
        Cholesky LS projection to PSD     -> G_hat      (the reported estimate)
        per-draw inflate + average        -> G_draws    (K x K x D, spread)
```

Three properties that make the assembly valid:

1. **Blocks are on a common scale.** stanimal standardises per trait
   (`Y_std[n,k] = (Y[n,k] - y_mean[k]) / y_sd[k]`, computed from that trait's
   observed values only) and back-transforms with
   `G = quad_form_diag(G_std, y_sd)`, i.e. `G_ij = y_sd_i * G_std_ij * y_sd_j`.
   `y_sd[k]` depends only on trait k's observed data, so **it is identical in
   every block that contains trait k**. Assembled entries are therefore on the
   same original scale regardless of which block produced them. (Verified in
   `../stanimal/package/inst/stan/animalModelG.stan` lines 55-93, 233-234.)
   Pin this with a test in M1 — it is the assumption the whole method rests on.

2. **A's eigen-decomposition is shared.** Every block needs the same
   `Q diag(lambda) Q'`. See section 5; this is the one place a custom Stan file
   is genuinely required.

3. **The projection is monotone, the averaging is linear.** The mean of the
   per-draw averages equals the average of the block means, so the point estimate
   is the same whether you assemble means or draws first. Project the **mean**
   (nonlinear, must not be applied draw-by-draw if you care about smoothness of
   the spread) and report spread from the unprojected draw ensemble.

---

## 5. The Stan delta (verified against stanimal)

### Facts established by reading the package and by experiment

1. **`A` is used only for its eigen-decomposition.** In
   `../stanimal/package/inst/stan/animalModelG.stan` it appears only in
   `transformed data`:

   ```stan
   matrix[N, N] Q = eigenvectors_sym(A);
   vector[N] ev = eigenvalues_sym(A);
   for (n in 1:N) lambda[n] = fmax(ev[n], 0);
   ...
   Ys_observed = Q' * Y_std;
   Xs = Q' * X_std;
   ```

   `model` and `generated quantities` never touch `A`. The patch is a data-block
   swap, not a fork.

2. **The cost saved is `O(N^3)` per chain per block fit.** `transformed data` runs
   once per chain, not per iteration, so the overhead is `4 * b` decompositions
   per run:

   | N | per chain | x 4 chains x 30 blocks | vs. sampling time |
   |---|---|---|---|
   | 1,000 | ~2 s | ~4 min | ~1-3% |
   | 4,000 | ~60 s | ~2 h | ~50% |
   | 10,000 | ~10 min | ~20 h | dominates |

   i.e. negligible in the package's stated niche (hundreds to a few thousand
   individuals), essential above it. Do not block M1 on this.

3. **`model_file=` cannot be used to inject `Q`/`lambda`.** `lmm_animal()` builds
   a fixed `stan_data` list containing `A = A` and passes it to
   `cstan(data = stan_data, ...)`. `...` is forwarded to the sampler, so passing
   `data =` there collides with the explicit `data =` ("formal argument matched
   by multiple actual arguments"). There is no hook.

4. **Extra data keys are ignored; missing declared keys are fatal.** Tested with a
   toy model: `cmdstanr` accepted `data = list(N = 10, y = ..., junk_not_in_model = 42)`
   and sampled normally. The reverse — a model declaring `Q` that is never passed
   it — fails at startup. So `model_file=` can replace the *model* but never its
   *inputs*, which rules it out.

5. **`model_file` is already incompatible with the factor structure.**
   `lmm_animal()` stops with "`model_file` cannot be combined with
   `genetic_structure = "factor"`" (R/lmm_animal.R:251-252). If indri wants
   horseshoe blocks later, it cannot go through that path anyway.

6. **Conventions that must match exactly when precomputing in R:**
   - Stan's `eigenvectors_sym` returns eigenvalues **ascending**; R's `eigen()`
     returns them **descending**. Reverse `values` and the columns of `vectors`.
   - stanimal floors at zero: `pmax(values, 0)`.

### Routes considered

| Route | Verdict |
|---|---|
| **A. indri owns the Stan file and the sampler call** | **Chosen.** No upstream change, works at any N. Cost: indri takes on the R-side data prep (validation, sample matching, missing-data mask, `x_mean_weight`). |
| B. Add a `q_lambda = list(Q=, lambda=)` hook to stanimal | Cleanest for indri, but needs an upstream PR. Keep as an option; do not depend on it. |
| C. Pass `A = diag(sort(eigenvalues))` so Stan's decomposition is free | **Rejected.** Requires pre-rotating Y/X, but `Y_std` is computed *before* rotation on the individual axis, so pre-rotation changes per-trait standardisation; and the missing-data mask is per (individual, trait), which rotation destroys. |

### The patch

`inst/stan/indriBlock.stan` = copy of `animalModelG.stan` with:

```stan
// data
- matrix[N, N] A;              // known relationship matrix
+ matrix[N, N] Q;              // eigenvectors of A, ascending eigenvalues
+ vector[N] lambda;            // eigenvalues of A, floored at 0

// transformed data
- vector[N] lambda;
- matrix[N, N] Q = eigenvectors_sym(A);
- vector[N] ev = eigenvalues_sym(A);
- for (n in 1:N) lambda[n] = fmax(ev[n], 0);
```

Everything else — standardisation, gap parameters, likelihood, generated
quantities — stays byte-identical. Start with the unstructured variance
parametrisation only; add `H2`/factor variants later if needed.

R side, computed **once per run**:

```r
qr <- eigen(A, symmetric = TRUE)
Q  <- qr$vectors[, rev(seq_len(ncol(qr$vectors))), drop = FALSE]
lambda <- pmax(rev(qr$values), 0)
```

---

## 6. API

```r
fit_G(Y, X = NULL, A,
      s = NULL,                     # subset size; default floor(sqrt(K)) clamped
      min_cooccurrence = 3,         # c
      subsets = NULL,               # precomputed b x s matrix, for reproducibility
      block = "stanimal",           # "stanimal" | "stanimal_qr" (M2)
      assembly = "ls",              # "ls" | "average"
      block_args = list(),          # forwarded: iter, chains, cores, priors
      seed = NULL)
```

Returns an `indri_fit`:

| element | contents |
|---|---|
| `G` | K x K point estimate (LS projection of the averaged blocks) |
| `corrG`, `sigma_G` | derived summaries |
| `G_draws` | K x K x D assembled ensemble (spread, not a posterior) |
| `G_bar` | pre-projection average, for diagnosis |
| `subsets` | b x s index matrix |
| `coverage` | K x K pair counts, plus which pairs are below threshold |
| `blocks` | the raw block fits (E included, for debugging only) |
| `diagnostics` | per-block Rhat / divergences / ESS, projection status |

Exported internals, each independently testable:
`indri_subsets()`, `assemble_G()`, `project_psd()`.

Methods: `print()`, `summary()`, `plot()` (correlation heatmap with coverage
overlay).

---

## 7. Package skeleton

```
package/
  DESCRIPTION        Imports: stanimal, posterior, cmdstanr | rstan;
                     Suggests: testthat, loo
  NAMESPACE
  R/fit_G.R          entry point and orchestration
  R/subsets.R        set-cover sampler (port of limmbo's multiple_set_covers_all)
  R/assemble.R       inflate / average / draw-wise assembly
  R/project_psd.R    Cholesky-parametrised LS fit (limmbo's BFGS objective)
  R/blocks.R         block engine abstraction ("stanimal", later "stanimal_qr")
  R/methods.R        print / summary / plot
  inst/stan/indriBlock.stan        (M2)
  tests/testthat/    test-subsets, test-assemble, test-project, test-fit-G,
                     test-scale, test-qr-equivalence
examples/            simulation_recovery.R, ratones_G.R
```

Reference implementation to port for `subsets.R` and `project_psd.R`:
`../limmbo/limmbo/utils/utils.py` (`multiple_set_covers_all`, `inflate_matrix`)
and `../limmbo/limmbo/core/vdbootstrap.py` (`__getBootstrapResults`,
`__fit_bootstrap_results`, `__rss`, `__RSS_grad_compute`).

---

## 8. Milestones

| # | Deliverable | Acceptance test |
|---|---|---|
| **M0** | Skeleton + `indri_subsets()` | Reproduces limmbo's `multiple_set_covers_all` for a grid of `(K, s, c)` at the same seed; every pair covered `>= c`; errors cleanly on `s < 2`, `c < 1`, `s > K` |
| **M1** | `fit_G()` on stock `stanimal::lmm_animal()`, average + LS projection | (a) a 2-trait block equals the full 2-trait `lmm_animal()` within MC error; (b) on a 12-trait simulation with known G, median per-entry relative error (correlation scale) within 2x the full 12-trait stanimal fit; (c) entries estimated by different blocks agree within MC error — **the scale-consistency test**; (d) runs end-to-end on the shipped `ratones` data |
| **M2** | `indriBlock.stan` + `block = "stanimal_qr"` | **Deterministic** equivalence: same data, fixed parameter vector, `log_prob()` matches the stock model to machine precision. Not a posterior comparison |
| **M3** | Uncertainty + diagnostics | Per-block Rhat/divergence table, pair-coverage report, `G_draws`; simulation study of interval coverage, documented in the man page |
| **M4** | Performance | Block caching keyed by (subset, hyperparameters); parallel schedule over blocks; target wall clock <= 2x the equivalent number of independent `lmm_animal()` fits |

M1 is the number that decides whether the project continues. If assembled G is
not within ~2x the full model's error, the estimator is not worth the machinery.

---

## 9. Validation plan

- **Simulation**: `stanimal::random_covariance(method = "factor")` for known G;
  sweep K in {8, 12, 20}, N in {300, 1000}, s in {3, 4, 5}, c in {3, 5}.
- **Metrics**: relative Frobenius error on the correlation scale (never absolute —
  entries can be ~1e-3), per-entry error against pair coverage, correlation of
  assembled `corrG` with truth, wall clock, block-level Rhat.
- **Reference points**: full stanimal fit (K small enough), and LiMMBo's own
  estimate where the Python stack is available.
- **Interval coverage** (M3): fraction of true entries inside the 95% band from
  `G_draws`, reported per coverage count. Expect this to be optimistic; the
  number is the deliverable, not a pass/fail gate.
- **Do not** claim posterior calibration. Claim what is measured.

---

## 10. Open decisions

1. **Block engine.** MVP: `stanimal::lmm_animal()` (pure R, handles gaps, slower
   per block). Later: a REML engine matching the paper. *Default: stanimal only.*
2. **Priors inside blocks.** stanimal's LKJ default vs
   `genetic_structure = "factor"` (horseshoe). This changes *what* G is
   estimated, not just speed — a factor prior shrinks blocks toward low rank and
   the assembly inherits it. Note that the factor structure is incompatible with
   `model_file` today, so a factor path implies owning the Stan file from M1
   onward. *Default: LKJ, factor exposed as an option later.*
3. **Assembly.** LS projection (matches LiMMBo, guarantees PSD) vs plain
   averaging (`Cg_average` in the paper). *Default: both computed, LS reported.*
4. **Uncertainty semantics.** Draw-wise ensemble vs point-only vs full block
   bootstrap. *Default: draw-wise, documented as a spread.*
5. **Dependency on stanimal**: `Imports:` (clean, version-coupled) vs vendoring
   its Stan files (independent, forks the priors). *Default: `Imports:` for the
   stock engine, vendored file only for `indriBlock.stan`.*
6. **License**: GPL-3, to match stanimal.
7. **Remote**: `gh repo create diogro/indri --private` — not done; needs a
   decision on public/private and org.

---

## 11. Risks

- **Coverage-starved pairs.** With small `s`, some pairs get few estimates and
  their entries are the noisy ones. The `coverage` field exists so this is
  visible rather than silent; consider refusing to fit when any pair falls below
  `c`.
- **Cost is `b` MCMC fits.** At K = 40, s = 5, expect tens of fits. Parallelism
  over blocks is a day-one requirement, not an optimization.
- **The uncertainty story is the weak point** and the one a reviewer will attack.
  Measure it (M3) rather than defending it.
- **Missing data interacts with subsetting.** A gap-heavy trait is weakly
  identified in a 4-trait block. Either document, or default to complete-case
  rows per block. `residual_structure = "diagonal"` with many gaps is explicitly
  warned against in stanimal's own docs.
- **Vendoring drift.** `indriBlock.stan` will diverge from `animalModelG.stan` as
  stanimal evolves. The deterministic `log_prob()` test in M2 is the guard; re-run
  it whenever the stanimal version changes.

---

## 12. Environment (verified on this machine)

- R 4.6.0; `cmdstanr`, `rstan`, `posterior`, `testthat`, `loo` installed.
- CmdStan 2.27.0 / 2.28.1 / 2.31.0 available.
- Python/LiMMBo reference implementation in `../limmbo`, runnable through pixi
  (`.pixi/envs/default`, Python 3.6 + limix 1.0.18) — this is the fastest way to
  generate reference subsets and reference G estimates.
- Repo: `/Users/diogro/Documents/GitHub/indri`, branch `main`, no commits.

---

## 13. References

- LiMMBo paper: <https://www.biorxiv.org/content/10.1101/255497v1>
- LiMMBo source: `../limmbo`
- stanimal source and interface guide: `../stanimal`,
  `../stanimal/AGENTS.md`
- Regularized horseshoe: Piironen & Vehtari (2017)
- LKJ prior: Lewandowski, Kurowicka & Joe (2009)
