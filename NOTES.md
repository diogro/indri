# indri — design notes and implementation plan

Status: the package is loadable and `R CMD check`-clean, with the subset sampler,
the assembly, the PSD projection (analytical gradient), the engine registry and
all three engines — `stanimal` (posterior draws), `sommer` and `asreml` (REML).
160+ test assertions pass; the asreml tests skip where no licence is available.

`examples/six_traits.R` runs the whole pipeline on a known 6 x 6 covariance
(about three minutes) and compares it entry by entry with the truth.

stanimal 0.8.0 takes the eigendecomposition of `A` as model data
(`lmm_animal(eigen_A =)`); indri decomposes `A` once per `fit_G()` call and the
`stanimal` engine passes it to every block (§5).

Scope decision of record: **G only**, with **pluggable block engines**.

For the current state of the work, what has been measured and what to do
next, read [STATUS.md](STATUS.md). This file is the design.

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

`indri` takes the subsetting idea from the first and the block fits from the
second — but the block fit is a **plugin**, not a hard dependency. LiMMBo calls
limix for its blocks; indri should be able to call **stanimal**, **sommer** or
**ASReml-R**, because the right engine depends on N, on whether the data have
gaps, on whether the machine has an ASReml licence, and on whether the user wants
a posterior or a likelihood.

Name: *indri* is a lemur — same genus as the animal-ish names the parent projects
already use, and short enough to be a valid R package name. Checked free on CRAN,
PyPI and GitHub.

---

## 2. Scope

In scope:

- `G` (K x K additive genetic covariance), on a documented response scale.
- **Pluggable block engines**, at least three:
  - `"stanimal"` — Bayesian, NUTS, returns posterior draws. Default.
  - `"sommer"` — REML in pure R, CRAN-installable, returns point estimates plus
    delta-method standard errors.
  - `"asreml"` — REML, licensed, sparse inverse of A, fastest per block.
- Prior/engine settings that affect `G`, and the diagnostics needed to compare
  them (pair coverage, per-block convergence, cross-engine agreement).
- Entry weights in the assembly (equal vs inverse-variance).

Out of scope, for now:

- Reporting `E`, `P`, `h2` (E is estimated inside blocks as a nuisance, not
  summarised).
- Association testing, breeding values, `size_shape()`, imputation API,
  LOO/PPC/model comparison.
- Any frequentist CLI parity with `limmbo`.

---

## 3. Why this must be two-stage

This is the load-bearing design argument, and it is independent of which engine
fits the blocks. Worth re-reading before changing the architecture.

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

Therefore a faithful translation must also decentralise `E`, which forces a
two-stage structure:

```
stage 1   b independent small fits (any engine)   ->  G_j, E_j, + uncertainty
stage 2   assemble those blocks into one K x K    ->  G, spread
```

Stage 2 is where LiMMBo's averaging and PSD projection live. Stage 1 is where the
engine works. The engine is interchangeable; stage 2 is not.

### What the output is, and what it is not

`G` (the point estimate) is a least-squares projection of block estimates onto
the PSD cone — the same estimand LiMMBo targets.

The spread (`G_draws` for `"stanimal"`, `G_se`/`G_draws` for the REML engines) is
an **assembled spread, not a joint posterior**. Blocks share traits, so block
estimates are dependent; combining them is a heuristic that carries block-level
uncertainty through, not a posterior over the joint matrix. Document this in the
man page and in anything published. Measuring what it actually covers
(section 10) is the honest way to justify it.

How the block draws are combined matters more than anything else here. Blocks
are fitted to the *same* individuals, so their estimation errors are nearly
perfectly correlated, while their chains are independent. Averaging draw `d`
across the `c` blocks that cover an entry therefore divides the spread by about
`sqrt(c)` without removing any real uncertainty — measured at 40% coverage for a
nominal 95% band. The **mixture** assembly (each entry of draw `d` taken from
one randomly chosen covering block) keeps the block spread and measured 89%,
against 80% for the full model (§10). Mixed draws are pointwise: entries are
drawn independently, so a mixed matrix need not be PSD, and they support
per-entry intervals, not functionals of the whole matrix.

---

## 4. Architecture

```
Y, X, A
  |
  +-> indri_subsets(K, s, c)          b x s index matrix, every pair covered >= c
  |
  +-> block fits (parallel over blocks), one engine per fit
  |     engine$fit(Y[, I_j], X, A)    ->  G_j, E_j, draws | se
  |
  +-> assemble_G()
        inflate + weighted average        -> G_bar      (K x K)
        Cholesky LS projection to PSD     -> G_hat      (the reported estimate)
        per-draw / per-se propagation     -> G_draws / G_se
```

Four properties that make the assembly valid:

1. **Blocks are on a common scale — by construction in indri.** Each engine
   treats scale differently (`stanimal` standardises internally and back-transforms
   with `G = quad_form_diag(G_std, y_sd)`; `sommer` and `ASReml` do not
   standardise at all). So indri standardises **once, globally**, before
   subsetting: `Y_std[, k] = (Y[, k] - mean_k) / sd_k` using all observed values
   of trait k, and keeps `y_sd` on the object. Every block, from every engine,
   then sees the same numbers and estimates the same G on the same scale; the
   reported `G` is on the standardised scale and can be back-transformed by the
   caller with `quad_form_diag(G, y_sd)`. This replaces an earlier design that
   relied on stanimal's per-block standardisation — that worked, but it made
   the estimator's validity depend on an internal of one engine.

2. **A is prepared once and shared.** `stanimal` wants dense `A`; `sommer` wants
   `solve(A)` marked as an inverse; `ASReml` wants a sparse three-column inverse.
   Compute each form once per run, not per block.

3. **`A`'s eigen-decomposition is shared** for the Stan engine. See section 5;
   this is the one place a custom Stan file is genuinely required. The REML
   engines solve their own systems and are unaffected.

   `cores` is a *budget* for the whole fit, not a worker count. Each engine
   brings its own parallelism — stanimal defaults to four chains running on four
   cores — so the budget is split between blocks and chains:

   ```
   per_block_cores    = engine_args$cores  %||%  min(chains, cores)
   blocks_in_parallel = max(1, floor(cores / per_block_cores))
   ```

   `chains` is never changed, because it is a statistical setting; only the
   number of processes running at once changes. `cores = 4` with the default
   engine arguments is one block with four chains in parallel; with
   `engine_args = list(chains = 1)` it is four blocks at once.

4. **The projection is monotone, the averaging is linear.** The mean of the
   per-draw averages equals the average of the block means, so the point estimate
   is the same whether you assemble draws or means first. Project the **mean**
   (nonlinear; do not project draw-by-draw if you care about smoothness of the
   spread) and report spread from the unprojected ensemble. That holds for the
   *average* of the draws; the spread itself must come from a mixture, not an
   average (section 3).

---

## 5. The Stan delta (verified against stanimal)

Only relevant to the `"stanimal"` engine.

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

   Negligible in the niche, essential above it. Do not block M1 on this. (The
   REML engines have their own N-side cost; ASReml with a sparse inverse of A is
   the only engine here that avoids the dense `O(N^3)`/`O(N^2)` entirely.)

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
   `genetic_structure = "factor"`" (`../stanimal/package/R/lmm_animal.R:251-252`).
   If indri wants horseshoe blocks later, it cannot go through that path anyway.

6. **Conventions that must match exactly when precomputing in R:**
   - Stan's `eigenvectors_sym` returns eigenvalues **ascending**; R's `eigen()`
     returns them **descending**. Reverse `values` and the columns of `vectors`.
   - stanimal floors at zero: `pmax(values, 0)`.

### Routes considered

| Route | Verdict |
|---|---|
| **A. indri owns the Stan file and the sampler call** | Originally chosen, **not what was built** (see Prototype status). No upstream change, works at any N. Cost: indri takes on the R-side data prep for its own model variant (validation, sample matching, missing-data mask, `x_mean_weight`). |
| B. Add a `q_lambda = list(Q=, lambda=)` hook to stanimal | Cleanest for indri, but needs an upstream PR. Keep as an option; do not depend on it. |
| C. Pass `A = diag(sort(eigenvalues))` so Stan's decomposition is free | **Rejected.** Requires pre-rotating Y/X, but `Y_std` is computed *before* rotation on the individual axis, so pre-rotation changes standardisation; and the missing-data mask is per (individual, trait), which rotation destroys. |

### Prototype status

What was actually built is route B, and generalised: on stanimal's
`indri-support` branch, now stanimal 0.8.0, *every* animal model takes `Q` and `lambda`
as data, `lmm_animal()` decomposes `A` once per fit in R (ascending, Stan's
order, so fits are bit-identical to 0.7.4), and `eigen_A` accepts a
precomputed pair for any structure or parametrization. The earlier duplicated
`animalModelGQR.stan` is gone, which removes the vendoring-drift risk (§12).
indri builds the pair once per `fit_G()` call (`indri_a_forms(need = "eigen")`)
and the stanimal engine passes it when the installed stanimal accepts
`eigen_A`. Measured saving: ~60 s per block at N = 4000, ~7 s at N = 2000,
negligible below 1000 (STATUS.md). Released as stanimal 0.8.0 (`master`
`9e0d831`).

The paragraphs below describe the first prototype and are kept for the record.

`lmm_animal(eigen_A = list(Q=, lambda=))` is implemented on the `indri-support`
branch of `../stanimal`, together with `animalModelGQR.stan`. Validated by
reconstructing `A` from the pair, and by asserting that the two models produce
**bit-identical draws** at the same seed and start — the trajectories can only
agree if the log density and gradient do. One thing the tests settled: the
eigenvalues and eigenvectors need only be consistently *paired*, not sorted, so
indri does not have to reverse R's descending order for correctness.

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

## 6. Block engines

### The contract

Every engine implements one function with one signature:

```r
engine$fit(Ysub, X, A_forms, control) -> list(
  G       = s x s numeric,           # required
  E       = s x s numeric,           # required (nuisance)
  draws   = s x s x D array | NULL,  # at least one of draws/se is required
  se      = s x s numeric | NULL,
  converged = logical,
  loglik  = numeric | NA,
  time    = numeric,
  warnings = character()
)
```

`A_forms` is the prepared list described in section 4 (`A`, `A_inv`,
`A_sparse_inv`); engines declare which they need so indri builds only what is
used. Each engine also declares capabilities:

| engine | estimator | uncertainty | A form | missing traits | cost/block | licence |
|---|---|---|---|---|---|---|
| `stanimal` | Bayesian posterior mean (NUTS) | `draws` | dense `A` | yes (gap parameters) | highest | GPL-3, CRAN |
| `sommer` | REML (`mmes`/`mmer`) | `se` (delta method) | dense `A^-1` | no — drops incomplete rows | low | GPL, CRAN |
| `asreml` | REML (AI) | `se` (from `$vcov`) | sparse `A^-1` | no | lowest | commercial |

### Rules the engines impose on the design

1. **One engine per fit.** Recorded on the object (`attr(fit, "engine")`). Mixing
   engines across blocks is an error, not a warning: the estimands differ (REML
   point estimate vs posterior mean under a prior), so mixed blocks would be
   averaged on unequal footing.
2. **The estimand depends on the engine.** With weak priors and enough N the
   three agree; with a few hundred individuals they do not, and `stanimal` will
   shrink more. Any cross-engine comparison must say which estimand it is
   comparing. This belongs in the man page, not just here.
3. **Missing data is a capability, not an assumption.** `stanimal` handles gaps
   natively; the REML engines silently drop incomplete rows. `missing =
   c("error", "complete_block")` decides the behaviour, defaulting to `"error"`
   when the engine cannot handle gaps and any gap falls in the block. Dropping
   rows per block also breaks the "same data in every block" property, so it must
   be recorded and reported.
4. **Weights.** Because every engine reports a per-entry uncertainty, the
   assembly can weight entries by inverse variance instead of LiMMBo's equal
   weights: `weights = c("none", "inverse_variance")`. Default `"none"` for
   fidelity to the paper, with `"inverse_variance"` as the obvious improvement to
   measure in M4.

### Verified call forms

Taken from `../stanimal/examples/package_comparison.R`, which fits the same
bivariate animal model with all six engines and is the reference for how each is
wired. Port these, don't reinvent them.

**sommer** (long format, `mmes`; term 1 is the animal/G term, term 2 the
residual/E term):

```r
Ainv <- solve(A); attr(Ainv, "inverse") <- TRUE
long <- data.frame(id = factor(rep(ids, K)),
                   trait = factor(rep(traits, each = N)),
                   value = as.vector(Y),
                   record = rep(seq_len(N), K), <covariates repeated per trait>)
fit <- sommer::mmes(
  value ~ 0 + trait + trait:x1,
  random = ~ vsm(usm(trait), ism(id), Gu = Ainv),
  rcov   = ~ vsm(usm(trait), ism(record)),
  data = long, control = mmes_control(diagnostics = list(verbose = FALSE)))
g  <- sommer::covmatrix_mmes(fit, 1, se = TRUE)
G  <- g$covariance
se <- g$covariance.se
```

`sommer::mmer()` is the wide-format alternative (`cbind(y1, y2) ~ ...`); `mmes`
is what the benchmark uses and is the safer default. `covmatrix_mmes()` exposes
`$covariance`, `$correlation`, `$covariance.se`, `$correlation.se`
(delta method from the fit's `theta_se`).

**ASReml-R** (wide format; bivariate fits need a second pass):

```r
Ainv <- <sparse three-column inverse of A, class "ginv">
fit <- asreml::asreml(
  fixed  = cbind(y1, y2) ~ trait + trait:x1,
  random = ~ vm(id, Ainv):us(trait, init = c(1, 0, 1)),
  rcov   = ~ units:us(trait, init = c(1, 0, 1)),
  data = dat)
fit <- update(fit)
vc <- summary(fit)$varcomp            # lower triangle, row by row
G  <- lower2mat(vc[grep("^vm\\(", rownames(vc)), 1])
```

Helpers to port from that file (`package_comparison.R:98-122`): `lower2mat()`
(lower triangle to symmetric matrix) and `sparse_ginv()` (three-column inverse
with `class "ginv"`). Note the benchmark builds the inverse from `solve(A)`
because `asreml::ainverse()` "fails intermittently on some pedigrees" — that
defeats the sparsity advantage. For the package: accept a user-supplied sparse
inverse (`nadiv::ainv(pedigree)` or `Matrix::solve(Matrix(A, sparse = TRUE))`),
fall back to the dense route with a warning.

### G extraction is engine-specific and easy to get wrong

| engine | how G comes out | format |
|---|---|---|
| `stanimal` | `fit_draws(fit)$G` / `draws_mean(fit, "G")` | K x K, original scale |
| `sommer` | `covmatrix_mmes(fit, 1)$covariance` | K x K |
| `asreml` | `summary(fit)$varcomp`, `"^vm\\("` rows, lower triangle | vector -> `lower2mat()` |

Written as one paragraph per engine in `R/engines/*.R`, each with its own test
against a known simulation (section 10).

Implementation notes. `sommer` runs `mmes()` single-threaded
(`OMP_NUM_THREADS = 1`): its solver multithreads even for bivariate blocks and
oversubscribes alongside block-level parallelism. `asreml` (4.2) reports the
standard errors in `summary(fit)$varcomp$std.error` rather than `$vcov`, and its
default "sigma" parameterisation still extracts cleanly with the `"^vm\\("` and
`"^units.*trait!trait"` patterns above. The PSD projection supplies an analytical
gradient; the numerical gradient over `K (K + 1) / 2` parameters is what made
large-K fits stall, not the engines.

---

## 7. API

```r
fit_G(Y, X = NULL, A,
      s = NULL,                     # subset size; engine default (stanimal
                                    #   min(floor(sqrt(K)), 4); REML engines 2)
      min_cooccurrence = NULL,      # c; engine default (stanimal 3, REML 1)
      subsets = NULL,               # precomputed b x s matrix, for reproducibility
      engine = c("stanimal", "sommer", "asreml"),
      engine_args = list(),         # forwarded to the engine (iter, chains, priors, ...)
      assembly = "ls",              # "ls" | "average"
      weights = "none",             # "none" | "inverse_variance"
      missing = "error",            # "error" | "complete_block"
      seed = NULL)
```

Returns an `indri_fit`:

| element | contents |
|---|---|
| `G` | K x K point estimate (LS projection of the averaged blocks) |
| `corrG`, `sigma_G` | derived summaries |
| `G_draws` / `G_se` | assembled spread, per engine capability |
| `G_bar` | pre-projection average, for diagnosis |
| `engine`, `engine_args` | what ran |
| `subsets` | b x s index matrix |
| `coverage` | K x K pair counts, plus pairs below threshold |
| `y_sd` | per-trait scale, for back-transformation |
| `blocks` | raw block fits (E included, debugging only) |
| `diagnostics` | per-block convergence (Rhat/divergences, or REML convergence), projection status, dropped rows |

Exported internals, each independently testable: `indri_subsets()`,
`assemble_G()`, `project_psd()`, and the engine constructors.

Methods: `print()`, `summary()`, `plot()` (correlation heatmap with coverage
overlay), `compare_engines()` (sugar over the same data, section 10).

---

## 8. Package skeleton

```
package/
  DESCRIPTION        Imports: posterior;
                     Suggests: stanimal, cmdstanr, rstan, sommer, asreml, nadiv
                     (engines are optional: stanimal is not a hard dependency once
                      sommer/asreml exist, but one of them must be installed)
  NAMESPACE
  R/fit_G.R              entry point and orchestration
  R/subsets.R            set-cover sampler (port of limmbo's multiple_set_covers_all)
  R/assemble.R           inflate / weighted average / draw and se propagation
  R/project_psd.R        Cholesky-parametrised LS fit (limmbo's BFGS objective)
  R/a_forms.R            dense A, dense inverse, sparse inverse, built once
  R/engines.R            registry, capability checks, dispatch
  R/engines/stanimal.R   block fit via stanimal::lmm_animal()
  R/engines/sommer.R     block fit via sommer::mmes()
  R/engines/asreml.R     block fit via asreml::asreml()
  R/methods.R            print / summary / plot / compare_engines
  (no Stan file: M3 went through stanimal's `eigen_A`, see §5)
  tests/testthat/        test-subsets, test-assemble, test-project, test-scale,
                         test-engine-stanimal, test-engine-sommer,
                         test-engine-asreml, test-cross-engine, test-qr-equivalence
examples/                simulation_recovery.R, engine_comparison.R, ratones_G.R
```

Reference implementations to port: `../limmbo/limmbo/utils/utils.py`
(`multiple_set_covers_all`, `inflate_matrix`) and
`../limmbo/limmbo/core/vdbootstrap.py` (`__getBootstrapResults`,
`__fit_bootstrap_results`, `__rss`, `__RSS_grad_compute`).

---

## 9. Milestones

Engines come early: they are cheap to add, and three of them make every later
decision measurable (a bug in the assembly shows up as disagreement between
engines).

| # | Deliverable | Acceptance test |
|---|---|---|
| **M0** | Skeleton + `indri_subsets()` — **done** | Every pair covered `>= c`; valid 1-based indices; deterministic given a seed; caller's RNG state restored; errors cleanly on `s < 2`, `c < 1`, `s > K`. Note: bit-for-bit reproduction of limmbo's `multiple_set_covers_all` is **not** achievable — its permutation comes from numpy's global generator, so LiMMBo's own subsets are not reproducible run to run. indri implements the same construction with a seeded generator, which is a deliberate divergence. Second divergence: LiMMBo's construction covers each pair only `ceiling(c / 2)` times by design and misses `c` in about half of settings (always for s = 2); indri tops up with greedy subsets until the guarantee holds (see REVIEW.md §1) |
| **M1** | `fit_G()` on the `stanimal` engine + global standardisation + LS projection — **done, gate passed** | (a) a 2-trait block equals the full 2-trait `lmm_animal()` within MC error, after undoing the global scale; (b) on a 12-trait simulation with known G, median per-entry relative error (correlation scale) within 2x the full 12-trait stanimal fit; (c) entries estimated by different blocks agree within MC error — **the scale-consistency test**; (d) runs end-to-end on the shipped `ratones` data |
| **M2** | `sommer` and `asreml` engines behind the contract | Each engine, on the same 3-trait block and the same data, gives G within REML tolerance of the other and of `stanimal` (MC error); `"asreml"` errors informatively when the licence is absent; no engine is imported at package load |
| **M3** | Share the eigendecomposition of `A` across stanimal blocks — **done** | stanimal takes `Q`/`lambda` as data in every animal model; a block fitted with the shared pair is bit-identical to a direct `lmm_animal()` fit; stock stanimal fits are bit-identical to 0.7.4; saving measured (~60 s per block at N = 4000). The original acceptance test (deterministic `log_prob()` equivalence) is subsumed by the bit-identical draws |
| **M4** | Uncertainty + weights + diagnostics | **mostly done**: per-block convergence table (Rhat/divergences or REML convergence), pair-coverage report, `G_draws`; equal vs inverse-variance weights compared on simulation (a wash — equal stays default); interval coverage measured: averaged draws 0.40, mixture 0.89, full model 0.80 at nominal 0.95. Left: `corrG_draws` and making `spread = "mixture"` the default |
| **M5** | Performance | Block caching keyed by (subset, engine, hyperparameters); parallel schedule over blocks; target wall clock <= 2x the equivalent number of independent block fits with that engine |

M1 is the number that decides whether the project continues: if assembled G is
not within ~2x the full model's error, the estimator is not worth the machinery.
M2 decides which engine is the default for which N.

**M1 result** (12 traits, N = 300, s = 4, 46 six-trait subsets, min_cooccurrence
3, three seeds, `results/validation_12_traits_*.csv`):

| seed | indri median per-entry abs corr error | full 12-trait model | ratio |
|---|---|---|---|
| 1 | 0.0895 | 0.0912 | 0.98 |
| 2 | 0.0829 | 0.0964 | 0.86 |
| 3 | 0.1079 | 0.1557 | 0.69 |
| median | 0.0895 | 0.0964 | **0.93** |

The subset estimator is as accurate as fitting all twelve traits at once, and on
these seeds slightly more so, at 4.6x the wall clock. That ratio is the whole
argument for the method turning on K: at 12 traits the full model is still
cheap, and the subsetting only pays once the full model's `O(P^5)` term bites.
The next measurement worth making is the same study at K = 30-40, where the full
model becomes the slow side.

---

## 10. Validation plan

- **Simulation**: `stanimal::random_covariance(method = "factor")` for known G;
  sweep K in {8, 12, 20}, N in {300, 1000}, s in {3, 4, 5}, c in {3, 5}.
- **Metrics**: relative Frobenius error on the correlation scale (never absolute —
  entries can be ~1e-3), per-entry error against pair coverage, correlation of
  assembled `corrG` with truth, wall clock, block-level convergence.
- **Cross-engine**: the same subset, the same standardised data, three engines;
  report pairwise relative differences in G. This is the cheapest guard against a
  wiring bug in any engine and the only way to know whether engine choice changes
  the answer at a given N.
- **Reference points**: full stanimal fit (K small enough), and LiMMBo's own
  estimate where the Python stack is available.
- **Interval coverage** (M4): fraction of true entries inside the 95% band,
  reported per coverage count and per engine. Expect this to be optimistic; the
  number is the deliverable, not a pass/fail gate. *Measured
  (`examples/uncertainty_coverage.R`, 12 traits, three seeds): indri's assembled
  band covers the true correlation in ~43% of entries at the nominal 95% level
  (per seed 53%, 48%, 29%), against ~81% for the full model (85%, 82%, 77%);
  indri's mean width is ~0.17 vs ~0.46 for the full model, and its coverage
  falls as the pair's coverage count rises.* *Diagnosed and fixed (six seeds,
  396 pairs): the narrowing is the averaging (§3). On the same fits the mixture
  assembly covers 0.89 (width 0.52), the average scaled by sqrt(coverage) 0.90,
  the plain average 0.40 (width 0.16), the full model 0.80 (width 0.49); the
  mixture's coverage is flat across pair-coverage bands. Still below nominal —
  report the measured number, not "95%".*
- **Weights** (M4): same simulated data, `weights = "none"` vs
  `"inverse_variance"`, reported on median per-entry error and correlation of
  correlations. Measured on the 12-trait seed-1 benchmark: inverse-variance
  (median |err| 0.0894, q90 0.209) does not beat equal weights (0.0866, 0.213)
  — a wash within Monte Carlo noise. Equal weights stay the default.
- **Do not** claim posterior calibration. Claim what is measured.

---

## 11. Open decisions

1. **Default engine.** `"stanimal"` gives uncertainty and handles gaps;
   `"asreml"`/`"sommer"` are far faster (REML) but drop incomplete rows, and
   `"asreml"` needs a licence. *All three are implemented and cross-validated;
   the default stays `"stanimal"`, with the REML engines for the N-heavy case.*
2. **Priors inside blocks.** stanimal's LKJ default vs
   `genetic_structure = "factor"` (horseshoe). This changes *what* G is
   estimated, not just speed. The factor structure is incompatible with
   `model_file` today, so a factor path implies owning the Stan file from M1
   onward. *Default: LKJ; factor later, and only for the Stan engine.*
3. **Assembly.** LS projection (matches LiMMBo, guarantees PSD) vs plain
   averaging (`Cg_average` in the paper). *Default: both computed, LS reported.*
4. **Weights.** Equal (paper fidelity) vs inverse-variance. Inverse-variance now
   weights the projection objective as well as the average
   (`sum_j || W_j * (C[I_j,I_j] - G_j) ||^2`), and the stanimal engine's weight
   is the per-entry posterior precision from its draws (REML engines keep their
   `se`). *Default: equal — measured as a wash at N = 300 (see §10), so the
   weighting remains an option rather than a recommendation.*
5. **Uncertainty per engine.** `stanimal` -> draws; `sommer`/`asreml` -> delta
   method SEs propagated as a Gaussian ensemble. The independence caveat in
   propagating a per-entry `se` matrix (entries are correlated) must be stated.
   *Default: propagate as documented, or `"none"` if the user prefers point
   estimates. For draws, assemble by mixture, not by average: averaged draws
   cover 0.40 at nominal 0.95, mixed draws 0.89 (full model 0.80); see §3 and
   §10. Present mixed bands as pointwise intervals with their measured
   coverage. The REML engines' Gaussian ensemble has the same averaging trap
   and is unmeasured.*
6. **Missing data.** `"error"` vs `"complete_block"` default. *Default: `"error"`.*
6b. **Non-converged REML blocks.** *Decided (2026-11-07, STATUS.md): reported
    in `diagnostics`, not refitted or dropped; they are uncertified optima,
    not bad estimates, and no refit tried (iterations, starts, larger block)
    improved them.*
6a. **Block size per engine.** *Decided (2026-11-07, STATUS.md), revised
    after the sommer ordering fix: REML engines use s = 3 (4-trait blocks)
    with c = 3 distinct contexts per pair; stanimal keeps the sqrt(K) rule.
    Subsets are never repeated, and c is capped by the block width. (The
    first decision, pairs once each, rested on a sommer sweep corrupted by
    the ordering bug.) asreml defaults to two workers, one licence seat each.*
6c. **REML uncertainty.** *Decided: a joint normal approximation per block
    from the REML sampling covariance, assembled by mixture; measured 0.96-0.98
    coverage at nominal 0.95.*
7. **stanimal as dependency.** `Imports:` for the default engine vs `Suggests:`
   with a runtime check. Now that three engines exist, `Suggests:` is more
   defensible. *Default: `Suggests:`, checked at call time.*
8. **Licence**: GPL-3, to match stanimal.
9. **Remote**: *done* — `diogro/indri`, private.
10. **`spread` default.** `"average"` (today) vs `"mixture"`. *Recommended:
    `"mixture"`, once `corrG_draws` exists, since that is what was measured.*

---

## 12. Risks

- **Engine heterogeneity.** Three estimands (posterior mean, REML, AI-REML),
  three uncertainty formats, three A representations. The contract and the
  cross-engine test exist to contain this; the failure mode is a plausible-looking
  G that is subtly on the wrong scale or from the wrong term.
- **Coverage-starved pairs.** With small `s`, some pairs get few estimates and
  their entries are the noisy ones. The `coverage` field exists so this is
  visible rather than silent; consider refusing to fit when any pair falls below
  `c`.
- **Cost is `b` block fits.** At K = 40, s = 5, expect tens of fits. With a REML
  engine this is minutes; with `stanimal` it is hours. Parallelism over blocks is
  a day-one requirement, and it is implemented (see below).
- **Do not fork R around cmdstan.** `parallel::mclapply()` looks like the obvious
  way to run blocks concurrently and silently destroys runs: the forked children
  conflict with the `processx` supervisor cmdstanr uses, the run dies mid-fit
  with no R error, and the only clue is that the log stops between chain
  messages. PSOCK workers (`parallel::makeCluster`, separate processes) run the
  same 23-block job in 0.8 min against 2.6 min sequentially, with identical
  estimates.
- **Licensed dependency.** ASReml-R cannot be in CI. Its tests must skip
  gracefully, and the engine must fail with an actionable message.
- **The uncertainty story is the weak point** and the one a reviewer will attack.
  Measure it (M4) rather than defending it. Measured: the mixture band reaches
  0.89 at nominal 0.95. The remaining gap is shared with the full model, which
  points at prior shrinkage and short chains rather than the assembly.
- **Missing data interacts with subsetting**, differently per engine. A gap-heavy
  trait is weakly identified in a 4-trait block under `stanimal`, and silently
  loses rows under the REML engines.
- **Vendoring drift.** *Retired*: indri ships no Stan file, and stanimal no
  longer carries a duplicated model for `eigen_A` (§5).

---

## 13. Environment (verified on this machine)

- R 4.6.0; `cmdstanr`, `rstan`, `posterior`, `testthat`, `loo` installed.
- CmdStan 2.27.0 / 2.28.1 / 2.31.0 available.
- Block engines available: `stanimal` (source at `../stanimal`; install with
  `devtools::install_local("../stanimal/package")` if not installed),
  `sommer` 4.5.12, `asreml` 4.2.0.482 (licensed; licence checks out),
  `nadiv` 2.18.0. `MCMCglmm` 2.36 also present but not an indri engine.
- Python/LiMMBo reference implementation in `../limmbo`, runnable through pixi
  (`.pixi/envs/default`, Python 3.6 + limix 1.0.18) — fastest way to generate
  reference subsets and reference G estimates.
- Repo: `/Users/diogro/Documents/GitHub/indri`, branch `main`.
- **Long CPU-bound runs get killed silently in this environment** after roughly
  15 CPU-minutes, with exit code 1 and nothing in the log between the last chain
  messages. A single-threaded run therefore survives about 15 minutes, and a run
  using six workers survives about two and a half. The validation study is
  chunked one seed per invocation for this reason, not because of a bug in
  indri; each seed costs about 8 CPU-minutes.

---

## 14. References

- LiMMBo paper: <https://www.biorxiv.org/content/10.1101/255497v1>
- LiMMBo source: `../limmbo`
- stanimal source and interface guide: `../stanimal`, `../stanimal/AGENTS.md`
- Engine comparison script with working call forms for all engines:
  `../stanimal/examples/package_comparison.R` (sommer ~line 212, ASReml ~line
  272, `lower2mat`/`sparse_ginv` at ~98-122)
- Regularized horseshoe: Piironen & Vehtari (2017)
- LKJ prior: Lewandowski, Kurowicka & Joe (2009)
