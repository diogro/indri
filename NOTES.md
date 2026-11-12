# indri — design notes

This is the design: what the estimator is, why it is built this way, the engine
contract and the decisions behind the defaults. What has been measured, what
changed when, and what to do next live in [STATUS.md](STATUS.md); coding rules
for the package are in [AGENTS.md](AGENTS.md).

Scope decision of record: **G only**, with **pluggable block engines**.

---

## 1. What indri is

`indri` estimates the additive genetic covariance matrix **G** (K x K) from a
phenotype matrix **Y** (N x K), a fixed-effect design **X** (N x J) and a known
relationship matrix **A** (N x N), for trait counts where one multivariate
animal model is too slow or too poorly conditioned to fit.

It combines two parents:

- **LiMMBo** (Python, `../limmbo`) fits REML on many small trait subsets and
  reassembles the blocks into one matrix. Cheap and assumption-light, but with
  no uncertainty.
- **stanimal** (R, `../stanimal`) fits full Bayesian animal models in Stan.
  Honest uncertainty, but expensive in the number of traits.

indri takes the subsetting from the first and the block fits from the second,
with the block fit as a **plugin**: `"stanimal"` (Bayesian), `"sommer"` (REML,
CRAN) or `"asreml"` (REML, licensed). The right engine depends on N, on gaps in
the data, on whether a licence is available and on whether a posterior or a
likelihood is wanted.

Name: *indri* is a lemur, in the spirit of the parent projects' names, and free
on CRAN, PyPI and GitHub.

## 2. Scope

In scope:

- `G` in the units of `Y`, with `corrG` and per-entry uncertainty from every
  engine.
- Pluggable block engines behind one contract, with the diagnostics needed to
  compare them (pair coverage, per-block convergence, cross-engine agreement).
- Assembly options: least-squares PSD projection or plain average, equal or
  inverse-variance weights.

- `E`, assembled from the blocks like `G` (point estimate only, kept
  positive definite), because the association scans need it.
- Association scans with the assembled `G` and `E` (`indri_gwas()`):
  multi- and single-trait, mixed or linear model, optional dominance. The
  design and its decisions are in [GWAS_PLAN.md](GWAS_PLAN.md).

Out of scope for now: reporting `P` or `h2`, uncertainty on `E`, breeding
values, imputation, model comparison, and any CLI parity with `limmbo`.

## 3. Why the estimator must be two-stage

The load-bearing argument, independent of the engine; re-read it before
changing the architecture.

LiMMBo's speed does not come from fitting fewer columns. It comes from
**decentralising E**: each subset gets its own small residual covariance and its
own small optimisation, so the cost is `b * O(N s^4 + s^5)` with `s` small and
fixed.

The tempting alternative, one global `E` with the likelihood evaluated on
subsets (a composite likelihood), buys nothing. Covering every pair needs
`b ~ K^2 / s^2` subsets, so the total is `sum_j O(N s^2 + s^3) ~ O(N K^2)`: the
cost of the full model, with a noisier objective and the same parameters.

So a faithful translation decentralises `E` too, which forces two stages:

```
stage 1   b independent small fits (any engine)   ->  G_j, E_j, + uncertainty
stage 2   assemble the blocks into one K x K       ->  G, corrG, draws
```

Stage 2 (averaging, PSD projection, draw assembly) is indri's; stage 1 is the
engine's. The engine is interchangeable; stage 2 is not.

**What the output is.** `G` is a least-squares PSD fit to the block estimates,
the estimand LiMMBo targets. The draws (`G_draws`, `corrG_draws`) are an
*assembled* spread, not a posterior over the joint matrix: blocks share traits
and data, so they are dependent. Their coverage is measured, not assumed (§8).

## 4. The pipeline

```
Y, X, A
  |  standardise each trait once (all observed values)      -> Y_std, y_sd
  |  prepare A once: dense, inverse, sparse inverse, eigen  -> a_forms
  +-> indri_subsets(K, s, c)        distinct subsets, every pair in >= c of them
  +-> block fits, parallel over blocks (one engine per fit)
  |     engine$fit(Y_std[, I_j], X, a_forms)  -> G_j, E_j, draws | se + vcov
  +-> REML only: normal draws of each block's G from its vcov
  +-> assembly, on the standardised scale
  |     weighted average                    -> G_bar
  |     least-squares PSD projection        -> G
  |     mixture of block draws              -> G_draws, corrG_draws
  +-> G, G_bar, G_draws, sigma_G to the units of Y (entry i,j times y_sd_i y_sd_j)
```

Properties the design relies on:

1. **A common scale, by construction.** Engines treat scale differently
   (stanimal standardises internally; sommer and asreml do not), so indri
   standardises once, globally, before subsetting, and every engine sees the
   same numbers. Fitting, weights, averaging, the projection and every draw run
   on that scale. Only the assembled covariances are converted to the units of
   `Y` at the end, and not re-projected: the unweighted projection objective
   depends on the units. `corrG` is scale-free; raw `blocks` stay standardised.
2. **A is prepared once.** Each engine declares the form it needs (`dense`,
   `inverse`, `sparse_inverse`, `eigen`), and each form is built once per call,
   not per block. The stanimal engine shares one eigendecomposition of `A`
   across all blocks (§6).
3. **`cores` is a budget, not a worker count.** A stanimal block brings its own
   chains, so `per_block_cores = engine_args$cores %||% min(chains, cores)` and
   `blocks_in_parallel = floor(cores / per_block_cores)`. `chains` is a
   statistical setting and is never changed. Blocks run in PSOCK workers.
4. **Average linearly, project once.** The point estimate projects the average
   of the block means; draws are never projected one by one.

## 5. Subsets

`indri_subsets(K, s, min_cooccurrence)` is LiMMBo's set-cover construction
(`multiple_set_covers_all`): trait pairs are put in a seeded random order and
covered in windows of `s - 1` pairs, so a subset holds up to `2 (s - 1)`
traits (four at `s = 3`), not `s`. Three deliberate departures from LiMMBo:

- **Seeded and reproducible**, and the caller's random state is restored
  (LiMMBo draws its order from numpy's global generator).
- **The coverage guarantee holds.** LiMMBo's passes cover each pair only
  `ceiling(c / 2)` times and leave the rest to chance; indri tops up with
  greedily built subsets until every pair reaches `c`.
- **Distinct contexts.** A subset is never repeated, and `c` counts distinct
  subsets: a pair estimated twice in the same trait context adds nothing (for
  a REML engine it is the same number). `c` is capped at the number of
  distinct subsets of that width containing a pair (one, for pairs).

Defaults: `s = min(floor(sqrt(K)), 3)`, at least 2 and kept narrow enough that
no block covers every trait; `c = 3` for stanimal and `c = 1` for the REML
engines (§9).

## 6. Block engines

### The contract

```r
engine$fit(Y, X, a_forms, control) -> list(
  G, E,          # s x s, required; on the standardised scale, in column order
  draws,         # s x s x D posterior draws (stanimal) or NULL
  se,            # s x s standard errors (REML) or NULL
  vcov,          # joint sampling covariance of the lower-triangle G entries
                 #   (REML), in (1,1), (2,1), (2,2), (3,1), ... order, or NULL
  converged,     # TRUE/FALSE (REML) or max Rhat (stanimal)
  loglik, time   # where available
)
```

Each engine declares `capabilities`: `draws`, `se`, `vcov`, `a_form`,
`missing` (whether it handles gaps in `Y`), and optionally `default_s`,
`default_cooccurrence`, `default_cores` and `default_chains`. Capabilities are
decided in the main process; `fit()` only uses what `a_forms` carries (PSOCK
workers load the installed package).

| engine | estimator | uncertainty | A form | gaps in Y | cost per block |
|---|---|---|---|---|---|
| `stanimal` | posterior mean (NUTS), `lkj_prior = 1` | posterior draws | eigendecomposition (shared) | modelled | highest |
| `sommer` | REML (`mmes`, AI with EM recovery) | `vcov` from its natural-scale parameters | dense `A^-1` | missing records skipped | low |
| `asreml` | REML (AI) | `vcov` = inverse AI matrix | sparse `A^-1` | partly observed rows kept | lowest; licence seats |
| `wombat` | REML (PX-EM then AI, WOMBAT program) | `vcov` = inverse AI matrix (`AvInfoCovs`) | sparse `A^-1` (`.gin` file) | missing records skipped; unidentified residual covariances held at zero | low; Linux/Windows binary |

Rules:

1. **One engine per fit.** The estimands differ (posterior mean under a prior
   vs REML), so blocks from different engines are never averaged together.
2. **Output in the column order of `Y`.** An engine that returns matrices in
   another order (sommer sorts factor levels alphabetically) silently puts
   variances and covariances on the wrong traits.
3. **One design for every engine.** `X` must have full column rank and be
   able to represent a constant (an intercept, or the indicators of every level
   of a factor), as stanimal requires: the responses are centred, and without
   a constant the REML engines would fit different models.
   `stanimal::animal_design()` builds `Y`, `X` and `A` from a formula.
4. **Gaps are handled by every engine.** stanimal integrates them out; in
   sommer's long format a gap is one missing record (`naMethodY =
   "exclude"` drops that row only), and asreml keeps partly observed rows
   (`na.method(y = "include")`), so an individual contributes every trait it
   has. The two REML engines agree exactly on gapped data. `missing =
   "model"` (default) passes gaps through; `"complete_block"` drops incomplete
   individuals once, before subsetting, and costs information (15% of values
   missing at random left 311 of 500 individuals and 29% wider standard
   errors, against 7% when the gaps are kept).

### Engine specifics

- **stanimal** calls `lmm_animal()` with `lkj_prior = 1` unless `engine_args`
  says otherwise, and passes one eigendecomposition of `A` as `eigen_A`
  (stanimal >= 0.8.0 takes `Q` and `lambda` as data; older versions decompose
  `A` themselves). Ascending eigenvalue order keeps a block bit-identical to a
  direct `lmm_animal()` fit. Each block writes CmdStan output to its own
  directory, deleted after the draws are read.
- **sommer** fits `mmes()` in long format: `value ~ 0 + trait (+ trait:x)`,
  `random = ~ vsm(usm(trait), ism(id), Gu = Ainv)`,
  `rcov = ~ vsm(usm(trait), ism(record))`, with bare `vsm`/`usm`/`ism` bound in
  the formula environment, trait levels in column order, `OMP_NUM_THREADS = 1`,
  and `covmatrix_mmes()` for G and E. `vcov` comes from
  `attr(fit$covParNative, "vcov")`, matched to traits by parameter name.
  Starting values go through `usm(trait, theta =)`.
- **asreml** fits in wide format, `cbind(...) ~ trait (+ trait:x)`,
  `random = ~ vm(id, Ainv):us(trait, init =)`, `rcov = ~ units:us(trait)`,
  then `update()`; G and E are read from `summary(fit)$varcomp` (`"^vm\\("`,
  `"^units.*trait!trait"`), `vcov` from `fit$ai`. `engine_args` takes `init`,
  `maxit` (50) and `updates` (up to 5, only while unconverged): asreml's own
  13 iterations stopped far from the optimum once `X` had covariates. Every
  process needs a licence seat.
- **wombat** writes `data.dat` (trait number, individual, `mu`, covariables,
  response; one record per observed trait, sorted by individual then trait),
  `animal.gin` (lower triangle of `A^-1` as column <= row, after a line with
  log |A|) with `animal.codes`, and `indri.par`
  (`FIX mu`, `COV x(1)` for the columns of `X` beyond the constant,
  `RAN animal GIN`, `SUBJ animal`), runs `wombat --batch --threadsN` in a
  temporary directory per block, and reads `BestPoint` (header, then E and G
  as upper triangles), `AvInfoCovs` (information for E then G, inverted and
  reordered to vech order) and `Iterates` (its "seems to have converged"
  verdict). A covariance whose starting value is zero stays at zero, so
  identified covariances start at `0.25 / s` and residual covariances of
  traits never observed together start, and stay, at zero. WOMBAT runs under
  `ulimit -s unlimited` (the Linux build segfaults entering AI-REML at about
  20 traits with an 8 MB stack), a non-zero exit is an error even when a
  `BestPoint` from an EM round exists, and `OVER MAXTRAIT` / `OVER MAXCOV`
  lines lift the compiled limits (42 traits, 19 covariables).
  `engine_args$wombat_rank = q` fits G at reduced rank (`ANAL MUV PC`,
  `VAR animal s q`); that also needs `animal.matrix` (A itself, inverted
  from the sparse inverse), and its log L depends on log |A| once per
  component, which is why the GIN header carries the true value.
- **Single-trait blocks** are fitted by both REML engines with formulas that
  drop the `trait` factor; `init = list(G =, E =)` sets starting covariances
  for both.

## 7. Assembly and projection

- **Average**: each entry is the (optionally inverse-variance weighted) mean
  over the subsets that contain both traits. The average of overlapping blocks
  is generally not PSD.
- **Projection** (`assembly = "ls"`, default): the PSD matrix `C = L L'`
  closest to all blocks at once, `sum_j || W_j * (C[I_j, I_j] - G_j) ||^2`,
  over the free entries of `L` by L-BFGS-B with an **analytic gradient** (a
  numerical one costs `K (K + 1) / 2` objective evaluations per step).
  `assembly = "average"` keeps the average, nudged onto the PSD cone.
- **Weights**: `"inverse_variance"` weights both the average and the
  projection, by `se` (REML) or posterior variance (stanimal). Measured as a
  wash; equal weights stay the default.
- **Alternatives measured and not adopted**: Higham's nearest correlation
  matrix, eigenvalue clipping and the weighted projection; none beat the
  current projection on all of slope, accuracy, the strongest pairs and
  integration (STATUS.md).

## 8. Uncertainty

Every engine ends with draws for each block: posterior draws from stanimal,
and for the REML engines `n_draws` (1000) joint normal draws from each block's
`vcov`, drawn in the main process from a seed derived from `seed`.

**Assemble by mixture, never by average.** Blocks are fitted to the same
individuals, so their errors are nearly perfectly correlated, while their
chains are independent; averaging draw `d` over the `c` blocks covering an
entry shrinks its spread by about `sqrt(c)` (measured: 40% coverage at nominal
95%). The mixture (`spread = "mixture"`, default) takes each entry of each draw
from one randomly chosen covering block. `corrG_draws` mixes each block's own
correlation draws, so every entry is a correlation.

Mixed draws are **pointwise**: entries are drawn independently, so a draw need
not be PSD. They support per-entry intervals, not functions of the whole matrix
(eigenvalues, integration, responses to selection). Measured coverage at
nominal 95%: about 0.96-0.98 for the REML engines (conservative) and 0.95-0.96
for stanimal at the current defaults (STATUS.md). Report the measured number,
not "95%".

## 9. Defaults and why

| setting | default | evidence (STATUS.md) |
|---|---|---|
| block size `s` | `min(floor(sqrt(K)), 3)`: blocks of up to four traits | at K = 35, 6-trait stanimal blocks shrank correlations and undercovered; REML gained nothing from larger blocks |
| `min_cooccurrence` | 3 (stanimal), 1 (REML) | for REML, three distinct contexts matched one in accuracy and coverage at twice the blocks |
| stanimal `lkj_prior` | 1 (stanimal's own default is 2) | LKJ 2 shrank correlations, undercovered and understated integration more |
| `spread` | `"mixture"` | averaging is overconfident by about `sqrt(c)` |
| `weights` | `"none"` | inverse-variance was a wash |
| `assembly` | `"ls"` | guarantees a PSD `G`; no measured alternative is better on every metric |
| asreml workers | 2, or 1 if asreml is already loaded | the licence allows two processes on one machine |
| non-converged REML blocks | reported, not refitted | they are uncertified optima, as accurate as the rest; no refit helped |

A method that understates integration is preferred to one that overstates it.

## 10. API

```r
fit_G(Y, X = NULL, A, s = NULL, min_cooccurrence = NULL, subsets = NULL,
      engine = "stanimal", engine_args = list(),
      assembly = c("ls", "average"), weights = c("none", "inverse_variance"),
      missing = c("model", "complete_block"), spread = c("mixture", "average"),
      cores = NULL, n_draws = 1000L, seed = NULL,
      progress = getOption("indri.progress", interactive()))
```

| element | contents |
|---|---|
| `G`, `corrG`, `sigma_G` | the estimate, in the units of `Y` (`corrG` scale-free) |
| `G_bar` | average before projection, units of `Y` |
| `G_draws`, `corrG_draws` | assembled draws, K x K x D, pointwise |
| `subsets`, `s`, `min_cooccurrence`, `coverage`, `counts` | the design and pair coverage |
| `blocks` | raw block fits, standardised scale, named by trait |
| `diagnostics` | `converged`, `max_rhat`, `n_unconverged`, `unconverged` (blocks and traits) |
| `y_mean`, `y_sd` | the internal standardisation |
| `engine`, `engine_args`, `assembly`, `weights`, `spread`, `parallel`, `call` | what ran |

Exported helpers, each testable on its own: `indri_subsets()`,
`indri_assemble()`, `indri_assemble_draws()`, `indri_project_psd()`,
`indri_inflate()`, `indri_a_forms()`, `indri_engine()`, `indri_engines()`.
Methods: `print()`, `summary()`, `plot()`.

## 11. Package layout

```
package/R/fit_G.R            entry point: standardisation, defaults, parallel
                             blocks, REML draws, assembly, units, diagnostics,
                             progress bar
package/R/subsets.R          set-cover sampler, top-up, distinct subsets
package/R/assemble.R         inflate, weighted average, draw assembly (mixture)
package/R/project_psd.R      least-squares PSD projection, analytic gradient
package/R/a_forms.R          dense, inverse, sparse inverse and eigen forms of A
package/R/engines.R          registry, shared helpers (init checks, normal draws)
package/R/engine-*.R         stanimal, sommer, asreml
package/R/methods.R          print / summary / plot
package/tests/testthat/      the spec
examples/                    studies and worked examples (examples/README.md)
```

DESCRIPTION: `Imports: graphics, grDevices, parallel, posterior, stats, utils`;
`Suggests: cmdstanr, rstan, sommer, asreml, nadiv, stanimal, testthat`. No
engine is loaded at package load; each is checked at call time.

## 12. Decisions

Decided (with the evidence in STATUS.md):

1. **Default engine: stanimal.** It gives a posterior; the REML
   engines are the fast choice for large N.
2. **Engines are `Suggests:`**, checked at call time.
3. **Assembly: least-squares projection**, average available.
4. **Weights: equal**, inverse-variance available.
5. **Uncertainty: mixture of block draws** for every engine, REML through a
   normal approximation from `vcov`.
6. **Block size, coverage and LKJ** as in §9.
7. **Non-converged REML blocks are reported**, not refitted or dropped.
8. **`G` in the units of `Y`** (indri 0.2.0); estimation stays on the
   standardised scale.
9. **Licence: GPL-3**, to match stanimal. **Remote**: `diogro/indri`, private.

Open:

- **Missing data across blocks.** Each engine handles gaps within a block;
  how gaps that differ between traits affect the assembly, and stanimal's
  imputation against REML on gapped data, have not been measured.
- **Factor (horseshoe) prior for stanimal blocks.** It changes what G is
  estimated; not exposed yet.
- **A REML full model as reference** was dropped: full multi-trait REML fits are
  too slow and rarely converge without curated starting values.

## 13. Risks and caveats

- **Shrinkage and integration.** Every engine understates the strongest
  correlations, so integration (eigenvalue variance of `corrG`) is understated
  too: by about 17% for REML and 26% for stanimal at K = 35. The loss is in the
  block estimates (prior for stanimal; noise and the |r| <= 1 boundary for
  REML), not in the projection.
- **The uncertainty is pointwise and assembled.** Say so wherever intervals are
  reported.
- **Engine heterogeneity.** Three estimands, three uncertainty sources, four A
  forms; the contract and the cross-engine comparison exist to catch a
  plausible-looking G that is on the wrong scale, in the wrong order or from
  the wrong term.
- **Licensed dependency.** asreml cannot run in CI; its tests skip without a
  licence, and the engine fails with an actionable message.
- **Missing data**: a trait with many gaps is weakly estimated in every block
  it enters, and pair coverage does not see it.

## 14. References

- LiMMBo: Meyer, Casale, Stegle and Birney (2018),
  <https://www.biorxiv.org/content/10.1101/255497v1>; source in `../limmbo`.
- stanimal: `../stanimal`, interface guide `../stanimal/AGENTS.md`; call forms
  for every engine in `../stanimal/examples/package_comparison.R`.
- LKJ prior: Lewandowski, Kurowicka and Joe (2009).
- Nearest correlation matrix: Higham (2002).
- Regularized horseshoe: Piironen and Vehtari (2017).
