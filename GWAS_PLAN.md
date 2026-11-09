# Plan: association mapping in indri

Status: in progress on the `gwas` branch. Written 2026-11-08.

Goal: port LiMMBo's association-testing layer (`limmbo/core/gwas.py` and the
`runAssociation` pipeline) to indri, so that the `G` indri estimates can drive a
multi-trait mixed-model genome scan. This reverses one line of NOTES.md §2
("association testing" is out of scope); update NOTES.md when work starts.

## 1. What LiMMBo does

LiMMBo's `GWAS` class wraps LIMIX 1.x. Everything it offers, and whether it
comes over:

| LiMMBo feature | LIMIX call | indri |
|---|---|---|
| multi-trait "any effect" LMM: one test per SNP, P degrees of freedom | `qtl_test_lmm_kronecker(Asnps = I, Acovs = I, K1c = Cg, K2c = Cn, K1r = K)` | **port** (the core) |
| single-trait LMM, one test per SNP and trait | `qtl_test_lmm(K = K)` | port |
| multi-trait LM: `Cg = 1e-9 I`, `Cn = cov(Y)`, `K = I` | same Kronecker call | port (falls out of the LMM code) |
| single-trait LM | `qtl_test_lmm(K = NULL)` | port |
| principal components as covariates | appended to the covariates | port, as columns of `X` |
| per-SNP effect sizes, P per SNP | `getBetaSNP()` | port |
| single-trait p-value adjustment: Bonferroni, or effective number of tests (Galwey 2009) | `effectiveTests()`: `(sum sqrt(lambda))^2 / sum(lambda)` from the Spearman correlation of `Y` | port |
| empirical p-values by permuting genotype rows | loop over the scan | port, with a caveat (§5) |
| empirical FDR: 10 permutations, threshold at rank `fdr * S * 10` | loop over the scan | port |
| Manhattan and QQ plots | `limix.plot`, `limmbo.plot.manhattan` | port, base graphics |
| CSV output: p-values, betas, adjusted and empirical p | `saveAssociationResults()` | `as.data.frame()`; users call `write.csv()` |
| REML `Cg`/`Cn` when none are given | `limix.mtset.MTSet` | **replace** with `fit_G()` |
| CLI, file readers, sample-ID matching | `parser.py`, `ReadData`, `InputData` | not ported: NOTES.md keeps CLI parity out of scope; R users pass matrices |

## 2. What indri is missing

- **`E` is not assembled.** The multi-trait test needs `Cg` *and* `Cn`. Every
  engine already returns a block `E` (`fit$blocks[[j]]$E`), so assembling it is
  the same average-and-project path as `G`, without draws. `fit_G()` gains `E`
  and `E_bar` in the units of `Y`. `E` must be positive definite for the
  whitening, so the projection floors its eigenvalues. LiMMBo added `1e-4 I`
  to both matrices.
- **No genotype handling.** The scan takes an `N x M` numeric matrix with
  sample row names, plus an optional data frame of SNP information
  (`chr`, `pos`, `id`). Missing genotypes are mean-imputed per SNP, as in
  `stanimal::genomic_relationship()`.

## 3. The computation

Model, for `N` individuals, `P` traits and one SNP `g`:

```
Y = X B + g beta' + U + Eps,   vec(U) ~ N(0, Cg (x) A),   vec(Eps) ~ N(0, Cn (x) I)
```

With `Cg` and `Cn` fixed, as in LiMMBo, the covariance diagonalises once for
the whole scan:

1. `A = Q diag(lambda) Q'`, decomposed once. indri already has this through
   the `"eigen"` A form in `a_forms.R`.
2. `Cn = L L'`, and `L^-1 Cg L^-T = W diag(d) W'`. Let `T = L^-T W`.
3. Rotate: `Y* = Q' Y T`, `X* = Q' X`, `G* = Q' G`, where `G` holds the
   genotypes. Entry `(i, k)` of `Y*` is then independent, with variance
   `lambda_i d_k + 1`.
4. Under the rotation the fixed effects stay free (`B T` and `T' beta` are
   reparameterisations). So the any-effect test splits into **P independent
   weighted regressions**: column `k` of `Y*` on `[X*, g*]` with weights
   `1 / (lambda_i d_k + delta)`.
5. **`G` and `E` are shapes, not literal covariances.** The covariance is
   `sigma^2 (G (x) A + delta E (x) I)`. The ratio `delta` is fitted once under
   the null by REML, unless the caller passes it. `sigma^2` is profiled under
   each model. The statistic is the likelihood ratio
   `NP log(q0 / q1) ~ chi^2_P`, where `q` is the generalised residual sum of
   squares. The per-trait effects are `beta = T'^-1 beta_tilde`.
   - Re-fitting `delta` keeps the test calibrated when `G` and `E` are off
     scale. limmbo's shipped `Cg`/`Cn` gave `delta = 1.35` on a 250-sample
     subset, and taking them literally moved a p-value from 0.12 to 0.38.
   - With indri's own `G` and `E`, fitted on the same data, `delta` should
     be near 1. It is reported, and a value far from 1 flags a mismatch.

Cost:

- The rotation `Q' G` is `O(N^2 M)` and dominates. Rotate in chunks of SNPs.
- Each column then costs `O(N J)` per SNP. Residualise `X*` once per column
  (Frisch-Waugh), so a chunk is a few matrix products.
- Nothing else is `O(N^2)` per SNP.

This uses the same rotate-once idea as stanimal's `animalModelG.stan`.

- **Single-trait LMM:** the same rotation with `P = 1`, covariance
  `sigma^2 (A + delta I)`. Each trait's `delta` is fitted by REML under the
  null and held fixed across SNPs (the P3D/EMMAX approximation). The option
  `single_trait_variance = "fit"` takes `delta = E_kk / G_kk` from the indri
  fit instead.
- **LM variants:** the multi-trait LM is `sigma^2 E (x) I`, with `E` by
  default the covariance of the residuals after `X`. limmbo used `cov(Y)`
  of the raw phenotypes, whose shape includes the covariates' effect. The
  single-trait LM is OLS.

## 4. Proposed interface

Names follow the package (`fit_G`, `indri_*`):

```r
scan <- indri_gwas(Y, genotypes, fit = NULL, G = NULL, E = NULL, A = NULL,
                   X = NULL, snp_info = NULL,
                   mode = c("multi_trait", "single_trait"),
                   model = c("lmm", "lm"),
                   single_trait_variance = c("reml", "fit"), delta = NULL,
                   dominance = FALSE,
                   adjust = c("none", "bonferroni", "effective"),
                   chunk_size = 1000L)
indri_empirical_p(scan, n_permutations = 1000L, seed = NULL, cores = 1L)
indri_empirical_fdr(scan, fdr = 0.05, n_permutations = 10L, seed = NULL)
indri_effective_tests(Y)
plot(scan, type = c("manhattan", "qq"))
as.data.frame(scan)        # one row per SNP (and trait), for write.csv()
```

- **Where `G` and `E` come from:** `fit` (an `indri_fit`) or explicit
  `G`/`E`, matched to the columns of `Y` by name. One of them is required for
  the multi-trait LMM. The single-trait LMM needs neither unless
  `single_trait_variance = "fit"`.
- **`A` is required for the LMMs**, as in `fit_G()`. An indri fit does not
  keep it, and storing it would make every fit `N x N`. A scan on a genomic relationship
  matrix with `G` estimated on a pedigree `A` is the caller's choice.
- **`delta`** fixes the null variance ratio instead of estimating it: a
  scalar for the multi-trait test, one value per trait for the single-trait
  test.
- **`X` is checked** with `.indri_check_design()`. Principal components are
  extra columns of `X`.
- **Requirements:** `Y` must be complete in v1, since a missing entry breaks
  the Kronecker rotation. Sample IDs are matched by row names across `Y`,
  `A` and `genotypes`, the way `fit_G()` matches `A`.
- **Return value:** an `indri_gwas` object holding `pvalues` (M, or M x P),
  `a` (additive effects, M x P), `delta`, `pvalues_adjust` and `n_tests` for single-trait
  scans with an adjustment, `snp_info`, `call`, `mode` and `model`.
- **Genotypes:** missing values are mean-imputed per SNP. A monomorphic SNP
  gets `NA` p-values and effects, with a warning.

## 5. Caveats to settle before porting

- **Permutation under relatedness.** LiMMBo permutes genotype rows across all
  samples. With related individuals the samples are not exchangeable, so the
  permuted null is not the null of the model. Port it as LiMMBo had it,
  document the caveat, and leave a better null for later, e.g. permuting in
  the rotated space or a parametric bootstrap. Each permutation re-rotates
  the genotypes (`O(N^2 M)`). The permutations are independent and pure R, so
  they can fork (AGENTS.md rule 1 only forbids forking around cmdstanr).
- **Empirical FDR** in LiMMBo uses exactly 10 permutations and the rank
  `fdr * S * 10`. Keep the procedure and expose `n_permutations`.
- **Proximal contamination / LOCO** is not in LiMMBo. Leave it out of v1 and
  note it for later.
- **Uncertainty in `G` and `E`** is not propagated, as in LiMMBo. The scan
  conditions on the point estimates.

## 6. Milestones

Each milestone ends with tests passing (`devtools::test("package")`).

| | work | done when |
|---|---|---|
| **A0** | `E` and `E_bar` assembled in `fit_G()`, in the units of `Y`, eigenvalues floored; NOTES.md §2 updated | **done**: `test-fit-E.R`, `E` matches the full-model REML `E` to 1e-3 |
| **A1** | rotation core: multi- and single-trait LMM and LM | **done**: `test-gwas-dense.R` matches dense GLS (`helper-gwas.R`) to 1e-8, and the REML `delta` to 1e-5 |
| **A2** | limmbo reference fixtures (§7) | **done**: `test-gwas-limmbo.R` matches limmbo given its `delta` |
| **A3** | adjustment, empirical p, empirical FDR | **partly**: Bonferroni and effective adjustment done and tested; empirical p and empirical FDR not started |
| **A4** | `indri_gwas()` interface, `fit =` path, sample matching, `as.data.frame()`, `print`/`summary`, Manhattan and QQ plots | **mostly**: interface, `fit =`, matching, `print`, `as.data.frame` and `examples/gwas.R` done; no `summary()` or `plot()` method yet (the example draws its own Manhattan and QQ) |
| **A5** | calibration and power studies (§8), STATUS.md entry | **not started**; `examples/gwas.R` gives a first look (λ_GC 0.98-1.05 on 3 seeds) |

## 7. Validation against LiMMBo

`limmbo/io/test/data/` has 1000 samples, 10 traits, 20 SNPs, covariates,
PCs, a relationship matrix, and `Cg`/`Cn`. The simulation is described in
`datadescription.md`.

Done. `fixtures/limmbo/make_limmbo_fixtures.py` runs limmbo's `GWAS` class
on the first 250 samples, with the shipped `Cg`/`Cn`, and writes every input
and output. Run it from the limmbo checkout with
`PYTHONPATH=. pixi run python <script> <outdir>`. The relationship matrix
(limmbo's `K.csv`) is stored as `A.rds`, an
xz-compressed `.rds`, so the fixtures take 520 KB.

What LIMIX does, found from these fixtures and its source:

- **Both tests re-fit an overall variance ratio under the null by ML,**
  `V = Cg (x) K + delta Cn (x) I` (or `K + delta I`). They hold it fixed
  under the alternative and profile `sigma^2`.
- **The search runs over log delta in [-5, 5].** It refines the grid optimum
  for the multi-trait test but not for the single-trait test (100 grid points
  only). Five of the 10 traits sit at the upper bound.
- **limmbo's reader appends an intercept column to the covariates.**

The tests pass LIMIX's fitted `delta` through `indri_gwas(delta =)`, so they
check the scan and not LIMIX's optimiser. LIMIX is a test reference only and
is not mentioned in user documentation.

## 8. Studies

They run on a machine with spare cores, synced through git, as in AGENTS.md
"Running studies":

- `examples/gwas_calibration.R` simulates G and E from
  `stanimal::random_covariance()` on the sibship design, with null SNPs. It
  reports λ_GC and type-I error at 1e-2 to 1e-4 for each mode, with `Cg`/`Cn`
  from indri against the truth.
- `examples/gwas_power.R` adds SNPs affecting 1, a few, or all traits. It
  reports power of the multi-trait test against single-trait with Bonferroni
  or effective adjustment. This is LiMMBo's paper comparison, now using
  indri's `G`.

## 9. Files

| file | content |
|---|---|
| `package/R/gwas.R` | `indri_gwas()`, rotation core, chunking |
| `package/R/gwas_null.R` | `indri_effective_tests()`, `indri_empirical_p()`, `indri_empirical_fdr()` |
| `package/R/methods.R` | `print`, `summary`, `plot`, `as.data.frame` for `indri_gwas` |
| `package/R/fit_G.R`, `assemble.R` | `E` assembly |
| `package/tests/testthat/test-gwas-*.R` | brute-force GLS, LiMMBo fixtures, adjustment, permutations, interface |
| `examples/gwas.R`, `gwas_calibration.R`, `gwas_power.R` | worked example and the two studies |

No new hard dependencies. A PLINK reader (`BEDMatrix` or `genio`) can be a
`Suggests:` later, since v1 takes a matrix.

## 10. Decisions (2026-11-08)

1. **Name:** `indri_gwas()`.
2. **Single-trait variance ratios:** both, through an argument
   (`single_trait_variance = c("reml", "fit")`); per-trait REML is the default,
   for parity with LIMIX and limmbo.
3. **Permutations:** port limmbo's genotype permutation as is and document
   that it ignores relatedness; a relatedness-aware null comes later.
4. **Input:** matrices only in v1; a PLINK reader can be a `Suggests:` later.
5. **LIMIX stays out of user documentation.** It is a test reference only;
   the multi-trait LMM re-fits `delta` and `sigma^2` and reports `delta`.
6. **`delta` is fitted by REML** for both tests (LIMIX used ML).
7. **Multi-trait LM residual shape:** the covariance of the residuals after
   `X`, not limmbo's `cov(Y)`.
8. **Effective number of tests:** from the residuals after `X`
   (`indri_effective_tests(Y, X)`), not limmbo's raw phenotypes.
9. **Dominance (added 2026-11-08):** `dominance = TRUE` adds a heterozygote
   indicator per SNP next to its 0/1/2 dosage. `pvalues` becomes the joint
   additive + dominance test (2 df per trait), `pvalues_dominance` tests
   dominance given additive (1 df per trait), and `a`/`d` are the joint
   model's additive and dominance effects (Falconer's symbols). SNPs seen in only two genotype classes keep
   the additive-only test and are `FALSE` in `dominance_estimable`.
   `test-gwas-dominance.R` checks it against dense nested GLS.

Work happens on the `gwas` branch, tests first.
