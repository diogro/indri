# indri

Estimates the **additive genetic covariance matrix G** (traits × traits) for
large trait sets, without fitting one giant multivariate animal model.

It does this LiMMBo's way: split the traits into many small, overlapping
subsets, fit a small animal model on each, and reassemble the estimates into a
single positive semi-definite matrix. The per-subset fit is pluggable, so the
same estimator runs on a Bayesian sampler or on a REML engine, and every engine
comes with per-entry uncertainty. See the [LiMMBo paper](#reference).

The assembled `G`, with the residual covariance `E` (the residual
covariance of `Y` after `X` minus `G`), then drives **genome-wide association scans**: one multi-trait test per SNP
across all the traits, or trait-by-trait tests, in a linear mixed model, with
optional dominance effects (`indri_gwas()`, [below](#association-scans)).

> **Status (0.3.0):** working and validated on simulations of up to 35 traits;
> association scans calibrated at 10 and 30 traits.
> The assembled `G` is as accurate as the full model where the full model can
> be fitted; the 95% bands cover the true correlations about 95-97% of the
> time. See [STATUS.md](STATUS.md) for the numbers.

## Install

```r
remotes::install_github("diogro/indri")
```

You also need at least one block engine:

| Engine | Package | Needs | Uncertainty |
|---|---|---|---|
| `"stanimal"` (default) | [`stanimal`](https://github.com/diogro/stanimal) (>= 0.8.0 recommended) | CmdStan | posterior draws |
| `"sommer"` | [`sommer`](https://cran.r-project.org/package=sommer) | | REML sampling distribution |
| `"asreml"` | [ASReml-R](https://vsni.co.uk/software/asreml-r/) | a licence (one seat per worker) | REML sampling distribution |

Every engine accepts gaps (`NA`) in `Y` and uses each individual's observed
traits. stanimal gives a posterior; the REML engines are much faster for large
N.

## Quick start

```r
library(indri)
set.seed(1)

# A known 6-trait genetic covariance on a sibship design.
G_true <- stanimal::random_covariance(
  n_traits = 6, method = "partial", min_var = 0.8, max_var = 1.4
)
E_true <- stanimal::random_covariance(
  n_traits = 6, method = "lkj", lkj_prior = 2, min_var = 0.5, max_var = 0.9
)
A <- kronecker(diag(20), matrix(0.5, 8, 8))
diag(A) <- 1
Y <- t(chol(A)) %*% matrix(stats::rnorm(nrow(A) * 6), nrow(A), 6) %*% chol(G_true) +
     matrix(stats::rnorm(nrow(A) * 6), nrow(A), 6) %*% chol(E_true)
colnames(Y) <- paste0("trait", 1:6)

fit <- fit_G(
  Y, A = A, engine = "stanimal",
  engine_args = list(chains = 1, iter = 600, refresh = 0, show_messages = FALSE)
)

fit                                            # summary print
fit$G                                          # G, in the units of Y
fit$corrG                                      # genetic correlations
apply(fit$corrG_draws, 1:2, quantile, 0.975)   # per-entry upper 95% bound
```

Fixed effects go in a design matrix with one row per individual, which must
include an intercept (or every level of a factor). stanimal's
`animal_design()` builds all three inputs from a formula, a data frame and a
pedigree:

```r
d <- stanimal::animal_design(cbind(t1, t2, t3) ~ sex + age, data = dat,
                             pedigree = ped)
fit <- fit_G(d$Y, X = d$X, A = d$A)
```

## What you get back

| slot | what it is |
|---|---|
| `G`, `corrG`, `sigma_G` | the assembled estimate (`G` and `sigma_G` in the units of `Y`) |
| `G_draws`, `corrG_draws` | assembled draws (traits × traits × draws), for per-entry intervals |
| `G_bar` | the average of the blocks before the PSD projection |
| `E`, `E_bar` | the residual covariance for association scans (residual covariance of `Y` after `X` minus `G`, kept well conditioned; units of `Y`), and the average of the block residual covariances |
| `subsets`, `coverage` | the trait subsets fitted, and how many covered each pair |
| `diagnostics` | per-block convergence; uncertified blocks listed by trait |
| `blocks` | the raw block fits, on the internal standardised scale |
| `y_mean`, `y_sd` | that standardisation |

`print()`, `summary()` and `plot()` are provided. A progress bar of the block
fits shows in interactive sessions (`progress = TRUE` forces it).

## Association scans

`indri_gwas()` tests SNPs against the traits. The multi-trait test is the
reason to estimate `G` and `E` for many traits at once: it tests each SNP
against all of them jointly, in a mixed model whose covariance is `G` and `E`
with the relationship matrix `A`.

```r
fit <- fit_G(Y, X = X, A = A)                  # G and E, once
scan <- indri_gwas(Y, genotypes, fit = fit, A = A, X = X,
                   snp_info = snp_info)        # one test per SNP, all traits
scan                                           # summary print
summary(scan)                                  # top SNPs, genomic inflation
plot(scan)                                     # Manhattan plot (ggplot2)
plot(scan, type = "qq")                        # QQ plot
head(as.data.frame(scan))                      # one row per SNP, for write.csv()

# trait by trait, adjusted for the effective number of traits
single <- indri_gwas(Y, genotypes, A = A, X = X, mode = "single_trait",
                     adjust = "effective")

# additive and dominance effects
dom <- indri_gwas(Y, genotypes, fit = fit, A = A, X = X, dominance = TRUE)
```

`genotypes` is an individuals × SNPs matrix (typically 0/1/2); samples are
matched to `Y`, `A` and `X` by row name. `snp_info` is an optional data frame
(for example `snp`, `chrom`, `pos`) carried into the results.

| `mode`, `model` | test per SNP | covariance |
|---|---|---|
| `"multi_trait"`, `"lmm"` (default) | joint, `P` degrees of freedom | `G ⊗ A + δ E ⊗ I` |
| `"multi_trait"`, `"lm"` | joint, `P` degrees of freedom | `E ⊗ I`, `E` the residual covariance after `X` |
| `"single_trait"`, `"lmm"` | one per trait | `A + δ I` per trait, `δ` by REML or from the fit |
| `"single_trait"`, `"lm"` | one per trait | ordinary least squares |

Every test is a likelihood ratio with the overall scale profiled. The
covariance is diagonalised once, so a scan costs one `N × N` rotation of the
genotypes plus `O(N P)` per SNP: thousands of SNPs take seconds.

- **`delta`.** `G` and `E` set the shape of the covariance; the ratio `δ`
  between them is re-fitted once under the null (no SNP) and held fixed for
  every SNP. With a `G` and `E` from the same data it comes out near 1, and
  it is returned (`scan$delta`) so that you can check. A value far from 1
  says the covariances do not describe the data being scanned; the test stays
  calibrated either way. Pass `delta =` to fix it instead.
- **Dominance.** `dominance = TRUE` adds a heterozygote indicator next to
  each SNP's dosage: `pvalues` becomes the joint additive + dominance test,
  `pvalues_dominance` tests dominance given the additive effect, and `a` and
  `d` are the additive and dominance effects (Falconer's symbols). It needs
  genotypes coded 0/1/2 and all three genotypes in numbers: a SNP seen in
  only two genotype classes keeps the additive test and gets `NA` dominance
  results (`dominance_estimable`).
- **Genotype by environment.** `gxe = "sex"` (any columns of `X`) adds each
  SNP's interaction with those environments: `pvalues` becomes the joint
  test of the SNP and its interactions, `pvalues_gxe` tests the
  interactions given the main effect, and `gxe` holds their effects. GxE in
  the background (a polygenic G x E, or residual variance differing between
  environments) is not modelled by default: pass a `fit_G_env()` fit
  (experimental) as `fit` for a background with its own `G` in each
  environment and between environments.
- **Permutations.** `n_permutations` reruns the scan with the genotypes
  permuted across individuals: `pvalues_empirical` per SNP, and
  `qvalues_empirical`, the permutation estimate of the false discovery rate.
  With `fdr = 0.05`, `fdr_threshold` is the p-value at or below which SNPs are
  discoveries at that rate. Permuting across relatives ignores their
  relatedness, so treat these as a check on the analytic p-values.
- **Adjusting single-trait p-values** for the number of traits:
  `adjust = "bonferroni"` multiplies by the number of traits, `"effective"`
  by the effective number of independent traits (Galwey 2009,
  `indri_effective_tests()`), counted on the residuals after `X`.
- **Inputs.** `Y` must be complete. Missing genotypes are replaced by the
  SNP's mean; monomorphic SNPs get `NA` with a warning. The mixed models need
  `A`; the multi-trait one needs `G` and `E`, from `fit` or given directly.

| slot | what it is |
|---|---|
| `pvalues` | one per SNP (multi-trait) or SNP × trait (single-trait) |
| `a` | additive effects, SNP × trait, in the units of `Y` per allele |
| `delta` | the null variance ratio (one, or one per trait) |
| `pvalues_dominance`, `d`, `dominance_estimable` | with `dominance = TRUE` |
| `pvalues_gxe`, `gxe`, `gxe_estimable` | with `gxe` |
| `pvalues_adjust`, `n_tests` | with `adjust` |
| `pvalues_empirical`, `qvalues_empirical`, `fdr_threshold` | with `n_permutations` (and `fdr`) |
| `snp_info` | the SNP information, or the SNP names |

`indri_locus_G(scan, genotypes = genotypes)` turns the loci's `a` and `d`
into the additive and dominance covariance matrices they imply, with the
allele frequencies and linkage disequilibrium of the genotypes (Melo,
Marroig and Wolf 2019): how the mapped loci build up the genetic
covariance between traits.

**Data.** The package ships an F3 intercross of the LG/J and SM/J mouse
strains: `lgsm_f3` (1611 mice; growth (weekly gains), raw weekly weight,
mandible area and necropsy traits, listed in `lgsm_traits`), `lgsm_genotypes` (353 markers coded additive `A`,
dominance `D` and imprinting `I`; `A + 1` is the 0/1/2 dosage
`indri_gwas()` takes), `lgsm_markers` (the map), and `lgsm_pedigree` and
`lgsm_A` (the F3 relationship matrix; the F2 parents are unrelated founders,
as the P strains are inbred and the F1 all heterozygous).

[`examples/gwas.R`](examples/gwas.R) runs the whole pipeline on simulated
sibs in under a minute, with four planted SNPs showing where each test is
strongest, permutations, and the plots.

## Caveats worth knowing

- **Intervals are per entry.** The draws are assembled from dependent blocks:
  use them for intervals on individual covariances and correlations, not for
  functions of the whole matrix (eigenvalues, integration). Measured coverage
  at nominal 95%: 0.96-0.98 for the REML engines, about 0.95 for stanimal.
- **Strong correlations are shrunk a little**, so integration is understated
  rather than overstated: the eigenvalue variance of `corrG` came out at about
  75-85% of the truth's at 35 traits.
- **Blocks are small by default**: up to four traits each, every pair in three
  distinct blocks for stanimal (fitted with `lkj_prior = 1`) and in at least one
  for the REML engines. Larger blocks shrank stanimal's correlations and gave
  the REML engines nothing for their cost.
- **asreml needs a licence seat per worker.** It runs two workers by default; a
  standard licence allows two processes on one machine at a time. Use
  `cores = 1` with a single seat.
- **Units.** `G` is in the units of `Y` since 0.2.0; fits saved by 0.1.0 hold a
  standardised `G`.
- **Association scans** condition on the point estimates of `G` and `E`;
  their uncertainty is not propagated. Fits made before `E` was assembled
  need refitting to drive the multi-trait test. Genotypes are read as
  matrices only (no PLINK reader).

## Repository map

- [`NOTES.md`](NOTES.md): design, engine contract and the reasons for the
  defaults.
- [`STATUS.md`](STATUS.md): what has been measured and what comes next.
- [`AGENTS.md`](AGENTS.md): guidance for coding agents working on the package.
- [`examples/`](examples/README.md): worked examples, including the genome
  scan in `gwas.R`.
- [`GWAS_PLAN.md`](GWAS_PLAN.md): the association-testing port from LiMMBo,
  its decisions and what remains.

## Scope

In scope: `G` and `E` from a phenotype matrix `Y`, a design `X` and a known
relationship matrix `A`, and association scans that use them. Out of scope
for now: `P` and `h2` reporting, breeding values, imputation, model
comparison.

## Reference

Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
a simple, scalable approach for linear mixed models in high-dimensional genetic
association studies. *bioRxiv*, preprint.
<https://www.biorxiv.org/content/10.1101/255497v1>
