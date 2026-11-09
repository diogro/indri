# indri (development version)

## New features

* `indri_gwas()` runs genome-wide association scans with the covariances
  indri estimates: one multi-trait "any effect" test per SNP across all the
  traits, or one test per SNP and trait, in a linear mixed model with the
  relationship matrix `A` or in a linear model. The covariance is
  diagonalised once, so a scan costs one `N x N` rotation of the genotypes
  plus `O(N P)` per SNP. Every test is a likelihood ratio with the overall
  scale profiled.
  - The multi-trait mixed model uses `G` and `E` as the shape of the
    covariance and re-fits one ratio `delta` between them by REML under the
    null, so the test stays calibrated if they are off in scale. `delta` is
    returned (near 1 when `G` and `E` describe the data) and can be fixed
    with `delta =`.
  - Single-trait mixed models fit each trait's variance ratio by REML, or
    take it from the fit with `single_trait_variance = "fit"`.
  - The multi-trait linear model uses the residual covariance of `Y` after
    `X`.
  - `dominance = TRUE` adds a heterozygote indicator for each SNP:
    `pvalues` becomes the joint additive + dominance test,
    `pvalues_dominance` tests dominance given the additive effect, and `a`
    and `d` hold the two effects. SNPs seen in only two genotype classes
    keep the additive test (`dominance_estimable`).
  - `adjust = "bonferroni"` or `"effective"` corrects single-trait p-values
    for the number of traits.
  - Samples are matched by row name across `Y`, `genotypes`, `A` and `X`;
    missing genotypes are mean-imputed; monomorphic SNPs get `NA`.
  - `n_permutations` reruns the scan with permuted genotypes:
    `pvalues_empirical`, and `qvalues_empirical`, a permutation estimate of
    the false discovery rate; `fdr` gives the p-value threshold at that rate
    (`fdr_threshold`). LiMMBo's "empirical FDR" was the `fdr`-quantile of the
    permuted p-values, an uncorrected threshold; indri estimates the FDR.
  - `print()`, `summary()` (top SNPs, genomic inflation), `as.data.frame()`,
    and ggplot2 Manhattan and QQ plots (`plot()`); `snp_info` carries SNP
    positions into the results.
  - Checked against dense GLS and against LiMMBo's results on its own test
    data.
* `indri_effective_tests()` gives Galwey's (2009) effective number of
  independent traits, from the residuals after `X` when a design is given.
* `fit_G()` returns `E` and `E_bar`: the residual covariance, assembled from
  the block fits like `G` and kept positive definite, in the units of `Y`.
  Fits made by earlier versions have no `E` and cannot drive the
  multi-trait scan.
* `examples/gwas.R` runs the whole pipeline on simulated sibs in seconds.
* Data from an F3 intercross of the LG/J and SM/J mouse strains:
  `lgsm_f3` (1611 mice with pedigree, rearing, and four trait classes:
  weekly weight gains, raw weekly weights, skull areas and necropsy
  weights; `lgsm_traits` lists them), `lgsm_genotypes` (353 markers in the additive, dominance and
  imprinting design, for the F3 mice and their F2 parents) and
  `lgsm_markers` (the map). Provenance still to be documented.

# indri 0.2.0

## Breaking changes

* `G`, `G_bar`, `G_draws` and `sigma_G` are returned in the units of `Y`.
  Fits saved by 0.1.0 hold a standardised `G`.

## New features

* The REML engines (sommer, asreml) report uncertainty: normal draws from
  each block's sampling covariance, assembled like stanimal's posterior
  draws.
* `spread = "mixture"` is the default way to assemble draws, and
  `corrG_draws` mixes each block's own correlation draws. Averaging the
  block draws had covered the true correlations 40% of the time at nominal
  95%; the mixture brought that to nominal.
* Non-converged REML blocks are reported in `fit$diagnostics` instead of
  trusting a single convergence flag.
* REML engines take starting values and fit single-trait blocks.
* stanimal blocks share one eigendecomposition of `A` (needs stanimal
  >= 0.8.0).

## Defaults

* Subsets of up to four traits (`s = 3`) for every engine; each pair in 3
  distinct blocks for stanimal and in at least 1 for the REML engines.
* The stanimal engine fits blocks with `lkj_prior = 1`.
* asreml runs two workers by default (one licence seat each).

## Bug fixes

* sommer returned blocks in alphabetical trait order, which put variances
  and covariances on the wrong traits; blocks now follow the columns of `Y`.

# indri 0.1.0

* First version: trait subsets with guaranteed pair coverage
  (`indri_subsets()`), block assembly (`indri_assemble()`,
  `indri_assemble_draws()`), a least-squares projection onto the PSD cone
  (`indri_project_psd()`), pluggable block engines (stanimal, sommer,
  asreml) and `fit_G()`.
