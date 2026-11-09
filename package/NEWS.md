# indri (development version)

* `fit_G_env()` (experimental) estimates the genetic covariance of traits
  within and between the groups of a discrete environment (sex, diet,
  site), fitting each trait in each environment as its own "character
  state" with `fit_G()` (sommer by default). `G` holds the
  within-environment blocks and the genetic covariance between
  environments (`corr_between`: each trait's cross-environment genetic
  correlation); `E` is block diagonal, as no individual is observed in two
  environments. The first step towards modelling polygenic GxE in the
  association scans: pass the result to `indri_gwas()` as `fit` and the
  scan uses an environment-specific background covariance (each pair of
  individuals' genetic covariance depends on both their environments). The
  scan rotates by a generalized eigendecomposition per connected component
  of `A`, so it stays fast for family data; it matches a dense GLS, and
  reduces to the usual scan when the background has no GxE.

* `indri_gwas(gxe = )` fits genotype-by-environment interactions: each SNP's
  dosage also enters multiplied by the named columns of `X` (sex, diet, a
  covariate), whose main effects are in the model already. `pvalues` is the
  joint test of the SNP's main effect and interactions, `pvalues_gxe` the
  interactions given the main effect, with the interaction effects in
  `gxe` (SNPs by traits by environments). Works for multi- and single-trait
  scans, LMM or LM, with dominance (as a main effect) and permutations;
  SNPs whose interactions cannot be estimated keep the main-effect test
  (`gxe_estimable`). `plot(scan, which = "gxe")` draws the interaction
  tests. GxE in the background covariance is not modelled.

* `indri_locus_G()` builds the additive (`G`) and dominance (`D`) covariance
  matrices implied by per-locus additive
  and dominance effects, allele frequencies and linkage disequilibrium
  (Melo, Marroig and Wolf 2019, Eqs. 7-14). It takes `indri_gwas()`'s `a`
  and `d` directly, with frequencies and composite LD from the genotypes,
  or effects from other software: a data frame with a column of locus
  names, traits in rows or in columns (`traits_in`), loci matched to the genotypes
  (and to `d`, `p`, `lambda`) by name, and errors for 0/1/2 codings and
  transposed effects. The help page sets out the effect-allele, dominance
  and unit conventions the effects must follow.

# indri 0.3.0

Association scans, the residual covariance, and data from an LG/SM
intercross. `fit$E` is new; see below for how it is computed.

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
* `fit_G()` returns `E`, the residual covariance in the units of `Y`: the
  covariance of the residuals of `Y` on `X` minus the assembled `G` (scaled
  for the relatedness in `A`), with eigenvalues floored at 1e-3 of the
  largest. `E_bar` is the average of the block residual covariances. The
  block average is not used for `E`: at 30 traits it was strongly
  indefinite, its projection had eigenvalues near zero, and the
  multi-trait test ran at 66 times its nominal rate at 1e-4.
  Fits made by earlier versions have no `E` and cannot drive the
  multi-trait scan.
* `examples/gwas.R` runs the whole pipeline on simulated sibs in seconds.
* Data from an F3 intercross of the LG/J and SM/J mouse strains:
  `lgsm_f3` (1611 mice with pedigree, rearing, and four trait classes:
  weekly weight gains, raw weekly weights, mandible areas and necropsy
  weights; `lgsm_traits` lists them), `lgsm_genotypes` (353 markers in the additive, dominance and
  imprinting design, for the F3 mice and their F2 parents), `lgsm_markers`
  (the map), `lgsm_pedigree` and `lgsm_A` (the F3 relationship matrix, with
  the F2 parents as unrelated founders: the P strains are inbred and the F1
  all heterozygous). The four sires missing in the source were inferred
  from the markers by Mendelian exclusion (`sire_inferred`). Provenance still to be documented.

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
