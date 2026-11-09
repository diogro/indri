# Plan: four vignettes on the LG/SM growth data

Status: proposal. Written 2026-11-09.

Four vignettes, each a worked analysis of the growth traits of the LG/SM F3
intercross shipped with indri (`lgsm_f3`, `lgsm_genotypes`, `lgsm_A`):

1. **G-matrix estimation** with the three engines, and how to choose the LKJ
   prior for stanimal.
2. **Genome scans** for additive and dominance effects.
3. **G from locus effects** (`indri_locus_G()`, Melo, Marroig and Wolf 2019).
4. **Genotype-by-environment scans** with sex as the environment.

They share one data preparation and build on each other: 1 fits the G and
E that 2 uses; 2 finds the loci whose effects 3 turns into a G; 4 extends 1
and 2 to environment-specific covariances.

## 1. Decisions before writing

1. **Which growth traits.** `lgsm_traits$growth` holds the 7 weekly gains
   of the original analyses (weeks 1-8); `lgsm_traits$weight` holds the 10
   raw weekly weights. The data reach week 10, so the gains could run to 9
   (`growth12` to `growth910`). **Recommendation: 9 weekly gains.** Gains
   are less correlated with each other than cumulative weights (|r| up to
   0.28 among the 7), so the G matrix has more to show, and 9 traits give
   36 pairs: a real test of the subset assembly without long fits. This
   needs `growth89` and `growth910` added to `lgsm_f3` (`data-raw/lgsm.R`)
   and `lgsm_traits$growth` widened, or a separate `growth9` set to keep the
   paper's 7. Alternative: the 10 raw weights (simpler, but G of cumulative
   weights is dominated by one size axis).
2. **Maternal and litter effects.** Full sibs share a dam, and growth to
   weaning depends on her. indri's models have one random effect (`A`), so
   a common-litter effect inflates `G` (the paper used a dam random
   effect). 632 mice were cross-fostered (nurse other than dam), which
   separates prenatal from nursing environment in part. **Recommendation:**
   fixed effects for sex, cohort, litter size at birth and at weaning in
   `X` (the original analyses' fixed effects); say plainly in vignettes 1
   and 4 that `G` here includes shared-litter variance; and list a second
   random effect (litter or nurse) as a feature request in STATUS.md.
3. **Precomputed or live.** Fitting nine traits on 1552 mice with stanimal
   takes far longer than a vignette build should (sommer and asreml take
   seconds to minutes). **Recommendation:** knitr's "precomputed vignette"
   pattern: each vignette has a source `vignettes/<name>.Rmd.orig` run on
   meyer (or locally) to produce `vignettes/<name>.Rmd` with outputs and
   figures in place; `R CMD build` then only renders Markdown. A
   `vignettes/precompute.R` script rebuilds all four. Alternative: ship
   saved fits in `inst/extdata` and run the light code live.
4. **Matrix comparison tools.** Vignettes 1 and 3 compare matrices (random
   skewers, Krzanowski subspace similarity, Mantel). evolqg has them, but is
   not a dependency. **Recommendation:** short helpers written in the
   vignettes (random skewers is ten lines), or evolqg in `Suggests:` with a
   fallback. No new hard dependency either way.

## 2. Shared setup (all four vignettes)

```r
library(indri)
traits <- lgsm_traits$growth                 # after decision 1
d <- lgsm_f3[lgsm_f3$ID %in% rownames(lgsm_genotypes$A), ]
d <- d[stats::complete.cases(d[c(traits, "LSB", "LSW")]), ]
Y <- as.matrix(d[traits]); rownames(Y) <- d$ID
X <- model.matrix(~ SEX + factor(COHORT) + LSB + LSW, d); rownames(X) <- d$ID
A <- lgsm_A[rownames(Y), rownames(Y)]
dosage <- lgsm_genotypes$A[rownames(Y), ] + 1     # 0/1/2, as indri_gwas() takes
info <- transform(lgsm_markers, pos = bp)        # or cM (rs3684333 has no bp)
```

About 1548 mice (1552 genotyped, 4 without litter sizes), 9 traits, 353
markers, 200 full-sib families of about 8, linked by 9 sires with two
mates. A short section in vignette 1 introduces the cross (P strains
inbred, F1 all heterozygous, F2 parents unrelated founders, F3 sibships),
the marker coding (A/D/I) and `lgsm_A`; the others point back to it.

## 3. Vignette 1: estimating G (`indri-G.Rmd`)

Audience: users who want G for many traits. Questions answered: what
indri does, which engine to pick, which prior.

1. **The idea in one paragraph:** trait subsets, block fits, assembly,
   projection (link to the LiMMBo paper and NOTES.md).
2. **The data:** the cross, the nine gains, the design `X`, `A`.
3. **Three engines on the same data:**
   - `fit_G(Y, X, A, engine = "sommer")`: seconds; REML; uncertainty from
     each block's sampling covariance.
   - `engine = "asreml"`: same model, needs a licence; agreement with
     sommer (expected within ~0.01, STATUS.md).
   - `engine = "stanimal"`: posterior draws; the slow one; `cores`,
     `engine_args = list(chains, iter)`.
   - A table of time, and a figure comparing the three `corrG` (scatter of
     pairs against each other, and heat maps).
4. **Reading a fit:** `print`, `summary`, `plot`, `corrG`, intervals from
   `corrG_draws`, `diagnostics` (unconverged blocks, Rhat), coverage of
   pairs, `G_bar` against `G` (what the projection changed), `E`.
5. **Choosing the LKJ prior for stanimal** (the section asked for):
   - What `lkj_prior` does in a block (prior on block correlations;
     uniform at 1, shrinking towards 0 above 1), and why indri's default is
     1 (simulations: 4-trait blocks at LKJ(1) gave nominal coverage and the
     least shrinkage; LKJ(2) undercovered, STATUS.md "Uncertainty").
   - When to raise it: blocks that do not mix (Rhat >= 1.05, divergences;
     indri's warning says so), many zero eigenvalues in `corrG`, few
     individuals per trait. The ratones data (STATUS.md) as the cautionary
     case: longer chains did not help, a stronger prior did.
   - On these data: fit LKJ 1 and 2 (and 4?), compare `corrG`, the number
     of flagged blocks, eigenvalues, and the REML fits; a figure of `corrG`
     entries under each prior against sommer.
   - The recommendation as a short rule of thumb, with its evidence.
6. **G against P:** the phenotypic covariance after `X`; similarity
   (random skewers, correlation of entries); integration.
7. **Caveats:** shared-litter variance in `G` (decision 2); intervals are
   per entry; what `E` is.

Runtime: sommer and asreml fits live-sized; stanimal precomputed.

## 4. Vignette 2: genome scans (`indri-gwas.Rmd`)

Audience: users mapping loci for many traits.

1. **Setup:** G and E from vignette 1 (`fit` from sommer, saved), `A`,
   dosages, map.
2. **Multi-trait scan:** `indri_gwas(Y, dosage, fit = fit, A = A, X = X,
   snp_info = info)`; `delta` and what it says; `summary` (top markers,
   genomic inflation); `plot` Manhattan and QQ; the 9-df joint test.
3. **Single-trait scans:** `mode = "single_trait"`, `adjust = "effective"`
   (and `indri_effective_tests()` on these traits); comparing which markers
   each finds; when each test has the power (the patterns in STATUS.md's
   power study).
4. **Thresholds:** 353 markers, Bonferroni and permutations
   (`n_permutations`, `fdr`, q-values), with the caveat that permuting
   across relatives ignores relatedness.
5. **Dominance:** `dominance = TRUE`; `pvalues_dominance`; `a` and `d` at
   the top markers; an F3 has all three genotypes at about 1/4, 1/2, 1/4, so
   dominance is well estimated here (unlike in outbred samples).
6. **Effects:** the `a` and `d` vectors across the nine weeks at the main
   peaks, plotted as growth profiles (early versus late growth, as in the
   paper).
7. **Comparison with the original analysis:** the 32 QTL of Melo, Marroig
   and Wolf (2019) used interval mapping and seven gains; which of their
   chromosomes indri's single-marker scan recovers.

Runtime: scans take seconds; permutations a minute or two (precomputed if
long).

## 5. Vignette 3: G from locus effects (`indri-locus-G.Rmd`)

Audience: users connecting mapped loci to the genetic covariance.

1. **The formulas** (Eqs. 7-14 of the paper) in words: average effects,
   2pq and LD terms, additive and dominance matrices.
2. **Choosing loci:** the significant markers of vignette 2, one per peak
   (marginal effects of linked markers overlap); the LD term
   (`ld = TRUE`, composite LD from the genotypes) handles what remains.
3. **`indri_locus_G(scan, genotypes = dosage)`** on those loci: `G`, `D`,
   `alpha`.
4. **Comparisons:**
   - locus `G` against the polygenic `G` of vignette 1 (random skewers,
     entries, eigenvectors): how much of the genetic covariance structure
     the mapped loci reproduce;
   - `D` against `G`: is dominance variance organised like additive;
   - with and without LD.
5. **Effects from other software:** the conventions section of the help
   page in practice: a data frame with a SNP column, `traits_in = "rows"`,
   checking the effect allele.
6. **Caveats:** single-marker effects are marginal and inflated at the
   selected peaks (winner's curse); a joint multi-locus model (the paper's
   Bayesian regression) gives better effects.

Runtime: seconds.

## 6. Vignette 4: genotype-by-environment scans (`indri-gxe.Rmd`)

Audience: users asking whether loci act differently between environments.
Sex is the environment: 191 of the 200 families have both sexes, so
cross-sex genetic covariance is estimable.

1. **Does growth have GxE in the background?** `fit_G_env(Y, env = SEX,
   X, A)`: within-sex `G`, cross-sex genetic correlations per week
   (`corr_between`), and what values below 1 mean.
2. **Why the background matters:** the calibration result (STATUS.md):
   unmodelled polygenic G x sex inflates interaction tests; the
   environment-specific background removes it (pending the current
   meyer run).
3. **The scan:** `indri_gwas(..., fit = fit_env, gxe = "SEXM")` (the
   column name from `model.matrix`); `pvalues_gxe`, `gxe` effects;
   Manhattan of interactions; the same scan with the usual background for
   comparison (genomic inflation of each).
4. **Reading interactions:** sex-specific effects at the top markers
   (`a` for females, `a + gxe` for males), plotted across weeks.
5. **Cost:** the environment-specific scan is about 7 times the usual one
   on family data; why (no single rotation; one decomposition per family
   cluster).
6. **Caveats:** `fit_G_env()` is experimental; cross-sex covariance rests
   on opposite-sex sibs; dominance is a main effect only.

Runtime: `fit_G_env()` with sommer on 18 character states and the scans:
minutes; precomputed.

## 7. Infrastructure

- `DESCRIPTION`: `Suggests: knitr, rmarkdown` (ggplot2 already),
  `VignetteBuilder: knitr`.
- `vignettes/`: four `.Rmd.orig` sources, the generated `.Rmd`, figures
  (`vignettes/figures/`), `precompute.R`; `.Rbuildignore` the `.orig` and
  the script.
- Saved fits for reuse between vignettes (`vignettes/cache/*.rds`,
  build-ignored), so vignette 2 does not refit vignette 1.
- Each vignette states its runtime and which engines it needs, and skips
  gracefully (precomputed output shown) when sommer, asreml or CmdStan is
  missing.
- `R CMD check` must stay at Status OK with the vignettes built.

## 8. Order of work

1. Decisions 1-4 above; extend the data if 9 gains (with a test).
2. Infrastructure and the shared setup chunk.
3. Vignette 1 (sommer and asreml first, then stanimal and the LKJ runs on
   meyer, with BLAS capped).
4. Vignette 2, from vignette 1's saved fit.
5. Vignette 3.
6. Vignette 4, after the GxE calibration with the environment-specific
   background is in (running now).
7. README: a "Vignettes" section; pkgdown later, if wanted.

## 9. Open questions

- 9 weekly gains, the paper's 7, or the 10 raw weights?
- Fixed effects as in the original analyses (sex, cohort, litter sizes),
  and accept shared-litter variance in `G`, or wait for a litter random
  effect?
- Precomputed `.Rmd.orig` vignettes, or saved fits with live code?
- evolqg in `Suggests:`, or helpers in the vignettes?
