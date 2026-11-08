# indri_gwas() against limmbo (LIMIX 1.x) on limmbo's own test data. The
# fixtures and the script that made them are in fixtures/limmbo/.
#
# LIMIX fits the null variance ratio its own way (a grid search over
# log delta in [-5, 5]; refined for the multi-trait test, not for the
# single-trait one), so these tests pass LIMIX's fitted delta and check the
# scan itself. How indri estimates delta is tested in test-gwas-dense.R.

test_that("the multi-trait LMM reproduces limmbo", {
  f <- limmbo_fixture()
  scan <- indri_gwas(f$Y, f$geno, G = f$Cg, E = f$Cn, A = f$A, X = f$X,
                     mode = "multi_trait", model = "lmm",
                     delta = f$delta$multi_lmm)
  expect_rel_equal(scan$pvalues, f$p$multi_lmm, 1e-6)
  expect_equal(unname(scan$betas), unname(f$beta$multi_lmm), tolerance = 1e-6)
  expect_equal(names(scan$pvalues), rownames(f$beta$multi_lmm))
  expect_equal(colnames(scan$betas), colnames(f$Y))
})

test_that("the multi-trait LM reproduces limmbo given its residual shape", {
  # limmbo used cov(Y) of the phenotypes as passed, with a 1e-9 genetic term
  # that LIMIX's delta pushes against its bound: agreement to ~1e-6.
  f <- limmbo_fixture()
  scan <- indri_gwas(f$Y, f$geno, E = f$Cn_lm, X = f$X,
                     mode = "multi_trait", model = "lm")
  expect_rel_equal(scan$pvalues, f$p$multi_lm, 1e-5)
  expect_equal(unname(scan$betas), unname(f$beta$multi_lm), tolerance = 1e-8)
})

test_that("the single-trait LMM reproduces limmbo", {
  f <- limmbo_fixture()
  scan <- indri_gwas(f$Y, f$geno, A = f$A, X = f$X, mode = "single_trait",
                     model = "lmm", delta = f$delta$single_lmm)
  expect_rel_equal(scan$pvalues, f$p$single_lmm, 1e-8)
  expect_equal(unname(scan$betas), unname(f$beta$single_lmm),
               tolerance = 1e-8)
})

test_that("the single-trait LM reproduces limmbo", {
  f <- limmbo_fixture()
  scan <- indri_gwas(f$Y, f$geno, X = f$X, mode = "single_trait",
                     model = "lm")
  expect_rel_equal(scan$pvalues, f$p$single_lm, 1e-8)
  expect_equal(unname(scan$betas), unname(f$beta$single_lm), tolerance = 1e-8)
})

test_that("the effective number of tests and the adjustments match limmbo", {
  f <- limmbo_fixture()
  # limmbo counted from the raw phenotypes; indri_gwas() uses the residuals
  # after X, which indri_effective_tests() gives with X.
  expect_equal(indri_effective_tests(f$Y), f$effective_tests, tolerance = 1e-10)
  bonf <- indri_gwas(f$Y, f$geno, A = f$A, X = f$X, mode = "single_trait",
                     delta = f$delta$single_lmm, adjust = "bonferroni")
  expect_equal(unname(bonf$pvalues_adjust), unname(f$padjust$lmm_bonferroni),
               tolerance = 1e-8)
  eff <- indri_gwas(f$Y, f$geno, A = f$A, X = f$X, mode = "single_trait",
                    delta = f$delta$single_lmm, adjust = "effective")
  expect_equal(eff$n_tests, indri_effective_tests(f$Y, f$X), tolerance = 1e-12)
  expect_equal(unname(eff$pvalues_adjust),
               unname(pmin(f$p$single_lmm * eff$n_tests, 1)), tolerance = 1e-8)
})

test_that("limmbo's SNP information comes through", {
  f <- limmbo_fixture()
  scan <- indri_gwas(f$Y, f$geno, X = f$X, model = "lm",
                     snp_info = f$snp_info)
  expect_equal(scan$snp_info$snp, colnames(f$geno))
  out <- as.data.frame(scan)
  expect_equal(nrow(out), ncol(f$geno))
  expect_true(all(c("snp", "chrom", "pos", "pvalue") %in% names(out)))
})
