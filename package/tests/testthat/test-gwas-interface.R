# The contract of indri_gwas() around the computation: sample and trait
# matching, input checks, genotype handling, the fit path and the outputs.

scan_lmm <- function(d, ...) {
  indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X, ...)
}

test_that("samples are matched by name across Y, genotypes, A and X", {
  d <- simulate_gwas_data(seed = 31)
  ref <- scan_lmm(d, delta = 1)
  o <- sample(d$N)
  shuffled <- indri_gwas(d$Y, d$geno[o, ], G = d$G, E = d$E,
                         A = d$A[rev(o), rev(o)], X = d$X[o, ], delta = 1)
  expect_equal(shuffled$pvalues, ref$pvalues, tolerance = 1e-10)
  expect_equal(shuffled$a, ref$a, tolerance = 1e-10)

  bad <- d$geno
  rownames(bad)[1] <- "nobody"
  expect_error(indri_gwas(d$Y, bad, G = d$G, E = d$E, A = d$A, X = d$X),
               "genotypes")
})

test_that("unnamed inputs are matched by position and must agree in size", {
  d <- simulate_gwas_data(seed = 32)
  named <- scan_lmm(d, delta = 1)
  plain <- indri_gwas(unname(d$Y), unname(d$geno), G = unname(d$G),
                      E = unname(d$E), A = unname(d$A), X = unname(d$X),
                      delta = 1)
  expect_equal(unname(plain$pvalues), unname(named$pvalues), tolerance = 1e-10)
  expect_error(indri_gwas(unname(d$Y), unname(d$geno)[-1, ], G = d$G,
                          E = d$E, A = unname(d$A), X = unname(d$X)),
               "rows")
})

test_that("G and E are matched to the traits of Y by name", {
  d <- simulate_gwas_data(seed = 33)
  ref <- scan_lmm(d, delta = 1)
  o <- c(3L, 1L, 2L)
  moved <- indri_gwas(d$Y, d$geno, G = d$G[o, o], E = d$E[o, o], A = d$A,
                      X = d$X, delta = 1)
  expect_equal(moved$pvalues, ref$pvalues, tolerance = 1e-10)
  expect_equal(colnames(moved$a), colnames(d$Y))
  wrong <- d$G
  dimnames(wrong) <- list(letters[1:3], letters[1:3])
  expect_error(indri_gwas(d$Y, d$geno, G = wrong, E = d$E, A = d$A, X = d$X),
               "traits")
})

test_that("an indri fit supplies G and E", {
  d <- simulate_gwas_data(seed = 34)
  fit <- structure(list(G = d$G, E = d$E), class = "indri_fit")
  via_fit <- indri_gwas(d$Y, d$geno, fit = fit, A = d$A, X = d$X, delta = 1)
  expect_equal(via_fit$pvalues, scan_lmm(d, delta = 1)$pvalues,
               tolerance = 1e-12)
  expect_error(indri_gwas(d$Y, d$geno, fit = fit, G = d$G, A = d$A, X = d$X),
               "either")
  # a fit made before E was assembled cannot drive the multi-trait test
  old <- structure(list(G = d$G), class = "indri_fit")
  expect_error(indri_gwas(d$Y, d$geno, fit = old, A = d$A, X = d$X), "E")
})

test_that("what each model needs is checked", {
  d <- simulate_gwas_data(seed = 35)
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = d$E, X = d$X), "A")
  expect_error(indri_gwas(d$Y, d$geno, A = d$A, X = d$X), "G")
  expect_error(indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                          single_trait_variance = "fit"), "G")
  # the single-trait REML path needs neither G nor E
  expect_silent(indri_gwas(d$Y, d$geno, A = d$A, X = d$X,
                           mode = "single_trait"))
  expect_error(scan_lmm(d, delta = -1), "delta")
  expect_error(indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                          delta = c(1, 2)), "delta")
  expect_error(scan_lmm(d, adjust = "bonferroni"), "single")
  not_pd <- d$E
  not_pd[1, 2] <- not_pd[2, 1] <- 5
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = not_pd, A = d$A, X = d$X),
               "positive definite")
})

test_that("Y must be complete and X must be a usable design", {
  d <- simulate_gwas_data(seed = 36)
  gap <- d$Y
  gap[2, 1] <- NA
  expect_error(indri_gwas(gap, d$geno, G = d$G, E = d$E, A = d$A, X = d$X),
               "missing")
  rank_deficient <- cbind(d$X, twice = 2 * d$X[, "x"])
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A,
                          X = rank_deficient), "rank")
  # the default design is an intercept
  default <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, delta = 1)
  explicit <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, delta = 1,
                         X = matrix(1, d$N, 1, dimnames = list(rownames(d$Y),
                                                               "(Intercept)")))
  expect_equal(default$pvalues, explicit$pvalues, tolerance = 1e-12)
})

test_that("missing genotypes are mean-imputed and monomorphic SNPs give NA", {
  d <- simulate_gwas_data(seed = 37)
  geno <- d$geno
  geno[c(1, 5), 2] <- NA
  imputed <- geno
  imputed[c(1, 5), 2] <- mean(geno[, 2], na.rm = TRUE)
  expect_equal(scan_lmm(modifyList(d, list(geno = geno)), delta = 1)$pvalues,
               scan_lmm(modifyList(d, list(geno = imputed)), delta = 1)$pvalues,
               tolerance = 1e-12)

  geno <- d$geno
  geno[, 3] <- 1
  expect_warning(scan <- scan_lmm(modifyList(d, list(geno = geno)), delta = 1),
                 "monomorphic")
  expect_true(is.na(scan$pvalues[3]))
  expect_true(all(is.na(scan$a[3, ])))
  expect_false(anyNA(scan$pvalues[-3]))
})

test_that("SNP information must describe the genotype columns", {
  d <- simulate_gwas_data(seed = 38)
  info <- data.frame(snp = colnames(d$geno), chrom = 1L,
                     pos = seq_len(d$M) * 100L)
  scan <- scan_lmm(d, delta = 1, snp_info = info)
  expect_equal(scan$snp_info, info)
  expect_equal(scan_lmm(d, delta = 1)$snp_info$snp, colnames(d$geno))
  expect_error(scan_lmm(d, delta = 1, snp_info = info[-1, ]), "snp_info")
})

test_that("results are named and convert to a data frame", {
  d <- simulate_gwas_data(seed = 39)
  multi <- scan_lmm(d, delta = 1)
  expect_named(multi$pvalues, colnames(d$geno))
  expect_equal(dimnames(multi$a), list(colnames(d$geno), colnames(d$Y)))
  expect_equal(multi$mode, "multi_trait")
  expect_equal(multi$model, "lmm")
  df <- as.data.frame(multi)
  expect_equal(nrow(df), d$M)
  expect_equal(df$pvalue, unname(multi$pvalues))
  expect_true(all(paste0("a_", colnames(d$Y)) %in% names(df)))

  single <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                       adjust = "bonferroni")
  long <- as.data.frame(single)
  expect_equal(nrow(long), d$M * d$P)
  expect_true(all(c("snp", "trait", "pvalue", "a", "pvalue_adjust") %in%
                    names(long)))
  expect_equal(long$pvalue[long$trait == "tail"],
               unname(single$pvalues[, "tail"]))
  expect_output(print(multi), "multi-trait")
  expect_output(print(single), "single-trait")
})

test_that("Bonferroni and effective adjustments multiply and cap at 1", {
  d <- simulate_gwas_data(seed = 40)
  raw <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait")
  expect_null(raw$pvalues_adjust)
  bonf <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                     adjust = "bonferroni")
  expect_equal(bonf$n_tests, d$P)
  expect_equal(bonf$pvalues_adjust, pmin(raw$pvalues * d$P, 1))
  eff <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                    adjust = "effective")
  expect_equal(eff$n_tests, indri_effective_tests(d$Y, d$X))
  expect_equal(eff$pvalues_adjust, pmin(raw$pvalues * eff$n_tests, 1))
})

test_that("the effective number of tests runs from 1 to the number of traits", {
  set.seed(41)
  z <- stats::rnorm(200)
  same <- cbind(z, z + 1e-8 * stats::rnorm(200), z * 2)
  expect_equal(indri_effective_tests(same), 1, tolerance = 1e-4)
  independent <- matrix(stats::rnorm(5000 * 4), 5000, 4)
  expect_equal(indri_effective_tests(independent), 4, tolerance = 0.05)
})

test_that("a covariate shared by the traits does not lower the effective count", {
  set.seed(42)
  n <- 2000
  sex <- stats::rbinom(n, 1, 0.5)
  Y <- matrix(stats::rnorm(n * 3), n, 3) + 3 * sex
  X <- cbind(intercept = 1, sex = sex)
  expect_lt(indri_effective_tests(Y), 2.5)  # measured 2.33 (correlations 0.7)
  expect_equal(indri_effective_tests(Y, X), 3, tolerance = 0.05)
})
