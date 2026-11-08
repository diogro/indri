# indri_gwas() against the dense reference in helper-gwas.R: the rotation
# must give the same statistics, p-values and effects as the explicit GLS.

test_that("the multi-trait LMM matches dense GLS at a given delta", {
  d <- simulate_gwas_data(seed = 11)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     mode = "multi_trait", model = "lmm", delta = 0.7)
  expect_s3_class(scan, "indri_gwas")
  expect_equal(scan$delta, 0.7)
  for (j in seq_len(d$M)) {
    ref <- dense_multi(d$Y, d$X, d$geno[, j], E = d$E, G = d$G, A = d$A,
                       delta = 0.7)
    expect_rel_equal(scan$pvalues[j], ref$p, 1e-8, label = paste("p", j))
    expect_equal(unname(scan$a[j, ]), unname(ref$beta), tolerance = 1e-8,
                 label = paste("beta", j))
  }
})

test_that("the multi-trait LMM estimates delta by REML under the null", {
  d <- simulate_gwas_data(seed = 12)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     mode = "multi_trait", model = "lmm")
  ref <- dense_multi(d$Y, d$X, d$geno[, 1], E = d$E, G = d$G, A = d$A)
  expect_equal(scan$delta, ref$delta, tolerance = 1e-5)
  # and the scan uses that delta for every SNP
  for (j in seq_len(d$M)) {
    at <- dense_multi(d$Y, d$X, d$geno[, j], E = d$E, G = d$G, A = d$A,
                      delta = scan$delta)
    expect_rel_equal(scan$pvalues[j], at$p, 1e-8, label = paste("p", j))
  }
})

test_that("delta is near 1 when G and E are the matrices that made the data", {
  d <- simulate_gwas_data(n_fam = 40, size = 5, causal = 0L, seed = 13)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X)
  expect_gt(scan$delta, 0.5)
  expect_lt(scan$delta, 2)
})

test_that("the multi-trait LM uses the residual covariance after X", {
  d <- simulate_gwas_data(seed = 14)
  scan <- indri_gwas(d$Y, d$geno, X = d$X, mode = "multi_trait", model = "lm")
  resid <- stats::lm.fit(d$X, d$Y)$residuals
  E_lm <- crossprod(resid) / (nrow(resid) - ncol(d$X))
  expect_null(scan$delta)
  for (j in seq_len(d$M)) {
    ref <- dense_multi(d$Y, d$X, d$geno[, j], E = E_lm)
    expect_rel_equal(scan$pvalues[j], ref$p, 1e-8, label = paste("p", j))
    expect_equal(unname(scan$a[j, ]), unname(ref$beta), tolerance = 1e-8)
  }
  # an explicit E replaces the residual covariance
  other <- indri_gwas(d$Y, d$geno, E = d$E, X = d$X, model = "lm")
  ref <- dense_multi(d$Y, d$X, d$geno[, 2], E = d$E)
  expect_rel_equal(other$pvalues[2], ref$p, 1e-8)
})

test_that("the single-trait LMM fits delta per trait by REML", {
  d <- simulate_gwas_data(seed = 15)
  scan <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                     model = "lmm")
  expect_equal(dim(scan$pvalues), c(d$M, d$P))
  expect_named(scan$delta, colnames(d$Y))
  for (k in seq_len(d$P)) {
    ref1 <- dense_single(d$Y[, k], d$X, d$geno[, 1], A = d$A)
    expect_equal(unname(scan$delta[k]), ref1$delta, tolerance = 1e-5)
    for (j in seq_len(d$M)) {
      ref <- dense_single(d$Y[, k], d$X, d$geno[, j], A = d$A,
                          delta = scan$delta[[k]])
      expect_rel_equal(scan$pvalues[j, k], ref$p, 1e-8,
                       label = paste("p", j, k))
      expect_equal(scan$a[j, k], ref$beta, tolerance = 1e-8)
    }
  }
})

test_that("single_trait_variance = \"fit\" takes delta from E and G", {
  d <- simulate_gwas_data(seed = 16)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     mode = "single_trait", single_trait_variance = "fit")
  expect_equal(unname(scan$delta), unname(diag(d$E) / diag(d$G)))
  ref <- dense_single(d$Y[, 2], d$X, d$geno[, 3], A = d$A,
                      delta = d$E[2, 2] / d$G[2, 2])
  expect_rel_equal(scan$pvalues[3, 2], ref$p, 1e-8)
})

test_that("a given delta overrides the estimate in the single-trait LMM", {
  d <- simulate_gwas_data(seed = 17)
  given <- c(wing = 0.5, tail = 2, beak = 1)
  scan <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                     delta = given)
  expect_equal(scan$delta, given)
  ref <- dense_single(d$Y[, "tail"], d$X, d$geno[, 4], A = d$A, delta = 2)
  expect_rel_equal(scan$pvalues[4, "tail"], ref$p, 1e-8)
})

test_that("the single-trait LM is ordinary least squares", {
  d <- simulate_gwas_data(seed = 18)
  scan <- indri_gwas(d$Y, d$geno, X = d$X, mode = "single_trait",
                     model = "lm")
  for (k in seq_len(d$P)) for (j in seq_len(d$M)) {
    ref <- dense_single(d$Y[, k], d$X, d$geno[, j])
    expect_rel_equal(scan$pvalues[j, k], ref$p, 1e-8)
    ols <- stats::lm.fit(cbind(d$X, d$geno[, j]), d$Y[, k])$coefficients
    expect_equal(scan$a[j, k], unname(utils::tail(ols, 1L)),
                 tolerance = 1e-8)
  }
})

test_that("the result does not depend on the chunk size", {
  d <- simulate_gwas_data(M = 13L, seed = 19)
  a <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                  chunk_size = 1000L)
  b <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                  chunk_size = 4L)
  expect_equal(a$pvalues, b$pvalues, tolerance = 1e-12)
  expect_equal(a$a, b$a, tolerance = 1e-12)
})

test_that("a causal SNP is found by the multi-trait test", {
  d <- simulate_gwas_data(n_fam = 30, size = 5, causal = 1L, effect = 0.5,
                          M = 10L, seed = 20)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X)
  expect_equal(unname(which.min(scan$pvalues)), 1L)
  # measured 1.6e-6 at this seed
  expect_lt(scan$pvalues[[1]], 1e-4)
})
