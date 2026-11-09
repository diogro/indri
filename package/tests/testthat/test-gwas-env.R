# indri_gwas(fit = <fit_G_env result>): scans with an environment-specific
# background covariance, against a dense covariance built entry by entry.

test_that("the multi-trait scan matches dense GLS with an environment-specific G", {
  d <- simulate_env_gwas(seed = 121)
  scan <- indri_gwas(d$Y, d$geno, fit = d$fit, A = d$A, X = d$X, delta = 0.8)
  expect_equal(scan$background, "environment")
  for (j in seq_len(d$M)) {
    g <- d$geno[, j]
    ref <- dense_env_multi(d$Y, d$X, cbind(d$X, g), d$A, d$env, d$Gs, d$Es,
                           delta = 0.8)
    expect_rel_equal(scan$pvalues[j], ref$p, 1e-8, label = paste("p", j))
    expect_equal(unname(scan$a[j, ]), ref$beta[1, ], tolerance = 1e-8)
  }
})

test_that("delta is fitted by REML on the environment-specific covariance", {
  d <- simulate_env_gwas(seed = 122)
  scan <- indri_gwas(d$Y, d$geno, fit = d$fit, A = d$A, X = d$X)
  expect_equal(scan$delta, dense_env_delta(d$Y, d$X, d$A, d$env, d$Gs, d$Es),
               tolerance = 1e-5)
  ref <- dense_env_multi(d$Y, d$X, cbind(d$X, d$geno[, 2]), d$A, d$env, d$Gs,
                         d$Es, delta = scan$delta)
  expect_rel_equal(scan$pvalues[2], ref$p, 1e-6)
})

test_that("gxe and dominance work with the environment-specific covariance", {
  d <- simulate_env_gwas(seed = 123, M = 6L)
  scan <- indri_gwas(d$Y, d$geno, fit = d$fit, A = d$A, X = d$X, delta = 1,
                     gxe = "male", dominance = TRUE)
  for (j in which(scan$dominance_estimable & scan$gxe_estimable)) {
    g <- d$geno[, j]
    h <- heterozygous(g)
    ge <- g * d$X[, "male"]
    full <- cbind(d$X, g, h, ge)
    joint <- dense_env_multi(d$Y, d$X, full, d$A, d$env, d$Gs, d$Es)
    inter <- dense_env_multi(d$Y, cbind(d$X, g, h), full, d$A, d$env, d$Gs,
                             d$Es)
    expect_rel_equal(scan$pvalues[j], joint$p, 1e-8)
    expect_rel_equal(scan$pvalues_gxe[j], inter$p, 1e-8)
    expect_equal(unname(scan$gxe[j, , "male"]), joint$beta[3, ],
                 tolerance = 1e-8)
  }
})

test_that("relatives across families use the dense path, with the same answer", {
  d <- simulate_env_gwas(seed = 124)
  # link two families: A is no longer block diagonal by family
  A <- d$A
  A["id1", "id6"] <- A["id6", "id1"] <- 0.25
  scan <- indri_gwas(d$Y, d$geno, fit = d$fit, A = A, X = d$X, delta = 1)
  ref <- dense_env_multi(d$Y, d$X, cbind(d$X, d$geno[, 3]), A, d$env, d$Gs,
                         d$Es)
  expect_rel_equal(scan$pvalues[3], ref$p, 1e-8)
  expect_lt(scan$n_components, nrow(A) / 5)   # the sibships, joined in one
})

test_that("single-trait scans use each trait's environment-specific G", {
  d <- simulate_env_gwas(seed = 125)
  scan <- indri_gwas(d$Y, d$geno, fit = d$fit, A = d$A, X = d$X,
                     mode = "single_trait", delta = c(wing = 0.7, tail = 1.2))
  for (k in seq_len(d$P)) {
    trait <- colnames(d$Y)[k]
    s <- paste(trait, c("f", "m"), sep = ":")
    for (j in 1:2) {
      ref <- dense_env_multi(d$Y[, k, drop = FALSE], d$X,
                             cbind(d$X, d$geno[, j]), d$A, d$env,
                             d$Gs[s, s, drop = FALSE],
                             d$Es[s, s, drop = FALSE],
                             delta = c(0.7, 1.2)[k])
      expect_rel_equal(scan$pvalues[j, k], ref$p, 1e-8)
    }
  }
})

test_that("an environment fit is matched to the scan and checked", {
  d <- simulate_env_gwas(seed = 126)
  # individuals matched by name
  o <- sample(d$N)
  moved <- indri_gwas(d$Y[o, ], d$geno, fit = d$fit, A = d$A, X = d$X,
                      delta = 1)
  expect_equal(moved$pvalues,
               indri_gwas(d$Y, d$geno, fit = d$fit, A = d$A, X = d$X,
                          delta = 1)$pvalues, tolerance = 1e-10)
  missing_env <- d$fit
  missing_env$env <- missing_env$env[-1]
  expect_error(indri_gwas(d$Y, d$geno, fit = missing_env, A = d$A, X = d$X),
               "environment")
  other <- d$Y
  colnames(other) <- c("a", "b")
  expect_error(indri_gwas(other, d$geno, fit = d$fit, A = d$A, X = d$X),
               "traits")
  expect_error(indri_gwas(d$Y, d$geno, fit = d$fit, X = d$X, model = "lm"),
               "lmm")
  expect_error(indri_gwas(d$Y, d$geno, fit = d$fit, X = d$X), "A")
})

test_that("with no GxE in the background it is the Kronecker scan", {
  # cross-sex genetic correlation 1 and the same E in both sexes: the
  # environment-specific covariance is G (x) A + E (x) I
  d <- simulate_env_gwas(seed = 127, r_between = 1)
  G0 <- d$Gs[1:2, 1:2]
  E0 <- d$Es[1:2, 1:2]
  dimnames(G0) <- dimnames(E0) <- list(colnames(d$Y), colnames(d$Y))
  env_scan <- indri_gwas(d$Y, d$geno, fit = d$fit, A = d$A, X = d$X,
                         delta = 0.9, gxe = "male")
  kron <- indri_gwas(d$Y, d$geno, G = G0, E = E0, A = d$A, X = d$X,
                     delta = 0.9, gxe = "male")
  expect_equal(env_scan$pvalues, kron$pvalues, tolerance = 1e-8)
  expect_equal(env_scan$pvalues_gxe, kron$pvalues_gxe, tolerance = 1e-8)
  expect_equal(env_scan$a, kron$a, tolerance = 1e-8)
})
