# Scans with extra random effects (indri_gwas(random =)): the low-rank
# (Woodbury) path against the rotation path and dense GLS.

sim_random <- function(seed = 1, N = 72, P = 2, L = 12) {
  set.seed(seed)
  fam <- rep(seq_len(L), each = N / L)
  A <- 0.5 * outer(fam, fam, "==") + diag(0.5, N)
  tr <- paste0("t", P:1)                      # names that do not sort in column order
  G <- 0.5 * diag(P) + 0.3; E <- diag(P) - 0.1 * (P > 1); C <- 0.4 * diag(P) - 0.15 * (P > 1)
  diag(E) <- 1; diag(C) <- 0.4
  dimnames(G) <- dimnames(E) <- dimnames(C) <- list(tr, tr)
  nurse <- sample(fam)
  dam <- sample(rep(seq_len(N / 4), 4))
  V <- kronecker(G, A) + kronecker(C, tcrossprod(incidence(nurse))) +
    kronecker(0.3 * C, tcrossprod(incidence(dam))) + kronecker(E, diag(N))
  Y <- matrix(t(chol(V)) %*% stats::rnorm(N * P), N, P, dimnames = list(NULL, tr))
  X <- cbind(1, stats::rnorm(N))
  g <- matrix(stats::rbinom(N * 4, 2, 0.4), N, 4, dimnames = list(NULL, paste0("s", 1:4)))
  H <- array(0, c(N, 3, 3), dimnames = list(NULL, c("f1", "f2", "f3"), paste0("p", 1:3)))
  for (m in 1:3) for (i in seq_len(N)) for (a in sample(3, 2, TRUE)) H[i, a, m] <- H[i, a, m] + 1
  list(Y = Y, X = X, A = A, G = G, E = E, C = C, nurse = nurse, dam = dam, g = g, H = H)
}

test_that("a zero random covariance reproduces the rotation path", {
  s <- sim_random()
  zero <- list(nurse = s$C * 0)
  for (mode in c("multi_trait", "single_trait")) {
    r0 <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X, mode = mode)
    rz <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X, mode = mode,
                     random = list(nurse = s$nurse), random_cov = zero,
                     delta = if (mode == "multi_trait") r0$delta else rep(1, 2))
    if (mode == "multi_trait") {
      expect_equal(rz$pvalues, r0$pvalues, tolerance = 1e-8)
      expect_equal(rz$a, r0$a, tolerance = 1e-8)
    }
  }
  # REML delta agrees too
  r0 <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X)
  rz <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X,
                   random = list(nurse = s$nurse), random_cov = zero)
  expect_equal(rz$delta, r0$delta, tolerance = 1e-5)
  expect_equal(rz$pvalues, r0$pvalues, tolerance = 1e-4)
  h0 <- indri_gwas_haplotypes(s$Y, s$H, G = s$G, E = s$E, A = s$A, X = s$X, delta = 0.8)
  hz <- indri_gwas_haplotypes(s$Y, s$H, G = s$G, E = s$E, A = s$A, X = s$X, delta = 0.8,
                              random = list(nurse = s$nurse), random_cov = zero)
  expect_equal(hz$pvalues, h0$pvalues, tolerance = 1e-8)
  expect_equal(hz$effects, h0$effects, tolerance = 1e-8)
})

test_that("SNP scans with a random term match dense GLS", {
  s <- sim_random(2)
  Z <- list(incidence(s$nurse)); Cs <- list(s$C)
  r <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X,
                  random = list(nurse = s$nurse), random_cov = list(nurse = s$C))
  for (j in seq_len(ncol(s$g))) {
    d <- dense_random(s$Y, s$X, s$g[, j, drop = FALSE], s$G, s$E, s$A, Z, Cs)
    expect_equal(unname(r$pvalues[j]), d$p, tolerance = 1e-6)
    expect_equal(unname(r$a[j, ]), d$beta[1, ], tolerance = 1e-6)
  }
  expect_equal(r$delta, d$delta, tolerance = 1e-5)
  expect_equal(r$random_terms, "nurse")
})

test_that("single-trait scans with a random term match dense GLS", {
  s <- sim_random(3)
  r <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X, mode = "single_trait",
                  random = list(nurse = s$nurse), random_cov = list(nurse = s$C))
  rf <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X, mode = "single_trait",
                   single_trait_variance = "fit",
                   random = list(nurse = s$nurse), random_cov = list(nurse = s$C))
  expect_equal(unname(rf$delta), c(1, 1))
  Z <- list(incidence(s$nurse))
  for (c in 1:2) for (j in 1:2) {
    one <- function(M) M[c, c, drop = FALSE]
    d <- dense_random(s$Y[, c, drop = FALSE], s$X, s$g[, j, drop = FALSE], one(s$G), one(s$E),
                      s$A, Z, list(one(s$C)))
    expect_equal(unname(r$pvalues[j, c]), d$p, tolerance = 1e-6)
    expect_equal(unname(r$delta[c]), d$delta, tolerance = 1e-5)
    dfit <- dense_random(s$Y[, c, drop = FALSE], s$X, s$g[, j, drop = FALSE], one(s$G), one(s$E),
                         s$A, Z, list(one(s$C)), delta = 1)
    expect_equal(unname(rf$pvalues[j, c]), dfit$p, tolerance = 1e-6)
  }
})

test_that("haplotype scans with two random terms match dense GLS", {
  s <- sim_random(4)
  Cd <- 0.3 * s$C
  r <- indri_gwas_haplotypes(s$Y, s$H, G = s$G, E = s$E, A = s$A, X = s$X,
                             random = list(nurse = s$nurse, dam = s$dam),
                             random_cov = list(nurse = s$C, dam = Cd))
  kept <- setdiff(dimnames(s$H)[[2]], r$baseline)
  for (m in 1:3) {
    d <- dense_random(s$Y, s$X, s$H[, kept, m], s$G, s$E, s$A,
                      list(incidence(s$nurse), incidence(s$dam)), list(s$C, Cd))
    expect_equal(unname(r$pvalues[m]), d$p, tolerance = 1e-6)
    eff <- rbind(0, d$beta)                     # baseline then kept, by trait
    rownames(eff) <- c(r$baseline, kept)
    eff <- sweep(eff[dimnames(s$H)[[2]], ], 2, colMeans(eff))
    expect_equal(unname(r$effects[m, , ]), unname(eff), tolerance = 1e-6)
  }
  expect_equal(r$df, 2L * 2L)
})

test_that("the covariances come from fit$random by default", {
  s <- sim_random(5)
  fit <- list(G = s$G, E = s$E, random = list(nurse = s$C))
  a <- indri_gwas(s$Y, s$g, fit = fit, A = s$A, X = s$X, random = list(nurse = s$nurse))
  b <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X,
                  random = list(nurse = s$nurse), random_cov = list(nurse = s$C))
  expect_equal(a$pvalues, b$pvalues)
})

test_that("missing levels are unshared and monomorphic SNPs get NA", {
  s <- sim_random(6)
  nurse <- s$nurse; nurse[1:3] <- NA
  g <- s$g; g[, 2] <- 1
  r <- suppressWarnings(indri_gwas(s$Y, g, G = s$G, E = s$E, A = s$A, X = s$X,
                                   random = list(nurse = nurse), random_cov = list(nurse = s$C),
                                   delta = 1))
  Zn <- incidence(ifelse(is.na(nurse), paste0("u", seq_along(nurse)), nurse))
  d <- dense_random(s$Y, s$X, g[, 1, drop = FALSE], s$G, s$E, s$A, list(Zn), list(s$C), delta = 1)
  expect_equal(unname(r$pvalues[1]), d$p, tolerance = 1e-6)
  expect_true(is.na(r$pvalues[2]))
})

test_that("random effects are checked", {
  s <- sim_random(7)
  run <- function(...) indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X, ...)
  expect_error(run(random = list(nurse = s$nurse)), "covariance of each term")
  expect_error(run(random = list(s$nurse), random_cov = list(nurse = s$C)), "named list")
  expect_error(run(random = list(nurse = s$nurse[-1]), random_cov = list(nurse = s$C)), "values")
  expect_error(run(random = list(nurse = s$nurse), random_cov = list(dam = s$C)), "no covariance")
  expect_error(run(random = list(nurse = s$nurse), random_cov = list(nurse = -s$C)), "semi-definite")
  expect_error(run(random = list(nurse = s$nurse), random_cov = list(nurse = s$C), dominance = TRUE),
               "not available")
  expect_error(indri_gwas(s$Y, s$g, X = s$X, model = "lm", random = list(nurse = s$nurse),
                          random_cov = list(nurse = s$C)), "lmm")
  expect_error(indri_gwas(s$Y, s$g, A = s$A, X = s$X, mode = "single_trait",
                          random = list(nurse = s$nurse), random_cov = list(nurse = s$C)),
               "need `G` and `E`")
})

test_that("indri_gls() matches dense GLS, with and without random terms", {
  s <- sim_random(8)
  Xf <- cbind(s$X, s$g[, 1:2])
  Z <- incidence(s$nurse)
  for (with_random in c(TRUE, FALSE)) {
    r <- if (with_random) {
      indri_gls(s$Y, Xf, G = s$G, E = s$E, A = s$A, random = list(nurse = s$nurse),
                random_cov = list(nurse = s$C), delta = 0.9)
    } else indri_gls(s$Y, Xf, G = s$G, E = s$E, A = s$A, delta = 0.9)
    V <- kronecker(s$G, s$A) + 0.9 * kronecker(s$E, diag(nrow(s$Y))) +
      if (with_random) kronecker(s$C, tcrossprod(Z)) else 0
    d <- dense_gls(c(s$Y), kronecker(diag(2), Xf), solve(V))
    expect_equal(c(r$coefficients), d$beta, tolerance = 1e-8)
    expect_equal(unname(r$vcov), unname(solve(d$A)), tolerance = 1e-8)
  }
  # REML delta as in the scans
  r <- indri_gls(s$Y, s$X, G = s$G, E = s$E, A = s$A, random = list(nurse = s$nurse),
                 random_cov = list(nurse = s$C))
  sc <- indri_gwas(s$Y, s$g, G = s$G, E = s$E, A = s$A, X = s$X,
                   random = list(nurse = s$nurse), random_cov = list(nurse = s$C))
  expect_equal(r$delta, sc$delta, tolerance = 1e-6)
})
