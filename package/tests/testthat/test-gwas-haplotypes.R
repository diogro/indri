# indri_gwas_haplotypes(): founder dosages tested jointly at each position,
# checked against dense GLS references (helper-gwas.R).

sim_hap <- function(N = 70, Fn = 4, M = 6, P = 2, seed = 3) {
  set.seed(seed)
  founders <- paste0("F", seq_len(Fn))
  H <- array(0, c(N, Fn, M), dimnames = list(paste0("i", seq_len(N)), founders, paste0("pos", seq_len(M))))
  for (m in seq_len(M)) for (i in seq_len(N)) {
    a <- sample(Fn, 2, replace = TRUE)
    H[i, a[1], m] <- H[i, a[1], m] + 1; H[i, a[2], m] <- H[i, a[2], m] + 1
  }
  Z <- matrix(rnorm(N * 30), N, 30)
  A <- tcrossprod(scale(Z)) / 30 + diag(0.05, N)
  dimnames(A) <- list(dimnames(H)[[1]], dimnames(H)[[1]])
  G <- matrix(c(1, 0.4, 0.4, 0.8), 2)[seq_len(P), seq_len(P), drop = FALSE]
  E <- matrix(c(1, -0.2, -0.2, 1.2), 2)[seq_len(P), seq_len(P), drop = FALSE]
  U <- t(chol(A)) %*% matrix(rnorm(N * P), N, P) %*% chol(G)
  Y <- U + matrix(rnorm(N * P), N, P) %*% chol(E)
  Y <- Y + H[, , 2] %*% matrix(c(0.8, 0, -0.5, 0, 0.6, 0.2, -0.3, 0)[seq_len(Fn * P)], Fn, P)
  dimnames(Y) <- list(dimnames(H)[[1]], paste0("t", seq_len(P)))
  names_tr <- colnames(Y); dimnames(G) <- dimnames(E) <- list(names_tr, names_tr)
  X <- cbind(1, rnorm(N))
  list(Y = Y, H = H, A = A, G = G, E = E, X = X)
}

# dense reference: GLS with the founder dosages (all but `drop`) for every trait
dense_hap_multi <- function(Y, X, Hm, G, E, A, delta) {
  N <- nrow(Y); P <- ncol(Y)
  y <- c(Y)
  D0 <- kronecker(diag(P), X)
  D1 <- cbind(D0, kronecker(diag(P), Hm))
  Vi <- solve(kronecker(G, A) + delta * kronecker(E, diag(N)))
  f0 <- dense_gls(y, D0, Vi); f1 <- dense_gls(y, D1, Vi)
  stat <- N * P * log(f0$q / f1$q)
  list(p = stats::pchisq(stat, P * ncol(Hm), lower.tail = FALSE),
       beta = matrix(utils::tail(f1$beta, P * ncol(Hm)), ncol(Hm), P))
}
dense_hap_single <- function(y, X, Hm, A, delta) {
  N <- length(y)
  Vi <- solve(A + delta * diag(N))
  f0 <- dense_gls(y, X, Vi); f1 <- dense_gls(y, cbind(X, Hm), Vi)
  stats::pchisq(N * log(f0$q / f1$q), ncol(Hm), lower.tail = FALSE)
}

test_that("the multi-trait haplotype test matches a dense GLS", {
  d <- sim_hap()
  s <- indri_gwas_haplotypes(d$Y, d$H, G = d$G, E = d$E, A = d$A, X = d$X, delta = 0.9,
                             baseline = "F1")
  expect_s3_class(s, "indri_gwas")
  expect_equal(s$df, 2L * 3L)
  for (m in 1:6) {
    ref <- dense_hap_multi(d$Y, d$X, d$H[, -1, m], d$G, d$E, d$A, 0.9)
    expect_equal(unname(s$pvalues[m]), ref$p, tolerance = 1e-6)
    eff <- rbind(0, ref$beta); eff <- sweep(eff, 2, colMeans(eff))
    expect_equal(unname(s$effects[m, , ]), unname(eff), tolerance = 1e-6)
  }
  expect_lt(s$pvalues[2], 1e-4)          # the simulated QTL
})

test_that("the single-trait haplotype test matches a dense GLS", {
  d <- sim_hap()
  s <- indri_gwas_haplotypes(d$Y, d$H, A = d$A, X = d$X, mode = "single_trait",
                             delta = c(t1 = 1.1, t2 = 0.7), baseline = "F2")
  expect_equal(dim(s$pvalues), c(6L, 2L))
  for (m in c(1, 2, 5)) for (k in 1:2)
    expect_equal(unname(s$pvalues[m, k]),
                 dense_hap_single(d$Y[, k], d$X, d$H[, -2, m], d$A, c(1.1, 0.7)[k]),
                 tolerance = 1e-6)
})

test_that("the baseline founder changes neither the test nor the effects", {
  d <- sim_hap()
  a <- indri_gwas_haplotypes(d$Y, d$H, G = d$G, E = d$E, A = d$A, X = d$X, delta = 1, baseline = "F1")
  b <- indri_gwas_haplotypes(d$Y, d$H, G = d$G, E = d$E, A = d$A, X = d$X, delta = 1, baseline = "F3")
  expect_equal(a$pvalues, b$pvalues, tolerance = 1e-8)
  expect_equal(a$effects, b$effects, tolerance = 1e-8)
})

test_that("delta is fitted by REML when not given, and the methods work", {
  d <- sim_hap()
  s <- indri_gwas_haplotypes(d$Y, d$H, G = d$G, E = d$E, A = d$A, X = d$X)
  expect_gt(s$delta, 0)
  expect_s3_class(summary(s), "summary.indri_gwas")
  expect_output(print(s), "founder haplotypes")
  expect_equal(nrow(as.data.frame(s)), 6L)
})

test_that("an absent founder at a position gives NA, and inputs are checked", {
  d <- sim_hap()
  H <- d$H; H[, "F4", 3] <- 0; H[, "F1", 3] <- 2 - rowSums(H[, -1, 3])
  s <- indri_gwas_haplotypes(d$Y, H, G = d$G, E = d$E, A = d$A, X = d$X, delta = 1, baseline = "F1")
  expect_true(is.na(s$pvalues[3]))
  expect_false(anyNA(s$pvalues[-3]))
  expect_error(indri_gwas_haplotypes(d$Y, d$H[, , 1], A = d$A), "array")
  expect_error(indri_gwas_haplotypes(d$Y, d$H, A = d$A), "fit")
  expect_error(indri_gwas_haplotypes(d$Y, d$H, G = d$G, E = d$E, A = d$A, baseline = "nope"), "baseline")
})
