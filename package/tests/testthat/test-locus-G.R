# indri_locus_G(): additive and dominance covariance matrices from per-locus
# effects (Melo, Marroig and Wolf 2019, Eqs. 7-14).

test_that("one locus gives 2pq alpha alpha' and (2pq)^2 d d'", {
  a <- matrix(c(0.5, -0.2, 0.1), 1, dimnames = list("rs1", c("x", "y", "z")))
  d <- matrix(c(0.3, 0.1, -0.4), 1)
  out <- indri_locus_G(a, d, p = 0.3)
  alpha <- a + (0.7 - 0.3) * d
  expect_equal(unname(out$G), unname(2 * 0.3 * 0.7 * crossprod(alpha)))
  expect_equal(unname(out$D), (2 * 0.3 * 0.7)^2 * crossprod(d))
  expect_equal(out$full_sib, out$G / 2 + out$D / 4)
  expect_equal(dimnames(out$G), list(c("x", "y", "z"), c("x", "y", "z")))
  expect_equal(unname(out$alpha), unname(alpha))
})

test_that("linkage disequilibrium enters as in Eqs. 9-12", {
  set.seed(1)
  M <- 4L; P <- 3L
  a <- matrix(stats::rnorm(M * P), M, P)
  d <- matrix(stats::rnorm(M * P), M, P)
  p <- c(0.2, 0.5, 0.4, 0.7); q <- 1 - p
  lambda <- matrix(0, M, M)
  lambda[1, 2] <- lambda[2, 1] <- 0.05
  lambda[3, 4] <- lambda[4, 3] <- -0.03
  out <- indri_locus_G(a, d, p = p, lambda = lambda)
  alpha <- a + (q - p) * d
  G <- D <- matrix(0, P, P)
  for (k in seq_len(P)) for (l in seq_len(P)) {
    G[k, l] <- sum(2 * p * q * alpha[, k] * alpha[, l])
    D[k, l] <- sum((2 * p * q)^2 * d[, k] * d[, l])
    for (j in seq_len(M)) for (x in seq_len(M)) if (j != x) {
      G[k, l] <- G[k, l] + 2 * lambda[j, x] * alpha[j, k] * alpha[x, l]
      D[k, l] <- D[k, l] + 4 * lambda[j, x]^2 * d[j, k] * d[x, l]
    }
  }
  expect_equal(unname(out$G), G)
  expect_equal(unname(out$D), D)
  # without LD only the per-locus sums remain
  none <- indri_locus_G(a, d, p = p, lambda = lambda, ld = FALSE)
  expect_equal(unname(none$G), unname(crossprod(alpha, 2 * p * q * alpha)))
  expect_equal(unname(none$D), unname(crossprod(d, (2 * p * q)^2 * d)))
})

test_that("frequencies and LD come from a 0/1/2 genotype matrix", {
  set.seed(2)
  N <- 400L
  g1 <- stats::rbinom(N, 2L, 0.4)
  g2 <- ifelse(stats::runif(N) < 0.8, g1, stats::rbinom(N, 2L, 0.4))  # linked
  g3 <- stats::rbinom(N, 2L, 0.6)
  geno <- cbind(rs1 = g1, rs2 = g2, rs3 = g3)
  a <- matrix(c(0.4, 0.2, -0.3, 0.1, 0, 0.5), 3,
              dimnames = list(c("rs1", "rs2", "rs3"), c("x", "y")))
  out <- indri_locus_G(a, genotypes = geno)
  p <- colMeans(geno) / 2
  expect_equal(unname(out$p), unname(p))
  lambda <- stats::cov(geno) / 2           # composite LD, Cov(dosages) / 2
  diag(lambda) <- 0
  S <- 2 * lambda
  diag(S) <- 2 * p * (1 - p)
  expect_equal(unname(out$G), unname(crossprod(a, S %*% a)))
  # loci are matched to the genotype columns by name
  shuffled <- indri_locus_G(a, genotypes = geno[, c(3, 1, 2)])
  expect_equal(shuffled$G, out$G)
  # no dominance effects, no dominance matrix
  expect_equal(unname(out$D), matrix(0, 2, 2))
})

test_that("G + D is the covariance of the genotypic values (unlinked loci)", {
  # Independent loci in Hardy-Weinberg equilibrium: the variance of the
  # genotypic values is V_A + V_D exactly, in expectation.
  set.seed(3)
  N <- 200000L
  p <- c(0.3, 0.5, 0.15)
  geno <- vapply(p, function(f) stats::rbinom(N, 2L, f), integer(N))
  a <- matrix(c(0.5, -0.3, 0.2, 0.1, 0.4, -0.2), 3)
  d <- matrix(c(0.4, 0.2, -0.3, -0.1, 0.3, 0.6), 3)
  values <- (geno - 1) %*% a + (geno == 1) %*% d
  out <- indri_locus_G(a, d, p = p)
  expect_equal(unname(out$G + out$D), unname(stats::cov(values)),
               tolerance = 0.02)
})

test_that("an indri_gwas scan supplies a and d", {
  d <- simulate_gwas_data(M = 6L, seed = 91)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 1, dominance = TRUE)
  from_scan <- indri_locus_G(scan, genotypes = d$geno)
  ok <- scan$dominance_estimable
  direct <- indri_locus_G(scan$a[ok, ], scan$d[ok, ],
                          genotypes = d$geno[, ok])
  expect_equal(from_scan$G, direct$G)
  expect_equal(from_scan$D, direct$D)
  expect_equal(from_scan$loci, colnames(d$geno)[ok])
  # an additive scan gives G only
  additive <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                         delta = 1)
  expect_equal(unname(indri_locus_G(additive, genotypes = d$geno)$D),
               matrix(0, d$P, d$P))
})

test_that("loci without effects are dropped, and inputs are checked", {
  a <- matrix(c(0.5, NA, 0.2, 0.1, NA, 0.3), 3,
              dimnames = list(c("rs1", "rs2", "rs3"), NULL))
  expect_message(out <- indri_locus_G(a, p = c(0.3, 0.4, 0.5)), "1 locus")
  expect_equal(out$loci, c("rs1", "rs3"))
  expect_error(indri_locus_G(a[1:2, ], p = c(0.3, 0.4, 0.5)), "p")
  expect_error(indri_locus_G(matrix(1, 2, 2)), "genotypes")
  expect_error(indri_locus_G(matrix(1, 2, 2), p = c(0.3, 1.2)), "p")
  expect_error(indri_locus_G(matrix(1, 2, 2), d = matrix(1, 3, 2),
                             p = c(0.3, 0.4)), "d")
})
