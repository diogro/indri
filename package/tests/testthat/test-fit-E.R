# fit_G() returns the residual covariance E for indri_gwas(): the
# covariance of the residuals of Y on X minus the assembled G (scaled for
# the relatedness in A), floored to stay well conditioned. E_bar is the
# average of the block E, kept for reference. The block average is not used
# for E: at 30 traits its projection had eigenvalues near zero, and the
# multi-trait test ran at 66 times its nominal rate at 1e-4 (STATUS.md).

residual_E <- function(Y_std, X, G_std, A, floor = 1e-3) {
  R <- stats::lm.fit(X, Y_std)$residuals
  M <- diag(nrow(X)) - X %*% solve(crossprod(X), t(X))
  c_A <- sum(diag(M %*% A)) / (nrow(X) - ncol(X))
  E <- crossprod(R) / (nrow(X) - ncol(X)) - c_A * G_std
  e <- eigen(E, symmetric = TRUE)
  e$vectors %*% (pmax(e$values, floor * e$values[1]) * t(e$vectors))
}

test_that("E is the residual covariance minus G, in the units of Y", {
  skip_if_not_installed("sommer")
  data <- simulate_block_data(n_fam = 30, size = 5, K = 4L)
  Y <- sweep(data$Y, 2L, c(1, 50, 0.02, 3), "*") + 100
  colnames(Y) <- c("wing", "tail", "beak", "leg")
  X <- cbind(intercept = 1, size = stats::rnorm(nrow(Y)))
  subsets <- rbind(c(2L, 4L, 1L), c(1L, 3L, 4L), c(2L, 3L, 1L))
  fit <- fit_G(Y, X = X, A = data$A, subsets = subsets, engine = "sommer",
               cores = 1L, n_draws = 30L, seed = 7L)

  D <- outer(fit$y_sd, fit$y_sd)
  Y_std <- sweep(sweep(Y, 2L, fit$y_mean, "-"), 2L, fit$y_sd, "/")
  G_std <- unname(fit$G / D)
  expect_equal(unname(fit$E), unname(residual_E(Y_std, X, G_std, data$A) * D),
               tolerance = 1e-8)
  # E_bar is the plain average of the block E
  E_blocks <- lapply(fit$blocks, function(b) list(G = unname(b$E)))
  avg <- indri_assemble(E_blocks, fit$subsets, 4L)$G
  expect_equal(unname(fit$E_bar), unname(avg * D), tolerance = 1e-10)
  expect_equal(dimnames(fit$E), list(colnames(Y), colnames(Y)))
  expect_equal(dimnames(fit$E_bar), list(colnames(Y), colnames(Y)))
  expect_equal(fit$E, t(fit$E), tolerance = 1e-10)
  expect_gt(min_eigen(fit$E), 0)
  # the blocks keep their standardised E
  expect_lt(max(diag(fit$blocks[[1]]$E)), 2)
})

test_that("E agrees with the full model's residual covariance", {
  skip_if_not_installed("sommer")
  # On this balanced sibship design REML reduces to the ANOVA estimators,
  # which split the residual covariance into c G + E exactly: measured gap
  # under 5e-4 over seeds 1-6.
  data <- simulate_block_data(n_fam = 60, size = 5, K = 4L, seed = 3)
  subsets <- rbind(c(2L, 4L, 1L), c(1L, 3L, 4L), c(2L, 3L, 1L))
  fit <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "sommer",
               cores = 1L, n_draws = 30L, seed = 7L)
  full <- fit_G(data$Y, A = data$A, subsets = matrix(1:4, 1L),
                assembly = "average", engine = "sommer", cores = 1L,
                n_draws = 30L, seed = 7L)
  # a single block covering every trait: its E is the full model's
  expect_lt(max(abs(fit$E - full$E_bar)) / mean(diag(full$E_bar)), 2e-3)
})

test_that("E is floored to stay well conditioned", {
  set.seed(5)
  N <- 60
  Y <- matrix(stats::rnorm(N * 3), N, 3)
  X <- matrix(1, N, 1)
  A <- diag(N)
  # a G larger than the residual covariance leaves nothing for E
  E <- .indri_residual_E(Y, X, diag(3) * 5, A)
  ev <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), 0)
  expect_gte(min(ev) / max(ev), 1e-3 - 1e-12)
  # with A = I and a small G the residual covariance splits exactly
  G <- diag(3) * 0.1
  P <- crossprod(scale(Y, scale = FALSE)) / (N - 1)
  expect_equal(.indri_residual_E(Y, X, G, A), unname(P - G), tolerance = 1e-10)
})

test_that("E is computed from the observed values when Y has gaps", {
  set.seed(6)
  N <- 80
  Y <- matrix(stats::rnorm(N * 3), N, 3) %*% chol(matrix(c(1, .5, .2,
                                                           .5, 1, .3,
                                                           .2, .3, 1), 3))
  Y[c(3, 9, 20), 1] <- NA
  Y[c(5, 40), 3] <- NA
  X <- cbind(1, stats::rnorm(N))
  E <- .indri_residual_E(Y, X, diag(3) * 0.1, diag(N))
  expect_false(anyNA(E))
  expect_gt(min(eigen(E, symmetric = TRUE, only.values = TRUE)$values), 0)
  # the (2, 3) entry uses the rows both traits have
  rows <- !is.na(Y[, 3])
  r2 <- stats::lm.fit(X, Y[, 2])$residuals[rows]
  r3 <- stats::lm.fit(X[rows, ], Y[rows, 3])$residuals
  expect_equal(E[2, 3], sum(r2 * r3) / (sum(rows) - 2), tolerance = 1e-10)
})
