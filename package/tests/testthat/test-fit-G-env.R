# fit_G_env(): genetic covariance within and between the groups of a
# discrete environment, from character states fitted with fit_G().

sim_env <- function(seed, n_fam = 100, size = 6, P = 3, r_between = 0.5) {
  set.seed(seed)
  N <- n_fam * size
  A <- kronecker(diag(n_fam), matrix(0.5, size, size))
  diag(A) <- 1
  ids <- paste0("id", seq_len(N))
  dimnames(A) <- list(ids, ids)
  sex <- stats::rbinom(N, 1L, 0.5)
  G0 <- matrix(0.3, P, P)
  diag(G0) <- 1
  Gs <- kronecker(matrix(c(1, r_between, r_between, 1), 2), G0)
  U <- t(chol(A)) %*% matrix(stats::rnorm(N * 2 * P), N, 2 * P) %*% chol(Gs)
  E0 <- matrix(0.2, P, P)
  diag(E0) <- 1
  noise <- matrix(stats::rnorm(N * P), N, P) %*% chol(E0)
  own <- ifelse(sex == 1, P, 0)
  Y <- t(vapply(seq_len(N), function(i) U[i, own[i] + seq_len(P)],
                numeric(P))) + noise
  dimnames(Y) <- list(ids, paste0("t", seq_len(P)))
  list(Y = Y, A = A, sex = ifelse(sex == 1, "m", "f"), N = N, P = P)
}

fit_env <- function(d, ...) {
  fit_G_env(d$Y, env = d$sex, A = d$A, cores = 1L, seed = 1L, n_draws = 30L,
            progress = FALSE, ...)
}

test_that("character states give G within and between environments", {
  skip_if_not_installed("sommer")
  d <- sim_env(2, n_fam = 150)
  fit <- fit_env(d)
  expect_s3_class(fit, "indri_fit_env")
  states <- c("t1:f", "t2:f", "t3:f", "t1:m", "t2:m", "t3:m")
  expect_equal(fit$states, states)
  expect_equal(dimnames(fit$G), list(states, states))
  expect_equal(fit$G, fit$fit$G)
  expect_equal(dimnames(fit$corr_between), list(c("t1", "t2", "t3"), "f-m"))
  expect_equal(fit$corr_between["t2", "f-m"],
               fit$G["t2:f", "t2:m"] / sqrt(fit$G["t2:f", "t2:f"] *
                                              fit$G["t2:m", "t2:m"]))
  # true cross-sex genetic correlation 0.5. At this seed the mean estimate
  # is 0.495 (trait-paired subsets); over seeds 1-8 the means ran 0.31-0.54
  # (individual traits 0.18-0.79): the cross-environment part rests on
  # relatives of the other sex only.
  expect_lt(abs(mean(fit$corr_between) - 0.5), 0.1)
  expect_output(print(fit), "between environments")
})

test_that("E is block diagonal: residual covariance minus G within each group", {
  skip_if_not_installed("sommer")
  d <- sim_env(3)
  fit <- fit_env(d)
  f <- c("t1:f", "t2:f", "t3:f")
  m <- c("t1:m", "t2:m", "t3:m")
  expect_true(all(fit$E[f, m] == 0))
  expect_true(all(fit$E[m, f] == 0))
  rows <- d$sex == "m"
  X <- matrix(1, d$N, 1L)
  ref <- .indri_residual_E(d$Y[rows, ], X[rows, , drop = FALSE],
                           unname(fit$G[m, m]), d$A[rows, rows])
  expect_equal(unname(fit$E[m, m]), unname(ref))
  expect_gt(min(eigen(fit$E[f, f], symmetric = TRUE,
                      only.values = TRUE)$values), 0)
})

test_that("env can name a column of X, whose main effect is then dropped", {
  skip_if_not_installed("sommer")
  d <- sim_env(4, n_fam = 40)
  X <- cbind(intercept = 1, male = as.numeric(d$sex == "m"),
             age = stats::rnorm(d$N))
  rownames(X) <- rownames(d$Y)
  by_name <- fit_G_env(d$Y, env = "male", X = X, A = d$A, cores = 1L,
                       seed = 1L, n_draws = 30L, progress = FALSE)
  by_vector <- fit_G_env(d$Y, env = X[, "male"], X = X[, c(1, 3)], A = d$A,
                         cores = 1L, seed = 1L, n_draws = 30L,
                         progress = FALSE)
  expect_equal(by_name$G, by_vector$G)
  expect_equal(by_name$levels, c("0", "1"))
})

test_that("inputs are checked", {
  d <- sim_env(5, n_fam = 20)
  expect_error(fit_G_env(d$Y, env = d$sex[-1], A = d$A), "one value per row")
  na <- d$sex
  na[3] <- NA
  expect_error(fit_G_env(d$Y, env = na, A = d$A), "missing")
  expect_error(fit_G_env(d$Y, env = rep("f", d$N), A = d$A), "two levels")
  expect_error(fit_G_env(d$Y, env = "diet", A = d$A), "column")
  # a covariate constant within an environment duplicates its main effect
  X <- cbind(intercept = 1, male = as.numeric(d$sex == "m"))
  rownames(X) <- rownames(d$Y)
  expect_error(fit_G_env(d$Y, env = d$sex, X = X, A = d$A), "rank deficient")
})

test_that("the default subsets pair traits across every environment", {
  sub <- .indri_env_subsets(4L, 2L)
  expect_equal(dim(sub), c(6L, 4L))
  expect_equal(sub[1, ], c(1L, 2L, 5L, 6L))
  expect_equal(sub[6, ], c(3L, 4L, 7L, 8L))
  expect_equal(ncol(.indri_env_subsets(3L, 3L)), 6L)
  expect_null(.indri_env_subsets(2L, 2L))
  # every state pair is covered, and the cross-environment entries of a
  # trait in every block that holds it
  cover <- indri_pair_coverage(sub, 8L)
  expect_true(all(cover[upper.tri(cover)] > 0))
  expect_equal(diag(cover[1:4, 5:8]), rep(3L, 4))
})

test_that("fit_G_env passes trait-paired subsets unless subsets are chosen", {
  skip_if_not_installed("sommer")
  d <- sim_env(5, n_fam = 40)
  fit <- fit_env(d)
  expect_equal(fit$fit$subsets, .indri_env_subsets(3L, 2L))
  own <- suppressWarnings(fit_env(d, s = 2L))  # 2-trait blocks may stall
  expect_false(identical(dim(own$fit$subsets), dim(fit$fit$subsets)))
})
