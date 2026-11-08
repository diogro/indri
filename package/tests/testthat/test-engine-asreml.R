test_that("the asreml engine returns a covariance and standard errors", {
  skip_if_no_asreml()
  data <- simulate_block_data(n_fam = 40, size = 5, K = 2L)
  engine <- indri_engine("asreml")
  expect_equal(engine$name, "asreml")
  expect_equal(engine$capabilities$a_form, "sparse_inverse")
  expect_false(engine$capabilities$draws)
  expect_true(engine$capabilities$se)
  expect_true(engine$capabilities$missing)   # gaps are single missing records

  forms <- indri_a_forms(data$A, need = engine$capabilities$a_form)
  expect_true("sparse_inverse" %in% names(forms))

  fit <- engine$fit(data$Y, matrix(1, data$N, 1L), forms, list())

  expect_equal(dim(fit$G), c(data$K, data$K))
  expect_equal(fit$G, t(fit$G), tolerance = 1e-10)
  expect_true(all(is.finite(fit$G)))
  expect_true(all(diag(fit$G) >= -1e-8))

  expect_equal(dim(fit$E), c(data$K, data$K))
  expect_equal(fit$E, t(fit$E), tolerance = 1e-10)
  expect_true(all(is.finite(fit$E)))

  expect_equal(dim(fit$se), c(data$K, data$K))
  expect_equal(fit$se, t(fit$se), tolerance = 1e-10)
  expect_true(all(diag(fit$se) >= 0))

  expect_null(fit$draws)
  expect_true(isTRUE(fit$converged))
})

test_that("fit_G runs end to end with the asreml engine", {
  skip_if_no_asreml()
  data <- simulate_block_data(n_fam = 40, size = 5, K = 3L)

  fit <- fit_G(data$Y, A = data$A, subsets = t(utils::combn(3L, 2L)),
               engine = "asreml")

  expect_s3_class(fit, "indri_fit")
  expect_equal(fit$engine, "asreml")
  expect_equal(dim(fit$G), c(3L, 3L))
  expect_equal(fit$G, t(fit$G), tolerance = 1e-8)
  expect_gt(min_eigen(fit$G), -1e-8)
  # REML blocks now yield draws from their sampling covariance
  expect_equal(dim(fit$G_draws)[1:2], c(3L, 3L))
  expect_true(all(abs(fit$corrG_draws) <= 1))
  expect_equal(unname(rownames(fit$G)), unname(colnames(data$Y)))
})

test_that("asreml start values come from init when given", {
  expect_equal(.indri_asreml_us_init(2L), c(1, 0, 1))
  m <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_equal(.indri_asreml_us_init(2L, m), c(2, 0.5, 1))
})

test_that("the asreml engine fits a single trait", {
  skip_if_no_asreml()
  data <- simulate_block_data(n_fam = 20, size = 5, K = 2L)
  Y <- scale(data$Y)
  forms <- indri_a_forms(data$A, need = "sparse_inverse")
  X <- matrix(1, nrow(Y), 1, dimnames = list(rownames(Y), "intercept"))
  one <- .indri_asreml_fit(Y[, 1, drop = FALSE], X, forms)
  expect_equal(dim(one$G), c(1L, 1L))
  expect_gt(one$G[1, 1], 0)
  expect_equal(one$G[1, 1] + one$E[1, 1], 1, tolerance = 0.25)
})

test_that("asreml blocks carry the inverse AI matrix for their G entries", {
  skip_if_no_asreml()
  data <- simulate_block_data(n_fam = 30, size = 5, K = 2L)
  Y <- scale(data$Y)
  X <- matrix(1, nrow(Y), 1, dimnames = list(rownames(Y), "intercept"))
  out <- .indri_asreml_fit(Y, X, indri_a_forms(data$A, need = "sparse_inverse"))
  expect_equal(dim(out$vcov), c(3L, 3L))
  expect_equal(sqrt(diag(out$vcov)), out$se[.indri_vech_index(2L)],
               tolerance = 1e-6)
})

test_that("asreml converges with covariates in X and agrees with sommer", {
  skip_if_no_asreml()
  skip_if_not_installed("sommer")
  # Regression: asreml's own budget (13 iterations, one update) stopped far
  # from the optimum once X had covariates, G off by 0.12 where sommer was
  # right.
  set.seed(3)
  A <- kronecker(diag(120), matrix(0.5, 5, 5)); diag(A) <- 1
  N <- nrow(A); ids <- paste0("id", seq_len(N))
  G <- matrix(c(1, 0.5, 0.3, 0.5, 0.8, -0.2, 0.3, -0.2, 0.9), 3)
  x <- rnorm(N, 50, 10); grp <- sample(1:3, N, TRUE)
  Y <- outer(x, c(0.3, -0.2, 0.15)) +
    rbind(c(0, 0, 0), c(3, -2, 1), c(-2, 4, 2))[grp, ] +
    t(chol(A)) %*% matrix(rnorm(N * 3), N) %*% chol(G) +
    matrix(rnorm(N * 3), N) %*% diag(sqrt(c(0.6, 0.7, 0.5)))
  Y <- scale(Y); dimnames(Y) <- list(ids, c("t1", "t2", "t3"))
  X <- cbind(intercept = 1, x = x, gb = grp == 2, gc = grp == 3)
  rownames(X) <- ids
  fa <- .indri_asreml_fit(Y, X, indri_a_forms(A, need = "sparse_inverse"))
  fs <- .indri_sommer_fit(Y, X, indri_a_forms(A, need = "inverse"))
  expect_true(fa$converged)
  expect_equal(unname(fa$G), unname(fs$G), tolerance = 1e-3)
})
