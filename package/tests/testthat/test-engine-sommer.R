test_that("the sommer engine returns a covariance and standard errors", {
  skip_if_no_sommer()
  data <- simulate_block_data(n_fam = 20, size = 4, K = 3L)
  engine <- indri_engine("sommer")
  expect_equal(engine$name, "sommer")
  expect_equal(engine$capabilities$a_form, "inverse")
  expect_false(engine$capabilities$draws)
  expect_true(engine$capabilities$se)
  expect_false(engine$capabilities$missing)

  forms <- indri_a_forms(data$A, need = engine$capabilities$a_form)
  expect_true("inverse" %in% names(forms))

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
  expect_true(all(is.finite(fit$se)))
  expect_true(all(diag(fit$se) >= 0))

  expect_null(fit$draws)
  expect_true(isTRUE(fit$converged))
})

test_that("fit_G runs end to end with the sommer engine", {
  skip_if_no_sommer()
  data <- simulate_block_data(n_fam = 8, size = 5, K = 3L)

  fit <- fit_G(data$Y, A = data$A, subsets = t(utils::combn(3L, 2L)),
               engine = "sommer")

  expect_s3_class(fit, "indri_fit")
  expect_equal(fit$engine, "sommer")
  expect_equal(dim(fit$G), c(3L, 3L))
  expect_equal(fit$G, t(fit$G), tolerance = 1e-8)
  expect_gt(min_eigen(fit$G), -1e-8)
  expect_null(fit$G_draws)
  expect_equal(unname(rownames(fit$G)), unname(colnames(data$Y)))
})

test_that("starting covariances are validated and reach the REML fit", {
  expect_equal(.indri_check_init(NULL, 2L), list(G = NULL, E = NULL))
  expect_equal(.indri_check_init(list(G = diag(2)), 2L),
               list(G = diag(2), E = NULL))
  expect_error(.indri_check_init(list(G = diag(3)), 2L), "positive-definite 2 x 2")
  expect_error(.indri_check_init(list(E = matrix(c(1, 2, 2, 1), 2)), 2L),
               "positive-definite")
  expect_error(.indri_check_init(diag(2), 2L), "must be a list")

  skip_if_not_installed("sommer")
  data <- simulate_block_data(n_fam = 20, size = 5, K = 2L)
  Y <- scale(data$Y)
  forms <- indri_a_forms(data$A, need = "inverse")
  X <- matrix(1, nrow(Y), 1, dimnames = list(rownames(Y), "intercept"))
  plain <- .indri_sommer_fit(Y, X, forms)
  started <- .indri_sommer_fit(Y, X, forms, list(init = list(
    G = plain$G, E = plain$E)))
  # starting at the optimum, the fit stays there
  expect_equal(started$G, plain$G, tolerance = 1e-3)
})

test_that("the sommer engine fits a single trait", {
  skip_if_not_installed("sommer")
  data <- simulate_block_data(n_fam = 20, size = 5, K = 2L)
  Y <- scale(data$Y)
  forms <- indri_a_forms(data$A, need = "inverse")
  X <- matrix(1, nrow(Y), 1, dimnames = list(rownames(Y), "intercept"))
  one <- .indri_sommer_fit(Y[, 1, drop = FALSE], X, forms)
  expect_equal(dim(one$G), c(1L, 1L))
  expect_equal(dim(one$E), c(1L, 1L))
  expect_gt(one$G[1, 1], 0)
  # on standardised data G + E is about the phenotypic variance, 1
  expect_equal(one$G[1, 1] + one$E[1, 1], 1, tolerance = 0.25)
})
