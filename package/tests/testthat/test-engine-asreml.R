test_that("the asreml engine returns a covariance and standard errors", {
  skip_if_no_asreml()
  data <- simulate_block_data(n_fam = 40, size = 5, K = 2L)
  engine <- indri_engine("asreml")
  expect_equal(engine$name, "asreml")
  expect_equal(engine$capabilities$a_form, "sparse_inverse")
  expect_false(engine$capabilities$draws)
  expect_true(engine$capabilities$se)
  expect_false(engine$capabilities$missing)

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
