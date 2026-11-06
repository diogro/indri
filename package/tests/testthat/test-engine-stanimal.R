test_that("the stanimal engine returns a covariance and its draws", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 6, size = 5)
  engine <- indri_engine("stanimal")
  forms <- indri_a_forms(data$A, need = engine$capabilities$a_form)

  fit <- engine$fit(data$Y, matrix(1, data$N, 1L), forms,
                    list(chains = 1L, iter = 200L, cores = 1L,
                         seed = 1L, refresh = 0, show_messages = FALSE))

  expect_equal(dim(fit$G), c(data$K, data$K))
  expect_equal(fit$G, t(fit$G), tolerance = 1e-10)
  expect_equal(dim(fit$draws)[1:2], c(data$K, data$K))
  expect_gt(dim(fit$draws)[3], 50L)
  expect_true(all(is.finite(fit$G)))
})

test_that("engines that are designed but not implemented say so", {
  expect_error(indri_engine("sommer"), "not implemented")
  expect_error(indri_engine("asreml"), "not implemented")
  expect_error(indri_engine("nonsense"), "Unknown engine")
  expect_equal(indri_engines(), "stanimal")
})
