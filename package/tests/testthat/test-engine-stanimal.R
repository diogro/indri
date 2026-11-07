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

test_that("unknown engines are refused and the registry lists the implemented ones", {
  expect_error(indri_engine("nonsense"), "Unknown engine")
  expect_equal(indri_engines(), c("stanimal", "sommer", "asreml"))
})

test_that("stanimal blocks write to their own output directory and clean up", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 6, size = 5, K = 3L)
  before <- list.files(tempdir(), pattern = "^indri-block-")
  # Same random state for both fits, as in tests that reseed between fits:
  # cmdstanr would give both the same output path in one minute.
  block <- function() {
    set.seed(42)
    fit_G(data$Y, A = data$A, subsets = t(utils::combn(3L, 2L)), engine = "stanimal",
          engine_args = list(chains = 1L, iter = 100L, seed = 1L, init = 0.5,
                             refresh = 0, show_messages = FALSE))
  }
  first <- block()
  second <- block()
  invisible(gc())
  expect_identical(first$G, second$G)
  expect_identical(list.files(tempdir(), pattern = "^indri-block-"), before)
})
