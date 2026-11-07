test_that("fit_G assembles a usable covariance matrix end to end", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 8, size = 5, K = 4L)

  subsets <- t(utils::combn(4L, 2L))
  fit <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "stanimal",
               engine_args = list(chains = 1L, iter = 200L, cores = 1L,
                                  seed = 1L, refresh = 0, show_messages = FALSE))

  expect_s3_class(fit, "indri_fit")
  expect_equal(dim(fit$G), c(4L, 4L))
  expect_equal(fit$G, t(fit$G), tolerance = 1e-8)
  expect_equal(fit$corrG, t(fit$corrG), tolerance = 1e-8)
  expect_equal(unname(diag(fit$corrG)), rep(1, 4L), tolerance = 1e-8)
  expect_gt(min_eigen(fit$G), -1e-8)
  expect_equal(dim(fit$G_draws), c(4L, 4L, dim(fit$G_draws)[3]))
  expect_equal(length(fit$y_sd), 4L)
  expect_true(all(fit$counts[upper.tri(fit$counts)] >= 1))
  expect_true(all(fit$coverage[upper.tri(fit$coverage)] >= 1))
  expect_equal(unname(rownames(fit$G)), unname(colnames(data$Y)))
})

test_that("the reported G is on the standardised scale and back-transforms", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 8, size = 5, K = 3L)
  fit <- fit_G(data$Y, A = data$A,
               subsets = t(utils::combn(3L, 2L)), engine = "stanimal",
               engine_args = list(chains = 1L, iter = 200L, cores = 1L,
                                  seed = 2L, refresh = 0, show_messages = FALSE))

  original_scale <- diag(fit$y_sd) %*% fit$G %*% diag(fit$y_sd)
  expect_equal(unname(diag(original_scale)), unname(fit$y_sd^2 * diag(fit$G)),
               tolerance = 1e-10)
  expect_gt(min_eigen(original_scale), -1e-8)
})

test_that("fit_G refuses what it cannot do", {
  data <- simulate_block_data(n_fam = 4, size = 5, K = 3L)

  expect_error(fit_G(data$Y[, 1, drop = FALSE], A = data$A,
                     subsets = matrix(c(1L, 1L), 1)),
               "at least two traits")
  expect_error(fit_G(data$Y, A = data$A[1:5, 1:5],
                     subsets = t(utils::combn(3L, 2L))),
               "individuals")
  expect_error(fit_G(data$Y, A = data$A,
                     subsets = matrix(rep(1:3, 2), 2, byrow = TRUE)),
               "nothing to reassemble")
})

test_that("diagnostics reads logical and Rhat convergence", {
  # REML engines report a logical flag; one FALSE must sink the summary.
  d <- .indri_diagnostics(list(list(converged = TRUE),
                               list(converged = TRUE),
                               list(converged = FALSE)))
  expect_false(d$converged)
  expect_true(all(is.na(d$rhat)))
  expect_false(is.finite(d$max_rhat))

  expect_true(.indri_diagnostics(list(list(converged = TRUE),
                                       list(converged = TRUE)))$converged)

  # Sampler engines report max Rhat per block; the 1.05 threshold applies.
  expect_true(.indri_diagnostics(list(list(converged = 1.01),
                                       list(converged = 1.02)))$converged)
  expect_false(.indri_diagnostics(list(list(converged = 1.06),
                                        list(converged = 1.02)))$converged)
  expect_equal(.indri_diagnostics(list(list(converged = 1.02),
                                        list(converged = 1.01)))$max_rhat, 1.02)

  # Unknown convergence (e.g. an rstan backend with no Rhat) stays optimistic,
  # matching the original behaviour.
  expect_true(.indri_diagnostics(list(list(converged = NA_real_)))$converged)
})

test_that("print and summary work on a fit", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 6, size = 5, K = 3L)
  fit <- fit_G(data$Y, A = data$A,
               subsets = t(utils::combn(3L, 2L)), engine = "stanimal",
               engine_args = list(chains = 1L, iter = 200L, cores = 1L,
                                  seed = 3L, refresh = 0, show_messages = FALSE))

  expect_output(print(fit), "indri fit")
  s <- summary(fit)
  expect_s3_class(s, "summary.indri_fit")
  expect_output(print(s), "indri fit summary")
})
