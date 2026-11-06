# M1 acceptance checks (a) and (c): the blocks must be on one common scale, and
# two blocks that share traits must agree about them.
#
# The whole estimator rests on the first of these. indri standardises every trait
# once, globally, before the subsets are drawn, so that a trait enters every
# block it belongs to as the same numbers. If that were not true, the average
# over overlapping blocks would be mixing estimates on different scales, and
# nothing downstream would notice.

# The standardisation indri applies, recomputed here so that the comparison does
# not have to reach into the fit.
standardise_like_indri <- function(Y) {
  y_mean <- colMeans(Y, na.rm = TRUE)
  y_sd <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  sweep(sweep(Y, 2L, y_mean, "-"), 2L, y_sd, "/")
}

test_that("a block is exactly the model fitted on the standardised columns", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 10, size = 5, K = 3L)
  subsets <- matrix(c(1L, 2L), nrow = 1L)
  engine_args <- list(chains = 1L, iter = 200L, cores = 1L, seed = 4L,
                      init = 0.5, refresh = 0, show_messages = FALSE)

  fit <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "stanimal",
               engine_args = engine_args)

  Y_std <- standardise_like_indri(data$Y)
  direct <- stanimal::lmm_animal(Y_std[, 1:2, drop = FALSE],
                                 matrix(1, data$N, 1L), data$A,
                                 chains = 1L, iter = 200L, cores = 1L,
                                 seed = 4L, init = 0.5, refresh = 0,
                                 show_messages = FALSE)

  # Same data, same seed, same start: the block must be the same model, which
  # means identical draws, not merely close ones.
  expect_identical(posterior::as_draws_matrix(fit$blocks[[1L]]$draws),
                   posterior::as_draws_matrix(stanimal::fit_draws(direct)$G))
})

test_that("the fit records the scale that undoes the standardisation", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 10, size = 5, K = 3L)
  subsets <- matrix(c(1L, 2L), nrow = 1L)

  fit <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "stanimal",
               engine_args = list(chains = 1L, iter = 300L, cores = 1L,
                                  seed = 5L, init = 0.5, refresh = 0,
                                  show_messages = FALSE))

  # the recorded scale is the data's, trait by trait
  expect_equal(unname(fit$y_sd), unname(apply(data$Y, 2L, stats::sd)))
  expect_equal(unname(fit$y_mean), unname(colMeans(data$Y)))

  # G is reported on the standardised scale, and undoing the scaling cannot
  # change a correlation
  block <- fit$blocks[[1L]]$G
  back <- diag(fit$y_sd[1:2]) %*% block %*% diag(fit$y_sd[1:2])
  expect_equal(unname(stats::cov2cor(back)), unname(stats::cov2cor(block)),
               tolerance = 1e-10)
  # and the back-transformed variances are the standardised ones scaled by the
  # squared trait standard deviations
  expect_equal(unname(diag(back)),
               unname(fit$y_sd[1:2]^2 * diag(block)), tolerance = 1e-10)
})

test_that("two blocks that share traits agree about them", {
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 15, size = 5, K = 4L)
  # traits 1,2,3 and 2,3,4: the pair 2,3 is estimated twice, by different models
  subsets <- rbind(c(1L, 2L, 3L), c(2L, 3L, 4L))

  fit <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "stanimal",
               engine_args = list(chains = 1L, iter = 500L, cores = 1L,
                                  seed = 6L, init = 0.5, refresh = 0,
                                  show_messages = FALSE))

  first <- fit$blocks[[1L]]
  second <- fit$blocks[[2L]]

  # position of traits 2 and 3 in each block
  shared_first <- first$G[2:3, 2:3]
  shared_second <- second$G[1:2, 1:2]
  sd_first <- apply(first$draws[2:3, 2:3, , drop = FALSE], c(1L, 2L), stats::sd)
  sd_second <- apply(second$draws[1:2, 1:2, , drop = FALSE], c(1L, 2L),
                     stats::sd)

  # One posterior standard deviation each, so the difference between two
  # estimates of the same quantity should sit well inside a few of them. Loose
  # enough not to be flaky, tight enough to catch a scale or index mismatch.
  tolerance <- 4 * sqrt(sd_first^2 + sd_second^2)
  expect_true(all(abs(shared_first - shared_second) <= tolerance))
  # and they should not be trivially equal (e.g. both zero)
  expect_gt(max(abs(shared_first)), 0.05)
})

test_that("shared traits are identical numbers in every block", {
  data <- simulate_block_data(n_fam = 4, size = 5, K = 4L)
  Y_std <- standardise_like_indri(data$Y)

  # what each block receives for trait 2 and 3
  expect_identical(Y_std[, 2:3, drop = FALSE],
                   Y_std[, c(2L, 3L), drop = FALSE])
  expect_equal(unname(apply(Y_std, 2L, stats::sd)), rep(1, 4L),
               tolerance = 1e-12)
})
