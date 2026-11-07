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
  # correlation draws: mixed per block, so every entry is a correlation
  expect_equal(fit$spread, "mixture")
  expect_equal(dim(fit$corrG_draws), dim(fit$G_draws))
  expect_true(all(abs(fit$corrG_draws) <= 1 + 1e-12))
  expect_true(all(fit$corrG_draws[cbind(1:4, 1:4, 1L)] == 1))
  expect_equal(dimnames(fit$corrG_draws)[1:2], dimnames(fit$G)) 
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

test_that("the default subset size grows with K but caps at 4", {
  expect_equal(.indri_default_s(3L), 2L)
  expect_equal(.indri_default_s(4L), 2L)
  expect_equal(.indri_default_s(12L), 3L)
  expect_equal(.indri_default_s(16L), 4L)
  expect_equal(.indri_default_s(30L), 4L)
  expect_equal(.indri_default_s(50L), 4L)
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

test_that("subset size, co-occurrence and workers default by engine", {
  stanimal <- .indri_engine_defaults(indri_engine("stanimal"), 16L)
  expect_equal(stanimal, list(s = 4L, min_cooccurrence = 3L))
  expect_equal(.indri_engine_defaults(indri_engine("stanimal"), 9L)$s, 3L)

  # REML: pairs, once each, whatever K is
  for (name in c("sommer", "asreml")) {
    for (K in c(3L, 16L, 50L)) {
      expect_equal(.indri_engine_defaults(indri_engine(name), K),
                   list(s = 2L, min_cooccurrence = 1L), label = paste(name, K))
    }
  }
  # one asreml licence seat per worker; a standard licence has two
  expect_equal(indri_engine("asreml")$capabilities$default_cores, 2L)
  # ...unless this session already holds one
  expected <- if ("asreml" %in% loadedNamespaces()) 1L else 2L
  expect_equal(.indri_default_cores(indri_engine("asreml")), expected)
  expect_equal(.indri_default_cores(indri_engine("sommer")), 1L)
  expect_null(indri_engine("sommer")$capabilities$default_cores)
  expect_null(indri_engine("stanimal")$capabilities$default_cores)

  # pairs once each: every pair covered, no repeats beyond the window wrap
  sub <- indri_subsets(16L, 2L, min_cooccurrence = 1L, seed = 1L)
  cov <- attr(sub, "coverage")
  expect_gte(min(cov[upper.tri(cov)]), 1L)
  expect_lte(nrow(sub), choose(16, 2) + 1L)
})
