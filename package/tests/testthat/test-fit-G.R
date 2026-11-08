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

test_that("the standardised-to-original conversion keeps NA, NULL, draws and names", {
  y_sd <- c(a = 2, b = 10, c = 0.5)
  G_std <- matrix(c(1, 0.5, NA, 0.5, 1, 0.2, NA, 0.2, 1), 3,
                  dimnames = list(names(y_sd), names(y_sd)))
  G <- .indri_to_original(G_std, y_sd)
  expect_equal(G, G_std * outer(y_sd, y_sd))
  expect_true(is.na(G[1, 3]) && is.na(G[3, 1]))
  expect_equal(dimnames(G), dimnames(G_std))
  expect_null(.indri_to_original(NULL, y_sd))
  draws <- array(rep(c(G_std), 4) * rep(1:4, each = 9), c(3, 3, 4),
                 dimnames = list(names(y_sd), names(y_sd), NULL))
  back <- .indri_to_original(draws, y_sd)
  expect_equal(dim(back), dim(draws))
  expect_equal(dimnames(back), dimnames(draws))
  for (d in 1:4) expect_equal(back[, , d], draws[, , d] * outer(y_sd, y_sd))
})

test_that("assembled covariances come back in the units of Y", {
  skip_if_not_installed("sommer")
  data <- simulate_block_data(n_fam = 30, size = 5, K = 4L)
  # very unequal scales, named traits that do not sort in column order
  Y <- sweep(data$Y, 2L, c(1, 50, 0.02, 3), "*") + 100
  colnames(Y) <- c("wing", "tail", "beak", "leg")
  X <- cbind(intercept = 1, size = rnorm(nrow(Y)))
  subsets <- rbind(c(2L, 4L, 1L), c(1L, 3L, 4L), c(2L, 3L, 1L))

  recompute <- function(fit, assembly, weights, spread) {
    blocks <- fit$blocks
    K <- ncol(Y)
    avg <- indri_assemble(blocks, fit$subsets, K, weights = weights)$G
    G_std <- if (assembly == "ls") {
      indri_project_psd(blocks, fit$subsets, K, init = avg,
                        weights = weights)$matrix
    } else {
      m <- min(eigen((avg + t(avg)) / 2, only.values = TRUE)$values)
      if (m < 0) avg + (-m + 1e-8) * diag(K) else avg
    }
    draws <- indri_assemble_draws(lapply(blocks, `[[`, "draws"), fit$subsets,
                                  K, method = spread, seed = 7L)
    list(G = G_std, G_bar = avg, draws = draws)
  }
  for (assembly in c("ls", "average")) for (weights in c("none", "inverse_variance"))
    for (spread in c("mixture", "average")) {
      info <- paste(assembly, weights, spread)
      fit <- fit_G(Y, X = X, A = data$A, subsets = subsets, engine = "sommer",
                   assembly = assembly, weights = weights, spread = spread,
                   cores = 1L, n_draws = 30L, seed = 7L)
      std <- recompute(fit, assembly, weights, spread)
      D <- outer(fit$y_sd, fit$y_sd)
      expect_equal(unname(fit$G), unname(std$G * D), tolerance = 1e-10,
                   label = info)
      expect_equal(unname(fit$G_bar), unname(std$G_bar * D), tolerance = 1e-10,
                   label = info)
      for (d in c(1L, 30L))
        expect_equal(unname(fit$G_draws[, , d]), unname(std$draws[, , d] * D),
                     tolerance = 1e-10, label = info)
      expect_equal(unname(fit$corrG), unname(stats::cov2cor(std$G)),
                   tolerance = 1e-12, label = info)
      expect_equal(unname(fit$sigma_G), unname(sqrt(diag(fit$G))), label = info)
      expect_equal(dimnames(fit$G), list(colnames(Y), colnames(Y)))
      expect_equal(fit$G, t(fit$G), tolerance = 1e-10)
      if (assembly == "ls") expect_gt(min_eigen(fit$G), -1e-6 * max(diag(fit$G)))
      # the raw blocks stay on the standardised scale
      expect_lt(max(diag(fit$blocks[[1]]$G) + diag(fit$blocks[[1]]$E)), 2)
    }
  # the scale is the data's, so the tail variance is ~2500 times the wing's
  expect_gt(fit$G["tail", "tail"] / fit$G["wing", "wing"], 100)

  # the complete-rows missing-data path converts with the scales of the rows kept
  Y_gap <- Y
  Y_gap[c(3, 17), 2] <- NA
  fit <- fit_G(Y_gap, X = X, A = data$A, subsets = subsets, engine = "sommer",
               missing = "complete_block", cores = 1L, n_draws = 30L,
               seed = 7L)
  kept <- stats::complete.cases(Y_gap)
  expect_equal(unname(fit$y_sd), unname(apply(Y_gap[kept, ], 2L, sd)))
  std <- recompute(fit, "ls", "none", "mixture")
  expect_equal(unname(fit$G), unname(std$G * outer(fit$y_sd, fit$y_sd)),
               tolerance = 1e-10)
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
  expect_equal(.indri_default_s(16L), 3L)
  expect_equal(.indri_default_s(30L), 3L)
  expect_equal(.indri_default_s(50L), 3L)
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
  # no block could tell: unknown, not converged
  expect_true(is.na(.indri_diagnostics(list(list(converged = NA_real_)))$converged))
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
  expect_equal(stanimal, list(s = 3L, min_cooccurrence = 3L))
  expect_equal(.indri_engine_defaults(indri_engine("stanimal"), 9L)$s, 3L)

  # REML: 4-trait blocks, each pair in at least one (s capped at K)
  for (name in c("sommer", "asreml")) {
    for (K in c(16L, 50L)) {
      expect_equal(.indri_engine_defaults(indri_engine(name), K),
                   list(s = 3L, min_cooccurrence = 1L), label = paste(name, K))
    }
    # blocks must stay narrower than K: pairs for 3 or 4 traits
    expect_equal(.indri_engine_defaults(indri_engine(name), 3L)$s, 2L)
    expect_equal(.indri_engine_defaults(indri_engine(name), 4L)$s, 2L)
    expect_equal(.indri_engine_defaults(indri_engine(name), 5L)$s, 3L)
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

test_that("diagnostics list the blocks that did not converge, by trait", {
  blocks <- list(list(converged = TRUE), list(converged = FALSE),
                 list(converged = 1.2), list(converged = NA))
  subsets <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L), c(1L, 2L))
  d <- .indri_diagnostics(blocks, subsets, c("a", "b", "c"))
  expect_false(d$converged)
  expect_equal(d$n_unconverged, 2L)
  expect_equal(d$n_unknown, 1L)
  expect_equal(d$unconverged$block, c(2L, 3L))
  expect_equal(d$unconverged$traits, c("a, c", "b, c"))
  expect_equal(d$unconverged$rhat, c(NA, 1.2))
})

test_that("only stanimal blocks that did not mix are worth a warning", {
  reml <- .indri_diagnostics(list(list(converged = FALSE), list(converged = TRUE)),
                             rbind(1:2, 2:3))
  expect_null(.indri_rhat_warning(reml, 2L))
  stan <- .indri_diagnostics(list(list(converged = 1.2), list(converged = 1.01)),
                             rbind(1:2, 2:3))
  expect_match(.indri_rhat_warning(stan, 2L), "1 of 2 blocks have Rhat >= 1.05")
})

test_that("the stanimal engine fits blocks with lkj_prior = 1 unless told otherwise", {
  body <- paste(deparse(.indri_stanimal_fit), collapse = " ")
  expect_match(body, "lkj_prior = 1", fixed = TRUE)
})

test_that("the progress bar counts blocks on stderr, and is silent when off", {
  out <- capture.output(type = "message", {
    bar <- .indri_progress(5L, TRUE)
    for (i in 1:5) bar$tick(1L)
    bar$close()
  })
  text <- paste(out, collapse = "")
  expect_match(text, "fitting blocks", fixed = TRUE)
  expect_match(text, "5/5", fixed = TRUE)
  silent <- capture.output(type = "message", {
    bar <- .indri_progress(5L, FALSE)
    bar$tick(5L)
    bar$close()
  })
  expect_length(silent, 0L)
  expect_false(.indri_progress(5L, FALSE)$enabled)
})

test_that("a fit with the progress bar matches one without, in order", {
  skip_if_not_installed("sommer")
  data <- simulate_block_data(n_fam = 20, size = 5, K = 5L)
  for (cores in c(1L, 2L)) {
    quiet <- fit_G(data$Y, A = data$A, engine = "sommer", cores = cores,
                   n_draws = 20L, seed = 4L, progress = FALSE)
    shown <- capture.output(type = "message",
      loud <- fit_G(data$Y, A = data$A, engine = "sommer", cores = cores,
                    n_draws = 20L, seed = 4L, progress = TRUE))
    expect_match(paste(shown, collapse = ""), "fitting blocks", fixed = TRUE)
    expect_identical(loud$G, quiet$G, label = paste("cores", cores))
    expect_identical(loud$corrG_draws, quiet$corrG_draws)
    expect_identical(lapply(loud$blocks, `[[`, "G"),
                     lapply(quiet$blocks, `[[`, "G"))
  }
})
