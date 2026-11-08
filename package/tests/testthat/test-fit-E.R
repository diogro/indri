# fit_G() assembles the residual covariance E from the block fits, the same
# way as G but without draws: the association tests need both.

test_that("E is assembled from the block E, in the units of Y", {
  skip_if_not_installed("sommer")
  data <- simulate_block_data(n_fam = 30, size = 5, K = 4L)
  Y <- sweep(data$Y, 2L, c(1, 50, 0.02, 3), "*") + 100
  colnames(Y) <- c("wing", "tail", "beak", "leg")
  subsets <- rbind(c(2L, 4L, 1L), c(1L, 3L, 4L), c(2L, 3L, 1L))
  fit <- fit_G(Y, A = data$A, subsets = subsets, engine = "sommer",
               cores = 1L, n_draws = 30L, seed = 7L)

  E_blocks <- lapply(fit$blocks, function(b) list(G = unname(b$E)))
  avg <- indri_assemble(E_blocks, fit$subsets, 4L)$G
  E_std <- indri_project_psd(E_blocks, fit$subsets, 4L, init = avg)$matrix
  D <- outer(fit$y_sd, fit$y_sd)
  expect_equal(unname(fit$E_bar), unname(avg * D), tolerance = 1e-10)
  expect_equal(unname(fit$E), unname(E_std * D), tolerance = 1e-8)
  expect_equal(dimnames(fit$E), list(colnames(Y), colnames(Y)))
  expect_equal(fit$E, t(fit$E), tolerance = 1e-10)
  # positive definite, as the whitening in indri_gwas() needs
  expect_gt(min_eigen(fit$E), 0)
  # the blocks keep their standardised E
  expect_lt(max(diag(fit$blocks[[1]]$E)), 2)
})

test_that("E agrees with the full model's residual covariance", {
  skip_if_not_installed("sommer")
  # The sibship design separates G from E poorly, so E is compared with the
  # full four-trait REML fit, not with the simulated truth: both were up to
  # 0.5 from it on the diagonal. Measured gap over seeds 1-6: <= 1e-4.
  data <- simulate_block_data(n_fam = 60, size = 5, K = 4L, seed = 3)
  subsets <- rbind(c(2L, 4L, 1L), c(1L, 3L, 4L), c(2L, 3L, 1L))
  fit <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "sommer",
               cores = 1L, n_draws = 30L, seed = 7L)
  full <- fit_G(data$Y, A = data$A, subsets = matrix(1:4, 1L),
                assembly = "average", engine = "sommer", cores = 1L,
                n_draws = 30L, seed = 7L)
  expect_lt(max(abs(fit$E - full$E)) / mean(diag(full$E)), 1e-3)
})

test_that("E is floored to stay positive definite", {
  # blocks whose E average to a singular matrix still give a usable E
  E <- matrix(c(1, 1, 1, 1), 2)
  blocks <- list(list(G = E), list(G = E))
  out <- .indri_assemble_E(blocks, rbind(c(1L, 2L), c(1L, 2L)), 2L)
  expect_gt(min_eigen(out$E), 0)
  expect_equal(out$E_bar, E)
})
