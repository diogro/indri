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

test_that("E recovers the residual covariance that made the data", {
  skip_if_not_installed("sommer")
  # simulate_block_data() draws residuals with sd 0.8, independent across traits
  data <- simulate_block_data(n_fam = 60, size = 5, K = 4L, seed = 3)
  fit <- fit_G(data$Y, A = data$A, engine = "sommer", cores = 1L,
               n_draws = 30L, seed = 7L)
  expect_lt(max(abs(diag(fit$E) - 0.64)), 0.2)
  off <- fit$E[upper.tri(fit$E)]
  expect_lt(max(abs(off)), 0.15)
})

test_that("E is floored to stay positive definite", {
  # blocks whose E average to a singular matrix still give a usable E
  E <- matrix(c(1, 1, 1, 1), 2)
  blocks <- list(list(G = E), list(G = E))
  out <- .indri_assemble_E(blocks, rbind(c(1L, 2L), c(1L, 2L)), 2L)
  expect_gt(min_eigen(out$E), 0)
  expect_equal(out$E_bar, E)
})
