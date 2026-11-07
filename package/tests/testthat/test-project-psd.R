test_that("exact blocks are recovered by the projection", {
  G <- known_G()
  subsets <- covering_subsets()
  blocks <- lapply(blocks_from(G, subsets), function(g) list(G = g))

  fit <- indri_project_psd(blocks, subsets)

  expect_equal(fit$matrix, G, tolerance = 1e-5)
  expect_equal(fit$convergence, 0L)
  expect_gt(min_eigen(fit$matrix), -1e-8)
})

test_that("the projection returns a positive semi-definite matrix from an indefinite average", {
  # a frustrated triangle: three pairwise estimates that cannot all hold, so
  # their average has a negative eigenvalue
  subsets <- t(utils::combn(3L, 2L))
  blocks <- list(
    list(G = matrix(c(1, 0.9, 0.9, 1), 2)),
    list(G = matrix(c(1, -0.9, -0.9, 1), 2)),
    list(G = matrix(c(1, 0.9, 0.9, 1), 2))
  )

  assembled <- indri_assemble(blocks, subsets, K = 3L)
  expect_lt(min_eigen(assembled$G), 0)

  fit <- indri_project_psd(blocks, subsets, K = 3L)
  expect_gt(min_eigen(fit$matrix), -1e-6)
  # and it did not simply zero the offending correlations away
  expect_gt(sum(abs(fit$matrix[upper.tri(fit$matrix)])), 0)
})

test_that("the fit is closer to the blocks than the average it started from", {
  G <- known_G()
  subsets <- covering_subsets()
  # perturb the blocks so that their average is not PSD but the truth is
  set.seed(1)
  blocks <- lapply(blocks_from(G, subsets), function(g) {
    list(G = g + matrix(stats::rnorm(length(g), sd = 0.2), nrow(g)))
  })
  assembled <- indri_assemble(blocks, subsets)
  loss <- function(C) {
    sum(vapply(seq_along(blocks), function(i) {
      idx <- subsets[i, ]
      sum((C[idx, idx, drop = FALSE] - blocks[[i]]$G)^2)
    }, numeric(1)))
  }

  fit <- indri_project_psd(blocks, subsets, init = assembled$G)

  # The average need not be a covariance matrix, so it is not a fair
  # competitor: the projection minimises over PSD matrices only. Compare
  # against the PSD matrix the optimiser starts from instead.
  start <- t(chol(assembled$G + 1e-8 * diag(4)))
  expect_lte(loss(fit$matrix), loss(start %*% t(start)) + 1e-8)
})

test_that("the analytical gradient matches finite differences", {
  G <- known_G()
  subsets <- covering_subsets()
  set.seed(2)
  blocks <- lapply(blocks_from(G, subsets), function(g) {
    list(G = g + matrix(stats::rnorm(length(g), sd = 0.1), nrow(g)))
  })
  K <- 4L
  parts <- .indri_psd_parts(blocks, subsets, K)
  theta <- stats::runif(nrow(parts$lt), 0.1, 1)

  ana <- .indri_psd_gradient(theta, K, parts)
  eps <- 1e-6
  num <- vapply(seq_along(theta), function(k) {
    up <- theta; up[k] <- up[k] + eps
    dn <- theta; dn[k] <- dn[k] - eps
    (.indri_psd_objective(up, K, parts) - .indri_psd_objective(dn, K, parts)) /
      (2 * eps)
  }, numeric(1))
  expect_equal(ana, num, tolerance = 1e-4)
})

test_that("missing entries in the average do not break the start", {
  subsets <- matrix(c(1L, 2L, 3L, 4L), nrow = 2, byrow = TRUE)
  blocks <- list(list(G = diag(2)), list(G = diag(2) * 2))
  fit <- indri_project_psd(blocks, subsets, K = 4L)
  expect_false(anyNA(fit$matrix))
  expect_gt(min_eigen(fit$matrix), -1e-8)
})
