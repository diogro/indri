test_that("inflate places the block and leaves the rest alone", {
  small <- matrix(1:4, 2, 2)
  big <- indri_inflate(small, index = c(1L, 3L), K = 3L, fill = 0)
  expect_equal(big[c(1, 3), c(1, 3)], small)
  expect_equal(big[2, ], c(0, 0, 0))
  expect_true(all(is.na(indri_inflate(small, c(1L, 3L), 3L)[2, ])))
})

test_that("exact blocks average back to the matrix they came from", {
  G <- known_G()
  subsets <- covering_subsets()
  blocks <- lapply(blocks_from(G, subsets), function(g) list(G = g))

  assembled <- indri_assemble(blocks, subsets)

  expect_equal(assembled$G, G, tolerance = 1e-12)
  expect_true(all(assembled$counts[upper.tri(assembled$counts)] >= 1))
  expect_equal(diag(assembled$counts), tabulate(as.vector(subsets), nbins = 4L))
})

test_that("uncovered entries stay missing rather than silently zero", {
  G <- known_G()
  subsets <- matrix(c(1L, 2L, 3L, 4L), nrow = 2, byrow = TRUE)
  blocks <- lapply(blocks_from(G, subsets), function(g) list(G = g))

  assembled <- indri_assemble(blocks, subsets)

  expect_false(is.na(assembled$G[1, 2]))
  expect_true(is.na(assembled$G[1, 3]))
  expect_equal(assembled$counts[1, 3], 0L)
})

test_that("inverse-variance weights favour the more precise block", {
  subsets <- matrix(c(1L, 2L, 1L, 2L), nrow = 2, byrow = TRUE)
  blocks <- list(
    list(G = matrix(c(1, 0.9, 0.9, 1), 2), se = matrix(0.01, 2, 2)),
    list(G = matrix(c(1, 0.1, 0.1, 1), 2), se = matrix(1, 2, 2))
  )

  equal <- indri_assemble(blocks, subsets, weights = "none")
  weighted <- indri_assemble(blocks, subsets, weights = "inverse_variance")

  expect_equal(equal$G[1, 2], 0.5)
  # the precise block says 0.9, the vague one says 0.1: weighting should land
  # near the precise estimate
  expect_gt(weighted$G[1, 2], 0.89)
  expect_lt(weighted$G[1, 2], 0.91)
})

test_that("draws are averaged the same way the means are", {
  G <- known_G()
  subsets <- covering_subsets()
  D <- 5L
  draws <- lapply(blocks_from(G, subsets), function(g) {
    array(rep(g, D), dim = c(nrow(g), ncol(g), D))
  })

  assembled <- indri_assemble_draws(draws, subsets)

  expect_equal(dim(assembled), c(4L, 4L, D))
  expect_equal(assembled[, , 1L], G, tolerance = 1e-12)
})
