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

test_that("inverse-variance weights fall back to posterior precision for draws", {
  subsets <- matrix(c(1L, 2L, 1L, 2L), nrow = 2, byrow = TRUE)
  make_draws <- function(g, sd, D = 400L) {
    out <- array(0, c(2L, 2L, D))
    for (i in 1:2) for (j in 1:2) out[i, j, ] <- stats::rnorm(D, g[i, j], sd)
    out
  }
  set.seed(42)
  blocks <- list(
    list(G = matrix(c(1, 0.9, 0.9, 1), 2),
         draws = make_draws(matrix(c(1, 0.9, 0.9, 1), 2), 0.01)),
    list(G = matrix(c(1, 0.1, 0.1, 1), 2),
         draws = make_draws(matrix(c(1, 0.1, 0.1, 1), 2), 1))
  )

  weighted <- indri_assemble(blocks, subsets, weights = "inverse_variance")

  # the tight-posterior block says 0.9, the loose one says 0.1
  expect_gt(weighted$G[1, 2], 0.85)
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

test_that("mixture draws copy a covering block's value entry by entry", {
  subsets <- rbind(c(1L, 2L, 3L), c(2L, 3L, 4L))
  D <- 200L
  # block 1 says 1 everywhere, block 2 says 2 everywhere
  draws <- list(array(1, c(3L, 3L, D)), array(2, c(3L, 3L, D)))

  mixed <- indri_assemble_draws(draws, subsets, K = 4L, method = "mixture",
                                seed = 1L)

  expect_equal(dim(mixed), c(4L, 4L, D))
  expect_true(all(mixed[1, 1, ] == 1))           # only block 1 covers (1, 1)
  expect_true(all(mixed[4, 4, ] == 2))           # only block 2 covers (4, 4)
  expect_true(all(is.na(mixed[1, 4, ])))         # nobody covers (1, 4)
  shared <- mixed[2, 3, ]                        # both blocks cover (2, 3)
  expect_true(all(shared %in% c(1, 2)))
  expect_true(mean(shared == 2) > 0.3 && mean(shared == 2) < 0.7)
  for (d in seq_len(D)) expect_equal(mixed[, , d], t(mixed[, , d]))
})

test_that("mixture draws are reproducible and leave the caller's stream alone", {
  subsets <- rbind(c(1L, 2L), c(1L, 2L))
  draws <- list(array(stats::rnorm(40), c(2L, 2L, 10L)),
                array(stats::rnorm(40), c(2L, 2L, 10L)))
  set.seed(5); before <- runif(1)
  set.seed(5)
  a <- indri_assemble_draws(draws, subsets, method = "mixture", seed = 9L)
  after <- runif(1)
  b <- indri_assemble_draws(draws, subsets, method = "mixture", seed = 9L)
  expect_identical(a, b)
  expect_equal(before, after)
})

test_that("mixing keeps the block spread that averaging divides away", {
  # c blocks that all saw the same data: their draws share one error term, the
  # situation in which averaging over independent chains is overconfident.
  set.seed(3)
  n_blocks <- 9L
  D <- 4000L
  subsets <- matrix(rep(1:2, n_blocks), n_blocks, byrow = TRUE)
  draws <- lapply(seq_len(n_blocks), function(j) {
    array(stats::rnorm(4L * D), c(2L, 2L, D))   # block posterior sd 1
  })

  averaged <- indri_assemble_draws(draws, subsets, method = "average")
  mixed <- indri_assemble_draws(draws, subsets, method = "mixture", seed = 1L)

  expect_equal(stats::sd(averaged[1, 2, ]), 1 / sqrt(n_blocks), tolerance = 0.1)
  expect_equal(stats::sd(mixed[1, 2, ]), 1, tolerance = 0.1)
})
