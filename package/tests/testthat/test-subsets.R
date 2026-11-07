test_that("subsets are valid indices with every pair covered", {
  subsets <- indri_subsets(K = 8, s = 3, min_cooccurrence = 3, seed = 1)

  expect_true(is.matrix(subsets))
  expect_type(subsets, "integer")
  expect_true(all(subsets >= 1 & subsets <= 8))
  expect_true(all(apply(subsets, 1L, function(row) !anyDuplicated(row))))
  expect_gte(ncol(subsets), 2L)

  coverage <- attr(subsets, "coverage")
  expect_equal(dim(coverage), c(8L, 8L))
  expect_true(all(coverage[upper.tri(coverage)] >= 3))
  # the diagonal counts how often a trait was sampled, not how many subsets exist
  expect_true(all(diag(coverage) >= 1))
  expect_true(all(diag(coverage) <= nrow(subsets)))
})

test_that("every pair reaches min_cooccurrence across seeds and settings", {
  # The LiMMBo construction alone misses this for about half of these
  # settings (always for s = 2), so a single seed is not a test.
  for (K in c(4, 6, 8, 12, 20, 30)) for (s in 2:4) for (c in 1:5) {
    if (s > K) next
    for (seed in 1:5) {
      subsets <- indri_subsets(K, s, min_cooccurrence = c, seed = seed)
      coverage <- attr(subsets, "coverage")
      info <- sprintf("K = %d, s = %d, c = %d, seed = %d", K, s, c, seed)
      # distinct subsets only, and the cap a block width imposes
      expect_false(anyDuplicated(apply(subsets, 1L, paste, collapse = "-")) > 0,
                   label = info)
      cap <- choose(K - 2, ncol(subsets) - 2)
      expect_gte(min(coverage[upper.tri(coverage)]), min(c, cap), label = info)
      expect_equal(attr(subsets, "min_cooccurrence"), min(c, cap), label = info)
      expect_true(all(apply(subsets, 1L, function(r) !anyDuplicated(r))),
                  label = info)
      expect_identical(coverage, indri_pair_coverage(subsets, K), label = info)
    }
  }
})

test_that("subsets are reproducible and do not disturb the caller's stream", {
  set.seed(99)
  before <- runif(1)

  set.seed(99)
  a <- indri_subsets(K = 6, s = 3, seed = 7)
  after <- runif(1)
  b <- indri_subsets(K = 6, s = 3, seed = 7)
  c <- indri_subsets(K = 6, s = 3, seed = 8)

  expect_equal(before, after)
  expect_identical(a, b)
  expect_false(identical(a, c))
})

test_that("the subset width follows from the pair windows, not from s", {
  subsets <- indri_subsets(K = 10, s = 3, min_cooccurrence = 1, seed = 2)
  expect_equal(ncol(subsets), 4L)

  subsets <- indri_subsets(K = 10, s = 4, min_cooccurrence = 1, seed = 2)
  expect_equal(ncol(subsets), 6L)
})

test_that("a larger minimum co-occurrence buys more subsets", {
  few <- indri_subsets(K = 6, s = 3, min_cooccurrence = 1, seed = 3)
  many <- indri_subsets(K = 6, s = 3, min_cooccurrence = 6, seed = 3)
  expect_gt(nrow(many), nrow(few))
})

test_that("degenerate inputs are refused", {
  expect_error(indri_subsets(K = 1, s = 2), "K")
  expect_error(indri_subsets(K = 4, s = 1), "s")
  expect_error(indri_subsets(K = 4, s = 5), "cannot exceed")
  expect_error(indri_subsets(K = 4, s = 2, min_cooccurrence = 0),
               "min_cooccurrence")
  expect_error(indri_subsets(K = 4.5, s = 2), "whole number")
})

test_that("coverage can be recomputed from a subset matrix", {
  subsets <- indri_subsets(K = 5, s = 2, min_cooccurrence = 2, seed = 4)
  expect_identical(indri_pair_coverage(subsets), attr(subsets, "coverage"))
})
