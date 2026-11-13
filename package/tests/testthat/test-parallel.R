# The core budget, and the parallel block path it drives.
#
# The budget exists because the engines bring their own parallelism: stanimal
# defaults to four chains and runs them on four cores, so "cores = 4" at the
# fit level must not mean four blocks at once on top of that.

test_that("the core budget is never exceeded", {
  engine <- indri_engine("stanimal")
  cases <- list(
    list(cores = 1L, args = list(), blocks = 1L, per_block = 1L),
    list(cores = 4L, args = list(), blocks = 1L, per_block = 4L),
    list(cores = 4L, args = list(chains = 1L), blocks = 4L, per_block = 1L),
    list(cores = 8L, args = list(chains = 2L), blocks = 4L, per_block = 2L),
    list(cores = 6L, args = list(chains = 2L), blocks = 3L, per_block = 2L),
    list(cores = 2L, args = list(chains = 4L), blocks = 1L, per_block = 2L),
    list(cores = 8L, args = list(chains = 4L, cores = 2L), blocks = 4L,
         per_block = 2L)
  )
  for (case in cases) {
    plan <- indri:::.indri_core_plan(case$cores, case$args, engine, 23L)
    expect_equal(plan$blocks_in_parallel, case$blocks)
    expect_equal(plan$per_block_cores, case$per_block)
    expect_lte(plan$blocks_in_parallel * plan$per_block_cores, case$cores)
    # the statistical setting is passed through untouched
    expect_equal(plan$chains, case$args$chains %||% 4L)
  }
})

test_that("a budget smaller than the chains runs the chains in turn", {
  engine <- indri_engine("stanimal")
  plan <- indri:::.indri_core_plan(1L, list(), engine, 10L)
  expect_equal(plan$chains, 4L)
  expect_equal(plan$per_block_cores, 1L)
  expect_equal(plan$blocks_in_parallel, 1L)
})

test_that("the number of blocks caps the number of workers", {
  engine <- indri_engine("stanimal")
  plan <- indri:::.indri_core_plan(32L, list(chains = 1L), engine, 3L)
  expect_equal(plan$blocks_in_parallel, 3L)
})

test_that("parallel blocks give the same answer as sequential ones", {
  # R CMD check allows two worker processes; this test needs four
  skip_if(nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")), "check limits cores")
  skip_if_no_cmdstan()
  skip_if_no_stanimal()
  data <- simulate_block_data(n_fam = 8, size = 5, K = 4L)
  subsets <- t(utils::combn(4L, 2L))
  engine_args <- list(chains = 1L, iter = 200L, seed = 1L, init = 0.5,
                      refresh = 0, show_messages = FALSE)

  sequential <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "stanimal",
                      engine_args = engine_args, cores = 1L)
  parallel <- fit_G(data$Y, A = data$A, subsets = subsets, engine = "stanimal",
                    engine_args = engine_args, cores = 4L)

  expect_equal(parallel$parallel$blocks_in_parallel, 4L)
  # the same seed gives the same draws per block, so the assembled matrix must
  # match exactly, not approximately
  expect_equal(parallel$G, sequential$G, tolerance = 1e-12)
  expect_equal(parallel$corrG, sequential$corrG, tolerance = 1e-12)
})
