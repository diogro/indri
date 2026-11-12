# The file writers and readers run everywhere; fits need the WOMBAT binary
# (Linux or Windows) and skip without it.

test_that("the wombat engine is registered with REML capabilities", {
  expect_true("wombat" %in% indri_engines())
  engine <- indri_engine("wombat")
  expect_equal(engine$name, "wombat")
  expect_equal(engine$capabilities$a_form, "sparse_inverse")
  expect_false(engine$capabilities$draws)
  expect_true(engine$capabilities$vcov)
  expect_true(engine$capabilities$missing)
})

test_that("a missing WOMBAT binary is a clear error", {
  expect_error(.indri_wombat_binary(tempfile()), "WOMBAT program")
})

test_that("covariables are the columns of X beyond the constant", {
  x <- c(0.3, -1, 2, 0.5)
  expect_equal(ncol(.indri_wombat_covariables(cbind(1, x))), 1L)
  expect_equal(.indri_wombat_covariables(cbind(1, x))[, 1], x)
  # indicators of every level span the constant: one of them is dropped
  f <- cbind(a = c(1, 1, 0, 0), b = c(0, 0, 1, 1))
  expect_equal(ncol(.indri_wombat_covariables(f)), 1L)
  expect_equal(ncol(.indri_wombat_covariables(matrix(1, 4, 1))), 0L)
})

test_that("data records are sorted by individual, then trait, without gaps", {
  Y <- matrix(c(1, NA, 3, 4, 5, NA), 3, 2)
  path <- tempfile()
  .indri_wombat_write_data(path, Y, matrix(0, 3, 0))
  rows <- read.table(path)
  expect_equal(rows[[1]], c(1, 2, 2, 1))          # trait
  expect_equal(rows[[2]], c(1, 1, 2, 3))          # individual
  expect_equal(rows[[4]], c(1, 4, 5, 3))          # response
})

test_that("the GIN file holds column <= row and a .codes line per level", {
  A <- kronecker(diag(2), matrix(c(1, 0.5, 0.5, 1), 2))
  gin <- tempfile(); codes <- tempfile()
  .indri_wombat_write_gin(gin, codes, indri_a_forms(A, "sparse_inverse")$
                            sparse_inverse, 4L)
  # the header is log |A| (not of its inverse)
  expect_equal(as.numeric(readLines(gin, n = 1L)),
               as.numeric(determinant(A)$modulus), tolerance = 1e-10)
  x <- read.table(gin, skip = 1)
  expect_true(all(x[[1]] <= x[[2]]))
  M <- matrix(0, 4, 4)
  M[cbind(x[[2]], x[[1]])] <- x[[3]]
  M[cbind(x[[1]], x[[2]])] <- x[[3]]
  expect_equal(M, solve(A), tolerance = 1e-12)
  expect_equal(nrow(read.table(codes)), 4L)
})

test_that("BestPoint and AvInfoCovs are read into indri's orders", {
  dir <- tempfile(); dir.create(dir)
  # E then G, upper triangles row by row, after a header line
  writeLines(c("  -10.5  6  100", "1 0.1 2", "3 0.4 5"),
             file.path(dir, "BestPoint"))
  point <- .indri_wombat_read_bestpoint(file.path(dir, "BestPoint"), 2L)
  expect_equal(point$loglik, -10.5)
  expect_equal(point$E, matrix(c(1, 0.1, 0.1, 2), 2))
  expect_equal(point$G, matrix(c(3, 0.4, 0.4, 5), 2))

  # a 3-trait information matrix: residual (6) then G (6); its inverse's G
  # block must come back in .indri_vech_index() order
  set.seed(1)
  m <- crossprod(matrix(rnorm(144), 12)) + diag(12)
  index <- which(upper.tri(m, diag = TRUE), arr.ind = TRUE)
  writeLines(c("-10.5", sprintf("%d %d %.17g", index[, 1], index[, 2],
                                m[index])), file.path(dir, "AvInfoCovs"))
  v <- .indri_wombat_read_vcov(file.path(dir, "AvInfoCovs"), 3L, -10.5)
  inverse <- solve(m)[7:12, 7:12]
  # upper row-wise: G11 G12 G13 G22 G23 G33; vech: G11 G21 G22 G31 G32 G33
  order <- c(1, 2, 4, 3, 5, 6)
  expect_equal(v, unname(inverse[order, order]), tolerance = 1e-10)
  # an information matrix for another point is not used
  expect_null(.indri_wombat_read_vcov(file.path(dir, "AvInfoCovs"), 3L, -9))
})

test_that("the parameter file lifts WOMBAT's limits only when needed", {
  path <- tempfile()
  .indri_wombat_write_par(path, 3L, 10L, 0L, diag(3), diag(3))
  expect_false(any(grepl("^OVER", readLines(path))))
  .indri_wombat_write_par(path, 60L, 10L, 0L, diag(60), diag(60))
  lines <- readLines(path)
  expect_true("OVER MAXTRAIT 60" %in% lines)
  # before the analysis type and the data layout
  expect_lt(match("OVER MAXTRAIT 60", lines), grep("^ANAL", lines))
})

test_that("a reduced rank writes a PC analysis with a full-rank residual", {
  path <- tempfile()
  .indri_wombat_write_par(path, 4L, 10L, 0L, diag(4), diag(4), rank = 2L)
  lines <- readLines(path)
  expect_true(all(c("ANAL MUV PC 4", "VAR residual 4 4", "VAR animal 4 2") %in%
                    lines))
  .indri_wombat_write_par(path, 4L, 10L, 0L, diag(4), diag(4))
  expect_true(all(c("ANAL MUV 4", "VAR residual 4", "VAR animal 4") %in%
                    readLines(path)))
})

test_that("a crashed WOMBAT is an error even when it left a BestPoint", {
  skip_on_os("windows")
  # a stand-in that writes a BestPoint from an "EM round", then crashes
  fake <- tempfile()
  writeLines(c("#!/bin/sh",
               "printf ' -10 6 100\\n1 0 1\\n1 0 1\\n' > BestPoint",
               "echo 'forrtl: severe (174): SIGSEGV'", "exit 174"), fake)
  Sys.chmod(fake, "755")
  Y <- matrix(rnorm(8), 4, 2)
  A <- diag(4)
  expect_error(.indri_wombat_fit(Y, matrix(1, 4, 1),
                                 indri_a_forms(A, "sparse_inverse"),
                                 list(wombat = fake)),
               "status 174")
})

test_that("convergence comes from WOMBAT's own verdict", {
  path <- tempfile()
  writeLines(c(" AI    1   -234.477   0.2525",
               " *** AIREML seems to have converged: Change in log L < 5E-04"),
             path)
  expect_true(.indri_wombat_converged(path))
  writeLines(" AI    1   -234.477   0.2525", path)
  expect_false(.indri_wombat_converged(path))
  expect_false(.indri_wombat_converged(tempfile()))
})

test_that("the wombat engine agrees with sommer, with covariates and gaps", {
  skip_if_no_wombat()
  skip_if_no_sommer()
  data <- simulate_block_data(n_fam = 50, size = 6, K = 3L, seed = 2)
  set.seed(3)
  x <- rnorm(data$N)
  Y <- scale(data$Y + outer(x, c(0.5, -0.3, 0.2)))
  Y[sample(length(Y), round(0.1 * length(Y)))] <- NA
  # names that do not sort in column order
  colnames(Y) <- c("trait2", "trait10", "trait1")
  X <- cbind("(Intercept)" = 1, x = x)
  rownames(data$A) <- colnames(data$A) <- rownames(Y)

  w <- .indri_wombat_fit(Y, X, indri_a_forms(data$A, "sparse_inverse"),
                         list(wombat = getOption("indri.wombat")))
  s <- .indri_sommer_fit(Y, X, indri_a_forms(data$A, "inverse"))
  expect_true(w$converged)
  expect_equal(rownames(w$G), colnames(Y))
  expect_equal(unname(w$G), unname(s$G), tolerance = 1e-3)
  expect_equal(unname(w$E), unname(s$E), tolerance = 1e-3)
  expect_equal(w$vcov, s$vcov, tolerance = 1e-3)
  expect_equal(unname(w$se), unname(s$se), tolerance = 1e-3)
})

test_that("the wombat engine fits a single trait", {
  skip_if_no_wombat()
  data <- simulate_block_data(n_fam = 30, size = 5, K = 1L)
  Y <- scale(data$Y)
  out <- indri_engine("wombat")$fit(Y, matrix(1, data$N, 1L),
                                    indri_a_forms(data$A, "sparse_inverse"))
  expect_equal(dim(out$G), c(1L, 1L))
  expect_gt(out$G[1, 1], 0)
  expect_true(out$converged)
})

test_that("residual covariances of traits never recorded together stay zero", {
  skip_if_no_wombat()
  # two traits in two groups: character states never observed together
  data <- simulate_block_data(n_fam = 60, size = 6, K = 2L, seed = 4)
  group <- rep(1:2, length.out = data$N)
  Y <- matrix(NA_real_, data$N, 4,
              dimnames = list(NULL, c("a:1", "b:1", "a:2", "b:2")))
  Y[group == 1, 1:2] <- data$Y[group == 1, ]
  Y[group == 2, 3:4] <- data$Y[group == 2, ]
  out <- indri_engine("wombat")$fit(scale(Y), matrix(1, data$N, 1L),
                                    indri_a_forms(data$A, "sparse_inverse"))
  expect_equal(unname(out$E[1:2, 3:4]), matrix(0, 2, 2))
  expect_true(all(out$E[1:2, 1:2] != 0))
  # the genetic covariances between groups are estimated
  expect_true(all(out$G[1:2, 3:4] != 0))
  expect_false(is.null(out$vcov))
})

test_that("fit_G runs end to end with the wombat engine in parallel", {
  skip_if_no_wombat()
  data <- simulate_block_data(n_fam = 40, size = 5, K = 4L)
  fit <- fit_G(data$Y, A = data$A, engine = "wombat", cores = 2L,
               progress = FALSE)
  expect_s3_class(fit, "indri_fit")
  expect_equal(fit$engine, "wombat")
  expect_gt(min_eigen(fit$G), -1e-8)
  expect_equal(unname(rownames(fit$G)), colnames(data$Y))
})

test_that("the .matrix file holds the lower triangle of A itself", {
  A <- kronecker(diag(2), matrix(c(1, 0.5, 0.5, 1), 2))
  path <- tempfile()
  .indri_wombat_write_matrix(path, indri_a_forms(A, "sparse_inverse")$
                               sparse_inverse, 4L)
  x <- read.table(path)
  expect_true(all(x[[1]] <= x[[2]]))
  M <- matrix(0, 4, 4)
  M[cbind(x[[2]], x[[1]])] <- x[[3]]
  M[cbind(x[[1]], x[[2]])] <- x[[3]]
  expect_equal(M, A, tolerance = 1e-12)
})
