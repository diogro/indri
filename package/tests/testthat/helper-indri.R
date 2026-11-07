skip_if_no_cmdstan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  ok <- tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) FALSE)
  if (!ok) testthat::skip("CmdStan is not installed")
}

skip_if_no_stanimal <- function() {
  testthat::skip_if_not_installed("stanimal")
}

skip_if_no_sommer <- function() {
  testthat::skip_if_not_installed("sommer")
}

skip_if_no_asreml <- function() {
  testthat::skip_if_not_installed("asreml")
  licenced <- tryCatch({
    suppressWarnings(asreml::asreml(
      fixed = y ~ 1, data = data.frame(y = c(1, 2))))
    TRUE
  }, error = function(e) FALSE)
  if (!licenced) testthat::skip("asreml is not licensed")
}

# A 4-trait covariance matrix with a clear structure, and the trait subsets that
# cover every pair. Blocks that are exact submatrices of this matrix are what a
# perfect engine would return, which is how the assembly steps are tested
# without fitting anything.
known_G <- function() {
  loadings <- matrix(c(1.0, 0.8, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.7), nrow = 4)
  loadings %*% t(loadings) + diag(c(0.3, 0.2, 0.4, 0.25))
}

covering_subsets <- function(K = 4L) {
  pairs <- t(utils::combn(K, 2L))
  pairs
}

blocks_from <- function(G, subsets) {
  lapply(seq_len(nrow(subsets)), function(i) G[subsets[i, ], subsets[i, ]])
}

min_eigen <- function(m) {
  min(eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values)
}

# Small balanced design: families of full sibs, so A is PSD without needing a
# pedigree, and a known genetic covariance to sample from.
simulate_block_data <- function(n_fam = 8, size = 5, K = 3L, G0 = NULL,
                                seed = 1) {
  set.seed(seed)
  K <- if (is.null(G0)) as.integer(K) else nrow(G0)
  if (is.null(G0)) {
    G0 <- matrix(0.3, K, K)
    diag(G0) <- 1
  }
  A <- kronecker(diag(n_fam), matrix(0.5, size, size))
  diag(A) <- 1
  N <- n_fam * size
  genetic <- t(chol(A)) %*% matrix(stats::rnorm(N * K), N, K) %*% chol(G0)
  Y <- genetic + matrix(stats::rnorm(N * K, sd = 0.8), N, K)
  colnames(Y) <- paste0("trait", seq_len(K))
  rownames(Y) <- paste0("id", seq_len(N))
  list(Y = Y, A = A, G0 = G0, N = N, K = K)
}
