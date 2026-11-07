# M3 timing: what does sharing the eigendecomposition of A save per block?
#
# Fits one six-trait block (the size indri uses at s = 4) at increasing N with
# short chains, so the fixed per-fit costs are visible, three ways:
#   old     stanimal decomposes A in Stan, once per chain (stanimal <= 0.7.4)
#   new     stanimal decomposes A once per fit, in R
#   shared  a precomputed decomposition is passed as eigen_A (what indri does)
#
# `old` needs a library holding stanimal master; the other two need stanimal
# with `eigen_A`. Each mode runs in its own process, so the two can coexist:
#   INDRI_TIMING_MODE=old INDRI_TIMING_LIB=/path/to/master/lib Rscript examples/eigen_A_timing.R
#   INDRI_TIMING_MODE=new    Rscript examples/eigen_A_timing.R
#   INDRI_TIMING_MODE=shared Rscript examples/eigen_A_timing.R
#
# Environment: INDRI_TIMING_N (500,1000,2000,4000), INDRI_TIMING_ITER (200),
# INDRI_TIMING_REPS (2). Writes results/eigen_A_timing_<mode>.csv (ignored).

mode <- Sys.getenv("INDRI_TIMING_MODE", "shared")
lib <- Sys.getenv("INDRI_TIMING_LIB", "")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
library(stanimal)
if (mode %in% c("new", "shared") &&
    !"eigen_A" %in% names(formals(stanimal::lmm_animal)))
  stop("mode `", mode, "` needs a stanimal with `eigen_A`")
if (mode == "old" && "eigen_A" %in% names(formals(stanimal::lmm_animal)))
  stop("mode `old` needs stanimal master; set INDRI_TIMING_LIB")

Ns <- as.integer(strsplit(Sys.getenv("INDRI_TIMING_N", "500,1000,2000,4000"), ",")[[1]])
iter <- as.integer(Sys.getenv("INDRI_TIMING_ITER", "200"))
reps <- as.integer(Sys.getenv("INDRI_TIMING_REPS", "2"))
K <- 6L
dir.create("results", showWarnings = FALSE)

rows <- list()
for (N in Ns) {
  set.seed(N)
  size <- 5L
  A <- kronecker(diag(N %/% size), matrix(0.5, size, size))
  diag(A) <- 1
  G <- stanimal::random_covariance(K, method = "lkj", min_var = 0.4, max_var = 0.6)
  E <- stanimal::random_covariance(K, method = "lkj", min_var = 0.4, max_var = 0.6)
  Y <- t(chol(A)) %*% matrix(rnorm(N * K), N, K) %*% chol(G) +
    matrix(rnorm(N * K), N, K) %*% chol(E)

  decompose <- system.time({
    e <- eigen(A, symmetric = TRUE)
    eigen_A <- list(Q = e$vectors, lambda = e$values)
  })[["elapsed"]]

  for (r in seq_len(reps)) {
    args <- list(Y = Y, X = matrix(1, N, 1L), A = A, chains = 4L, cores = 4L,
                 iter = iter, seed = r, init = 0.5, refresh = 0,
                 show_messages = FALSE)
    if (mode == "shared") args$eigen_A <- eigen_A
    elapsed <- system.time(fit <- do.call(stanimal::lmm_animal, args))[["elapsed"]]
    sampling <- max(fit$time()$chains$total)
    rows[[length(rows) + 1L]] <- data.frame(
      mode = mode, N = N, rep = r, iter = iter, wall = elapsed,
      slowest_chain = sampling, overhead = elapsed - sampling,
      r_eigen = decompose)
    message(sprintf("%-6s N = %5d rep %d: %6.1f s wall, %6.1f s slowest chain",
                    mode, N, r, elapsed, sampling))
    rm(fit); invisible(gc(FALSE))
  }
}

out <- do.call(rbind, rows)
utils::write.csv(out, sprintf("results/eigen_A_timing_%s.csv", mode),
                 row.names = FALSE)
print(aggregate(cbind(wall, slowest_chain, overhead, r_eigen) ~ mode + N,
                data = out, FUN = median))
