# Block-size and co-occurrence sweep with a known genetic truth.
#
#   Rscript examples/block_size_sweep.R
#
# Fits the same simulated data once per (s, c) combination, so configurations
# are compared on identical data; run several seeds as separate processes and
# compare within seed.
#
#   INDRI_SWEEP_ENGINE  block engine: stanimal, sommer or asreml (stanimal)
#   INDRI_SWEEP_K       traits (50)
#   INDRI_SWEEP_SEED    simulation seed (1)
#   INDRI_SWEEP_CORES   process budget passed to fit_G() (128)
#   INDRI_SWEEP_ITER    stanimal iterations per chain (300)
#   INDRI_SWEEP_CHAINS  stanimal chains per block (2)
#   INDRI_SWEEP_S       comma-separated s values to sweep (2,3,4,5,6)
#   INDRI_SWEEP_C       comma-separated min_cooccurrence values (3)
#   INDRI_SWEEP_METHOD  truth method for random_covariance (partial)
#   INDRI_SWEEP_LKJ     lkj_prior when METHOD = "lkj" (2)
#
# Writes results/block_size_sweep_<engine>_K<K>_seed<seed>.csv (git-ignored).

library(indri)

engine <- Sys.getenv("INDRI_SWEEP_ENGINE", "stanimal")
K <- as.integer(Sys.getenv("INDRI_SWEEP_K", "50"))
seed <- as.integer(Sys.getenv("INDRI_SWEEP_SEED", "1"))
cores <- as.integer(Sys.getenv("INDRI_SWEEP_CORES", "128"))
iter <- as.integer(Sys.getenv("INDRI_SWEEP_ITER", "300"))
chains <- as.integer(Sys.getenv("INDRI_SWEEP_CHAINS", "2"))
s_values <- as.integer(strsplit(Sys.getenv("INDRI_SWEEP_S", "2,3,4,5,6"), ",")[[1]])
c_values <- as.integer(strsplit(Sys.getenv("INDRI_SWEEP_C", "3"), ",")[[1]])
method <- Sys.getenv("INDRI_SWEEP_METHOD", "partial")
lkj_prior <- as.numeric(Sys.getenv("INDRI_SWEEP_LKJ", "2"))
n_fam <- 60L
size <- 5L
# sommer multithreads each fit through OpenMP; let `cores` drive parallelism
Sys.setenv(OMP_NUM_THREADS = "1")
dir.create("results", showWarnings = FALSE)

set.seed(seed)
truth_args <- list(n_traits = K, method = method, min_var = 0.8, max_var = 1.4,
                   trait_names = paste0("trait", seq_len(K)))
if (method == "lkj") truth_args$lkj_prior <- lkj_prior
truth <- do.call(stanimal::random_covariance, truth_args)
resid <- stanimal::random_covariance(
  n_traits = K, method = "lkj", lkj_prior = 2, min_var = 0.5, max_var = 0.9)
A <- kronecker(diag(n_fam), matrix(0.5, size, size)); diag(A) <- 1
N <- nrow(A)
genetic <- t(chol(A)) %*% matrix(stats::rnorm(N * K), N, K) %*% chol(truth)
Y <- genetic + matrix(stats::rnorm(N * K), N, K) %*% chol(resid)
rownames(Y) <- paste0("id", seq_len(N))
colnames(Y) <- colnames(truth)
corr_true <- stats::cov2cor(truth)
pairs <- upper.tri(corr_true)

engine_args <- if (engine == "stanimal") {
  list(chains = chains, iter = iter, init = 0.5, refresh = 0,
       show_messages = FALSE)
} else list()

rows <- list()
for (s in s_values) {
  for (c in c_values) {
    message("== ", engine, ", K = ", K, ", seed ", seed, ", s = ", s,
            ", c = ", c, " ==")
    started <- Sys.time()
    fit <- fit_G(Y, A = A, s = s, min_cooccurrence = c, engine = engine,
                 engine_args = engine_args, cores = cores,
                 seed = 1000L + seed)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    error <- fit$corrG[pairs] - corr_true[pairs]
    converged <- vapply(fit$blocks, function(b) {
      v <- b$converged
      if (is.logical(v)) isTRUE(v) else is.numeric(v) && isTRUE(v < 1.05)
    }, logical(1))
    rows[[length(rows) + 1L]] <- data.frame(
      engine = engine, K = K, seed = seed, s = s, c = c,
      block_size = ncol(fit$subsets), n_blocks = nrow(fit$subsets),
      distinct_blocks = length(unique(apply(fit$subsets, 1L, paste,
                                            collapse = "-"))),
      converged = mean(converged),
      time_s = round(elapsed, 2),
      block_time_s = round(sum(vapply(fit$blocks, function(b)
        if (is.null(b$time)) NA_real_ else b$time, numeric(1))), 2),
      median_abs = stats::median(abs(error)),
      q90_abs = as.numeric(stats::quantile(abs(error), 0.9)),
      relative_frobenius = norm(fit$corrG - corr_true, "F") /
        norm(corr_true, "F"),
      corr_of_corr = stats::cor(fit$corrG[pairs], corr_true[pairs]),
      min_eig = min(eigen(fit$G, symmetric = TRUE, only.values = TRUE)$values))
    print(rows[[length(rows)]], row.names = FALSE)
  }
}

out <- do.call(rbind, rows)
file <- file.path("results", sprintf("block_size_sweep_%s_K%d_seed%d.csv",
                                     engine, K, seed))
utils::write.csv(out, file, row.names = FALSE)
cat("\n== block-size sweep ==\n")
print(out, row.names = FALSE)
cat("\nwritten: ", file, "\n", sep = "")
