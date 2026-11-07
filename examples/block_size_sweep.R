# stanimal block-size sweep at K = 50 with a "partial" (near-boundary) truth.
#
#   Rscript examples/block_size_sweep.R
#
#   INDRI_SWEEP_K      traits (50)
#   INDRI_SWEEP_CORES  Stan process budget (128)
#   INDRI_SWEEP_ITER   iterations per chain (300)
#   INDRI_SWEEP_CHAINS chains per block (2)
#   INDRI_SWEEP_S      comma-separated s values to sweep (2,3,4,5,6)

library(indri)

K <- as.integer(Sys.getenv("INDRI_SWEEP_K", "50"))
cores <- as.integer(Sys.getenv("INDRI_SWEEP_CORES", "128"))
iter <- as.integer(Sys.getenv("INDRI_SWEEP_ITER", "300"))
chains <- as.integer(Sys.getenv("INDRI_SWEEP_CHAINS", "2"))
s_values <- as.integer(strsplit(Sys.getenv("INDRI_SWEEP_S", "2,3,4,5,6"), ",")[[1]])
seed <- 1L
n_fam <- 60L
size <- 5L
min_cooccurrence <- 3L

set.seed(seed)
truth <- stanimal::random_covariance(
  n_traits = K, method = "partial", min_var = 0.8, max_var = 1.4,
  trait_names = paste0("trait", seq_len(K)))
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

ctrl <- list(chains = chains, iter = iter, init = 0.5,
             refresh = 0, show_messages = FALSE)

rows <- list()
for (s in s_values) {
  message("== s = ", s, " ==")
  started <- Sys.time()
  fit <- fit_G(Y, A = A, s = s, min_cooccurrence = min_cooccurrence,
               engine = "stanimal", engine_args = ctrl, cores = cores,
               seed = 1000L + seed)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  rows[[length(rows) + 1L]] <- data.frame(
    s = s, block_size = ncol(fit$subsets), n_blocks = nrow(fit$subsets),
    time_min = round(elapsed, 2),
    median_abs = round(stats::median(abs(fit$corrG[pairs] - corr_true[pairs])), 4),
    corr_of_corr = round(stats::cor(fit$corrG[pairs], corr_true[pairs]), 4),
    min_eig = round(min(eigen(fit$G, symmetric = TRUE, only.values = TRUE)$values), 4))
  print(rows[[length(rows)]], row.names = FALSE)
}

out <- do.call(rbind, rows)
utils::write.csv(out, file.path("results", "block_size_sweep.csv"), row.names = FALSE)
cat("\n== block-size sweep ==\n")
print(out, row.names = FALSE)
cat("\nwritten: results/block_size_sweep.csv\n")
