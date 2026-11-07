# Cross-engine comparison on simulated data.
#
# Simulates a known K-trait genetic covariance and estimates it with indri on
# two block engines:
#
#   1. indri with stanimal blocks (posterior mean per block),
#   2. indri with sommer blocks (REML point estimate per block).
#
# Every estimate is compared with the truth on the correlation scale, and the
# engines are compared with each other on the same trait pairs.
#
# Run from the repository root:
#   Rscript examples/compare_engines.R
#
# Needs stanimal, sommer and CmdStan. Environment overrides:
#   INDRI_ENGINE_K       number of traits (10)
#   INDRI_ENGINE_NFAM    families (60)
#   INDRI_ENGINE_SIZE    sibship size (5)
#   INDRI_ENGINE_S       subset size parameter for stanimal blocks (4)
#   INDRI_ENGINE_S_SOMMER subset size parameter for sommer blocks (2)
#   INDRI_ENGINE_MINCOOC minimum pair coverage (3)
#   INDRI_ENGINE_SEED    simulation seed (1)
#   INDRI_ENGINE_CORES   Stan process budget for the stanimal blocks (64)
#   INDRI_ENGINE_CORES_SOMMER number of sommer blocks to run in parallel (16)
#   INDRI_ENGINE_ITER    stanimal iterations per chain (400)
#   INDRI_ENGINE_CHAINS  chains per stanimal block (2)
#   INDRI_ENGINE_SOMMER  1 to run the sommer engine (1)
#   INDRI_ENGINE_OMP_THREADS  OpenMP threads per sommer block (1)

library(indri)

if (!requireNamespace("stanimal", quietly = TRUE))
  stop("This study needs stanimal")
if (!requireNamespace("sommer", quietly = TRUE))
  stop("This study needs sommer")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("This study needs cmdstanr and a CmdStan installation")

## Settings ----------------------------------------------------------------

K <- as.integer(Sys.getenv("INDRI_ENGINE_K", "10"))
n_fam <- as.integer(Sys.getenv("INDRI_ENGINE_NFAM", "60"))
sibship_size <- as.integer(Sys.getenv("INDRI_ENGINE_SIZE", "5"))
s <- as.integer(Sys.getenv("INDRI_ENGINE_S", "4"))
s_sommer <- as.integer(Sys.getenv("INDRI_ENGINE_S_SOMMER", "2"))
min_cooccurrence <- as.integer(Sys.getenv("INDRI_ENGINE_MINCOOC", "3"))
seed <- as.integer(Sys.getenv("INDRI_ENGINE_SEED", "1"))
cores <- as.integer(Sys.getenv("INDRI_ENGINE_CORES", "64"))
cores_sommer <- as.integer(Sys.getenv("INDRI_ENGINE_CORES_SOMMER", "16"))
iter <- as.integer(Sys.getenv("INDRI_ENGINE_ITER", "400"))
chains <- as.integer(Sys.getenv("INDRI_ENGINE_CHAINS", "2"))
run_sommer <- Sys.getenv("INDRI_ENGINE_SOMMER", "1") == "1"

stanimal_control <- list(chains = chains, iter = iter, init = 0.5,
                         refresh = 0, show_messages = FALSE)
# sommer's mmes solver multithreads through OpenMP even for tiny bivariate
# blocks, which oversubscribes the machine when combined with parallel blocks
# and slows each fit down (measured ~2.6x). Pin it to one thread per block and
# let `cores` control the parallelism instead.
Sys.setenv(OMP_NUM_THREADS = Sys.getenv("INDRI_ENGINE_OMP_THREADS", "1"))
dir.create("results", showWarnings = FALSE)

metrics <- function(estimate, truth) {
  pairs <- upper.tri(truth)
  error <- estimate[pairs] - truth[pairs]
  list(
    median_abs = stats::median(abs(error)),
    q90_abs = as.numeric(stats::quantile(abs(error), 0.9, names = FALSE)),
    max_abs = max(abs(error)),
    relative_frobenius = norm(estimate - truth, "F") / norm(truth, "F"),
    corr_of_corr = stats::cor(estimate[pairs], truth[pairs])
  )
}

simulate <- function(seed) {
  set.seed(seed)
  truth <- stanimal::random_covariance(
    n_traits = K, method = "partial", min_var = 0.8, max_var = 1.4,
    trait_names = paste0("trait", seq_len(K))
  )
  residual <- stanimal::random_covariance(
    n_traits = K, method = "lkj", lkj_prior = 2, min_var = 0.5, max_var = 0.9
  )
  A <- kronecker(diag(n_fam), matrix(0.5, sibship_size, sibship_size))
  diag(A) <- 1
  N <- nrow(A)
  genetic <- t(chol(A)) %*% matrix(stats::rnorm(N * K), N, K) %*% chol(truth)
  Y <- genetic + matrix(stats::rnorm(N * K), N, K) %*% chol(residual)
  rownames(Y) <- paste0("id", seq_len(N))
  colnames(Y) <- colnames(truth)
  list(Y = Y, A = A, truth = truth, corr_true = stats::cov2cor(truth), N = N)
}

## Simulate ----------------------------------------------------------------

message("== simulating ", K, " traits, N = ", n_fam * sibship_size, " ==")
data <- simulate(seed)
pairs <- upper.tri(data$corr_true)

## indri + sommer -----------------------------------------------------------

if (run_sommer) {
  message("== indri, sommer engine (", cores_sommer, " parallel blocks) ==")
  started <- Sys.time()
  fit_sommer <- fit_G(data$Y, A = data$A, s = s_sommer,
                      min_cooccurrence = min_cooccurrence, engine = "sommer",
                      cores = cores_sommer, seed = 1000L + seed)
  time_sommer <- as.numeric(difftime(Sys.time(), started, units = "mins"))
} else {
  fit_sommer <- NULL
  time_sommer <- NA_real_
}

## indri + stanimal ----------------------------------------------------------

message("== indri, stanimal engine (", cores, " process budget) ==")
started <- Sys.time()
fit_stanimal <- fit_G(data$Y, A = data$A, s = s,
                      min_cooccurrence = min_cooccurrence, engine = "stanimal",
                      engine_args = stanimal_control, cores = cores,
                      seed = 2000L + seed)
time_stanimal <- as.numeric(difftime(Sys.time(), started, units = "mins"))

## Compare -------------------------------------------------------------------

estimates <- list(
  indri_sommer = if (!is.null(fit_sommer)) fit_sommer$corrG else NULL,
  indri_stanimal = fit_stanimal$corrG
)
times <- c(
  indri_sommer = time_sommer,
  indri_stanimal = time_stanimal
)

summary_rows <- list()
block_s <- c(indri_sommer = s_sommer, indri_stanimal = s)
for (method in names(estimates)) {
  estimate <- estimates[[method]]
  if (is.null(estimate)) next
  m <- metrics(estimate, data$corr_true)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    method = method,
    s = block_s[[method]],
    m,
    time_min = times[[method]],
    stringsAsFactors = FALSE
  )
}
summary_table <- do.call(rbind, summary_rows)

cat("\n== accuracy vs truth (correlation scale) ==\n")
print(summary_table, row.names = FALSE)

cat("\n== cross-engine agreement (correlation of per-pair estimates) ==\n")
agree <- matrix(NA_real_, length(estimates), length(estimates),
                dimnames = list(names(estimates), names(estimates)))
for (i in names(estimates)) {
  for (j in names(estimates)) {
    if (is.null(estimates[[i]]) || is.null(estimates[[j]])) next
    agree[i, j] <- stats::cor(estimates[[i]][pairs], estimates[[j]][pairs])
  }
}
print(round(agree, 3))

cat("\n== truth vs estimate, per pair ==\n")
entry <- data.frame(
  pair = paste(colnames(data$corr_true)[row(data$corr_true)[pairs]],
               colnames(data$corr_true)[col(data$corr_true)[pairs]], sep = "-"),
  truth = data$corr_true[pairs]
)
for (method in names(estimates)) {
  if (!is.null(estimates[[method]])) entry[[method]] <- estimates[[method]][pairs]
}
print(entry, row.names = FALSE)

utils::write.csv(summary_table,
                 file.path("results", sprintf("compare_engines_summary_seed%s.csv", seed)),
                 row.names = FALSE)
utils::write.csv(entry,
                 file.path("results", sprintf("compare_engines_entries_seed%s.csv", seed)),
                 row.names = FALSE)
cat("\nwritten: results/compare_engines_*_seed", seed, ".csv\n", sep = "")
