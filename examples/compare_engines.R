# Cross-engine comparison on simulated data.
#
# Simulates a known K-trait genetic covariance and estimates it three ways:
#
#   1. indri with stanimal blocks (posterior mean per block),
#   2. indri with sommer blocks (REML point estimate per block),
#   3. the full K-trait stanimal model as the reference.
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
#   INDRI_ENGINE_S       subset size parameter (4)
#   INDRI_ENGINE_MINCOOC minimum pair coverage (3)
#   INDRI_ENGINE_SEED    simulation seed (1)
#   INDRI_ENGINE_CORES   Stan process budget for the stanimal blocks (64)
#   INDRI_ENGINE_ITER    stanimal iterations per chain (400)
#   INDRI_ENGINE_CHAINS  chains per stanimal block (2)
#   INDRI_ENGINE_FULL    1 to run the full stanimal reference (1)

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
min_cooccurrence <- as.integer(Sys.getenv("INDRI_ENGINE_MINCOOC", "3"))
seed <- as.integer(Sys.getenv("INDRI_ENGINE_SEED", "1"))
cores <- as.integer(Sys.getenv("INDRI_ENGINE_CORES", "64"))
iter <- as.integer(Sys.getenv("INDRI_ENGINE_ITER", "400"))
chains <- as.integer(Sys.getenv("INDRI_ENGINE_CHAINS", "2"))
run_full <- Sys.getenv("INDRI_ENGINE_FULL", "1") == "1"

stanimal_control <- list(chains = chains, iter = iter, init = 0.5,
                         refresh = 0, show_messages = FALSE)
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

message("== indri, sommer engine ==")
started <- Sys.time()
fit_sommer <- fit_G(data$Y, A = data$A, s = s,
                    min_cooccurrence = min_cooccurrence, engine = "sommer",
                    seed = 1000L + seed)
time_sommer <- as.numeric(difftime(Sys.time(), started, units = "mins"))

## indri + stanimal ----------------------------------------------------------

message("== indri, stanimal engine (", cores, " process budget) ==")
started <- Sys.time()
fit_stanimal <- fit_G(data$Y, A = data$A, s = s,
                      min_cooccurrence = min_cooccurrence, engine = "stanimal",
                      engine_args = stanimal_control, cores = cores,
                      seed = 2000L + seed)
time_stanimal <- as.numeric(difftime(Sys.time(), started, units = "mins"))

## full stanimal reference ----------------------------------------------------

if (run_full) {
  message("== full ", K, "-trait stanimal model ==")
  started <- Sys.time()
  full <- stanimal::lmm_animal(data$Y, matrix(1, data$N, 1L), data$A,
                               chains = 2L, iter = 500L, cores = 2L,
                               seed = 3000L + seed, refresh = 0,
                               show_messages = FALSE)
  time_full <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  corr_full <- stanimal::draws_mean(full, "corrG")
} else {
  corr_full <- NULL
  time_full <- NA_real_
}

## Compare -------------------------------------------------------------------

estimates <- list(
  indri_sommer = fit_sommer$corrG,
  indri_stanimal = fit_stanimal$corrG,
  full_stanimal = corr_full
)
times <- c(
  indri_sommer = time_sommer,
  indri_stanimal = time_stanimal,
  full_stanimal = time_full
)

summary_rows <- list()
for (method in names(estimates)) {
  estimate <- estimates[[method]]
  if (is.null(estimate)) next
  m <- metrics(estimate, data$corr_true)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    method = method,
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
