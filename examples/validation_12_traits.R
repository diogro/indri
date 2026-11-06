# M1 acceptance study: is the assembled 12-trait G as good as fitting the whole
# model at once?
#
# For each seed this simulates a known 12-trait genetic covariance, estimates it
# two ways (indri on trait subsets, and stanimal on all twelve traits at once),
# and compares both with the truth on the correlation scale. The gate from
# NOTES.md is that indri's median per-entry error stays within 2x the full
# model's.
#
# Run from the repository root:
#   Rscript examples/validation_12_traits.R
#
# Takes roughly half an hour. Requires stanimal and CmdStan.

library(indri)

if (!requireNamespace("stanimal", quietly = TRUE))
  stop("This study needs stanimal")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("This study needs cmdstanr and a CmdStan installation")

## Settings ----------------------------------------------------------------

K <- 12L
n_fam <- 60L
sibship_size <- 5L
s <- 4L
min_cooccurrence <- 3L
seeds <- 1:3
block_control <- list(chains = 2L, iter = 800L, init = 0.5,
                      refresh = 0, show_messages = FALSE)
# A budget, not a worker count: with two chains per block, six cores means three
# blocks at a time.
cores <- max(1L, min(6L, parallel::detectCores() - 1L))
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

summary_rows <- list()
entry_rows <- list()

for (seed in seeds) {
  message("\n== seed ", seed, " ==")
  data <- simulate(seed)
  pairs <- upper.tri(data$corr_true)

  started <- Sys.time()
  fit <- fit_G(data$Y, A = data$A, s = s,
               min_cooccurrence = min_cooccurrence, engine = "stanimal",
               engine_args = block_control, cores = cores, seed = 1000L + seed)
  time_indri <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  started <- Sys.time()
  full <- stanimal::lmm_animal(data$Y, matrix(1, data$N, 1L), data$A,
                               chains = 2L, iter = 800L, cores = 2L,
                               seed = 2000L + seed, refresh = 0,
                               show_messages = FALSE)
  time_full <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  corr_full <- stanimal::draws_mean(full, "corrG")

  for (method in c("indri", "full")) {
    estimate <- if (method == "indri") fit$corrG else corr_full
    m <- metrics(estimate, data$corr_true)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      seed = seed, method = method, m,
      time_min = if (method == "indri") time_indri else time_full,
      subsets = if (method == "indri") nrow(fit$subsets) else NA_integer_,
      coverage_min = if (method == "indri") min(fit$coverage[pairs]) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  entry_rows[[length(entry_rows) + 1L]] <- data.frame(
    seed = seed,
    pair = paste(colnames(data$corr_true)[row(data$corr_true)[pairs]],
                 colnames(data$corr_true)[col(data$corr_true)[pairs]], sep = "-"),
    truth = data$corr_true[pairs],
    indri = fit$corrG[pairs],
    full = corr_full[pairs],
    coverage = fit$coverage[pairs],
    stringsAsFactors = FALSE
  )

  rm(full); invisible(gc(FALSE))
}

# One seed with inverse-variance weights, to see whether using the block
# uncertainty helps at this sample size.
message("\n== weights comparison (seed 1) ==")
data <- simulate(1L)
fit_weighted <- fit_G(data$Y, A = data$A, s = s,
                      min_cooccurrence = min_cooccurrence, engine = "stanimal",
                      engine_args = block_control,
                      weights = "inverse_variance", cores = cores, seed = 1001L)
m <- metrics(fit_weighted$corrG, data$corr_true)
summary_rows[[length(summary_rows) + 1L]] <- data.frame(
  seed = 1L, method = "indri_weighted", m, time_min = NA_real_, subsets = NA_integer_,
  coverage_min = NA_integer_, stringsAsFactors = FALSE
)

summary_table <- do.call(rbind, summary_rows)
entries <- do.call(rbind, entry_rows)
utils::write.csv(summary_table, "results/validation_12_traits_summary.csv",
                 row.names = FALSE)
utils::write.csv(entries, "results/validation_12_traits_entries.csv",
                 row.names = FALSE)

cat("\n== summary ==\n")
print(summary_table, row.names = FALSE)

indri_median <- stats::median(summary_table$median_abs[summary_table$method == "indri"])
full_median <- stats::median(summary_table$median_abs[summary_table$method == "full"])
cat("\nmedian over seeds of the median per-entry |corrG error|\n")
cat("  indri: ", round(indri_median, 4), "\n", sep = "")
cat("  full:  ", round(full_median, 4), "\n", sep = "")
cat("  ratio: ", round(indri_median / full_median, 2),
    " (gate: at most 2)\n", sep = "")

cat("\nper-entry error against coverage:\n")
coverage_bands <- cut(entries$coverage, breaks = c(0, 3, 6, 9, Inf),
                      labels = c("1-3", "4-6", "7-9", "10+"))
print(tapply(abs(entries$indri - entries$truth), coverage_bands, stats::median))

cat("\nwritten: results/validation_12_traits_summary.csv\n")
cat("written: results/validation_12_traits_entries.csv\n")
