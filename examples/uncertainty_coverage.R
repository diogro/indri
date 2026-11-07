# M4 uncertainty study: does the assembled block uncertainty cover the truth?
#
# Simulates a known 12-trait genetic covariance, fits it two ways — indri on
# trait subsets (stanimal engine) and the full 12-trait model — and reports the
# fraction of true correlation entries inside each estimate's 95% band, per
# pair-coverage count for indri.
#
# The indri band comes from `G_draws`, which is an assembled spread across
# dependent blocks, NOT a joint posterior (subsets share traits). This study is
# the measurement that says how far from calibrated that spread is, compared
# with the full model's own posterior.
#
# Run from the repository root:
#   Rscript examples/uncertainty_coverage.R
#
# Requires stanimal and CmdStan. Environment overrides:
#   INDRI_COVERAGE_SEEDS   seeds to simulate (1,2,3)
#   INDRI_COVERAGE_CORES   block-parallelism budget for indri (1)
#   INDRI_COVERAGE_ITER    indri iterations per block chain (400)
#   INDRI_COVERAGE_CI      interval level (0.95)

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
seeds <- as.integer(strsplit(Sys.getenv("INDRI_COVERAGE_SEEDS", "1,2,3"), ",")[[1]])
cores <- as.integer(Sys.getenv("INDRI_COVERAGE_CORES", "1"))
iter <- as.integer(Sys.getenv("INDRI_COVERAGE_ITER", "400"))
ci <- as.numeric(Sys.getenv("INDRI_COVERAGE_CI", "0.95"))

block_control <- list(chains = 1L, iter = iter, init = 0.5,
                      refresh = 0, show_messages = FALSE)
dir.create("results", showWarnings = FALSE)

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
  list(Y = Y, A = A, corr_true = stats::cov2cor(truth), N = N)
}

# Posterior draws of G -> draws of the correlation matrix, one slice at a time.
corr_draws <- function(G_draws) {
  out <- array(0, dim(G_draws))
  for (d in seq_len(dim(G_draws)[3L]))
    out[, , d] <- stats::cov2cor(G_draws[, , d])
  out
}

# Per-entry (lo, hi) quantiles over the draw axis of a K x K x D array.
quantile_bands <- function(draws, ci) {
  lo <- (1 - ci) / 2
  hi <- 1 - lo
  list(
    lo = apply(draws, c(1L, 2L), stats::quantile, probs = lo, names = FALSE),
    hi = apply(draws, c(1L, 2L), stats::quantile, probs = hi, names = FALSE)
  )
}

# Fraction of true entries inside the band, overall and per pair-coverage band.
coverage_rows <- function(lo, hi, corr_true, coverage = NULL) {
  pairs <- upper.tri(corr_true)
  inside <- corr_true[pairs] >= lo[pairs] & corr_true[pairs] <= hi[pairs]
  width <- hi[pairs] - lo[pairs]

  if (is.null(coverage)) {
    return(data.frame(
      band = "all", n_pairs = sum(pairs),
      coverage = mean(inside), mean_width = mean(width),
      stringsAsFactors = FALSE
    ))
  }

  bands <- cut(coverage[pairs], breaks = c(0, 3, 6, 9, Inf),
               labels = c("1-3", "4-6", "7-9", "10+"))
  data.frame(
    band = levels(bands),
    n_pairs = as.integer(table(bands)),
    coverage = as.numeric(tapply(inside, bands, mean)),
    mean_width = as.numeric(tapply(width, bands, mean)),
    stringsAsFactors = FALSE
  )
}

rows <- list()

for (seed in seeds) {
  message("\n== seed ", seed, " ==")
  data <- simulate(seed)

  fit <- fit_G(data$Y, A = data$A, s = s,
               min_cooccurrence = min_cooccurrence, engine = "stanimal",
               engine_args = block_control, cores = cores, seed = 1000L + seed)
  indri <- quantile_bands(corr_draws(fit$G_draws), ci)
  indri_rows <- coverage_rows(indri$lo, indri$hi, data$corr_true, fit$coverage)
  indri_rows$seed <- seed
  indri_rows$method <- "indri"

  full <- stanimal::lmm_animal(data$Y, matrix(1, data$N, 1L), data$A,
                               chains = 2L, iter = 500L, cores = 2L,
                               seed = 2000L + seed, refresh = 0,
                               show_messages = FALSE)
  full_draws <- stanimal::fit_draws(full)$corrG
  full_bands <- quantile_bands(full_draws, ci)
  full_rows <- coverage_rows(full_bands$lo, full_bands$hi, data$corr_true)
  full_rows$seed <- seed
  full_rows$method <- "full"

  rows[[length(rows) + 1L]] <- rbind(indri_rows, full_rows)
  rm(full); invisible(gc(FALSE))
}

summary <- do.call(rbind, rows)
summary <- summary[, c("seed", "method", "band", "n_pairs", "coverage", "mean_width")]
utils::write.csv(
  summary,
  sprintf("results/uncertainty_coverage_summary_seed%s.csv",
          paste(unique(summary$seed), collapse = "-")),
  row.names = FALSE
)

cat("\n== 95% coverage of the true correlation, indri vs full model ==\n")
print(summary, row.names = FALSE)
cat("\nwritten: results/uncertainty_coverage_summary_seed*.csv\n")
