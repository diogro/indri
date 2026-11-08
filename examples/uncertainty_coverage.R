# Uncertainty study: does the assembled block uncertainty cover the truth?
#
# Simulates a known K-trait genetic covariance (12 by default), fits it with
# indri on any engine and, optionally, the full stanimal model, and reports the
# fraction of true correlation entries inside each estimate's 95% band (per
# pair-coverage count for indri), the band width, the slope of estimate on
# truth, and the eigenvalue variance of corrG against the truth's.
#
# For stanimal the band comes from the block posteriors, assembled three ways
# from the same fit, on the correlation scale (each block's draws -> correlation
# draws):
#   indri_average   draw d averaged over covering blocks (what G_draws did by
#                   default); shrinks the spread by ~1/sqrt(coverage)
#   indri_mixture   each entry of draw d taken from one random covering block
#   indri_inflated  the average, with deviations from its mean scaled by
#                   sqrt(coverage): what the average would give if the block
#                   errors were perfectly correlated
# None is a joint posterior (subsets share traits and data). The fits are saved
# to results/uncertainty_coverage_fits_seed<N>.rds (git-ignored) so further
# arms can be scored without refitting.
#
# Run from the repository root:
#   Rscript examples/uncertainty_coverage.R
#
# Requires stanimal and CmdStan. Environment overrides:
#   INDRI_COVERAGE_SEEDS   seeds to simulate (1,2,3)
#   INDRI_COVERAGE_CORES   block-parallelism budget for indri (1)
#   INDRI_COVERAGE_ITER    indri iterations per block chain (400)
#   INDRI_COVERAGE_CI      interval level (0.95)
#   INDRI_COVERAGE_ENGINE  stanimal, sommer or asreml (stanimal). REML engines
#                          use their defaults (s = 3, c = 1) and are scored
#                          on fit$corrG_draws, drawn from each pair's REML
#                          sampling covariance; with one block per pair there
#                          is nothing to average or mix.
#   INDRI_COVERAGE_FULL    1 to fit the full stanimal model as a reference
#                          (1 for stanimal, 0 otherwise; same seeds, same data)
#   INDRI_COVERAGE_DRAWS   draws per block for the REML engines (1000)
#   INDRI_COVERAGE_LKJ     stanimal lkj_prior for the blocks (stanimal default 2)
#   INDRI_COVERAGE_S_STAN  stanimal subset size parameter (4)
#   INDRI_COVERAGE_K       traits (12)
#   INDRI_COVERAGE_NFAM    families of 5 full sibs; N = 5 x NFAM (60, N = 300)
#   INDRI_COVERAGE_S, INDRI_COVERAGE_C
#                          REML subset size and distinct contexts per pair
#                          (package defaults: 3 and 3)

library(indri)

if (!requireNamespace("stanimal", quietly = TRUE))
  stop("This study needs stanimal")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("This study needs cmdstanr and a CmdStan installation")

## Settings ----------------------------------------------------------------

K <- as.integer(Sys.getenv("INDRI_COVERAGE_K", "12"))
n_fam <- as.integer(Sys.getenv("INDRI_COVERAGE_NFAM", "60"))   # N = 5 n_fam
sibship_size <- 5L
s <- as.integer(Sys.getenv("INDRI_COVERAGE_S_STAN", "4"))     # stanimal blocks
min_cooccurrence <- 3L
lkj <- Sys.getenv("INDRI_COVERAGE_LKJ", "")                   # stanimal lkj_prior
lkj <- if (nzchar(lkj)) as.numeric(lkj) else NULL
# REML block settings; empty means the package default (s = 3, c = 3)
reml_s <- Sys.getenv("INDRI_COVERAGE_S", "")
reml_c <- Sys.getenv("INDRI_COVERAGE_C", "")
reml_s <- if (nzchar(reml_s)) as.integer(reml_s) else NULL
reml_c <- if (nzchar(reml_c)) as.integer(reml_c) else NULL
seeds <- as.integer(strsplit(Sys.getenv("INDRI_COVERAGE_SEEDS", "1,2,3"), ",")[[1]])
cores <- as.integer(Sys.getenv("INDRI_COVERAGE_CORES", "1"))
iter <- as.integer(Sys.getenv("INDRI_COVERAGE_ITER", "400"))
ci <- as.numeric(Sys.getenv("INDRI_COVERAGE_CI", "0.95"))
engine <- Sys.getenv("INDRI_COVERAGE_ENGINE", "stanimal")
reml <- engine != "stanimal"
run_full <- Sys.getenv("INDRI_COVERAGE_FULL", if (reml) "0" else "1") == "1"
n_draws <- as.integer(Sys.getenv("INDRI_COVERAGE_DRAWS", "1000"))
# File tag: engine, and any setting that differs from the original study
tag <- paste0(
  if (reml) paste0(engine, "_") else "",
  if (K != 12L || n_fam != 60L) sprintf("K%d_N%d_", K, n_fam * sibship_size) else "",
  if (!is.null(reml_s)) sprintf("s%d_", reml_s) else "",
  if (!is.null(reml_c)) sprintf("c%d_", reml_c) else "",
  if (!reml && s != 4L) sprintf("sstan%d_", s) else "",
  if (!reml && !is.null(lkj)) sprintf("lkj%s_", format(lkj)) else "")
Sys.setenv(OMP_NUM_THREADS = "1")

block_control <- list(chains = 1L, iter = iter, init = 0.5,
                      refresh = 0, show_messages = FALSE)
if (!is.null(lkj)) block_control$lkj_prior <- lkj
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

# Scale each entry's deviations from its mean by sqrt(coverage).
inflate_spread <- function(draws, coverage) {
  centre <- apply(draws, c(1L, 2L), mean)
  scale <- sqrt(pmax(coverage, 1))
  out <- draws
  for (d in seq_len(dim(draws)[3L]))
    out[, , d] <- centre + (draws[, , d] - centre) * scale
  out
}

rows <- list()
ev_rows <- list()
# Integration: relative eigenvalue variance (evolqg::CalcEigenVar) of the
# genetic correlation matrix, truth vs estimate. For the indri point estimate
# only: mixed draws are pointwise, not valid matrices, so they give no
# interval for a whole-matrix quantity. The full model's posterior draws are
# valid matrices and give one.
eigenvar <- function(m) {
  if (!requireNamespace("evolqg", quietly = TRUE)) return(NA_real_)
  evolqg::CalcEigenVar((m + t(m)) / 2)
}

for (seed in seeds) {
  message("\n== seed ", seed, " ==")
  data <- simulate(seed)

  if (reml) {
    fit <- fit_G(data$Y, A = data$A, engine = engine, s = reml_s,
                 min_cooccurrence = reml_c, cores = cores,
                 n_draws = n_draws, seed = 1000L + seed)
    arms <- list(fit$corrG_draws)
    names(arms) <- paste0("indri_", engine)
  } else {
    fit <- fit_G(data$Y, A = data$A, s = s,
                 min_cooccurrence = min_cooccurrence, engine = "stanimal",
                 engine_args = block_control, cores = cores,
                 seed = 1000L + seed)
    block_corr <- lapply(fit$blocks, function(b) corr_draws(b$draws))
    averaged <- indri_assemble_draws(block_corr, fit$subsets, K)
    arms <- list(
      indri_average = averaged,
      indri_mixture = indri_assemble_draws(block_corr, fit$subsets, K,
                                           method = "mixture", seed = seed),
      indri_inflated = inflate_spread(averaged, fit$coverage)
    )
  }
  indri_rows <- do.call(rbind, lapply(names(arms), function(arm) {
    band <- quantile_bands(arms[[arm]], ci)
    out <- coverage_rows(band$lo, band$hi, data$corr_true, fit$coverage)
    out$method <- arm
    out
  }))
  indri_rows$seed <- seed
  # shrinkage: slope of the point estimate on the truth (1 = unbiased scale)
  pairs <- upper.tri(data$corr_true)
  indri_rows$slope <- unname(coef(stats::lm(fit$corrG[pairs] ~
                                              data$corr_true[pairs]))[2])

  full_draws <- NULL
  full_rows <- NULL
  if (run_full) {
    full <- stanimal::lmm_animal(data$Y, matrix(1, data$N, 1L), data$A,
                                 chains = 2L, iter = 500L, cores = 2L,
                                 seed = 2000L + seed, refresh = 0,
                                 show_messages = FALSE)
    full_draws <- stanimal::fit_draws(full)$corrG
    full_bands <- quantile_bands(full_draws, ci)
    full_rows <- coverage_rows(full_bands$lo, full_bands$hi, data$corr_true)
    full_rows$seed <- seed
    full_rows$method <- "full"
    full_mean <- apply(full_draws, c(1L, 2L), mean)
    full_rows$slope <- unname(coef(stats::lm(full_mean[pairs] ~
                                               data$corr_true[pairs]))[2])
    rm(full)
  }

  rows[[length(rows) + 1L]] <- rbind(indri_rows, full_rows)
  ev_true <- eigenvar(data$corr_true)
  ev_rows[[length(ev_rows) + 1L]] <- data.frame(
    seed = seed, method = if (reml) paste0("indri_", engine) else "indri",
    ev_true = ev_true, ev_est = eigenvar(fit$corrG),
    ev_lo = NA_real_, ev_hi = NA_real_)
  if (run_full) {
    full_ev <- apply(full_draws, 3L, eigenvar)
    full_ev_mean <- eigenvar(apply(full_draws, c(1L, 2L), mean))
    ev_rows[[length(ev_rows) + 1L]] <- data.frame(
      seed = seed, method = "full", ev_true = ev_true, ev_est = full_ev_mean,
      ev_lo = unname(stats::quantile(full_ev, 0.025)),
      ev_hi = unname(stats::quantile(full_ev, 0.975)))
  }
  saveRDS(list(data = data, subsets = fit$subsets, coverage = fit$coverage,
               G = fit$G, corrG = fit$corrG,
               block_draws = lapply(fit$blocks, `[[`, "draws"),
               block_vcov = lapply(fit$blocks, `[[`, "vcov"),
               diagnostics = fit$diagnostics, full_corr_draws = full_draws),
          sprintf("results/uncertainty_coverage_fits_%sseed%d.rds", tag, seed))
  invisible(gc(FALSE))
}

ev_table <- do.call(rbind, ev_rows)
utils::write.csv(
  ev_table,
  sprintf("results/uncertainty_coverage_eigenvar_%sseed%s.csv", tag,
          paste(unique(ev_table$seed), collapse = "-")),
  row.names = FALSE
)
print(ev_table, row.names = FALSE)

summary <- do.call(rbind, rows)
summary <- summary[, c("seed", "method", "band", "n_pairs", "coverage", "mean_width",
                       "slope")]
utils::write.csv(
  summary,
  sprintf("results/uncertainty_coverage_summary_%sseed%s.csv", tag,
          paste(unique(summary$seed), collapse = "-")),
  row.names = FALSE
)

cat("\n== 95% coverage of the true correlation, indri arms vs full model ==\n")
print(summary, row.names = FALSE)
cat("\nwritten: results/uncertainty_coverage_summary_seed*.csv\n")
