# Does shrinking the assembled genetic correlations help? Schafer and
# Strimmer's (2005) analytic shrinkage, applied after fit_G()'s projection,
# against three targets, scored against the simulated truth.
#
# Run from the repository root, with BLAS and OpenMP limited to one thread:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript examples/shrinkage_G.R
#
# Needs sommer and stanimal. Settings:
#   INDRI_SHRINK_K       traits                             (default 35)
#   INDRI_SHRINK_N_FAM   families of 5 full sibs            (default 60)
#   INDRI_SHRINK_SEEDS   replicates, as "from:to"           (default 1:20)
#   INDRI_SHRINK_G       G generator, "partial" or "lkj"    (default partial)
#   INDRI_SHRINK_CORES   replicates run at once (forked)    (default 1)
#
# For each replicate the assembled corrG r (off-diagonal entries) is shrunk
# towards a target t, r_s = lambda t + (1 - lambda) r, with the Ledoit-Wolf
# intensity lambda = sum(var(r)) / sum((r - t)^2) clamped to [0, 1], var(r)
# from the mixture draws (corrG_draws, one variance per entry). Targets:
#   zero        no genetic correlation (Schafer and Strimmer's target D)
#   phenotypic  the phenotypic correlations (residuals of Y on X): G
#               proportional to P, Cheverud's conjecture
#   constant    the mean genetic correlation (Ledoit and Wolf 2004)
# Each target is positive definite, so r_s is a correlation matrix. Scores
# against the true genetic correlations: median absolute error, error on the
# strongest 10% (by |truth|), slope of estimate on truth (1 = no shrinkage),
# relative Frobenius error, and integration (eigenvalue variance of the
# correlation matrix over the truth's; STATUS.md: indri understates it).
#
# Writes results/shrinkage_G_K<K>[_lkj].csv and its _summary.csv.

Sys.setenv(OMP_NUM_THREADS = "1")
library(indri)
source("examples/gwas_simulate.R")

K <- as.integer(Sys.getenv("INDRI_SHRINK_K", "35"))
n_fam <- as.integer(Sys.getenv("INDRI_SHRINK_N_FAM", "60"))
seeds <- eval(parse(text = Sys.getenv("INDRI_SHRINK_SEEDS", "1:20")))
G_method <- Sys.getenv("INDRI_SHRINK_G", "partial")
cores <- as.integer(Sys.getenv("INDRI_SHRINK_CORES", "1"))

upper <- function(m) m[upper.tri(m)]
eigen_var <- function(r) stats::var(eigen(r, symmetric = TRUE,
                                          only.values = TRUE)$values)

shrink <- function(r, v, target) {
  d <- upper(r) - upper(target)
  lambda <- min(1, max(0, sum(upper(v)) / sum(d^2)))
  out <- lambda * target + (1 - lambda) * r
  diag(out) <- 1
  list(r = out, lambda = lambda)
}

score <- function(est, truth) {
  e <- upper(est); t <- upper(truth)
  strong <- abs(t) >= stats::quantile(abs(t), 0.9)
  c(median_err = stats::median(abs(e - t)),
    strong_err = mean(abs(e - t)[strong]),
    slope = unname(stats::coef(stats::lm(e ~ t))[2]),
    rel_frobenius = sqrt(sum((e - t)^2) / sum(t^2)),
    integration = eigen_var(est) / eigen_var(truth),
    min_eigen = min(eigen(est, symmetric = TRUE, only.values = TRUE)$values))
}

one_replicate <- function(seed) {
  sim <- gwas_simulate_population(P = K, n_fam = n_fam, M = 1L, seed = seed,
                                  G_method = G_method)
  fit <- fit_G(sim$Y, X = sim$X, A = sim$A, engine = "sommer", cores = 1L,
               n_draws = 500L, seed = seed, progress = FALSE)
  r <- fit$corrG
  v <- apply(fit$corrG_draws, c(1L, 2L), stats::var)
  truth <- stats::cov2cor(sim$G)
  resid <- stats::lm.fit(sim$X, sim$Y)$residuals
  targets <- list(
    zero = diag(K),
    phenotypic = stats::cor(resid),
    constant = {
      m <- matrix(mean(upper(r)), K, K)
      diag(m) <- 1
      m
    })
  rows <- list(data.frame(seed = seed, K = K, method = "none", lambda = 0,
                          t(score(r, truth))))
  for (name in names(targets)) {
    s <- shrink(r, v, targets[[name]])
    rows[[length(rows) + 1L]] <- data.frame(seed = seed, K = K, method = name,
                                            lambda = s$lambda,
                                            t(score(s$r, truth)))
  }
  # the target itself, for reference
  rows[[length(rows) + 1L]] <- data.frame(seed = seed, K = K,
                                          method = "phenotypic_only",
                                          lambda = 1,
                                          t(score(targets$phenotypic, truth)))
  do.call(rbind, rows)
}

message(sprintf("shrinkage of G: K = %d, N = %d, G %s, %d replicates on %d cores",
                K, 5L * n_fam, G_method, length(seeds), cores))
results <- parallel::mclapply(seeds, function(s)
  tryCatch(one_replicate(s), error = function(e) {
    message("seed ", s, " failed: ", conditionMessage(e))
    NULL
  }), mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
tag <- if (G_method == "partial") "" else paste0("_", G_method)
file <- sprintf("results/shrinkage_G_K%d%s.csv", K, tag)
utils::write.csv(results, file, row.names = FALSE)
metrics <- c("lambda", "median_err", "strong_err", "slope", "rel_frobenius",
             "integration", "min_eigen")
summary <- stats::aggregate(results[metrics], list(method = results$method),
                            mean)
summary <- summary[match(c("none", "zero", "phenotypic", "constant",
                           "phenotypic_only"), summary$method), ]
utils::write.csv(summary, sub("\\.csv$", "_summary.csv", file),
                 row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
cat("\nwrote", file, "\n")
