# Calibration of indri_gwas() under the null: genomic inflation and type-I
# error of every test, on simulated sibs with no causal SNPs.
#
# Run from the repository root, with BLAS and OpenMP limited to one thread
# (they must be set before R starts; each forked replicate is meant to use
# one core, and a multithreaded BLAS otherwise starts hundreds of threads
# per process):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript examples/gwas_calibration.R
#
# Needs sommer and stanimal (for random_covariance()). Settings:
#   INDRI_GWAS_STUDY_P       traits                          (default 10)
#   INDRI_GWAS_STUDY_SEEDS   replicates, as "from:to"         (default 1:20)
#   INDRI_GWAS_STUDY_N_FAM   families of 5 full sibs          (default 100)
#   INDRI_GWAS_STUDY_M       null SNPs per replicate          (default 5000)
#   INDRI_GWAS_STUDY_CORES   replicates run at once (forked)  (default 1)
#
# Writes results/gwas_calibration_P<P>.csv (one row per replicate and test)
# and results/gwas_calibration_P<P>_summary.csv.

Sys.setenv(OMP_NUM_THREADS = "1")
if (!identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "1"))
  warning("OPENBLAS_NUM_THREADS is not 1: with a multithreaded BLAS every ",
          "replicate starts many threads. Set it before starting R.",
          call. = FALSE, immediate. = TRUE)
library(indri)
source("examples/gwas_simulate.R")

P <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_P", "10"))
seeds <- eval(parse(text = Sys.getenv("INDRI_GWAS_STUDY_SEEDS", "1:20")))
n_fam <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_N_FAM", "100"))
M <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_M", "5000"))
cores <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_CORES", "1"))
alphas <- c(1e-2, 1e-3, 1e-4)

one_replicate <- function(seed) {
  sim <- gwas_simulate_population(P = P, n_fam = n_fam, M = M, seed = seed)
  Y <- sim$Y
  started <- Sys.time()
  fit <- fit_G(Y, X = sim$X, A = sim$A, engine = "sommer", cores = 1L,
               seed = seed, progress = FALSE)
  fit_secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  # Candidate E: the residual covariance of Y after X minus the fitted G,
  # eigenvalues floored at 1e-3 of the largest (the assembled E can be
  # nearly singular at many traits; see STATUS.md).
  resid <- stats::lm.fit(sim$X, Y)$residuals
  E_res <- crossprod(resid) / (nrow(Y) - ncol(sim$X)) - fit$G
  e <- eigen((E_res + t(E_res)) / 2, symmetric = TRUE)
  E_res <- e$vectors %*% (pmax(e$values, 1e-3 * e$values[1]) * t(e$vectors))
  dimnames(E_res) <- dimnames(fit$G)
  scans <- list(
    multi_lmm_fit = function()
      indri_gwas(Y, sim$genotypes, fit = fit, A = sim$A, X = sim$X),
    multi_lmm_fit_residual_E = function()
      indri_gwas(Y, sim$genotypes, G = fit$G, E = E_res, A = sim$A, X = sim$X),
    multi_lmm_truth = function()
      indri_gwas(Y, sim$genotypes, G = sim$G, E = sim$E, A = sim$A, X = sim$X),
    multi_lm = function()
      indri_gwas(Y, sim$genotypes, X = sim$X, model = "lm"),
    single_lmm = function()
      indri_gwas(Y, sim$genotypes, A = sim$A, X = sim$X,
                 mode = "single_trait"),
    single_lm = function()
      indri_gwas(Y, sim$genotypes, X = sim$X, mode = "single_trait",
                 model = "lm")
  )
  rows <- lapply(names(scans), function(test) {
    started <- Sys.time()
    scan <- scans[[test]]()
    secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    p <- c(scan$pvalues)
    p <- p[!is.na(p)]
    df_test <- if (scan$mode == "multi_trait") P else 1L
    lambda <- stats::qchisq(stats::median(p), df_test, lower.tail = FALSE) /
      stats::qchisq(0.5, df_test, lower.tail = FALSE)
    rate <- vapply(alphas, function(a) mean(p <= a), numeric(1))
    data.frame(seed = seed, P = P, N = nrow(Y), M = M, test = test,
               lambda_gc = lambda,
               t1e_1e2 = rate[1], t1e_1e3 = rate[2], t1e_1e4 = rate[3],
               n_tests = length(p),
               delta = if (length(scan$delta) == 1L) scan$delta else NA_real_,
               delta_median = if (length(scan$delta) > 1L)
                 stats::median(scan$delta) else NA_real_,
               scan_secs = secs, fit_secs = fit_secs)
  })
  do.call(rbind, rows)
}

message(sprintf("gwas calibration: P = %d, N = %d, M = %d, %d replicates on %d cores",
                P, 5L * n_fam, M, length(seeds), cores))
results <- parallel::mclapply(seeds, function(s)
  tryCatch(one_replicate(s), error = function(e) {
    message("seed ", s, " failed: ", conditionMessage(e))
    NULL
  }), mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
file <- sprintf("results/gwas_calibration_P%d.csv", P)
utils::write.csv(results, file, row.names = FALSE)

summary <- do.call(rbind, lapply(split(results, results$test), function(r)
  data.frame(test = r$test[1], replicates = nrow(r),
             lambda_gc = mean(r$lambda_gc),
             lambda_gc_range = sprintf("%.3f-%.3f", min(r$lambda_gc),
                                       max(r$lambda_gc)),
             t1e_1e2 = mean(r$t1e_1e2), t1e_1e3 = mean(r$t1e_1e3),
             t1e_1e4 = mean(r$t1e_1e4),
             delta = mean(r$delta), scan_secs = mean(r$scan_secs))))
summary_file <- sprintf("results/gwas_calibration_P%d_summary.csv", P)
utils::write.csv(summary, summary_file, row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
cat("\nwrote", file, "and", summary_file, "\n")
