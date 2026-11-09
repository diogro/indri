# Calibration of the genotype-by-environment tests of indri_gwas(gxe = )
# when the background has GxE: null SNPs, sex as the environment, three
# scenarios, each scanned with the usual background (G and E from fit_G())
# and with an environment-specific one (fit_G_env(), tests "*_env").
#   none          no background GxE
#   residual      residual variance twice as large in males
#   polygenic     genetic values differ between the sexes (cross-sex genetic
#                 correlation 0.5 for every trait)
#
# Run from the repository root, with BLAS and OpenMP limited to one thread:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript examples/gwas_gxe_calibration.R
#
# Needs sommer and stanimal. Settings:
#   INDRI_GWAS_STUDY_P       traits                           (default 10)
#   INDRI_GWAS_STUDY_SEEDS   replicates, as "from:to"          (default 1:20)
#   INDRI_GWAS_STUDY_N_FAM   families of 5 full sibs           (default 100)
#   INDRI_GWAS_STUDY_M       null SNPs per replicate           (default 3000)
#   INDRI_GWAS_STUDY_G       G generator, "partial" or "lkj"   (default lkj)
#   INDRI_GWAS_STUDY_CORES   replicates run at once (forked)   (default 1)
#
# Writes results/gwas_gxe_calibration_P<P>.csv and its _summary.csv.

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
M <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_M", "3000"))
G_method <- Sys.getenv("INDRI_GWAS_STUDY_G", "lkj")
cores <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_CORES", "1"))
alphas <- c(1e-2, 1e-3, 1e-4)
scenarios <- c("none", "residual", "polygenic")

# Background phenotypes with the scenario's GxE, on the population's
# design: sex is X[, "sex"].
background <- function(pop, scenario) {
  N <- pop$N
  sex <- pop$X[, "sex"]
  L_A <- t(chol(pop$A))
  draw_U <- function() L_A %*% matrix(stats::rnorm(N * P), N, P) %*% chol(pop$G)
  U <- if (scenario == "polygenic") {
    # shared and sex-specific parts with cross-sex genetic correlation 0.5
    shared <- draw_U()
    sqrt(0.5) * shared + sqrt(0.5) * ifelse(sex == 1, 1, 0) * draw_U() +
      sqrt(0.5) * ifelse(sex == 0, 1, 0) * draw_U()
  } else {
    draw_U()
  }
  noise <- matrix(stats::rnorm(N * P), N, P) %*% chol(pop$E)
  if (scenario == "residual") noise <- noise * ifelse(sex == 1, sqrt(2), 1)
  Y <- 10 + outer(sex, rep(0.5, P)) + U + noise
  dimnames(Y) <- dimnames(pop$Y)
  Y
}

one_replicate <- function(seed) {
  pop <- gwas_simulate_population(P = P, n_fam = n_fam, M = M, seed = seed,
                                  G_method = G_method)
  rows <- list()
  for (scenario in scenarios) {
    Y <- background(pop, scenario)
    fit <- fit_G(Y, X = pop$X, A = pop$A, engine = "sommer", cores = 1L,
                 seed = seed, progress = FALSE)
    fit_env <- fit_G_env(Y, env = "sex", X = pop$X, A = pop$A, cores = 1L,
                         seed = seed, progress = FALSE)
    scans <- list(
      multi_lmm = indri_gwas(Y, pop$genotypes, fit = fit, A = pop$A,
                             X = pop$X, gxe = "sex"),
      single_lmm = indri_gwas(Y, pop$genotypes, A = pop$A, X = pop$X,
                              mode = "single_trait", gxe = "sex"),
      multi_lmm_env = indri_gwas(Y, pop$genotypes, fit = fit_env, A = pop$A,
                                 X = pop$X, gxe = "sex"),
      single_lmm_env = indri_gwas(Y, pop$genotypes, fit = fit_env,
                                  A = pop$A, X = pop$X,
                                  mode = "single_trait", gxe = "sex"))
    for (test in names(scans)) for (which in c("gxe", "joint")) {
      s <- scans[[test]]
      p <- c(if (which == "gxe") s$pvalues_gxe else s$pvalues)
      p <- p[!is.na(p)]
      df <- (if (s$mode == "multi_trait") P else 1L) *
        (if (which == "gxe") 1L else 2L)
      lambda <- stats::qchisq(stats::median(p), df, lower.tail = FALSE) /
        stats::qchisq(0.5, df, lower.tail = FALSE)
      rate <- vapply(alphas, function(a) mean(p <= a), numeric(1))
      rows[[length(rows) + 1L]] <- data.frame(
        seed = seed, P = P, N = pop$N, scenario = scenario, test = test,
        which = which, lambda_gc = lambda, t1e_1e2 = rate[1],
        t1e_1e3 = rate[2], t1e_1e4 = rate[3], n_tests = length(p),
        r_between = mean(fit_env$corr_between))
    }
  }
  do.call(rbind, rows)
}

message(sprintf("gxe calibration: P = %d, N = %d, M = %d, %d replicates on %d cores",
                P, 5L * n_fam, M, length(seeds), cores))
results <- parallel::mclapply(seeds, function(s)
  tryCatch(one_replicate(s), error = function(e) {
    message("seed ", s, " failed: ", conditionMessage(e))
    NULL
  }), mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
file <- sprintf("results/gwas_gxe_calibration_P%d.csv", P)
utils::write.csv(results, file, row.names = FALSE)
cells <- split(results, list(results$scenario, results$test, results$which),
               drop = TRUE)
summary <- do.call(rbind, lapply(cells, function(r)
  data.frame(scenario = r$scenario[1], test = r$test[1], which = r$which[1],
             replicates = nrow(r), lambda_gc = mean(r$lambda_gc),
             t1e_1e2 = mean(r$t1e_1e2), t1e_1e3 = mean(r$t1e_1e3),
             t1e_1e4 = mean(r$t1e_1e4), r_between = mean(r$r_between))))
summary <- summary[order(summary$which, summary$test, summary$scenario), ]
utils::write.csv(summary, sub("\\.csv$", "_summary.csv", file),
                 row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
cat("\nwrote", file, "\n")
