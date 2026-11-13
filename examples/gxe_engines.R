# Why is fit_G_env() conservative? Three estimates of the character-state G
# on the simulations of gwas_gxe_calibration.R (sex as the environment), each
# then used for the environment-specific genotype-by-sex scan:
#   sommer_blocks  fit_G_env() as it is: sommer blocks, cross-sex residual
#                  covariances estimated freely, blocks assembled
#   wombat_blocks  the same blocks fitted by WOMBAT, which holds the
#                  unidentified cross-sex residual covariances at zero
#   wombat_full    one WOMBAT fit of all 2P character states: no subsets, no
#                  assembly or projection
# Scenarios: "none" (cross-sex genetic correlation 1) and "polygenic" (0.5).
#
# Run from the repository root, with BLAS and OpenMP limited to one thread:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript examples/gxe_engines.R
#
# Needs sommer, stanimal and the WOMBAT program (options(indri.wombat =) or
# the PATH). Settings:
#   INDRI_GXE_P        traits                                  (default 10)
#   INDRI_GXE_SEEDS    replicates, as "from:to"                (default 1:20)
#   INDRI_GXE_N_FAM    families of 5 full sibs                 (default 100)
#   INDRI_GXE_M        null SNPs per replicate                 (default 3000)
#   INDRI_GXE_FITS     estimates to run, comma separated       (default all)
#   INDRI_GXE_CORES    replicates run at once (forked)         (default 1)
#
# Writes results/gxe_engines_P<P>.csv and its _summary.csv.

Sys.setenv(OMP_NUM_THREADS = "1")
if (!identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "1"))
  warning("OPENBLAS_NUM_THREADS is not 1: with a multithreaded BLAS every ",
          "replicate starts many threads. Set it before starting R.",
          call. = FALSE, immediate. = TRUE)
library(indri)
source("examples/gwas_simulate.R")

P <- as.integer(Sys.getenv("INDRI_GXE_P", "10"))
seeds <- eval(parse(text = Sys.getenv("INDRI_GXE_SEEDS", "1:20")))
n_fam <- as.integer(Sys.getenv("INDRI_GXE_N_FAM", "100"))
M <- as.integer(Sys.getenv("INDRI_GXE_M", "3000"))
fits <- strsplit(Sys.getenv("INDRI_GXE_FITS",
                            "sommer_blocks,wombat_blocks,wombat_full"),
                 ",")[[1]]
cores <- as.integer(Sys.getenv("INDRI_GXE_CORES", "1"))
alphas <- c(1e-2, 1e-3, 1e-4)
scenarios <- c(none = 1, polygenic = 0.5)   # true cross-sex r_G

# Background phenotypes with cross-sex genetic correlation r (sex is
# X[, "sex"]): a shared part and sex-specific parts, as in
# gwas_gxe_calibration.R.
background <- function(pop, r) {
  N <- pop$N
  sex <- pop$X[, "sex"]
  L_A <- t(chol(pop$A))
  draw_U <- function() L_A %*% matrix(stats::rnorm(N * P), N, P) %*% chol(pop$G)
  U <- sqrt(r) * draw_U() +
    sqrt(1 - r) * (ifelse(sex == 1, 1, 0) * draw_U() +
                   ifelse(sex == 0, 1, 0) * draw_U())
  noise <- matrix(stats::rnorm(N * P), N, P) %*% chol(pop$E)
  Y <- 10 + outer(sex, rep(0.5, P)) + U + noise
  dimnames(Y) <- dimnames(pop$Y)
  Y
}

fit_env <- function(how, Y, pop, seed) {
  args <- list(Y, env = "sex", X = pop$X, A = pop$A, cores = 1L, seed = seed,
               progress = FALSE)
  switch(how,
    sommer_blocks = do.call(fit_G_env, c(args, engine = "sommer")),
    wombat_blocks = do.call(fit_G_env, c(args, engine = "wombat")),
    wombat_full = do.call(fit_G_env, c(args, engine = "wombat",
                                       list(subsets = matrix(seq_len(2L * P), 1L),
                                            assembly = "average"))))
}

one_replicate <- function(seed) {
  pop <- gwas_simulate_population(P = P, n_fam = n_fam, M = M, seed = seed,
                                  G_method = "lkj")
  rows <- list()
  for (scenario in names(scenarios)) {
    r <- scenarios[[scenario]]
    set.seed(seed + 1000L)
    Y <- background(pop, r)
    # true character-state G, states ordered trait:0 ..., trait:1 ...
    truth <- kronecker(matrix(c(1, r, r, 1), 2), pop$G)
    for (how in fits) {
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(fit_env(how, Y, pop, seed), error = function(e) e)
      elapsed <- proc.time()[["elapsed"]] - started
      if (inherits(fit, "error")) {
        message("seed ", seed, " ", scenario, " ", how, ": ",
                conditionMessage(fit))
        next
      }
      levels <- fit$levels
      order <- as.vector(outer(colnames(Y), levels, paste, sep = ":"))
      G <- fit$G[order, order]
      between <- G[seq_len(P), P + seq_len(P)]
      within <- G[seq_len(P), seq_len(P)]
      scan <- indri_gwas(Y, pop$genotypes, fit = fit, A = pop$A, X = pop$X,
                         gxe = "sex")
      for (which in c("gxe", "joint")) {
        p <- c(if (which == "gxe") scan$pvalues_gxe else scan$pvalues)
        p <- p[!is.na(p)]
        df <- P * (if (which == "gxe") 1L else 2L)
        lambda <- stats::qchisq(stats::median(p), df, lower.tail = FALSE) /
          stats::qchisq(0.5, df, lower.tail = FALSE)
        rate <- vapply(alphas, function(a) mean(p <= a), numeric(1))
        rows[[length(rows) + 1L]] <- data.frame(
          seed = seed, P = P, N = pop$N, scenario = scenario, r_true = r,
          fit = how, which = which, lambda_gc = lambda, t1e_1e2 = rate[1],
          t1e_1e3 = rate[2], t1e_1e4 = rate[3],
          r_between = mean(fit$corr_between),
          # the within-sex G and the between-sex covariance against the truth,
          # relative Frobenius errors
          err_within = norm(within - truth[seq_len(P), seq_len(P)], "F") /
            norm(truth[seq_len(P), seq_len(P)], "F"),
          err_between = norm(between - truth[seq_len(P), P + seq_len(P)], "F") /
            norm(truth[seq_len(P), P + seq_len(P)], "F"),
          # between-sex covariance relative to the truth (1 = unbiased)
          ratio_between = sum(diag(between)) /
            sum(diag(truth[seq_len(P), P + seq_len(P)])),
          unconverged = fit$fit$diagnostics$n_unconverged %||% NA_integer_,
          blocks = nrow(fit$fit$subsets), seconds = elapsed)
      }
    }
  }
  do.call(rbind, rows)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

message(sprintf("gxe engines: P = %d, N = %d, M = %d, %d replicates on %d cores; %s",
                P, 5L * n_fam, M, length(seeds), cores,
                paste(fits, collapse = ", ")))
results <- parallel::mclapply(seeds, function(s)
  tryCatch(one_replicate(s), error = function(e) {
    message("seed ", s, " failed: ", conditionMessage(e))
    NULL
  }), mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
tag <- Sys.getenv("INDRI_GXE_TAG", "")
file <- sprintf("results/gxe_engines_P%d%s.csv", P, tag)
utils::write.csv(results, file, row.names = FALSE)
cells <- split(results, list(results$scenario, results$fit, results$which),
               drop = TRUE)
summary <- do.call(rbind, lapply(cells, function(r)
  data.frame(scenario = r$scenario[1], fit = r$fit[1], which = r$which[1],
             replicates = nrow(r), r_true = r$r_true[1],
             r_between = mean(r$r_between),
             ratio_between = mean(r$ratio_between),
             err_within = mean(r$err_within), err_between = mean(r$err_between),
             lambda_gc = mean(r$lambda_gc), t1e_1e2 = mean(r$t1e_1e2),
             t1e_1e3 = mean(r$t1e_1e3), t1e_1e4 = mean(r$t1e_1e4),
             unconverged = mean(r$unconverged), seconds = mean(r$seconds))))
summary <- summary[order(summary$which, summary$scenario, summary$fit), ]
utils::write.csv(summary, sub("\\.csv$", "_summary.csv", file),
                 row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
cat("\nwrote", file, "\n")
