# Power of indri_gwas()'s tests against the pattern of a SNP's effects
# across traits: the multi-trait test against single-trait tests adjusted
# for the number of traits (LiMMBo's comparison, here with indri's G and E).
#
# Run from the repository root, with BLAS and OpenMP limited to one thread
# (they must be set before R starts; each forked replicate is meant to use
# one core, and a multithreaded BLAS otherwise starts hundreds of threads
# per process):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript examples/gwas_power.R
#
# Needs sommer and stanimal (for random_covariance()). Settings:
#   INDRI_GWAS_STUDY_P       traits                           (default 10)
#   INDRI_GWAS_STUDY_SEEDS   populations, as "from:to"         (default 1:20)
#   INDRI_GWAS_STUDY_N_FAM   families of 5 full sibs           (default 100)
#   INDRI_GWAS_STUDY_REPS    causal SNPs per population, per
#                            scenario and effect size          (default 10)
#   INDRI_GWAS_STUDY_H       effect sizes, comma-separated     (default
#                            0.001,0.0025,0.005)
#   INDRI_GWAS_STUDY_CORES   populations run at once (forked)  (default 1)
#
# For each population, G and E are fitted once on the background phenotypes
# (no causal SNP). Each replicate then adds one causal SNP to them, with
# effects beta = s * v on the traits: v a unit vector giving the pattern, s
# set so that the SNP explains a fraction h of the total phenotypic variance
# (h * trace(G + E)). Patterns:
#   one        a single trait
#   three      three traits, same sign
#   all_same   every trait, same sign
#   all_random every trait, random signs
#   against    the direction of least phenotypic variance (the smallest
#              eigenvector of G + E), where a joint test gains most
#   dominance  a heterozygote effect on three traits, no additive effect
#
# Writes results/gwas_power_P<P>.csv (one row per replicate and test) and
# results/gwas_power_P<P>_summary.csv (power at alpha = 1e-5 and 5e-8).

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
reps <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_REPS", "10"))
hs <- as.numeric(strsplit(Sys.getenv("INDRI_GWAS_STUDY_H", "0.001,0.0025,0.005"),
                          ",")[[1]])
cores <- as.integer(Sys.getenv("INDRI_GWAS_STUDY_CORES", "1"))
alphas <- c(1e-5, 5e-8)
scenarios <- c("one", "three", "all_same", "all_random", "against",
               "dominance")

pattern <- function(scenario, Ptot) {
  v <- switch(scenario,
    one = c(1, rep(0, P - 1L)),
    three = c(rep(1, 3L), rep(0, P - 3L)),
    all_same = rep(1, P),
    all_random = stats::rnorm(P),
    against = eigen(Ptot, symmetric = TRUE)$vectors[, P],
    dominance = c(rep(1, 3L), rep(0, P - 3L)))
  v / sqrt(sum(v^2))
}

one_population <- function(seed) {
  pop <- gwas_simulate_population(P = P, n_fam = n_fam, M = 1L, seed = seed)
  fit <- fit_G(pop$Y, X = pop$X, A = pop$A, engine = "sommer", cores = 1L,
               seed = seed, progress = FALSE)
  Ptot <- pop$G + pop$E
  rows <- list()
  for (scenario in scenarios) for (h in hs) for (r in seq_len(reps)) {
    maf <- stats::runif(1L, 0.2, 0.5)
    g <- gwas_simulate_genotypes(pop$n_fam, pop$size, maf)
    dimnames(g) <- list(rownames(pop$Y), "causal")
    code <- if (scenario == "dominance") as.numeric(g == 1) else c(g)
    v <- pattern(scenario, Ptot)
    s <- sqrt(h * sum(diag(Ptot)) / stats::var(code))
    Y <- pop$Y + outer(code, s * v)

    multi <- indri_gwas(Y, g, fit = fit, A = pop$A, X = pop$X)
    single <- indri_gwas(Y, g, A = pop$A, X = pop$X, mode = "single_trait",
                         adjust = "effective")
    p <- c(multi = unname(multi$pvalues),
           single_effective = min(single$pvalues_adjust),
           single_bonferroni = min(1, min(single$pvalues) * P))
    if (scenario == "dominance") {
      dom <- indri_gwas(Y, g, fit = fit, A = pop$A, X = pop$X,
                        dominance = TRUE)
      p <- c(p, multi_dominance_joint = unname(dom$pvalues),
             multi_dominance = unname(dom$pvalues_dominance))
    }
    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed, P = P, N = pop$N, scenario = scenario, h = h, rep = r,
      maf = maf, test = names(p), p = unname(p), delta = multi$delta,
      n_effective = single$n_tests)
  }
  do.call(rbind, rows)
}

message(sprintf("gwas power: P = %d, N = %d, %d populations x %d scenarios x %d sizes x %d reps on %d cores",
                P, 5L * n_fam, length(seeds), length(scenarios), length(hs),
                reps, cores))
results <- parallel::mclapply(seeds, function(s)
  tryCatch(one_population(s), error = function(e) {
    message("seed ", s, " failed: ", conditionMessage(e))
    NULL
  }), mc.cores = cores, mc.preschedule = FALSE)
results <- do.call(rbind, results)

dir.create("results", showWarnings = FALSE)
file <- sprintf("results/gwas_power_P%d.csv", P)
utils::write.csv(results, file, row.names = FALSE)

cells <- split(results, list(results$scenario, results$h, results$test),
               drop = TRUE)
summary <- do.call(rbind, lapply(cells, function(r)
  data.frame(scenario = r$scenario[1], h = r$h[1], test = r$test[1],
             replicates = nrow(r),
             power_1e5 = mean(r$p <= alphas[1]),
             power_5e8 = mean(r$p <= alphas[2]),
             median_log10p = stats::median(-log10(r$p)))))
summary <- summary[order(summary$scenario, summary$h, summary$test), ]
summary_file <- sprintf("results/gwas_power_P%d_summary.csv", P)
utils::write.csv(summary, summary_file, row.names = FALSE)
print(summary, row.names = FALSE, digits = 3)
cat("\nwrote", file, "and", summary_file, "\n")
