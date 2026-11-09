# A genome scan with the covariances indri estimates.
#
# Run from the repository root:
#   Rscript examples/gwas.R
#
# Needs sommer (the REML engine) and ggplot2; under a minute. Settings:
#   INDRI_GWAS_CORES  process budget for fit_G()     (default 2)
#   INDRI_GWAS_SEED   seed for the simulation          (default 2026)
#
# The data are simulated: 60 families of 5 full sibs with Mendelian
# genotypes at unlinked SNPs, 10 correlated traits, and four causal SNPs
# chosen to show where each test is strong:
#   rs_big3    large additive effect on 3 of the 10 traits
#   rs_pleio   small additive effects on all 10 traits, of opposite sign in
#              the two groups of genetically correlated traits
#   rs_dom     overdominance on 2 traits (the heterozygote stands out)
#   rs_one     large additive effect on a single trait

Sys.setenv(OMP_NUM_THREADS = "1")   # before sommer loads (see AGENTS.md)
library(indri)

if (!requireNamespace("sommer", quietly = TRUE))
  stop("This example needs sommer: install.packages('sommer')")
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("This example needs ggplot2: install.packages('ggplot2')")

cores <- as.integer(Sys.getenv("INDRI_GWAS_CORES", "2"))
seed <- as.integer(Sys.getenv("INDRI_GWAS_SEED", "2026"))
set.seed(seed)

## The population ----------------------------------------------------------

n_fam <- 60L
size <- 5L
N <- n_fam * size
P <- 10L
n_chrom <- 5L
M_per_chrom <- 400L
M <- n_chrom * M_per_chrom
traits <- paste0("trait", seq_len(P))
ids <- paste0("id", seq_len(N))

# Full sibs share half their alleles on average: the relationship matrix the
# genotypes are drawn under, and the one the analysis uses.
A <- kronecker(diag(n_fam), matrix(0.5, size, size))
diag(A) <- 1
dimnames(A) <- list(ids, ids)

# Parents draw two alleles per SNP from its frequency; each child takes one
# allele from each parent at random. Unlinked SNPs, so no LD between them.
maf <- stats::runif(M, 0.1, 0.5)
causal <- c(rs_big3 = 150L, rs_pleio = 700L, rs_dom = 1300L, rs_one = 1750L)
# Dominance needs all three genotypes in numbers: at frequency 0.5 a quarter
# of the individuals are minor homozygotes.
maf[causal["rs_dom"]] <- 0.5
genotypes <- matrix(0L, N, M)
for (f in seq_len(n_fam)) {
  parents <- array(stats::rbinom(2L * 2L * M, 1L, rep(maf, each = 4L)),
                   c(2L, 2L, M))            # parent x allele x SNP
  kids <- (f - 1L) * size + seq_len(size)
  for (k in kids) {
    pick <- function(parent)
      parents[cbind(parent, stats::rbinom(M, 1L, 0.5) + 1L, seq_len(M))]
    from_mum <- pick(1L)
    from_dad <- pick(2L)
    genotypes[k, ] <- from_mum + from_dad
  }
}
snp_info <- data.frame(snp = paste0("rs", seq_len(M)),
                       chrom = rep(seq_len(n_chrom), each = M_per_chrom),
                       pos = rep(seq_len(M_per_chrom), n_chrom) * 1e5)
snp_info$snp[causal] <- names(causal)
dimnames(genotypes) <- list(ids, snp_info$snp)

## Traits ------------------------------------------------------------------

# Two genetic factors and two residual ones, so the traits are correlated
# through both G and E.
loadings_G <- cbind(rep(c(0.8, 0.3), each = P / 2), rep(c(0.2, 0.7), P / 2))
loadings_E <- cbind(seq(0.6, -0.3, length.out = P), rep(0.4, P))
G_true <- loadings_G %*% t(loadings_G) + diag(0.3, P)
E_true <- loadings_E %*% t(loadings_E) + diag(0.5, P)
dimnames(G_true) <- dimnames(E_true) <- list(traits, traits)

sex <- stats::rbinom(N, 1L, 0.5)
X <- cbind(intercept = 1, sex = sex)
rownames(X) <- ids
U <- t(chol(A)) %*% matrix(stats::rnorm(N * P), N, P) %*% chol(G_true)
noise <- matrix(stats::rnorm(N * P), N, P) %*% chol(E_true)
Y <- 10 + outer(sex, rep(0.5, P)) + U + noise

g <- function(name) genotypes[, name]
Y[, 1:3] <- Y[, 1:3] + 0.55 * g("rs_big3")
# Against the genetic correlation (traits 1-5 and 6-10 are positively
# correlated through the first factor): the direction where a joint test
# gains most over the traits taken one at a time.
Y <- Y + outer(g("rs_pleio"), rep(c(0.25, -0.25), each = P / 2))
Y[, 4:5] <- Y[, 4:5] + 0.9 * (g("rs_dom") == 1)
Y[, 8] <- Y[, 8] + 0.7 * g("rs_one")
dimnames(Y) <- list(ids, traits)

## G and E from trait subsets ----------------------------------------------

# The association scan needs the genetic and residual covariance of all 10
# traits. fit_G() estimates them from small trait subsets and reassembles
# them; with a REML engine this takes seconds. (In a real study, fit G and E
# on the phenotypes once, then scan as many genotype sets as you like.)
started <- Sys.time()
fit <- fit_G(Y, X = X, A = A, engine = "sommer", cores = cores, seed = seed)
fit_time <- as.numeric(difftime(Sys.time(), started, units = "secs"))

## Scans --------------------------------------------------------------------

started <- Sys.time()
# One joint test per SNP across the 10 traits (10 degrees of freedom), with
# empirical p-values and an FDR threshold from 99 genotype permutations.
multi <- indri_gwas(Y, genotypes, fit = fit, A = A, X = X,
                    snp_info = snp_info, n_permutations = 99L, fdr = 0.05,
                    seed = seed, cores = cores)
# The same SNPs trait by trait, adjusted for the effective number of traits.
single <- indri_gwas(Y, genotypes, A = A, X = X, snp_info = snp_info,
                     mode = "single_trait", adjust = "effective")
# Additive and dominance together, and dominance given additive.
dominance <- indri_gwas(Y, genotypes, fit = fit, A = A, X = X,
                        snp_info = snp_info, dominance = TRUE)
scan_time <- as.numeric(difftime(Sys.time(), started, units = "secs"))

## Results ------------------------------------------------------------------

cat("\n== indri genome scan: ", N, " individuals, ", P, " traits, ", M,
    " SNPs ==\n", sep = "")
cat(sprintf("fit_G (sommer):   %.0f s\n", fit_time))
cat(sprintf("three scans:      %.1f s (99 permutations included)\n", scan_time))
cat(sprintf("null variance ratio delta (multi-trait): %.2f\n", multi$delta))
cat("  (near 1: the fitted G and E describe these data)\n")

# Genomic inflation on the null SNPs: the median chi-square over its
# expectation. Near 1 when the test is calibrated.
lambda_gc <- function(p, df) {
  stats::qchisq(stats::median(p, na.rm = TRUE), df, lower.tail = FALSE) /
    stats::qchisq(0.5, df, lower.tail = FALSE)
}
null <- setdiff(seq_len(M), causal)
best_single <- apply(single$pvalues_adjust, 1L, min)
cat(sprintf("genomic inflation, null SNPs: multi-trait %.3f, single-trait %.3f\n",
            lambda_gc(multi$pvalues[null], P),
            lambda_gc(single$pvalues[null, ], 1)))

threshold <- 0.05 / M
rank_of <- function(p) rank(p, ties.method = "min")[causal]
cat(sprintf("\ncausal SNPs (Bonferroni threshold over SNPs: %.1e)\n", threshold))
print(data.frame(
  snp = names(causal),
  multi_trait_p = signif(multi$pvalues[causal], 3),
  multi_rank = rank_of(multi$pvalues),
  best_single_trait_p = signif(best_single[causal], 3),
  single_rank = rank_of(best_single),
  dominance_p = signif(dominance$pvalues_dominance[causal], 3),
  row.names = NULL
), row.names = FALSE)
cat("\n  rs_pleio: no single trait carries much signal; the joint test adds it up\n")
cat("  rs_one: one trait only; the joint test spends 10 df on it\n")
cat("  rs_dom: the additive tests see little of it; the dominance test finds it\n")

cat("\nadditive and dominance effects at rs_dom (true d = 0.9 on traits 4-5):\n")
print(round(rbind(a = dominance$a["rs_dom", ], d = dominance$d["rs_dom", ]), 2))

cat("\nempirical p-values of the causal SNPs (99 permutations; the smallest",
    "possible is 0.01):\n")
print(signif(multi$pvalues_empirical[causal], 2))
cat("  (permuting genotypes across relatives ignores the relatedness; a check,",
    "not a replacement)\n\n")
print(summary(multi, n = 5))

## Figures -----------------------------------------------------------------

# plot() returns ggplot objects: a Manhattan plot by default, a QQ plot with
# type = "qq". Mark the planted SNPs.
dir.create("examples/figures", showWarnings = FALSE, recursive = TRUE)
save <- function(plot, name, height = 3)
  ggplot2::ggsave(file.path("examples/figures", name), plot, width = 9,
                  height = height, dpi = 120)
save(plot(multi, highlight = names(causal)), "gwas_manhattan_multi.png")
save(plot(single, highlight = names(causal)), "gwas_manhattan_single.png")
save(plot(dominance, which = "dominance", highlight = names(causal)),
     "gwas_manhattan_dominance.png")
save(plot(multi, type = "qq"), "gwas_qq_multi.png", height = 4)
cat("\nfigures: examples/figures/gwas_*.png\n")
