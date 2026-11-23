# Scans with a litter (nurse) effect: calibration and power of the three ways
# of handling it.
#
# Two generations of random mating from unrelated founders, so that the
# scanned families are related to each other (cousins), as in an advanced
# intercross; the scanned generation is families of full sibs, each reared
# together as one litter (a nurse level is one family, as when whole litters
# are cross-fostered or kept with their mother). Phenotypes: a polygenic
# background (pedigree A), a litter effect with its own trait covariance
# (negative between the traits, like catch-up growth), and noise; null SNPs
# are passed down the same pedigree.
# Per seed, G (and the litter covariance) are estimated with fit_G() (sommer)
# for each model, then three multi-trait scans:
#   ignored  the litter effect is not modelled (it is taken for relatedness)
#   fixed    one dummy per litter in X (a within-family test)
#   random   indri_gwas(random = list(nurse = litter)), the low-rank path
# Null SNPs give the type-I error and the genomic inflation; one planted QTL
# per seed gives the power.
#
# Settings: INDRI_GWAS_RANDOM_SEEDS (default 20), _N_FAM (150 scanned
# litters, from 75 founder pairs with 4 offspring each), _SIZE (6),
# _M (2000 null SNPs), _QTL (0.25, the QTL effect in SD of each trait per
# allele). Writes results/gwas_random_calibration.csv.
suppressPackageStartupMessages(library(indri))
source("examples/gwas_simulate.R")
env <- function(name, default) as.numeric(Sys.getenv(paste0("INDRI_GWAS_RANDOM_", name), default))
seeds <- seq_len(env("SEEDS", 20)); n_fam <- env("N_FAM", 150); size <- env("SIZE", 6)
M <- env("M", 2000); qtl <- env("QTL", 0.25)
P <- 2L; traits <- c("t1", "t2")
G <- matrix(c(0.6, 0.2, 0.2, 0.5), 2, dimnames = list(traits, traits))
C <- matrix(c(0.5, -0.25, -0.25, 0.4), 2, dimnames = list(traits, traits))
E <- matrix(c(1, 0.1, 0.1, 1), 2, dimnames = list(traits, traits))
options(indri.warn_rhat = FALSE)

# Mendelian genotypes down a pedigree at M unlinked SNPs (founders in HWE).
pedigree_genotypes <- function(ped, maf) {
  n <- nrow(ped); M <- length(maf)
  hap <- array(0L, c(n, 2L, M))
  for (i in seq_len(n)) {
    if (is.na(ped$dam[i])) {
      hap[i, , ] <- matrix(stats::rbinom(2L * M, 1L, rep(maf, each = 2L)), 2L)
    } else {
      for (h in 1:2) {
        par <- if (h == 1L) ped$dam[i] else ped$sire[i]
        hap[i, h, ] <- hap[cbind(par, stats::rbinom(M, 1L, 0.5) + 1L, seq_len(M))]
      }
    }
  }
  hap[, 1, ] + hap[, 2, ]
}

one_seed <- function(seed) {
  set.seed(seed)
  n_pairs0 <- n_fam / 2
  founders <- data.frame(id = seq_len(2 * n_pairs0), dam = NA, sire = NA)
  g1_dam <- rep(seq(1, 2 * n_pairs0, 2), each = 4); g1_sire <- g1_dam + 1L
  g1 <- data.frame(id = 2 * n_pairs0 + seq_along(g1_dam), dam = g1_dam, sire = g1_sire,
                   sex = rep(c(0L, 1L), length(g1_dam) / 2))
  # scanned litters: a female and a male of the first generation from different families
  fem <- sample(g1$id[g1$sex == 0L]); mal <- g1$id[g1$sex == 1L]
  mal <- mal[order(match(g1$dam[match(mal, g1$id)], unique(g1$dam)))]
  mal <- mal[c(2:length(mal), 1)]                       # never a full sib of its mate
  pairs <- data.frame(dam = fem, sire = sample(mal))
  fam_of <- function(x) g1$dam[match(x, g1$id)]
  while (any(fam_of(pairs$dam) == fam_of(pairs$sire))) pairs$sire <- sample(pairs$sire)
  start <- max(g1$id)
  g2 <- data.frame(id = start + seq_len(n_fam * size), dam = rep(pairs$dam, each = size),
                   sire = rep(pairs$sire, each = size))
  ped <- rbind(founders, g1[, 1:3], g2)
  Afull <- stanimal::makeA(data.frame(ped$id, ped$sire, ped$dam))
  scanned <- as.character(g2$id)
  N <- length(scanned)
  ids <- paste0("id", scanned)
  A <- Afull[scanned, scanned]; dimnames(A) <- list(ids, ids)
  litter <- factor(rep(seq_len(n_fam), each = size))
  sex <- stats::rbinom(N, 1L, 0.5)
  geno <- pedigree_genotypes(ped, stats::runif(M + 1L, 0.1, 0.5))[g2$id - 0L, , drop = FALSE]
  dimnames(geno) <- list(ids, c(paste0("null", seq_len(M)), "qtl"))
  g_q <- geno[, "qtl"] - mean(geno[, "qtl"])
  Y <- t(chol(A + diag(1e-8, N))) %*% matrix(stats::rnorm(N * P), N, P) %*% chol(G) +
    matrix(stats::rnorm(n_fam * P), n_fam, P)[as.integer(litter), ] %*% chol(C) +
    matrix(stats::rnorm(N * P), N, P) %*% chol(E) +
    outer(sex, rep(0.5, P)) + outer(g_q, rep(qtl, P))
  dimnames(Y) <- list(ids, traits)
  X <- cbind(intercept = 1, sex = sex); rownames(X) <- ids
  Xf <- cbind(X, stats::model.matrix(~ litter)[, -1]); rownames(Xf) <- ids
  fit_args <- list(engine = "sommer", subsets = matrix(1:2, 1), assembly = "average",
                   cores = 1L, progress = FALSE, seed = seed)
  f_ign <- do.call(fit_G, c(list(Y, X = X, A = A), fit_args))
  f_fix <- do.call(fit_G, c(list(Y, X = Xf, A = A), fit_args))
  f_rnd <- do.call(fit_G, c(list(Y, X = X, A = A, random = list(nurse = litter)), fit_args))
  scans <- list(
    ignored = indri_gwas(Y, geno, fit = f_ign, A = A, X = X),
    fixed = indri_gwas(Y, geno, fit = f_fix, A = A, X = Xf),
    random = indri_gwas(Y, geno, fit = f_rnd, A = A, X = X, random = list(nurse = litter)))
  fits <- list(ignored = f_ign, fixed = f_fix, random = f_rnd)
  do.call(rbind, lapply(names(scans), function(k) {
    p <- scans[[k]]$pvalues; p0 <- p[seq_len(M)]
    data.frame(seed = seed, model = k, delta = scans[[k]]$delta,
               rate_05 = mean(p0 < 0.05), rate_01 = mean(p0 < 0.01), rate_001 = mean(p0 < 0.001),
               lambda_gc = stats::median(stats::qchisq(p0, P, lower.tail = FALSE)) /
                 stats::qchisq(0.5, P),
               qtl_log10p = -log10(p[["qtl"]]), qtl_hit_001 = p[["qtl"]] < 1e-3,
               G_t1 = fits[[k]]$G[1, 1], rG = fits[[k]]$corrG[1, 2])
  }))
}

t0 <- Sys.time()
res <- do.call(rbind, lapply(seeds, one_seed))
dir.create("results", showWarnings = FALSE)
utils::write.csv(res, "results/gwas_random_calibration.csv", row.names = FALSE)
cat(sprintf("%d seeds, %d families of %d, %d null SNPs, QTL %.2f SD per allele, %.1f min\n",
            length(seeds), n_fam, size, M, qtl, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("true G[1,1] %.2f, rG %.2f; litter correlation %.2f\n", G[1, 1], stats::cov2cor(G)[1, 2],
            stats::cov2cor(C)[1, 2]))
summ <- aggregate(cbind(rate_05, rate_01, rate_001, lambda_gc, qtl_log10p, qtl_hit_001, G_t1, rG, delta) ~ model,
                  res, mean)
print(format(summ, digits = 3), row.names = FALSE)
