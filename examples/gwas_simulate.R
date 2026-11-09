# Shared simulation for the GWAS studies (gwas_calibration.R, gwas_power.R):
# families of full sibs with Mendelian genotypes at unlinked SNPs, and P
# traits with random genetic and residual covariances. Not part of the
# package; sourced by the studies.

# Genotypes of `size` full sibs in each of `n_fam` families at M unlinked
# SNPs with allele frequencies `maf`: parents draw two alleles each, every
# child takes one from each parent.
gwas_simulate_genotypes <- function(n_fam, size, maf) {
  M <- length(maf)
  geno <- matrix(0L, n_fam * size, M)
  for (f in seq_len(n_fam)) {
    parents <- array(stats::rbinom(4L * M, 1L, rep(maf, each = 4L)),
                     c(2L, 2L, M))
    pick <- function(parent)
      parents[cbind(parent, stats::rbinom(M, 1L, 0.5) + 1L, seq_len(M))]
    for (k in (f - 1L) * size + seq_len(size))
      geno[k, ] <- pick(1L) + pick(2L)
  }
  geno
}

# The relationship matrix of that design.
gwas_sibship_A <- function(n_fam, size) {
  A <- kronecker(diag(n_fam), matrix(0.5, size, size))
  diag(A) <- 1
  A
}

# A population: G, E, a design with sex, the background phenotypes (no
# causal SNP) and M null SNPs. `G_method` is passed to
# stanimal::random_covariance(): "partial" gives strong, nearly singular
# genetic correlation structures at many traits (smallest eigenvalue ~1e-8
# at 30 traits), "lkj" (LKJ(2)) weaker and better conditioned ones.
gwas_simulate_population <- function(P, n_fam, M, seed, size = 5L,
                                     G_method = "partial") {
  set.seed(seed)
  N <- n_fam * size
  traits <- paste0("trait", seq_len(P))
  ids <- paste0("id", seq_len(N))
  G <- if (G_method == "lkj") {
    stanimal::random_covariance(n_traits = P, method = "lkj", lkj_prior = 2,
                                min_var = 0.8, max_var = 1.4,
                                trait_names = traits)
  } else {
    stanimal::random_covariance(n_traits = P, method = "partial",
                                min_var = 0.8, max_var = 1.4,
                                trait_names = traits)
  }
  E <- stanimal::random_covariance(n_traits = P, method = "lkj",
                                   lkj_prior = 2, min_var = 0.8,
                                   max_var = 1.6, trait_names = traits)
  A <- gwas_sibship_A(n_fam, size)
  dimnames(A) <- list(ids, ids)
  X <- cbind(intercept = 1, sex = stats::rbinom(N, 1L, 0.5))
  rownames(X) <- ids
  U <- t(chol(A)) %*% matrix(stats::rnorm(N * P), N, P) %*% chol(G)
  noise <- matrix(stats::rnorm(N * P), N, P) %*% chol(E)
  Y <- 10 + outer(X[, "sex"], rep(0.5, P)) + U + noise
  dimnames(Y) <- list(ids, traits)
  geno <- gwas_simulate_genotypes(n_fam, size, stats::runif(M, 0.05, 0.5))
  dimnames(geno) <- list(ids, paste0("rs", seq_len(M)))
  list(Y = Y, X = X, A = A, G = G, E = E, genotypes = geno, N = N, P = P,
       n_fam = n_fam, size = size)
}

# Settings shared by the studies: the G generator and a tag for the output
# file names ("" for the default "partial" generator).
gwas_study_G_method <- function() Sys.getenv("INDRI_GWAS_STUDY_G", "partial")
gwas_study_tag <- function() {
  m <- gwas_study_G_method()
  extra <- Sys.getenv("INDRI_GWAS_STUDY_TAG", "")
  paste0(if (m == "partial") "" else paste0("_", m),
         if (nzchar(extra)) paste0("_", extra) else "")
}

# E from a shrunk residual covariance: Schafer and Strimmer's (2005)
# estimator (corpcor::cov.shrink: correlations towards 0, variances towards
# their median, intensities by the Ledoit-Wolf lemma) of the residuals of Y
# on X, minus c G, floored like fit_G()'s E. For the studies only.
gwas_shrunk_E <- function(Y, X, G, A) {
  R <- stats::lm.fit(X, Y)$residuals
  S <- unclass(corpcor::cov.shrink(R, verbose = FALSE))
  S <- S * (nrow(R) - 1) / (nrow(R) - ncol(X))
  M_A <- sum(diag(A)) - sum(diag(solve(crossprod(X), crossprod(X, A %*% X))))
  E <- S - (M_A / (nrow(R) - ncol(X))) * G
  e <- eigen((E + t(E)) / 2, symmetric = TRUE)
  E <- e$vectors %*% (pmax(e$values, 1e-3 * e$values[1]) * t(e$vectors))
  dimnames(E) <- dimnames(G)
  (E + t(E)) / 2
}
