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
# causal SNP) and M null SNPs.
gwas_simulate_population <- function(P, n_fam, M, seed, size = 5L) {
  set.seed(seed)
  N <- n_fam * size
  traits <- paste0("trait", seq_len(P))
  ids <- paste0("id", seq_len(N))
  G <- stanimal::random_covariance(n_traits = P, method = "partial",
                                   min_var = 0.8, max_var = 1.4,
                                   trait_names = traits)
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
