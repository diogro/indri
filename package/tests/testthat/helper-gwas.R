# Reference implementations for indri_gwas(), written the slow, obvious way:
# dense covariance matrices, explicit solves, no rotation. Small N only.
#
# Every test profiles an overall scale sigma^2, so the statistic is a
# likelihood ratio with sigma^2 at its ML value under each model:
#   multi-trait LMM   V = G (x) A + delta E (x) I       stat = NP log(q0 / q1)
#   multi-trait LM    V = E (x) I                       stat = NP log(q0 / q1)
#   single-trait LMM  V = A + delta I, per trait        stat =  N log(q0 / q1)
#   single-trait LM   V = I                             stat =  N log(q0 / q1)
# where q is the generalised residual sum of squares r' V^-1 r. vec() stacks
# the columns of Y, so the covariance of vec(Y) is G (x) A.

dense_gls <- function(y, D, Vi) {
  A <- crossprod(D, Vi %*% D)
  b <- solve(A, crossprod(D, Vi %*% y))
  r <- y - D %*% b
  list(beta = drop(b), q = drop(crossprod(r, Vi %*% r)), A = A)
}

# REML criterion for the scale-free variance ratio, sigma^2 profiled.
dense_reml <- function(log_delta, y, D, V_of) {
  V <- V_of(exp(log_delta))
  Vi <- solve(V)
  fit <- dense_gls(y, D, Vi)
  n <- length(y) - ncol(D)
  0.5 * (determinant(V)$modulus + determinant(fit$A)$modulus +
           n * log(fit$q / n))
}

dense_delta <- function(y, D, V_of) {
  exp(stats::optimize(dense_reml, c(-12, 12), y = y, D = D, V_of = V_of,
                      tol = 1e-10)$minimum)
}

# Multi-trait any-effect test. G = NULL gives the LM (V = E (x) I).
dense_multi <- function(Y, X, g, E, G = NULL, A = NULL, delta = NULL) {
  N <- nrow(Y); P <- ncol(Y)
  y <- c(Y)
  D0 <- kronecker(diag(P), X)
  D1 <- cbind(D0, kronecker(diag(P), matrix(g, ncol = 1L)))
  V_of <- if (is.null(G)) function(d) kronecker(E, diag(N))
          else function(d) kronecker(G, A) + d * kronecker(E, diag(N))
  if (!is.null(G) && is.null(delta)) delta <- dense_delta(y, D0, V_of)
  Vi <- solve(V_of(if (is.null(delta)) 1 else delta))
  f0 <- dense_gls(y, D0, Vi)
  f1 <- dense_gls(y, D1, Vi)
  stat <- N * P * log(f0$q / f1$q)
  list(p = stats::pchisq(stat, P, lower.tail = FALSE), stat = stat,
       beta = utils::tail(f1$beta, P), delta = delta)
}

# Single-trait test for one trait. A = NULL gives the LM.
dense_single <- function(y, X, g, A = NULL, delta = NULL) {
  N <- length(y)
  D1 <- cbind(X, g)
  V_of <- if (is.null(A)) function(d) diag(N)
          else function(d) A + d * diag(N)
  if (!is.null(A) && is.null(delta)) delta <- dense_delta(y, X, V_of)
  Vi <- solve(V_of(if (is.null(delta)) 1 else delta))
  f0 <- dense_gls(y, X, Vi)
  f1 <- dense_gls(y, D1, Vi)
  stat <- N * log(f0$q / f1$q)
  list(p = stats::pchisq(stat, 1, lower.tail = FALSE), stat = stat,
       beta = unname(utils::tail(f1$beta, 1L)), delta = delta)
}

# A small association data set: families of full sibs, P correlated traits,
# M SNPs of which `causal` affect every trait. Names on everything, with
# traits that do not sort in column order.
simulate_gwas_data <- function(n_fam = 8, size = 5, P = 3L, M = 6L,
                               causal = 1L, effect = 0.6, seed = 1) {
  set.seed(seed)
  N <- n_fam * size
  A <- kronecker(diag(n_fam), matrix(0.5, size, size))
  diag(A) <- 1
  G0 <- matrix(0.4, P, P); diag(G0) <- 1
  E0 <- matrix(-0.2, P, P); diag(E0) <- 1
  ids <- paste0("id", seq_len(N))
  traits <- c("wing", "tail", "beak", "leg", "neck")[seq_len(P)]
  snps <- paste0("rs", seq_len(M) * 7L)
  geno <- matrix(stats::rbinom(N * M, 2L, 0.3), N, M,
                 dimnames = list(ids, snps))
  X <- cbind(intercept = 1, x = stats::rnorm(N))
  rownames(X) <- ids
  U <- t(chol(A)) %*% matrix(stats::rnorm(N * P), N, P) %*% chol(G0)
  noise <- matrix(stats::rnorm(N * P), N, P) %*% chol(E0)
  Y <- X %*% matrix(c(5, 0.5), 2L, P) + U + noise
  if (causal > 0L)
    Y <- Y + geno[, seq_len(causal), drop = FALSE] %*%
      matrix(effect, causal, P)
  dimnames(Y) <- list(ids, traits)
  dimnames(A) <- list(ids, ids)
  dimnames(G0) <- dimnames(E0) <- list(traits, traits)
  list(Y = Y, X = X, A = A, geno = geno, G = G0, E = E0, N = N, P = P, M = M)
}

# limmbo's association results on its own test data (250 samples, 10 traits,
# 20 SNPs), made by fixtures/limmbo/make_limmbo_fixtures.py.
limmbo_fixture <- function() {
  dir <- testthat::test_path("fixtures", "limmbo")
  rd <- function(f) as.matrix(utils::read.csv(file.path(dir, f), row.names = 1,
                                              check.names = FALSE))
  ld <- utils::read.csv(file.path(dir, "log_delta.csv"),
                        stringsAsFactors = FALSE)
  delta <- function(test) {
    v <- exp(ld$log_delta[ld$test == test])
    if (length(v) > 1L) stats::setNames(v, ld$trait[ld$test == test]) else v
  }
  list(
    Y = rd("Y.csv"), X = rd("X.csv"), geno = rd("genotypes.csv"),
    A = readRDS(file.path(dir, "A.rds")), Cg = rd("Cg.csv"), Cn = rd("Cn.csv"),
    Cn_lm = rd("Cn_lm.csv"),
    snp_info = utils::read.csv(file.path(dir, "snp_info.csv"),
                               stringsAsFactors = FALSE),
    p = list(
      multi_lmm = utils::read.csv(file.path(dir, "multi_trait_lmm_pvalues.csv"),
                                  stringsAsFactors = FALSE)$pvalue,
      multi_lm = utils::read.csv(file.path(dir, "multi_trait_lm_pvalues.csv"),
                                 stringsAsFactors = FALSE)$pvalue,
      single_lmm = rd("single_trait_lmm_pvalues.csv"),
      single_lm = rd("single_trait_lm_pvalues.csv")),
    beta = list(
      multi_lmm = rd("multi_trait_lmm_betas.csv"),
      multi_lm = rd("multi_trait_lm_betas.csv"),
      single_lmm = rd("single_trait_lmm_betas.csv"),
      single_lm = rd("single_trait_lm_betas.csv")),
    padjust = list(
      lmm_bonferroni = rd("single_trait_lmm_padjust_bonferroni.csv"),
      lmm_effective = rd("single_trait_lmm_padjust_effective.csv")),
    effective_tests = utils::read.csv(
      file.path(dir, "effective_tests.csv"))$effective_tests,
    delta = list(multi_lmm = delta("multi_trait_lmm"),
                 single_lmm = delta("single_trait_lmm"))
  )
}

expect_rel_equal <- function(object, expected, tol, ...) {
  rel <- max(abs(unname(object) - unname(expected)) /
               pmax(abs(unname(expected)), 1e-300))
  testthat::expect_lt(rel, tol, ...)
}

# Likelihood ratio between two nested fixed-effect designs, for the
# dominance tests: X0 (null) inside X1 (alternative), the extra columns of X1
# being the SNP terms. Multi-trait: every trait gets every column. `delta`
# must be given (the scans fit it once, without any SNP). Returns the p-value
# and the effects of the extra columns, one row per column and one column
# per trait.
dense_nested_multi <- function(Y, X0, X1, E, G = NULL, A = NULL, delta = 1) {
  N <- nrow(Y); P <- ncol(Y)
  y <- c(Y)
  V <- if (is.null(G)) kronecker(E, diag(N))
       else kronecker(G, A) + delta * kronecker(E, diag(N))
  Vi <- solve(V)
  f0 <- dense_gls(y, kronecker(diag(P), X0), Vi)
  f1 <- dense_gls(y, kronecker(diag(P), X1), Vi)
  extra <- ncol(X1) - ncol(X0)
  stat <- N * P * log(f0$q / f1$q)
  b <- matrix(f1$beta, ncol(X1), P)
  list(p = stats::pchisq(stat, P * extra, lower.tail = FALSE),
       beta = b[ncol(X0) + seq_len(extra), , drop = FALSE])
}

dense_nested_single <- function(y, X0, X1, A = NULL, delta = 1) {
  N <- length(y)
  Vi <- solve(if (is.null(A)) diag(N) else A + delta * diag(N))
  f0 <- dense_gls(y, X0, Vi)
  f1 <- dense_gls(y, X1, Vi)
  extra <- ncol(X1) - ncol(X0)
  stat <- N * log(f0$q / f1$q)
  list(p = stats::pchisq(stat, extra, lower.tail = FALSE),
       beta = unname(utils::tail(f1$beta, extra)))
}

# Heterozygote indicator: the dominance coding of a 0/1/2 genotype.
heterozygous <- function(g) as.numeric(g == 1)

# Dense reference for scans with an environment-specific background: the
# covariance of vec(Y) (trait-major) is built entry by entry,
#   Cov(y_ik, y_jl) = A_ij Gs[(k, env_i), (l, env_j)] + [i == j] delta E_env_i[k, l],
# with Gs over the character states "trait:level" and Es block diagonal.
dense_env_V <- function(A, env, Gs, Es, traits, delta = 1) {
  N <- nrow(A); P <- length(traits)
  state <- function(k, i) paste(traits[k], env[i], sep = ":")
  V <- matrix(0, N * P, N * P)
  for (k in seq_len(P)) for (l in seq_len(P)) {
    rows <- (k - 1) * N + seq_len(N)
    cols <- (l - 1) * N + seq_len(N)
    g <- outer(seq_len(N), seq_len(N), function(i, j)
      A[cbind(i, j)] * Gs[cbind(state(k, i), state(l, j))])
    e <- diag(vapply(seq_len(N), function(i) Es[state(k, i), state(l, i)],
                     numeric(1)))
    V[rows, cols] <- g + delta * e
  }
  V
}

dense_env_multi <- function(Y, X0, X1, A, env, Gs, Es, delta = 1) {
  N <- nrow(Y); P <- ncol(Y)
  V <- dense_env_V(A, env, Gs, Es, colnames(Y), delta)
  Vi <- solve(V)
  y <- c(Y)
  f0 <- dense_gls(y, kronecker(diag(P), X0), Vi)
  f1 <- dense_gls(y, kronecker(diag(P), X1), Vi)
  extra <- ncol(X1) - ncol(X0)
  stat <- N * P * log(f0$q / f1$q)
  b <- matrix(f1$beta, ncol(X1), P)
  list(p = stats::pchisq(stat, P * extra, lower.tail = FALSE),
       beta = b[ncol(X0) + seq_len(extra), , drop = FALSE])
}

dense_env_delta <- function(Y, X, A, env, Gs, Es) {
  P <- ncol(Y)
  dense_delta(c(Y), kronecker(diag(P), X), function(d)
    dense_env_V(A, env, Gs, Es, colnames(Y), d))
}

# Sibs of both sexes with a sex-specific G (cross-sex genetic correlation
# r_between), and known character-state Gs and block-diagonal Es.
simulate_env_gwas <- function(n_fam = 12, size = 5, P = 2L, M = 5L,
                              r_between = 0.5, seed = 1) {
  set.seed(seed)
  N <- n_fam * size
  A <- kronecker(diag(n_fam), matrix(0.5, size, size))
  diag(A) <- 1
  ids <- paste0("id", seq_len(N))
  dimnames(A) <- list(ids, ids)
  env <- ifelse(stats::rbinom(N, 1L, 0.5) == 1L, "m", "f")
  names(env) <- ids
  traits <- c("wing", "tail", "beak")[seq_len(P)]
  states <- as.vector(outer(traits, c("f", "m"), paste, sep = ":"))
  G0 <- matrix(0.3, P, P); diag(G0) <- 1
  Gs <- kronecker(matrix(c(1, r_between, r_between, 1), 2), G0)
  E0 <- matrix(-0.2, P, P); diag(E0) <- 1
  Es <- kronecker(diag(2), E0)
  dimnames(Gs) <- dimnames(Es) <- list(states, states)
  # a symmetric square root, so that a singular Gs (r_between = 1) works
  eg <- eigen(Gs, symmetric = TRUE)
  root <- eg$vectors %*% (sqrt(pmax(eg$values, 0)) * t(eg$vectors))
  U <- t(chol(A)) %*% matrix(stats::rnorm(N * 2 * P), N, 2 * P) %*% root
  own <- ifelse(env == "m", P, 0)
  Y <- t(vapply(seq_len(N), function(i) U[i, own[i] + seq_len(P)],
                numeric(P))) + matrix(stats::rnorm(N * P), N, P) %*% chol(E0)
  dimnames(Y) <- list(ids, traits)
  X <- cbind(intercept = 1, male = as.numeric(env == "m"),
             x = stats::rnorm(N))
  rownames(X) <- ids
  geno <- matrix(stats::rbinom(N * M, 2L, 0.4), N, M,
                 dimnames = list(ids, paste0("rs", seq_len(M))))
  fit <- structure(list(G = Gs, E = Es, levels = c("f", "m"),
                        traits = traits, states = states,
                        env = factor(env, levels = c("f", "m"))),
                   class = "indri_fit_env")
  names(fit$env) <- ids
  list(Y = Y, X = X, A = A, env = env, Gs = Gs, Es = Es, geno = geno,
       fit = fit, N = N, P = P, M = M)
}
