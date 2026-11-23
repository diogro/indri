# Association scans with extra random effects (indri_gwas(random =),
# indri_gwas_haplotypes(random =)).
#
# The covariance of vec(Y) is, up to the profiled scale sigma^2,
#   V = G (x) A + sum_r C_r (x) Z_r Z_r' + delta E (x) I,
# with Z_r the incidence matrix of random term r (L_r levels) and C_r its trait
# covariance. G (x) A + delta E (x) I is whitened as in the other scans: the
# traits by T (E -> I, G -> diag(d)), the individuals by the eigenvectors Q of
# A, giving a diagonal covariance D[i, c] = lambda_i d_c + delta. The random
# terms then add a low-rank part, W W' with
#   W = [Lc_1 (x) Q'Z_1, Lc_2 (x) Q'Z_2, ...],   T' C_r T = Lc_r Lc_r',
# of width w = P sum_r L_r. By the Woodbury identity, with S = I + W' D^-1 W,
#   a' V^-1 b = a' D^-1 b - (R^-T W' D^-1 a)' (R^-T W' D^-1 b),  S = R'R,
#   log |V| = log |D| + log |S|,
# so every quadratic form of the GLS is the diagonal one minus a correction
# through the w x w Cholesky factor R. The scan never forms the NP x NP V: it
# costs the rotation of the markers, O(N w) per marker column, and O(w^2) per
# triangular solve. REML for delta repeats the null fit (one Cholesky of S)
# per evaluation.
#
# Columns of the designs are trait-specific, e_c (x) v for trait c and an
# N-vector v (rotated). For such a column, W' D^-1 (e_c (x) v) is, block by
# block, Lc_r[c, ]' (x) (Q'Z_r)' (v / D[, c]).

# Validate `random` (named list of factors, one value per individual) and
# match `random_cov` (named list of P x P matrices, in the units of Y) to it.
# Missing levels count as unshared, as in fit_G(random =). Returns a list of
# terms with the incidence matrix `Z` (N x L) and the covariance `C`.
.indri_gwas_random_terms <- function(random, random_cov, traits, N) {
  if (!is.list(random) || !length(random) || is.null(names(random)) ||
      any(names(random) == "") || anyDuplicated(names(random)))
    stop("`random` must be a named list of factors, one per random term.",
         call. = FALSE)
  if (is.null(random_cov))
    stop("`random` needs the covariance of each term: pass a `fit` from ",
         "fit_G(random = ) (its `random`), or `random_cov`.", call. = FALSE)
  if (!is.list(random_cov) || is.null(names(random_cov)))
    stop("`random_cov` must be a named list of trait covariance matrices, one ",
         "per term of `random`.", call. = FALSE)
  missing <- setdiff(names(random), names(random_cov))
  if (length(missing))
    stop("`random_cov` has no covariance for the random term(s) ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  lapply(stats::setNames(names(random), names(random)), function(nm) {
    f <- random[[nm]]
    if (length(f) != N)
      stop("`random$", nm, "` has ", length(f), " values; `Y` has ", N,
           " individuals.", call. = FALSE)
    f <- as.character(f)
    gone <- is.na(f)
    f[gone] <- paste0(".indri_unshared_", seq_len(sum(gone)))
    f <- factor(f, levels = unique(f))
    Z <- matrix(0, N, nlevels(f))
    Z[cbind(seq_len(N), as.integer(f))] <- 1
    C <- .indri_gwas_trait_matrix(as.matrix(random_cov[[nm]]), traits,
                                  paste0("random_cov$", nm))
    C <- (C + t(C)) / 2
    if (min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) <
        -1e-8 * max(abs(diag(C)), 1e-12))
      stop("`random_cov$", nm, "` must be positive semi-definite.", call. = FALSE)
    list(Z = Z, C = C, levels = nlevels(f))
  })
}

# A factor B with B B' = M for a positive semi-definite M.
.indri_psd_factor <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
}

# W' D^-1 (e_c (x) v) for the columns v of `V` (N x k), trait c: w x k.
.indri_lr_m <- function(V, dinv_c, c, Zt, Lc) {
  Vd <- V * dinv_c
  do.call(rbind, lapply(seq_along(Zt), function(r)
    kronecker(matrix(Lc[[r]][c, ], ncol = 1L), crossprod(Zt[[r]], Vd))))
}

# S = I + W' D^-1 W.
.indri_lr_S <- function(Dinv, Zt, Lc) {
  P <- ncol(Dinv)
  blocks <- lapply(seq_along(Zt), function(r) lapply(seq_along(Zt), function(s) {
    B <- 0
    for (c in seq_len(P))
      B <- B + kronecker(tcrossprod(Lc[[r]][c, ], Lc[[s]][c, ]),
                         crossprod(Zt[[r]], Zt[[s]] * Dinv[, c]))
    B
  }))
  S <- do.call(rbind, lapply(blocks, function(row) do.call(cbind, row)))
  S <- (S + t(S)) / 2
  diag(S) <- diag(S) + 1
  S
}

.indri_block_diag <- function(mats) {
  n <- vapply(mats, nrow, 0L)
  out <- matrix(0, sum(n), sum(n))
  at <- 0L
  for (m in mats) {
    i <- at + seq_len(nrow(m))
    out[i, i] <- m
    at <- at + nrow(m)
  }
  out
}

# The null model (no marker) for a given D: everything the per-marker tests
# reuse, and the REML criterion.
.indri_lr_null <- function(Ys, Xs, D, Zt, Lc) {
  P <- ncol(Ys); J <- ncol(Xs)
  Dinv <- 1 / D
  R <- chol(.indri_lr_S(Dinv, Zt, Lc))
  Mx <- do.call(cbind, lapply(seq_len(P), function(c) .indri_lr_m(Xs, Dinv[, c], c, Zt, Lc)))
  Ox <- backsolve(R, Mx, transpose = TRUE)
  my <- Reduce(`+`, lapply(seq_len(P), function(c)
    .indri_lr_m(Ys[, c, drop = FALSE], Dinv[, c], c, Zt, Lc)))
  oy <- drop(backsolve(R, my, transpose = TRUE))
  XX <- .indri_block_diag(lapply(seq_len(P), function(c) crossprod(Xs, Xs * Dinv[, c]))) -
    crossprod(Ox)
  XX <- (XX + t(XX)) / 2
  Xy <- unlist(lapply(seq_len(P), function(c) crossprod(Xs, Ys[, c] * Dinv[, c]))) -
    drop(crossprod(Ox, oy))
  yy <- sum(Ys^2 * Dinv) - sum(oy^2)
  RX <- chol(XX)
  beta <- backsolve(RX, backsolve(RX, Xy, transpose = TRUE))
  q0 <- yy - sum(Xy * beta)
  n <- length(Ys) - P * J
  list(Dinv = Dinv, R = R, Ox = Ox, oy = oy, beta = drop(beta), q0 = q0,
       XXinv = chol2inv(RX),
       reml = 0.5 * (sum(log(D)) + 2 * sum(log(diag(R))) +
                       2 * sum(log(diag(RX))) + n * log(q0 / n)))
}

# One scan: whitening of the traits, rotation, delta, null fit, then the
# nested test of each marker's terms. `terms_of(idx)` returns a list of K
# matrices (N x length(idx)), the marker terms of positions idx, unrotated.
# Returns stat-ready pieces: q0, q1 (per marker), the coefficients (markers x
# K x P, on the traits of Y) and delta.
.indri_gwas_random_scan <- function(Y, X, G, E, terms, rot, delta, terms_of,
                                    M, K, chunk_size, skip = rep(FALSE, M)) {
  N <- nrow(Y); P <- ncol(Y); J <- ncol(X)
  L <- .indri_gwas_chol(E, "E")
  Gw <- backsolve(L, t(backsolve(L, G, upper.tri = FALSE)), upper.tri = FALSE)
  eg <- eigen((Gw + t(Gw)) / 2, symmetric = TRUE)
  d <- pmax(eg$values, 0)
  Tm <- backsolve(t(L), eg$vectors)
  Ti <- solve(Tm)
  Ys <- .indri_gwas_rotate(rot$Q, Y %*% Tm)
  Xs <- .indri_gwas_rotate(rot$Q, X)
  Zt <- lapply(terms, function(tm) .indri_gwas_rotate(rot$Q, tm$Z))
  Lc <- lapply(terms, function(tm) .indri_psd_factor(crossprod(Tm, tm$C %*% Tm)))
  D_of <- function(dl) outer(rot$lambda, d) + dl
  if (is.null(delta))
    delta <- exp(stats::optimize(function(ld) .indri_lr_null(Ys, Xs, D_of(exp(ld)), Zt, Lc)$reml,
                                 c(-12, 12), tol = 1e-9)$minimum)
  null <- .indri_lr_null(Ys, Xs, D_of(delta), Zt, Lc)
  Dinv <- null$Dinv; R <- null$R

  q <- P * K                          # term columns, ordered (trait c, term k)
  ck <- expand.grid(k = seq_len(K), c = seq_len(P))
  q1 <- rep(NA_real_, M)
  coef <- array(NA_real_, c(M, K, P))
  for (idx in .indri_gwas_chunks(M, chunk_size)) {
    m <- length(idx)
    Tr <- lapply(terms_of(idx), function(t) .indri_gwas_rotate(rot$Q, t))
    Om <- vector("list", q); GTx <- vector("list", q); GTy <- vector("list", q)
    for (a in seq_len(q)) {
      c <- ck$c[a]; k <- ck$k[a]
      Om[[a]] <- backsolve(R, .indri_lr_m(Tr[[k]], Dinv[, c], c, Zt, Lc), transpose = TRUE)
      gx <- -crossprod(Om[[a]], null$Ox)                      # m x PJ
      cols <- (c - 1L) * J + seq_len(J)
      gx[, cols] <- gx[, cols] + crossprod(Tr[[k]] * Dinv[, c], Xs)
      GTx[[a]] <- gx
      GTy[[a]] <- drop(crossprod(Tr[[k]] * Dinv[, c], Ys[, c])) - drop(crossprod(Om[[a]], null$oy))
    }
    GTT <- array(0, c(m, q, q))
    for (a in seq_len(q)) for (b in a:q) {
      v <- -colSums(Om[[a]] * Om[[b]])
      if (ck$c[a] == ck$c[b])
        v <- v + colSums(Tr[[ck$k[a]]] * Tr[[ck$k[b]]] * Dinv[, ck$c[a]])
      GTT[, a, b] <- v
      GTT[, b, a] <- v
    }
    for (i in seq_len(m)) {
      if (skip[idx[i]]) next
      Hx <- t(vapply(GTx, function(g) g[i, ], numeric(P * J)))      # q x PJ
      if (q == 1L) Hx <- matrix(Hx, 1L)
      Aq <- matrix(GTT[i, , ], q, q) - Hx %*% null$XXinv %*% t(Hx)
      Aq <- (Aq + t(Aq)) / 2
      r <- vapply(GTy, `[`, 0, i) - drop(Hx %*% null$beta)
      ch <- suppressWarnings(chol(Aq, pivot = TRUE))
      if (attr(ch, "rank") < q ||
          min(diag(ch))^2 < 1e-10 * max(diag(Aq), 1e-300)) next
      b <- drop(solve(Aq, r))
      q1[idx[i]] <- null$q0 - sum(r * b)
      coef[idx[i], , ] <- matrix(b, K, P) %*% Ti                    # back to Y's traits
    }
  }
  list(q0 = null$q0, q1 = q1, coef = coef, delta = unname(delta),
       n = N * P)
}

# Dispatch for indri_gwas() and indri_gwas_haplotypes(): multi-trait, or one
# trait at a time (each with its own G, E and C_r entries).
.indri_gwas_random <- function(Y, X, G, E, terms, A, mode, delta,
                               single_trait_variance, terms_of, M, K,
                               chunk_size, skip) {
  P <- ncol(Y); traits <- colnames(Y)
  rot <- .indri_gwas_rotation(A)
  if (mode == "multi_trait") {
    s <- .indri_gwas_random_scan(Y, X, G, E, terms, rot, delta, terms_of, M, K,
                                 chunk_size, skip)
    stat <- s$n * log(s$q0 / s$q1)
    return(list(stat = stat, df = P * K, coef = s$coef, delta = s$delta))
  }
  stat <- matrix(NA_real_, M, P, dimnames = list(NULL, traits))
  coef <- array(NA_real_, c(M, K, P))
  dl <- stats::setNames(rep(NA_real_, P), traits)
  for (c in seq_len(P)) {
    if (G[c, c] <= 0)
      stop("Extra random effects in a single-trait scan need a positive genetic ",
           "variance for every trait (`", traits[c], "` has ", G[c, c], ").",
           call. = FALSE)
    one <- lapply(terms, function(tm) list(Z = tm$Z, C = tm$C[c, c, drop = FALSE]))
    d_c <- if (!is.null(delta)) delta[c] else if (single_trait_variance == "fit") 1 else NULL
    s <- .indri_gwas_random_scan(Y[, c, drop = FALSE], X, G[c, c, drop = FALSE],
                                 E[c, c, drop = FALSE], one, rot, d_c, terms_of, M, K,
                                 chunk_size, skip)
    stat[, c] <- s$n * log(s$q0 / s$q1)
    coef[, , c] <- s$coef[, , 1]
    dl[c] <- s$delta
  }
  list(stat = stat, df = K, coef = coef, delta = dl)
}

#' Generalised least squares with the scan covariance
#'
#' Fits the fixed effects `B` of `Y = X B + U + Eps` by generalised least
#' squares under the covariance of the multi-trait scans,
#' `V = G (x) A + sum_r C_r (x) Z_r Z_r' + delta E (x) I` (see
#' [indri_gwas()], sections Models and Extra random effects), with `delta`
#' fitted by REML. Useful for models with many marker terms at once, such as
#' the founder effects of all QTL estimated jointly, where the scan's
#' one-marker-at-a-time tests do not apply. It uses the same rotation and
#' low-rank update as the scans, so it never forms the `NP x NP` covariance.
#'
#' @param Y Numeric matrix, individuals by traits, with no missing values.
#' @param X Fixed-effects design (individuals by effects), full column rank.
#' @param G,E Genetic and residual covariances of the traits, in the units of
#'   `Y` (or `fit`).
#' @param A Relationship matrix of the individuals.
#' @param fit Optional `indri_fit` supplying `G`, `E` and `random_cov`.
#' @param random,random_cov Optional extra random effects, as in
#'   [indri_gwas()].
#' @param delta Optional fixed multiplier of `E`; by default fitted by REML.
#' @return A list: `coefficients` (effects by traits, in the units of `Y`),
#'   `vcov` (the covariance of `c(coefficients)`, which stacks the traits),
#'   `se` (effects by traits), `delta`, `sigma2` (the profiled scale; the
#'   covariance is `vcov` for `sigma2 = 1`, as the scans assume), `reml` (the
#'   REML criterion at `delta`, half of minus twice the restricted
#'   log-likelihood up to a constant) and `n`.
#' @seealso [indri_gwas()], [indri_gwas_haplotypes()]
#' @export
indri_gls <- function(Y, X, G = NULL, E = NULL, A, fit = NULL, random = NULL,
                      random_cov = NULL, delta = NULL) {
  Y <- .indri_gwas_matrix(Y, "Y")
  if (anyNA(Y)) stop("`Y` has missing values.", call. = FALSE)
  N <- nrow(Y); P <- ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("trait", seq_len(P))
  traits <- colnames(Y)
  X <- .indri_gwas_rows(.indri_gwas_matrix(X, "X"), rownames(Y), N, "X")
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  .indri_check_design(X)
  A <- .indri_gwas_matrix(A, "A")
  if (!is.null(fit)) {
    if (!is.null(G) || !is.null(E))
      stop("Pass either `fit` or `G` and `E`, not both.", call. = FALSE)
    G <- fit$G; E <- fit$E
    if (is.null(random_cov)) random_cov <- fit$random
  }
  if (is.null(G) || is.null(E)) stop("`indri_gls()` needs `G` and `E` (or `fit`).", call. = FALSE)
  G <- .indri_gwas_trait_matrix(G, traits, "G")
  E <- .indri_gwas_trait_matrix(E, traits, "E")
  if (!is.null(delta) && (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0))
    stop("`delta` must be one positive number.", call. = FALSE)
  terms <- if (is.null(random)) {
    # no extra term: one level with zero covariance leaves V unchanged
    list(.none = list(Z = matrix(1, N, 1L), C = G * 0))
  } else .indri_gwas_random_terms(random, random_cov, traits, N)
  rot <- .indri_gwas_rotation(A)
  L <- .indri_gwas_chol(E, "E")
  Gw <- backsolve(L, t(backsolve(L, G, upper.tri = FALSE)), upper.tri = FALSE)
  eg <- eigen((Gw + t(Gw)) / 2, symmetric = TRUE)
  d <- pmax(eg$values, 0)
  Tm <- backsolve(t(L), eg$vectors)
  Ti <- solve(Tm)
  Ys <- .indri_gwas_rotate(rot$Q, Y %*% Tm)
  Xs <- .indri_gwas_rotate(rot$Q, X)
  Zt <- lapply(terms, function(tm) .indri_gwas_rotate(rot$Q, tm$Z))
  Lc <- lapply(terms, function(tm) .indri_psd_factor(crossprod(Tm, tm$C %*% Tm)))
  D_of <- function(dl) outer(rot$lambda, d) + dl
  if (is.null(delta))
    delta <- exp(stats::optimize(function(ld) .indri_lr_null(Ys, Xs, D_of(exp(ld)), Zt, Lc)$reml,
                                 c(-12, 12), tol = 1e-9)$minimum)
  nf <- .indri_lr_null(Ys, Xs, D_of(delta), Zt, Lc)
  J <- ncol(X)
  B <- matrix(nf$beta, J, P) %*% Ti
  Kt <- kronecker(t(Ti), diag(J))                 # vec(B) = (T^-T (x) I) vec(B_t)
  V <- Kt %*% nf$XXinv %*% t(Kt)
  V <- (V + t(V)) / 2
  dimnames(B) <- list(colnames(X), traits)
  nm <- paste(rep(traits, each = J), rep(colnames(X), P), sep = ":")
  dimnames(V) <- list(nm, nm)
  n <- N * P - P * J
  list(coefficients = B, vcov = V,
       se = matrix(sqrt(pmax(diag(V), 0)), J, P, dimnames = dimnames(B)),
       delta = unname(delta), sigma2 = nf$q0 / n, reml = nf$reml, n = N)
}
