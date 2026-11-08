#' Genome-wide association scan with the assembled covariances
#'
#' Tests each SNP for association with the traits in `Y`, either jointly
#' (one multi-trait "any effect" test per SNP) or trait by trait, in a linear
#' mixed model with relationship matrix `K` or in a linear model. The
#' multi-trait mixed model uses the genetic and residual covariances `G` and
#' `E`, usually from [fit_G()], which is what makes it feasible for many
#' traits.
#'
#' @section Models:
#' For `N` individuals, `P` traits and one SNP `g`,
#' `Y = X B + g beta' + U + Eps`. Every test profiles an overall scale
#' `sigma^2` and is a likelihood-ratio test:
#'
#' * multi-trait LMM: `vec(U + Eps) ~ N(0, sigma^2 (G (x) K + delta E (x) I))`,
#'   `P` degrees of freedom. `G` and `E` set the shape; the ratio `delta` is
#'   fitted once by REML under the null (no SNP) and held fixed for every SNP.
#'   With `G` and `E` fitted to the same data it should be near 1, and is
#'   returned so that it can be checked.
#' * multi-trait LM: `sigma^2 E (x) I`, `P` degrees of freedom, with `E` by
#'   default the covariance of the residuals of `Y` on `X`.
#' * single-trait LMM: per trait, `sigma^2 (K + delta I)`, one degree of
#'   freedom, with `delta` fitted by REML under the null
#'   (`single_trait_variance = "reml"`) or taken as `E[k, k] / G[k, k]`
#'   (`"fit"`).
#' * single-trait LM: ordinary least squares.
#'
#' The covariance is diagonalised once (eigenvectors of `K`, then of `G`
#' relative to `E`), so the scan costs one `N x N` rotation of the genotypes
#' plus `O(N P)` per SNP.
#'
#' @param Y Numeric matrix, individuals by traits, with no missing values.
#' @param genotypes Numeric matrix, individuals by SNPs (typically 0/1/2).
#'   Missing values are replaced by the SNP's mean; monomorphic SNPs get `NA`.
#' @param fit An `indri_fit` from [fit_G()] supplying `G` and `E`.
#' @param G,E Genetic and residual covariance matrices of the traits, in the
#'   units of `Y`, instead of `fit`. Matched to the columns of `Y` by name
#'   when both are named. For the multi-trait LM, `E` replaces the residual
#'   covariance.
#' @param K Relationship matrix of the individuals; required for
#'   `model = "lmm"`.
#' @param X Fixed-effects design with full column rank that can represent a
#'   constant, as for [fit_G()]; defaults to an intercept. Principal
#'   components go here.
#' @param snp_info Optional data frame with one row per column of
#'   `genotypes` (for example `snp`, `chrom`, `pos`), carried into the
#'   results.
#' @param mode `"multi_trait"` for one joint test per SNP, `"single_trait"`
#'   for one test per SNP and trait.
#' @param model `"lmm"` (with `K`) or `"lm"`.
#' @param single_trait_variance For the single-trait LMM, `"reml"` fits each
#'   trait's variance ratio, `"fit"` takes it from `G` and `E`.
#' @param delta Optional fixed variance ratio for the LMMs, instead of the
#'   fitted one: a number for the multi-trait test, one value per trait for
#'   the single-trait test.
#' @param adjust Single-trait p-value adjustment for the number of traits:
#'   `"bonferroni"` multiplies by `P`, `"effective"` by
#'   [indri_effective_tests()] of the residuals after `X`; both cap at 1.
#' @param chunk_size Number of SNPs rotated at a time.
#'
#' @return An object of class `indri_gwas`: `pvalues` (a vector over SNPs,
#'   or SNPs by traits), `betas` (SNPs by traits, in the units of `Y` per
#'   unit of genotype), `delta`, `pvalues_adjust` and `n_tests` when
#'   adjusted, `snp_info`, `mode`, `model` and `call`.
#' @seealso [fit_G()], [indri_effective_tests()]
#' @references
#' Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
#' a simple, scalable approach for linear mixed models in high-dimensional
#' genetic association studies. bioRxiv, preprint.
#' \doi{10.1101/255497}
#' @export
indri_gwas <- function(Y, genotypes, fit = NULL, G = NULL, E = NULL, K = NULL,
                       X = NULL, snp_info = NULL,
                       mode = c("multi_trait", "single_trait"),
                       model = c("lmm", "lm"),
                       single_trait_variance = c("reml", "fit"),
                       delta = NULL,
                       adjust = c("none", "bonferroni", "effective"),
                       chunk_size = 1000L) {
  call <- match.call()
  mode <- match.arg(mode)
  model <- match.arg(model)
  single_trait_variance <- match.arg(single_trait_variance)
  adjust <- match.arg(adjust)
  chunk_size <- .indri_count(chunk_size, "chunk_size")

  Y <- .indri_gwas_matrix(Y, "Y")
  if (anyNA(Y))
    stop("`Y` has missing values; indri_gwas() needs complete phenotypes ",
         "(drop or impute the incomplete individuals first).", call. = FALSE)
  N <- nrow(Y)
  P <- ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("trait", seq_len(P))
  traits <- colnames(Y)

  genotypes <- .indri_gwas_matrix(genotypes, "genotypes")
  if (is.null(colnames(genotypes)))
    colnames(genotypes) <- paste0("snp", seq_len(ncol(genotypes)))
  genotypes <- .indri_gwas_rows(genotypes, rownames(Y), N, "genotypes")
  M <- ncol(genotypes)

  if (is.null(X)) {
    X <- matrix(1, N, 1L, dimnames = list(rownames(Y), "(Intercept)"))
  } else {
    X <- .indri_gwas_rows(.indri_gwas_matrix(X, "X"), rownames(Y), N, "X")
  }
  .indri_check_design(X)

  if (model == "lmm") {
    if (is.null(K))
      stop("`K` is required for model = \"lmm\".", call. = FALSE)
    K <- .indri_gwas_matrix(K, "K")
    if (is.null(rownames(K)) && !is.null(colnames(K))) rownames(K) <- colnames(K)
    if (is.null(colnames(K))) colnames(K) <- rownames(K)
    K <- .indri_gwas_rows(K, rownames(Y), N, "K")
    K <- t(.indri_gwas_rows(t(K), rownames(Y), N, "K"))
  }

  if (!is.null(snp_info)) {
    snp_info <- as.data.frame(snp_info)
    if (nrow(snp_info) != M)
      stop("`snp_info` must have one row per column of `genotypes` (", M,
           "), not ", nrow(snp_info), ".", call. = FALSE)
  } else {
    snp_info <- data.frame(snp = colnames(genotypes), stringsAsFactors = FALSE)
  }

  if (!is.null(fit)) {
    if (!is.null(G) || !is.null(E))
      stop("Pass either `fit` or `G` and `E`, not both.", call. = FALSE)
    G <- fit$G
    E <- fit$E
    if (is.null(E) && (model == "lmm" && (mode == "multi_trait" ||
                                          single_trait_variance == "fit")))
      stop("`fit` has no `E`: it was made before indri assembled the ",
           "residual covariance. Refit it with fit_G().", call. = FALSE)
  }
  needs_GE <- model == "lmm" &&
    (mode == "multi_trait" || single_trait_variance == "fit")
  if (needs_GE && (is.null(G) || is.null(E)))
    stop(if (mode == "multi_trait") "The multi-trait LMM" else
           "`single_trait_variance = \"fit\"`",
         " needs `G` and `E` (or an indri `fit`).", call. = FALSE)
  if (!is.null(G)) G <- .indri_gwas_trait_matrix(G, traits, "G")
  if (!is.null(E)) E <- .indri_gwas_trait_matrix(E, traits, "E")

  if (adjust != "none" && mode == "multi_trait")
    stop("`adjust` corrects single-trait p-values for the number of traits; ",
         "use it with mode = \"single_trait\".", call. = FALSE)
  if (!is.null(delta)) {
    if (model != "lmm")
      stop("`delta` applies to model = \"lmm\" only.", call. = FALSE)
    delta <- .indri_gwas_delta(delta, mode, traits)
  }

  genotypes <- .indri_gwas_impute(genotypes)
  monomorphic <- apply(genotypes, 2L, function(g) max(g) - min(g) < 1e-12)
  if (any(monomorphic))
    warning(sum(monomorphic), " monomorphic SNP(s) get NA results: ",
            paste(utils::head(colnames(genotypes)[monomorphic], 5L),
                  collapse = ", "),
            if (sum(monomorphic) > 5L) ", ...", call. = FALSE)

  scan <- if (mode == "multi_trait") {
    .indri_gwas_multi(Y, X, genotypes, G, E, K, model, delta, chunk_size,
                      monomorphic)
  } else {
    .indri_gwas_single(Y, X, genotypes, G, E, K, model, delta,
                       single_trait_variance, chunk_size, monomorphic)
  }

  out <- list(pvalues = scan$pvalues, betas = scan$betas, delta = scan$delta,
              pvalues_adjust = NULL, n_tests = NULL, snp_info = snp_info,
              mode = mode, model = model, n_individuals = N, traits = traits,
              call = call)
  if (adjust != "none") {
    out$n_tests <- if (adjust == "bonferroni") P else indri_effective_tests(Y, X)
    out$pvalues_adjust <- pmin(scan$pvalues * out$n_tests, 1)
  }
  structure(out, class = "indri_gwas")
}

#' Effective number of independent traits
#'
#' Galwey's (2009) effective number of tests, from the eigenvalues
#' `lambda` of the Spearman correlation matrix of the traits:
#' `(sum(sqrt(lambda)))^2 / sum(lambda)`. It is 1 for perfectly correlated
#' traits and `ncol(Y)` for independent ones, and is the factor
#' `adjust = "effective"` multiplies single-trait p-values by in
#' [indri_gwas()].
#'
#' With a design `X`, the correlation is that of the residuals of `Y` on `X`:
#' the single-trait tests see the traits after the fixed effects, and
#' covariates shared by the traits (sex, age, line) would otherwise make
#' them look more correlated than they are to the tests.
#'
#' @param Y Numeric matrix, individuals by traits.
#' @param X Optional fixed-effects design, as in [indri_gwas()].
#' @return A number between 1 and `ncol(Y)`.
#' @references Galwey, N. W. (2009). A new measure of the effective number of
#'   tests, a practical tool for comparing families of non-independent
#'   significance tests. Genetic Epidemiology 33, 559-568.
#'   \doi{10.1002/gepi.20408}
#' @export
indri_effective_tests <- function(Y, X = NULL) {
  Y <- as.matrix(Y)
  if (!is.null(X)) Y <- stats::lm.fit(as.matrix(X), Y)$residuals
  R <- stats::cor(Y, method = "spearman", use = "pairwise.complete.obs")
  lambda <- pmax(eigen(R, symmetric = TRUE, only.values = TRUE)$values, 0)
  sum(sqrt(lambda))^2 / sum(lambda)
}

# ---------------------------------------------------------------------------
# The scans. Both reduce to independent weighted regressions, one per column
# of a rotated response, with per-entry variances V[i, c] (see
# .indri_gwas_columns()).

.indri_gwas_multi <- function(Y, X, genotypes, G, E, K, model, delta,
                              chunk_size, monomorphic) {
  N <- nrow(Y); P <- ncol(Y)
  if (model == "lm" && is.null(E)) {
    resid <- stats::lm.fit(X, Y)$residuals
    E <- crossprod(resid) / (N - ncol(X))
  }
  L <- .indri_gwas_chol(E, "E")
  if (model == "lmm") {
    # E^-1/2 G E^-T/2 = W diag(d) W'; T = L^-T W whitens E and diagonalises G.
    Gw <- backsolve(L, t(backsolve(L, G, upper.tri = FALSE)), upper.tri = FALSE)
    eg <- eigen((Gw + t(Gw)) / 2, symmetric = TRUE)
    d <- pmax(eg$values, 0)
    Tm <- backsolve(t(L), eg$vectors)
    rot <- .indri_gwas_rotation(K)
  } else {
    d <- rep(0, P)
    Tm <- backsolve(t(L), diag(P))
    rot <- list(Q = NULL, lambda = rep(0, N))
  }
  Ys <- .indri_gwas_rotate(rot$Q, Y %*% Tm)
  Xs <- .indri_gwas_rotate(rot$Q, X)
  V_of <- function(dl) outer(rot$lambda, d) + if (model == "lmm") dl else 1
  if (model == "lmm" && is.null(delta))
    delta <- .indri_gwas_reml_delta(Ys, Xs, V_of)
  V <- V_of(if (model == "lmm") delta else 1)

  stat <- rep(NA_real_, ncol(genotypes))
  beta_tilde <- matrix(NA_real_, ncol(genotypes), P)
  null <- .indri_gwas_null(Ys, Xs, V)
  for (idx in .indri_gwas_chunks(ncol(genotypes), chunk_size)) {
    Gs <- .indri_gwas_rotate(rot$Q, genotypes[, idx, drop = FALSE])
    alt <- .indri_gwas_columns(null, Xs, Gs, V)
    stat[idx] <- N * P * log(sum(null$q0) / colSums(alt$q1))
    beta_tilde[idx, ] <- t(alt$beta)
  }
  stat[monomorphic] <- NA_real_
  beta_tilde[monomorphic, ] <- NA_real_
  # beta' T = beta_tilde', so beta' = beta_tilde' T^-1
  betas <- beta_tilde %*% solve(Tm)
  dimnames(betas) <- list(colnames(genotypes), colnames(Y))
  pvalues <- stats::pchisq(stat, P, lower.tail = FALSE)
  names(pvalues) <- colnames(genotypes)
  list(pvalues = pvalues, betas = betas,
       delta = if (model == "lmm") unname(delta) else NULL)
}

.indri_gwas_single <- function(Y, X, genotypes, G, E, K, model, delta,
                               single_trait_variance, chunk_size,
                               monomorphic) {
  N <- nrow(Y); P <- ncol(Y)
  rot <- if (model == "lmm") .indri_gwas_rotation(K)
         else list(Q = NULL, lambda = rep(0, N))
  Ys <- .indri_gwas_rotate(rot$Q, Y)
  Xs <- .indri_gwas_rotate(rot$Q, X)
  if (model == "lmm" && is.null(delta)) {
    delta <- if (single_trait_variance == "fit") {
      g <- diag(G)
      if (any(g <= 0))
        stop("`single_trait_variance = \"fit\"` needs a positive genetic ",
             "variance for every trait.", call. = FALSE)
      diag(E) / g
    } else {
      vapply(seq_len(P), function(k)
        .indri_gwas_reml_delta(Ys[, k, drop = FALSE], Xs,
                               function(dl) matrix(rot$lambda + dl, ncol = 1L)),
        numeric(1))
    }
    names(delta) <- colnames(Y)
  }
  V <- if (model == "lmm") outer(rot$lambda, rep(1, P)) +
         matrix(delta, N, P, byrow = TRUE)
       else matrix(1, N, P)

  pvalues <- betas <- matrix(NA_real_, ncol(genotypes), P,
                             dimnames = list(colnames(genotypes), colnames(Y)))
  null <- .indri_gwas_null(Ys, Xs, V)
  for (idx in .indri_gwas_chunks(ncol(genotypes), chunk_size)) {
    Gs <- .indri_gwas_rotate(rot$Q, genotypes[, idx, drop = FALSE])
    alt <- .indri_gwas_columns(null, Xs, Gs, V)
    stat <- N * log(null$q0 / alt$q1)
    pvalues[idx, ] <- t(stats::pchisq(stat, 1, lower.tail = FALSE))
    betas[idx, ] <- t(alt$beta)
  }
  pvalues[monomorphic, ] <- NA_real_
  betas[monomorphic, ] <- NA_real_
  list(pvalues = pvalues, betas = betas,
       delta = if (model == "lmm") delta else NULL)
}

# Weighted least squares of every column of Ys on Xs under the null, kept for
# the SNP tests: per column c, weights 1 / V[, c], the residuals r0 and their
# weighted sum of squares q0.
.indri_gwas_null <- function(Ys, Xs, V) {
  C <- ncol(Ys)
  W <- 1 / V
  r0 <- matrix(0, nrow(Ys), C)
  XtWX_inv <- vector("list", C)
  for (c in seq_len(C)) {
    A <- crossprod(Xs, Xs * W[, c])
    XtWX_inv[[c]] <- chol2inv(chol(A))
    b <- XtWX_inv[[c]] %*% crossprod(Xs, Ys[, c] * W[, c])
    r0[, c] <- Ys[, c] - Xs %*% b
  }
  list(W = W, r0 = r0, q0 = colSums(W * r0^2), XtWX_inv = XtWX_inv)
}

# Add each SNP to every column's regression (Frisch-Waugh: residualise the
# rotated SNP on Xs with the column's weights). Returns, per column (rows)
# and SNP (columns), the SNP effect and the residual sum of squares.
.indri_gwas_columns <- function(null, Xs, Gs, V) {
  C <- ncol(null$r0); m <- ncol(Gs)
  beta <- q1 <- matrix(NA_real_, C, m)
  for (c in seq_len(C)) {
    w <- null$W[, c]
    Gt <- Gs - Xs %*% (null$XtWX_inv[[c]] %*% crossprod(Xs, Gs * w))
    num <- colSums(Gt * (w * null$r0[, c]))
    den <- colSums(Gt^2 * w)
    beta[c, ] <- num / den
    q1[c, ] <- null$q0[c] - num^2 / den
  }
  list(beta = beta, q1 = q1)
}

# REML for the variance ratio, sigma^2 profiled: V_of(delta) gives the
# per-entry variances of the rotated response for that ratio.
.indri_gwas_reml_delta <- function(Ys, Xs, V_of) {
  n <- length(Ys) - ncol(Ys) * ncol(Xs)
  criterion <- function(log_delta) {
    V <- V_of(exp(log_delta))
    W <- 1 / V
    q <- 0
    logdet_A <- 0
    for (c in seq_len(ncol(Ys))) {
      A <- crossprod(Xs, Xs * W[, c])
      R <- chol(A)
      b <- backsolve(R, backsolve(R, crossprod(Xs, Ys[, c] * W[, c]),
                                  transpose = TRUE))
      r <- Ys[, c] - Xs %*% b
      q <- q + sum(W[, c] * r^2)
      logdet_A <- logdet_A + 2 * sum(log(diag(R)))
    }
    0.5 * (sum(log(V)) + logdet_A + n * log(q / n))
  }
  exp(stats::optimize(criterion, c(-12, 12), tol = 1e-9)$minimum)
}

.indri_gwas_rotation <- function(K) {
  e <- eigen((K + t(K)) / 2, symmetric = TRUE)
  if (min(e$values) < -1e-8 * max(abs(e$values)))
    stop("`K` must be positive semi-definite.", call. = FALSE)
  list(Q = e$vectors, lambda = pmax(e$values, 0))
}

.indri_gwas_rotate <- function(Q, x) if (is.null(Q)) unname(x) else crossprod(Q, x)

.indri_gwas_chunks <- function(M, chunk_size) {
  split(seq_len(M), ceiling(seq_len(M) / chunk_size))
}

# ---------------------------------------------------------------------------
# Input handling.

.indri_gwas_matrix <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x))
    stop("`", name, "` must be numeric.", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

# Put the rows of `x` in the order of the samples of Y: by name when both are
# named, otherwise by position.
.indri_gwas_rows <- function(x, ids, N, name) {
  if (!is.null(ids) && !is.null(rownames(x))) {
    missing <- setdiff(ids, rownames(x))
    if (length(missing))
      stop("`", name, "` has no rows for ", length(missing),
           " sample(s) of `Y`: ",
           paste(utils::head(missing, 5L), collapse = ", "),
           if (length(missing) > 5L) ", ...", call. = FALSE)
    return(x[ids, , drop = FALSE])
  }
  if (nrow(x) != N)
    stop("`", name, "` has ", nrow(x), " rows but `Y` has ", N,
         "; name the rows to match samples by name.", call. = FALSE)
  x
}

.indri_gwas_trait_matrix <- function(x, traits, name) {
  x <- .indri_gwas_matrix(x, name)
  P <- length(traits)
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    if (!setequal(rownames(x), traits) || !setequal(colnames(x), traits))
      stop("`", name, "` does not describe the traits of `Y` (",
           paste(utils::head(traits, 5L), collapse = ", "),
           if (P > 5L) ", ...", ").", call. = FALSE)
    x <- x[traits, traits, drop = FALSE]
  } else if (nrow(x) != P || ncol(x) != P) {
    stop("`", name, "` must be ", P, " x ", P, ", one row and column per ",
         "trait of `Y`.", call. = FALSE)
  }
  if (any(!is.finite(x)) || max(abs(x - t(x))) > 1e-8 * max(abs(x), 1))
    stop("`", name, "` must be a finite symmetric matrix.", call. = FALSE)
  x <- (x + t(x)) / 2
  dimnames(x) <- list(traits, traits)
  x
}

.indri_gwas_chol <- function(x, name) {
  L <- tryCatch(t(chol(x)), error = function(e) NULL)
  if (is.null(L))
    stop("`", name, "` must be positive definite.", call. = FALSE)
  L
}

.indri_gwas_delta <- function(delta, mode, traits) {
  if (!is.numeric(delta) || any(!is.finite(delta)) || any(delta <= 0))
    stop("`delta` must be positive and finite.", call. = FALSE)
  if (mode == "multi_trait") {
    if (length(delta) != 1L)
      stop("`delta` must be a single number for the multi-trait test.",
           call. = FALSE)
    return(unname(delta))
  }
  if (length(delta) != length(traits))
    stop("`delta` must have one value per trait (", length(traits), ").",
         call. = FALSE)
  if (!is.null(names(delta))) {
    if (!setequal(names(delta), traits))
      stop("The names of `delta` must be the traits of `Y`.", call. = FALSE)
    return(delta[traits])
  }
  stats::setNames(delta, traits)
}

.indri_gwas_impute <- function(genotypes) {
  if (!anyNA(genotypes)) return(genotypes)
  for (j in which(colSums(is.na(genotypes)) > 0)) {
    g <- genotypes[, j]
    g[is.na(g)] <- if (all(is.na(g))) 0 else mean(g, na.rm = TRUE)
    genotypes[, j] <- g
  }
  genotypes
}

# ---------------------------------------------------------------------------
# Methods.

#' @export
print.indri_gwas <- function(x, ...) {
  what <- if (x$mode == "multi_trait") "multi-trait" else "single-trait"
  cat("indri association scan:", what, toupper(x$model), "\n")
  M <- NROW(x$pvalues)
  cat(sprintf("  %d SNPs, %d traits, %d individuals\n", M, length(x$traits),
              x$n_individuals))
  if (!is.null(x$delta)) {
    if (length(x$delta) == 1L) cat(sprintf("  null variance ratio delta: %.3g\n",
                                           x$delta))
    else cat(sprintf("  null variance ratio delta: %.3g to %.3g over traits\n",
                     min(x$delta), max(x$delta)))
  }
  if (!is.null(x$n_tests))
    cat(sprintf("  p-values adjusted for %.3g tests\n", x$n_tests))
  if (all(is.na(x$pvalues))) {
    cat("  no SNP could be tested\n")
  } else {
    cat(sprintf("  smallest p-value: %.3g\n", min(x$pvalues, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
as.data.frame.indri_gwas <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  info <- x$snp_info
  if (!"snp" %in% names(info))
    info <- cbind(snp = rownames(x$betas), info, stringsAsFactors = FALSE)
  if (x$mode == "multi_trait") {
    betas <- as.data.frame(x$betas)
    names(betas) <- paste0("beta_", colnames(x$betas))
    out <- cbind(info, pvalue = unname(x$pvalues), betas)
  } else {
    out <- do.call(rbind, lapply(colnames(x$pvalues), function(k) {
      part <- cbind(info, trait = k, pvalue = unname(x$pvalues[, k]),
                    beta = unname(x$betas[, k]), stringsAsFactors = FALSE)
      if (!is.null(x$pvalues_adjust))
        part$pvalue_adjust <- unname(x$pvalues_adjust[, k])
      part
    }))
  }
  rownames(out) <- row.names
  out
}
