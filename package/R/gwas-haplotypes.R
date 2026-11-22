#' Association scan with founder haplotypes
#'
#' Tests each position for an association between the traits in `Y` and the
#' founder haplotypes of the individuals (multi-parent crosses: advanced
#' intercrosses, the Diversity Outbred, MAGIC, heterogeneous stocks), in the
#' same linear mixed models as [indri_gwas()]. At each position the dosages of
#' the founders (expected copies of each founder's haplotype, for example
#' from `qtl2::genoprob_to_alleleprob()` times 2) enter together, so the test
#' has `F - 1` degrees of freedom per trait for `F` founders.
#'
#' @param Y Numeric matrix, individuals by traits, with no missing values.
#' @param haplotypes Numeric array, individuals by founders by positions, of
#'   founder dosages (each individual's dosages sum to 2 at every position).
#'   The founder and position names are carried into the results.
#' @param fit,G,E,A,X,mode,model,single_trait_variance,delta As for
#'   [indri_gwas()].
#' @param snp_info Optional data frame with one row per position (for example
#'   `snp`, `chrom`, `pos`), carried into the results.
#' @param baseline Founder whose dosage is left out of the model (dosages sum
#'   to 2, so one is redundant with the intercept); by default the founder with
#'   the largest mean dosage. Effects are reported for every founder, centred
#'   on their mean, so the choice does not change them or the test.
#' @param chunk_size Number of positions rotated at a time.
#'
#' @return An object of class `indri_gwas`, so that `plot()`, `summary()` and
#'   `as.data.frame()` work on it: `pvalues` (a vector over positions, or
#'   positions by traits), `effects` (positions by founders by traits: founder
#'   effects on each trait, in the units of `Y` per copy of the founder's
#'   haplotype, centred on the mean founder), `a` (positions by traits: the
#'   spread of the founder effects, largest minus smallest), `df` (degrees of
#'   freedom of each test), `delta`, `snp_info`, `mode`, `model`, `founders`,
#'   `baseline`, `traits`, `n_individuals`, `haplotypes = TRUE` and `call`.
#'
#' @details
#' The model is that of [indri_gwas()] with the single dosage replaced by the
#' `F - 1` founder dosages: for the multi-trait LMM,
#' `vec(Y - X B - H Beta) ~ N(0, sigma^2 (G (x) A + delta E (x) I))`, and a
#' likelihood-ratio test of `Beta = 0` with `P (F - 1)` degrees of freedom.
#' The variance ratio `delta` is fitted once under the null by REML and held
#' fixed, as in [indri_gwas()]. A position where the founder dosages are
#' collinear in some rotated column (a founder absent, or two founders never
#' told apart) gets `NA`.
#'
#' @seealso [indri_gwas()] for SNP scans, [indri_locus_G()].
#' @export
indri_gwas_haplotypes <- function(Y, haplotypes, fit = NULL, G = NULL, E = NULL, A = NULL,
                                  X = NULL, snp_info = NULL,
                                  mode = c("multi_trait", "single_trait"),
                                  model = c("lmm", "lm"),
                                  single_trait_variance = c("reml", "fit"),
                                  delta = NULL, baseline = NULL, chunk_size = 500L) {
  call <- match.call()
  mode <- match.arg(mode)
  model <- match.arg(model)
  single_trait_variance <- match.arg(single_trait_variance)
  chunk_size <- .indri_count(chunk_size, "chunk_size")

  Y <- .indri_gwas_matrix(Y, "Y")
  if (anyNA(Y))
    stop("`Y` has missing values; drop or impute the incomplete individuals first.",
         call. = FALSE)
  N <- nrow(Y); P <- ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("trait", seq_len(P))
  traits <- colnames(Y)

  if (!is.array(haplotypes) || length(dim(haplotypes)) != 3L || !is.numeric(haplotypes))
    stop("`haplotypes` must be a numeric array, individuals by founders by positions.",
         call. = FALSE)
  if (dim(haplotypes)[1] != N)
    stop("`haplotypes` must have one row per row of `Y`.", call. = FALSE)
  if (!is.null(rownames(Y)) && !is.null(dimnames(haplotypes)[[1]]) &&
      !identical(dimnames(haplotypes)[[1]], rownames(Y)))
    stop("The individuals of `haplotypes` are not in the order of the rows of `Y`.",
         call. = FALSE)
  if (anyNA(haplotypes)) stop("`haplotypes` has missing values.", call. = FALSE)
  Fn <- dim(haplotypes)[2]; M <- dim(haplotypes)[3]
  if (Fn < 2L) stop("`haplotypes` needs at least two founders.", call. = FALSE)
  founders <- dimnames(haplotypes)[[2]] %||% paste0("founder", seq_len(Fn))
  positions <- dimnames(haplotypes)[[3]] %||% paste0("pos", seq_len(M))
  storage.mode(haplotypes) <- "double"

  if (is.null(baseline)) {
    baseline <- founders[which.max(apply(haplotypes, 2, mean))]
  } else if (!baseline %in% founders) {
    stop("`baseline` must be one of the founders.", call. = FALSE)
  }
  b0 <- match(baseline, founders)
  kept <- setdiff(seq_len(Fn), b0)

  if (is.null(X)) {
    X <- matrix(1, N, 1L, dimnames = list(rownames(Y), "(Intercept)"))
  } else {
    X <- .indri_gwas_rows(.indri_gwas_matrix(X, "X"), rownames(Y), N, "X")
  }
  .indri_check_design(X)
  if (model == "lmm") {
    if (is.null(A)) stop("`A` is required for model = \"lmm\".", call. = FALSE)
    A <- .indri_gwas_matrix(A, "A")
    if (is.null(rownames(A)) && !is.null(colnames(A))) rownames(A) <- colnames(A)
    if (is.null(colnames(A))) colnames(A) <- rownames(A)
    A <- .indri_gwas_rows(A, rownames(Y), N, "A")
    A <- t(.indri_gwas_rows(t(A), rownames(Y), N, "A"))
  }
  if (!is.null(snp_info)) {
    snp_info <- as.data.frame(snp_info)
    if (nrow(snp_info) != M)
      stop("`snp_info` must have one row per position (", M, ").", call. = FALSE)
    if (is.null(snp_info$snp)) snp_info$snp <- positions
  } else {
    snp_info <- data.frame(snp = positions, stringsAsFactors = FALSE)
  }
  if (!is.null(fit)) {
    if (!is.null(G) || !is.null(E))
      stop("Pass either `fit` or `G` and `E`, not both.", call. = FALSE)
    G <- fit$G; E <- fit$E
  }
  if (!is.null(G)) G <- G[traits, traits, drop = FALSE]
  if (!is.null(E)) E <- E[traits, traits, drop = FALSE]
  if (mode == "multi_trait" && model == "lmm" && (is.null(G) || is.null(E)))
    stop("The multi-trait LMM needs `fit`, or `G` and `E`.", call. = FALSE)
  if (!is.null(delta)) delta <- .indri_gwas_delta(delta, mode, traits)

  # ---- rotation and null model (as in indri_gwas) ----------------------------
  if (mode == "multi_trait") {
    if (model == "lm" && is.null(E)) {
      resid <- stats::lm.fit(X, Y)$residuals
      E <- crossprod(resid) / (N - ncol(X))
    }
    L <- .indri_gwas_chol(E, "E")
    if (model == "lmm") {
      Gw <- backsolve(L, t(backsolve(L, G, upper.tri = FALSE)), upper.tri = FALSE)
      eg <- eigen((Gw + t(Gw)) / 2, symmetric = TRUE)
      d <- pmax(eg$values, 0)
      Tm <- backsolve(t(L), eg$vectors)
      rot <- .indri_gwas_rotation(A)
    } else {
      d <- rep(0, P); Tm <- backsolve(t(L), diag(P))
      rot <- list(Q = NULL, lambda = rep(0, N))
    }
    Ys <- .indri_gwas_rotate(rot$Q, Y %*% Tm)
    Xs <- .indri_gwas_rotate(rot$Q, X)
    V_of <- function(dl) outer(rot$lambda, d) + if (model == "lmm") dl else 1
    if (model == "lmm" && is.null(delta)) delta <- .indri_gwas_reml_delta(Ys, Xs, V_of)
    V <- V_of(if (model == "lmm") delta else 1)
  } else {
    rot <- if (model == "lmm") .indri_gwas_rotation(A) else list(Q = NULL, lambda = rep(0, N))
    Ys <- .indri_gwas_rotate(rot$Q, Y)
    Xs <- .indri_gwas_rotate(rot$Q, X)
    if (model == "lmm" && is.null(delta)) {
      delta <- if (single_trait_variance == "fit") {
        if (is.null(G) || is.null(E) || any(diag(G) <= 0))
          stop("`single_trait_variance = \"fit\"` needs `G` and `E` with positive ",
               "genetic variances.", call. = FALSE)
        diag(E) / diag(G)
      } else {
        vapply(seq_len(P), function(k)
          .indri_gwas_reml_delta(Ys[, k, drop = FALSE], Xs,
                                 function(dl) matrix(rot$lambda + dl, ncol = 1L)), numeric(1))
      }
      names(delta) <- traits
    }
    V <- if (model == "lmm") outer(rot$lambda, rep(1, P)) + matrix(delta, N, P, byrow = TRUE)
         else matrix(1, N, P)
  }
  null <- .indri_gwas_null(Ys, Xs, V)

  # ---- scan ---------------------------------------------------------------------
  K <- length(kept)
  df <- if (mode == "multi_trait") P * K else K
  pvalues <- if (mode == "multi_trait") rep(NA_real_, M) else matrix(NA_real_, M, P)
  beta <- array(NA_real_, c(M, K, P))                 # baseline-coded, rotated columns
  for (idx in .indri_gwas_chunks(M, chunk_size)) {
    terms <- lapply(kept, function(f)
      .indri_gwas_rotate(rot$Q, matrix(haplotypes[, f, idx], N, length(idx))))
    fits <- .indri_gwas_nested(null, Xs, terms, V, list(all = seq_len(K)))$all
    if (mode == "multi_trait") {
      ok <- colSums(!fits$ok) == 0L
      q1 <- colSums(fits$q)
      stat <- N * P * log(sum(null$q0) / q1)
      pvalues[idx] <- ifelse(ok, stats::pchisq(stat, df, lower.tail = FALSE), NA_real_)
      for (c in seq_len(P)) beta[idx, , c] <- fits$beta[c, , ]
    } else {
      stat <- N * log(null$q0 / fits$q)               # columns x positions
      p <- stats::pchisq(stat, df, lower.tail = FALSE)
      p[!fits$ok] <- NA_real_
      pvalues[idx, ] <- t(p)
      for (c in seq_len(P)) beta[idx, , c] <- fits$beta[c, , ]
    }
  }

  # effects on the traits: back-transform the multi-trait columns, add the
  # baseline (0) and centre on the mean founder
  effects <- array(0, c(M, Fn, P), dimnames = list(positions, founders, traits))
  for (k in seq_len(K)) {
    bk <- matrix(beta[, k, ], M, P)
    effects[, kept[k], ] <- if (mode == "multi_trait") bk %*% solve(Tm) else bk
  }
  effects <- sweep(effects, c(1, 3), apply(effects, c(1, 3), mean))
  bad <- if (mode == "multi_trait") is.na(pvalues) else apply(is.na(pvalues), 1, all)
  effects[bad, , ] <- NA_real_
  spread <- apply(effects, c(1, 3), function(v) diff(range(v)))
  dimnames(spread) <- list(positions, traits)
  if (mode == "multi_trait") names(pvalues) <- positions
  else dimnames(pvalues) <- list(positions, traits)

  structure(list(pvalues = pvalues, a = spread, effects = effects, df = df,
                 delta = if (model == "lmm") unname(delta) else NULL,
                 pvalues_adjust = NULL, n_tests = NULL, snp_info = snp_info,
                 dominance = FALSE, gxe_terms = NULL, background = "kronecker",
                 mode = mode, model = model, n_individuals = N, traits = traits,
                 founders = founders, baseline = baseline, haplotypes = TRUE,
                 call = call),
            class = "indri_gwas")
}
