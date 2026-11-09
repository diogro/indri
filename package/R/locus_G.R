#' Genetic covariance matrices from locus effects
#'
#' The additive (`G`) and dominance (`D`) covariance matrices of a set of
#' traits implied by the additive and dominance effects of individual loci,
#' their allele frequencies and the linkage disequilibrium between them
#' (Melo, Marroig and Wolf 2019, Eqs. 7-14). Each locus contributes through
#' its vector of pleiotropic effects, so `G` and `D` show how the loci found
#' by [indri_gwas()], or by any other mapping method, build up the genetic
#' covariance between the traits.
#'
#' @details For locus `j` with allele frequency `p_j` (of the allele the
#' 0/1/2 dosage counts; `q_j = 1 - p_j`), additive effects `a_j` and
#' dominance effects `d_j` (vectors over the traits), the average effect is
#' `alpha_j = a_j + (q_j - p_j) d_j` (Eq. 7), and
#'
#' * `G = sum_j 2 p_j q_j alpha_j alpha_j' + sum_{j != x} 2 lambda_jx
#'   alpha_j alpha_x'` (Eqs. 9-10),
#' * `D = sum_j (2 p_j q_j)^2 d_j d_j' + sum_{j != x} 4 lambda_jx^2 d_j d_x'`
#'   (Eqs. 11-12),
#'
#' where `lambda_jx` is the linkage disequilibrium between loci `j` and `x`.
#' From an unphased genotype matrix it is the composite estimate
#' `lambda_jx = Cov(g_j, g_x) / 2` of the dosages.
#'
#' @section Conventions for effects from other software:
#' The formulas assume one convention for the effects; estimates from other
#' packages have to be put on it, or `G` and `D` are wrong without any error.
#'
#' * **Effect allele.** `a` is the effect of one copy of the allele that the
#'   0/1/2 genotypes count (half the difference between the two
#'   homozygotes), and `p` is that allele's frequency. GWAS software reports
#'   each beta for a named allele (for example `A1` in PLINK 2's `--glm`
#'   output, `ALLELE1` in regenie, `allele1` in GEMMA): the genotypes must
#'   count that allele, or that locus's `a` must be negated. Check this per
#'   locus; tools differ in which allele they report.
#' * **Dominance coding.** `d` is the deviation of the heterozygote from the
#'   midpoint of the homozygotes, which is what a regression on the dosage
#'   plus a heterozygote indicator (1 for genotype 1, 0 otherwise) estimates
#'   (as [indri_gwas()] does, and as the LL = +1, LS = 0, SS = -1 coding of
#'   [lgsm_genotypes] with its `D` column). Orthogonal or frequency-centred
#'   dominance codings (Cockerham's, for instance) estimate other quantities
#'   and must be converted first. Most GWAS software estimates no dominance:
#'   leave `d = NULL`.
#' * **Units.** Effects in the units of the traits give `G` and `D` in those
#'   units squared. Standardised betas give matrices on the standardised
#'   scale.
#' * **Shape.** Effects are loci by traits (one row per locus); a
#'   matrix with one row per trait needs `traits_in = "rows"`. A data
#'   frame may carry one non-numeric column of locus names.
#'
#' Loci are matched to the columns of `genotypes`, and to `d`, `p` and
#' `lambda`, by name whenever both sides are named. Matching to `genotypes`
#' by position warns, as does matching `d` by position when only one side
#' is named.
#'
#' @param a Additive effects, loci by traits (a matrix, a numeric data
#'   frame, or a data frame with one column of locus names; a vector for one
#'   trait), or an `indri_gwas` scan (whose `a`, and `d` with
#'   `dominance = TRUE`, are used). Loci without an estimate are dropped.
#' @param d Dominance effects, in the same layout as `a`; `NULL` for none.
#' @param genotypes Individuals by loci, 0/1/2 dosages of the allele the
#'   effects refer to, for the allele frequencies and the linkage
#'   disequilibrium. May have more loci than `a`.
#' @param p Allele frequencies, one per locus, instead of `genotypes`.
#' @param lambda Linkage disequilibrium between the loci (loci by loci; the
#'   diagonal is ignored), instead of the estimate from `genotypes`.
#' @param ld If `FALSE`, ignore linkage disequilibrium (only the per-locus
#'   sums).
#' @param traits_in Where the traits are in `a` and `d`: `"columns"` (the
#'   default: one row per locus, one column per trait) or `"rows"`.
#' @return A list with `G` and `D` (traits by traits), `alpha` (the average
#'   effects, loci by traits), `p` and `loci` (the loci used).
#' @references Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic
#'   perspective on multivariate variation, pleiotropy, and evolution.
#'   Journal of Heredity 110, 479-493. \doi{10.1093/jhered/esz011}
#' @seealso [indri_gwas()]
#' @examples
#' # one locus with effects on three traits
#' indri_locus_G(matrix(c(0.5, -0.2, 0.1), 1), d = matrix(c(0.3, 0, 0), 1),
#'               p = 0.4)$G
#' @export
indri_locus_G <- function(a, d = NULL, genotypes = NULL, p = NULL,
                          lambda = NULL, ld = TRUE,
                          traits_in = c("columns", "rows")) {
  traits_in <- match.arg(traits_in)
  if (inherits(a, "indri_gwas")) {
    scan <- a
    a <- scan$a
    d <- if (isTRUE(scan$dominance)) scan$d else NULL
  }
  a <- .indri_effects(a, "a", traits_in)
  d <- if (is.null(d)) NULL else .indri_effects(d, "d", traits_in)
  a_named <- !is.null(rownames(a))
  M <- nrow(a)

  if (!is.null(genotypes)) {
    genotypes <- .indri_locus_genotypes(genotypes)
    .indri_check_transposed(a, genotypes)
    genotypes <- .indri_match_loci(genotypes, if (a_named) rownames(a), M,
                                   "genotypes", columns = TRUE)
  } else if (is.null(p)) {
    stop("Give `genotypes` (0/1/2 dosages) or the allele frequencies `p`.",
         call. = FALSE)
  }
  if (!a_named) rownames(a) <- paste0("locus", seq_len(M))
  loci <- rownames(a)

  if (!is.null(d)) {
    if (ncol(d) != ncol(a))
      stop("`d` must have the traits of `a` (", ncol(a), ").", call. = FALSE)
    d <- .indri_match_loci(d, if (a_named) loci, M, "d")
    rownames(d) <- loci
  }

  # loci with no estimate (monomorphic, or dominance not estimable) go
  keep <- stats::complete.cases(a)
  if (!is.null(d)) keep <- keep & stats::complete.cases(d)
  if (!all(keep)) {
    message(sum(!keep), ngettext(sum(!keep), " locus", " loci"),
            " without estimated effects dropped: ",
            paste(utils::head(loci[!keep], 5L), collapse = ", "))
  }

  if (is.null(p)) {
    p <- colMeans(genotypes, na.rm = TRUE) / 2
  } else {
    p <- .indri_match_vector(p, if (a_named) loci, M, "p")
  }
  if (any(!is.finite(p)) || any(p < 0 | p > 1))
    stop("Allele frequencies `p` must lie in [0, 1].", call. = FALSE)

  if (ld && is.null(lambda) && !is.null(genotypes)) {
    lambda <- stats::cov(genotypes, use = "pairwise.complete.obs") / 2
  }
  if (!ld || is.null(lambda)) {
    lambda <- matrix(0, M, M)
  } else {
    lambda <- as.matrix(lambda)
    if (!identical(dim(lambda), c(M, M)))
      stop("`lambda` must be ", M, " x ", M, ", one row and column per locus.",
           call. = FALSE)
    if (a_named && !is.null(rownames(lambda)) && !is.null(colnames(lambda))) {
      lambda <- .indri_match_loci(lambda, loci, M, "lambda")
      lambda <- .indri_match_loci(lambda, loci, M, "lambda", columns = TRUE)
    }
  }

  a <- a[keep, , drop = FALSE]
  if (!is.null(d)) d <- d[keep, , drop = FALSE]
  p <- unname(p[keep])
  lambda <- unname(lambda[keep, keep, drop = FALSE])
  q <- 1 - p
  het <- 2 * p * q

  alpha <- if (is.null(d)) a else a + (q - p) * d
  S <- 2 * lambda
  diag(S) <- het
  G <- crossprod(alpha, S %*% alpha)
  D <- if (is.null(d)) {
    matrix(0, ncol(a), ncol(a))
  } else {
    Delta <- 4 * lambda^2
    diag(Delta) <- het^2
    crossprod(d, Delta %*% d)
  }
  traits <- colnames(a)
  G <- (G + t(G)) / 2
  D <- (D + t(D)) / 2
  dimnames(G) <- dimnames(D) <- list(traits, traits)
  list(G = G, D = D, alpha = alpha,
       p = stats::setNames(p, rownames(a)), loci = rownames(a))
}

# Effects as a numeric matrix, loci by traits. A data frame may carry one
# non-numeric column, taken as the locus names.
.indri_effects <- function(x, name, traits_in) {
  ids <- NULL
  if (is.data.frame(x)) {
    id <- !vapply(x, is.numeric, logical(1))
    if (sum(id) > 1L)
      stop("`", name, "` has more than one non-numeric column (",
           paste(names(x)[id], collapse = ", "), "); keep at most one, of ",
           "locus names, and numeric effects.", call. = FALSE)
    if (any(id)) ids <- as.character(x[[which(id)]])
    x <- as.matrix(x[!id])
    if (!is.null(ids)) rownames(x) <- ids
  }
  if (is.null(dim(x)))
    x <- matrix(x, ncol = 1L, dimnames = list(names(x), NULL))
  x <- as.matrix(x)
  if (!is.numeric(x))
    stop("`", name, "` must hold numeric effects, loci by traits.",
         call. = FALSE)
  if (traits_in == "rows") x <- t(x)
  storage.mode(x) <- "double"
  x
}

.indri_locus_genotypes <- function(genotypes) {
  genotypes <- as.matrix(genotypes)
  if (!is.numeric(genotypes))
    stop("`genotypes` must be numeric 0/1/2 dosages.", call. = FALSE)
  range <- suppressWarnings(range(genotypes, na.rm = TRUE))
  if (range[1] < 0 || range[2] > 2)
    stop("`genotypes` must be 0/1/2 dosages of the allele the effects refer ",
         "to (values run from ", signif(range[1], 3), " to ",
         signif(range[2], 3), "); a -1/0/1 coding needs + 1.", call. = FALSE)
  genotypes
}

# Effects given traits by loci: the column names of `a` are genotype columns
# and its row names are not, or (without names to go by) the trait count
# matches the genotype columns and the locus count does not.
.indri_check_transposed <- function(a, genotypes) {
  g_names <- colnames(genotypes)
  by_name <- !is.null(g_names) && !is.null(colnames(a)) &&
    all(colnames(a) %in% g_names) && !any(rownames(a) %in% g_names)
  by_size <- is.null(rownames(a)) && ncol(a) == ncol(genotypes) &&
    nrow(a) != ncol(genotypes)
  if (by_name || by_size)
    stop("The effects look transposed: they must be loci by traits (one row ",
         "per locus). Use t(), or `traits_in = \"rows\"`.",
         call. = FALSE)
}

# Rows (or, with `columns`, columns) of `x` put in the order of `loci`: by
# name when both are named (every locus must be there), by position
# otherwise. Position matching warns for genotypes, and for `d` when only
# one side is named.
.indri_match_loci <- function(x, loci, M, name, columns = FALSE) {
  have <- if (columns) colnames(x) else rownames(x)
  n <- if (columns) ncol(x) else nrow(x)
  what <- if (columns) "column" else "row"
  if (!is.null(loci) && !is.null(have)) {
    missing <- setdiff(loci, have)
    if (length(missing))
      stop("`", name, "` has no ", what, " for ", length(missing),
           ngettext(length(missing), " locus", " loci"), " of the effects: ",
           paste(utils::head(missing, 5L), collapse = ", "),
           if (length(missing) > 5L) ", ...", ". Loci are matched by name.",
           call. = FALSE)
    return(if (columns) x[, loci, drop = FALSE] else x[loci, , drop = FALSE])
  }
  if (n != M)
    stop("`", name, "` must have one ", what, " per locus of the effects (",
         M, "), not ", n, ".", call. = FALSE)
  one_named <- xor(is.null(loci), is.null(have))
  if (name == "genotypes" || (name == "d" && one_named))
    warning("The loci of the effects and of `", name, "` were matched by ",
            "position; name both (row names of the effects, ",
            if (columns) "column" else "row", " names of `", name,
            "`) to match them by name.", call. = FALSE)
  x
}

.indri_match_vector <- function(v, loci, M, name) {
  if (!is.null(loci) && !is.null(names(v))) {
    missing <- setdiff(loci, names(v))
    if (length(missing))
      stop("`", name, "` has no value for ", length(missing),
           ngettext(length(missing), " locus", " loci"), " of the effects: ",
           paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    return(v[loci])
  }
  if (length(v) != M)
    stop("`", name, "` must have one value per locus of the effects (", M,
         "), not ", length(v), ".", call. = FALSE)
  v
}
