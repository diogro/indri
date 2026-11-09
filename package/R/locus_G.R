#' Genetic covariance matrices from locus effects
#'
#' The additive (`G`) and dominance (`D`) covariance matrices of a set of
#' traits implied by the additive and dominance effects of individual loci,
#' their allele frequencies and the linkage disequilibrium between them
#' (Melo, Marroig and Wolf 2019, Eqs. 7-14). Each locus contributes through
#' its vector of pleiotropic effects, so `G` and `D` show how the loci found
#' by [indri_gwas()] build up the genetic covariance between the traits.
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
#' `lambda_jx = Cov(g_j, g_x) / 2` of the dosages. The covariance between full
#' sibs is `G / 2 + D / 4`.
#'
#' The effects are on the scale of [indri_gwas()]'s `a` and `d`: `a` is the
#' effect of one copy of the counted allele (half the difference between the
#' homozygotes) and `d` the deviation of the heterozygote from their
#' midpoint, as in the LL = +1, LS = 0, SS = -1 coding of [lgsm_genotypes].
#'
#' @param a Additive effects, loci by traits, or an `indri_gwas` scan (whose
#'   `a`, and `d` with `dominance = TRUE`, are used; loci without an
#'   additive or dominance estimate are dropped).
#' @param d Dominance effects, loci by traits; `NULL` for none.
#' @param genotypes Individuals by loci, coded 0/1/2, for the allele
#'   frequencies and the linkage disequilibrium; matched to the loci of `a`
#'   by column name when both are named.
#' @param p Allele frequencies, one per locus, instead of `genotypes`.
#' @param lambda Linkage disequilibrium between the loci (loci by loci; the
#'   diagonal is ignored), instead of the estimate from `genotypes`.
#' @param ld If `FALSE`, ignore linkage disequilibrium (only the per-locus
#'   sums).
#' @return A list with `G`, `D` and `full_sib` (traits by traits), `alpha`
#'   (the average effects, loci by traits), `p` and `loci` (the loci used).
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
                          lambda = NULL, ld = TRUE) {
  if (inherits(a, "indri_gwas")) {
    scan <- a
    a <- scan$a
    d <- if (isTRUE(scan$dominance)) scan$d else NULL
  }
  a <- as.matrix(a)
  M <- nrow(a)
  if (is.null(rownames(a))) rownames(a) <- paste0("locus", seq_len(M))
  if (!is.null(d)) {
    d <- as.matrix(d)
    if (!identical(dim(d), dim(a)))
      stop("`d` must have the dimensions of `a` (", M, " x ", ncol(a), ").",
           call. = FALSE)
  }

  if (!is.null(genotypes)) {
    genotypes <- as.matrix(genotypes)
    if (!is.null(colnames(genotypes)) && all(rownames(a) %in% colnames(genotypes))) {
      genotypes <- genotypes[, rownames(a), drop = FALSE]
    } else if (ncol(genotypes) != M) {
      stop("`genotypes` must have one column per locus of `a` (", M, ").",
           call. = FALSE)
    }
  } else if (is.null(p)) {
    stop("Give `genotypes` (0/1/2) or the allele frequencies `p`.",
         call. = FALSE)
  }

  # loci with no estimate (monomorphic, or dominance not estimable) go
  keep <- stats::complete.cases(a)
  if (!is.null(d)) keep <- keep & stats::complete.cases(d)
  if (!all(keep)) {
    message(sum(!keep), ngettext(sum(!keep), " locus", " loci"),
            " without estimated effects dropped: ",
            paste(utils::head(rownames(a)[!keep], 5L), collapse = ", "))
  }

  if (is.null(p)) {
    p <- colMeans(genotypes, na.rm = TRUE) / 2
  } else if (length(p) != M) {
    stop("`p` must have one allele frequency per locus of `a` (", M, ").",
         call. = FALSE)
  }
  if (any(!is.finite(p)) || any(p < 0 | p > 1))
    stop("Allele frequencies `p` must lie in [0, 1].", call. = FALSE)

  if (ld && is.null(lambda) && !is.null(genotypes)) {
    lambda <- stats::cov(genotypes, use = "pairwise.complete.obs") / 2
  }
  if (!ld || is.null(lambda)) lambda <- matrix(0, M, M)
  lambda <- as.matrix(lambda)
  if (!identical(dim(lambda), c(M, M)))
    stop("`lambda` must be ", M, " x ", M, ", one row and column per locus.",
         call. = FALSE)

  a <- a[keep, , drop = FALSE]
  if (!is.null(d)) d <- d[keep, , drop = FALSE]
  p <- p[keep]
  lambda <- lambda[keep, keep, drop = FALSE]
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
  list(G = G, D = D, full_sib = G / 2 + D / 4, alpha = alpha,
       p = stats::setNames(p, rownames(a)), loci = rownames(a))
}
