#' Relationship matrices in the forms the engines want
#'
#' The engines do not agree on how the relationship matrix arrives: a Bayesian
#' model takes it dense or as its eigendecomposition, sommer takes its inverse
#' marked as such, and ASReml takes a sparse three-column inverse. Build the
#' forms once per run so that a loop over many trait subsets pays for each one
#' only once.
#'
#' @param A Square relationship matrix.
#' @param need Which forms to build. `"dense"` is always returned.
#'
#' @return A list with some of `dense`, `inverse`, `sparse_inverse` and
#'   `eigen`. `eigen` is `list(Q, lambda)`, the eigenvectors and eigenvalues of
#'   `A` in ascending order, as `stanimal::lmm_animal(eigen_A = )` takes them
#'   (and in the order stanimal would compute them itself, so a block fitted
#'   with it is identical to one fitted without).
#'
#' @examples
#' # A relationship matrix for two pairs of full siblings.
#' A <- kronecker(diag(2), matrix(c(1, 0.5, 0.5, 1), 2))
#' forms <- indri_a_forms(A, need = c("inverse", "eigen"))
#' forms$inverse
#' forms$eigen$lambda
#'
#' # Reconstruct A from the shared eigendecomposition.
#' Q <- forms$eigen$Q
#' Q %*% diag(forms$eigen$lambda) %*% t(Q)
#'
#' @export
indri_a_forms <- function(A, need = c("dense", "inverse", "sparse_inverse",
                                     "eigen")) {
  need <- match.arg(need, several.ok = TRUE)
  A <- as.matrix(A)
  if (nrow(A) != ncol(A))
    stop("`A` must be a square matrix.", call. = FALSE)
  if (anyNA(A))
    stop("`A` contains missing values.", call. = FALSE)

  out <- list(dense = A)
  if ("inverse" %in% need) {
    inverse <- solve(A)
    attr(inverse, "inverse") <- TRUE
    out$inverse <- inverse
  }
  if ("sparse_inverse" %in% need) {
    out$sparse_inverse <- .indri_sparse_inverse(A)
  }
  if ("eigen" %in% need) {
    e <- eigen(A, symmetric = TRUE)
    ascending <- rev(seq_len(nrow(A)))
    out$eigen <- list(Q = e$vectors[, ascending, drop = FALSE],
                      lambda = e$values[ascending])
  }
  out
}

# Three-column sparse inverse in the layout asreml::vm() expects: row, column
# and value for the lower triangle, sorted row-major. Ported from the comparison
# script in stanimal's examples, which builds it this way instead of calling
# asreml::ainverse(). A pedigree-based inverse (nadiv::ainv) avoids the dense
# solve entirely and should be preferred when one is available.
.indri_sparse_inverse <- function(A) {
  inverse <- solve(A)
  small <- abs(inverse) < 1e-10 * max(abs(inverse))
  inverse[small] <- 0
  lower <- lower.tri(inverse, diag = TRUE) & inverse != 0
  index <- which(lower, arr.ind = TRUE)
  index <- index[order(index[, 1L], index[, 2L]), , drop = FALSE]
  out <- cbind(Row = index[, 1L], Column = index[, 2L],
               Ainverse = inverse[index])
  attr(out, "rowNames") <- rownames(A)
  attr(out, "geneticGroups") <- c(0, 0)
  class(out) <- c(class(out), "ginv")
  out
}
