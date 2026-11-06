#' Fit the assembled blocks to the closest positive semi-definite matrix
#'
#' The average of overlapping subset estimates is not a covariance matrix: it is
#' a patchwork, and nothing constrains it to be positive semi-definite. This fits
#' the full matrix that is closest to all the subset estimates at once,
#'
#'   `C* = argmin over C >= 0 of sum_j || C[I_j, I_j] - G_j ||_F^2`
#'
#' by writing `C = L L'` with `L` lower triangular, which makes every candidate
#' positive semi-definite by construction, and optimising the `K (K + 1) / 2`
#' free entries of `L` with L-BFGS-B. The average is used as the starting point.
#'
#' @param blocks List of subset fits, each with a `G` matrix.
#' @param subsets Integer matrix of trait indices, one subset per row.
#' @param K Number of traits; inferred from `subsets` when `NULL`.
#' @param init Starting matrix. Defaults to the average of the blocks.
#' @param control Optional list of controls for [stats::optim()].
#'
#' @return A list with `matrix`, the fitted PSD matrix, plus `convergence`,
#'   `value` and the raw `optim` result.
#' @export
indri_project_psd <- function(blocks, subsets, K = NULL, init = NULL,
                              control = list()) {
  if (is.null(K)) K <- max(subsets)
  if (is.null(init)) {
    init <- indri_assemble(blocks, subsets, K)$G
  }
  init <- as.matrix(init)
  if (anyNA(init)) {
    # entries no subset covered: leave them at zero, they are unconstrained by
    # the objective anyway
    init[is.na(init)] <- 0
  }
  ridge <- max(1e-8, 1e-8 * max(abs(diag(init))))
  start <- tryCatch(t(chol(init + ridge * diag(K))),
                    error = function(e) diag(sqrt(pmax(diag(init), 1e-8)), K))
  theta0 <- start[lower.tri(start, diag = TRUE)]

  objective <- function(theta) {
    candidate <- .indri_from_theta(theta, K)
    total <- 0
    for (i in seq_along(blocks)) {
      index <- subsets[i, ]
      residual <- candidate[index, index, drop = FALSE] - as.matrix(blocks[[i]]$G)
      total <- total + sum(residual^2)
    }
    total
  }

  fit <- stats::optim(theta0, objective, method = "L-BFGS-B",
                      control = utils::modifyList(list(maxit = 1000), control))
  list(matrix = .indri_from_theta(fit$par, K),
       convergence = fit$convergence,
       value = fit$value,
       optim = fit)
}

# Rebuild C = L L' from the free entries of the lower triangle of L.
.indri_from_theta <- function(theta, K) {
  L <- matrix(0, K, K)
  L[lower.tri(L, diag = TRUE)] <- theta
  L %*% t(L)
}

# Smallest eigenvalue of a symmetric matrix, for reporting whether a result is
# usable. Negative values below tolerance are rounding.
.indri_min_eigen <- function(m) {
  min(eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values)
}
