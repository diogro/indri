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

  parts <- .indri_psd_parts(blocks, subsets, K)
  theta0 <- start[parts$lt]

  fit <- stats::optim(
    theta0,
    fn = function(theta) .indri_psd_objective(theta, K, parts),
    gr = function(theta) .indri_psd_gradient(theta, K, parts),
    method = "L-BFGS-B",
    control = utils::modifyList(list(maxit = 1000), control)
  )
  list(matrix = .indri_from_theta(fit$par, K),
       convergence = fit$convergence,
       value = fit$value,
       optim = fit)
}

# Free entries of L: the lower triangle, column by column.
.indri_theta_map <- function(K) {
  which(lower.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
}

# Precompute the per-block pieces so the objective and its gradient do not
# allocate (or as.matrix()) inside the optimiser's hot loop. A numerical
# gradient would evaluate the objective once per free parameter (K (K + 1) / 2
# of them), which for K = 30 and hundreds of blocks is prohibitive.
.indri_psd_parts <- function(blocks, subsets, K) {
  lt <- .indri_theta_map(K)
  list(
    lt = lt,
    G = lapply(blocks, function(b) as.matrix(b$G)),
    idx = lapply(seq_len(nrow(subsets)), function(i) subsets[i, ]),
    rows = lapply(seq_len(K), function(r) which(lt[, 1L] == r))
  )
}

.indri_psd_objective <- function(theta, K, parts) {
  L <- matrix(0, K, K)
  L[parts$lt] <- theta
  C <- L %*% t(L)
  total <- 0
  for (j in seq_along(parts$G)) {
    residual <- C[parts$idx[[j]], parts$idx[[j]], drop = FALSE] - parts$G[[j]]
    total <- total + sum(residual^2)
  }
  total
}

# d/d theta_k of sum_j || C[I_j, I_j] - G_j ||^2, with C = L L'. For a free
# entry L[r, c] (r >= c) the derivative is
#   2 * sum_j ((R_j + R_j') %*% L[I_j, ])[pos(r), c],
# where R_j = C[I_j, I_j] - G_j and pos(r) is r's position in I_j. G_j is
# symmetric in normal use, in which case R_j + R_j' = 2 R_j and the factor is 4.
.indri_psd_gradient <- function(theta, K, parts) {
  L <- matrix(0, K, K)
  L[parts$lt] <- theta
  C <- L %*% t(L)
  grad <- numeric(length(theta))
  for (j in seq_along(parts$G)) {
    idx <- parts$idx[[j]]
    R <- C[idx, idx, drop = FALSE] - parts$G[[j]]
    M <- (R + t(R)) %*% L[idx, , drop = FALSE]   # length(idx) x K
    for (p in seq_along(idx)) {
      for (q in parts$rows[[idx[p]]]) {
        grad[q] <- grad[q] + 2 * M[p, parts$lt[q, 2L]]
      }
    }
  }
  grad
}

# Rebuild C = L L' from the free entries of the lower triangle of L.
.indri_from_theta <- function(theta, K) {
  L <- matrix(0, K, K)
  L[.indri_theta_map(K)] <- theta
  L %*% t(L)
}

# Smallest eigenvalue of a symmetric matrix, for reporting whether a result is
# usable. Negative values below tolerance are rounding.
.indri_min_eigen <- function(m) {
  min(eigen((m + t(m)) / 2, symmetric = TRUE, only.values = TRUE)$values)
}
