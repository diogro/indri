#' Place a subset covariance matrix into the full trait space
#'
#' @param small Square subset covariance matrix.
#' @param index Traits the subset covers, as indices into the full set.
#' @param K Number of traits in the full set.
#' @param fill Value for the entries the subset does not cover.
#'
#' @return A `K` by `K` matrix.
#' @export
indri_inflate <- function(small, index, K, fill = NA_real_) {
  small <- as.matrix(small)
  if (length(index) != nrow(small))
    stop("`index` must have one entry per row of `small`.", call. = FALSE)
  out <- matrix(fill, K, K)
  out[index, index] <- small
  out
}

#' Average the subset estimates into one matrix
#'
#' Every entry of the assembled matrix is a weighted mean over the subsets that
#' contain both of its traits. With `weights = "none"` all subsets count
#' equally, which is what LiMMBo does; `"inverse_variance"` weights each
#' estimate by the inverse of its reported uncertainty, which the engines all
#' provide in one form or another.
#'
#' This is the averaging step, not the estimate: the average of overlapping
#' blocks is generally not positive semi-definite. `indri_project_psd()` is what
#' turns it into a covariance matrix.
#'
#' @param blocks List of subset fits, each with a `G` matrix and optionally a
#'   `se` matrix of the same shape.
#' @param subsets Integer matrix of trait indices, one subset per row.
#' @param K Number of traits; inferred from `subsets` when `NULL`.
#' @param weights `"none"` or `"inverse_variance"`.
#'
#' @return A list with `G` (the averaged matrix) and `counts` (how many subsets
#'   contributed to each entry).
#' @export
indri_assemble <- function(blocks, subsets, K = NULL,
                           weights = c("none", "inverse_variance")) {
  weights <- match.arg(weights)
  if (length(blocks) != nrow(subsets))
    stop("`blocks` and `subsets` disagree on the number of subsets.",
         call. = FALSE)
  if (is.null(K)) K <- max(subsets)

  accumulated <- matrix(0, K, K)
  weight_sum <- matrix(0, K, K)
  counts <- matrix(0L, K, K)

  for (i in seq_along(blocks)) {
    index <- subsets[i, ]
    G <- blocks[[i]]$G
    w <- .indri_weights(blocks[[i]], weights)
    if (is.null(w)) w <- matrix(1, nrow(G), ncol(G))
    value <- indri_inflate(G * w, index, K, fill = 0)
    weight <- indri_inflate(w, index, K, fill = 0)
    accumulated <- accumulated + value
    weight_sum <- weight_sum + weight
    counts[index, index] <- counts[index, index] + 1L
  }

  assembled <- matrix(NA_real_, K, K)
  covered <- weight_sum > 0
  assembled[covered] <- accumulated[covered] / weight_sum[covered]
  list(G = assembled, counts = counts)
}

# Inverse-variance weights for one subset, or NULL when unweighted or when the
# engine reported nothing to weight by.
.indri_weights <- function(block, weights) {
  if (weights == "none") return(NULL)
  se <- block$se
  if (is.null(se)) return(NULL)
  se <- as.matrix(se)
  se[!is.finite(se) | se <= 0] <- NA_real_
  out <- 1 / se^2
  # A missing standard error should not silently drop the estimate: fall back to
  # weight 1 for that entry.
  out[is.na(out)] <- 1
  out
}

#' Average bootstrap draws across subsets
#'
#' Assembles the per-draw block posteriors into one array, drawing `d` from every
#' block and averaging the overlaps exactly as the means are averaged. The result
#' carries block-level uncertainty through; it is **not** a posterior over the
#' joint matrix, because subsets share traits and are therefore dependent. See
#' `NOTES.md`.
#'
#' @param draws List of `s` by `s` by `D` arrays, one per subset.
#' @param subsets Integer matrix of trait indices, one subset per row.
#'
#' @return A `K` by `K` by `D` array.
#' @export
indri_assemble_draws <- function(draws, subsets) {
  K <- max(subsets)
  D <- dim(draws[[1L]])[3L]
  accumulated <- array(0, c(K, K, D))
  counts <- matrix(0, K, K)

  for (i in seq_along(draws)) {
    index <- subsets[i, ]
    block <- draws[[i]]
    if (dim(block)[3L] != D)
      stop("every draw array must have the same number of draws.", call. = FALSE)
    big <- array(0, c(K, K, D))
    big[index, index, ] <- block
    accumulated <- accumulated + big
    counts[index, index] <- counts[index, index] + 1
  }

  counts[counts == 0] <- 1
  sweep(accumulated, c(1L, 2L), counts, "/")
}
