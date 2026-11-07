#' Place a subset covariance matrix into the full trait space
#'
#' @param small Square subset covariance matrix.
#' @param index Traits the subset covers, as indices into the full set.
#' @param K Number of traits in the full set.
#' @param fill Value for the entries the subset does not cover.
#'
#' @return A `K` by `K` matrix.
#'
#' @examples
#' block <- matrix(c(1, 0.4, 0.4, 2), nrow = 2)
#' indri_inflate(block, index = c(1, 3), K = 4)
#' indri_inflate(block, index = c(1, 3), K = 4, fill = 0)
#'
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
#' estimate by the inverse of its reported uncertainty. REML engines report a
#' per-entry standard error; the stanimal engine reports posterior draws, whose
#' per-entry variance supplies the weight.
#'
#' This is the averaging step, not the estimate: the average of overlapping
#' blocks is generally not positive semi-definite. `indri_project_psd()` is what
#' turns it into a covariance matrix.
#'
#' @param blocks List of subset fits, each with a `G` matrix and optionally a
#'   `se` matrix of the same shape.
#' @param subsets Integer matrix of trait indices, one subset per row.
#' @param K Number of traits; inferred from `subsets` when `NULL`.
#' @param weights `"none"` or `"inverse_variance"`. The latter needs each block
#'   to carry a `se` matrix or a `draws` array (used for its per-entry variance
#'   when `se` is absent); blocks that carry neither are treated as weight 1.
#'
#' @return A list with `G` (the averaged matrix) and `counts` (how many subsets
#'   contributed to each entry).
#'
#' @examples
#' # Two estimates of traits 1 and 2, plus blocks covering trait 3.
#' subsets <- rbind(c(1, 2), c(1, 2), c(1, 3), c(2, 3))
#' blocks <- list(
#'   list(G = matrix(c(1, 0.8, 0.8, 1), 2), se = matrix(0.1, 2, 2)),
#'   list(G = matrix(c(1, 0.2, 0.2, 1), 2), se = matrix(0.4, 2, 2)),
#'   list(G = matrix(c(1, 0.3, 0.3, 1), 2), se = matrix(0.2, 2, 2)),
#'   list(G = matrix(c(1, 0.4, 0.4, 1), 2), se = matrix(0.2, 2, 2))
#' )
#' equal <- indri_assemble(blocks, subsets, K = 3)
#' weighted <- indri_assemble(blocks, subsets, K = 3,
#'                            weights = "inverse_variance")
#' equal$G
#' equal$counts
#' weighted$G[1, 2]  # Closer to the more precise first block.
#'
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
# engine reported nothing to weight by. REML engines report `se`; the stanimal
# engine reports posterior draws instead, so its precision is the inverse
# posterior variance.
.indri_weights <- function(block, weights) {
  if (weights == "none") return(NULL)

  precision <- function(sd) {
    sd <- as.matrix(sd)
    sd[!is.finite(sd) | sd <= 0] <- NA_real_
    out <- 1 / sd^2
    # A missing standard error should not silently drop the estimate: fall back
    # to weight 1 for that entry.
    out[is.na(out)] <- 1
    out
  }

  if (!is.null(block$se)) return(precision(block$se))

  draws <- block$draws
  if (!is.null(draws)) return(precision(apply(draws, c(1L, 2L), stats::sd)))

  NULL
}

#' Assemble block posterior draws into one array
#'
#' Combines the per-block posterior draws into a `K` by `K` by `D` array. Either
#' way the result carries block-level uncertainty through; it is **not** a
#' posterior over the joint matrix, because subsets share traits and data and
#' are therefore dependent. See `NOTES.md`.
#'
#' `method = "average"` averages draw `d` across every block that covers an
#' entry, exactly as the means are averaged. Because the block chains are
#' independent, this shrinks the spread of an entry covered `c` times by about
#' `1 / sqrt(c)`, while the block estimates themselves, fitted to the same
#' individuals, have nearly perfectly correlated errors: the band is too narrow
#' and gets narrower the better a pair is covered.
#'
#' `method = "mixture"` instead takes each entry of draw `d` from one covering
#' block chosen at random, so an entry's draws are a mixture of its block
#' posteriors: the spread is the average within-block spread plus the spread
#' between blocks. Entries are drawn independently of one another, so a mixed
#' covariance draw need not be positive semi-definite; mix per-block
#' correlation draws (pass them as `draws`) when correlations are the target.
#'
#' @param draws List of `s` by `s` by `D` arrays, one per subset.
#' @param subsets Integer matrix of trait indices, one subset per row.
#' @param K Number of traits; inferred from `subsets` when `NULL`, which is only
#'   right when the subsets happen to reach the last trait.
#' @param method `"average"` or `"mixture"`; see Details.
#' @param seed Optional seed for the block choice in `"mixture"`. The caller's
#'   random state is left untouched.
#'
#' @return A `K` by `K` by `D` array. Entries no subset covers are 0 for
#'   `"average"` and `NA` for `"mixture"`.
#'
#' @examples
#' # Illustrative draws for the same two-trait subset from two block fits.
#' subsets <- rbind(c(1, 2), c(1, 2))
#' make_draws <- function(correlations) {
#'   array(unlist(lapply(correlations, function(r)
#'     matrix(c(1, r, r, 1), 2))), dim = c(2, 2, length(correlations)))
#' }
#' draws <- list(make_draws(c(0.1, 0.2, 0.3, 0.4)),
#'               make_draws(c(0.5, 0.6, 0.7, 0.8)))
#' averaged <- indri_assemble_draws(draws, subsets)
#' mixed <- indri_assemble_draws(draws, subsets, method = "mixture")
#' averaged[1, 2, ]
#' mixed[1, 2, ]
#' # Mixture entries are pointwise draws, not a joint PSD matrix posterior.
#'
#' @export
indri_assemble_draws <- function(draws, subsets, K = NULL,
                                 method = c("average", "mixture"),
                                 seed = NULL) {
  method <- match.arg(method)
  if (is.null(K)) K <- max(subsets)
  if (length(draws) != nrow(subsets))
    stop("`draws` and `subsets` disagree on the number of subsets.",
         call. = FALSE)
  D <- dim(draws[[1L]])[3L]
  if (method == "mixture")
    return(.indri_mix_draws(draws, subsets, K, D, seed))
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

# For every entry (a, b), a <= b, pick one covering block per draw and copy its
# value; mirror to (b, a).
.indri_mix_draws <- function(draws, subsets, K, D, seed) {
  for (block in draws)
    if (dim(block)[3L] != D)
      stop("every draw array must have the same number of draws.", call. = FALSE)
  restore <- .indri_use_seed(seed)
  on.exit(restore(), add = TRUE)

  out <- array(NA_real_, c(K, K, D))
  for (b in seq_len(K)) {
    for (a in seq_len(b)) {
      covering <- which(apply(subsets, 1L, function(r) a %in% r && b %in% r))
      if (!length(covering)) next
      values <- vapply(covering, function(j) {
        index <- subsets[j, ]
        draws[[j]][match(a, index), match(b, index), ]
      }, numeric(D))
      values <- matrix(values, nrow = D)          # D x n_covering
      pick <- sample.int(length(covering), D, replace = TRUE)
      out[a, b, ] <- out[b, a, ] <- values[cbind(seq_len(D), pick)]
    }
  }
  out
}
