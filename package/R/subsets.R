#' Trait subsets for the block estimator
#'
#' Draws the trait subsets on which the small models are fitted. This is the
#' set-cover construction used by LiMMBo from version 1.0 on: each subset is
#' built from a window of trait pairs, and the number of subsets is chosen so
#' that every pair is covered at least `min_cooccurrence` times.
#'
#' Three things differ from LiMMBo's implementation. The permutation that
#' orders the pairs is seeded here, so the result is reproducible (LiMMBo draws
#' it from the global numpy generator and is not), and the caller's random seed
#' is restored afterwards. And the coverage is guaranteed: LiMMBo's passes
#' cover each pair only `ceiling(min_cooccurrence / 2)` times by construction
#' and leave the rest to chance, which often falls short; indri adds greedily
#' built subsets of the same width until every pair reaches
#' `min_cooccurrence`.
#'
#' Note that the subset size is not exactly `s`. A subset covers `s - 1` trait
#' pairs, which is up to `2 * (s - 1)` distinct traits; the subsets therefore
#' come out larger than `s`, and the assembled blocks use
#' `ncol(subsets)` rather than `s`.
#'
#' @param K Number of traits.
#' @param s Subset size parameter, at least 2. Controls how many trait pairs
#'   each subset covers.
#' @param min_cooccurrence Minimum number of times every trait pair must be
#'   sampled.
#' @param seed Optional seed for the subset construction. The caller's random
#'   state is left untouched.
#'
#' @return An integer matrix with one distinct subset per row, holding 1-based
#'   indices of the traits in that subset. Every pair is in at least
#'   `min_cooccurrence` distinct subsets, capped at the number of distinct
#'   subsets of this width that contain it (one, for two-trait subsets); the
#'   coverage reached is stored in the `min_cooccurrence` attribute and the
#'   pair counts in the `coverage` attribute.
#'
#' @examples
#' subsets <- indri_subsets(K = 6, s = 3, min_cooccurrence = 3, seed = 42)
#' subsets
#' ncol(subsets)  # Actual block width; it can be larger than s.
#' coverage <- attr(subsets, "coverage")
#' coverage
#' min(coverage[upper.tri(coverage)])
#'
#' @export
indri_subsets <- function(K, s, min_cooccurrence = 3L, seed = NULL) {
  K <- .indri_count(K, "K", min = 2L)
  s <- .indri_count(s, "s", min = 2L)
  min_cooccurrence <- .indri_count(min_cooccurrence, "min_cooccurrence", min = 1L)
  if (s > K)
    stop("`s` (", s, ") cannot exceed the number of traits (", K, ").",
         call. = FALSE)

  restore <- .indri_use_seed(seed)
  on.exit(restore(), add = TRUE)

  n_pairs <- as.integer(K * (K - 1L) / 2L)
  pairs <- utils::combn(K, 2L)
  # Position `perm[c]` refers to the c-th trait pair, as in LiMMBo's
  # `number_to_tuple[order[counter]] <- (i, j)`.
  perm <- sample.int(n_pairs, n_pairs)
  pair_a <- integer(n_pairs)
  pair_b <- integer(n_pairs)
  pair_a[perm] <- pairs[1L, ]
  pair_b[perm] <- pairs[2L, ]

  inflated <- s - 1L
  goal <- n_pairs %/% inflated + 1L
  n_covers <- as.integer(ceiling(min_cooccurrence / 2))

  offsets <- seq_len(inflated) - 1L
  subsets <- vector("list", n_covers * goal)
  k <- 0L
  for (i in seq_len(n_covers)) {
    for (j in seq_len(goal)) {
      positions <- ((offsets + (i - 1L) + (j - 1L) * inflated) %% n_pairs) + 1L
      k <- k + 1L
      subsets[[k]] <- sort(unique(c(pair_a[positions], pair_b[positions])))
    }
  }

  # Subsets come out at slightly different sizes when a pair shares a trait with
  # another in the same window. LiMMBo pads the short ones with randomly chosen
  # traits, which is what makes the assembled matrix rectangular.
  width <- max(lengths(subsets))
  for (i in seq_along(subsets)) {
    need <- width - length(subsets[[i]])
    if (need > 0L) {
      available <- setdiff(seq_len(K), subsets[[i]])
      subsets[[i]] <- sort(c(subsets[[i]], .indri_sample(available, need)))
    }
  }

  out <- do.call(rbind, subsets)
  # Coverage counts distinct trait combinations: fitting the same block twice
  # estimates a pair in the same context twice, which adds nothing (and for a
  # deterministic REML engine returns the same number).
  out <- out[!duplicated(out), , drop = FALSE]
  coverage <- indri_pair_coverage(out, K)

  # A pair can only be in choose(K - 2, width - 2) distinct blocks of this
  # width (one, for pairs), so that caps the coverage any construction can
  # reach.
  target <- min(min_cooccurrence, choose(K - 2L, width - 2L))

  # The passes above cover every pair `n_covers = ceiling(c / 2)` times by
  # construction; LiMMBo relies on incidental co-membership for the rest, and
  # that often falls short (min coverage 2 for c = 3 is common). Top up with
  # greedily built, new subsets of the same width until every pair reaches
  # the target.
  out <- .indri_top_up(out, coverage, K, width, target)
  attr(out, "min_cooccurrence") <- target
  out
}

# Add subsets of `width` traits until every pair is covered at least
# `min_cooccurrence` times. Each new subset starts from one of the most
# under-covered pairs and grows by the trait that closes the most remaining
# deficits, so every addition covers at least one deficient pair and the loop
# terminates.
.indri_top_up <- function(subsets, coverage, K, width, min_cooccurrence) {
  extra <- list()
  keys <- apply(subsets, 1L, paste, collapse = "-")
  repeat {
    deficit <- min_cooccurrence - coverage
    deficit[!upper.tri(deficit)] <- 0L
    if (all(deficit <= 0L)) break
    deficit <- pmax(deficit, 0L)
    deficit <- deficit + t(deficit)

    worst <- which(deficit == max(deficit) & upper.tri(deficit), arr.ind = TRUE)
    members <- unname(worst[.indri_sample(seq_len(nrow(worst)), 1L), ])
    worst_pair <- members
    while (length(members) < width) {
      others <- setdiff(seq_len(K), members)
      gain <- colSums(deficit[members, others, drop = FALSE] > 0L)
      members <- c(members, .indri_sample(others[gain == max(gain)], 1L))
    }
    members <- sort(members)
    key <- paste(members, collapse = "-")
    if (key %in% keys) {
      # the greedy choice already exists: take the best block not yet fitted
      # that contains the most under-covered pair (one always exists, since
      # the target never exceeds the number of distinct blocks with the pair)
      members <- .indri_new_block(worst_pair, K, width, deficit, keys)
      key <- paste(members, collapse = "-")
    }
    keys <- c(keys, key)
    extra[[length(extra) + 1L]] <- members
    coverage[members, members] <- coverage[members, members] + 1L
  }
  if (length(extra)) subsets <- rbind(subsets, do.call(rbind, extra))
  storage.mode(subsets) <- "integer"
  attr(subsets, "coverage") <- coverage
  subsets
}

# The block of `width` traits containing `pair` that is not in `keys` and
# covers the most pairs still short of their target. Candidates are every
# completion when there are few, a random sample otherwise.
.indri_new_block <- function(pair, K, width, deficit, keys) {
  rest <- setdiff(seq_len(K), pair)
  need <- width - 2L
  if (need == 0L) return(sort(pair))
  n <- choose(length(rest), need)
  candidates <- if (n <= 5000) utils::combn(rest, need, simplify = FALSE) else
    lapply(seq_len(2000L), function(i) .indri_sample(rest, need))
  best <- NULL
  best_gain <- -1
  for (cand in candidates) {
    block <- sort(c(pair, cand))
    if (paste(block, collapse = "-") %in% keys) next
    gain <- sum(deficit[block, block][upper.tri(diag(width))] > 0L)
    if (gain > best_gain) {
      best <- block
      best_gain <- gain
    }
  }
  if (is.null(best))
    stop("no new block contains pair ", paste(pair, collapse = "-"),
         "; this is a bug in the top-up.", call. = FALSE)
  best
}

# sample() treats a single number n as 1:n; this never does.
.indri_sample <- function(x, size) {
  x[sample.int(length(x), size)]
}

# Number of times each trait pair is covered by a subset matrix.
indri_pair_coverage <- function(subsets, K = NULL) {
  if (is.null(K)) K <- max(subsets)
  coverage <- matrix(0L, K, K)
  for (i in seq_len(nrow(subsets))) {
    idx <- subsets[i, ]
    coverage[idx, idx] <- coverage[idx, idx] + 1L
  }
  coverage
}

.indri_count <- function(x, arg, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < min || x != trunc(x))
    stop("`", arg, "` must be a whole number of at least ", min, ".",
         call. = FALSE)
  as.integer(x)
}

# Seed the generator for a computation and hand back a function that restores
# whatever state the caller had, so that calling indri does not move the user's
# random stream along.
.indri_use_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  had_state <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_state <- if (had_state) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  function() {
    if (had_state) {
      assign(".Random.seed", old_state, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}
