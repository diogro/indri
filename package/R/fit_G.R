#' Estimate the additive genetic covariance matrix
#'
#' Fits a multivariate animal model on many small trait subsets and reassembles
#' the subset estimates into one positive semi-definite matrix. The estimator is
#' LiMMBo's: the speed comes from giving each subset its own nuisance residual
#' covariance, so the many small fits replace one large one.
#'
#' Traits are standardised once, using every observed value, before the subsets
#' are drawn. Every engine therefore sees the same numbers and estimates the same
#' matrix on the same scale, whatever it does internally. The returned `G` is on
#' that standardised scale; multiply it back with
#' `diag(y_sd) %*% G %*% diag(y_sd)` for the original units.
#'
#' @param Y Phenotype matrix, individuals in rows and traits in columns. Gaps are
#'   allowed when the engine can handle them; see `missing`.
#' @param X Fixed-effect design. Defaults to an intercept.
#' @param A Relationship matrix, individuals in the same order as `Y`. When both
#'   have row names they are matched by name.
#' @param s Subset size parameter. Controls how many trait pairs each subset
#'   covers, not the number of traits in a subset; see [indri_subsets()].
#'   Defaults to `min(floor(sqrt(K)), 4)`, clamped to at least 2.
#' @param min_cooccurrence Minimum number of times every trait pair must be
#'   sampled.
#' @param subsets Optional precomputed subset matrix, from [indri_subsets()].
#' @param engine Block engine: `"stanimal"` (posterior draws), `"sommer"` or
#'   `"asreml"` (REML point estimates).
#' @param engine_args List of arguments passed to the engine's fit function.
#' @param assembly `"ls"` fits the closest positive semi-definite matrix to all
#'   the subset estimates; `"average"` keeps the plain average, nudged onto the
#'   positive semi-definite cone.
#' @param weights `"none"` for LiMMBo's equal weights, or `"inverse_variance"`
#'   to weight each subset estimate by its reported uncertainty (per-entry
#'   standard error for REML engines, per-entry posterior variance for the
#'   draws-based stanimal engine). Weights apply to both the average and the
#'   projection.
#' @param missing What to do when `Y` has gaps and the engine cannot handle them:
#'   `"error"`, or `"complete_block"` to drop incomplete individuals once, before
#'   subsetting, so every block sees the same rows.
#' @param cores Total number of concurrent Stan processes for the whole fit,
#'   not the number of subsets to run at once. The budget is divided between
#'   blocks and chains: a block fit supplies `chains` processes by default (four
#'   for `"stanimal"`, its own default), so `cores = 4` with the default engine
#'   arguments runs one block at a time with four chains in parallel, while
#'   `cores = 4` with `engine_args = list(chains = 1)` runs four blocks at once.
#'   `chains` is never changed here; only how many run concurrently.
#' @param seed Optional seed for the subset construction. The engine's sampler
#'   seed is passed through `engine_args`.
#'
#' @return An object of class `indri_fit`. The estimate is `G`; `G_draws` carries
#'   the assembled block uncertainty, which is a spread across dependent blocks
#'   and **not** a posterior over the joint matrix.
#'
#' @export
fit_G <- function(Y, X = NULL, A, s = NULL, min_cooccurrence = 3L,
                  subsets = NULL, engine = "stanimal", engine_args = list(),
                  assembly = c("ls", "average"),
                  weights = c("none", "inverse_variance"),
                  missing = c("error", "complete_block"),
                  cores = 1L, seed = NULL) {
  call <- match.call()
  assembly <- match.arg(assembly)
  weights <- match.arg(weights)
  missing <- match.arg(missing)
  engine <- indri_engine(engine)

  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (ncol(Y) < 2L)
    stop("`Y` needs at least two traits: a subset estimate needs a covariance.",
         call. = FALSE)
  if (is.null(colnames(Y)))
    colnames(Y) <- paste0("trait", seq_len(ncol(Y)))
  if (is.null(rownames(Y)))
    rownames(Y) <- paste0("id", seq_len(nrow(Y)))

  A <- .indri_match_A(A, rownames(Y))
  X <- .indri_design(X, nrow(Y), rownames(Y))

  if (anyNA(Y)) {
    if (isTRUE(engine$capabilities$missing)) {
      # the engine models the gaps
    } else if (missing == "error") {
      stop("`Y` has gaps and the ", engine$name, " engine cannot handle them. ",
           "Use `missing = \"complete_block\"` to drop incomplete individuals.",
           call. = FALSE)
    } else {
      complete <- stats::complete.cases(Y)
      if (sum(complete) < 2L)
        stop("fewer than two individuals are complete across all traits.",
             call. = FALSE)
      Y <- Y[complete, , drop = FALSE]
      X <- X[complete, , drop = FALSE]
      A <- A[complete, complete, drop = FALSE]
    }
  }

  K <- ncol(Y)
  y_mean <- colMeans(Y, na.rm = TRUE)
  y_sd <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  if (any(!is.finite(y_sd)) || any(y_sd <= 0))
    stop("every trait needs at least two observed values and some variation.",
         call. = FALSE)
  Y_std <- sweep(sweep(Y, 2L, y_mean, "-"), 2L, y_sd, "/")

  if (is.null(subsets)) {
    if (is.null(s)) s <- .indri_default_s(K)
    subsets <- indri_subsets(K, s, min_cooccurrence = min_cooccurrence,
                             seed = seed)
  }
  subsets <- as.matrix(subsets)
  if (max(subsets) > K)
    stop("`subsets` refers to traits beyond the ", K, " in `Y`.", call. = FALSE)
  if (ncol(subsets) >= K && assembly == "ls")
    stop("every subset covers every trait, so there is nothing to reassemble; ",
         "use fewer traits per subset or `assembly = \"average\"`.",
         call. = FALSE)

  a_forms <- indri_a_forms(A, need = engine$capabilities$a_form)

  indices <- seq_len(nrow(subsets))
  plan <- .indri_core_plan(cores, engine_args, engine, length(indices))
  engine_args <- plan$engine_args
  blocks <- if (plan$blocks_in_parallel > 1L) {
    cluster <- parallel::makeCluster(plan$blocks_in_parallel)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::parLapply(cluster, indices, .indri_fit_block, Y = Y_std,
                        design = X, a_forms = a_forms, engine = engine,
                        engine_args = engine_args, subsets = subsets)
  } else {
    lapply(indices, .indri_fit_block, Y = Y_std, design = X,
           a_forms = a_forms, engine = engine, engine_args = engine_args,
           subsets = subsets)
  }
  failed <- vapply(blocks, function(b) inherits(b, "try-error") || !is.list(b),
                   logical(1))
  if (any(failed))
    stop("subset fits failed for ", sum(failed), " of ", nrow(subsets),
         " subsets; the first error is: ",
         as.character(blocks[[which(failed)[1L]]]), call. = FALSE)

  if (weights == "inverse_variance" &&
      !any(vapply(blocks, function(b) !is.null(b$se) || !is.null(b$draws),
                  logical(1)))) {
    warning("`weights = \"inverse_variance\"` but the ", engine$name,
            " engine reports neither standard errors nor posterior draws; ",
            "falling back to equal weights.", call. = FALSE)
  }

  assembled <- indri_assemble(blocks, subsets, K, weights = weights)
  G <- switch(assembly,
    ls = indri_project_psd(blocks, subsets, K, init = assembled$G,
                           weights = weights)$matrix,
    average = assembled$G
  )
  if (assembly == "average") {
    # the plain average is not a covariance matrix; nudge it onto the cone so
    # that everything downstream can rely on a usable matrix
    minimum <- .indri_min_eigen(G)
    if (minimum < 0) G <- G + (-minimum + 1e-8) * diag(K)
  }

  G_draws <- NULL
  if (isTRUE(engine$capabilities$draws)) {
    block_draws <- lapply(blocks, `[[`, "draws")
    if (all(!vapply(block_draws, is.null, logical(1L))))
      G_draws <- indri_assemble_draws(block_draws, subsets, K)
  }

  # The engine may or may not carry trait names through; naming is indri's job.
  trait_names <- colnames(Y)
  dimnames(G) <- list(trait_names, trait_names)
  dimnames(assembled$G) <- list(trait_names, trait_names)
  corrG <- stats::cov2cor(G)
  dimnames(corrG) <- list(trait_names, trait_names)
  if (!is.null(G_draws))
    dimnames(G_draws) <- list(trait_names, trait_names, NULL)

  structure(list(
    G = G,
    G_bar = assembled$G,
    G_draws = G_draws,
    corrG = corrG,
    sigma_G = sqrt(pmax(diag(G), 0)),
    counts = assembled$counts,
    coverage = indri_pair_coverage(subsets, K),
    subsets = subsets,
    engine = engine$name,
    engine_args = engine_args,
    assembly = assembly,
    weights = weights,
    y_mean = y_mean,
    y_sd = y_sd,
    y_names = colnames(Y),
    x_names = colnames(X),
    n_individuals = nrow(Y),
    parallel = plan[c("cores", "chains", "blocks_in_parallel",
                      "per_block_cores")],
    blocks = blocks,
    diagnostics = .indri_diagnostics(blocks),
    call = call
  ), class = "indri_fit")
}

.indri_diagnostics <- function(blocks) {
  values <- lapply(blocks, `[[`, "converged")

  # Rhat is only meaningful for engines that report it (stanimal returns the
  # block's max Rhat as a number); REML engines report a logical convergence
  # flag instead, which carries no Rhat.
  rhat <- vapply(values, function(v) {
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) as.numeric(v)
    else NA_real_
  }, numeric(1))

  # A block is converged when it says TRUE, or when its Rhat is below the usual
  # 1.05 threshold. NULL, length != 1 and NA mean "could not tell".
  block_ok <- vapply(values, function(v) {
    if (is.null(v) || length(v) != 1L) return(NA)
    if (is.logical(v)) return(isTRUE(v))
    if (is.numeric(v)) {
      if (!is.finite(v)) return(NA)
      return(v < 1.05)
    }
    NA
  }, logical(1))

  list(
    rhat = rhat,
    max_rhat = suppressWarnings(max(rhat, na.rm = TRUE)),
    converged = if (all(is.na(block_ok))) TRUE else all(block_ok, na.rm = TRUE)
  )
}

# Default subset size: grow with K like sqrt(K), but never past 4. Larger
# blocks stop buying accuracy (and get slower) at K = 50, so 4 is the cap.
.indri_default_s <- function(K) {
  max(2L, min(K - 1L, floor(sqrt(K)), 4L))
}

# Match A to the individuals in Y, by name when both have names.
.indri_match_A <- function(A, ids) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A))
    stop("`A` must be a square matrix.", call. = FALSE)
  if (nrow(A) != length(ids))
    stop("`A` has ", nrow(A), " rows but `Y` has ", length(ids),
         " individuals.", call. = FALSE)
  if (!is.null(rownames(A)) && !is.null(colnames(A))) {
    if (!setequal(rownames(A), ids))
      stop("the row names of `A` do not match the row names of `Y`.",
           call. = FALSE)
    A <- A[ids, ids, drop = FALSE]
  }
  A
}

.indri_design <- function(X, n, ids) {
  if (is.null(X)) {
    X <- matrix(1, n, 1L, dimnames = list(ids, "intercept"))
    return(X)
  }
  X <- as.matrix(X)
  if (nrow(X) != n)
    stop("`X` has ", nrow(X), " rows but `Y` has ", n, " individuals.",
         call. = FALSE)
  if (is.null(colnames(X)))
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (is.null(rownames(X))) rownames(X) <- ids
  storage.mode(X) <- "double"
  X
}

# Split a total process budget between blocks and the chains inside a block.
#
# The engines bring their own default number of chains (four for stanimal, which
# also runs them in parallel), so "four cores" cannot mean "four blocks at once"
# without multiplying processes by four. `chains` is a statistical setting and is
# left alone; only the concurrency changes.
.indri_core_plan <- function(cores, engine_args, engine, n_blocks) {
  cores <- max(1L, as.integer(cores %||% 1L))
  chains <- engine_args$chains %||% engine$capabilities$default_chains %||% 1L
  chains <- max(1L, as.integer(chains))

  per_block_cores <- engine_args$cores %||% min(chains, cores)
  per_block_cores <- max(1L, min(as.integer(per_block_cores), chains))

  blocks_in_parallel <- max(1L, cores %/% per_block_cores)
  blocks_in_parallel <- min(blocks_in_parallel, n_blocks)

  engine_args$cores <- per_block_cores
  list(engine_args = engine_args, cores = cores, chains = chains,
       blocks_in_parallel = blocks_in_parallel,
       per_block_cores = per_block_cores)
}

# Top-level so that a PSOCK worker does not have to be sent the whole of
# fit_G()'s frame.
.indri_fit_block <- function(i, Y, design, a_forms, engine, engine_args,
                             subsets) {
  engine$fit(Y[, subsets[i, ], drop = FALSE], design, a_forms, engine_args)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
