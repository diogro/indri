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
#'   The default depends on the engine: `min(floor(sqrt(K)), 4)`, clamped to at
#'   least 2, for `"stanimal"`, and 2 (blocks of two traits) for the REML
#'   engines, whose multi-trait fits converge poorly (sommer's badly) and were
#'   measured to gain no accuracy from larger blocks.
#' @param min_cooccurrence Minimum number of times every trait pair must be
#'   sampled. The default depends on the engine: 3 for `"stanimal"`, 1 for the
#'   REML engines. A REML fit of the same pair on the same data is the same
#'   number, so with two-trait blocks repeats add cost and nothing else.
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
#' @param cores Total number of concurrent processes for the whole fit,
#'   not the number of subsets to run at once. Defaults to 1, except for
#'   `"asreml"`, which defaults to 2: every worker checks out its own ASReml
#'   licence seat, and a standard licence allows two at once, on one machine at
#'   a time. When asreml is already loaded in the session (which usually holds
#'   a seat) the default drops to 1 and the blocks are fitted in-session. Use
#'   `cores = 1` with a single-seat licence; more workers than free seats stop
#'   with a licence error. The budget is divided between
#'   blocks and chains: a block fit supplies `chains` processes by default (four
#'   for `"stanimal"`, its own default), so `cores = 4` with the default engine
#'   arguments runs one block at a time with four chains in parallel, while
#'   `cores = 4` with `engine_args = list(chains = 1)` runs four blocks at once.
#'   `chains` is never changed here; only how many run concurrently.
#' @param spread How `G_draws` and `corrG_draws` combine the block posteriors.
#'   `"mixture"` (the default) takes each entry of each draw from one covering
#'   block at random; `"average"` averages draw `d` across the covering blocks.
#'   Blocks are fitted to the same individuals, so averaging their independent
#'   chains is overconfident by roughly the square root of a pair's coverage:
#'   on a 12-trait simulation the averaged 95% band covered the true
#'   correlation 40% of the time and the mixture 89% (the full model: 80%).
#'   See [indri_assemble_draws()].
#' @param seed Optional seed for the subset construction. The engine's sampler
#'   seed is passed through `engine_args`.
#'
#' @return An object of class `indri_fit`. The estimate is `G` (and `corrG`).
#'   For engines that return posterior draws, `G_draws` and `corrG_draws`
#'   (traits x traits x draws) carry the assembled block uncertainty; use
#'   `corrG_draws` for intervals on correlations, since it mixes each block's
#'   own correlation draws and stays within \[-1, 1\]. Mixed draws are
#'   pointwise: their entries are drawn independently, so a single draw need
#'   not be positive semi-definite, and they support per-entry intervals, not
#'   functions of the whole matrix. Neither is a posterior over the joint
#'   matrix.
#'
#' @export
fit_G <- function(Y, X = NULL, A, s = NULL, min_cooccurrence = NULL,
                  subsets = NULL, engine = "stanimal", engine_args = list(),
                  assembly = c("ls", "average"),
                  weights = c("none", "inverse_variance"),
                  missing = c("error", "complete_block"),
                  spread = c("mixture", "average"),
                  cores = NULL, seed = NULL) {
  call <- match.call()
  assembly <- match.arg(assembly)
  weights <- match.arg(weights)
  missing <- match.arg(missing)
  spread <- match.arg(spread)
  engine <- indri_engine(engine)
  if (is.null(cores)) cores <- .indri_default_cores(engine)

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
    defaults <- .indri_engine_defaults(engine, K)
    if (is.null(s)) s <- defaults$s
    if (is.null(min_cooccurrence)) min_cooccurrence <- defaults$min_cooccurrence
    subsets <- indri_subsets(K, s, min_cooccurrence = min_cooccurrence,
                             seed = seed)
  }
  subsets <- as.matrix(subsets)
  if (max(subsets) > K)
    stop("`subsets` refers to traits beyond the ", K, " in `Y`.", call. = FALSE)
  # indri_subsets() guarantees min_cooccurrence; hand-built subsets may leave
  # pairs out on purpose, but those entries of G then come from the PSD fill-in
  # alone, not from data, and the caller should know.
  uncovered <- sum(indri_pair_coverage(subsets, K)[upper.tri(diag(K))] == 0L)
  if (uncovered > 0L)
    warning(uncovered, " trait pair(s) are not in any subset; their entries ",
            "of `G` are not estimated from the data.", call. = FALSE)
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
    tryCatch(
      parallel::parLapply(cluster, indices, .indri_fit_block, Y = Y_std,
                          design = X, a_forms = a_forms, engine = engine,
                          engine_args = engine_args, subsets = subsets),
      error = function(e) {
        if (grepl("licen", conditionMessage(e), ignore.case = TRUE))
          stop("a block worker could not check out an ", engine$name,
               " licence (", conditionMessage(e), "). Each worker needs ",
               "its own seat: lower `cores` to the number of free seats ",
               "(`cores = 1` fits the blocks in this session).",
               call. = FALSE)
        stop(e)
      })
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
  corrG_draws <- NULL
  if (isTRUE(engine$capabilities$draws)) {
    block_draws <- lapply(blocks, `[[`, "draws")
    if (all(!vapply(block_draws, is.null, logical(1L)))) {
      G_draws <- indri_assemble_draws(block_draws, subsets, K,
                                      method = spread, seed = seed)
      # Correlations are mixed from each block's own correlation draws, as in
      # the coverage study, so every entry is a real correlation.
      block_corr <- lapply(block_draws, .indri_corr_draws)
      corrG_draws <- indri_assemble_draws(
        block_corr, subsets, K, method = spread,
        seed = if (is.null(seed)) NULL else seed + 1L)
      for (k in seq_len(K)) corrG_draws[k, k, ] <- 1
    }
  }

  # The engine may or may not carry trait names through; naming is indri's job.
  trait_names <- colnames(Y)
  dimnames(G) <- list(trait_names, trait_names)
  dimnames(assembled$G) <- list(trait_names, trait_names)
  corrG <- stats::cov2cor(G)
  dimnames(corrG) <- list(trait_names, trait_names)
  if (!is.null(G_draws)) {
    dimnames(G_draws) <- list(trait_names, trait_names, NULL)
    dimnames(corrG_draws) <- list(trait_names, trait_names, NULL)
  }

  structure(list(
    G = G,
    G_bar = assembled$G,
    G_draws = G_draws,
    corrG = corrG,
    corrG_draws = corrG_draws,
    sigma_G = sqrt(pmax(diag(G), 0)),
    counts = assembled$counts,
    coverage = indri_pair_coverage(subsets, K),
    subsets = subsets,
    s = s,
    min_cooccurrence = min_cooccurrence,
    engine = engine$name,
    engine_args = engine_args,
    assembly = assembly,
    weights = weights,
    spread = if (is.null(G_draws)) NA_character_ else spread,
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

# Workers when the caller does not say. asreml defaults to two, one licence
# seat each (a standard licence allows two processes, on one machine at a
# time). If asreml is already loaded here, this session has most likely
# checked out a seat itself, and two more workers would exceed the licence, so
# fit in-session instead.
.indri_default_cores <- function(engine) {
  cores <- engine$capabilities$default_cores %||% 1L
  if (identical(engine$name, "asreml") && "asreml" %in% loadedNamespaces())
    cores <- 1L
  cores
}

# Subset size and co-occurrence when the caller gives neither. stanimal keeps
# the sqrt(K) rule and c = 3. The REML engines use pairs, once each: sweeps at
# K = 16 and 30 (examples/block_size_sweep.R) found no accuracy gain from
# larger blocks for asreml, sommer's multi-trait fits mostly fail to converge,
# and a deterministic fit of the same pair is the same number every time.
.indri_engine_defaults <- function(engine, K) {
  caps <- engine$capabilities
  s <- caps$default_s %||% .indri_default_s(K)
  list(s = max(2L, min(as.integer(s), K)),
       min_cooccurrence = caps$default_cooccurrence %||% 3L)
}

# Correlation draws of one block: cov2cor() of every draw.
.indri_corr_draws <- function(draws) {
  out <- draws
  for (d in seq_len(dim(draws)[3L])) out[, , d] <- stats::cov2cor(draws[, , d])
  out
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
