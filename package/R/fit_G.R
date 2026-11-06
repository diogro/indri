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
#'   Defaults to `floor(sqrt(K))`, clamped to at least 2.
#' @param min_cooccurrence Minimum number of times every trait pair must be
#'   sampled.
#' @param subsets Optional precomputed subset matrix, from [indri_subsets()].
#' @param engine Block engine: `"stanimal"` today, with `"sommer"` and
#'   `"asreml"` designed and not yet implemented.
#' @param engine_args List of arguments passed to the engine's fit function.
#' @param assembly `"ls"` fits the closest positive semi-definite matrix to all
#'   the subset estimates; `"average"` keeps the plain average, nudged onto the
#'   positive semi-definite cone.
#' @param weights `"none"` for LiMMBo's equal weights, or `"inverse_variance"`
#'   to weight each subset estimate by its reported uncertainty.
#' @param missing What to do when `Y` has gaps and the engine cannot handle them:
#'   `"error"`, or `"complete_block"` to drop incomplete individuals once, before
#'   subsetting, so every block sees the same rows.
#' @param cores Number of worker processes used for the subset fits. Fork-based,
#'   so ignored on Windows.
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
    if (is.null(s)) s <- max(2L, min(K - 1L, floor(sqrt(K))))
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

  fit_one <- function(i) {
    engine$fit(Y_std[, subsets[i, ], drop = FALSE], X, a_forms, engine_args)
  }
  indices <- seq_len(nrow(subsets))
  cores <- max(1L, min(as.integer(cores), length(indices)))
  blocks <- if (cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(indices, fit_one, mc.cores = cores)
  } else {
    lapply(indices, fit_one)
  }
  failed <- vapply(blocks, function(b) inherits(b, "try-error") || !is.list(b),
                   logical(1))
  if (any(failed))
    stop("subset fits failed for ", sum(failed), " of ", nrow(subsets),
         " subsets; the first error is: ",
         as.character(blocks[[which(failed)[1L]]]), call. = FALSE)

  assembled <- indri_assemble(blocks, subsets, K, weights = weights)
  G <- switch(assembly,
    ls = indri_project_psd(blocks, subsets, K, init = assembled$G)$matrix,
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
      G_draws <- indri_assemble_draws(block_draws, subsets)
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
    coverage = indri_pair_coverage(subsets),
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
    blocks = blocks,
    diagnostics = .indri_diagnostics(blocks),
    call = call
  ), class = "indri_fit")
}

.indri_diagnostics <- function(blocks) {
  rhat <- vapply(blocks, function(b) {
    value <- b$converged
    if (is.null(value) || length(value) != 1L) NA_real_ else as.numeric(value)
  }, numeric(1))
  list(
    rhat = rhat,
    max_rhat = suppressWarnings(max(rhat, na.rm = TRUE)),
    converged = all(is.na(rhat)) || max(rhat, na.rm = TRUE) < 1.05
  )
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
