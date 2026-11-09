#' Environment-specific genetic covariance (experimental)
#'
#' Estimates the genetic covariance of a set of traits within and between
#' the groups of a discrete environment (sex, diet, site), treating each
#' trait in each environment as its own trait ("character states"). Each
#' individual is observed only in its own environment, so the genetic
#' covariance between environments is identified through relatives in
#' different environments, and the residual covariance between environments
#' is not identified at all: it is zero by construction.
#'
#' @details The response is reshaped to one column per trait and
#' environment (`"trait:level"`), each observed only in that environment's
#' individuals, and [fit_G()] fits it with its usual subsets, assembly and
#' engines. `G` is the assembled genetic covariance of the character states:
#' its diagonal blocks are the within-environment `G` and its off-diagonal
#' blocks the genetic covariance between environments (a cross-environment
#' genetic correlation below 1 is genotype-by-environment interaction in the
#' background). `E` is block diagonal: within each environment, the residual
#' covariance of the traits after `X` minus that environment's block of
#' `G` (scaled for the relatedness within it), floored like [fit_G()]'s `E`.
#'
#' The environment's main effect is carried by the character states' own
#' means, so a column of `X` named by `env` is dropped; every other column
#' of `X` must vary within each environment.
#'
#' This is a trial: the cross-environment covariances rest on relatives in
#' different environments only. In simulations (900 sibs of both sexes, 10
#' and 30 traits; STATUS.md) the cross-environment genetic correlations
#' came out **too low**: 0.67-0.77 where the truth was 1, 0.36-0.45 where it
#' was 0.5, more so with more traits. Scans with this covariance
#' ([indri_gwas()] with this fit) are then conservative for interaction
#' tests (genomic inflation 0.88-0.98), not anti-conservative, but lose
#' power. The cause is under investigation.
#'
#' @param Y Numeric matrix, individuals by traits.
#' @param env The environment of each individual: a vector or factor with
#'   one value per row of `Y` and at least two levels, or the name of a
#'   column of `X`.
#' @param X,A As for [fit_G()].
#' @param engine Block engine; `"sommer"` by default (the stanimal engine
#'   imputes every unobserved character state and is much slower here).
#' @param ... Passed to [fit_G()].
#' @return An object of class `indri_fit_env`: `G` and `E` (character states
#'   by character states, in the units of `Y`), `corr_between` (traits by
#'   pairs of environments: each trait's genetic correlation between two
#'   environments), `levels`, `traits`, `states` (the character-state
#'   names), `env` (each individual's environment, named by individual) and
#'   `fit` (the [fit_G()] fit of the character states). Pass it to
#'   [indri_gwas()] as `fit` to scan with this covariance.
#' @seealso [fit_G()], [indri_gwas()]
#' @export
fit_G_env <- function(Y, env, X = NULL, A, engine = "sommer", ...) {
  Y <- .indri_gwas_matrix(Y, "Y")
  N <- nrow(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("trait", seq_len(ncol(Y)))
  if (is.null(rownames(Y))) rownames(Y) <- paste0("id", seq_len(N))
  traits <- colnames(Y)
  if (is.null(X)) {
    X <- matrix(1, N, 1L, dimnames = list(rownames(Y), "(Intercept)"))
  } else {
    X <- .indri_gwas_matrix(X, "X")
    if (is.null(rownames(X))) rownames(X) <- rownames(Y)
  }

  if (is.character(env) && length(env) == 1L && N > 1L) {
    if (is.null(colnames(X)) || !env %in% colnames(X))
      stop("`env` names no column of `X`.", call. = FALSE)
    column <- env
    env <- X[, column]
    X <- X[, colnames(X) != column, drop = FALSE]
  }
  if (length(env) != N)
    stop("`env` must have one value per row of `Y` (", N, "), or name a ",
         "column of `X`.", call. = FALSE)
  if (anyNA(env))
    stop("`env` has missing values.", call. = FALSE)
  env <- factor(env)
  levels <- levels(env)
  if (length(levels) < 2L)
    stop("`env` must have at least two levels.", call. = FALSE)
  for (l in levels) {
    rows <- env == l
    if (qr(X[rows, , drop = FALSE])$rank < ncol(X))
      stop("`X` is rank deficient within environment ", l, ": every column ",
           "must vary within each environment (the environment's main ",
           "effect is in the character states' means already).",
           call. = FALSE)
  }

  # character states: one column per trait and environment
  states <- as.vector(outer(traits, levels, paste, sep = ":"))
  Ys <- matrix(NA_real_, N, length(states),
               dimnames = list(rownames(Y), states))
  for (l in levels) {
    rows <- env == l
    Ys[rows, paste(traits, l, sep = ":")] <- Y[rows, ]
  }
  fit <- fit_G(Ys, X = X, A = A, engine = engine, missing = "model", ...)
  G <- fit$G

  # block-diagonal E: within each environment, residual covariance minus
  # c G of that environment's block
  A <- .indri_match_A(A, rownames(Y))
  E <- matrix(0, length(states), length(states),
              dimnames = list(states, states))
  for (l in levels) {
    rows <- env == l
    s <- paste(traits, l, sep = ":")
    E[s, s] <- .indri_residual_E(Y[rows, , drop = FALSE],
                                 X[rows, , drop = FALSE],
                                 unname(G[s, s]),
                                 A[rows, rows, drop = FALSE])
  }

  # each trait's genetic correlation between pairs of environments
  pairs <- utils::combn(levels, 2L)
  corr_between <- matrix(NA_real_, length(traits), ncol(pairs),
                         dimnames = list(traits, apply(pairs, 2L, paste,
                                                       collapse = "-")))
  for (k in seq_along(traits)) for (i in seq_len(ncol(pairs))) {
    a <- paste(traits[k], pairs[1L, i], sep = ":")
    b <- paste(traits[k], pairs[2L, i], sep = ":")
    corr_between[k, i] <- G[a, b] / sqrt(G[a, a] * G[b, b])
  }

  structure(list(G = G, E = E, corr_between = corr_between, levels = levels,
                 traits = traits, states = states,
                 env = stats::setNames(env, rownames(Y)), fit = fit),
            class = "indri_fit_env")
}

#' @export
print.indri_fit_env <- function(x, ...) {
  cat("indri environment-specific G (experimental):", length(x$traits),
      "traits in", length(x$levels), "environments (",
      paste(x$levels, collapse = ", "), ")\n")
  cat("genetic correlation of each trait between environments:\n")
  print(round(x$corr_between, 3))
  invisible(x)
}
