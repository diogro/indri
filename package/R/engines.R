#' Block engines
#'
#' The per-subset fit is a plugin. An engine is a list with a `name`, a
#' `capabilities` list and a `fit` function with the signature
#' `fit(Y, X, a_forms, control)` returning at least `G` and `E` (on the
#' standardised scale, in the column order of `Y`), plus either posterior
#' `draws` (an `s` by `s` by `D` array; stanimal) or `se` and `vcov`, the joint
#' sampling covariance of the lower-triangle G entries (the REML engines), from
#' which [fit_G()] draws.
#'
#' LiMMBo delegates its subset fits to limix; indri delegates to stanimal,
#' sommer, ASReml-R or the WOMBAT program, because the right engine depends on
#' the number of individuals, on whether the data have gaps, on whether a
#' licence is available, and on whether the caller wants a posterior or a
#' likelihood. WOMBAT is a standalone REML program (Linux and Windows, free
#' for academic use); the `"wombat"` engine finds it through
#' `engine_args = list(wombat = <path>)`, `options(indri.wombat = <path>)` or
#' the PATH, `engine_args = list(wombat_dir = <dir>)` keeps a block's
#' WOMBAT files for inspection, and `engine_args = list(wombat_rank = q)`
#' fits G at reduced rank q (the residual stays full rank).
#'
#' @param name Engine name.
#'
#' @return An engine list.
#'
#' @examples
#' # Inspect the registry without fitting a model or loading an engine.
#' indri_engines()
#' engine <- indri_engine("sommer")
#' engine$name
#' engine$capabilities
#'
#' @export
indri_engine <- function(name = "stanimal") {
  if (!is.character(name) || length(name) != 1L || is.na(name))
    stop("`engine` must be a single engine name.", call. = FALSE)
  switch(name,
    stanimal = .indri_engine_stanimal(),
    sommer = .indri_engine_sommer(),
    asreml = .indri_engine_asreml(),
    wombat = .indri_engine_wombat(),
    stop("Unknown engine `", name, "`. Implemented: ",
         paste(indri_engines(), collapse = ", "), ".", call. = FALSE)
  )
}

#' @rdname indri_engine
#' @return `indri_engines()` returns the names of the implemented engines.
#' @export
indri_engines <- function() {
  c("stanimal", "sommer", "asreml", "wombat")
}

# Starting covariances for a block fit, `control$init = list(G =, E =)` on the
# standardised scale. Either may be NULL; a supplied one must be a symmetric
# positive-definite k x k matrix.
.indri_check_init <- function(init, k) {
  if (is.null(init)) return(list(G = NULL, E = NULL))
  if (!is.list(init))
    stop("`init` must be a list with elements `G` and/or `E`.", call. = FALSE)
  for (part in c("G", "E")) {
    m <- init[[part]]
    if (is.null(m)) next
    m <- as.matrix(m)
    if (!identical(dim(m), c(k, k)) || anyNA(m) ||
        max(abs(m - t(m))) > 1e-8 ||
        min(eigen((m + t(m)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values) <= 0)
      stop("`init$", part, "` must be a symmetric positive-definite ", k,
           " x ", k, " matrix.", call. = FALSE)
    init[[part]] <- unname((m + t(m)) / 2)
  }
  list(G = init$G, E = init$E)
}

# Lower-triangle positions of a k x k matrix, row by row: (1,1), (2,1), (2,2),
# (3,1), ... The order of the REML engines' `vcov`.
.indri_vech_index <- function(k) {
  cbind(rep(seq_len(k), seq_len(k)), sequence(seq_len(k)))
}

# Normal draws of a block's G from its point estimate and the sampling
# covariance of its lower-triangle entries (the REML engines' `vcov`): the
# asymptotic REML sampling distribution, drawn jointly so that the entries
# keep their correlation. Without a usable `vcov` the block contributes its
# point estimate only (no spread). A draw near a boundary can be indefinite;
# .indri_corr_draws() copes with that.
.indri_normal_draws <- function(G, vcov, n_draws) {
  k <- nrow(G)
  out <- array(G, c(k, k, n_draws))
  if (is.null(vcov)) return(out)
  index <- .indri_vech_index(k)
  e <- eigen((vcov + t(vcov)) / 2, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(index))
  x <- G[index] + L %*% matrix(stats::rnorm(nrow(index) * n_draws),
                               nrow(index), n_draws)
  for (d in seq_len(n_draws)) {
    m <- matrix(0, k, k)
    m[index] <- x[, d]
    m[index[, 2:1, drop = FALSE]] <- x[, d]
    out[, , d] <- m
  }
  out
}
