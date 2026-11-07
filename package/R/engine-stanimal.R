# stanimal engine: one multivariate animal model per trait subset, fitted with
# stanimal::lmm_animal(). Returns posterior draws, so the assembly can carry
# block-level uncertainty through to the full matrix.

.indri_engine_stanimal <- function() {
  list(
    name = "stanimal",
    capabilities = list(
      draws = TRUE,
      se = FALSE,
      # the eigendecomposition is shared by every block when stanimal can take
      # it, so the O(N^3) step runs once per fit_G() call instead of once per
      # block (and, before stanimal took it as data, once per chain)
      a_form = if (.indri_stanimal_takes_eigen()) c("dense", "eigen") else "dense",
      missing = TRUE,
      # what lmm_animal() uses when the caller does not say otherwise, so that
      # fit_G() can divide a core budget without oversubscribing
      default_chains = 4L
    ),
    fit = .indri_stanimal_fit
  )
}

.indri_stanimal_fit <- function(Y, X, a_forms, control = list()) {
  if (!requireNamespace("stanimal", quietly = TRUE))
    stop("The `stanimal` engine needs the stanimal package. Install it, or ",
         "choose a different engine.", call. = FALSE)

  args <- utils::modifyList(list(chains = 4L, iter = 3000L), control)
  # a_forms carries `eigen` only when the engine asked for it, which it does
  # only when stanimal takes it (decided once, in the main process: a PSOCK
  # worker may not see the same helpers)
  if (!is.null(a_forms$eigen)) args$eigen_A <- a_forms$eigen
  fit <- do.call(stanimal::lmm_animal,
                 c(list(Y = Y, X = X, A = a_forms$dense), args))

  draws <- stanimal::fit_draws(fit)
  G_draws <- draws$G
  E_draws <- draws$E

  # The CmdStanMCMC object is deliberately not kept: it carries the draws of
  # every parameter plus handles on its output files, and a fit keeps one per
  # subset. Everything the assembly needs is extracted here.
  list(
    G = apply(G_draws, c(1L, 2L), mean),
    E = apply(E_draws, c(1L, 2L), mean),
    draws = G_draws,
    se = NULL,
    converged = .indri_stanimal_rhat(fit),
    n_draws = dim(G_draws)[3L]
  )
}

# Max Rhat over the covariance entries, when the backend reports it. rstan fits
# go through stanimal's draw extraction, which drops the chain structure, so the
# number is only available for cmdstanr.
.indri_stanimal_rhat <- function(fit) {
  if (!inherits(fit, "CmdStanMCMC")) return(NA_real_)
  tryCatch({
    draws <- fit$draws(variables = "G")
    value <- max(posterior::rhat(draws), na.rm = TRUE)
    if (is.finite(value)) value else NA_real_
  }, error = function(e) NA_real_)
}

# stanimal accepts a precomputed eigendecomposition of A from the version that
# takes it as model data; older versions decompose A themselves.
.indri_stanimal_takes_eigen <- function() {
  requireNamespace("stanimal", quietly = TRUE) &&
    "eigen_A" %in% names(formals(stanimal::lmm_animal))
}
