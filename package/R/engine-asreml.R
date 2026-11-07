# asreml engine: one multivariate animal model per trait subset, fitted with
# asreml::asreml() in wide format. REML point estimates with standard errors
# from the variance-component table, for inverse-variance assembly weighting.
#
# ASReml-R is licensed and not on CRAN; the engine errors with asreml's own
# message when the package is missing or unlicensed.

.indri_engine_asreml <- function() {
  list(
    name = "asreml",
    capabilities = list(
      draws = FALSE,
      se = TRUE,
      a_form = "sparse_inverse",
      missing = FALSE,
      # pairs, once each; see .indri_engine_defaults()
      default_s = 2L,
      default_cooccurrence = 1L,
      # each worker checks out a licence seat; a standard licence has two
      default_cores = 2L
    ),
    fit = .indri_asreml_fit
  )
}

.indri_asreml_fit <- function(Y, X, a_forms, control = list()) {
  if (!requireNamespace("asreml", quietly = TRUE))
    stop("The `asreml` engine needs the asreml package. It is licensed and ",
         "not on CRAN; install it separately, or choose a different engine.",
         call. = FALSE)

  N <- nrow(Y)
  s <- ncol(Y)
  ids <- rownames(Y)
  traits <- colnames(Y)
  if (is.null(ids)) ids <- as.character(seq_len(N))
  if (is.null(traits)) traits <- paste0("trait", seq_len(s))
  traits <- make.names(traits, unique = TRUE)

  # asreml's vm() matches the inverse relationship matrix to the `id` levels
  # through its rowNames attribute.
  Ainv <- a_forms$sparse_inverse
  attr(Ainv, "rowNames") <- ids

  # Wide format: one column per trait. Fixed effects are trait-specific
  # intercepts (when X has an intercept column) plus trait-specific slopes.
  x_names <- make.names(colnames(X), unique = TRUE)
  constant <- apply(X, 2L, function(v) all(v == v[1L]))
  covars <- which(!constant)

  dat <- data.frame(id = factor(ids))
  if (length(covars)) {
    for (j in covars) dat[[x_names[j]]] <- X[, j]
  }
  for (k in seq_len(s)) dat[[traits[k]]] <- Y[, k]

  terms <- character(0)
  if (any(constant)) terms <- c(terms, "trait")
  if (length(covars)) terms <- c(terms, paste0("trait:", x_names[covars]))
  if (!length(terms)) {
    # no intercept: trait-specific slopes only, and no global intercept
    terms <- c("0 +", paste0("trait:", x_names))
  }

  fixed <- stats::as.formula(
    paste("cbind(", paste(traits, collapse = ", "), ") ~ ",
          paste(terms, collapse = " + ")),
    env = environment()
  )
  # Starting values: the identity unless `control$init = list(G =, E =)` gives
  # covariances on the standardised scale.
  init <- .indri_check_init(control$init, s)
  init_G <- .indri_asreml_us_init(s, init$G)
  init_E <- .indri_asreml_us_init(s, init$E)
  random <- stats::as.formula(
    paste("~ vm(id, Ainv):us(trait, init = c(",
          paste(init_G, collapse = ", "), "))"),
    env = environment()
  )
  rcov <- stats::as.formula(
    paste("~ units:us(trait, init = c(",
          paste(init_E, collapse = ", "), "))"),
    env = environment()
  )

  # control carries the core budget fit_G() injects, which asreml() does not
  # understand. `maxit` (AI iterations per pass, asreml's default 13) and
  # `updates` (extra passes, default 1) are the convergence controls.
  maxit <- control$maxit %||% 13L
  updates <- control$updates %||% 1L
  started <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(
    asreml::asreml(fixed = fixed, random = random, rcov = rcov, data = dat,
                   maxit = maxit)
  )
  # The first pass can stop early; the updates finish it.
  for (u in seq_len(updates)) {
    if (isTRUE(fit$converge) && u > 1L) break
    fit <- suppressWarnings(stats::update(fit, maxit = maxit))
  }
  elapsed <- proc.time()[["elapsed"]] - started

  vc <- summary(fit)$varcomp
  gi <- grep("^vm\\(", rownames(vc))
  ei <- grep("^units.*trait!trait", rownames(vc))

  list(
    G = .indri_lower2mat(vc[gi, "component"], s),
    E = .indri_lower2mat(vc[ei, "component"], s),
    draws = NULL,
    se = .indri_lower2mat(vc[gi, "std.error"], s),
    converged = isTRUE(fit$converge),
    loglik = as.numeric(summary(fit)$loglik),
    time = elapsed
  )
}

# us(trait, init = ...) start values on the lower triangle, in asreml's
# row-by-row order: the identity, or the lower triangle of `m`.
.indri_asreml_us_init <- function(k, m = NULL) {
  rows <- rep(seq_len(k), seq_len(k))
  cols <- sequence(seq_len(k))
  if (!is.null(m)) return(signif(m[cbind(rows, cols)], 8))
  out <- numeric(length(rows))
  out[rows == cols] <- 1
  out
}

# asreml reports us() parameters as the lower triangle, row by row.
.indri_lower2mat <- function(x, k) {
  m <- matrix(0, k, k)
  rows <- rep(seq_len(k), seq_len(k))
  cols <- sequence(seq_len(k))
  m[cbind(rows, cols)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}
