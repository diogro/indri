# sommer engine: one multivariate animal model per trait subset, fitted with
# sommer::mmes() in long format. REML estimates with standard errors and the
# joint sampling covariance of the G entries (`vcov`), from which fit_G()
# draws the block's uncertainty.

.indri_engine_sommer <- function() {
  list(
    name = "sommer",
    capabilities = list(
      draws = FALSE,
      se = TRUE,
      # joint sampling covariance of the G entries, from which fit_G() draws
      vcov = TRUE,
      a_form = "inverse",
      # gaps in Y are fine: a missing value is one missing observation, and the
      # REML fit uses everything else the individual has
      missing = TRUE,
      # 4-trait blocks, each pair in at least one; see
      # .indri_engine_defaults()
      default_s = 3L,
      default_cooccurrence = 1L
    ),
    fit = .indri_sommer_fit
  )
}

.indri_sommer_fit <- function(Y, X, a_forms, control = list()) {
  if (!requireNamespace("sommer", quietly = TRUE))
    stop("The `sommer` engine needs the sommer package. Install it, or ",
         "choose a different engine.", call. = FALSE)

  N <- nrow(Y)
  K <- ncol(Y)
  ids <- rownames(Y)
  traits <- colnames(Y)
  if (is.null(ids)) ids <- as.character(seq_len(N))
  if (is.null(traits)) traits <- paste0("trait", seq_len(K))

  # mmes() matches Gu against the `id` levels, so the inverse relationship
  # matrix has to carry the individual names.
  Ainv <- a_forms$inverse
  dimnames(Ainv) <- list(ids, ids)

  # Fixed effects: trait-specific intercepts (when X has an intercept column)
  # plus trait-specific slopes for every other column. Long format below.
  x_names <- make.names(colnames(X), unique = TRUE)
  constant <- apply(X, 2L, function(v) all(v == v[1L]))
  covars <- which(!constant)

  long <- data.frame(
    id = factor(rep(ids, K)),
    # levels in column order: factor() would sort them alphabetically, and
    # sommer returns its matrices in level order
    trait = factor(rep(traits, each = N), levels = traits),
    value = as.vector(Y),
    record = factor(rep(seq_len(N), K)),
    stringsAsFactors = FALSE
  )
  if (length(covars)) {
    for (j in covars) long[[x_names[j]]] <- rep(X[, j], K)
  }

  # A single trait is an ordinary univariate animal model: the formulas drop
  # the `trait` factor (and its usm() structure) altogether.
  single <- K == 1L
  terms <- character(0)
  if (any(constant)) terms <- c(terms, if (single) "1" else "trait")
  if (length(covars))
    terms <- c(terms, if (single) x_names[covars]
                      else paste0("trait:", x_names[covars]))
  if (!length(terms))
    stop("`X` has no columns; the sommer engine needs a design matrix.",
         call. = FALSE)

  fixed <- stats::as.formula(
    paste("value ~ 0 +", paste(terms, collapse = " + ")),
    env = environment()
  )

  # sommer evaluates the structure calls in the formula's environment, and its
  # internal evaluator re-emits them with bare names (vsm/usm/ism), so bind the
  # helpers there and let `Ainv` resolve through the parent frame.
  .env <- new.env(parent = environment())
  .env$vsm <- sommer::vsm
  .env$usm <- sommer::usm
  .env$ism <- sommer::ism
  # Optional starting covariances (`control$init = list(G =, E =)`, on the
  # standardised scale): sommer takes their shape through usm(theta =) and
  # finds the overall scale from the data.
  init <- .indri_check_init(control$init, K)
  .env$G0 <- init$G
  .env$E0 <- init$E
  if (single) {
    random <- stats::as.formula("~ vsm(ism(id), Gu = Ainv)", env = .env)
    rcov <- stats::as.formula("~ vsm(ism(record))", env = .env)
  } else {
    random <- stats::as.formula(paste0(
      "~ vsm(usm(trait", if (!is.null(init$G)) ", theta = G0", "), ism(id), ",
      "Gu = Ainv)"), env = .env)
    rcov <- stats::as.formula(paste0(
      "~ vsm(usm(trait", if (!is.null(init$E)) ", theta = E0", "), ",
      "ism(record))"), env = .env)
  }

  # .indri_core_plan() adds a `cores` budget that sommer's control object does
  # not understand; keep only the entries mmes_control() accepts.
  control_names <- names(formals(sommer::mmes_control))
  ctrl <- utils::modifyList(
    list(diagnostics = list(verbose = FALSE)),
    control[intersect(names(control), control_names)]
  )
  mmes_control <- do.call(sommer::mmes_control, ctrl)

  started <- proc.time()[["elapsed"]]
  fit <- sommer::mmes(fixed, random = random, rcov = rcov, data = long,
                      control = mmes_control)
  elapsed <- proc.time()[["elapsed"]] - started

  # covmatrix_mmes() term 1 is the animal effect (G), term 2 the residual (E).
  G <- sommer::covmatrix_mmes(fit, 1, se = TRUE)
  E <- sommer::covmatrix_mmes(fit, 2)

  vcov <- .indri_sommer_vcov(fit, traits, single)

  llik <- utils::tail(fit$llik, 1)
  if (length(llik) != 1L || !is.finite(llik)) llik <- NA_real_

  # Return every matrix in the column order of Y. sommer labels them by trait
  # level; reorder by name in case its order ever differs from the levels.
  in_order <- function(m) {
    m <- as.matrix(m)
    if (!single && !is.null(rownames(m)) && setequal(rownames(m), traits))
      m <- m[traits, traits, drop = FALSE]
    m
  }
  list(
    G = in_order(G$covariance),
    E = in_order(E$covariance),
    draws = NULL,
    se = in_order(G$covariance.se),
    vcov = vcov,
    converged = isTRUE(fit$convergence),
    loglik = as.numeric(llik),
    time = elapsed
  )
}

# Sampling covariance of the block's G entries, in .indri_vech_index() order,
# from sommer's natural-scale parameter table. Parameters are matched by name
# ("variance[a]", "covariance[b,a]"), so the result follows the columns of Y
# whatever order sommer keeps them in. NULL when sommer reports none.
.indri_sommer_vcov <- function(fit, traits, single) {
  native <- fit$covParNative
  V <- attr(native, "vcov")
  if (is.null(native) || is.null(V)) return(NULL)
  rows <- which(native$term == names(fit$theta)[1L])
  k <- length(traits)
  index <- .indri_vech_index(k)
  if (single) {
    if (length(rows) != 1L) return(NULL)
    out <- V[rows, rows, drop = FALSE]
  } else {
    parm <- native$parameter[rows]
    pos <- vapply(parm, function(p) {
      inside <- sub("^[a-z]+\\[(.*)\\]$", "\\1", p)
      who <- match(strsplit(inside, ",", fixed = TRUE)[[1L]], traits)
      if (length(who) == 1L) who <- c(who, who)
      if (anyNA(who) || length(who) != 2L) return(NA_integer_)
      which(index[, 1L] == max(who) & index[, 2L] == min(who))
    }, integer(1))
    if (anyNA(pos) || length(pos) != nrow(index) || anyDuplicated(pos))
      return(NULL)
    out <- matrix(NA_real_, nrow(index), nrow(index))
    out[pos, pos] <- V[rows, rows]
  }
  if (anyNA(out) || any(!is.finite(out))) return(NULL)
  unname(out)
}
