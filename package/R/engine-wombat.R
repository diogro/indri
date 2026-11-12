# wombat engine: one multivariate animal model per trait subset, fitted by
# Karin Meyer's WOMBAT, a standalone REML program (a binary, not an R
# package). Each block runs in its own temporary directory: indri writes the
# data, the inverse relationship matrix and a parameter file, runs WOMBAT, and
# reads back the estimates (BestPoint), the average information matrix for the
# covariance components (AvInfoCovs), from which come `se` and `vcov`, and the
# iterate log (Iterates), from which comes `converged`.
#
# WOMBAT builds exist for Linux and Windows only. The binary is found through
# `engine_args$wombat`, then `options(indri.wombat =)`, then the PATH.
#
# Gaps in Y are left out of the data file (WOMBAT reads one record per trait
# and individual). Residual covariances between traits that no individual has
# both of, such as one trait in two environments (fit_G_env()), are not
# identified; they start at zero, WOMBAT leaves them out of the estimation and
# they come back as zero.

.indri_engine_wombat <- function() {
  # The option is read here, in the main process: PSOCK workers do not
  # inherit options().
  option <- getOption("indri.wombat")
  list(
    name = "wombat",
    capabilities = list(
      draws = FALSE,
      se = TRUE,
      # joint sampling covariance of the G entries, from which fit_G() draws
      vcov = TRUE,
      a_form = "sparse_inverse",
      missing = TRUE,
      # 4-trait blocks, each pair in at least one, as for the other REML
      # engines; see .indri_engine_defaults()
      default_s = 3L,
      default_cooccurrence = 1L
    ),
    fit = function(Y, X, a_forms, control = list()) {
      if (is.null(control[["wombat"]])) control[["wombat"]] <- option
      .indri_wombat_fit(Y, X, a_forms, control)
    }
  )
}

# The WOMBAT executable, or an error naming where indri looked.
.indri_wombat_binary <- function(path = NULL) {
  path <- path %||% getOption("indri.wombat") %||% unname(Sys.which("wombat"))
  if (!nzchar(path) || !file.exists(path))
    stop("The `wombat` engine needs the WOMBAT program (Linux or Windows; ",
         "free for academic use from Karin Meyer). Put it on the PATH, or ",
         "give its location with `engine_args = list(wombat = <path>)` or ",
         "`options(indri.wombat = <path>)`.", call. = FALSE)
  normalizePath(path)
}

.indri_wombat_fit <- function(Y, X, a_forms, control = list()) {
  binary <- .indri_wombat_binary(control[["wombat"]])

  N <- nrow(Y)
  s <- ncol(Y)
  ids <- rownames(Y) %||% as.character(seq_len(N))
  Ainv <- a_forms$sparse_inverse
  if (!is.null(attr(Ainv, "rowNames")) &&
      !identical(as.character(attr(Ainv, "rowNames")), as.character(ids)))
    stop("The inverse relationship matrix is not in the row order of `Y`.",
         call. = FALSE)

  # `wombat_dir` keeps the block's files there (for inspection); by default
  # they go to a temporary directory that is removed afterwards
  if (is.null(control[["wombat_dir"]])) {
    dir <- tempfile("indri-wombat-")
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  } else {
    dir <- control[["wombat_dir"]]
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir <- normalizePath(dir)

  covars <- .indri_wombat_covariables(X)
  .indri_wombat_write_gin(file.path(dir, "animal.gin"),
                          file.path(dir, "animal.codes"), Ainv, N)
  .indri_wombat_write_data(file.path(dir, "data.dat"), Y, covars)

  init <- .indri_check_init(control[["init"]], s)
  observed <- !is.na(Y)
  together <- crossprod(observed) > 0
  # standardised traits: half the variance genetic, half residual. (WOMBAT
  # pauses to warn when a starting residual variance reaches the raw variance,
  # as the identity for both would.)
  start_G <- init$G %||% diag(0.5, s)
  start_E <- init$E %||% diag(0.5, s)
  # WOMBAT holds a covariance whose starting value is zero at zero (measured:
  # a bivariate residual covariance stayed at 0 from a zero start and matched
  # sommer from 0.05). Start the identified ones away from zero; 0.25 / s
  # keeps the default start positive definite at any s.
  nudge <- 0.25 / s
  off <- row(start_G) != col(start_G)
  start_G[off & start_G == 0] <- nudge
  start_E[off & together & start_E == 0] <- nudge
  # residual covariances of traits never recorded together are not
  # identified: zero, and so not estimated
  start_E[!together] <- 0
  # `wombat_rank`: reduced-rank G (WOMBAT's PC analysis), for G whose
  # REML estimate sits on the positive-definite boundary
  rank <- control[["wombat_rank"]]
  if (!is.null(rank)) {
    rank <- as.integer(rank)
    if (length(rank) != 1L || is.na(rank) || rank < 1L)
      stop("`wombat_rank` must be a positive integer.", call. = FALSE)
    if (rank >= s) rank <- NULL
  }
  if (!is.null(rank))
    .indri_wombat_write_matrix(file.path(dir, "animal.matrix"), Ainv, N)
  .indri_wombat_write_par(file.path(dir, "indri.par"), s, N, ncol(covars),
                          start_G, start_E, rank)

  threads <- max(1L, as.integer(control[["cores"]] %||% 1L))
  # --batch: no programmed pauses for interactive confirmation
  args <- c("--batch", paste0("--threads", threads), control[["wombat_args"]],
            "indri.par")
  started <- proc.time()[["elapsed"]]
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE, after = FALSE)
  # WOMBAT's OpenMP/MKL threads follow --threads; keep the environment from
  # asking for more
  env <- paste0(c("OMP_NUM_THREADS=", "MKL_NUM_THREADS="), 2L * threads)
  if (.Platform$OS.type == "unix") {
    # The Linux build keeps large arrays on the stack: with the usual 8 MB
    # limit it segfaults entering AI-REML once a block has about 20 traits
    # (measured at 20 character states). Lift the soft limit for its shell.
    command <- paste("ulimit -s unlimited 2>/dev/null || ulimit -s hard",
                     "2>/dev/null; exec", paste(shQuote(c(binary, args)),
                                                collapse = " "))
    log <- suppressWarnings(system2("/bin/sh", c("-c", shQuote(command)),
                                    stdout = TRUE, stderr = TRUE, env = env))
  } else {
    log <- suppressWarnings(system2(binary, args, stdout = TRUE,
                                    stderr = TRUE, env = env))
  }
  status <- attr(log, "status") %||% 0L
  setwd(old)
  elapsed <- proc.time()[["elapsed"]] - started

  # a crash can leave a BestPoint from the last EM round: never return it
  best <- file.path(dir, "BestPoint")
  if (status != 0L || !file.exists(best))
    stop("WOMBAT ", if (status != 0L) paste("stopped with status", status)
         else "produced no estimates", ". Its last output:\n",
         paste(utils::tail(log, 15L), collapse = "\n"), call. = FALSE)
  point <- .indri_wombat_read_bestpoint(best, s)
  vcov <- .indri_wombat_read_vcov(file.path(dir, "AvInfoCovs"), s,
                                  point$loglik)
  se <- matrix(NA_real_, s, s)
  if (!is.null(vcov)) {
    index <- .indri_vech_index(s)
    se[index] <- sqrt(pmax(diag(vcov), 0))
    se[index[, 2:1, drop = FALSE]] <- se[index]
  }

  traits <- colnames(Y) %||% paste0("trait", seq_len(s))
  names <- list(traits, traits)
  list(
    G = structure(point$G, dimnames = names),
    E = structure(point$E, dimnames = names),
    draws = NULL,
    se = structure(se, dimnames = names),
    vcov = vcov,
    converged = .indri_wombat_converged(file.path(dir, "Iterates")),
    loglik = point$loglik,
    time = elapsed
  )
}

# The columns of X that WOMBAT fits as covariables: indri's X always spans a
# constant (.indri_check_design()), which becomes WOMBAT's single-level fixed
# effect `mu`, and the covariables are the columns that add to it, chosen
# greedily so that mu plus covariables has full rank.
.indri_wombat_covariables <- function(X) {
  X <- as.matrix(X)
  keep <- integer(0)
  basis <- matrix(1, nrow(X), 1L)
  for (j in seq_len(ncol(X))) {
    trial <- cbind(basis, X[, j])
    if (qr(trial)$rank > ncol(basis)) {
      basis <- trial
      keep <- c(keep, j)
    }
  }
  X[, keep, drop = FALSE]
}

# animal.gin: the log determinant of A (0: a constant in log L, not needed
# here), then the lower triangle of A^-1 as `column row value` with
# column <= row. animal.codes: level i has code i, so that individuals with no
# record in this block keep their place in A^-1.
.indri_wombat_write_gin <- function(gin, codes, Ainv, N) {
  lines <- sprintf("%d %d %.15g", as.integer(Ainv[, 2L]),
                   as.integer(Ainv[, 1L]), as.numeric(Ainv[, 3L]))
  writeLines(c(sprintf("%.15g", .indri_wombat_logdet(Ainv, N)), lines), gin)
  writeLines(sprintf("%d %d", seq_len(N), seq_len(N)), codes)
}

# log |A|, the GIN header: a constant of the log likelihood in a full-rank
# fit, but in a reduced-rank fit it enters once per principal component, so
# log L is comparable across ranks only with the true value. Computed densely
# (dense Cholesky of the N x N inverse) up to 10000 individuals; beyond that
# 0, and log L is then defined up to a rank-dependent constant.
.indri_wombat_logdet <- function(Ainv, N) {
  if (N > 10000L) return(0)
  M <- .indri_wombat_dense(Ainv, N)
  -2 * sum(log(diag(chol(M))))
}

.indri_wombat_dense <- function(Ainv, N) {
  M <- matrix(0, N, N)
  M[cbind(as.integer(Ainv[, 1L]), as.integer(Ainv[, 2L]))] <- Ainv[, 3L]
  M[cbind(as.integer(Ainv[, 2L]), as.integer(Ainv[, 1L]))] <- Ainv[, 3L]
  M
}

# The relationship matrix itself (animal.matrix, next to animal.gin): WOMBAT's
# reduced-rank AI-REML with a GIN effect multiplies by it in every iterate and
# stops without it. indri hands the engine the sparse inverse only, so it is
# inverted here (dense, N x N), and the lower triangle's non-zero elements
# are written as `column row value`, without the GIN's determinant line.
.indri_wombat_write_matrix <- function(path, Ainv, N) {
  A <- chol2inv(chol(.indri_wombat_dense(Ainv, N)))
  index <- which(lower.tri(A, diag = TRUE) & abs(A) > 1e-14, arr.ind = TRUE)
  writeLines(sprintf("%d %d %.15g", index[, 2L], index[, 1L], A[index]), path)
}

# One record per trait and individual, sorted by individual then trait:
# trait number, individual, mu (always 1), covariables, response. Gaps are
# absent records.
.indri_wombat_write_data <- function(path, Y, covars) {
  s <- ncol(Y)
  rows <- which(!is.na(t(Y)), arr.ind = TRUE)   # trait within individual
  trait <- rows[, 1L]
  ind <- rows[, 2L]
  fields <- cbind(
    sprintf("%d", trait), sprintf("%d", ind), "1",
    if (ncol(covars))
      matrix(sprintf("%.12g", covars[ind, , drop = FALSE]), length(ind)),
    sprintf("%.12g", Y[cbind(ind, trait)]))
  writeLines(apply(fields, 1L, paste, collapse = " "), path)
}

# The parameter file. Every trait has the same model: the mean (`mu`), the
# covariables and an animal effect with covariance proportional to the
# supplied inverse (GIN). Starting values are upper triangles, one per line
# (WOMBAT's lines are limited to 78 characters).
.indri_wombat_write_par <- function(path, s, N, n_cov, G, E, rank = NULL) {
  # `rank` < s: a reduced-rank (PC) analysis, fitting the leading `rank`
  # principal components of G; the residual stays full rank
  pc <- !is.null(rank) && rank < s
  cov_names <- if (n_cov) paste0("x", seq_len(n_cov)) else character(0)
  columns <- function(k) c(
    sprintf("TR%d traitno %d", k, s), sprintf("TR%d animal %d", k, N),
    sprintf("TR%d mu 1", k), sprintf("TR%d %s", k, cov_names),
    sprintf("TR%d y%d", k, k))
  upper <- function(m) {
    index <- which(upper.tri(m, diag = TRUE), arr.ind = TRUE)
    index <- index[order(index[, 1L], index[, 2L]), , drop = FALSE]
    sprintf("  %.10g", m[index])
  }
  lines <- c(
    "COM indri block",
    # compiled-in limits (`wombat --limits`: 42 traits, 19 covariables) are
    # raised by OVERRIDE lines, which must come before the data layout
    if (s > 42L) sprintf("OVER MAXTRAIT %d", s),
    if (n_cov > 19L) sprintf("OVER MAXCOV %d", n_cov),
    sprintf("ANAL MUV %s%d", if (pc) "PC " else "", s),
    "DATA data.dat",
    unlist(lapply(seq_len(s), columns)),
    "END",
    "MODEL",
    "  FIX mu",
    if (n_cov) sprintf("  COV %s(1)", cov_names),
    "  RAN animal GIN",
    # WOMBAT takes the subject of a record from the first genetic (NRM)
    # effect; with GIN it has to be told, or records of one individual count
    # as different subjects and the residual covariances stay at zero
    "  SUBJ animal",
    sprintf("  TR y%d %d", seq_len(s), seq_len(s)),
    "END MOD",
    sprintf("VAR residual %d%s", s, if (pc) sprintf(" %d", s) else ""),
    upper(E),
    sprintf("VAR animal %d%s", s, if (pc) sprintf(" %d", rank) else ""),
    upper(G))
  writeLines(lines, path)
}

# BestPoint: a header line (log L, the number of parameters, the number of
# records), then each covariance matrix as its upper triangle, row by row,
# residual first.
.indri_wombat_read_bestpoint <- function(path, s) {
  lines <- readLines(path)
  loglik <- as.numeric(strsplit(trimws(lines[1L]), "\\s+")[[1L]][1L])
  x <- scan(text = lines[-1L], quiet = TRUE)
  m <- s * (s + 1L) / 2L
  if (length(x) != 2L * m)
    stop("WOMBAT's BestPoint has ", length(x), " covariance components; ",
         "expected ", 2L * m, ".", call. = FALSE)
  list(loglik = loglik,
       E = .indri_wombat_upper2mat(x[seq_len(m)], s),
       G = .indri_wombat_upper2mat(x[m + seq_len(m)], s))
}

.indri_wombat_upper2mat <- function(x, s) {
  m <- matrix(0, s, s)
  index <- which(upper.tri(m, diag = TRUE), arr.ind = TRUE)
  index <- index[order(index[, 1L], index[, 2L]), , drop = FALSE]
  m[index] <- x
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

# Sampling covariance of the G entries in .indri_vech_index() order: the
# inverse of the average information for the covariance components
# (AvInfoCovs: log L, then `row column value` for the upper triangle; residual
# components first, then G, each as an upper triangle row by row). Components
# WOMBAT did not estimate have zero rows and are dropped before inverting.
# NULL when the file is missing, belongs to another point or is singular.
.indri_wombat_read_vcov <- function(path, s, loglik) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    lines <- readLines(path)
    at <- as.numeric(strsplit(trimws(lines[1L]), "\\s+")[[1L]][1L])
    if (!is.finite(at) || abs(at - loglik) > 1e-4 * max(1, abs(loglik)))
      return(NULL)
    x <- matrix(scan(text = lines[-1L], quiet = TRUE), ncol = 3L,
                byrow = TRUE)
    m <- s * (s + 1L) / 2L
    p <- max(x[, 1:2], 2L * m)
    info <- matrix(0, p, p)
    info[x[, 1:2]] <- x[, 3L]
    info[x[, 2:1]] <- x[, 3L]
    estimated <- which(diag(info) > 0)
    inverse <- matrix(NA_real_, p, p)
    inverse[estimated, estimated] <- solve(info[estimated, estimated])
    # G's position j in the upper triangle (row by row) of the vech pair
    # (row, column) with row >= column, i.e. upper-triangle cell (column, row)
    upper <- which(upper.tri(diag(s), diag = TRUE), arr.ind = TRUE)
    upper <- upper[order(upper[, 1L], upper[, 2L]), , drop = FALSE]
    vech <- .indri_vech_index(s)
    order <- match(paste(vech[, 2L], vech[, 1L]),
                   paste(upper[, 1L], upper[, 2L]))
    v <- inverse[m + order, m + order, drop = FALSE]
    if (anyNA(v) || any(!is.finite(v))) NULL else unname(v)
  }, error = function(e) NULL)
}

# WOMBAT ends the Iterates file with its verdict, e.g. "AIREML seems to have
# converged: Change in log L < 0.5E-03" (its default AI criterion).
.indri_wombat_converged <- function(path) {
  if (!file.exists(path)) return(FALSE)
  any(grepl("seems to have converged", readLines(path), fixed = TRUE))
}
