#' Estimate the additive genetic covariance matrix
#'
#' Fits a multivariate animal model on many small trait subsets and reassembles
#' the subset estimates into one positive semi-definite matrix. The estimator is
#' LiMMBo's: the speed comes from giving each subset its own nuisance residual
#' covariance, so the many small fits replace one large one.
#'
#' Traits are standardised once, using every observed value, before the subsets
#' are drawn. Every engine therefore sees the same numbers and estimates the same
#' matrix on the same scale, whatever it does internally; the fixed effects in
#' `X` are handled inside each block fit. Assembly and the positive
#' semi-definite projection also run on that scale, and the assembled
#' covariances (`G`, `G_bar`, `G_draws`, `sigma_G`) are then returned in the
#' units of `Y`. The raw block fits in `blocks` stay on the standardised scale:
#' do not combine them directly with the assembled outputs (`y_sd` converts
#' between the two).
#'
#' @param Y Phenotype matrix, individuals in rows and traits in columns. Gaps
#'   (`NA`) are allowed with every engine; see `missing`.
#' @param X Fixed-effect design, one row per individual. Defaults to an
#'   intercept. It must have full column rank and be able to represent a
#'   constant (an intercept, or the indicators of every level of a factor).
#'   `stanimal::animal_design()` builds `Y`, `X` and `A` from a formula, a data
#'   frame and a pedigree: `d <- animal_design(cbind(t1, t2) ~ sex + age, dat,
#'   pedigree = ped); fit_G(d$Y, d$X, d$A)`.
#' @param A Relationship matrix, individuals in the same order as `Y`. When both
#'   have row names they are matched by name.
#' @param s Subset size parameter. Controls how many trait pairs each subset
#'   covers, not the number of traits in a subset; see [indri_subsets()].
#'   Defaults to `min(floor(sqrt(K)), 3)`, at least 2 and kept narrow enough
#'   that a block never covers every trait: blocks of up to four traits. Larger
#'   blocks shrank the stanimal correlations more (and undercovered) and gave
#'   the REML engines nothing for their extra cost.
#' @param min_cooccurrence Minimum number of distinct subsets every trait pair
#'   must appear in. The default depends on the engine: 3 for `"stanimal"`, 1
#'   for the REML engines, where three distinct contexts per pair measured no
#'   better than one (same accuracy and interval coverage at half the cost).
#'   Identical subsets are never repeated, since a pair estimated twice in the
#'   same trait context adds nothing, so the value is capped at what the block
#'   width allows (1 for two-trait blocks).
#' @param subsets Optional precomputed subset matrix, from [indri_subsets()].
#' @param engine Block engine: `"stanimal"` (Bayesian, posterior draws, models
#'   gaps in `Y`), `"sommer"`, `"asreml"` or `"wombat"` (REML, with
#'   uncertainty from each block's sampling covariance; much faster for large
#'   N). `"wombat"` runs the external WOMBAT program; see [indri_engine()].
#' @param engine_args List of arguments passed to the engine's fit function.
#' @param assembly `"ls"` fits the closest positive semi-definite matrix to all
#'   the subset estimates; `"average"` keeps the plain average, nudged onto the
#'   positive semi-definite cone.
#' @param weights `"none"` for LiMMBo's equal weights, or `"inverse_variance"`
#'   to weight each subset estimate by its reported uncertainty (per-entry
#'   standard error for REML engines, per-entry posterior variance for the
#'   draws-based stanimal engine). Weights apply to both the average and the
#'   projection.
#' @param missing What to do with gaps in `Y`. `"model"` (the default) keeps
#'   every observed value and lets the engine handle the gaps: stanimal
#'   integrates them out, and the REML engines use each individual's observed
#'   traits (a missing value is one missing record, not a missing individual).
#'   `"complete_block"` drops incomplete individuals once, before subsetting;
#'   it discards information (with 15% of values missing at random, a third of
#'   the individuals and about 20% wider standard errors).
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
#'   chains is overconfident by roughly the square root of a pair's coverage
#'   (40% coverage of the true correlations at nominal 95%, against 95-97% for
#'   the mixture at the default settings). See [indri_assemble_draws()].
#' @param n_draws Number of draws per block for the REML engines, which
#'   report a joint sampling covariance of each block's G entries (the inverse
#'   average-information matrix) instead of posterior draws: indri draws from
#'   that normal approximation and assembles the draws like stanimal's.
#'   Ignored for `"stanimal"`.
#' @param progress Show a progress bar of the block fits (blocks done, elapsed
#'   time, an estimate of the time left) on the console. Defaults to
#'   `getOption("indri.progress", interactive())`: on in an interactive
#'   session, off under `Rscript`. With parallel blocks the bar updates after
#'   each load-balanced round of about four blocks per worker. Engine output
#'   (for example CmdStan's own progress) is printed as before; silence it with
#'   `engine_args = list(refresh = 0, show_messages = FALSE)` for stanimal.
#' @param seed Optional seed for the subset construction. The engine's sampler
#'   seed is passed through `engine_args`.
#'
#' @return An object of class `indri_fit`. The estimate is `G`, in the units
#'   of `Y`, with `corrG` its correlation matrix; `G_bar` (the average before
#'   projection) and `sigma_G` are in the units of `Y` too. `E` is the
#'   residual covariance, for [indri_gwas()]: the covariance of the residuals
#'   of `Y` on `X` minus the assembled `G` (scaled for the relatedness in
#'   `A`), with eigenvalues floored at 1e-3 of the largest. `E_bar` is the
#'   average of the block residual covariances, for reference: at many traits
#'   it is not a covariance matrix, and a projection of it miscalibrates the
#'   association tests. `blocks`
#'   holds the raw standardised-scale block fits. `G_draws` and `corrG_draws`
#'   (traits x traits x draws) carry the assembled block uncertainty, from
#'   the posterior draws of `"stanimal"` blocks or from the normal sampling
#'   approximation of the REML blocks (see `n_draws`); use
#'   `corrG_draws` for intervals on correlations, since it mixes each block's
#'   own correlation draws and stays within \[-1, 1\]. Mixed draws are
#'   pointwise: their entries are drawn independently, so a single draw need
#'   not be positive semi-definite, and they support per-entry intervals, not
#'   functions of the whole matrix. Neither is a posterior over the joint
#'   matrix.
#'
#' @references
#' Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
#' a simple, scalable approach for linear mixed models in high-dimensional
#' genetic association studies. bioRxiv, preprint.
#' \doi{10.1101/255497}
#'
#' @examples
#' \donttest{
#' # This example fits six small Stan models and needs stanimal and CmdStan.
#' if (requireNamespace("stanimal", quietly = TRUE) &&
#'     requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.na(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   K <- 4L
#'   N <- 60L
#'   # Twelve families of five full siblings.
#'   A <- kronecker(diag(12), matrix(0.5, 5, 5))
#'   diag(A) <- 1
#'   G_true <- matrix(0.2, K, K)
#'   diag(G_true) <- 0.8
#'   E_true <- diag(0.6, K)
#'   X <- cbind(intercept = 1, group = rep(c(0, 1), length.out = N))
#'   genetic <- t(chol(A)) %*% matrix(rnorm(N * K), N, K) %*% chol(G_true)
#'   residual <- matrix(rnorm(N * K), N, K) %*% chol(E_true)
#'   Y <- genetic + residual + X[, "group"] %o% c(0.2, -0.2, 0.1, 0)
#'   colnames(Y) <- paste0("trait", seq_len(K))
#'
#'   # These short chains are for a quick example, not reliable inference.
#'   # For real work, use more chains and iterations and check diagnostics.
#'   fit <- fit_G(Y, X = X, A = A, subsets = t(utils::combn(K, 2)),
#'                engine = "stanimal", cores = 1,
#'                engine_args = list(chains = 2, iter = 400,
#'                                   refresh = 0, show_messages = FALSE))
#'   fit$corrG
#'   fit$coverage
#'   fit$diagnostics
#'
#'   # G is in the units of Y; the standardised-scale estimate, if needed:
#'   G_std <- diag(1 / fit$y_sd) %*% fit$G %*% diag(1 / fit$y_sd)
#'   fit$sigma_G
#'
#'   # Pointwise correlation spread, not a joint posterior interval.
#'   quantile(fit$corrG_draws[1, 2, ], c(0.025, 0.975))
#' }
#' }
#'
#' @export
fit_G <- function(Y, X = NULL, A, s = NULL, min_cooccurrence = NULL,
                  subsets = NULL, engine = "stanimal", engine_args = list(),
                  assembly = c("ls", "average"),
                  weights = c("none", "inverse_variance"),
                  missing = c("model", "complete_block"),
                  spread = c("mixture", "average"),
                  cores = NULL, n_draws = 1000L, seed = NULL,
                  progress = getOption("indri.progress", interactive())) {
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
    if (missing == "model" && isTRUE(engine$capabilities$missing)) {
      # the engine handles the gaps
    } else if (missing == "model") {
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
  bar <- .indri_progress(length(indices), isTRUE(progress))
  blocks <- if (plan$blocks_in_parallel > 1L) {
    cluster <- parallel::makeCluster(plan$blocks_in_parallel)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    tryCatch({
      if (bar$enabled) {
        # Progress needs results back as blocks finish, so send them in
        # load-balanced rounds of a few blocks per worker (a round waits for
        # its slowest block). Without a bar, one parLapply() call.
        size <- 4L * plan$blocks_in_parallel
        rounds <- split(indices, ceiling(seq_along(indices) / size))
        out <- vector("list", length(indices))
        for (r in rounds) {
          out[r] <- parallel::clusterApplyLB(
            cluster, r, .indri_fit_block, Y = Y_std, design = X,
            a_forms = a_forms, engine = engine, engine_args = engine_args,
            subsets = subsets)
          bar$tick(length(r))
        }
        out
      } else {
        parallel::parLapply(cluster, indices, .indri_fit_block, Y = Y_std,
                            design = X, a_forms = a_forms, engine = engine,
                            engine_args = engine_args, subsets = subsets)
      }
    },
      error = function(e) {
        bar$close()
        if (grepl("licen", conditionMessage(e), ignore.case = TRUE))
          stop("a block worker could not check out an ", engine$name,
               " licence (", conditionMessage(e), "). Each worker needs ",
               "its own seat: lower `cores` to the number of free seats ",
               "(`cores = 1` fits the blocks in this session).",
               call. = FALSE)
        stop(e)
      })
  } else {
    lapply(indices, function(i) {
      out <- .indri_fit_block(i, Y = Y_std, design = X, a_forms = a_forms,
                              engine = engine, engine_args = engine_args,
                              subsets = subsets)
      bar$tick(1L)
      out
    })
  }
  bar$close()
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

  # The residual covariance, for indri_gwas(). Not assembled from the block
  # E: their average can be strongly indefinite at many traits, and its
  # projection then has eigenvalues near zero that miscalibrate the
  # association tests (30 traits: 66x the nominal rate at 1e-4). E is the
  # residual covariance of Y after X minus the assembled G instead; the
  # block average is kept as E_bar.
  E_parts <- list(
    E = .indri_residual_E(Y_std, X, G, A),
    E_bar = indri_assemble(lapply(blocks, function(b) list(G = b$E)),
                           subsets, K)$G)

  # REML blocks carry a sampling covariance instead of draws: draw them here,
  # in the main process, so the result does not depend on how many workers
  # ran. Each block then holds draws like a stanimal block.
  if (!isTRUE(engine$capabilities$draws) && isTRUE(engine$capabilities$vcov)) {
    n_draws <- .indri_count(n_draws, "n_draws", min = 2L)
    restore <- .indri_use_seed(if (is.null(seed)) NULL else seed + 2L)
    for (i in seq_along(blocks))
      blocks[[i]]$draws <- .indri_normal_draws(blocks[[i]]$G, blocks[[i]]$vcov,
                                               n_draws)
    restore()
  }

  G_draws <- NULL
  corrG_draws <- NULL
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

  # The engine may or may not carry trait names through; naming is indri's job.
  trait_names <- colnames(Y)
  blocks <- lapply(seq_along(blocks), function(i)
    .indri_name_block(blocks[[i]], trait_names[subsets[i, ]]))
  diagnostics <- .indri_diagnostics(blocks, subsets, trait_names)
  # A stanimal block with Rhat >= 1.05 did not mix: that is a sampling problem
  # worth a warning. REML blocks that report no convergence are listed in
  # `diagnostics` only (see .indri_diagnostics()).
  message <- if (identical(engine$name, "stanimal"))
    .indri_convergence_warning(diagnostics, nrow(subsets),
                               engine_args$lkj_prior %||% 1)
  if (!is.null(message) && isTRUE(getOption("indri.warn_rhat", TRUE)))
    warning(message, call. = FALSE)
  # Everything above runs on the standardised scale: block fits, weights,
  # averaging, the PSD projection, the REML draws and the correlation draws.
  # Correlations come from the standardised G (they are dimensionless); the
  # covariance outputs are then returned in the units of Y. Not projected
  # again: the unweighted projection objective depends on the units, so a
  # projection on the original scale would be a different estimate.
  corrG <- stats::cov2cor(G)
  dimnames(corrG) <- list(trait_names, trait_names)
  G <- .indri_to_original(G, y_sd)
  G_bar <- .indri_to_original(assembled$G, y_sd)
  G_draws <- .indri_to_original(G_draws, y_sd)
  E <- .indri_to_original(E_parts$E, y_sd)
  E_bar <- .indri_to_original(E_parts$E_bar, y_sd)
  dimnames(G) <- dimnames(G_bar) <- list(trait_names, trait_names)
  dimnames(E) <- dimnames(E_bar) <- list(trait_names, trait_names)
  if (!is.null(G_draws)) {
    dimnames(G_draws) <- list(trait_names, trait_names, NULL)
    dimnames(corrG_draws) <- list(trait_names, trait_names, NULL)
  }

  structure(list(
    G = G,
    G_bar = G_bar,
    G_draws = G_draws,
    corrG = corrG,
    corrG_draws = corrG_draws,
    sigma_G = sqrt(pmax(diag(G), 0)),
    E = E,
    E_bar = E_bar,
    counts = .indri_named(assembled$counts, trait_names),
    coverage = .indri_named(indri_pair_coverage(subsets, K), trait_names),
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
    diagnostics = diagnostics,
    call = call
  ), class = "indri_fit")
}

# The residual covariance E from the phenotypes and the assembled G. The
# residuals R of Y on X have expected covariance E + c G, with
# c = tr(M A) / (n - J) and M the residual-maker of X (c = 1 when A = I, and
# close to 1 for most designs), so E = R'R / (n - J) - c G. Traits with gaps
# are regressed on their observed rows, and each covariance uses the rows
# both traits have. The eigenvalues are floored at 1e-3 of the largest, so E
# is positive definite and well conditioned for the whitening in
# indri_gwas().
.indri_residual_E <- function(Y, X, G, A) {
  N <- nrow(Y)
  J <- ncol(X)
  R <- matrix(NA_real_, N, ncol(Y))
  for (k in seq_len(ncol(Y))) {
    seen <- !is.na(Y[, k])
    R[seen, k] <- stats::lm.fit(X[seen, , drop = FALSE], Y[seen, k])$residuals
  }
  seen <- !is.na(R)
  R[!seen] <- 0
  pairs <- crossprod(seen * 1)
  P <- crossprod(R) / pmax(pairs - J, 1)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  c_A <- (sum(diag(A)) - sum(diag(XtX_inv %*% crossprod(X, A %*% X)))) /
    (N - J)
  E <- P - c_A * G
  e <- eigen((E + t(E)) / 2, symmetric = TRUE)
  floor <- 1e-3 * max(e$values[1L], .Machine$double.eps)
  E <- e$vectors %*% (pmax(e$values, floor) * t(e$vectors))
  (E + t(E)) / 2
}

.indri_diagnostics <- function(blocks, subsets = NULL, trait_names = NULL) {
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
    if (is.logical(v)) return(if (is.na(v)) NA else v)
    if (is.numeric(v)) {
      if (!is.finite(v)) return(NA)
      return(v < 1.05)
    }
    NA
  }, logical(1))

  # Which blocks said no. For the REML engines this is almost always "optimum
  # reached, convergence not certified": a convergence study (STATUS.md) found the
  # flagged pairs as accurate as the rest, at the same optimum from every
  # starting point, and no refit (more iterations, other starts, a larger
  # block) improved them. They are reported, not refitted or dropped.
  # Divergent transitions (stanimal) also mark a block as not converged.
  divergences <- vapply(blocks, function(b) {
    v <- b$divergences
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) as.integer(v)
    else NA_integer_
  }, integer(1))
  block_ok[which(divergences > 0L)] <- FALSE

  bad <- which(block_ok %in% FALSE)
  traits <- if (is.null(subsets)) rep(NA_character_, length(bad)) else
    vapply(bad, function(i) {
      idx <- subsets[i, ]
      paste(if (is.null(trait_names)) idx else trait_names[idx],
            collapse = ", ")
    }, character(1))

  list(
    rhat = rhat,
    max_rhat = suppressWarnings(max(rhat, na.rm = TRUE)),
    # NA when no block could tell (e.g. stanimal fits through rstan)
    converged = if (all(is.na(block_ok))) NA else all(block_ok, na.rm = TRUE),
    n_unconverged = length(bad),
    n_unknown = sum(is.na(block_ok)),
    n_divergent_blocks = sum(divergences > 0L, na.rm = TRUE),
    divergences = divergences,
    unconverged = data.frame(block = bad, traits = traits,
                             rhat = rhat[bad], divergences = divergences[bad],
                             stringsAsFactors = FALSE)
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
# the sqrt(K) rule; the REML engines use s = 3 (4-trait blocks), which the
# sweeps at K = 16 and 30 (STATUS.md) found slightly more
# accurate than pairs at the same or lower cost. stanimal asks for each pair
# in 3 distinct trait contexts, the REML engines for 1: at s = 3 three
# contexts gave the same accuracy and coverage as one (examples/
# STATUS.md, K = 12 and 30) at about twice the blocks.
# indri_subsets() never repeats a block.
.indri_engine_defaults <- function(engine, K) {
  caps <- engine$capabilities
  s <- caps$default_s %||% .indri_default_s(K)
  # a subset covers up to 2 (s - 1) traits; keep it narrower than K, or there
  # is nothing to reassemble (pairs at K = 3 or 4)
  s <- min(as.integer(s), max(2L, (K + 1L) %/% 2L))
  list(s = max(2L, min(s, K)),
       min_cooccurrence = caps$default_cooccurrence %||% 3L)
}

# The warning for stanimal blocks that did not converge, or NULL: which
# diagnostics flagged them (rank-normalised split-Rhat >= 1.05 on the block's
# G draws, from posterior::rhat(), and divergent transitions), and a
# suggestion. On real data whose block posteriors were poorly identified
# (ratones, 35 traits, 329 mice) three times longer chains changed nothing,
# while a stronger LKJ prior is what regularises the block correlations.
# Silence it with options(indri.warn_rhat = FALSE) (the test suite does: its
# fits are short on purpose).
.indri_convergence_warning <- function(diagnostics, n_blocks, lkj_prior = 1) {
  slow <- sum(diagnostics$unconverged$rhat >= 1.05, na.rm = TRUE)
  divergent <- diagnostics$n_divergent_blocks %||% 0L
  if (slow == 0L && divergent == 0L) return(NULL)
  found <- c(
    if (slow > 0L)
      sprintf("%d of %d have rank-normalised split-Rhat >= 1.05 on G (max %s)",
              slow, n_blocks, signif(diagnostics$max_rhat, 3)),
    if (divergent > 0L)
      sprintf("%d had divergent transitions (%d in total)", divergent,
              sum(diagnostics$divergences, na.rm = TRUE)))
  n_bad <- length(unique(diagnostics$unconverged$block))
  suggest <- if (lkj_prior < 2) 2 else lkj_prior + 1
  paste0(if (n_bad == 1L) "A block" else "Several blocks",
         " failed to converge: ", paste(found, collapse = "; "), ". ",
         "The likelihood may not identify the genetic correlations strongly ",
         "enough with these data. Consider increasing the LKJ prior ",
         "parameter, e.g. ",
         "`engine_args = list(lkj_prior = ", suggest, ")` (this fit used ",
         lkj_prior, "), which regularises the block correlations; longer ",
         "chains do not help when block posteriors are poorly identified. See ",
         "`fit$diagnostics$unconverged`.")
}

# A console progress bar for the block fits, on stderr: blocks done, elapsed
# time and a linear estimate of the time left. Disabled, it does nothing.
.indri_progress <- function(total, enabled) {
  noop <- function(...) invisible(NULL)
  if (!enabled || total < 1L)
    return(list(enabled = FALSE, tick = noop, close = noop))
  started <- proc.time()[["elapsed"]]
  done <- 0L
  open <- TRUE
  clock <- function(sec) {
    sec <- round(sec)
    if (sec < 60) sprintf("%ds", sec)
    else if (sec < 3600) sprintf("%dm %02ds", sec %/% 60, sec %% 60)
    else sprintf("%dh %02dm", sec %/% 3600, (sec %% 3600) %/% 60)
  }
  draw <- function() {
    elapsed <- proc.time()[["elapsed"]] - started
    width <- 30L
    filled <- round(width * done / total)
    left <- if (done > 0L && done < total)
      paste0(", ~", clock(elapsed / done * (total - done)), " left") else ""
    cat(sprintf("\r  fitting blocks [%s%s] %d/%d, %s%s   ",
                strrep("=", filled), strrep(" ", width - filled), done, total,
                clock(elapsed), left), file = stderr())
    utils::flush.console()
  }
  draw()
  list(
    enabled = TRUE,
    tick = function(n = 1L) {
      done <<- min(total, done + as.integer(n))
      draw()
    },
    close = function() {
      if (open) cat("\n", file = stderr())
      open <<- FALSE
      invisible(NULL)
    }
  )
}

# Standardised-scale covariance (a K x K matrix, or K x K x D draws) to the
# units of Y: entry (i, j) times y_sd[i] * y_sd[j]. NULL stays NULL and NA
# (uncovered) stays NA.
.indri_to_original <- function(x, y_sd) {
  if (is.null(x)) return(NULL)
  scale <- outer(y_sd, y_sd)
  if (length(dim(x)) == 3L) {
    out <- x
    for (d in seq_len(dim(x)[3L])) out[, , d] <- x[, , d] * scale
    return(out)
  }
  x * scale
}

# A K x K matrix with the trait names on both dimensions.
.indri_named <- function(m, names) {
  dimnames(m) <- list(names, names)
  m
}

# Name a block's matrices after its own traits, whatever the engine returned:
# G, E and se are traits x traits, draws traits x traits x draws, and vcov is
# indexed by the lower-triangle entries ("b,a" is the covariance of b with a).
.indri_name_block <- function(block, traits) {
  k <- length(traits)
  for (part in c("G", "E", "se")) {
    m <- block[[part]]
    if (is.matrix(m) && identical(dim(m), c(k, k)))
      block[[part]] <- .indri_named(m, traits)
  }
  if (!is.null(block$draws) && identical(dim(block$draws)[1:2], c(k, k)))
    dimnames(block$draws) <- list(traits, traits, NULL)
  if (is.matrix(block$vcov)) {
    index <- .indri_vech_index(k)
    if (nrow(block$vcov) == nrow(index)) {
      entry <- paste(traits[index[, 1L]], traits[index[, 2L]], sep = ",")
      dimnames(block$vcov) <- list(entry, entry)
    }
  }
  block
}

# Correlation draws of one block, draw by draw. A normal-approximation draw
# near a boundary can have a non-positive variance or an indefinite matrix;
# its variances are floored and its correlations clamped to [-1, 1], so every
# entry is a correlation (posterior draws are unaffected).
.indri_corr_draws <- function(draws) {
  out <- draws
  for (d in seq_len(dim(draws)[3L])) {
    m <- draws[, , d, drop = FALSE][, , 1L]
    s <- sqrt(pmax(diag(as.matrix(m)), 1e-12))
    r <- pmin(pmax(as.matrix(m) / outer(s, s), -1), 1)
    diag(r) <- 1
    out[, , d] <- r
  }
  out
}

# Default subset size: grow with K like sqrt(K), but never past 3 (4-trait
# blocks). At K = 35 the stanimal engine covered nominally with s = 3 and
# undercovered with s = 4 (0.96 vs 0.94 at LKJ 1, 0.94 vs 0.88 at LKJ 2), and
# was more accurate; the REML engines measured the same or better at s = 3.
.indri_default_s <- function(K) {
  max(2L, min(K - 1L, floor(sqrt(K)), 3L))
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
  if (ncol(X) == 0L) .indri_check_design(X)
  if (is.null(colnames(X)))
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (is.null(rownames(X))) rownames(X) <- ids
  storage.mode(X) <- "double"
  .indri_check_design(X)
  X
}

# The design every engine can use, as stanimal requires of its own: finite,
# full column rank, and able to represent a constant (an intercept column, or
# the 0/1 indicators of every level of a factor). The responses are centred,
# and without a constant the REML engines would fit different models (sommer
# has no intercept at all, asreml an implicit one).
.indri_check_design <- function(X) {
  hint <- paste0("Build `Y`, `X` and `A` from a formula with ",
                 "`stanimal::animal_design()`, or add an intercept column.")
  if (ncol(X) == 0L)
    stop("`X` has no columns. ", hint, call. = FALSE)
  if (any(!is.finite(X)))
    stop("`X` has missing or non-finite values; drop those individuals ",
         "first (`stanimal::animal_design()` does).", call. = FALSE)
  if (qr(X)$rank < ncol(X))
    stop("`X` does not have full column rank: some columns are linear ",
         "combinations of others (for example an intercept together with the ",
         "indicators of every level of a factor). ", hint, call. = FALSE)
  fitted <- stats::lm.fit(X, rep(1, nrow(X)))$fitted.values
  if (max(abs(fitted - 1)) > 1e-8)
    stop("`X` cannot represent a constant, so the trait means are not in the ",
         "model (the responses are centred). Include an intercept, or the ",
         "indicators of every level of a factor. ", hint, call. = FALSE)
  invisible(TRUE)
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
