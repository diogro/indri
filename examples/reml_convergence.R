# What are the non-converged REML pairs, and can a refit rescue them?
#
# Fits every trait pair once (s = 2, c = 1, the REML default), then for every
# pair whose fit reports no convergence:
#   1. describes it: estimated genetic correlation, heritabilities, error
#      against the truth, next to the same numbers for the converged pairs;
#   2. refits it three ways and records convergence and error:
#        iter    more iterations (sommer maxIters 300; asreml maxit 50 with
#                up to 5 update passes)
#        restart starting values at the non-converged estimate itself (an AI
#                step that stalled often moves on from a fresh start)
#        half    starting values G = E = P / 2, P the pair's phenotypic
#                covariance: h2 = 0.5 and the phenotypic correlation split
#                evenly
#        uni     starting variances from the two univariate fits, with the
#                phenotypic correlation split evenly between G and E
#        triple  the pair inside a 3-trait block, with the trait most
#                correlated with the pair (phenotypically) added
#
#   INDRI_CONV_ENGINE   sommer or asreml (sommer)
#   INDRI_CONV_K        traits (16)
#   INDRI_CONV_SEED     simulation seed (1); same simulation as
#                       examples/block_size_sweep.R
#   INDRI_CONV_CORES    workers for the initial fit (8; asreml: 2)
#
# Writes results/reml_convergence_<engine>_K<K>_seed<seed>.csv (git-ignored).

library(indri)

engine_name <- Sys.getenv("INDRI_CONV_ENGINE", "sommer")
K <- as.integer(Sys.getenv("INDRI_CONV_K", "16"))
seed <- as.integer(Sys.getenv("INDRI_CONV_SEED", "1"))
cores <- as.integer(Sys.getenv("INDRI_CONV_CORES",
                               if (engine_name == "asreml") "2" else "8"))
Sys.setenv(OMP_NUM_THREADS = "1")
dir.create("results", showWarnings = FALSE)

## Simulate, as block_size_sweep.R does -------------------------------------

n_fam <- 60L; size <- 5L
set.seed(seed)
truth <- stanimal::random_covariance(
  n_traits = K, method = "partial", min_var = 0.8, max_var = 1.4,
  trait_names = paste0("trait", seq_len(K)))
resid <- stanimal::random_covariance(
  n_traits = K, method = "lkj", lkj_prior = 2, min_var = 0.5, max_var = 0.9)
A <- kronecker(diag(n_fam), matrix(0.5, size, size)); diag(A) <- 1
N <- nrow(A)
Y <- t(chol(A)) %*% matrix(stats::rnorm(N * K), N, K) %*% chol(truth) +
  matrix(stats::rnorm(N * K), N, K) %*% chol(resid)
rownames(Y) <- paste0("id", seq_len(N)); colnames(Y) <- colnames(truth)
corr_true <- stats::cov2cor(truth)
h2_true <- diag(truth) / (diag(truth) + diag(resid))

## Fit every pair once --------------------------------------------------------

fit <- fit_G(Y, A = A, engine = engine_name, cores = cores, seed = 1000L + seed)
stopifnot(ncol(fit$subsets) == 2L)

engine <- indri_engine(engine_name)
Y_std <- sweep(sweep(Y, 2L, colMeans(Y), "-"), 2L, apply(Y, 2L, stats::sd), "/")
X <- matrix(1, N, 1L, dimnames = list(rownames(Y), "intercept"))
a_forms <- indri_a_forms(A, need = engine$capabilities$a_form)
P_corr <- stats::cor(Y_std)

describe <- function(block, idx) {
  G <- block$G; E <- block$E
  i <- idx[1]; j <- idx[2]
  a <- match(i, idx); b <- match(j, idx)
  r <- G[a, b] / sqrt(G[a, a] * G[b, b])
  data.frame(est_corr = r, abs_err = abs(r - corr_true[i, j]),
             h2_i = G[a, a] / (G[a, a] + E[a, a]),
             h2_j = G[b, b] / (G[b, b] + E[b, b]),
             converged = isTRUE(block$converged))
}

refit <- function(idx, control = list()) {
  out <- tryCatch(engine$fit(Y_std[, idx, drop = FALSE], X, a_forms, control),
                  error = function(e) NULL)
  out
}

# A start must be positive definite: pull a near-singular estimate inwards.
shrink <- function(m) {
  r <- stats::cov2cor(m)[1, 2]
  if (!is.finite(r) || abs(r) > 0.95) {
    r <- sign(r) * 0.95
    d <- sqrt(pmax(diag(m), 0.05))
    m <- diag(d) %*% matrix(c(1, r, r, 1), 2) %*% diag(d)
  }
  m
}
restart_values <- function(block) list(G = shrink(block$G), E = shrink(block$E))
univariate <- lapply(seq_len(K), function(k) refit(k))
uni_values <- function(idx) {
  u <- univariate[idx]
  if (any(vapply(u, is.null, logical(1)))) return(NULL)
  g <- pmax(vapply(u, function(x) x$G[1, 1], numeric(1)), 0.05)
  e <- pmax(vapply(u, function(x) x$E[1, 1], numeric(1)), 0.05)
  r <- max(min(P_corr[idx[1], idx[2]], 0.9), -0.9)
  cov2 <- function(v) matrix(c(v[1], r * sqrt(v[1] * v[2]),
                               r * sqrt(v[1] * v[2]), v[2]), 2)
  list(G = cov2(g), E = cov2(e))
}
half_values <- function(idx) {
  P <- stats::cov(Y_std[, idx])
  list(G = P / 2, E = P / 2)
}

more_iterations <- if (engine_name == "sommer") {
  list(convergence = list(maxIters = 300L))
} else list(maxit = 50L, updates = 5L)

rows <- list()
for (b in seq_len(nrow(fit$subsets))) {
  idx <- fit$subsets[b, ]
  base <- describe(fit$blocks[[b]], idx)
  rows[[length(rows) + 1L]] <- cbind(pair = paste(idx, collapse = "-"),
    truth = corr_true[idx[1], idx[2]], h2_true_i = h2_true[idx[1]],
    h2_true_j = h2_true[idx[2]], strategy = "base", base)
  if (base$converged) next

  message("pair ", paste(idx, collapse = "-"), " did not converge; refitting")
  attempts <- list(iter = list(idx = idx, control = more_iterations),
                   restart = list(idx = idx, control = list(
                     init = restart_values(fit$blocks[[b]]))),
                   half = list(idx = idx, control = list(init = half_values(idx))),
                   uni = list(idx = idx, control = list(init = uni_values(idx))))
  third <- setdiff(seq_len(K), idx)
  third <- third[which.max(pmax(abs(P_corr[idx[1], third]), abs(P_corr[idx[2], third])))]
  attempts$triple <- list(idx = sort(c(idx, third)), control = list())
  for (name in names(attempts)) {
    a <- attempts[[name]]
    if (name == "uni" && is.null(a$control$init)) next
    out <- refit(a$idx, a$control)
    if (is.null(out)) next
    keep <- match(idx, a$idx)
    sub <- list(G = out$G[keep, keep], E = out$E[keep, keep],
                converged = out$converged)
    rows[[length(rows) + 1L]] <- cbind(pair = paste(idx, collapse = "-"),
      truth = corr_true[idx[1], idx[2]], h2_true_i = h2_true[idx[1]],
      h2_true_j = h2_true[idx[2]], strategy = name, describe(sub, idx))
  }
}

out <- do.call(rbind, rows)
out$engine <- engine_name; out$K <- K; out$seed <- seed
file <- file.path("results", sprintf("reml_convergence_%s_K%d_seed%d.csv",
                                     engine_name, K, seed))
utils::write.csv(out, file, row.names = FALSE)
base <- out[out$strategy == "base", ]
cat(sprintf("\n%s K = %d seed %d: %d of %d pairs did not converge\n",
            engine_name, K, seed, sum(!base$converged), nrow(base)))
failed <- out[out$pair %in% base$pair[!base$converged], ]
if (nrow(failed))
  print(aggregate(cbind(converged, abs_err) ~ strategy, failed, mean))
cat("converged pairs, median abs err:", round(median(base$abs_err[base$converged]), 3), "\n")
cat("written: ", file, "\n", sep = "")
