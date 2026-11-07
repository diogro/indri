#' @export
print.indri_fit <- function(x, ...) {
  K <- nrow(x$G)
  cat("<indri fit>\n")
  cat("  traits:            ", K, "\n", sep = "")
  cat("  individuals:       ", x$n_individuals, "\n", sep = "")
  cat("  engine:            ", x$engine, "\n", sep = "")
  cat("  subsets:           ", nrow(x$subsets), " of up to ",
      ncol(x$subsets), " traits\n", sep = "")
  cat("  assembly:          ", x$assembly, " (", x$weights, " weights)\n",
      sep = "")
  cat("  pair coverage:     ", min(x$coverage[upper.tri(x$coverage)]),
      " to ", max(x$coverage[upper.tri(x$coverage)]), "\n", sep = "")
  cat("  min eigenvalue:    ", signif(.indri_min_eigen(x$G), 4), "\n", sep = "")
  if (!is.null(x$diagnostics$max_rhat) && is.finite(x$diagnostics$max_rhat))
    cat("  max Rhat (G):      ", signif(x$diagnostics$max_rhat, 4), "\n", sep = "")
  n_bad <- x$diagnostics$n_unconverged
  if (!is.null(n_bad) && n_bad > 0L) {
    reml <- x$engine %in% c("sommer", "asreml")
    cat("  ", if (reml) "not certified:     " else "not converged:     ",
        n_bad, " of ", nrow(x$subsets), " blocks",
        if (reml) " (at the optimum as far as measured; see $diagnostics)"
        else " (see $diagnostics$unconverged)", "\n", sep = "")
  }
  cat("\nG (first entries):\n")
  show <- seq_len(min(4L, K))
  print(round(x$G[show, show, drop = FALSE], 3))
  invisible(x)
}

#' @export
summary.indri_fit <- function(object, ...) {
  out <- list(
    engine = object$engine,
    traits = nrow(object$G),
    subsets = nrow(object$subsets),
    subset_width = ncol(object$subsets),
    coverage_min = min(object$coverage[upper.tri(object$coverage)]),
    coverage_max = max(object$coverage[upper.tri(object$coverage)]),
    min_eigenvalue = .indri_min_eigen(object$G),
    sigma_G = object$sigma_G,
    corrG = object$corrG,
    diagnostics = object$diagnostics
  )
  class(out) <- "summary.indri_fit"
  out
}

#' @export
print.summary.indri_fit <- function(x, ...) {
  cat("<indri fit summary>\n")
  cat("  engine:        ", x$engine, "\n", sep = "")
  cat("  traits:        ", x$traits, "\n", sep = "")
  cat("  subsets:       ", x$subsets, " x ", x$subset_width, "\n", sep = "")
  cat("  coverage:      ", x$coverage_min, " to ", x$coverage_max, "\n", sep = "")
  cat("  min eigenvalue:", signif(x$min_eigenvalue, 4), "\n", sep = "")
  cat("\nsigma_G:\n")
  print(round(x$sigma_G, 3))
  invisible(x)
}

#' @export
plot.indri_fit <- function(x, ...) {
  K <- nrow(x$corrG)
  graphics::image(seq_len(K), seq_len(K), x$corrG[, K:1],
                  zlim = c(-1, 1), col = grDevices::hcl.colors(11, "Blue-Red 3"),
                  xlab = "", ylab = "", axes = FALSE, ...)
  graphics::title(main = "corrG")
  invisible(x)
}
