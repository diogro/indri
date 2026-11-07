# Estimating a 6-trait additive genetic covariance matrix from trait subsets.
#
# Run from the repository root:
#   Rscript examples/six_traits.R
#
# Needs stanimal and CmdStan. The 6 x 6 truth comes from
# stanimal::random_covariance(), which draws a random correlation matrix and
# scales it by per-trait variances.

library(indri)

if (!requireNamespace("stanimal", quietly = TRUE))
  stop("This example needs stanimal: install.packages('stanimal')")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("This example needs cmdstanr and a CmdStan installation.")

set.seed(2024)
K <- 6L
trait_names <- paste0("trait", seq_len(K))

## The truth ---------------------------------------------------------------

# Genetic covariance. "partial" draws partial correlations uniformly and turns
# them into a correlation matrix, so the truth has both signs and a spread of
# magnitudes, which makes the recovery visible.
G_true <- stanimal::random_covariance(
  n_traits = K, method = "partial",
  min_var = 0.8, max_var = 1.4, trait_names = trait_names
)

# Residual covariance, independent of the genetic one.
E_true <- stanimal::random_covariance(
  n_traits = K, method = "lkj", lkj_prior = 2,
  min_var = 0.5, max_var = 0.9, trait_names = trait_names
)

## Simulate a cohort -------------------------------------------------------

# 25 families of 8 full sibs: a relationship matrix with real structure that
# does not need a pedigree.
n_fam <- 25L
sibship_size <- 8L
A <- kronecker(diag(n_fam), matrix(0.5, sibship_size, sibship_size))
diag(A) <- 1
N <- nrow(A)

# Y = genetic + residual, with breeding values integrated over the relatedness
# structure: rows of the genetic term are correlated through A.
genetic <- t(chol(A)) %*% matrix(stats::rnorm(N * K), N, K) %*% chol(G_true)
Y <- genetic + matrix(stats::rnorm(N * K), N, K) %*% chol(E_true)
rownames(Y) <- paste0("id", seq_len(N))
colnames(Y) <- trait_names

## Fit ---------------------------------------------------------------------

# A note on what to expect. The block fits are Bayesian, and stanimal puts an
# LKJ(2) prior on the trait correlations, which pulls them towards zero; at 200
# individuals that prior still has weight, so the reported correlations sit a
# little closer to zero than the truth. This is the estimator, not the
# subsetting: stanimal's own full 6-trait model shrinks the same way. Widen the
# prior with `engine_args = list(lkj_prior = 1)` if you want it weaker.

# `s = 3` makes each subset cover 3 - 1 = 2 trait pairs, which works out to 4
# traits per subset. With min_cooccurrence = 3 every pair of the 6 traits is
# estimated at least three times.
started <- Sys.time()
fit <- fit_G(
  Y, A = A, s = 3L, min_cooccurrence = 3L, engine = "stanimal",
  engine_args = list(chains = 2L, iter = 1000L, cores = 2L, seed = 1L,
                     refresh = 0, show_messages = FALSE),
  seed = 2024L
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

## Compare with the truth --------------------------------------------------

corr_true <- stats::cov2cor(G_true)
corr_hat <- fit$corrG
pairs <- upper.tri(corr_true)

relative_error <- norm(corr_hat - corr_true, "F") / norm(corr_true, "F")

cat("\n== indri, 6 traits ==\n")
cat("individuals:            ", N, "\n")
cat("subsets:                ", nrow(fit$subsets), " of ", ncol(fit$subsets),
    " traits\n", sep = "")
cat("pair coverage:          ", min(fit$coverage[pairs]), " to ",
    max(fit$coverage[pairs]), "\n", sep = "")
cat("fit time:               ", round(elapsed, 1), " min\n", sep = "")
cat("min eigenvalue of G:    ", signif(min(eigen(fit$G, symmetric = TRUE,
                                                only.values = TRUE)$values), 4),
    "\n", sep = "")
cat("relative error, corrG:  ", round(relative_error, 3), "\n", sep = "")

cat("\ngenetic correlations, truth vs estimate:\n")
table_out <- data.frame(
  pair = paste(trait_names[row(corr_true)[pairs]], trait_names[col(corr_true)[pairs]],
               sep = "-"),
  true = round(corr_true[pairs], 2),
  indri = round(corr_hat[pairs], 2),
  difference = round((corr_hat - corr_true)[pairs], 2)
)
print(table_out, row.names = FALSE)

# Variances: fit$G is in the units of the traits, so it compares directly.
G_original <- fit$G
cat("\ngenetic variances, truth vs estimate:\n")
print(data.frame(trait = trait_names,
                 true = round(diag(G_true), 2),
                 indri = round(diag(G_original), 2)), row.names = FALSE)

## Figure ------------------------------------------------------------------

dir.create("examples/figures", showWarnings = FALSE, recursive = TRUE)
grDevices::png("examples/figures/six_traits_corrG.png", width = 1200,
               height = 420, res = 110)
graphics::par(mfrow = c(1, 3), mar = c(6, 6, 3, 1))

draw_corr <- function(m, main) {
  k <- nrow(m)
  graphics::image(seq_len(k), seq_len(k), m[, k:1], zlim = c(-1, 1),
                  col = grDevices::hcl.colors(21, "Blue-Red 3"),
                  axes = FALSE, xlab = "", ylab = "", main = main)
  graphics::axis(1, at = seq_len(k), labels = colnames(m), las = 2, cex.axis = 0.7)
  graphics::axis(2, at = seq_len(k), labels = rev(colnames(m)), las = 2,
                 cex.axis = 0.7)
  graphics::box()
}

draw_corr(corr_true, "true corrG")
draw_corr(corr_hat, "indri corrG")
graphics::plot(corr_true[pairs], corr_hat[pairs], pch = 19,
               xlim = c(-1, 1), ylim = c(-1, 1),
               xlab = "true", ylab = "indri",
               main = sprintf("genetic correlations (relative error %.2f)",
                              relative_error))
graphics::abline(0, 1, lty = 2)
grDevices::dev.off()

cat("\nfigure: examples/figures/six_traits_corrG.png\n")
