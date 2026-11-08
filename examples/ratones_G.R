# indri on the ratones mice: estimate the full additive genetic covariance
# matrix G and compare it with the phenotypic covariance matrix P.
#
# The data ships with stanimal (ratones_morpho and ratones_pedigree): 329 mice,
# 35 skull traits, five lines. AGE, SEX and LIN are fitted as fixed effects; the
# relationship matrix is built from the pedigree with stanimal::makeA().
#
# Run from the repository root:
#   Rscript examples/ratones_G.R
#
# Needs stanimal and CmdStan. Environment overrides:
#   INDRI_RATONES_CORES   Stan process budget (64)
#   INDRI_RATONES_ITER    iterations per chain (1000)
#   INDRI_RATONES_CHAINS  chains per block (4)
#   INDRI_RATONES_SEED    seed (1)

library(indri)

if (!requireNamespace("stanimal", quietly = TRUE))
  stop("This example needs stanimal")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  stop("This example needs cmdstanr and a CmdStan installation")

cores <- as.integer(Sys.getenv("INDRI_RATONES_CORES", "64"))
iter <- as.integer(Sys.getenv("INDRI_RATONES_ITER", "1000"))
chains <- as.integer(Sys.getenv("INDRI_RATONES_CHAINS", "4"))
seed <- as.integer(Sys.getenv("INDRI_RATONES_SEED", "1"))

## Data ----------------------------------------------------------------------

data(ratones_morpho, package = "stanimal")
data(ratones_pedigree, package = "stanimal")

meta <- c(".id", "ID", "SEX", "GER", "LIN", "MADRE", "PADRE", "AGE", "P49",
          "line", "original_line", "selection")
traits <- setdiff(names(ratones_morpho), meta)
Y <- as.matrix(log(ratones_morpho[, traits]))

# Fixed effects: age, sex and line.
X <- stats::model.matrix(~ AGE + SEX + LIN, data = ratones_morpho)

# Additive relationship matrix from the pedigree, aligned to the phenotyped
# mice. stanimal::makeA() sorts the pedigree itself.
A <- stanimal::makeA(ratones_pedigree)[rownames(ratones_morpho),
                                       rownames(ratones_morpho)]

message("N = ", nrow(Y), " individuals, ", ncol(Y), " traits")

## Fit ------------------------------------------------------------------------

# The default subset size caps at 4, so a 35-trait fit runs as ~400 six-trait
# blocks. The stanimal engine returns posterior draws per block.
fit <- fit_G(
  Y, X = X, A = A, engine = "stanimal",
  engine_args = list(chains = chains, iter = iter, init = 0.5,
                     refresh = 0, show_messages = FALSE),
  cores = cores, seed = seed
)

corrG <- fit$corrG                     # additive genetic correlations
# Phenotypic correlation of the same size-corrected traits, so that G and P
# are comparable.
resid <- stats::residuals(stats::lm(Y ~ 0 + X))
corrP <- stats::cor(resid)
pairs <- upper.tri(corrP)

cat("\n== G (indri) vs P (phenotypic) ==\n")
cat("correlation of correlations: ", round(stats::cor(corrG[pairs], corrP[pairs]), 3), "\n", sep = "")
cat("median |corrG - corrP|:       ", round(stats::median(abs(corrG[pairs] - corrP[pairs])), 3), "\n", sep = "")
cat("relative Frobenius distance:  ", round(norm(corrG - corrP, "F") / norm(corrP, "F"), 3), "\n", sep = "")
cat("min eigenvalue of G:          ", round(min(eigen(fit$G, symmetric = TRUE, only.values = TRUE)$values), 4), "\n", sep = "")

dir.create("results", showWarnings = FALSE)
utils::write.csv(corrG, "results/ratones_corrG.csv")
utils::write.csv(corrP, "results/ratones_corrP.csv")
cat("\nwritten: results/ratones_corrG.csv, results/ratones_corrP.csv\n")

## Plot -----------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Plotting needs ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE))
  stop("Plotting needs patchwork")

pG <- stanimal::plot_corr(corrG, "corrG  (genetic)")
pP <- stanimal::plot_corr(corrP, "corrP  (phenotypic)")
both <- patchwork::wrap_plots(pG, pP, ncol = 2)

ggplot2::ggsave("results/ratones_corrG_vs_P.png", both,
                width = 13, height = 7, dpi = 150)
cat("written: results/ratones_corrG_vs_P.png\n")
