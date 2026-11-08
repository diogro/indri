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
# Needs stanimal (for the data) and the chosen engine. Environment overrides:
#   INDRI_RATONES_ENGINE  stanimal, sommer or asreml (stanimal)
#   INDRI_RATONES_CORES   process budget (64)
#   INDRI_RATONES_ITER    stanimal iterations per chain (1000; stanimal's own
#                         default is 3000, use it when the blocks do not mix)
#   INDRI_RATONES_CHAINS  stanimal chains per block (4)
#   INDRI_RATONES_SEED    seed (1)
#   INDRI_RATONES_TAG     suffix for the output files (none), so several runs
#                         can sit side by side

library(indri)

if (!requireNamespace("stanimal", quietly = TRUE))
  stop("This example needs stanimal")
engine <- Sys.getenv("INDRI_RATONES_ENGINE", "stanimal")
if (engine == "stanimal" && !requireNamespace("cmdstanr", quietly = TRUE))
  stop("The stanimal engine needs cmdstanr and a CmdStan installation")
tag <- Sys.getenv("INDRI_RATONES_TAG", "")
if (nzchar(tag)) tag <- paste0("_", tag)
Sys.setenv(OMP_NUM_THREADS = "1")   # sommer: let `cores` drive parallelism

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

# With the defaults a 35-trait fit runs as about 600 four-trait blocks for
# stanimal (each pair in three distinct ones) and about 200 for the REML
# engines (each pair in at least one).
engine_args <- if (engine == "stanimal") {
  list(chains = chains, iter = iter, init = 0.5, refresh = 0,
       show_messages = FALSE)
} else list()
fit <- fit_G(Y, X = X, A = A, engine = engine, engine_args = engine_args,
             cores = cores, seed = seed)

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
ev <- eigen(corrG, symmetric = TRUE, only.values = TRUE)$values
cat("eigenvalues of corrG below 1e-6: ", sum(ev < 1e-6), " of ", length(ev), "\n", sep = "")
cat("blocks not converged/certified:  ", fit$diagnostics$n_unconverged, " of ",
    nrow(fit$subsets), if (is.finite(fit$diagnostics$max_rhat))
      paste0(" (max Rhat ", signif(fit$diagnostics$max_rhat, 3), ")"), "\n", sep = "")

dir.create("results", showWarnings = FALSE)
out <- function(name) file.path("results", paste0("ratones_", name, tag, ".", "csv"))
utils::write.csv(corrG, out("corrG"))
utils::write.csv(corrP, out("corrP"))
saveRDS(list(G = fit$G, corrG = corrG, corrG_draws = fit$corrG_draws,
             diagnostics = fit$diagnostics, subsets = fit$subsets,
             engine = engine, iter = iter),
        file.path("results", paste0("ratones_fit", tag, ".rds")))
cat("\nwritten: ", out("corrG"), ", ", out("corrP"), "\n", sep = "")

## Plot -----------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Plotting needs ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE))
  stop("Plotting needs patchwork")

pG <- stanimal::plot_corr(corrG, "corrG  (genetic)")
pP <- stanimal::plot_corr(corrP, "corrP  (phenotypic)")
both <- patchwork::wrap_plots(pG, pP, ncol = 2)

png_file <- file.path("results", paste0("ratones_corrG_vs_P", tag, ".png"))
ggplot2::ggsave(png_file, both, width = 13, height = 7, dpi = 150)
cat("written: ", png_file, "\n", sep = "")
