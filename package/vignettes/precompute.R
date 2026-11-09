# Precompute the vignettes: run each .Rmd.orig (fits, scans, figures) and
# write the .Rmd that R CMD build renders. Run from package/vignettes, with
# indri installed and BLAS limited to one thread:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     INDRI_VIGNETTE_CORES=64 Rscript precompute.R [names...]
# Vignettes reuse earlier ones' fits through cache/, so run them in order.
names <- commandArgs(trailingOnly = TRUE)
if (!length(names)) names <- c("indri-G", "indri-gwas", "indri-locus-G",
                               "indri-gxe")
for (name in names) {
  message("== ", name)
  started <- Sys.time()
  knitr::knit(paste0(name, ".Rmd.orig"), paste0(name, ".Rmd"),
              envir = new.env())
  message(name, ": ", round(as.numeric(difftime(Sys.time(), started,
                                                 units = "mins")), 1), " min")
}
