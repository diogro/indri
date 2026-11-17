# Render the precomputed vignettes to self-contained HTML in docs/ (the
# repository's GitHub Pages site). Run from package/vignettes after
# precompute.R:
#   Rscript render_site.R
# docs/index.html is the hand-written landing page and is not touched here:
# when adding a vignette, add it to `vignettes` below and give it a card in
# docs/index.html.
site <- normalizePath(file.path("..", "..", "docs"), mustWork = FALSE)
dir.create(site, showWarnings = FALSE)
vignettes <- c(`indri-G` = "Estimating G for many traits: engines and priors",
               `indri-gwas` = "Genome scans for many traits: additive and dominance effects",
               `indri-locus-G` = "G from locus effects",
               `indri-gxe` = "Genotype-by-environment scans: loci acting differently in each sex")
for (name in names(vignettes)) {
  rmarkdown::render(paste0(name, ".Rmd"), output_dir = site, quiet = TRUE,
                    envir = new.env())
}
message("rendered ", length(vignettes), " vignettes into ", site)
