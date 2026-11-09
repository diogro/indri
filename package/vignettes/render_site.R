# Render the precomputed vignettes to self-contained HTML in docs/ (the
# repository's GitHub Pages site), with an index page. Run from
# package/vignettes after precompute.R:
#   Rscript render_site.R
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
items <- paste0('  <li><a href="', names(vignettes), '.html">', vignettes,
                "</a></li>", collapse = "\n")
writeLines(c(
  "<!DOCTYPE html>", '<html lang="en">', "<head>", '<meta charset="utf-8">',
  "<title>indri vignettes</title>",
  "<style>body { font-family: sans-serif; max-width: 42em; margin: 3em auto; line-height: 1.5; padding: 0 1em; }</style>",
  "</head>", "<body>", "<h1>indri vignettes</h1>",
  "<p>Worked analyses of the nine weekly weight gains of the LG/SM F3 mice that ship with <a href=\"https://github.com/diogro/indri\">indri</a>.</p>",
  "<ul>", items, "</ul>", "</body>", "</html>"), file.path(site, "index.html"))
message("rendered ", length(vignettes), " vignettes into ", site)
