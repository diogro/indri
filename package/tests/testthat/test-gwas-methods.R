# plot() (ggplot2 Manhattan and QQ plots) and summary() for indri_gwas.

scan_with_map <- function(dominance = FALSE, ...) {
  d <- simulate_gwas_data(M = 12L, causal = 1L, effect = 0.7, seed = 81)
  info <- data.frame(snp = colnames(d$geno), chrom = rep(1:3, each = 4),
                     pos = rep(c(10, 20, 30, 40), 3) * 1e5)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     snp_info = info, dominance = dominance, ...)
  list(d = d, scan = scan, info = info)
}

test_that("the Manhattan plot draws -log10 p along the genome", {
  skip_if_not_installed("ggplot2")
  s <- scan_with_map()
  p <- plot(s$scan)
  expect_s3_class(p, "ggplot")
  points <- p$data
  expect_equal(nrow(points), 12L)
  expect_equal(points$logp, -log10(unname(s$scan$pvalues)))
  # chromosomes laid end to end, positions increasing within each
  expect_true(all(diff(points$x) > 0))
  expect_equal(as.character(points$chrom), as.character(s$info$chrom))
  # the Bonferroni line over SNPs by default
  hline <- Filter(function(l) inherits(l$geom, "GeomHline"), p$layers)
  expect_length(hline, 1L)
  expect_equal(hline[[1]]$data$yintercept, -log10(0.05 / 12))
  built <- ggplot2::ggplot_build(p)
  expect_gt(length(built$data), 0L)
})

test_that("without positions the Manhattan plot uses the SNP order", {
  skip_if_not_installed("ggplot2")
  d <- simulate_gwas_data(M = 7L, seed = 82)
  scan <- indri_gwas(d$Y, d$geno, X = d$X, model = "lm")
  p <- plot(scan, threshold = 1e-3)
  expect_equal(p$data$x, seq_len(7))
  hline <- Filter(function(l) inherits(l$geom, "GeomHline"), p$layers)
  expect_equal(hline[[1]]$data$yintercept, 3)
})

test_that("SNPs can be highlighted by name", {
  skip_if_not_installed("ggplot2")
  s <- scan_with_map()
  p <- plot(s$scan, highlight = colnames(s$d$geno)[c(1, 5)])
  expect_equal(sum(p$data$highlight), 2L)
  expect_error(plot(s$scan, highlight = "nope"), "highlight")
})

test_that("the QQ plot compares sorted p-values with their expectation", {
  skip_if_not_installed("ggplot2")
  s <- scan_with_map()
  p <- plot(s$scan, type = "qq")
  expect_s3_class(p, "ggplot")
  expected <- -log10(stats::ppoints(12))
  expect_equal(p$data$expected, expected)
  expect_equal(p$data$observed, -log10(sort(unname(s$scan$pvalues))))
  expect_gt(length(ggplot2::ggplot_build(p)$data), 0L)
})

test_that("single-trait scans plot one trait or the best per SNP", {
  skip_if_not_installed("ggplot2")
  d <- simulate_gwas_data(M = 6L, seed = 83)
  scan <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                     adjust = "bonferroni")
  best <- plot(scan)
  expect_equal(best$data$logp,
               -log10(unname(apply(scan$pvalues_adjust, 1L, min))))
  tail <- plot(scan, trait = "tail")
  expect_equal(tail$data$logp, -log10(unname(scan$pvalues[, "tail"])))
  expect_error(plot(scan, trait = "nope"), "trait")
  qq <- plot(scan, type = "qq")
  expect_equal(nrow(qq$data), length(scan$pvalues))
})

test_that("dominance and empirical p-values can be plotted", {
  skip_if_not_installed("ggplot2")
  s <- scan_with_map(dominance = TRUE, n_permutations = 9L, seed = 1L)
  dom <- plot(s$scan, which = "dominance")
  ok <- !is.na(s$scan$pvalues_dominance)
  expect_equal(dom$data$logp, -log10(unname(s$scan$pvalues_dominance[ok])))
  emp <- plot(s$scan, which = "empirical")
  expect_equal(emp$data$logp, -log10(unname(s$scan$pvalues_empirical)))
  plain <- scan_with_map()
  expect_error(plot(plain$scan, which = "dominance"), "dominance")
  expect_error(plot(plain$scan, which = "empirical"), "n_permutations")
})

test_that("summary reports the top SNPs and genomic inflation", {
  s <- scan_with_map()
  sm <- summary(s$scan, n = 3)
  expect_s3_class(sm, "summary.indri_gwas")
  expect_equal(nrow(sm$top), 3L)
  expect_equal(sm$top$pvalue, unname(sort(s$scan$pvalues)[1:3]))
  expect_equal(sm$top$snp[1], colnames(s$d$geno)[1])
  P <- s$d$P
  expected_lambda <- stats::qchisq(stats::median(s$scan$pvalues), P,
                                   lower.tail = FALSE) /
    stats::qchisq(0.5, P, lower.tail = FALSE)
  expect_equal(sm$lambda_gc, expected_lambda)
  expect_equal(sm$n_tested, 12L)
  expect_output(print(sm), "inflation")

  single <- indri_gwas(s$d$Y, s$d$geno, A = s$d$A, X = s$d$X,
                       mode = "single_trait")
  ss <- summary(single)
  expect_named(ss$lambda_gc, colnames(s$d$Y))
  expect_true(all(c("snp", "trait", "pvalue") %in% names(ss$top)))
  expect_output(print(ss), "single-trait")
})
