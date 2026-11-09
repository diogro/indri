# Empirical p-values, q-values and an FDR threshold, from scans of the same
# data with the genotype rows permuted across individuals (as LiMMBo does).
# LiMMBo's own "empirical FDR" was the fdr-quantile of the permuted
# p-values, an uncorrected threshold; indri estimates the FDR instead.
#
# The permutations are drawn in the main process before any scan:
# with `seed`, set.seed(seed) and then `n_permutations` calls to
# sample.int(N). Each permuted scan keeps the null variance ratio of the
# observed one (it does not depend on the genotypes).

permuted_pvalues <- function(d, perms, ..., which = "pvalues") {
  lapply(seq_len(ncol(perms)), function(i) {
    geno <- d$geno[perms[, i], , drop = FALSE]
    rownames(geno) <- rownames(d$geno)
    indri_gwas(d$Y, geno, ...)[[which]]
  })
}

draw_perms <- function(seed, n, N) {
  set.seed(seed)
  vapply(seq_len(n), function(i) sample.int(N), integer(N))
}

test_that("empirical p-values count permuted p-values at least as small", {
  d <- simulate_gwas_data(M = 8L, seed = 61)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     n_permutations = 30L, seed = 4L)
  expect_equal(scan$n_permutations, 30L)
  perms <- draw_perms(4L, 30L, d$N)
  null_p <- permuted_pvalues(d, perms, G = d$G, E = d$E, A = d$A, X = d$X,
                             delta = scan$delta)
  count <- Reduce(`+`, lapply(null_p, function(p) p <= scan$pvalues))
  expect_equal(scan$pvalues_empirical, (count + 1) / 31)
  expect_named(scan$pvalues_empirical, colnames(d$geno))
  expect_true(all(scan$pvalues_empirical >= 1 / 31 &
                    scan$pvalues_empirical <= 1))
})

test_that("single-trait scans get one empirical p-value per SNP and trait", {
  d <- simulate_gwas_data(M = 5L, seed = 62)
  scan <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                     n_permutations = 20L, seed = 9L)
  expect_equal(dim(scan$pvalues_empirical), c(d$M, d$P))
  perms <- draw_perms(9L, 20L, d$N)
  null_p <- permuted_pvalues(d, perms, A = d$A, X = d$X,
                             mode = "single_trait", delta = scan$delta)
  count <- Reduce(`+`, lapply(null_p, function(p) p <= scan$pvalues))
  expect_equal(scan$pvalues_empirical, (count + 1) / 21)
})

# The permutation estimate of the false discovery rate at a threshold t:
# the mean number of permuted p-values <= t per permutation over the number
# of observed p-values <= t. q-values take the smallest estimate over
# thresholds at or above each observed value.
expected_fdr <- function(observed, null_list, fdr) {
  n <- length(null_list)
  pooled <- unlist(null_list)
  at <- vapply(observed, function(t)
    min(1, sum(pooled <= t) / n / sum(observed <= t)), numeric(1))
  q <- vapply(observed, function(t) min(at[observed >= t]), numeric(1))
  passing <- observed[at <= fdr]
  list(q = q, threshold = if (length(passing)) max(passing) else NA_real_)
}

test_that("q-values and the FDR threshold come from the permuted scans", {
  d <- simulate_gwas_data(M = 10L, causal = 2L, effect = 0.5, seed = 63)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     n_permutations = 12L, fdr = 0.1, seed = 5L)
  perms <- draw_perms(5L, 12L, d$N)
  null_p <- permuted_pvalues(d, perms, G = d$G, E = d$E, A = d$A, X = d$X,
                             delta = scan$delta)
  ref <- expected_fdr(unname(scan$pvalues), null_p, 0.1)
  expect_equal(unname(scan$qvalues_empirical), ref$q)
  expect_named(scan$qvalues_empirical, colnames(d$geno))
  expect_equal(scan$fdr, 0.1)
  expect_equal(scan$fdr_threshold, ref$threshold)
})

test_that("single-trait q-values use each SNP's best adjusted p-value", {
  d <- simulate_gwas_data(M = 6L, causal = 1L, effect = 0.6, seed = 64)
  scan <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                     adjust = "effective", n_permutations = 10L, fdr = 0.2,
                     seed = 6L)
  perms <- draw_perms(6L, 10L, d$N)
  best <- function(p) apply(p, 1L, min)
  null_p <- permuted_pvalues(d, perms, A = d$A, X = d$X,
                             mode = "single_trait", adjust = "effective",
                             delta = scan$delta, which = "pvalues_adjust")
  ref <- expected_fdr(unname(best(scan$pvalues_adjust)), lapply(null_p, best),
                      0.2)
  expect_equal(unname(scan$qvalues_empirical), ref$q)
  expect_equal(scan$fdr_threshold, ref$threshold)
})

test_that("a strong causal SNP passes the FDR threshold", {
  d <- simulate_gwas_data(n_fam = 20, size = 5, causal = 1L, effect = 0.8,
                          M = 20L, seed = 71)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     n_permutations = 30L, fdr = 0.05, seed = 3L)
  expect_lte(scan$qvalues_empirical[[1]], 0.05)
  expect_gte(scan$fdr_threshold, scan$pvalues[[1]])
  # q-values are monotone in the p-values
  o <- order(scan$pvalues)
  expect_true(all(diff(scan$qvalues_empirical[o]) >= -1e-12))
})

test_that("fdr alone runs limmbo's ten permutations", {
  d <- simulate_gwas_data(M = 4L, seed = 65)
  scan <- indri_gwas(d$Y, d$geno, X = d$X, model = "lm", fdr = 0.25,
                     seed = 1L)
  expect_equal(scan$n_permutations, 10L)
  expect_length(scan$qvalues_empirical, 4L)
  expect_true("fdr_threshold" %in% names(scan))
})

test_that("permutations are reproducible and do not depend on cores", {
  skip_on_os("windows")
  d <- simulate_gwas_data(M = 6L, seed = 66)
  run <- function(cores, seed = 3L)
    indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
               n_permutations = 16L, fdr = 0.1, seed = seed, cores = cores)
  one <- run(1L)
  expect_equal(run(2L)$pvalues_empirical, one$pvalues_empirical)
  expect_equal(run(2L)$qvalues_empirical, one$qvalues_empirical)
  expect_equal(run(2L)$fdr_threshold, one$fdr_threshold)
  expect_equal(run(1L)$pvalues_empirical, one$pvalues_empirical)
  # and the caller's random stream is left alone
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  run(1L)
  expect_equal(stats::runif(1), before)
})

test_that("a strong causal SNP reaches the smallest empirical p-value", {
  d <- simulate_gwas_data(n_fam = 20, size = 5, causal = 1L, effect = 0.8,
                          M = 6L, seed = 67)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     n_permutations = 40L, seed = 2L)
  expect_equal(scan$pvalues_empirical[[1]], 1 / 41)
})

test_that("permutations follow monomorphic SNPs and dominance", {
  d <- simulate_gwas_data(M = 5L, seed = 68)
  d$geno[, 2] <- 1
  expect_warning(scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A,
                                    X = d$X, n_permutations = 10L, seed = 1L),
                 "monomorphic")
  expect_true(is.na(scan$pvalues_empirical[[2]]))
  expect_false(anyNA(scan$pvalues_empirical[-2]))

  d <- simulate_gwas_data(M = 5L, seed = 69)
  dom <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                    dominance = TRUE, n_permutations = 10L, seed = 1L)
  perms <- draw_perms(1L, 10L, d$N)
  null_p <- permuted_pvalues(d, perms, G = d$G, E = d$E, A = d$A, X = d$X,
                             delta = dom$delta, dominance = TRUE)
  count <- Reduce(`+`, lapply(null_p, function(p) p <= dom$pvalues))
  expect_equal(dom$pvalues_empirical, (count + 1) / 11)
})

test_that("permutation arguments are checked, and absent by default", {
  d <- simulate_gwas_data(M = 3L, seed = 70)
  plain <- indri_gwas(d$Y, d$geno, X = d$X, model = "lm")
  expect_null(plain$pvalues_empirical)
  expect_null(plain$fdr_threshold)
  expect_false("pvalue_empirical" %in% names(as.data.frame(plain)))
  expect_error(indri_gwas(d$Y, d$geno, X = d$X, model = "lm",
                          n_permutations = -1L), "n_permutations")
  expect_error(indri_gwas(d$Y, d$geno, X = d$X, model = "lm", fdr = 1.5),
               "fdr")
  with_p <- indri_gwas(d$Y, d$geno, X = d$X, model = "lm",
                       n_permutations = 5L, seed = 1L)
  df <- as.data.frame(with_p)
  expect_equal(df$pvalue_empirical, unname(with_p$pvalues_empirical))
})
