# indri_gwas(dominance = TRUE): each SNP enters as its additive dosage and a
# heterozygote indicator. `pvalues` is the joint test of both, and
# `pvalues_dominance` tests dominance with the additive term in the null.

test_that("the multi-trait LMM tests additive and dominance jointly", {
  d <- simulate_gwas_data(seed = 51)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 0.8, dominance = TRUE)
  expect_true(scan$dominance)
  # SNPs seen in only two genotype classes make the dense joint model
  # singular; they are tested separately below
  testable <- which(scan$dominance_estimable)
  expect_gte(length(testable), 4L)
  for (j in testable) {
    a <- d$geno[, j]
    h <- heterozygous(a)
    joint <- dense_nested_multi(d$Y, d$X, cbind(d$X, a, h), E = d$E,
                                G = d$G, A = d$A, delta = 0.8)
    dom <- dense_nested_multi(d$Y, cbind(d$X, a), cbind(d$X, a, h), E = d$E,
                              G = d$G, A = d$A, delta = 0.8)
    expect_rel_equal(scan$pvalues[j], joint$p, 1e-8, label = paste("joint", j))
    expect_rel_equal(scan$pvalues_dominance[j], dom$p, 1e-8,
                     label = paste("dominance", j))
    expect_equal(unname(scan$a[j, ]), joint$beta[1, ], tolerance = 1e-8)
    expect_equal(unname(scan$d[j, ]), joint$beta[2, ],
                 tolerance = 1e-8)
  }
})

test_that("delta is fitted without the SNP, as without dominance", {
  d <- simulate_gwas_data(seed = 52)
  plain <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X)
  dom <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                    dominance = TRUE)
  expect_equal(dom$delta, plain$delta)
})

test_that("the multi-trait LM tests additive and dominance jointly", {
  d <- simulate_gwas_data(seed = 53)
  scan <- indri_gwas(d$Y, d$geno, E = d$E, X = d$X, model = "lm",
                     dominance = TRUE)
  a <- d$geno[, 2]
  h <- heterozygous(a)
  joint <- dense_nested_multi(d$Y, d$X, cbind(d$X, a, h), E = d$E)
  dom <- dense_nested_multi(d$Y, cbind(d$X, a), cbind(d$X, a, h), E = d$E)
  expect_rel_equal(scan$pvalues[2], joint$p, 1e-8)
  expect_rel_equal(scan$pvalues_dominance[2], dom$p, 1e-8)
})

test_that("the single-trait LMM and LM test additive and dominance jointly", {
  d <- simulate_gwas_data(seed = 54)
  delta <- c(wing = 0.6, tail = 1.5, beak = 1)
  lmm <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                    delta = delta, dominance = TRUE)
  lm <- indri_gwas(d$Y, d$geno, X = d$X, mode = "single_trait",
                   model = "lm", dominance = TRUE)
  expect_gte(sum(lmm$dominance_estimable), 4L)
  for (k in seq_len(d$P)) for (j in which(lmm$dominance_estimable)) {
    a <- d$geno[, j]
    h <- heterozygous(a)
    joint <- dense_nested_single(d$Y[, k], d$X, cbind(d$X, a, h), A = d$A,
                                 delta = delta[[k]])
    dom <- dense_nested_single(d$Y[, k], cbind(d$X, a), cbind(d$X, a, h),
                               A = d$A, delta = delta[[k]])
    expect_rel_equal(lmm$pvalues[j, k], joint$p, 1e-8)
    expect_rel_equal(lmm$pvalues_dominance[j, k], dom$p, 1e-8)
    expect_equal(c(lmm$a[j, k], lmm$d[j, k]), joint$beta,
                 tolerance = 1e-8)
    ols <- dense_nested_single(d$Y[, k], d$X, cbind(d$X, a, h))
    expect_rel_equal(lm$pvalues[j, k], ols$p, 1e-8)
  }
})

test_that("a SNP without all three genotypes keeps an additive-only test", {
  d <- simulate_gwas_data(seed = 55)
  geno <- d$geno
  geno[geno[, 1] == 2, 1] <- 1          # no minor homozygotes: h equals a
  d$geno <- geno
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 1, dominance = TRUE)
  additive <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                         delta = 1)
  expect_false(scan$dominance_estimable[[1]])
  expect_true(all(scan$dominance_estimable[-1]))
  expect_named(scan$dominance_estimable, colnames(d$geno))
  expect_equal(scan$pvalues[[1]], additive$pvalues[[1]], tolerance = 1e-10)
  expect_equal(scan$a[1, ], additive$a[1, ], tolerance = 1e-10)
  expect_true(is.na(scan$pvalues_dominance[[1]]))
  expect_true(all(is.na(scan$d[1, ])))
})

test_that("missing genotypes are imputed in both codings", {
  d <- simulate_gwas_data(seed = 56)
  geno <- d$geno
  geno[c(2, 9), 3] <- NA
  d$geno <- geno
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 1, dominance = TRUE)
  seen <- !is.na(geno[, 3])
  a <- ifelse(seen, geno[, 3], mean(geno[seen, 3]))
  h <- ifelse(seen, heterozygous(geno[, 3]), mean(heterozygous(geno[seen, 3])))
  joint <- dense_nested_multi(d$Y, d$X, cbind(d$X, a, h), E = d$E, G = d$G,
                              A = d$A, delta = 1)
  expect_rel_equal(scan$pvalues[3], joint$p, 1e-8)
})

test_that("dominance needs genotypes coded 0, 1 and 2", {
  d <- simulate_gwas_data(seed = 57)
  dosage <- d$geno + 0.1
  expect_error(indri_gwas(d$Y, dosage, G = d$G, E = d$E, A = d$A, X = d$X,
                          dominance = TRUE), "0, 1 and 2")
  # without dominance, dosages are fine
  expect_silent(indri_gwas(d$Y, dosage, G = d$G, E = d$E, A = d$A, X = d$X))
})

test_that("a dominant SNP is found by the dominance test", {
  set.seed(58)
  d <- simulate_gwas_data(n_fam = 40, size = 5, causal = 0L, M = 8L,
                          seed = 58)
  d$geno[, 1] <- stats::rbinom(d$N, 2L, 0.5)
  d$Y <- d$Y + heterozygous(d$geno[, 1]) * 0.8   # overdominance, no additive
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     dominance = TRUE)
  expect_equal(unname(which.min(scan$pvalues_dominance)), 1L)
  expect_lt(scan$pvalues_dominance[[1]], 1e-4)
  additive <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X)
  expect_gt(additive$pvalues[[1]], scan$pvalues[[1]])
})

test_that("dominance results reach the data frame and the print", {
  d <- simulate_gwas_data(seed = 59)
  multi <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                      delta = 1, dominance = TRUE)
  df <- as.data.frame(multi)
  expect_equal(df$pvalue_dominance, unname(multi$pvalues_dominance))
  expect_true(all(paste0("d_", colnames(d$Y)) %in% names(df)))
  expect_equal(df$dominance_estimable, unname(multi$dominance_estimable))
  expect_output(print(multi), "dominance")

  single <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                       dominance = TRUE)
  long <- as.data.frame(single)
  expect_true(all(c("pvalue_dominance", "d") %in% names(long)))
  expect_equal(long$d[long$trait == "beak"],
               unname(single$d[, "beak"]))

  plain <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait")
  expect_false(plain$dominance)
  expect_null(plain$pvalues_dominance)
  expect_false("pvalue_dominance" %in% names(as.data.frame(plain)))
})
