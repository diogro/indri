# indri_gwas(gxe = ): each SNP also interacts with environments named among
# the columns of X. `pvalues` is the joint test of the SNP and its
# interactions (and dominance), `pvalues_gxe` the interactions given the
# SNP's main effect (and dominance).

gxe_data <- function(seed, M = 6L, ...) {
  d <- simulate_gwas_data(M = M, seed = seed, ...)
  set.seed(seed + 1000)
  d$X <- cbind(d$X, sex = stats::rbinom(d$N, 1L, 0.5))
  d
}

test_that("the multi-trait LMM tests the SNP and its interaction", {
  d <- gxe_data(101)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 0.8, gxe = "x")
  expect_equal(scan$gxe_terms, "x")
  for (j in seq_len(d$M)) {
    g <- d$geno[, j]
    ge <- g * d$X[, "x"]
    joint <- dense_nested_multi(d$Y, d$X, cbind(d$X, g, ge), E = d$E,
                                G = d$G, A = d$A, delta = 0.8)
    inter <- dense_nested_multi(d$Y, cbind(d$X, g), cbind(d$X, g, ge),
                                E = d$E, G = d$G, A = d$A, delta = 0.8)
    expect_rel_equal(scan$pvalues[j], joint$p, 1e-8, label = paste("joint", j))
    expect_rel_equal(scan$pvalues_gxe[j], inter$p, 1e-8,
                     label = paste("gxe", j))
    expect_equal(unname(scan$a[j, ]), joint$beta[1, ], tolerance = 1e-8)
    expect_equal(unname(scan$gxe[j, , "x"]), joint$beta[2, ], tolerance = 1e-8)
  }
  expect_equal(dim(scan$gxe), c(d$M, d$P, 1L))
  expect_equal(dimnames(scan$gxe), list(colnames(d$geno), colnames(d$Y), "x"))
})

test_that("several environments enter together", {
  d <- gxe_data(102)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 1, gxe = c("x", "sex"))
  g <- d$geno[, 3]
  ge <- cbind(g * d$X[, "x"], g * d$X[, "sex"])
  joint <- dense_nested_multi(d$Y, d$X, cbind(d$X, g, ge), E = d$E, G = d$G,
                              A = d$A, delta = 1)
  inter <- dense_nested_multi(d$Y, cbind(d$X, g), cbind(d$X, g, ge),
                              E = d$E, G = d$G, A = d$A, delta = 1)
  expect_rel_equal(scan$pvalues[3], joint$p, 1e-8)
  expect_rel_equal(scan$pvalues_gxe[3], inter$p, 1e-8)
  expect_equal(unname(scan$gxe[3, , "sex"]), joint$beta[3, ], tolerance = 1e-8)
})

test_that("the multi-trait LM and the single-trait tests take gxe", {
  d <- gxe_data(103)
  lm <- indri_gwas(d$Y, d$geno, E = d$E, X = d$X, model = "lm", gxe = "sex")
  g <- d$geno[, 2]
  ge <- g * d$X[, "sex"]
  expect_rel_equal(lm$pvalues_gxe[2],
                   dense_nested_multi(d$Y, cbind(d$X, g), cbind(d$X, g, ge),
                                      E = d$E)$p, 1e-8)
  delta <- c(wing = 0.6, tail = 1.5, beak = 1)
  single <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                       delta = delta, gxe = "x")
  for (k in seq_len(d$P)) for (j in seq_len(d$M)) {
    g <- d$geno[, j]
    ge <- g * d$X[, "x"]
    joint <- dense_nested_single(d$Y[, k], d$X, cbind(d$X, g, ge), A = d$A,
                                 delta = delta[[k]])
    inter <- dense_nested_single(d$Y[, k], cbind(d$X, g), cbind(d$X, g, ge),
                                 A = d$A, delta = delta[[k]])
    expect_rel_equal(single$pvalues[j, k], joint$p, 1e-8)
    expect_rel_equal(single$pvalues_gxe[j, k], inter$p, 1e-8)
    expect_equal(c(single$a[j, k], single$gxe[j, k, "x"]), joint$beta,
                 tolerance = 1e-8)
  }
  expect_equal(dim(single$pvalues_gxe), c(d$M, d$P))
})

test_that("dominance and gxe together: GxE on the additive term only", {
  d <- gxe_data(104)
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 1, gxe = "x", dominance = TRUE)
  for (j in which(scan$dominance_estimable)) {
    g <- d$geno[, j]
    h <- heterozygous(g)
    ge <- g * d$X[, "x"]
    full <- cbind(d$X, g, h, ge)
    joint <- dense_nested_multi(d$Y, d$X, full, E = d$E, G = d$G, A = d$A,
                                delta = 1)
    inter <- dense_nested_multi(d$Y, cbind(d$X, g, h), full, E = d$E,
                                G = d$G, A = d$A, delta = 1)
    dom <- dense_nested_multi(d$Y, cbind(d$X, g, ge), full, E = d$E,
                              G = d$G, A = d$A, delta = 1)
    expect_rel_equal(scan$pvalues[j], joint$p, 1e-8)
    expect_rel_equal(scan$pvalues_gxe[j], inter$p, 1e-8)
    expect_rel_equal(scan$pvalues_dominance[j], dom$p, 1e-8)
    expect_equal(unname(scan$d[j, ]), joint$beta[2, ], tolerance = 1e-8)
    expect_equal(unname(scan$gxe[j, , "x"]), joint$beta[3, ],
                 tolerance = 1e-8)
  }
})

test_that("an interaction that cannot be estimated keeps the main-effect test", {
  d <- gxe_data(105)
  # carriers of SNP 1 are all males: its interaction with sex is its dosage
  d$geno[d$X[, "sex"] == 0, 1] <- 0L
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     delta = 1, gxe = "sex")
  plain <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                      delta = 1)
  expect_false(scan$gxe_estimable[[1]])
  expect_true(all(scan$gxe_estimable[-1]))
  expect_named(scan$gxe_estimable, colnames(d$geno))
  expect_equal(scan$pvalues[[1]], plain$pvalues[[1]], tolerance = 1e-10)
  expect_equal(scan$a[1, ], plain$a[1, ], tolerance = 1e-10)
  expect_true(is.na(scan$pvalues_gxe[[1]]))
  expect_true(all(is.na(scan$gxe[1, , ])))
})

test_that("gxe names columns of X that can interact", {
  d <- gxe_data(106)
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                          gxe = "diet"), "diet")
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                          gxe = "intercept"), "constant")
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A,
                          X = unname(d$X), gxe = "x"), "names")
  expect_error(indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                          gxe = 2), "names")
})

test_that("a SNP whose effect depends on the environment is found", {
  d <- gxe_data(107, M = 10L, causal = 0L, n_fam = 30, size = 5)
  sex <- d$X[, "sex"]
  d$Y <- d$Y + outer(d$geno[, 1] * (2 * sex - 1), rep(0.5, d$P))
  scan <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                     gxe = "sex")
  expect_equal(unname(which.min(scan$pvalues_gxe)), 1L)
  expect_lt(scan$pvalues_gxe[[1]], 1e-4)
  plain <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X)
  expect_gt(plain$pvalues[[1]], scan$pvalues[[1]])
})

test_that("gxe results reach the outputs, plots and permutations", {
  d <- gxe_data(108)
  multi <- indri_gwas(d$Y, d$geno, G = d$G, E = d$E, A = d$A, X = d$X,
                      delta = 1, gxe = "sex", n_permutations = 5L, seed = 1L)
  df <- as.data.frame(multi)
  expect_equal(df$pvalue_gxe, unname(multi$pvalues_gxe))
  expect_true(all(paste0("gxe_sex_", colnames(d$Y)) %in% names(df)))
  expect_equal(df$gxe_estimable, unname(multi$gxe_estimable))
  expect_output(print(multi), "sex")
  expect_length(multi$pvalues_empirical, d$M)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(multi, which = "gxe"), "ggplot")
  single <- indri_gwas(d$Y, d$geno, A = d$A, X = d$X, mode = "single_trait",
                       gxe = "sex")
  long <- as.data.frame(single)
  expect_true(all(c("pvalue_gxe", "gxe_sex") %in% names(long)))
  plain <- indri_gwas(d$Y, d$geno, X = d$X, model = "lm")
  expect_null(plain$pvalues_gxe)
  expect_error(plot(plain, which = "gxe"), "gxe")
})
