# indri_gemma(): GEMMA from R. The runs need the GEMMA program
# (options(indri.gemma = <path>) or `gemma` on the PATH); the relationship
# matrices and the argument checks do not.

gemma_path <- function() {
  p <- getOption("indri.gemma") %||% unname(Sys.which("gemma"))
  if (nzchar(p) && file.exists(p)) p else NULL
}
skip_if_no_gemma <- function() {
  if (is.null(gemma_path())) skip("GEMMA is not available")
}

sim_gemma_data <- function(N = 240, M = 240, P = 2, seed = 1) {
  set.seed(seed)
  chrom <- rep(c("1", "2", "3"), each = M / 3)
  # three founder haplotypes per chromosome segment, so SNPs are correlated
  geno <- matrix(rbinom(N * M, 2, rep(runif(M, 0.15, 0.85), each = N)), N, M,
                 dimnames = list(paste0("id", seq_len(N)), paste0("snp", seq_len(M))))
  causal <- round(M * c(1, 10, 20) / 24)
  B <- matrix(c(1, 0.8, 0, 0.8, -0.8, 0.8), 3, P)[, seq_len(P), drop = FALSE]
  K <- indri_gemma_kinship(geno)
  L <- t(chol(K + diag(1e-6, N)))
  U <- L %*% matrix(rnorm(N * P), N, P) %*% chol(matrix(c(1, 0.5, 0.5, 1), 2)[seq_len(P), seq_len(P)])
  Y <- scale(geno[, causal], scale = FALSE) %*% B + U + matrix(rnorm(N * P), N, P)
  dimnames(Y) <- list(rownames(geno), paste0("t", seq_len(P)))
  X <- cbind(intercept = 1, sex = rep(0:1, length.out = N))
  Y[, 1] <- Y[, 1] + 0.5 * X[, "sex"]
  list(Y = Y, geno = geno, X = X, causal = colnames(geno)[causal],
       snp_info = data.frame(snp = colnames(geno), chrom = chrom, pos = seq_len(M) * 1e5))
}

test_that("indri_gemma_kinship() is GEMMA's centred relationship matrix", {
  d <- sim_gemma_data()
  Z <- scale(d$geno, scale = FALSE)
  expect_equal(unname(indri_gemma_kinship(d$geno)), unname(tcrossprod(Z) / ncol(Z)))
})

test_that("LOCO matrices equal the matrix of the other chromosomes' SNPs", {
  d <- sim_gemma_data()
  loco <- indri_gemma_kinship(d$geno, d$snp_info$chrom)
  expect_named(loco, c("1", "2", "3"))
  other <- d$snp_info$chrom != "2"
  expect_equal(loco[["2"]], indri_gemma_kinship(d$geno[, other]))
})

test_that("monomorphic SNPs do not enter the relationship matrix", {
  d <- sim_gemma_data()
  g <- cbind(d$geno, mono = 1)
  expect_equal(indri_gemma_kinship(g), indri_gemma_kinship(d$geno))
})

test_that("indri_gemma() checks its inputs before running GEMMA", {
  d <- sim_gemma_data()
  expect_error(indri_gemma(d$Y, d$geno, gemma = tempfile()), "needs the GEMMA program")
  skip_if_no_gemma()
  Y <- d$Y; Y[1, 1] <- NA
  expect_error(indri_gemma(Y, d$geno), "missing values")
  expect_error(indri_gemma(d$Y, d$geno, X = d$X[, "sex", drop = FALSE]), "constant")
  expect_error(indri_gemma(d$Y, d$geno, loco = TRUE), "chrom")
  expect_error(indri_gemma(d$Y[, 1, drop = FALSE], d$geno), "at least two traits")
})

test_that("multi-trait GEMMA agrees with indri_gwas() given GEMMA's covariances", {
  skip_if_no_gemma()
  d <- sim_gemma_data()
  K <- indri_gemma_kinship(d$geno)
  gm <- indri_gemma(d$Y, d$geno, A = K, X = d$X, snp_info = d$snp_info)
  expect_s3_class(gm, "indri_gwas")
  expect_equal(gm$engine, "gemma")
  expect_equal(dim(gm$a), c(ncol(d$geno), 2L))
  nl <- gm$null[[1]]
  expect_true(isSymmetric(nl$Vg)); expect_true(isSymmetric(nl$Ve))
  ig <- indri_gwas(d$Y, d$geno, G = nl$Vg, E = nl$Ve, A = K, X = d$X,
                   snp_info = d$snp_info, delta = 1)
  ok <- !is.na(gm$pvalues)
  # GEMMA's Wald test vs indri's likelihood ratio test, same covariances
  expect_gt(cor(-log10(gm$pvalues[ok]), -log10(ig$pvalues[ok])), 0.99)
  expect_equal(unname(gm$a[ok, ]), unname(ig$a[ok, ]), tolerance = 0.02)
  # the causal SNPs are found
  expect_true(all(-log10(gm$pvalues[d$causal]) > 4))
})

test_that("single-trait GEMMA returns one column per trait", {
  skip_if_no_gemma()
  d <- sim_gemma_data()
  gs <- indri_gemma(d$Y, d$geno, X = d$X, snp_info = d$snp_info, mode = "single_trait")
  expect_equal(dim(gs$pvalues), c(ncol(d$geno), 2L))
  expect_equal(colnames(gs$pvalues), colnames(d$Y))
  expect_length(gs$null, 2L)
  K <- indri_gemma_kinship(d$geno)
  ig <- indri_gwas(d$Y, d$geno, A = K, X = d$X, mode = "single_trait",
                   snp_info = d$snp_info)
  ok <- !is.na(gs$pvalues[, 1])
  expect_gt(cor(-log10(gs$pvalues[ok, 1]), -log10(ig$pvalues[ok, 1])), 0.99)
})

test_that("LOCO runs give one null model per chromosome and work with the methods", {
  skip_if_no_gemma()
  d <- sim_gemma_data()
  gl <- indri_gemma(d$Y, d$geno, X = d$X, snp_info = d$snp_info, loco = TRUE, cores = 2)
  expect_named(gl$null, c("1", "2", "3"))
  expect_false(anyNA(gl$pvalues))
  s <- summary(gl)
  expect_s3_class(s, "summary.indri_gwas")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(gl), "ggplot")
  lg <- indri_locus_G(gl$a[d$causal, ], genotypes = d$geno[, d$causal])
  expect_equal(dim(lg$G), c(2L, 2L))
})

test_that("keep_dir keeps one directory per run", {
  skip_if_no_gemma()
  d <- sim_gemma_data(M = 60)
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE))
  indri_gemma(d$Y, d$geno, X = d$X, snp_info = d$snp_info, loco = TRUE, keep_dir = dir)
  expect_equal(sort(list.files(dir)), paste0("run_", 1:3))
  expect_true(file.exists(file.path(dir, "run_1", "output", "result.assoc.txt")))
})
