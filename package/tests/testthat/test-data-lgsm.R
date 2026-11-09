# The LG/SM F3 data sets: shape, coding, and that they line up with each
# other and with indri_gwas().

test_that("the phenotypes hold every trait class", {
  expect_equal(nrow(lgsm_f3), 1611L)
  expect_false(anyDuplicated(lgsm_f3$ID) > 0)
  expect_named(lgsm_traits, c("growth", "weight", "area", "necropsy"))
  expect_true(all(unlist(lgsm_traits) %in% names(lgsm_f3)))
  expect_false(anyNA(lgsm_f3[lgsm_traits$weight]))
  # growth is the weekly gain in weight
  expect_length(lgsm_traits$growth, 7L)
  expect_equal(lgsm_f3$growth34, lgsm_f3$WEEK4 - lgsm_f3$WEEK3)
  expect_equal(unname(as.matrix(lgsm_f3[lgsm_traits$growth])),
               unname(as.matrix(lgsm_f3[paste0("WEEK", 2:8)]) -
                        as.matrix(lgsm_f3[paste0("WEEK", 1:7)])))
  expect_equal(sum(!stats::complete.cases(lgsm_f3[lgsm_traits$area])), 95L)
  expect_true(all(lgsm_f3$SEX %in% c("F", "M")))
  expect_equal(length(unique(lgsm_f3$FAMILY)), 201L)
})

test_that("the genotypes follow the additive, dominance, imprinting design", {
  A <- lgsm_genotypes$A
  expect_named(lgsm_genotypes, c("A", "D", "I"))
  expect_equal(dim(A), c(1934L, 353L))
  for (m in lgsm_genotypes) expect_identical(dimnames(m), dimnames(A))
  expect_true(all(A %in% c(-1L, 0L, 1L)))
  expect_identical(lgsm_genotypes$D, (A == 0L) * 1L)
  expect_true(all((lgsm_genotypes$I != 0L) == (A == 0L)))
  # F3 mice and their F2 parents
  f3 <- rownames(A) %in% lgsm_f3$ID
  expect_equal(sum(f3), 1552L)
  expect_true(all(rownames(A)[!f3] %in% c(lgsm_f3$Sire, lgsm_f3$Dam)))
})

test_that("the map describes the genotype columns", {
  expect_equal(lgsm_markers$snp, colnames(lgsm_genotypes$A))
  expect_equal(as.vector(table(lgsm_markers$chrom)),
               c(31, 28, 26, 24, 19, 22, 18, 13, 19, 18, 21, 19, 18, 14, 16,
                 13, 13, 13, 8))
  expect_true(all(tapply(lgsm_markers$cM, lgsm_markers$chrom,
                         function(x) all(diff(x) >= 0))))
  expect_equal(lgsm_markers$snp[is.na(lgsm_markers$bp)], "rs3684333")
})

test_that("the pedigree has the F2 as unrelated founders", {
  f2 <- lgsm_pedigree[lgsm_pedigree$generation == "F2", ]
  f3 <- lgsm_pedigree[lgsm_pedigree$generation == "F3", ]
  expect_equal(nrow(f2), 391L)
  expect_true(all(is.na(f2$Sire) & is.na(f2$Dam)))
  expect_equal(f3$ID, lgsm_f3$ID)
  expect_true(all(c(f3$Sire, f3$Dam) %in% c(f2$ID, NA)))
  expect_equal(sum(is.na(f3$Sire)), 4L)
})

test_that("A is the F3 relationship matrix from that pedigree", {
  expect_equal(dim(lgsm_A), c(1611L, 1611L))
  expect_equal(rownames(lgsm_A), as.character(lgsm_f3$ID))
  expect_true(isSymmetric(lgsm_A))
  expect_true(all(diag(lgsm_A) == 1))
  expect_true(all(lgsm_A[upper.tri(lgsm_A)] %in% c(0, 0.25, 0.5)))
  id <- function(i) as.character(lgsm_f3$ID[i])
  # full sibs, half sibs through a sire with two mates, unrelated
  sibs <- which(lgsm_f3$FAMILY == lgsm_f3$FAMILY[10] & !is.na(lgsm_f3$FAMILY))
  expect_equal(lgsm_A[id(sibs[1]), id(sibs[2])], 0.5)
  pairs <- unique(lgsm_f3[!is.na(lgsm_f3$Sire), c("Sire", "Dam")])
  sire <- pairs$Sire[duplicated(pairs$Sire)][1]
  dams <- pairs$Dam[pairs$Sire == sire]
  a <- which(lgsm_f3$Dam == dams[1])[1]
  b <- which(lgsm_f3$Dam == dams[2])[1]
  expect_equal(lgsm_A[id(a), id(b)], 0.25)
  other <- which(lgsm_f3$Sire != lgsm_f3$Sire[a] & lgsm_f3$Dam != lgsm_f3$Dam[a])[1]
  expect_equal(lgsm_A[id(a), id(other)], 0)
  # the three mice of family 199 (sire unrecorded) are full sibs
  f199 <- id(which(lgsm_f3$FAMILY == 199))
  expect_true(all(lgsm_A[f199, f199][upper.tri(diag(3))] == 0.5))
  expect_gt(min(eigen(lgsm_A, symmetric = TRUE, only.values = TRUE)$values),
            0.49)
})

test_that("the data drive indri_gwas(), with dominance", {
  keep <- stats::complete.cases(lgsm_f3[lgsm_traits$necropsy]) &
    lgsm_f3$ID %in% rownames(lgsm_genotypes$A)
  d <- lgsm_f3[keep, ]
  Y <- as.matrix(d[lgsm_traits$necropsy])
  rownames(Y) <- d$ID
  X <- cbind(intercept = 1, male = as.numeric(d$SEX == "M"))
  rownames(X) <- d$ID
  chr1 <- lgsm_markers$chrom == 1
  info <- transform(lgsm_markers[chr1, ], pos = cM)
  scan <- indri_gwas(Y, lgsm_genotypes$A[, chr1] + 1, X = X, model = "lm",
                     snp_info = info, dominance = TRUE)
  lmm <- indri_gwas(Y, lgsm_genotypes$A[, chr1] + 1, A = lgsm_A, X = X,
                    mode = "single_trait", snp_info = info)
  expect_false(anyNA(lmm$pvalues))
  expect_equal(length(scan$pvalues), sum(chr1))
  expect_true(all(scan$dominance_estimable))
  expect_false(anyNA(scan$pvalues))
})
