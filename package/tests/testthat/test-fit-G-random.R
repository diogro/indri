# fit_G(random = ): extra random effects (nurse, dam) next to G.

sim_nurse <- function(n_fam = 60, per_fam = 8, seed = 21) {
  set.seed(seed)
  N <- n_fam * per_fam
  fam <- rep(seq_len(n_fam), each = per_fam)
  A <- outer(fam, fam, function(a, b) ifelse(a == b, 0.5, 0)); diag(A) <- 1
  ids <- paste0("i", seq_len(N)); dimnames(A) <- list(ids, ids)
  # cross-fostering: half of each family goes to another family's nurse
  nurse <- fam
  swap <- rep(c(TRUE, FALSE), length.out = N)
  nurse[swap] <- ((fam[swap] + 17) %% n_fam) + 1
  G <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3)
  C <- matrix(c(1.2, 0.8, 0, 0.8, 1, 0, 0, 0, 0.05), 3)       # nurse effects on traits 1-2
  E <- diag(3)
  U <- t(chol(A)) %*% matrix(rnorm(N * 3), N, 3) %*% chol(G)
  W <- matrix(rnorm(n_fam * 3), n_fam, 3) %*% chol(C)
  Y <- U + W[nurse, ] + matrix(rnorm(N * 3), N, 3) %*% chol(E)
  dimnames(Y) <- list(ids, paste0("t", 1:3))
  list(Y = Y, A = A, nurse = nurse, G = G, C = C, E = E)
}

test_that("sommer recovers G and the nurse covariance", {
  skip_if_not_installed("sommer")
  d <- sim_nurse()
  sub <- t(utils::combn(3, 2))
  fit <- fit_G(d$Y, A = d$A, engine = "sommer", subsets = sub, random = list(nurse = d$nurse),
               cores = 1, progress = FALSE)
  expect_named(fit$random, "nurse")
  expect_equal(dim(fit$random$nurse), c(3L, 3L))
  expect_equal(unname(fit$random$nurse[1:2, 1:2]), d$C[1:2, 1:2], tolerance = 0.35)
  expect_equal(unname(diag(fit$G)), diag(d$G), tolerance = 0.35)
  expect_equal(unname(diag(fit$E)), diag(d$E), tolerance = 0.3)
  # leaving the nurse out pushes its variance into G (sibs share nurses)
  naive <- fit_G(d$Y, A = d$A, engine = "sommer", subsets = sub, cores = 1, progress = FALSE)
  expect_gt(naive$G[1, 2], fit$G[1, 2])
  expect_null(naive$random)
})

test_that("random terms are checked", {
  d <- sim_nurse(n_fam = 10, per_fam = 4)
  sub <- t(utils::combn(3, 2))
  expect_error(fit_G(d$Y, A = d$A, engine = "stanimal", subsets = sub,
                     random = list(nurse = d$nurse)), "sommer or asreml")
  skip_if_not_installed("sommer")
  expect_error(fit_G(d$Y, A = d$A, engine = "sommer", subsets = sub,
                     random = list(d$nurse)), "named list")
  expect_error(fit_G(d$Y, A = d$A, engine = "sommer", subsets = sub,
                     random = list(nurse = d$nurse[-1])), "individuals")
})

test_that("asreml gives the same nurse covariance as sommer", {
  skip_if_no_asreml()
  d <- sim_nurse()
  sub <- t(utils::combn(3, 2))
  a <- fit_G(d$Y, A = d$A, engine = "asreml", subsets = sub, random = list(nurse = d$nurse),
             cores = 1, progress = FALSE)
  s <- fit_G(d$Y, A = d$A, engine = "sommer", subsets = sub, random = list(nurse = d$nurse),
             cores = 1, progress = FALSE)
  expect_equal(a$random$nurse, s$random$nurse, tolerance = 0.05)
  expect_equal(a$G, s$G, tolerance = 0.05)
})
