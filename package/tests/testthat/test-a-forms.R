test_that("the dense form is the matrix itself", {
  A <- diag(3)
  forms <- indri_a_forms(A, need = "dense")
  expect_equal(forms$dense, A)
  expect_null(forms$inverse)
})

test_that("the inverse form is marked for sommer", {
  A <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  forms <- indri_a_forms(A, need = "inverse")
  expect_equal(unclass(forms$inverse), unclass(solve(A)),
               ignore_attr = TRUE)
  expect_true(attr(forms$inverse, "inverse"))
})

test_that("the sparse inverse is a three-column ginv of the lower triangle", {
  A <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  forms <- indri_a_forms(A, need = "sparse_inverse")
  sparse <- forms$sparse_inverse

  expect_true("ginv" %in% class(sparse))
  expect_equal(colnames(sparse), c("Row", "Column", "Ainverse"))
  expect_true(all(sparse[, "Row"] >= sparse[, "Column"]))

  dense <- solve(A)
  for (i in seq_len(nrow(sparse))) {
    expect_equal(unname(sparse[i, "Ainverse"]),
                 unname(dense[sparse[i, "Row"], sparse[i, "Column"]]))
  }
  expect_equal(attr(sparse, "geneticGroups"), c(0, 0))
})

test_that("bad relationship matrices are refused", {
  expect_error(indri_a_forms(matrix(1:6, 2, 3), need = "dense"), "square")
  expect_error(indri_a_forms(matrix(c(1, NA, NA, 1), 2), need = "dense"),
               "missing")
})

test_that("the eigen form is the ascending decomposition stanimal takes", {
  A <- matrix(c(2, 0.5, 0.2, 0.5, 1, 0.1, 0.2, 0.1, 1.5), 3, 3)
  forms <- indri_a_forms(A, need = c("dense", "eigen"))
  e <- forms$eigen
  expect_false(is.unsorted(e$lambda))
  expect_equal(e$Q %*% (e$lambda * t(e$Q)), A, tolerance = 1e-12)
  expect_equal(crossprod(e$Q), diag(3), tolerance = 1e-12)
  expect_null(indri_a_forms(A, need = "dense")$eigen)
})

test_that("the stanimal engine shares the decomposition when stanimal takes it", {
  skip_if_not_installed("stanimal")
  engine <- indri_engine("stanimal")
  takes <- "eigen_A" %in% names(formals(stanimal::lmm_animal))
  expect_identical("eigen" %in% engine$capabilities$a_form, takes)
})
