test_that("prepare_matrix() builds a peptide x protein matrix of ones", {
  data <- data.frame(protein_id = c("A", "B", "B"), peptide_seq = c("p1", "p1", "p2"))

  mat <- as.matrix(prepare_matrix(data))

  expect_equal(mat, matrix(c(1, 0, 1, 1), 2, dimnames = list(c("p1", "p2"), c("A", "B"))))
})

test_that("prepare_matrix() gives the same matrix for a tibble and a data frame", {
  skip_if_not_installed("dplyr")
  data <- data.frame(protein_id = c("A", "B", "B"), peptide_seq = c("p1", "p1", "p2"))

  expect_equal(prepare_matrix(dplyr::as_tibble(data)), prepare_matrix(data))
})

test_that("prepare_matrix() weighting = 'inverse' splits a shared peptide between its proteins", {
  data <- data.frame(protein_id = c("A", "B", "B"), peptide_seq = c("p1", "p1", "p2"))

  mat <- as.matrix(prepare_matrix(data, weighting = "inverse"))

  expect_equal(mat, matrix(c(0.5, 0, 0.5, 1), 2, dimnames = list(c("p1", "p2"), c("A", "B"))))
})

test_that("prepare_matrix() uses numeric weights, one per row", {
  data <- data.frame(protein_id = c("A", "B"), peptide_seq = c("p1", "p1"))

  mat <- as.matrix(prepare_matrix(data, weighting = c(2, 3)))

  expect_equal(mat, matrix(c(2, 3), 1, dimnames = list("p1", c("A", "B"))))
})

test_that("prepare_matrix() keeps or sums repeated pairs according to use.last.ij", {
  data <- data.frame(protein_id = c("A", "A"), peptide_seq = c("p1", "p1"))

  expect_equal(as.matrix(prepare_matrix(data))[1, 1], 1)
  expect_equal(as.matrix(prepare_matrix(data, use.last.ij = FALSE))[1, 1], 2)
})

test_that("prepare_matrix() keeps or sums repeated pairs with numeric weights too", {
  data <- data.frame(protein_id = c("A", "A", "B"), peptide_seq = c("p1", "p1", "p1"))

  expect_equal(as.matrix(prepare_matrix(data, weighting = c(2, 3, 4)))[1, ], c(A = 3, B = 4))
  expect_equal(as.matrix(prepare_matrix(data, weighting = c(2, 3, 4), use.last.ij = FALSE))[1, ], c(A = 5, B = 4))
})

test_that("prepare_matrix() drops unused factor levels", {
  data(protpepmetashort, package = "prozor", envir = environment())
  first_rows <- protpepmetashort[1:50, ]
  # the first peptide level is unused, so no error but an empty first row without the fix
  without_first <- protpepmetashort[protpepmetashort$peptide_seq != levels(protpepmetashort$peptide_seq)[1], ]

  mat <- prepare_matrix(first_rows)
  mat_without_first <- prepare_matrix(without_first)

  expect_equal(dim(mat), c(length(unique(first_rows$peptide_seq)), length(unique(first_rows$protein_id))))
  expect_setequal(rownames(mat), as.character(first_rows$peptide_seq))
  expect_setequal(rownames(mat_without_first), as.character(without_first$peptide_seq))
  expect_setequal(colnames(mat_without_first), as.character(without_first$protein_id))
  expect_true(all(Matrix::rowSums(mat_without_first) > 0))
  expect_true(all(Matrix::colSums(mat_without_first) > 0))
})

test_that("prepare_matrix() rejects unknown weightings", {
  data <- data.frame(protein_id = c("A", "B"), peptide_seq = c("p1", "p1"))

  expect_error(prepare_matrix(data, weighting = "AA"), "`weighting` must be NULL.*not a character of length 1")
})
