test_that("prepare_matrix() builds a peptide x protein matrix of ones", {
  data <- data.frame(proteinID = c("A", "B", "B"), strippedSequence = c("p1", "p1", "p2"))

  mat <- as.matrix(prepare_matrix(data))

  expect_equal(mat, matrix(c(1, 0, 1, 1), 2, dimnames = list(c("p1", "p2"), c("A", "B"))))
})

test_that("prepare_matrix() gives the same matrix for a tibble and a data frame", {
  data <- data.frame(proteinID = c("A", "B", "B"), strippedSequence = c("p1", "p1", "p2"))

  expect_equal(prepare_matrix(dplyr::as_tibble(data)), prepare_matrix(data))
})

test_that("prepare_matrix() weighting = 'inverse' splits a shared peptide between its proteins", {
  data <- data.frame(proteinID = c("A", "B", "B"), strippedSequence = c("p1", "p1", "p2"))

  mat <- as.matrix(prepare_matrix(data, weighting = "inverse"))

  expect_equal(mat, matrix(c(0.5, 0, 0.5, 1), 2, dimnames = list(c("p1", "p2"), c("A", "B"))))
})

test_that("prepare_matrix() uses numeric weights, one per row", {
  data <- data.frame(proteinID = c("A", "B"), strippedSequence = c("p1", "p1"))

  mat <- as.matrix(prepare_matrix(data, weighting = c(2, 3)))

  expect_equal(mat, matrix(c(2, 3), 1, dimnames = list("p1", c("A", "B"))))
})

test_that("prepare_matrix() keeps or sums repeated pairs according to use.last.ij", {
  data <- data.frame(proteinID = c("A", "A"), strippedSequence = c("p1", "p1"))

  expect_equal(as.matrix(prepare_matrix(data))[1, 1], 1)
  expect_equal(as.matrix(prepare_matrix(data, use.last.ij = FALSE))[1, 1], 2)
})

test_that("prepare_matrix() rejects unknown weightings", {
  data <- data.frame(proteinID = c("A", "B"), strippedSequence = c("p1", "p1"))

  expect_error(prepare_matrix(data, weighting = "AA"), "`weighting` must be NULL.*not a character of length 1")
})
