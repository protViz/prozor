test_that("prepareMatrix() builds a peptide x protein matrix of ones", {
  data <- data.frame(proteinID = c("A", "B", "B"), strippedSequence = c("p1", "p1", "p2"))

  mat <- as.matrix(prepareMatrix(data))

  expect_equal(mat, matrix(c(1, 0, 1, 1), 2, dimnames = list(c("p1", "p2"), c("A", "B"))))
})

test_that("prepareMatrix() weighting = 'inverse' splits a shared peptide between its proteins", {
  data <- data.frame(proteinID = c("A", "B", "B"), strippedSequence = c("p1", "p1", "p2"))

  mat <- as.matrix(prepareMatrix(data, weighting = "inverse"))

  expect_equal(mat, matrix(c(0.5, 0, 0.5, 1), 2, dimnames = list(c("p1", "p2"), c("A", "B"))))
})

test_that("prepareMatrix() uses numeric weights, one per row", {
  data <- data.frame(proteinID = c("A", "B"), strippedSequence = c("p1", "p1"))

  mat <- as.matrix(prepareMatrix(data, weighting = c(2, 3)))

  expect_equal(mat, matrix(c(2, 3), 1, dimnames = list("p1", c("A", "B"))))
})

test_that("prepareMatrix() keeps or sums repeated pairs according to use.last.ij", {
  data <- data.frame(proteinID = c("A", "A"), strippedSequence = c("p1", "p1"))

  expect_equal(as.matrix(prepareMatrix(data))[1, 1], 1)
  expect_equal(as.matrix(prepareMatrix(data, use.last.ij = FALSE))[1, 1], 2)
})

test_that("prepareMatrix() rejects unknown weightings", {
  data <- data.frame(proteinID = c("A", "B"), strippedSequence = c("p1", "p1"))

  expect_error(prepareMatrix(data, weighting = "AA"), "weighting must be NULL")
})
