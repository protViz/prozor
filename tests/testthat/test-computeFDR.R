test_that("computeFDR() matches hand-computed FDR and q-values", {
  res <- computeFDR(score = c(4, 3, 2, 1), decoy_hit = c(FALSE, FALSE, TRUE, FALSE), larger_better = TRUE)

  expect_equal(res$score, c(4, 3, 2, 1))
  expect_equal(res$SimpleFDR, c(0, 0, 1 / 2, 1 / 3))
  expect_equal(res$qValue_SimpleFDR, c(0, 0, 1 / 3, 1 / 3))
  expect_equal(res$FPR, c(0, 0, 2 / 3, 2 / 4))
  expect_equal(res$qValue_FPR, c(0, 0, 1 / 2, 1 / 2))
})

test_that("computeFDR() sorts by score in the direction of larger_better", {
  data(fdrSample, package = "prozor", envir = environment())
  decoy <- grepl("REV_", fdrSample$proteinID)

  smaller <- computeFDR(fdrSample$score, decoy, larger_better = FALSE)
  larger <- computeFDR(fdrSample$score2, decoy, larger_better = TRUE)

  expect_false(is.unsorted(smaller$score))
  expect_false(is.unsorted(rev(larger$score)))
  for (res in list(smaller, larger)) {
    expect_length(res$score, nrow(fdrSample))
    expect_equal(res$decoy_hit, decoy[res$order])
    expect_false(is.unsorted(res$qValue_FPR))
    expect_false(is.unsorted(res$qValue_SimpleFDR))
    expect_true(all(res$qValue_FPR >= 0))
    expect_true(all(res$qValue_FPR <= res$FPR))
  }
})

test_that("computeFDRwithID() keeps IDs aligned with the sorted scores", {
  data(fdrSample, package = "prozor", envir = environment())

  res <- computeFDRwithID(fdrSample$score, fdrSample$proteinID, larger_better = FALSE)

  expect_equal(res$ID, fdrSample$proteinID[res$order])
  expect_equal(res$score, fdrSample$score[res$order])
  expect_equal(res$decoy_hit, grepl("REV_", res$ID))
})

test_that("predictScoreFDR() returns a score cut-off and rejects unknown methods", {
  data(fdrSample, package = "prozor", envir = environment())
  res <- computeFDRwithID(fdrSample$score2, fdrSample$proteinID, larger_better = TRUE)

  cutoff <- predictScoreFDR(res, qValue = 5)

  expect_true(is.finite(cutoff))
  expect_true(all(res$qValue_SimpleFDR[res$score >= cutoff] <= 0.05))
  expect_error(predictScoreFDR(res, method = "unknown"), "no such method")
})
