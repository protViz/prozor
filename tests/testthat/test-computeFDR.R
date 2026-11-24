test_that("computeFDR() matches hand-computed FDR and q-values", {
  res <- computeFDR(score = c(4, 3, 2, 1), decoy_hit = c(FALSE, FALSE, TRUE, FALSE), larger_better = TRUE)

  expect_true(res$larger_better)
  expect_equal(res$order, 1:4)
  expect_equal(res$score, c(4, 3, 2, 1))
  expect_equal(res$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(res$SimpleFDR, c(0, 0, 1 / 2, 1 / 3))
  expect_equal(res$qValue_SimpleFDR, c(0, 0, 1 / 3, 1 / 3))
  expect_equal(res$FPR, c(0, 0, 2 / 3, 2 / 4))
  expect_equal(res$qValue_FPR, c(0, 0, 1 / 2, 1 / 2))
})

test_that("computeFDR() ranks small scores first when larger_better = FALSE", {
  res <- computeFDR(score = c(2, 1, 4, 3), decoy_hit = c(FALSE, FALSE, FALSE, TRUE), larger_better = FALSE)

  expect_false(res$larger_better)
  expect_equal(res$order, c(2, 1, 4, 3))
  expect_equal(res$score, c(1, 2, 3, 4))
  expect_equal(res$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(res$qValue_SimpleFDR, c(0, 0, 1 / 3, 1 / 3))
})

test_that("computeFDR() ranks decoys first within tied scores, whatever the input order", {
  decoy_first <- computeFDR(c(2, 2, 1), c(TRUE, FALSE, FALSE))
  decoy_second <- computeFDR(c(2, 2, 1), c(FALSE, TRUE, FALSE))
  smaller <- computeFDR(c(1, 1, 2), c(FALSE, TRUE, FALSE), larger_better = FALSE)

  expect_equal(decoy_first$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_second$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_first$qValue_SimpleFDR, c(0.5, 0.5, 0.5))
  expect_equal(decoy_second$qValue_SimpleFDR, c(0.5, 0.5, 0.5))
  expect_equal(smaller$decoy_hit, c(TRUE, FALSE, FALSE))
})

test_that("computeFDR() on fdrSample: sorted scores and monotone q-values", {
  data(fdrSample, package = "prozor", envir = environment())
  decoy <- grepl("^REV_", fdrSample$proteinID)

  smaller <- computeFDR(fdrSample$score, decoy, larger_better = FALSE)
  larger <- computeFDR(fdrSample$score2, decoy, larger_better = TRUE)

  expect_false(is.unsorted(smaller$score))
  expect_false(is.unsorted(rev(larger$score)))
  expect_equal(smaller$decoy_hit, decoy[smaller$order])
  expect_equal(larger$decoy_hit, decoy[larger$order])
  expect_false(is.unsorted(smaller$qValue_SimpleFDR))
  expect_false(is.unsorted(larger$qValue_SimpleFDR))
  expect_true(all(smaller$qValue_FPR <= smaller$FPR))
  expect_true(all(larger$qValue_FPR <= larger$FPR))
})

test_that("computeFDRwithID() reorders IDs with the scores and uses the decoy pattern", {
  res <- computeFDRwithID(c(1, 3, 2), c("a", "DEC_b", "c"), decoy = "DEC_")

  expect_equal(res$score, c(3, 2, 1))
  expect_equal(res$ID, c("DEC_b", "c", "a"))
  expect_equal(res$decoy_hit, c(TRUE, FALSE, FALSE))
})

test_that("computeFDRwithID() counts only IDs starting with REV_ as decoys by default", {
  res <- computeFDRwithID(c(2, 1), c("REV_sp|P1|X_HUMAN", "sp|P2|PREV_HUMAN"))

  expect_equal(res$decoy_hit, c(TRUE, FALSE))
})

test_that("predictScoreFDR() returns the least strict score within the FDR", {
  res <- computeFDR(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))
  smaller <- computeFDR(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  # SimpleFDR q-values 0, 0, 1/3, 1/3; FPR q-values 0, 0, 1/2, 1/2
  expect_equal(predictScoreFDR(res, qValue = 20), 3)
  expect_equal(predictScoreFDR(res, qValue = 40), 1)
  expect_equal(predictScoreFDR(res, qValue = 40, method = "FPR"), 3)
  expect_equal(predictScoreFDR(smaller, qValue = 20), 2)
})

test_that("predictScoreFDR() errors when no score reaches the FDR", {
  res <- computeFDR(c(2, 1), c(TRUE, FALSE))

  expect_error(predictScoreFDR(res, qValue = 1), "no score reaches a q-value of 1 %")
})

test_that("predictScoreFDR() rejects unknown methods", {
  res <- computeFDR(c(2, 1), c(FALSE, FALSE))

  expect_error(predictScoreFDR(res, method = "unknown"), "no such method")
})

test_that("plotFDR() restores the graphics parameters", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  graphics::par(mar = c(1, 2, 3, 4))
  res <- computeFDR(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))

  plotFDR(res)

  expect_equal(graphics::par("mar"), c(1, 2, 3, 4))
})
