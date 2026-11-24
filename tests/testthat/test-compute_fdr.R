test_that("compute_fdr() matches hand-computed FDR and q-values", {
  res <- compute_fdr(score = c(4, 3, 2, 1), decoy_hit = c(FALSE, FALSE, TRUE, FALSE), larger_better = TRUE)

  expect_true(res$larger_better)
  expect_equal(res$order, 1:4)
  expect_equal(res$score, c(4, 3, 2, 1))
  expect_equal(res$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(res$simple_fdr, c(0, 0, 1 / 2, 1 / 3))
  expect_equal(res$q_value_simple_fdr, c(0, 0, 1 / 3, 1 / 3))
  expect_equal(res$fdr_concatenated, c(0, 0, 2 / 3, 2 / 4))
  expect_equal(res$q_value_fdr_concatenated, c(0, 0, 1 / 2, 1 / 2))
})

test_that("compute_fdr() ranks small scores first when larger_better = FALSE", {
  res <- compute_fdr(score = c(2, 1, 4, 3), decoy_hit = c(FALSE, FALSE, FALSE, TRUE), larger_better = FALSE)

  expect_false(res$larger_better)
  expect_equal(res$order, c(2, 1, 4, 3))
  expect_equal(res$score, c(1, 2, 3, 4))
  expect_equal(res$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(res$q_value_simple_fdr, c(0, 0, 1 / 3, 1 / 3))
})

test_that("compute_fdr() ranks decoys first within tied scores, whatever the input order", {
  decoy_first <- compute_fdr(c(2, 2, 1), c(TRUE, FALSE, FALSE))
  decoy_second <- compute_fdr(c(2, 2, 1), c(FALSE, TRUE, FALSE))
  smaller <- compute_fdr(c(1, 1, 2), c(FALSE, TRUE, FALSE), larger_better = FALSE)

  expect_equal(decoy_first$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_second$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_first$q_value_simple_fdr, c(0.5, 0.5, 0.5))
  expect_equal(decoy_second$q_value_simple_fdr, c(0.5, 0.5, 0.5))
  expect_equal(smaller$decoy_hit, c(TRUE, FALSE, FALSE))
})

test_that("compute_fdr() on fdr_sample: sorted scores and monotone q-values", {
  data(fdr_sample, package = "prozor", envir = environment())
  decoy <- grepl("^REV_", fdr_sample$protein_id)

  smaller <- compute_fdr(fdr_sample$score, decoy, larger_better = FALSE)
  larger <- compute_fdr(fdr_sample$score2, decoy, larger_better = TRUE)

  expect_false(is.unsorted(smaller$score))
  expect_false(is.unsorted(rev(larger$score)))
  expect_equal(smaller$decoy_hit, decoy[smaller$order])
  expect_equal(larger$decoy_hit, decoy[larger$order])
  expect_false(is.unsorted(smaller$q_value_simple_fdr))
  expect_false(is.unsorted(larger$q_value_simple_fdr))
  expect_true(all(smaller$q_value_fdr_concatenated <= smaller$fdr_concatenated))
  expect_true(all(larger$q_value_fdr_concatenated <= larger$fdr_concatenated))
})

test_that("compute_fdr() stops for NA, non-logical decoy flags and unequal lengths", {
  score <- c(4, 3, 2, 1)
  decoy <- c(FALSE, FALSE, TRUE, FALSE)

  expect_error(compute_fdr(score, c(FALSE, TRUE, NA, FALSE)), "`decoy_hit` must not contain NA; found 1\\.")
  expect_error(compute_fdr(c(4, NA, NA, 1), decoy), "`score` must not contain NA; found 2\\.")
  expect_error(
    compute_fdr(score, c("FALSE", "FALSE", "TRUE", "FALSE")),
    "`decoy_hit` must be logical, not character\\."
  )
  expect_error(compute_fdr(score, decoy[1:3]), "`score` and `decoy_hit` must have the same length, not 4 and 3\\.")
  expect_error(compute_fdr_with_id(score, c("a", "REV_b")), "must have the same length, not 4 and 2\\.")
})

test_that("compute_fdr_with_id() reorders IDs with the scores and uses the decoy pattern", {
  res <- compute_fdr_with_id(c(1, 3, 2), c("a", "DEC_b", "c"), decoy = "DEC_")

  expect_equal(res$score, c(3, 2, 1))
  expect_equal(res$id, c("DEC_b", "c", "a"))
  expect_equal(res$decoy_hit, c(TRUE, FALSE, FALSE))
})

test_that("compute_fdr_with_id() counts only IDs starting with REV_ as decoys by default", {
  res <- compute_fdr_with_id(c(2, 1), c("REV_sp|P1|X_HUMAN", "sp|P2|PREV_HUMAN"))

  expect_equal(res$decoy_hit, c(TRUE, FALSE))
})

test_that("predict_score_fdr() returns the least strict score within the FDR", {
  res <- compute_fdr(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))
  smaller <- compute_fdr(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  # simple_fdr q-values 0, 0, 1/3, 1/3; fdr_concatenated q-values 0, 0, 1/2, 1/2
  expect_equal(predict_score_fdr(res, q_value = 0.2), 3)
  expect_equal(predict_score_fdr(res, q_value = 0.4), 1)
  expect_equal(predict_score_fdr(res, q_value = 0.4, method = "fdr_concatenated"), 3)
  expect_equal(predict_score_fdr(smaller, q_value = 0.2), 2)
})

test_that("predict_score_fdr() takes q_value as a fraction, like the q-values of compute_fdr()", {
  data(fdr_sample, package = "prozor", envir = environment())
  fdr <- compute_fdr_with_id(fdr_sample$score2, fdr_sample$protein_id, larger_better = TRUE)

  cutoff <- predict_score_fdr(fdr, q_value = 0.01)

  # 1 % FDR; read as percent, 0.01 gave a 0.01 % cut-off with 18339 hits
  expect_equal(sum(fdr$score >= cutoff), 34876)
  expect_error(predict_score_fdr(fdr, q_value = 5), "`q_value` must be a fraction between 0 and 1 .*, not 5\\.")
})

test_that("predict_score_fdr() errors when no score reaches the FDR", {
  res <- compute_fdr(c(2, 1), c(TRUE, FALSE))

  expect_error(predict_score_fdr(res, q_value = 0.01), "No score reaches a q-value of 0.01")
})

test_that("predict_score_fdr() rejects unknown methods", {
  res <- compute_fdr(c(2, 1), c(FALSE, FALSE))

  expect_error(
    predict_score_fdr(res, method = "unknown"),
    "`method` must be \"simple_fdr\" or \"fdr_concatenated\", not \"unknown\""
  )
})
