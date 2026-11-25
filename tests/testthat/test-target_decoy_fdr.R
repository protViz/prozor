test_that("TargetDecoyFDR matches hand-computed FDR and q-values", {
  fdr <- TargetDecoyFDR$new(score = c(4, 3, 2, 1), decoy_hit = c(FALSE, FALSE, TRUE, FALSE), larger_better = TRUE)
  tab <- fdr$table

  expect_true(fdr$larger_better)
  expect_identical(fdr$n, 4L)
  expect_s3_class(tab, "data.frame")
  expect_named(
    tab,
    c(
      "order",
      "score",
      "decoy_hit",
      "fdr_concatenated",
      "simple_fdr",
      "q_value_fdr_concatenated",
      "q_value_simple_fdr"
    )
  )
  expect_equal(tab$order, 1:4)
  expect_equal(tab$score, c(4, 3, 2, 1))
  expect_equal(tab$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(tab$simple_fdr, c(0, 0, 1 / 2, 1 / 3))
  expect_equal(tab$q_value_simple_fdr, c(0, 0, 1 / 3, 1 / 3))
  expect_equal(tab$fdr_concatenated, c(0, 0, 2 / 3, 2 / 4))
  expect_equal(tab$q_value_fdr_concatenated, c(0, 0, 1 / 2, 1 / 2))
})

test_that("TargetDecoyFDR ranks small scores first when larger_better = FALSE", {
  fdr <- TargetDecoyFDR$new(score = c(2, 1, 4, 3), decoy_hit = c(FALSE, FALSE, FALSE, TRUE), larger_better = FALSE)

  expect_false(fdr$larger_better)
  expect_equal(fdr$table$order, c(2, 1, 4, 3))
  expect_equal(fdr$table$score, c(1, 2, 3, 4))
  expect_equal(fdr$table$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(fdr$table$q_value_simple_fdr, c(0, 0, 1 / 3, 1 / 3))
})

test_that("TargetDecoyFDR ranks decoys first within tied scores, whatever the input order", {
  decoy_first <- TargetDecoyFDR$new(c(2, 2, 1), c(TRUE, FALSE, FALSE))
  decoy_second <- TargetDecoyFDR$new(c(2, 2, 1), c(FALSE, TRUE, FALSE))
  smaller <- TargetDecoyFDR$new(c(1, 1, 2), c(FALSE, TRUE, FALSE), larger_better = FALSE)

  expect_equal(decoy_first$table$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_second$table$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_first$table$q_value_simple_fdr, c(0.5, 0.5, 0.5))
  expect_equal(decoy_second$table$q_value_simple_fdr, c(0.5, 0.5, 0.5))
  expect_equal(smaller$table$decoy_hit, c(TRUE, FALSE, FALSE))
})

test_that("TargetDecoyFDR on fdr_sample: sorted scores and monotone q-values", {
  data(fdr_sample, package = "prozor", envir = environment())
  decoy <- grepl("^REV_", fdr_sample$protein_id)

  smaller <- TargetDecoyFDR$new(fdr_sample$score, decoy, larger_better = FALSE)$table
  larger <- TargetDecoyFDR$new(fdr_sample$score2, decoy, larger_better = TRUE)$table

  expect_false(is.unsorted(smaller$score))
  expect_false(is.unsorted(rev(larger$score)))
  expect_equal(smaller$decoy_hit, decoy[smaller$order])
  expect_equal(larger$decoy_hit, decoy[larger$order])
  expect_false(is.unsorted(smaller$q_value_simple_fdr))
  expect_false(is.unsorted(larger$q_value_simple_fdr))
  expect_true(all(smaller$q_value_fdr_concatenated <= smaller$fdr_concatenated))
  expect_true(all(larger$q_value_fdr_concatenated <= larger$fdr_concatenated))
})

test_that("TargetDecoyFDR stops for NA, non-logical decoy flags and unequal lengths", {
  score <- c(4, 3, 2, 1)
  decoy <- c(FALSE, FALSE, TRUE, FALSE)

  expect_error(TargetDecoyFDR$new(score, c(FALSE, TRUE, NA, FALSE)), "`decoy_hit` must not contain NA; found 1\\.")
  expect_error(TargetDecoyFDR$new(c(4, NA, NA, 1), decoy), "`score` must not contain NA; found 2\\.")
  expect_error(
    TargetDecoyFDR$new(score, c("FALSE", "FALSE", "TRUE", "FALSE")),
    "`decoy_hit` must be logical, not character\\."
  )
  expect_error(
    TargetDecoyFDR$new(score, decoy[1:3]),
    "`score` and `decoy_hit` must have the same length, not 4 and 3\\."
  )
  expect_error(
    TargetDecoyFDR$new(score, decoy, id = c("a", "REV_b")),
    "`score` and `id` must have the same length, not 4 and 2\\."
  )
})

test_that("TargetDecoyFDR puts the IDs, sorted with the scores, after `order`", {
  id <- c("a", "DEC_b", "c")
  fdr <- TargetDecoyFDR$new(c(1, 3, 2), decoy_hit = grepl("^DEC_", id), id = id)

  expect_identical(names(fdr$table)[1:3], c("order", "id", "score"))
  expect_equal(fdr$table$score, c(3, 2, 1))
  expect_equal(fdr$table$id, c("DEC_b", "c", "a"))
  expect_equal(fdr$table$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(fdr$table$id, id[fdr$table$order])
  expect_false("id" %in% names(TargetDecoyFDR$new(c(1, 3, 2), decoy_hit = grepl("^DEC_", id))$table))
})

test_that("TargetDecoyFDR keeps factor IDs and has the rank as row names", {
  data(fdr_sample, package = "prozor", envir = environment())
  fdr <- TargetDecoyFDR$new(
    fdr_sample$score2,
    decoy_hit = grepl("^REV_", fdr_sample$protein_id),
    id = fdr_sample$protein_id
  )

  expect_identical(fdr$table$id, fdr_sample$protein_id[fdr$table$order])
  expect_identical(rownames(TargetDecoyFDR$new(c(a = 1, b = 2), c(FALSE, FALSE))$table), c("1", "2"))
})

test_that("the decoy pattern ^REV_ counts only IDs starting with REV_ as decoys", {
  id <- c("REV_sp|P1|X_HUMAN", "sp|P2|PREV_HUMAN")
  fdr <- TargetDecoyFDR$new(c(2, 1), decoy_hit = grepl("^REV_", id), id = id)

  expect_equal(fdr$table$decoy_hit, c(TRUE, FALSE))
})

test_that("TargetDecoyFDR fields are read-only", {
  fdr <- TargetDecoyFDR$new(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))
  tab <- fdr$table

  expect_error(fdr$table <- tab[1:2, ], "`table` of a TargetDecoyFDR object is read-only")
  expect_error(fdr$larger_better <- FALSE, "`larger_better` of a TargetDecoyFDR object is read-only")
  expect_error(fdr$n <- 2, "`n` of a TargetDecoyFDR object is read-only")
  expect_identical(fdr$table, tab)
  expect_true(fdr$larger_better)
  expect_identical(fdr$n, 4L)
})

test_that("TargetDecoyFDR$score_at() returns the least strict score within the FDR", {
  fdr <- TargetDecoyFDR$new(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))
  smaller <- TargetDecoyFDR$new(c(1, 2, 3, 4), c(FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  # simple_fdr q-values 0, 0, 1/3, 1/3; fdr_concatenated q-values 0, 0, 1/2, 1/2
  expect_equal(fdr$score_at(q_value = 0.2), 3)
  expect_equal(fdr$score_at(q_value = 0.4), 1)
  expect_equal(fdr$score_at(q_value = 0.4, method = "fdr_concatenated"), 3)
  expect_equal(smaller$score_at(q_value = 0.2), 2)
})

test_that("TargetDecoyFDR$score_at() takes q_value as a fraction, like the q-values in table", {
  data(fdr_sample, package = "prozor", envir = environment())
  fdr <- TargetDecoyFDR$new(fdr_sample$score2, decoy_hit = grepl("^REV_", fdr_sample$protein_id))

  cutoff <- fdr$score_at(q_value = 0.01)

  # 1 % FDR; read as percent, 0.01 gave a 0.01 % cut-off with 18339 hits
  expect_equal(sum(fdr$table$score >= cutoff), 34876)
  expect_error(fdr$score_at(q_value = 5), "`q_value` must be a fraction between 0 and 1 .*, not 5\\.")
  expect_error(fdr$score_at(q_value = NA), "`q_value` must be a fraction between 0 and 1 .*, not NA\\.")
})

test_that("TargetDecoyFDR$score_at() errors when no score reaches the FDR", {
  fdr <- TargetDecoyFDR$new(c(2, 1), c(TRUE, FALSE))

  expect_error(fdr$score_at(q_value = 0.01), "No score reaches a q-value of 0.01 \\(method \"simple_fdr\"\\)")
})

test_that("TargetDecoyFDR$score_at() rejects unknown methods", {
  fdr <- TargetDecoyFDR$new(c(2, 1), c(FALSE, FALSE))

  expect_error(
    fdr$score_at(method = "unknown"),
    "`method` must be one of \"simple_fdr\", \"fdr_concatenated\"; not \"unknown\""
  )
})

test_that("TargetDecoyFDR$print() shows hits, decoys, larger_better and the scores at 1 % and 5 %", {
  # simple_fdr q-values 0, 0, 0, 0, 1/4: the score at 0.01 and 0.05 is 2
  fdr <- TargetDecoyFDR$new(c(5, 4, 3, 2, 1), c(FALSE, FALSE, FALSE, FALSE, TRUE))
  none <- TargetDecoyFDR$new(c(1, 2), c(TRUE, FALSE), larger_better = FALSE)

  expect_output(
    expect_invisible(print(fdr)),
    "<TargetDecoyFDR> hits: 5, decoys: 1, larger_better: TRUE\nscore at simple_fdr q-value 0.01: 2, at 0.05: 2",
    fixed = TRUE
  )
  expect_output(
    print(none),
    "<TargetDecoyFDR> hits: 2, decoys: 1, larger_better: FALSE\nscore at simple_fdr q-value 0.01: none, at 0.05: none",
    fixed = TRUE
  )
})
