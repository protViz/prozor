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
      "fdr_simple",
      "q_value_fdr_concatenated",
      "q_value_fdr_simple"
    )
  )
  expect_equal(tab$order, 1:4)
  expect_equal(tab$score, c(4, 3, 2, 1))
  expect_equal(tab$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(tab$fdr_simple, c(0, 0, 1 / 2, 1 / 3))
  expect_equal(tab$q_value_fdr_simple, c(0, 0, 1 / 3, 1 / 3))
  expect_equal(tab$fdr_concatenated, c(0, 0, 2 / 3, 2 / 4))
  expect_equal(tab$q_value_fdr_concatenated, c(0, 0, 1 / 2, 1 / 2))
})

test_that("TargetDecoyFDR ranks small scores first when larger_better = FALSE", {
  fdr <- TargetDecoyFDR$new(score = c(2, 1, 4, 3), decoy_hit = c(FALSE, FALSE, FALSE, TRUE), larger_better = FALSE)

  expect_false(fdr$larger_better)
  expect_equal(fdr$table$order, c(2, 1, 4, 3))
  expect_equal(fdr$table$score, c(1, 2, 3, 4))
  expect_equal(fdr$table$decoy_hit, c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(fdr$table$q_value_fdr_simple, c(0, 0, 1 / 3, 1 / 3))
})

test_that("TargetDecoyFDR ranks decoys first within tied scores, whatever the input order", {
  decoy_first <- TargetDecoyFDR$new(c(2, 2, 1), c(TRUE, FALSE, FALSE))
  decoy_second <- TargetDecoyFDR$new(c(2, 2, 1), c(FALSE, TRUE, FALSE))
  smaller <- TargetDecoyFDR$new(c(1, 1, 2), c(FALSE, TRUE, FALSE), larger_better = FALSE)

  expect_equal(decoy_first$table$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_second$table$decoy_hit, c(TRUE, FALSE, FALSE))
  expect_equal(decoy_first$table$q_value_fdr_simple, c(0.5, 0.5, 0.5))
  expect_equal(decoy_second$table$q_value_fdr_simple, c(0.5, 0.5, 0.5))
  expect_equal(smaller$table$decoy_hit, c(TRUE, FALSE, FALSE))
})

test_that("TargetDecoyFDR computes the FDR per score: tied hits share the counts at the end of the tie", {
  # 100 targets at 10; 10 decoys and 1 target at 5: accepting 5 gives 10 decoys / 101 targets
  fdr <- TargetDecoyFDR$new(c(rep(10, 100), rep(5, 11)), c(rep(FALSE, 100), rep(TRUE, 10), FALSE))

  expect_equal(fdr$score_at(q_value = 0.01), 10)
  expect_equal(unique(fdr$table$fdr_simple[fdr$table$score == 5]), 10 / 101)
})

test_that("TargetDecoyFDR gives one q-value per score on tied bcd scores", {
  decoy <- grepl("^REV_", bcd$protein)
  fdr <- TargetDecoyFDR$new(round(bcd$xcorr, 1), decoy)
  tab <- fdr$table

  n_q_simple <- tapply(tab$q_value_fdr_simple, tab$score, function(q) length(unique(q)))
  n_q_concatenated <- tapply(tab$q_value_fdr_concatenated, tab$score, function(q) length(unique(q)))
  expect_true(all(n_q_simple == 1))
  expect_true(all(n_q_concatenated == 1))
  # per row, the cut-offs were 1.7 (1.07 % realized) and 1.3 (5.04 %)
  expect_equal(fdr$score_at(q_value = 0.01), 1.8)
  expect_equal(fdr$score_at(q_value = 0.05), 1.4)
})

test_that("TargetDecoyFDR on fdr_sample: sorted scores and monotone q-values", {
  decoy <- grepl("^REV_", fdr_sample$protein_id)

  smaller <- TargetDecoyFDR$new(fdr_sample$score, decoy, larger_better = FALSE)$table
  larger <- TargetDecoyFDR$new(fdr_sample$neg_log10_score, decoy, larger_better = TRUE)$table

  expect_false(is.unsorted(smaller$score))
  expect_false(is.unsorted(rev(larger$score)))
  expect_equal(smaller$decoy_hit, decoy[smaller$order])
  expect_equal(larger$decoy_hit, decoy[larger$order])
  expect_false(is.unsorted(smaller$q_value_fdr_simple))
  expect_false(is.unsorted(larger$q_value_fdr_simple))
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

test_that("TargetDecoyFDR stops for zero hits", {
  expect_error(TargetDecoyFDR$new(numeric(), logical()), "`score` and `decoy_hit` must not be empty\\.")
})

test_that("TargetDecoyFDR stops for scores that are not numeric", {
  # sorted as text, the scores would be ranked 9, 100, 10
  expect_error(
    TargetDecoyFDR$new(c("10", "9", "100"), c(FALSE, TRUE, FALSE)),
    "`score` must be numeric, not character\\."
  )
  expect_error(TargetDecoyFDR$new(factor(c(10, 9)), c(FALSE, TRUE)), "`score` must be numeric, not factor\\.")
})

test_that("TargetDecoyFDR stops unless larger_better is TRUE or FALSE", {
  score <- c(4, 3, 2, 1)
  decoy <- c(FALSE, FALSE, TRUE, FALSE)

  expect_error(TargetDecoyFDR$new(score, decoy, larger_better = NA), "`larger_better` must be TRUE or FALSE, not NA\\.")
  expect_error(
    TargetDecoyFDR$new(score, decoy, larger_better = "no"),
    "`larger_better` must be TRUE or FALSE, not \"no\"\\."
  )
  expect_error(
    TargetDecoyFDR$new(score, decoy, larger_better = c(TRUE, FALSE)),
    "`larger_better` must be TRUE or FALSE, not c\\(TRUE, FALSE\\)\\."
  )
})

test_that("TargetDecoyFDR caps the FDR estimates and q-values at 1 when the best hit is a decoy", {
  fdr <- TargetDecoyFDR$new(c(3, 2, 1), c(TRUE, FALSE, TRUE))
  tab <- fdr$table

  # uncapped: fdr_simple Inf, 1, 2 and fdr_concatenated 2, 1, 4/3, with q-values up to 2
  expect_equal(tab$fdr_simple, c(1, 1, 1))
  expect_equal(tab$fdr_concatenated, c(1, 1, 1))
  expect_equal(tab$q_value_fdr_simple, c(1, 1, 1))
  expect_equal(tab$q_value_fdr_concatenated, c(1, 1, 1))
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
  fdr <- TargetDecoyFDR$new(
    fdr_sample$neg_log10_score,
    decoy_hit = grepl("^REV_", fdr_sample$protein_id),
    id = fdr_sample$protein_id
  )

  expect_identical(fdr$table$id, fdr_sample$protein_id[fdr$table$order])
  expect_identical(rownames(TargetDecoyFDR$new(c(a = 1, b = 2), c(FALSE, FALSE))$table), c("1", "2"))
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

  # fdr_simple q-values 0, 0, 1/3, 1/3; fdr_concatenated q-values 0, 0, 1/2, 1/2
  expect_equal(fdr$score_at(q_value = 0.2), 3)
  expect_equal(fdr$score_at(q_value = 0.4), 1)
  expect_equal(fdr$score_at(q_value = 0.4, method = "fdr_simple"), 1)
  expect_equal(fdr$score_at(q_value = 0.4, method = "fdr_concatenated"), 3)
  expect_equal(smaller$score_at(q_value = 0.2), 2)
  # q-values equal to `q_value` are accepted; 0 and 1 are valid
  expect_equal(fdr$score_at(q_value = 1 / 3), 1)
  expect_equal(fdr$score_at(q_value = 0), 3)
  expect_equal(fdr$score_at(q_value = 1), 1)
})

test_that("TargetDecoyFDR$score_at() takes q_value as a fraction, like the q-values in table", {
  fdr <- TargetDecoyFDR$new(fdr_sample$neg_log10_score, decoy_hit = grepl("^REV_", fdr_sample$protein_id))

  cutoff <- fdr$score_at(q_value = 0.01)

  # 1 % FDR; read as percent, 0.01 gave a 0.01 % cut-off with 18339 hits
  expect_equal(sum(fdr$table$score >= cutoff), 34876)
})

test_that("TargetDecoyFDR$score_at() stops for a q_value that is not a fraction", {
  fdr <- TargetDecoyFDR$new(fdr_sample$neg_log10_score, decoy_hit = grepl("^REV_", fdr_sample$protein_id))

  expect_error(fdr$score_at(q_value = 5), "`q_value` must be a fraction between 0 and 1 .*, not 5\\.")
  expect_error(fdr$score_at(q_value = NA), "`q_value` must be a fraction between 0 and 1 .*, not NA\\.")
})

test_that("TargetDecoyFDR$score_at() errors when no score reaches the FDR", {
  fdr <- TargetDecoyFDR$new(c(2, 1), c(TRUE, FALSE))

  expect_error(fdr$score_at(q_value = 0.01), "No score reaches a q-value of 0.01 \\(method \"fdr_simple\"\\)")
})

test_that("TargetDecoyFDR$score_at() rejects unknown methods", {
  fdr <- TargetDecoyFDR$new(c(2, 1), c(FALSE, FALSE))

  expect_error(
    fdr$score_at(method = "unknown"),
    "`method` must be one of \"fdr_simple\", \"fdr_concatenated\", not \"unknown\""
  )
})

test_that("TargetDecoyFDR$print() shows hits, decoys, larger_better and the scores at 1 % and 5 %", {
  # scores 40 to 1: 33 targets, 1 decoy, 1 target, 4 decoys, 1 target; fdr_simple q-values 0 for scores 40 to 8,
  # 1/34 at 7 and 6, then 2/34 and more; fdr_concatenated q-values 2/35 at 7 and 6
  fdr <- TargetDecoyFDR$new(40:1, rep(c(FALSE, TRUE, FALSE, TRUE, FALSE), c(33, 1, 1, 4, 1)))
  none <- TargetDecoyFDR$new(c(1, 2), c(TRUE, FALSE), larger_better = FALSE)

  expect_output(
    expect_invisible(print(fdr)),
    "<TargetDecoyFDR> hits: 40, decoys: 5, larger_better: TRUE\nscore at fdr_simple q-value 0.01: 8, at 0.05: 6",
    fixed = TRUE
  )
  expect_equal(fdr$score_at(q_value = 0.05, method = "fdr_concatenated"), 8)
  expect_equal(fdr$score_at(q_value = 0.2), 1)
  expect_output(
    print(none),
    "<TargetDecoyFDR> hits: 2, decoys: 1, larger_better: FALSE\nscore at fdr_simple q-value 0.01: none, at 0.05: none",
    fixed = TRUE
  )
})

test_that("TargetDecoyFDR$plot() draws the score histogram and the q-value curves in two panels", {
  fdr <- TargetDecoyFDR$new(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))

  p <- fdr$plot(bins = 4)

  expect_s3_class(p, "ggplot")
  hits <- ggplot2::layer_data(p, 1)
  curves <- ggplot2::layer_data(p, 2)
  expect_equal(sum(hits$count), 4)
  expect_true(all(hits$PANEL == 1))
  # targets blue, decoys orange
  expect_setequal(unique(hits$fill[hits$count > 0]), c("#0072B2", "#D55E00"))
  expect_true(all(curves$PANEL == 2))
  # one curve per estimate: fdr_concatenated black, fdr_simple orange
  expect_setequal(unique(curves$colour), c("#000000", "#E69F00"))
  # q-values in percent: fdr_simple 0, 0, 1/3, 1/3 and fdr_concatenated 0, 0, 1/2, 1/2
  expect_setequal(round(curves$y, 2), c(0, 33.33, 50))
})

test_that("TargetDecoyFDR$plot() uses a log score axis when smaller scores are better", {
  scores <- c(1e-5, 1e-3, 1e-2, 1e-1)
  fdr <- TargetDecoyFDR$new(scores, c(FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  p <- fdr$plot()

  expect_equal(sort(unique(ggplot2::layer_data(p, 2)$x)), log10(scores))
})

test_that("TargetDecoyFDR$plot() keeps zero scores, on a linear axis, when smaller scores are better", {
  scores <- c(0, 0, 1e-3, 1e-2, 1e-1)
  fdr <- TargetDecoyFDR$new(scores, c(FALSE, FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  p <- fdr$plot(bins = 5)

  expect_no_warning(hits <- ggplot2::layer_data(p, 1))
  expect_equal(sum(hits$count), 5)
  expect_equal(sort(unique(ggplot2::layer_data(p, 2)$x)), unique(scores))
})
