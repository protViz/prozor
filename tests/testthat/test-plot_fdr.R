test_that("plot_fdr() draws the score histogram and the q-value curves in two panels", {
  fdr <- compute_fdr(c(4, 3, 2, 1), c(FALSE, FALSE, TRUE, FALSE))

  p <- plot_fdr(fdr, bins = 4)

  expect_s3_class(p, "ggplot")
  hits <- ggplot2::layer_data(p, 1)
  curves <- ggplot2::layer_data(p, 2)
  expect_equal(sum(hits$count), 4)
  expect_true(all(hits$PANEL == 1))
  expect_setequal(unique(hits$fill[hits$count > 0]), c("#0072B2", "#D55E00"))
  expect_true(all(curves$PANEL == 2))
  # q-values in percent: simple_fdr 0, 0, 1/3, 1/3 and fdr_concatenated 0, 0, 1/2, 1/2
  expect_setequal(round(curves$y, 2), c(0, 33.33, 50))
})

test_that("plot_fdr() uses a log score axis when smaller scores are better", {
  scores <- c(1e-5, 1e-3, 1e-2, 1e-1)
  fdr <- compute_fdr(scores, c(FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  p <- plot_fdr(fdr)

  expect_equal(sort(unique(ggplot2::layer_data(p, 2)$x)), log10(scores))
})

test_that("plot_fdr() keeps zero scores, on a linear axis, when smaller scores are better", {
  scores <- c(0, 0, 1e-3, 1e-2, 1e-1)
  fdr <- compute_fdr(scores, c(FALSE, FALSE, FALSE, TRUE, FALSE), larger_better = FALSE)

  p <- plot_fdr(fdr, bins = 5)

  expect_no_warning(hits <- ggplot2::layer_data(p, 1))
  expect_equal(sum(hits$count), 5)
  expect_equal(sort(unique(ggplot2::layer_data(p, 2)$x)), unique(scores))
})
