test_that("Cdsw$new(method = \"constant\") gives n_bins equal-width windows over the mass range", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "constant")

  expect_identical(cdsw$method, "constant")
  expect_length(cdsw$breaks, 11)
  expect_lt(cdsw$breaks[1], min(masses))
  expect_gt(cdsw$breaks[11], max(masses))
  expect_equal(diff(cdsw$breaks), rep(diff(cdsw$breaks)[1], 10), tolerance = 0.01)
})

test_that("Cdsw$new(method = \"constant\") encloses the masses for any digits", {
  min_masses <- seq(400.001, 400.099, by = 0.001)

  first_breaks <- vapply(
    min_masses,
    function(min_mass) Cdsw$new(c(min_mass, 450, 499.99), n_bins = 2, method = "constant", digits = 3)$breaks[1],
    numeric(1)
  )

  expect_true(all(first_breaks < min_masses))
  cdsw <- Cdsw$new(c(400.008, 450, 499.99), n_bins = 2, method = "constant", digits = 3)
  expect_equal(sum(ggplot2::layer_data(cdsw$plot(), 1)$ymax), 3)
  expect_equal(sum(cdsw$error()$normalized_counts > 0), 2)
})

test_that("Cdsw$new() rounds the breaks to digits", {
  for (method in c("mixture", "quantile", "constant")) {
    breaks <- Cdsw$new(masses, method = method, digits = 1)$breaks
    expect_equal(breaks, round(breaks, 1), label = paste(method, "breaks"))
  }
})

test_that("Cdsw$new() uses the mixture method by default and keeps window widths within the limits", {
  cdsw <- Cdsw$new(masses, max_window = 100, min_window = 5)
  breaks <- cdsw$breaks

  expect_identical(cdsw$method, "mixture")
  expect_length(breaks, 26)
  # the limits hold for the rounded and padded breaks
  expect_true(all(diff(breaks) >= 5))
  expect_true(all(diff(breaks) <= 100))
  # the smallest share of uniform data that meets the limits: the widest window is at the limit, up to the rounding
  expect_equal(max(diff(breaks)), 100, tolerance = 0.01 / 100)
})

test_that("Cdsw$new(method = \"mixture\") keeps the precursor counts per window close to equal", {
  cdsw <- Cdsw$new(masses, max_window = 100)

  counts <- cdsw$as_table(overlap = 0)$n_precursors
  expect_lt(stats::sd(counts) / mean(counts), 0.2)
})

test_that("Cdsw$new(method = \"mixture\") is reproducible and does not use the random number generator", {
  withr::local_seed(1)
  seed <- .Random.seed

  first <- Cdsw$new(masses, max_window = 100)$breaks

  expect_identical(.Random.seed, seed)
  expect_identical(Cdsw$new(masses, max_window = 100)$breaks, first)
})

test_that("Cdsw$new(method = \"mixture\") meets limits close to the equal width", {
  equal_width <- diff(range(masses)) / 25

  breaks <- Cdsw$new(masses, max_window = 1.02 * equal_width, min_window = 0.98 * equal_width)$breaks

  expect_true(all(diff(breaks) >= 0.98 * equal_width))
  expect_true(all(diff(breaks) <= 1.02 * equal_width))
})

test_that("Cdsw$new(method = \"mixture\") gives the quantile breaks when they meet the limits", {
  breaks <- Cdsw$new(masses, max_window = 1000, min_window = 1)$breaks

  expect_identical(breaks, Cdsw$new(masses, method = "quantile")$breaks)
})

test_that("Cdsw$new() pads the outer breaks so that the windows enclose the masses, for every method", {
  # unpadded, the first quantile break would be 400.01, above the smallest mass
  ms1 <- c(400.006, 450, 499.994)

  for (method in c("mixture", "quantile", "constant")) {
    cdsw <- Cdsw$new(ms1, n_bins = 2, method = method)
    expect_lt(cdsw$breaks[1], min(ms1), label = paste(method, "first break"))
    expect_gt(cdsw$breaks[3], max(ms1), label = paste(method, "last break"))
    expect_equal(sum(cdsw$as_table(overlap = 0)$n_precursors), 3, label = paste(method, "precursors in the windows"))
  }
})

test_that(".mixture_quantiles() caps a break at the next mass", {
  # without the cap, the middle break of this mixture with 1 % uniform data would be 170
  expect_equal(.mixture_quantiles(c(0, 1, 10), 2, 0.01), c(0, 1, 10))
})

test_that("Cdsw$new(method = \"mixture\") errors when equal-width windows violate the limits", {
  # 25 equal windows over 350-1200 m/z are 34 wide
  expect_error(Cdsw$new(masses, max_window = 10), "No window layout satisfies the limits")
  expect_no_error(Cdsw$new(masses, max_window = 10, method = "quantile"))
})

test_that("Cdsw$new() stops for invalid arguments, naming the argument and the value", {
  expect_error(Cdsw$new(masses, method = "linear"), "`method` must be one of \"mixture\", .*, not \"linear\"")
  expect_error(Cdsw$new(as.character(masses)), "`masses` must be a non-empty numeric vector, not character")
  expect_error(Cdsw$new(numeric()), "`masses` must be a non-empty numeric vector")
  expect_error(Cdsw$new(c(400, NA, 500)), "`masses` must not contain NA, NaN or infinite values; found 1")
  expect_error(Cdsw$new(masses, n_bins = 0), "`n_bins` must be a single positive whole number, not 0")
  expect_error(Cdsw$new(masses, n_bins = 2.5), "`n_bins` must be a single positive whole number, not 2.5")
  expect_error(Cdsw$new(masses, max_window = NA), "`max_window` must be a single positive number, not NA")
  expect_error(Cdsw$new(masses, min_window = -1), "`min_window` must be a single non-negative number, not -1")
  expect_error(Cdsw$new(masses, digits = Inf), "`digits` must be a single non-negative whole number, not Inf")
  expect_error(
    Cdsw$new(masses, digits = c(1, 2)),
    "`digits` must be a single non-negative whole number, not c\\(1, 2\\)"
  )
})

test_that("Cdsw fields are read-only", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "quantile")
  breaks <- cdsw$breaks

  expect_error(cdsw$breaks <- seq(400, 1400, by = 100), "`breaks` of a Cdsw object is read-only")
  expect_error(cdsw$masses <- 500, "`masses` of a Cdsw object is read-only")
  expect_error(cdsw$n_bins <- 5, "`n_bins` of a Cdsw object is read-only")
  expect_error(cdsw$method <- "constant", "`method` of a Cdsw object is read-only")
  expect_identical(cdsw$breaks, breaks)
  expect_identical(cdsw$n_bins, 10)
  expect_identical(cdsw$masses, masses)
})

test_that("Cdsw results do not depend on the order of the method calls", {
  first <- Cdsw$new(masses, max_window = 100)
  second <- Cdsw$new(masses, max_window = 100)

  table_first <- first$as_table()
  error_first <- first$error()
  optimized_first <- first$optimize_windows(min_count = 10)
  first$plot_breaks()
  first$plot_cycle()

  second$plot_breaks()
  optimized_second <- second$optimize_windows(min_count = 10)
  second$plot_cycle(overlap = 3)
  error_second <- second$error()
  table_second <- second$as_table()

  expect_identical(table_second, table_first)
  expect_identical(error_second, error_first)
  expect_identical(optimized_second, optimized_first)
  expect_identical(first$breaks, second$breaks)
})

test_that("Cdsw$clone() returns an independent object with the same state", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "quantile")

  cloned <- cdsw$clone()

  expect_identical(cloned$breaks, cdsw$breaks)
  expect_identical(cloned$masses, cdsw$masses)
  expect_identical(cloned$n_bins, cdsw$n_bins)
  expect_identical(cloned$method, cdsw$method)
  expect_identical(cloned$as_table(), cdsw$as_table())
  # the fields of the clone are read-only too
  expect_error(cloned$breaks <- seq(300, 1300, by = 100), "`breaks` of a Cdsw object is read-only")
  expect_identical(cloned$breaks, cdsw$breaks)
})

test_that("Cdsw$print() shows the masses, windows, method and window widths", {
  cdsw <- Cdsw$new(c(405, 410, 420, 430, 435), n_bins = 2, method = "constant", digits = 0)

  expect_output(
    expect_invisible(print(cdsw)),
    "<Cdsw> 5 masses in 2 windows, method \"constant\"\nwindow widths 16 to 16 m/z"
  )
})

test_that("Cdsw$as_table() returns one row per window", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "constant")

  tab <- cdsw$as_table(overlap = 1)

  expect_named(tab, c("from", "to", "mid", "width", "n_precursors"))
  expect_equal(nrow(tab), 10)
  expect_equal(tab$from, cdsw$breaks[1:10] - 0.5)
  expect_equal(tab$to, cdsw$breaks[2:11] + 0.5)
  expect_equal(tab$n_precursors[1], sum(masses > tab$from[1] & masses <= tab$to[1]))
})

test_that("Cdsw$as_table(), $optimize_windows() and $plot_cycle() stop for an invalid overlap", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "constant")

  expect_error(cdsw$as_table(overlap = "1"), "`overlap` must be a single non-negative number, not \"1\"")
  expect_error(cdsw$optimize_windows(overlap = -1), "`overlap` must be a single non-negative number, not -1")
  expect_error(cdsw$plot_cycle(overlap = NA), "`overlap` must be a single non-negative number, not NA")
})

test_that("Cdsw$as_table(), $plot() and $error() count the masses in (from, to]", {
  # digits = 0 gives the breaks 404, 420 and 436, so the mass 420 lies on the middle break
  cdsw <- Cdsw$new(c(405, 410, 420, 430, 435), n_bins = 2, method = "constant", digits = 0)

  expect_equal(cdsw$breaks, c(404, 420, 436))
  expect_equal(cdsw$as_table(overlap = 0)$n_precursors, c(3, 2))
  expect_equal(cdsw$as_table(overlap = 20)$n_precursors, c(4, 3))
  expect_equal(ggplot2::layer_data(cdsw$plot(), 1)$ymax, c(3, 2))
  expect_equal(cdsw$error()$normalized_counts, c(3, 2) / sqrt(13))
})

test_that("Cdsw$error() gives both distances and the normalized observed and equal counts", {
  cdsw <- Cdsw$new(c(405, 410, 420, 430, 435), n_bins = 2, method = "constant", digits = 0)

  err <- cdsw$error()

  expect_named(err, c("manhattan", "euclidean", "normalized_counts", "equal_counts"))
  # counts 3 and 2, normalized to length 1
  observed <- c(3, 2) / sqrt(13)
  expect_equal(err$equal_counts, c(1, 1) / sqrt(2))
  expect_equal(err$manhattan, mean(1 - observed / max(observed)))
  expect_equal(err$euclidean, sqrt(sum((observed - err$equal_counts)^2)) / 2)
})

test_that("Cdsw$plot() draws one bar per window with the precursor counts and a line at equal counts", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "constant")

  p <- cdsw$plot()

  expect_s3_class(p, "ggplot")
  bars <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(bars), 10)
  expect_equal(bars$xmin, cdsw$breaks[1:10])
  expect_equal(sum(bars$ymax), length(masses))
  expect_equal(ggplot2::layer_data(p, 2)$yintercept, length(masses) / 10)
})

test_that("Cdsw$plot_breaks() compares the breaks with constant and quantile breaks", {
  cdsw <- Cdsw$new(masses, n_bins = 10, max_window = 100)

  p <- cdsw$plot_breaks()

  expect_s3_class(p, "ggplot")
  points <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(points), 3 * 11)
  expect_length(unique(points$colour), 3)
  expect_equal(
    points$x,
    c(
      Cdsw$new(masses, n_bins = 10, method = "constant")$breaks,
      Cdsw$new(masses, n_bins = 10, method = "quantile")$breaks,
      cdsw$breaks
    )
  )
})

test_that("Cdsw$plot_cycle() draws one rectangle per window", {
  cdsw <- Cdsw$new(masses, n_bins = 10, method = "constant")

  p <- cdsw$plot_cycle(overlap = 1)

  expect_s3_class(p, "ggplot")
  rects <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(rects), 10)
  expect_equal(rects$ymin, cdsw$as_table(overlap = 1)$from)
  # alternating fills, starting with the second colour at step 1
  expect_equal(rects$fill, rep(c("#E69F00", "#56B4E9"), 5))
})

test_that("Cdsw$plot_cycle() passes overlap on to as_table()", {
  cdsw <- Cdsw$new(masses, max_window = 100)

  rects <- ggplot2::layer_data(cdsw$plot_cycle(overlap = 3), 1)

  expect_equal(rects$ymin, cdsw$as_table(overlap = 3)$from)
})

test_that("Cdsw$optimize_windows() passes overlap, bin_width and min_count on", {
  cdsw <- Cdsw$new(masses, max_window = 100)

  res <- cdsw$optimize_windows(overlap = 3, bin_width = 0.5, min_count = 10)

  expect_identical(res, readjust_windows(cdsw$as_table(overlap = 3), masses, bin_width = 0.5, min_count = 10))
})

test_that("readjust_windows() moves boundaries but keeps one window per row", {
  cdsw <- Cdsw$new(masses, max_window = 100)
  tab <- cdsw$as_table()

  res <- readjust_windows(tab, masses, min_count = 10)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("from", "to", "mid", "width", "n_precursors"))
  expect_equal(nrow(res), nrow(tab))
  expect_true(all(res$to > res$from))
  expect_true(all(abs(res$from - tab$from) < 1))
  expect_equal(res$n_precursors[1], sum(masses > res$from[1] & masses <= res$to[1]))
  expect_identical(cdsw$optimize_windows(min_count = 10), res)
})

test_that("readjust_windows() moves starts down and ends up to the nearest empty bin, also the first and last", {
  # 0.1 m/z bins with 10 masses each, except the empty bins at 399.95, 450.05 and 500.05
  ms1 <- seq(400.005, 500, by = 0.01)
  ms1 <- ms1[ms1 <= 450 | ms1 > 450.1]
  windows <- data.frame(from = c(400, 450.1), to = c(450, 500))

  res <- readjust_windows(windows, ms1, min_count = 5)

  expect_equal(res$from, c(399.95, 450.05))
  expect_equal(res$to, c(450.05, 500.05))
  expect_equal(res$n_precursors, c(sum(ms1 <= 450), sum(ms1 > 450.1)))
})

test_that("readjust_windows() keeps a boundary that is already on an empty-bin mid", {
  # bins of 1 m/z from 399 to 404: the empty bins have the mids 399.5, 401.5 and 403.5
  windows <- data.frame(from = c(399.5, 401.5), to = c(401.5, 403.5))

  res <- readjust_windows(windows, c(400.2, 400.7, 402.5), bin_width = 1, min_count = 1)

  expect_equal(res$from, c(399.5, 401.5))
  expect_equal(res$to, c(401.5, 403.5))
  expect_equal(res$n_precursors, c(2, 1))
})

test_that("readjust_windows() returns its own result unchanged", {
  once <- readjust_windows(Cdsw$new(masses, max_window = 100)$as_table(), masses, min_count = 10)

  expect_identical(readjust_windows(once, masses, min_count = 10), once)
})

test_that("readjust_windows() keeps the windows when no bin is empty", {
  ms1 <- c(405, 410, 420, 430, 435)
  windows <- data.frame(from = c(400, 420), to = c(420, 440))

  res <- readjust_windows(windows, ms1, min_count = 0)

  expect_equal(res$from, windows$from)
  expect_equal(res$to, windows$to)
  expect_equal(res$n_precursors, c(3, 2))
})

test_that("readjust_windows() bins the masses in bins of 0.01 m/z", {
  ms1 <- c(400.12, 512.37, 700.5, 1000.04)
  windows <- data.frame(from = c(400, 700), to = c(700, 1000.1))

  res <- readjust_windows(windows, ms1, bin_width = 0.01, min_count = 1)

  expect_equal(res$from, c(400, 699.995))
  expect_equal(res$to, c(700.005, 1000.1))
  expect_equal(res$n_precursors, c(2, 2))
})

test_that("readjust_windows() runs wherever the smallest mass lies in its bin of 0.01 m/z", {
  expect_no_error(
    for (min_mass in seq(400.001, 400.099, by = 0.001)) {
      readjust_windows(list(from = 500, to = 600), c(min_mass, 500, 700), bin_width = 0.01, min_count = 1)
    }
  )
})

test_that("readjust_windows() moves the boundaries to empty bins of 1 m/z", {
  # bins of 1 m/z from 399 to 404: the empty bins have the mids 399.5, 401.5 and 403.5
  res <- readjust_windows(list(from = 401, to = 402), c(400.2, 400.7, 402.5), bin_width = 1, min_count = 1)

  expect_equal(c(res$from, res$to, res$n_precursors), c(399.5, 403.5, 3))
})

test_that("readjust_windows() stops for a bin_width that is not a positive number", {
  windows <- data.frame(from = c(400, 700), to = c(700, 1000.1))

  expect_error(
    readjust_windows(windows, c(400.12, 512.37), bin_width = 0),
    "`bin_width` must be a single positive number, not 0"
  )
})

test_that("readjust_windows() stops for invalid masses and windows without `from` or `to`", {
  windows <- data.frame(from = c(400, 450), to = c(450, 500))

  expect_error(
    readjust_windows(windows, c(420, NA)),
    "`masses` must not contain NA, NaN or infinite values; found 1\\."
  )
  expect_error(readjust_windows(windows, "420"), "`masses` must be a non-empty numeric vector, not character")
  expect_error(
    readjust_windows(windows["from"], 420),
    "`windows` must have the columns `from`, `to`; missing: `to`\\."
  )
  expect_error(readjust_windows(c(400, 450), 420), "`windows` must be a data frame or list, not numeric\\.")
})

test_that("readjust_windows() stops for a min_count that is not a non-negative number", {
  windows <- data.frame(from = c(400, 450), to = c(450, 500))

  expect_error(
    readjust_windows(windows, 420, min_count = "a"),
    "`min_count` must be a single non-negative number, not \"a\""
  )
  expect_error(
    readjust_windows(windows, 420, min_count = -1),
    "`min_count` must be a single non-negative number, not -1"
  )
})

test_that("plot_windows() draws the masses with a start and an end line per window", {
  windows <- data.frame(from = c(400, 500), to = c(510, 610))

  p <- plot_windows(windows, masses)

  expect_s3_class(p, "ggplot")
  lines <- ggplot2::layer_data(p, 2)
  expect_setequal(lines$xintercept, c(400, 500, 510, 610))
  bars <- ggplot2::layer_data(p, 1)
  expect_equal(sum(bars$count), length(masses))
  expect_equal(bars$xmax - bars$xmin, rep(0.1, nrow(bars)))
  wide <- ggplot2::layer_data(plot_windows(windows, masses, bin_width = 5), 1)
  expect_equal(wide$xmax - wide$xmin, rep(5, nrow(wide)))
})

test_that("plot_windows() stops for invalid masses and windows without `from` or `to`", {
  windows <- data.frame(from = c(400, 500), to = c(510, 610))

  expect_error(plot_windows(windows, c(420, Inf)), "`masses` must not contain NA, NaN or infinite values; found 1\\.")
  expect_error(plot_windows(windows, numeric()), "`masses` must be a non-empty numeric vector, not numeric of length 0")
  expect_error(plot_windows(list(start = 400, to = 510), 420), "missing: `from`\\.")
})
