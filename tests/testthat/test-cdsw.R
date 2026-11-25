test_that("Cdsw$new(method = \"constant\") gives nbins equal-width windows over the mass range", {
  data(masses, package = "prozor", envir = environment())

  cdsw <- Cdsw$new(masses, nbins = 10, method = "constant")

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
    function(min_mass) Cdsw$new(c(min_mass, 450, 499.99), nbins = 2, method = "constant", digits = 3)$breaks[1],
    numeric(1)
  )

  expect_true(all(first_breaks < min_masses))
  cdsw <- Cdsw$new(c(400.008, 450, 499.99), nbins = 2, method = "constant", digits = 3)
  expect_equal(sum(ggplot2::layer_data(cdsw$plot(), 1)$ymax), 3)
  expect_equal(sum(cdsw$error()$counts > 0), 2)
})

test_that("Cdsw$new() rounds the breaks to digits", {
  data(masses, package = "prozor", envir = environment())

  for (method in c("sampling", "quantile", "constant")) {
    breaks <- Cdsw$new(masses, method = method, digits = 1)$breaks
    expect_equal(breaks, round(breaks, 1))
  }
})

test_that("Cdsw$new() samples by default and keeps window widths within the limits", {
  data(masses, package = "prozor", envir = environment())

  cdsw <- Cdsw$new(masses, maxwindow = 100, minwindow = 5)
  breaks <- cdsw$breaks

  expect_identical(cdsw$method, "sampling")
  expect_length(breaks, 26)
  # rounding to digits = 2 moves a break by up to 0.005; the outer breaks are padded by up to 0.02
  expect_true(all(diff(breaks) >= 5))
  expect_true(all(diff(breaks) <= 100 + 0.025))
  # the smallest share of uniform data that meets the limits: the widest window is at the limit
  expect_equal(max(diff(breaks)), 100, tolerance = 0.025 / 100)
})

test_that("Cdsw$new(method = \"sampling\") keeps the precursor counts per window close to equal", {
  data(masses, package = "prozor", envir = environment())

  cdsw <- Cdsw$new(masses, maxwindow = 100)

  counts <- cdsw$as_table(overlap = 0)$counts
  expect_lt(stats::sd(counts) / mean(counts), 0.2)
})

test_that("Cdsw$new(method = \"sampling\") is reproducible and does not use the random number generator", {
  data(masses, package = "prozor", envir = environment())
  withr::local_seed(1)
  seed <- .Random.seed

  first <- Cdsw$new(masses, maxwindow = 100)$breaks

  expect_identical(.Random.seed, seed)
  expect_identical(Cdsw$new(masses, maxwindow = 100)$breaks, first)
})

test_that("Cdsw$new(method = \"sampling\") meets limits close to the equal width", {
  data(masses, package = "prozor", envir = environment())
  equal_width <- diff(range(masses)) / 25

  breaks <- Cdsw$new(masses, maxwindow = 1.02 * equal_width, minwindow = 0.98 * equal_width)$breaks

  expect_true(all(diff(breaks) >= 0.98 * equal_width - 0.01))
  expect_true(all(diff(breaks) <= 1.02 * equal_width + 0.025))
})

test_that("Cdsw$new(method = \"sampling\") gives the quantile breaks when they meet the limits", {
  data(masses, package = "prozor", envir = environment())

  breaks <- Cdsw$new(masses, maxwindow = 1000, minwindow = 1)$breaks

  expect_identical(breaks, Cdsw$new(masses, method = "quantile")$breaks)
})

test_that("Cdsw$new(method = \"sampling\") errors when equal-width windows violate the limits", {
  data(masses, package = "prozor", envir = environment())

  # 25 equal windows over 350-1200 m/z are 34 wide
  expect_error(Cdsw$new(masses, maxwindow = 10), "No window layout satisfies the limits")
  expect_no_error(Cdsw$new(masses, maxwindow = 10, method = "quantile"))
})

test_that("Cdsw$new() stops for invalid arguments, naming the argument and the value", {
  data(masses, package = "prozor", envir = environment())

  expect_error(Cdsw$new(masses, method = "linear"), "`method` must be one of \"sampling\", .*; not \"linear\"")
  expect_error(Cdsw$new(as.character(masses)), "`masses` must be a non-empty numeric vector, not character")
  expect_error(Cdsw$new(numeric()), "`masses` must be a non-empty numeric vector")
  expect_error(Cdsw$new(c(400, NA, 500)), "`masses` must not contain NA, NaN or infinite values; found 1")
  expect_error(Cdsw$new(masses, nbins = 0), "`nbins` must be a single positive whole number, not 0")
  expect_error(Cdsw$new(masses, nbins = 2.5), "`nbins` must be a single positive whole number, not 2.5")
  expect_error(Cdsw$new(masses, maxwindow = NA), "`maxwindow` must be a single positive number, not NA")
  expect_error(Cdsw$new(masses, minwindow = -1), "`minwindow` must be a single non-negative number, not -1")
  expect_error(
    Cdsw$new(masses, digits = c(1, 2)),
    "`digits` must be a single non-negative whole number, not c\\(1, 2\\)"
  )
})

test_that("Cdsw fields are read-only", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, nbins = 10, method = "quantile")
  breaks <- cdsw$breaks

  expect_error(cdsw$breaks <- seq(400, 1400, by = 100), "`breaks` of a Cdsw object is read-only")
  expect_error(cdsw$masses <- 500, "`masses` of a Cdsw object is read-only")
  expect_error(cdsw$nbins <- 5, "`nbins` of a Cdsw object is read-only")
  expect_error(cdsw$method <- "constant", "`method` of a Cdsw object is read-only")
  expect_identical(cdsw$breaks, breaks)
  expect_identical(cdsw$nbins, 10)
  expect_identical(cdsw$masses, masses)
})

test_that("Cdsw results do not depend on the order of the method calls", {
  data(masses, package = "prozor", envir = environment())
  first <- Cdsw$new(masses, maxwindow = 100)
  second <- Cdsw$new(masses, maxwindow = 100)

  table_first <- first$as_table()
  error_first <- first$error()
  optimized_first <- first$optimize_windows(maxbin = 10)
  first$plot_breaks()
  first$show_cycle()

  second$plot_breaks()
  optimized_second <- second$optimize_windows(maxbin = 10)
  second$show_cycle(overlap = 3)
  error_second <- second$error()
  table_second <- second$as_table()

  expect_identical(table_second, table_first)
  expect_identical(error_second, error_first)
  expect_identical(optimized_second, optimized_first)
  expect_identical(first$breaks, second$breaks)
})

test_that("Cdsw$clone() returns an independent object with the same state", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, nbins = 10, method = "quantile")

  cloned <- cdsw$clone()

  expect_false(identical(cloned, cdsw))
  expect_identical(cloned$breaks, cdsw$breaks)
  expect_identical(cloned$masses, cdsw$masses)
  expect_identical(cloned$nbins, cdsw$nbins)
  expect_identical(cloned$method, cdsw$method)
  expect_identical(cloned$as_table(), cdsw$as_table())
  # the state is private and read-only; changing it in the clone through R6 internals leaves the original as it was
  cloned$.__enclos_env__$private$.breaks <- seq(300, 1300, by = 100)
  expect_identical(cdsw$breaks, Cdsw$new(masses, nbins = 10, method = "quantile")$breaks)
})

test_that("Cdsw$print() shows the masses, windows, method and window widths", {
  cdsw <- Cdsw$new(c(405, 410, 420, 430, 435), nbins = 2, method = "constant", digits = 0)

  expect_output(
    expect_invisible(print(cdsw)),
    "<Cdsw> 5 masses in 2 windows, method \"constant\"\nwindow widths 16 to 16 m/z"
  )
})

test_that("Cdsw$as_table() returns one row per window", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, nbins = 10, method = "constant")

  tab <- cdsw$as_table(overlap = 1)

  expect_named(tab, c("from", "to", "mid", "width", "counts"))
  expect_equal(nrow(tab), 10)
  expect_equal(tab$from, cdsw$breaks[1:10] - 0.5)
  expect_equal(tab$to, cdsw$breaks[2:11] + 0.5)
  expect_equal(tab$counts[1], sum(masses > tab$from[1] & masses <= tab$to[1]))
})

test_that("Cdsw$as_table(), $plot() and $error() count the masses in (from, to]", {
  # digits = 0 gives the breaks 404, 420 and 436, so the mass 420 lies on the middle break
  cdsw <- Cdsw$new(c(405, 410, 420, 430, 435), nbins = 2, method = "constant", digits = 0)

  expect_equal(cdsw$breaks, c(404, 420, 436))
  expect_equal(cdsw$as_table(overlap = 0)$counts, c(3, 2))
  expect_equal(cdsw$as_table(overlap = 20)$counts, c(4, 3))
  expect_equal(ggplot2::layer_data(cdsw$plot(), 1)$ymax, c(3, 2))
  expect_equal(cdsw$error()$counts, c(3, 2) / sqrt(13))
})

test_that("Cdsw$plot() draws one bar per window with the precursor counts", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, nbins = 10, method = "constant")

  p <- cdsw$plot()

  expect_s3_class(p, "ggplot")
  bars <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(bars), 10)
  expect_equal(bars$xmin, cdsw$breaks[1:10])
  expect_equal(sum(bars$ymax), length(masses))
})

test_that("Cdsw$plot_breaks() compares the breaks with constant and quantile breaks", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, nbins = 10, maxwindow = 100)

  p <- cdsw$plot_breaks()

  expect_s3_class(p, "ggplot")
  points <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(points), 3 * 11)
  expect_length(unique(points$colour), 3)
  expect_equal(
    points$x,
    c(
      Cdsw$new(masses, nbins = 10, method = "constant")$breaks,
      Cdsw$new(masses, nbins = 10, method = "quantile")$breaks,
      cdsw$breaks
    )
  )
})

test_that("Cdsw$show_cycle() draws one rectangle per window", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, nbins = 10, method = "constant")

  p <- cdsw$show_cycle(overlap = 1)

  expect_s3_class(p, "ggplot")
  rects <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(rects), 10)
  expect_equal(rects$ymin, cdsw$as_table(overlap = 1)$from)
})

test_that("readjust_windows() moves boundaries but keeps one window per row", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw$new(masses, maxwindow = 100)
  tab <- cdsw$as_table()

  res <- readjust_windows(tab, masses, maxbin = 10)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("from", "to", "mid", "width", "counts"))
  expect_equal(nrow(res), nrow(tab))
  expect_true(all(res$to > res$from))
  expect_true(all(abs(res$from - tab$from) < 1))
  expect_equal(res$counts[1], sum(masses > res$from[1] & masses <= res$to[1]))
  expect_identical(cdsw$optimize_windows(maxbin = 10, overlap = 1), res)
})

test_that("readjust_windows() moves starts down and ends up to the nearest empty bin, also the first and last", {
  # 0.1 m/z bins with 10 masses each, except the empty bins at 399.95, 450.05 and 500.05
  ms1 <- seq(400.005, 500, by = 0.01)
  ms1 <- ms1[ms1 <= 450 | ms1 > 450.1]
  wind <- data.frame(from = c(400, 450.1), to = c(450, 500))

  res <- readjust_windows(wind, ms1, maxbin = 5)

  expect_equal(res$from, c(399.95, 450.05))
  expect_equal(res$to, c(450.05, 500.05))
  expect_equal(res$counts, c(sum(ms1 <= 450), sum(ms1 > 450.1)))
})

test_that("readjust_windows() keeps the windows when no bin is empty", {
  ms1 <- c(405, 410, 420, 430, 435)
  wind <- data.frame(from = c(400, 420), to = c(420, 440))

  res <- readjust_windows(wind, ms1, maxbin = 0)

  expect_equal(res$from, wind$from)
  expect_equal(res$to, wind$to)
  expect_equal(res$counts, c(3, 2))
})

test_that("readjust_windows() bins the masses for any bin width", {
  ms1 <- c(400.12, 512.37, 700.5, 1000.04)
  wind <- data.frame(from = c(400, 700), to = c(700, 1000.1))

  res <- readjust_windows(wind, ms1, bin_width = 0.01, maxbin = 1)

  expect_equal(res$from, c(400, 699.995))
  expect_equal(res$to, c(700.005, 1000.1))
  expect_equal(res$counts, c(2, 2))
  expect_no_error(
    for (min_mass in seq(400.001, 400.099, by = 0.001)) {
      readjust_windows(list(from = 500, to = 600), c(min_mass, 500, 700), bin_width = 0.01, maxbin = 1)
    }
  )
  # bins of 1 m/z from 399 to 404: the empty bins have the mids 399.5, 401.5 and 403.5
  res <- readjust_windows(list(from = 401, to = 402), c(400.2, 400.7, 402.5), bin_width = 1, maxbin = 1)
  expect_equal(c(res$from, res$to, res$counts), c(399.5, 403.5, 3))
  expect_error(readjust_windows(wind, ms1, bin_width = 0), "`bin_width` must be a single positive number, not 0")
})

test_that("plot_windows() draws the masses with a start and an end line per window", {
  data(masses, package = "prozor", envir = environment())
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
