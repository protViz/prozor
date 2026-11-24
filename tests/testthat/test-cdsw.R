test_that("Cdsw starts with nbins equal-width windows over the mass range", {
  data(masses, package = "prozor", envir = environment())

  cdsw <- Cdsw(masses, nbins = 10)

  expect_length(cdsw$breaks, 11)
  expect_lt(cdsw$breaks[1], min(masses))
  expect_gt(cdsw$breaks[11], max(masses))
  expect_equal(diff(cdsw$breaks), rep(diff(cdsw$breaks)[1], 10), tolerance = 0.01)
})

test_that("Cdsw$constant_breaks() encloses the masses for any digits", {
  min_masses <- seq(400.001, 400.099, by = 0.001)

  first_breaks <- vapply(
    min_masses,
    function(min_mass) Cdsw(c(min_mass, 450, 499.99), nbins = 2)$constant_breaks(digits = 3)[1],
    numeric(1)
  )

  expect_true(all(first_breaks < min_masses))
  cdsw <- Cdsw(c(400.008, 450, 499.99), nbins = 2)
  cdsw$constant_breaks(digits = 3)
  expect_equal(sum(ggplot2::layer_data(cdsw$plot(), 1)$ymax), 3)
  expect_equal(sum(cdsw$error()$counts > 0), 2)
})

test_that("Cdsw$sampling_breaks() keeps window widths within the limits", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)

  breaks <- cdsw$sampling_breaks(maxwindow = 100, minwindow = 5)

  expect_length(breaks, 26)
  expect_equal(cdsw$breaks, breaks)
  # rounding to digits = 2 moves a break by up to 0.005; the outer breaks are padded by up to 0.02
  expect_true(all(diff(breaks) >= 5))
  expect_true(all(diff(breaks) <= 100 + 0.025))
  # the smallest share of uniform data that meets the limits: the widest window is at the limit
  expect_equal(max(diff(breaks)), 100, tolerance = 0.025 / 100)
})

test_that("Cdsw$sampling_breaks() keeps the precursor counts per window close to equal", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)

  cdsw$sampling_breaks(maxwindow = 100)

  counts <- cdsw$as_table(overlap = 0)$counts
  expect_lt(stats::sd(counts) / mean(counts), 0.2)
})

test_that("Cdsw$sampling_breaks() is reproducible and does not use the random number generator", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)
  withr::local_seed(1)
  seed <- .Random.seed

  first <- cdsw$sampling_breaks(maxwindow = 100)

  expect_identical(.Random.seed, seed)
  expect_identical(cdsw$sampling_breaks(maxwindow = 100), first)
})

test_that("Cdsw$sampling_breaks() meets limits close to the equal width", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)
  equal_width <- diff(range(masses)) / 25

  breaks <- cdsw$sampling_breaks(maxwindow = 1.02 * equal_width, minwindow = 0.98 * equal_width)

  expect_true(all(diff(breaks) >= 0.98 * equal_width - 0.01))
  expect_true(all(diff(breaks) <= 1.02 * equal_width + 0.025))
})

test_that("Cdsw$sampling_breaks() gives the quantile breaks when they meet the limits", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)

  breaks <- cdsw$sampling_breaks(maxwindow = 1000, minwindow = 1)

  expect_identical(breaks, cdsw$quantile_breaks())
})

test_that("Cdsw$sampling_breaks() errors when equal-width windows violate the limits", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)

  # 25 equal windows over 350-1200 m/z are 34 wide
  expect_error(cdsw$sampling_breaks(maxwindow = 10), "No window layout satisfies the limits")
})

test_that("Cdsw$copy() returns an independent object with the same breaks", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)
  cdsw$quantile_breaks()

  copied <- cdsw$copy()

  expect_identical(copied$breaks, cdsw$breaks)
  expect_identical(copied$masses, cdsw$masses)
  expect_identical(copied$nbins, cdsw$nbins)
  copied$constant_breaks()
  expect_identical(cdsw$breaks, Cdsw(masses, nbins = 10)$quantile_breaks())
})

test_that("Cdsw$as_table() returns one row per window", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)

  tab <- cdsw$as_table(overlap = 1)

  expect_named(tab, c("from", "to", "mid", "width", "counts"))
  expect_equal(nrow(tab), 10)
  expect_equal(tab$from, cdsw$breaks[1:10] - 0.5)
  expect_equal(tab$to, cdsw$breaks[2:11] + 0.5)
  expect_equal(tab$counts[1], sum(masses > tab$from[1] & masses <= tab$to[1]))
})

test_that("Cdsw$as_table(), $plot() and $error() count the masses in (from, to]", {
  cdsw <- Cdsw(c(405, 410, 420, 430, 435), nbins = 2)
  cdsw$breaks <- c(400, 420, 440)

  expect_equal(cdsw$as_table(overlap = 0)$counts, c(3, 2))
  expect_equal(cdsw$as_table(overlap = 20)$counts, c(4, 3))
  expect_equal(ggplot2::layer_data(cdsw$plot(), 1)$ymax, c(3, 2))
  expect_equal(cdsw$error()$counts, c(3, 2) / sqrt(13))
})

test_that("Cdsw$plot() draws one bar per window with the precursor counts", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)

  p <- cdsw$plot()

  expect_s3_class(p, "ggplot")
  bars <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(bars), 10)
  expect_equal(bars$xmin, cdsw$breaks[1:10])
  expect_equal(sum(bars$ymax), length(masses))
})

test_that("Cdsw$plot_breaks() compares three break layouts", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)

  p <- cdsw$plot_breaks()

  expect_s3_class(p, "ggplot")
  points <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(points), 3 * 11)
  expect_length(unique(points$colour), 3)
})

test_that("Cdsw$show_cycle() draws one rectangle per window", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)

  p <- cdsw$show_cycle(overlap = 1)

  expect_s3_class(p, "ggplot")
  rects <- ggplot2::layer_data(p, 1)
  expect_equal(nrow(rects), 10)
  expect_equal(rects$ymin, cdsw$as_table(overlap = 1)$from)
})

test_that("readjust_windows() moves boundaries but keeps one window per row", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)
  cdsw$sampling_breaks(maxwindow = 100)
  tab <- cdsw$as_table()

  res <- readjust_windows(tab, masses, maxbin = 10)

  expect_named(res, c("from", "to", "mid", "width", "counts"))
  expect_length(res$from, nrow(tab))
  expect_true(all(res$to > res$from))
  expect_true(all(abs(res$from - tab$from) < 1))
  expect_equal(res$counts[1], sum(masses > res$from[1] & masses <= res$to[1]))
  expect_equal(cdsw$optimize_windows(maxbin = 10, overlap = 1), data.frame(res))
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

test_that("readjust_windows() bins the masses for any digits", {
  ms1 <- c(400.12, 512.37, 700.5, 1000.04)
  wind <- data.frame(from = c(400, 700), to = c(700, 1000.1))

  res <- readjust_windows(wind, ms1, digits = 2, maxbin = 1)

  expect_equal(res$from, c(400, 699.995))
  expect_equal(res$to, c(700.005, 1000.1))
  expect_equal(res$counts, c(2, 2))
  expect_no_error(
    for (min_mass in seq(400.001, 400.099, by = 0.001)) {
      readjust_windows(list(from = 500, to = 600), c(min_mass, 500, 700), digits = 2, maxbin = 1)
    }
  )
})

test_that("plot_windows() draws the masses with a start and an end line per window", {
  data(masses, package = "prozor", envir = environment())
  windows <- data.frame(from = c(400, 500), to = c(510, 610))

  p <- plot_windows(windows, masses)

  expect_s3_class(p, "ggplot")
  lines <- ggplot2::layer_data(p, 2)
  expect_setequal(lines$xintercept, c(400, 500, 510, 610))
  expect_equal(sum(ggplot2::layer_data(p, 1)$count), length(masses))
})
