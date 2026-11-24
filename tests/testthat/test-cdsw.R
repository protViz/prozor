test_that("Cdsw starts with nbins equal-width windows over the mass range", {
  data(masses, package = "prozor", envir = environment())

  cdsw <- Cdsw(masses, nbins = 10)

  expect_length(cdsw$breaks, 11)
  expect_lt(cdsw$breaks[1], min(masses))
  expect_gt(cdsw$breaks[11], max(masses))
  expect_equal(diff(cdsw$breaks), rep(diff(cdsw$breaks)[1], 10), tolerance = 0.01)
})

test_that("Cdsw$sampling_breaks() keeps window widths within the limits", {
  data(masses, package = "prozor", envir = environment())
  withr::local_seed(1)
  cdsw <- Cdsw(masses)

  breaks <- cdsw$sampling_breaks(maxwindow = 100, minwindow = 5)

  expect_length(breaks, 26)
  expect_equal(cdsw$breaks, breaks)
  # .adjust_breaks() widens the outer windows by 1 / 10^digits
  expect_true(all(diff(breaks) >= 5))
  expect_true(all(diff(breaks) <= 100 + 0.02))
})

test_that("Cdsw$sampling_breaks() errors when equal-width windows violate the limits", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)

  # 25 equal windows over 350-1200 m/z are 34 wide
  expect_error(cdsw$sampling_breaks(maxwindow = 10), "No window layout satisfies the limits")
})

test_that("Cdsw$as_table() returns one row per window", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)

  tab <- cdsw$as_table(overlap = 1)

  expect_named(tab, c("from", "to", "mid", "width", "counts"))
  expect_equal(nrow(tab), 10)
  expect_equal(tab$from, cdsw$breaks[1:10] - 0.5)
  expect_equal(tab$to, cdsw$breaks[2:11] + 0.5)
  expect_equal(tab$counts[1], sum(masses > tab$from[1] & masses < tab$to[1]))
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
  withr::local_seed(1)
  cdsw <- Cdsw(masses)
  cdsw$sampling_breaks(maxwindow = 100)
  tab <- cdsw$as_table()

  res <- readjust_windows(tab, masses, maxbin = 10)

  expect_named(res, c("from", "to", "mid", "width", "counts"))
  expect_length(res$from, nrow(tab))
  expect_true(all(res$to > res$from))
  expect_true(all(abs(res$from - tab$from) < 1))
  expect_equal(res$counts[1], sum(masses >= res$from[1] & masses <= res$to[1]))
  expect_equal(cdsw$optimize_windows(maxbin = 10, overlap = 1), data.frame(res))
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
