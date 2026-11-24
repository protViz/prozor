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
  # .adjustBreaks() widens the outer windows by 1 / 10^digits
  expect_true(all(diff(breaks) >= 5))
  expect_true(all(diff(breaks) <= 100 + 0.02))
})

test_that("Cdsw$sampling_breaks() errors when equal-width windows violate the limits", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses)

  # 25 equal windows over 350-1200 m/z are 34 wide
  expect_error(cdsw$sampling_breaks(maxwindow = 10), "no windows satisfy the constraints")
})

test_that("Cdsw$asTable() returns one row per window", {
  data(masses, package = "prozor", envir = environment())
  cdsw <- Cdsw(masses, nbins = 10)

  tab <- cdsw$asTable(overlap = 1)

  expect_named(tab, c("from", "to", "mid", "width", "counts"))
  expect_equal(nrow(tab), 10)
  expect_equal(tab$from, cdsw$breaks[1:10] - 0.5)
  expect_equal(tab$to, cdsw$breaks[2:11] + 0.5)
  expect_equal(tab$counts[1], sum(masses > tab$from[1] & masses < tab$to[1]))
})

test_that("readjustWindows() moves boundaries but keeps one window per row", {
  data(masses, package = "prozor", envir = environment())
  withr::local_seed(1)
  cdsw <- Cdsw(masses)
  cdsw$sampling_breaks(maxwindow = 100)
  tab <- cdsw$asTable()

  res <- readjustWindows(tab, masses, maxbin = 10)

  expect_named(res, c("from", "to", "mid", "width", "counts"))
  expect_length(res$from, nrow(tab))
  expect_true(all(res$to > res$from))
  expect_true(all(abs(res$from - tab$from) < 1))
  expect_equal(res$counts[1], sum(masses >= res$from[1] & masses <= res$to[1]))
})
