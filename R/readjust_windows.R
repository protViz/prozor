# Moves each window start down to the nearest empty bin below it; starts without an empty bin below stay.
# `empty` holds the sorted mids of the empty bins.
.move_starts <- function(starts, empty) {
  below <- findInterval(starts, empty, left.open = TRUE)
  starts[below > 0] <- empty[below[below > 0]]
  starts
}

# Moves each window end up to the nearest empty bin above it; ends without an empty bin above stay
.move_ends <- function(ends, empty) {
  above <- findInterval(ends, empty) + 1L
  moved <- above <= length(empty)
  ends[moved] <- empty[above[moved]]
  ends
}

#' Move DIA window boundaries to regions with few precursors
#'
#' @param wind A data frame or list with the columns `from` and `to`, e.g. from `Cdsw$as_table()`.
#' @param ms1data Precursor masses (m/z).
#' @param bin_width Width of the bins (m/z) into which the masses are binned.
#' @param maxbin Bins with fewer than `maxbin` precursors count as empty.
#' @return A data frame with one row per window and the columns `from`, `to`, `mid`, `width` and `counts` (precursors
#'   with `from < m/z <= to`).
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw$new(masses, maxwindow = 100)
#' windows <- cdsw$as_table()
#' head(readjust_windows(windows, masses, maxbin = 10))
readjust_windows <- function(wind, ms1data, bin_width = 0.1, maxbin = 15) {
  .check_number(bin_width, "bin_width", positive = TRUE, whole = FALSE)
  mass_range <- .padded_range(ms1data, 1 / bin_width)
  breaks <- seq(mass_range[1], mass_range[2], by = bin_width)
  reshist <- graphics::hist(ms1data, breaks = breaks, plot = FALSE)
  empty <- reshist$mids[reshist$counts < maxbin]
  new_from <- .move_starts(wind$from, empty)
  new_to <- .move_ends(wind$to, empty)
  counts <- .window_counts(ms1data, new_from, new_to)
  data.frame(from = new_from, to = new_to, mid = (new_from + new_to) * 0.5, width = new_to - new_from, counts = counts)
}

#' Plot precursor masses with DIA window boundaries
#'
#' @param windows A data frame or list with the columns `from` and `to`, e.g. from `Cdsw$as_table()` or
#'   [readjust_windows()].
#' @param ms1data Precursor masses (m/z).
#' @param bin_width Width of the histogram bins (m/z).
#' @return A ggplot object.
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw$new(masses, maxwindow = 100)
#' windows <- cdsw$optimize_windows(maxbin = 10)
#' plot_windows(windows, masses)
#' # zoom into one window
#' start <- windows$from[5]
#' plot_windows(windows, masses) + ggplot2::coord_cartesian(xlim = c(start - 3, start + 3))
plot_windows <- function(windows, ms1data, bin_width = 0.1) {
  .check_number(bin_width, "bin_width", positive = TRUE, whole = FALSE)
  boundaries <- data.frame(
    mz = c(windows$from, windows$to),
    boundary = rep(c("window start", "window end"), each = length(windows$from))
  )
  ggplot2::ggplot(data.frame(mz = ms1data), ggplot2::aes(x = .data$mz)) +
    ggplot2::geom_histogram(binwidth = bin_width, fill = "grey40") +
    ggplot2::geom_vline(
      data = boundaries,
      ggplot2::aes(xintercept = .data$mz, colour = .data$boundary, linetype = .data$boundary)
    ) +
    ggplot2::scale_colour_manual(values = c("window start" = "#0072B2", "window end" = "#D55E00"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("window start" = "solid", "window end" = "dashed"), name = NULL) +
    ggplot2::labs(x = "precursor m/z", y = "precursors per bin")
}
