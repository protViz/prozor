# Moves each window start down to the nearest empty-bin mid at or below it, so a start already on such a mid stays;
# starts without an empty bin below stay too. `empty` holds the sorted mids of the empty bins.
.move_starts <- function(starts, empty) {
  below <- findInterval(starts, empty)
  starts[below > 0] <- empty[below[below > 0]]
  starts
}

# Moves each window end up to the nearest empty-bin mid at or above it, so an end already on such a mid stays; ends
# without an empty bin above stay too
.move_ends <- function(ends, empty) {
  above <- findInterval(ends, empty, left.open = TRUE) + 1L
  moved <- above <= length(empty)
  ends[moved] <- empty[above[moved]]
  ends
}

# Stops unless `windows` is a data frame or list with the columns `from` and `to`
.check_windows <- function(windows) {
  if (!is.list(windows)) {
    stop("`windows` must be a data frame or list, not ", class(windows)[1], ".", call. = FALSE)
  }
  .check_columns(windows, c("from", "to"), "windows")
}

#' Move DIA window boundaries to regions with few precursors
#'
#' Moves each window start down, and each window end up, to the mid of the nearest empty bin. A boundary already on
#' such a mid stays, so `readjust_windows()` returns its own result unchanged.
#'
#' @param windows A data frame or list with the columns `from` and `to`, e.g. from `Cdsw$as_table()`.
#' @param masses Precursor masses (m/z): a numeric vector of finite values.
#' @param bin_width Width of the bins (m/z) into which the masses are binned.
#' @param min_count Bins with fewer than `min_count` precursors count as empty; a number, at least 0.
#' @return A data frame with one row per window and the columns `from`, `to`, `mid`, `width` and `n_precursors`
#'   (precursors with `from < m/z <= to`).
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw$new(masses, max_window = 100)
#' windows <- cdsw$as_table()
#' head(readjust_windows(windows, masses, min_count = 10))
readjust_windows <- function(windows, masses, bin_width = 0.1, min_count = 15) {
  .check_windows(windows)
  .check_masses(masses)
  .check_number(bin_width, "bin_width", positive = TRUE, whole = FALSE)
  .check_number(min_count, "min_count", positive = FALSE, whole = FALSE)
  mass_range <- .padded_range(masses, 1 / bin_width)
  breaks <- seq(mass_range[1], mass_range[2], by = bin_width)
  reshist <- graphics::hist(masses, breaks = breaks, plot = FALSE)
  empty <- reshist$mids[reshist$counts < min_count]
  new_from <- .move_starts(windows$from, empty)
  new_to <- .move_ends(windows$to, empty)
  n_precursors <- .window_counts(masses, new_from, new_to)
  data.frame(
    from = new_from,
    to = new_to,
    mid = (new_from + new_to) * 0.5,
    width = new_to - new_from,
    n_precursors = n_precursors
  )
}

#' Plot precursor masses with DIA window boundaries
#'
#' @param windows A data frame or list with the columns `from` and `to`, e.g. from `Cdsw$as_table()` or
#'   [readjust_windows()].
#' @inheritParams readjust_windows
#' @param bin_width Width of the histogram bins (m/z).
#' @return A ggplot object.
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw$new(masses, max_window = 100)
#' windows <- cdsw$optimize_windows(min_count = 10)
#' plot_windows(windows, masses)
#' # zoom into one window
#' start <- windows$from[5]
#' plot_windows(windows, masses) + ggplot2::coord_cartesian(xlim = c(start - 3, start + 3))
plot_windows <- function(windows, masses, bin_width = 0.1) {
  .check_windows(windows)
  .check_masses(masses)
  .check_number(bin_width, "bin_width", positive = TRUE, whole = FALSE)
  boundaries <- data.frame(
    mz = c(windows$from, windows$to),
    boundary = rep(c("window start", "window end"), each = length(windows$from))
  )
  ggplot2::ggplot(data.frame(mz = masses), ggplot2::aes(x = .data$mz)) +
    ggplot2::geom_histogram(binwidth = bin_width, fill = "grey40") +
    ggplot2::geom_vline(
      data = boundaries,
      ggplot2::aes(xintercept = .data$mz, colour = .data$boundary, linetype = .data$boundary)
    ) +
    ggplot2::scale_colour_manual(values = c("window start" = "#0072B2", "window end" = "#D55E00"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("window start" = "solid", "window end" = "dashed"), name = NULL) +
    ggplot2::labs(x = "precursor m/z", y = "precursors per bin")
}
