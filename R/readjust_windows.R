# Moves each window boundary to the nearest bin with few peaks: starts move down, ends move up.
.make_new_from_to <- function(boundaries, empty, is_from = TRUE) {
  vapply(
    boundaries,
    function(boundary) {
      idx <- which.min(abs(boundary - empty))
      if (is_from) {
        if (idx == 1) {
          boundary
        } else if (empty[idx] < boundary) {
          empty[idx]
        } else {
          empty[idx - 1]
        }
      } else {
        if (idx == length(empty)) {
          boundary
        } else if (empty[idx] > boundary) {
          empty[idx]
        } else {
          empty[idx + 1]
        }
      }
    },
    numeric(1)
  )
}

#' Move DIA window boundaries to regions with few precursors
#'
#' @param wind A data frame or list with the columns `from` and `to`, e.g. from `Cdsw$as_table()`.
#' @param ms1data Precursor masses (m/z).
#' @param digits Mass accuracy: the masses are binned into bins of width `1 / 10^digits`.
#' @param maxbin Bins with fewer than `maxbin` precursors count as empty.
#' @return A list with the elements `from`, `to`, `mid`, `width` and `counts` (precursors per window).
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw(masses)
#' breaks <- cdsw$sampling_breaks(maxwindow = 100)
#' windows <- cdsw$as_table()
#' head(data.frame(readjust_windows(windows, masses, maxbin = 10)))
readjust_windows <- function(wind, ms1data, digits = 1, maxbin = 15) {
  breaks <- seq(
    round(min(ms1data) - 1 / 10^digits, digits = 1),
    round(max(ms1data) + 1 / 10^digits, digits = 1),
    by = 1 / 10^digits
  )
  reshist <- graphics::hist(ms1data, breaks = breaks, plot = FALSE)
  empty <- reshist$mids[reshist$counts < maxbin]
  new_from <- .make_new_from_to(wind$from, empty)
  new_to <- .make_new_from_to(wind$to, empty, is_from = FALSE)
  counts <- vapply(
    seq_along(new_from),
    function(i) sum(ms1data >= new_from[i] & ms1data <= new_to[i]),
    integer(1)
  )
  list(from = new_from, to = new_to, mid = (new_from + new_to) * 0.5, width = new_to - new_from, counts = counts)
}

#' Plot precursor masses with DIA window boundaries
#'
#' @param windows A data frame or list with the columns `from` and `to`, e.g. from `Cdsw$as_table()` or
#'   [readjust_windows()].
#' @param ms1data Precursor masses (m/z).
#' @param digits The histogram bin width is `1 / 10^digits` m/z.
#' @return A ggplot object.
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw(masses)
#' breaks <- cdsw$sampling_breaks(maxwindow = 100)
#' windows <- readjust_windows(cdsw$as_table(), masses, maxbin = 10)
#' plot_windows(windows, masses)
#' # zoom into one window
#' start <- windows$from[5]
#' plot_windows(windows, masses) + ggplot2::coord_cartesian(xlim = c(start - 3, start + 3))
plot_windows <- function(windows, ms1data, digits = 1) {
  boundaries <- data.frame(
    mz = c(windows$from, windows$to),
    boundary = rep(c("window start", "window end"), each = length(windows$from))
  )
  ggplot2::ggplot(data.frame(mz = ms1data), ggplot2::aes(x = .data$mz)) +
    ggplot2::geom_histogram(binwidth = 1 / 10^digits, fill = "grey40") +
    ggplot2::geom_vline(
      data = boundaries,
      ggplot2::aes(xintercept = .data$mz, colour = .data$boundary, linetype = .data$boundary)
    ) +
    ggplot2::scale_colour_manual(values = c("window start" = "#0072B2", "window end" = "#D55E00"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("window start" = "solid", "window end" = "dashed"), name = NULL) +
    ggplot2::labs(x = "precursor m/z", y = "precursors per bin")
}
