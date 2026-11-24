.adjust_breaks <- function(breaks, digits) {
  breaks[1] <- breaks[1] - 1 / (10^digits)
  breaks[length(breaks)] <- breaks[length(breaks)] + 1 / (10^digits)
  round(breaks, digits = digits)
}


#' Deviation of a window layout from equal numbers of MS1 masses per window
#'
#' @param splits The window boundaries.
#' @param data The MS1 masses.
#' @return A list with `score1` (Manhattan distance), `score2` (Euclidean distance), `counts` (observed counts)
#'   and `optimumN` (optimal counts), both counts normalised.
#' @keywords internal
objectiveMS1Function <- function(splits, data) {
  counts <- graphics::hist(data, breaks = splits, plot = FALSE)$counts
  optimum_n <- rep(length(data) / (length(splits) - 1), length(splits) - 1)
  counts <- counts / sqrt(sum(counts^2))
  optimum_n <- optimum_n / sqrt(sum(optimum_n^2))
  list(
    score1 = sum(abs(counts / max(abs(counts)) - optimum_n / max(abs(optimum_n)))) / length(counts),
    score2 = sqrt(sum((counts - optimum_n)^2)) / length(counts),
    counts = counts,
    optimumN = optimum_n
  )
}

# TRUE if all window widths lie within [minwindow, maxwindow]
hardconstrain <- function(splits, minwindow = 5, maxwindow = 50) {
  widths <- diff(splits)
  all(widths >= minwindow & widths <= maxwindow)
}


#' Compute dynamic SWATH (DIA) windows
#'
#' A reference class that splits the MS1 mass range into windows with about equal numbers of precursors,
#' subject to minimum and maximum window widths.
#'
#' @field masses MS1 masses.
#' @field breaks The window boundaries.
#' @field nbins Number of windows.
#' @field digits Mass accuracy of the result.
#' @import methods
#' @importFrom stats quantile runif
#' @include readjust_windows.R
#' @family DIA windows
#' @export Cdsw
#' @exportClass Cdsw
#' @examples
#' data(masses)
#' cdsw <- Cdsw(masses)
#' breaks <- cdsw$sampling_breaks(maxwindow = 100)
#' cdsw$plot_breaks()
#' cdsw$plot()
#' cdsw$as_table()
#' cdsw$optimize_windows()
#' cdsw$show_cycle()
Cdsw <- setRefClass(
  "Cdsw",
  fields = list(
    masses = "numeric",
    breaks = "numeric",
    nbins = "numeric",
    digits = "numeric"
  ),
  methods = list(
    initialize = function(masses, nbins = 25, digits = 1) {
      "Create the object from a vector of MS1 masses, with nbins windows of equal width and mass accuracy digits."
      .self$masses <- masses
      .self$nbins <- nbins
      .self$digits <- digits
      constant_breaks()
    },
    constant_breaks = function(digits = 2) {
      "Equal-width windows; returns the breaks invisibly."
      minmass <- round(min(.self$masses) - 1 / 10^digits, digits = 2)
      maxmass <- round(max(.self$masses) + 1 / 10^digits, digits = 2)
      .self$breaks <- round(seq(minmass, maxmass, length.out = nbins + 1), digits = digits)
      invisible(.self$breaks)
    },
    quantile_breaks = function(digits = 2) {
      "Windows with equal numbers of MS1 masses; widths may violate the limits of sampling_breaks()."
      qqs <- stats::quantile(.self$masses, seq(0, 1, length.out = nbins + 1))
      .self$breaks <- .adjust_breaks(qqs, digits = digits)
      invisible(.self$breaks)
    },
    sampling_breaks = function(maxwindow = 150, minwindow = 5, digits = 2) {
      "Quantile windows, mixed with uniform data until all widths lie within [minwindow, maxwindow].
      Stops if even equal-width windows violate the limits; returns the breaks invisibly."
      qqs <- stats::quantile(masses, probs = seq(0, 1, by = 1 / nbins))
      unif <- seq(min(masses), max(masses), length.out = nbins + 1)
      if (!hardconstrain(unif, minwindow = minwindow, maxwindow = maxwindow)) {
        stop(
          "No window layout satisfies the limits: ",
          nbins,
          " windows of equal width are ",
          round(diff(unif)[1], digits = 2),
          " wide, outside [minwindow, maxwindow] = [",
          minwindow,
          ", ",
          maxwindow,
          "]."
        )
      }
      mixeddata <- masses
      while (!hardconstrain(qqs, minwindow, maxwindow)) {
        uniformdata <- stats::runif(round(length(masses) / 20), min = min(masses), max = max(masses))
        mixeddata <- c(mixeddata, uniformdata)
        qqs <- stats::quantile(mixeddata, probs = seq(0, 1, by = 1 / nbins))
      }
      .self$breaks <- .adjust_breaks(qqs, digits = digits)
      invisible(.self$breaks)
    },
    plot_breaks = function() {
      "ggplot of the current breaks against equal-width and quantile breaks."
      index <- seq(0, nbins)
      layouts <- c("equal width", "quantile", "current")
      breaks_long <- data.frame(
        mz = c(
          seq(min(.self$masses), max(.self$masses), length.out = nbins + 1),
          stats::quantile(.self$masses, seq(0, 1, length.out = nbins + 1), names = FALSE),
          .self$breaks
        ),
        index = rep(index, 3),
        layout = factor(rep(layouts, each = nbins + 1), levels = layouts)
      )
      ggplot2::ggplot(breaks_long, ggplot2::aes(x = .data$mz, y = .data$index, colour = .data$layout)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 1) +
        ggplot2::scale_colour_manual(
          values = c("equal width" = "#999999", "quantile" = "#E69F00", "current" = "#0072B2"),
          name = "breaks"
        ) +
        ggplot2::labs(x = "precursor m/z", y = "break index")
    },
    plot = function() {
      "ggplot of the precursors per window; the dashed line marks equal counts."
      counts <- graphics::hist(.self$masses, breaks = .self$breaks, plot = FALSE)$counts
      windows <- data.frame(from = utils::head(.self$breaks, -1), to = utils::tail(.self$breaks, -1), count = counts)
      ggplot2::ggplot(windows) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin = .data$from, xmax = .data$to, ymin = 0, ymax = .data$count),
          fill = "#56B4E9",
          colour = "white"
        ) +
        ggplot2::geom_hline(yintercept = length(.self$masses) / nbins, linetype = "dashed") +
        ggplot2::labs(x = "precursor m/z", y = "precursors per window")
    },
    as_table = function(overlap = 1) {
      "Windows as a data.frame with start (from), end (to), centre (mid), width and precursor count; overlap in m/z."
      q <- .self$breaks
      idx <- seq_len(length(q) - 1)
      from <- q[idx] - overlap * 0.5
      to <- q[idx + 1] + overlap * 0.5
      width <- to - from
      counts <- vapply(seq_along(from), function(i) sum(.self$masses > from[i] & .self$masses < to[i]), integer(1))
      data.frame(from = from, to = to, mid = from + width * 0.5, width = width, counts = counts)
    },
    error = function() {
      "Deviation of the windows from equal counts: a list, see objectiveMS1Function()."
      objectiveMS1Function(.self$breaks, .self$masses)
    },
    optimize_windows = function(digits = 1, maxbin = 15, overlap = 0) {
      "Move window starts and ends to regions with few precursor masses; returns a data.frame."
      data.frame(readjust_windows(as_table(overlap = overlap), masses, digits = digits, maxbin = maxbin))
    },
    show_cycle = function(overlap = 1) {
      "ggplot of the DIA cycle, one window per step; overlap in m/z."
      windows <- .self$as_table(overlap = overlap)
      windows$step <- seq_len(nrow(windows))
      windows$alternate <- factor(windows$step %% 2)
      ggplot2::ggplot(windows) +
        ggplot2::geom_rect(
          ggplot2::aes(
            xmin = .data$step,
            xmax = .data$step + 1,
            ymin = .data$from,
            ymax = .data$to,
            fill = .data$alternate
          ),
          colour = "black",
          alpha = 0.6,
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(values = c("0" = "#56B4E9", "1" = "#E69F00")) +
        ggplot2::labs(x = "window", y = "m/z", title = "DIA cycle")
    }
  )
)
