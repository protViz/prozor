# Range of x rounded outward to `digits` decimals and padded by 1 / 10^digits, so that it encloses x
.padded_range <- function(x, digits) {
  scale <- 10^digits
  c(floor(min(x) * scale) - 1, ceiling(max(x) * scale) + 1) / scale
}

# Breaks rounded to `digits`; the outer breaks are padded so that they enclose the masses
.adjust_breaks <- function(breaks, digits) {
  outer <- .padded_range(breaks, digits)
  breaks <- round(breaks, digits = digits)
  breaks[c(1, length(breaks))] <- outer
  breaks
}


#' Deviation of a window layout from equal numbers of MS1 masses per window
#'
#' @param splits The window boundaries.
#' @param data The MS1 masses.
#' @return A list with `score1` (Manhattan distance), `score2` (Euclidean distance), `counts` (observed counts)
#'   and `optimum_n` (optimal counts), both counts normalized.
#' @keywords internal
objective_ms1 <- function(splits, data) {
  counts <- .window_counts(data, utils::head(splits, -1), utils::tail(splits, -1))
  counts <- counts / sqrt(sum(counts^2))
  optimum_n <- rep(1 / sqrt(length(counts)), length(counts))
  list(
    score1 = mean(1 - counts / max(counts)),
    score2 = sqrt(sum((counts - optimum_n)^2)) / length(counts),
    counts = counts,
    optimum_n = optimum_n
  )
}

# TRUE if all window widths lie within [minwindow, maxwindow]
.hard_constrain <- function(splits, minwindow, maxwindow) {
  widths <- diff(splits)
  all(widths >= minwindow & widths <= maxwindow)
}

# Number of masses x with from < x <= to, for each window (the intervals of hist())
.window_counts <- function(x, from, to) {
  sorted <- sort(x)
  findInterval(to, sorted) - findInterval(from, sorted)
}

# Breaks of nbins windows with equal numbers of masses
.quantile_breaks <- function(x, nbins) {
  stats::quantile(x, seq(0, 1, length.out = nbins + 1), names = FALSE)
}

# Breaks at the quantiles of the mixture (1 - share) * masses + share * uniform over the mass range: the quantile breaks
# of the masses plus equally spaced points that make up `share` of the data, in the limit of infinitely many points.
# `sorted` holds the sorted masses; 0 < share <= 1.
.mixture_breaks <- function(sorted, nbins, share) {
  n <- length(sorted)
  span <- sorted[n] - sorted[1]
  # the mixture CDF jumps by (1 - share) / n at each mass and rises with slope share / span in between
  at <- c(sorted[1], sorted)
  cdf <- c(0, (1 - share) * seq_len(n) / n + share * (sorted - sorted[1]) / span)
  probs <- seq_len(nbins - 1) / nbins
  k <- findInterval(probs, cdf)
  inner <- pmin(at[k] + (probs - cdf[k]) * span / share, at[k + 1])
  c(sorted[1], inner, sorted[n])
}

# Breaks with all widths within [minwindow, maxwindow]: the quantile breaks if they fit, otherwise the mixture breaks
# with the smallest share of uniform data that fits, found by bisection to 1e-9. `equal` (share 1) must fit; the
# search keeps the last layout that fit, so it returns one even if fitting is not monotone in the share.
.sampling_breaks <- function(x, nbins, equal, minwindow, maxwindow) {
  quantiles <- .quantile_breaks(x, nbins)
  if (.hard_constrain(quantiles, minwindow = minwindow, maxwindow = maxwindow)) {
    return(quantiles)
  }
  sorted <- sort(x)
  fitting <- equal
  lower <- 0
  upper <- 1
  while (upper - lower > 1e-9) {
    share <- (lower + upper) / 2
    breaks <- .mixture_breaks(sorted, nbins, share)
    if (.hard_constrain(breaks, minwindow = minwindow, maxwindow = maxwindow)) {
      fitting <- breaks
      upper <- share
    } else {
      lower <- share
    }
  }
  fitting
}


#' Compute dynamic SWATH (DIA) windows
#'
#' A reference class that splits the MS1 mass range into windows with about equal numbers of precursors,
#' subject to minimum and maximum window widths.
#'
#' @field masses MS1 masses.
#' @field breaks The window boundaries.
#' @field nbins Number of windows.
#' @import methods
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
    nbins = "numeric"
  ),
  methods = list(
    initialize = function(masses, nbins = 25) {
      "Create the object from a vector of MS1 masses, with nbins windows of equal width."
      .self$masses <- masses
      .self$nbins <- nbins
      .self$constant_breaks()
    },
    copy = function() {
      "A copy of the object with the same masses, nbins and breaks."
      copied <- Cdsw(.self$masses, nbins = .self$nbins)
      copied$breaks <- .self$breaks
      copied
    },
    constant_breaks = function(digits = 2) {
      "Equal-width windows; returns the breaks invisibly."
      mass_range <- .padded_range(.self$masses, digits)
      .self$breaks <- round(seq(mass_range[1], mass_range[2], length.out = .self$nbins + 1), digits = digits)
      invisible(.self$breaks)
    },
    quantile_breaks = function(digits = 2) {
      "Windows with equal numbers of MS1 masses; widths may violate the limits of sampling_breaks()."
      .self$breaks <- .adjust_breaks(.quantile_breaks(.self$masses, .self$nbins), digits = digits)
      invisible(.self$breaks)
    },
    sampling_breaks = function(maxwindow = 150, minwindow = 5, digits = 2) {
      "Quantile breaks of the masses mixed with uniformly spread data, with the smallest share of uniform data
      that keeps all widths within [minwindow, maxwindow]. Stops if even equal-width windows violate the limits;
      returns the breaks invisibly."
      equal <- seq(min(.self$masses), max(.self$masses), length.out = .self$nbins + 1)
      if (!.hard_constrain(equal, minwindow = minwindow, maxwindow = maxwindow)) {
        stop(
          "No window layout satisfies the limits: ",
          .self$nbins,
          " windows of equal width are ",
          round(diff(equal)[1], digits = 2),
          " wide, outside [minwindow, maxwindow] = [",
          minwindow,
          ", ",
          maxwindow,
          "]."
        )
      }
      mixed <- .sampling_breaks(.self$masses, .self$nbins, equal, minwindow = minwindow, maxwindow = maxwindow)
      .self$breaks <- .adjust_breaks(mixed, digits = digits)
      invisible(.self$breaks)
    },
    plot_breaks = function() {
      "ggplot of the current breaks against equal-width and quantile breaks."
      index <- seq(0, .self$nbins)
      layouts <- c("equal width", "quantile", "current")
      breaks_long <- data.frame(
        mz = c(
          seq(min(.self$masses), max(.self$masses), length.out = .self$nbins + 1),
          .quantile_breaks(.self$masses, .self$nbins),
          .self$breaks
        ),
        index = rep(index, 3),
        layout = factor(rep(layouts, each = .self$nbins + 1), levels = layouts)
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
      "ggplot of the precursors per window (from < m/z <= to); the dashed line marks equal counts."
      from <- utils::head(.self$breaks, -1)
      to <- utils::tail(.self$breaks, -1)
      windows <- data.frame(from = from, to = to, count = .window_counts(.self$masses, from, to))
      ggplot2::ggplot(windows) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin = .data$from, xmax = .data$to, ymin = 0, ymax = .data$count),
          fill = "#56B4E9",
          colour = "white"
        ) +
        ggplot2::geom_hline(yintercept = length(.self$masses) / .self$nbins, linetype = "dashed") +
        ggplot2::labs(x = "precursor m/z", y = "precursors per window")
    },
    as_table = function(overlap = 1) {
      "Windows as a data.frame with start (from), end (to), center (mid), width and the number of precursors with
      from < m/z <= to (counts); overlap in m/z."
      from <- utils::head(.self$breaks, -1) - overlap * 0.5
      to <- utils::tail(.self$breaks, -1) + overlap * 0.5
      width <- to - from
      counts <- .window_counts(.self$masses, from, to)
      data.frame(from = from, to = to, mid = from + width * 0.5, width = width, counts = counts)
    },
    error = function() {
      "Deviation of the windows from equal counts: a list, see objective_ms1()."
      objective_ms1(.self$breaks, .self$masses)
    },
    optimize_windows = function(digits = 1, maxbin = 15, overlap = 0) {
      "Move window starts and ends to regions with few precursor masses; returns a data.frame."
      data.frame(readjust_windows(.self$as_table(overlap = overlap), .self$masses, digits = digits, maxbin = maxbin))
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
