# Range of x rounded outward to multiples of 1 / scale and padded by 1 / scale, so that it encloses x
.padded_range <- function(x, scale) {
  c(floor(min(x) * scale) - 1, ceiling(max(x) * scale) + 1) / scale
}

# Breaks rounded to `digits`; the outer breaks are padded so that they enclose the masses
.adjust_breaks <- function(breaks, digits) {
  padded <- .padded_range(breaks, 10^digits)
  breaks <- round(breaks, digits = digits)
  breaks[c(1, length(breaks))] <- padded
  breaks
}


#' Deviation of a window layout from equal numbers of MS1 masses per window
#'
#' @param splits The window boundaries.
#' @param data The MS1 masses.
#' @return A list with `manhattan` (Manhattan distance), `euclidean` (Euclidean distance), `normalized_counts`
#'   (observed counts) and `equal_counts` (the counts of windows with equal numbers of masses), both counts
#'   normalized.
#' @noRd
.objective_ms1 <- function(splits, data) {
  counts <- .window_counts(data, utils::head(splits, -1), utils::tail(splits, -1))
  normalized_counts <- counts / sqrt(sum(counts^2))
  equal_counts <- rep(1 / sqrt(length(counts)), length(counts))
  list(
    manhattan = mean(1 - normalized_counts / max(normalized_counts)),
    euclidean = sqrt(sum((normalized_counts - equal_counts)^2)) / length(counts),
    normalized_counts = normalized_counts,
    equal_counts = equal_counts
  )
}

# TRUE if all window widths lie within [min_window, max_window]
.hard_constrain <- function(splits, min_window, max_window) {
  widths <- diff(splits)
  all(widths >= min_window & widths <= max_window)
}

# Number of masses x with from < x <= to, for each window (the intervals of hist())
.window_counts <- function(x, from, to) {
  sorted <- sort(x)
  findInterval(to, sorted) - findInterval(from, sorted)
}

# Unrounded breaks of n_bins windows with equal numbers of masses
.mass_quantiles <- function(x, n_bins) {
  stats::quantile(x, seq(0, 1, length.out = n_bins + 1), names = FALSE)
}

# Breaks at the quantiles of the mixture (1 - share) * masses + share * uniform over the mass range: the quantile breaks
# of the masses plus equally spaced points that make up `share` of the data, in the limit of infinitely many points.
# `sorted` holds the sorted masses; 0 < share <= 1.
.mixture_quantiles <- function(sorted, n_bins, share) {
  n <- length(sorted)
  span <- sorted[n] - sorted[1]
  # the mixture CDF jumps by (1 - share) / n at each mass and rises with slope share / span in between
  at <- c(sorted[1], sorted)
  cdf <- c(0, (1 - share) * seq_len(n) / n + share * (sorted - sorted[1]) / span)
  probs <- seq_len(n_bins - 1) / n_bins
  k <- findInterval(probs, cdf)
  inner <- pmin(at[k] + (probs - cdf[k]) * span / share, at[k + 1])
  c(sorted[1], inner, sorted[n])
}

# Breaks rounded to `digits` with all widths within [min_window, max_window]: the quantile breaks if they fit,
# otherwise the mixture breaks with the smallest share of uniform data that fits, found by bisection to 1e-9. The
# limits are checked on the rounded, padded breaks. `equal` (share 1, rounded) must fit; the search keeps the last
# layout that fit, so it returns one even if fitting is not monotone in the share.
.fitting_mixture_breaks <- function(x, n_bins, equal, min_window, max_window, digits) {
  quantiles <- .adjust_breaks(.mass_quantiles(x, n_bins), digits = digits)
  if (.hard_constrain(quantiles, min_window = min_window, max_window = max_window)) {
    return(quantiles)
  }
  sorted <- sort(x)
  fitting <- equal
  lower <- 0
  upper <- 1
  while (upper - lower > 1e-9) {
    share <- (lower + upper) / 2
    breaks <- .adjust_breaks(.mixture_quantiles(sorted, n_bins, share), digits = digits)
    if (.hard_constrain(breaks, min_window = min_window, max_window = max_window)) {
      fitting <- breaks
      upper <- share
    } else {
      lower <- share
    }
  }
  fitting
}

# Breaks of n_bins windows of equal width, rounded to `digits`, enclosing the masses
.constant_breaks <- function(masses, n_bins, digits) {
  mass_range <- .padded_range(masses, 10^digits)
  round(seq(mass_range[1], mass_range[2], length.out = n_bins + 1), digits = digits)
}

# Breaks of n_bins windows with equal numbers of masses, rounded to `digits`, enclosing the masses
.quantile_breaks <- function(masses, n_bins, digits) {
  .adjust_breaks(.mass_quantiles(masses, n_bins), digits = digits)
}

# Quantile breaks of the masses mixed with uniformly spread data, rounded to `digits`, with the smallest share of
# uniform data that keeps all widths of the rounded breaks within [min_window, max_window]. Stops if even equal-width
# windows violate the limits.
.mixture_breaks <- function(masses, n_bins, max_window, min_window, digits) {
  equal <- .adjust_breaks(seq(min(masses), max(masses), length.out = n_bins + 1), digits = digits)
  if (!.hard_constrain(equal, min_window = min_window, max_window = max_window)) {
    widths <- round(range(diff(equal)), digits = digits)
    stop(
      "No window layout satisfies the limits: ",
      n_bins,
      " windows of equal width, with the breaks rounded to `digits`, are ",
      widths[1],
      " to ",
      widths[2],
      " wide, not all within [`min_window`, `max_window`] = [",
      min_window,
      ", ",
      max_window,
      "].",
      call. = FALSE
    )
  }
  .fitting_mixture_breaks(masses, n_bins, equal, min_window = min_window, max_window = max_window, digits = digits)
}


#' Compute dynamic SWATH (DIA) windows
#'
#' @description
#' An R6 class that splits the MS1 mass range into `n_bins` windows with about equal numbers of precursors. The window
#' boundaries (breaks) are computed once, when the object is created, with one of three methods:
#' * `"mixture"` (default): the quantile breaks of the masses mixed with uniformly spread data, with the smallest share
#'   of uniform data that keeps all window widths within `[min_window, max_window]`; the quantile breaks themselves if
#'   they fit. The limits hold for the rounded and padded breaks. Stops if even windows of equal width violate the
#'   limits.
#' * `"quantile"`: windows with equal numbers of masses, of any width.
#' * `"constant"`: windows of equal width.
#'
#' The breaks are rounded to `digits` decimals, and the outer breaks are rounded outward and padded by `1 / 10^digits`,
#' so that the windows enclose all masses.
#'
#' The object cannot be changed: its fields are read-only and its methods only read them, so their results do not
#' depend on the order in which they are called. To compare methods or settings, create one object for each.
#' `clone()` copies an object.
#'
#' @return `Cdsw$new()` returns a `Cdsw` R6 object; see Methods.
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw$new(masses, max_window = 100)
#' cdsw
#' cdsw$plot_breaks()
#' cdsw$plot()
#' head(cdsw$as_table())
#' head(cdsw$optimize_windows())
#' cdsw$plot_cycle()
#' # one object per method
#' quantile_windows <- Cdsw$new(masses, method = "quantile")
#' c(mixture = cdsw$error()$manhattan, quantile = quantile_windows$error()$manhattan)
Cdsw <- R6::R6Class(
  "Cdsw",
  public = list(
    #' @description
    #' Compute the breaks of `n_bins` windows over the masses with `method`.
    #' @param masses MS1 masses (m/z): a numeric vector of finite values.
    #' @param n_bins Number of windows.
    #' @param method How to place the breaks: `"mixture"` (default), `"quantile"` or `"constant"`; see Description.
    #' @param max_window,min_window Largest and smallest window width in m/z; used by `method = "mixture"` only.
    #' @param digits Number of decimal places of the breaks.
    #' @return A new `Cdsw` object.
    initialize = function(
      masses,
      n_bins = 25,
      method = c("mixture", "quantile", "constant"),
      max_window = 150,
      min_window = 5,
      digits = 2
    ) {
      method <- .match_method(method, c("mixture", "quantile", "constant"), arg = "method")
      .check_masses(masses)
      .check_number(n_bins, "n_bins", positive = TRUE, whole = TRUE)
      .check_number(max_window, "max_window", positive = TRUE, whole = FALSE)
      .check_number(min_window, "min_window", positive = FALSE, whole = FALSE)
      .check_number(digits, "digits", positive = FALSE, whole = TRUE)
      private$.masses <- masses
      private$.n_bins <- n_bins
      private$.method <- method
      private$.digits <- digits
      private$.breaks <- switch(
        method,
        mixture = .mixture_breaks(masses, n_bins, max_window = max_window, min_window = min_window, digits = digits),
        quantile = .quantile_breaks(masses, n_bins, digits = digits),
        constant = .constant_breaks(masses, n_bins, digits = digits)
      )
    },
    #' @description
    #' Print the number of masses and windows, the method and the range of window widths.
    #' @param ... Not used.
    #' @return The object, invisibly.
    print = function(...) {
      widths <- formatC(range(diff(private$.breaks)), format = "f", digits = private$.digits)
      n_masses <- format(length(private$.masses), big.mark = ",")
      cat(sprintf("<Cdsw> %s masses in %s windows, method \"%s\"\n", n_masses, private$.n_bins, private$.method))
      cat(sprintf("window widths %s to %s m/z\n", widths[1], widths[2]))
      invisible(self)
    },
    #' @description
    #' The windows as a data frame.
    #' @param overlap Overlap of neighboring windows in m/z, at least 0: each window is widened by `overlap / 2` on both
    #'   sides.
    #' @return A data frame with one row per window: start `from`, end `to`, center `mid`, `width` and the number of
    #'   precursors with `from < m/z <= to` (`n_precursors`).
    as_table = function(overlap = 1) {
      .check_number(overlap, "overlap", positive = FALSE, whole = FALSE)
      from <- utils::head(private$.breaks, -1) - overlap * 0.5
      to <- utils::tail(private$.breaks, -1) + overlap * 0.5
      width <- to - from
      n_precursors <- .window_counts(private$.masses, from, to)
      data.frame(from = from, to = to, mid = from + width * 0.5, width = width, n_precursors = n_precursors)
    },
    #' @description
    #' Deviation of the windows from equal numbers of precursors.
    #' @return A list with `manhattan` (Manhattan distance), `euclidean` (Euclidean distance), `normalized_counts`
    #'   (precursors per window) and `equal_counts` (the counts of windows with equal numbers of precursors), both
    #'   counts normalized to length 1.
    error = function() {
      .objective_ms1(private$.breaks, private$.masses)
    },
    #' @description
    #' Move window starts and ends to regions with few precursors, with [readjust_windows()].
    #' @param bin_width Width of the bins (m/z) into which the masses are binned.
    #' @param min_count Bins with fewer than `min_count` precursors count as empty.
    #' @param overlap Overlap of neighboring windows in m/z, as in `as_table()`.
    #' @return A data frame with the columns of `as_table()`.
    optimize_windows = function(bin_width = 0.1, min_count = 15, overlap = 1) {
      readjust_windows(self$as_table(overlap = overlap), private$.masses, bin_width = bin_width, min_count = min_count)
    },
    #' @description
    #' Plot the precursors per window (`from < m/z <= to`); the dashed line marks equal counts.
    #' @return A ggplot object.
    plot = function() {
      ggplot2::ggplot(self$as_table(overlap = 0)) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin = .data$from, xmax = .data$to, ymin = 0, ymax = .data$n_precursors),
          fill = "#56B4E9",
          colour = "white"
        ) +
        ggplot2::geom_hline(yintercept = length(private$.masses) / private$.n_bins, linetype = "dashed") +
        ggplot2::labs(x = "precursor m/z", y = "precursors per window")
    },
    #' @description
    #' Plot the breaks of this object against the constant and quantile breaks of the same masses.
    #' @return A ggplot object.
    plot_breaks = function() {
      masses <- private$.masses
      n_bins <- private$.n_bins
      layouts <- c("constant", "quantile", "this object")
      breaks_long <- data.frame(
        mz = c(
          .constant_breaks(masses, n_bins, digits = private$.digits),
          .quantile_breaks(masses, n_bins, digits = private$.digits),
          private$.breaks
        ),
        index = rep(seq(0, n_bins), 3),
        layout = factor(rep(layouts, each = n_bins + 1), levels = layouts)
      )
      ggplot2::ggplot(breaks_long, ggplot2::aes(x = .data$mz, y = .data$index, colour = .data$layout)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 1) +
        ggplot2::scale_colour_manual(
          values = c("constant" = "#999999", "quantile" = "#E69F00", "this object" = "#0072B2"),
          name = "breaks"
        ) +
        ggplot2::labs(x = "precursor m/z", y = "break index")
    },
    #' @description
    #' Plot the DIA cycle, one window per step.
    #' @param overlap Overlap of neighboring windows in m/z, as in `as_table()`.
    #' @return A ggplot object.
    plot_cycle = function(overlap = 1) {
      windows <- self$as_table(overlap = overlap)
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
  ),
  active = list(
    #' @field masses MS1 masses (read-only).
    masses = function(value) .read_only(value, "masses", "Cdsw", private$.masses),
    #' @field n_bins Number of windows (read-only).
    n_bins = function(value) .read_only(value, "n_bins", "Cdsw", private$.n_bins),
    #' @field method Method that placed the breaks (read-only).
    method = function(value) .read_only(value, "method", "Cdsw", private$.method),
    #' @field breaks Window boundaries, `n_bins + 1` values (read-only).
    breaks = function(value) .read_only(value, "breaks", "Cdsw", private$.breaks)
  ),
  private = list(
    .masses = NULL,
    .n_bins = NULL,
    .method = NULL,
    .digits = NULL,
    .breaks = NULL
  )
)
