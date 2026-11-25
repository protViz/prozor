# Range of x rounded outward to multiples of 1 / scale and padded by 1 / scale, so that it encloses x
.padded_range <- function(x, scale) {
  c(floor(min(x) * scale) - 1, ceiling(max(x) * scale) + 1) / scale
}

# Breaks rounded to `digits`; the outer breaks are padded so that they enclose the masses
.adjust_breaks <- function(breaks, digits) {
  outer <- .padded_range(breaks, 10^digits)
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

# Unrounded breaks of nbins windows with equal numbers of masses
.mass_quantiles <- function(x, nbins) {
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

# Unrounded breaks with all widths within [minwindow, maxwindow]: the quantile breaks if they fit, otherwise the mixture
# breaks with the smallest share of uniform data that fits, found by bisection to 1e-9. `equal` (share 1) must fit; the
# search keeps the last layout that fit, so it returns one even if fitting is not monotone in the share.
.fitting_mixture_breaks <- function(x, nbins, equal, minwindow, maxwindow) {
  quantiles <- .mass_quantiles(x, nbins)
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

# Breaks of nbins windows of equal width, rounded to `digits`, enclosing the masses
.constant_breaks <- function(masses, nbins, digits) {
  mass_range <- .padded_range(masses, 10^digits)
  round(seq(mass_range[1], mass_range[2], length.out = nbins + 1), digits = digits)
}

# Breaks of nbins windows with equal numbers of masses, rounded to `digits`, enclosing the masses
.quantile_breaks <- function(masses, nbins, digits) {
  .adjust_breaks(.mass_quantiles(masses, nbins), digits = digits)
}

# Quantile breaks of the masses mixed with uniformly spread data, with the smallest share of uniform data that keeps
# all widths within [minwindow, maxwindow], rounded to `digits`. Stops if even equal-width windows violate the limits.
.sampling_breaks <- function(masses, nbins, maxwindow, minwindow, digits) {
  equal <- seq(min(masses), max(masses), length.out = nbins + 1)
  if (!.hard_constrain(equal, minwindow = minwindow, maxwindow = maxwindow)) {
    stop(
      "No window layout satisfies the limits: ",
      nbins,
      " windows of equal width are ",
      round(diff(equal)[1], digits = 2),
      " wide, outside [minwindow, maxwindow] = [",
      minwindow,
      ", ",
      maxwindow,
      "].",
      call. = FALSE
    )
  }
  mixed <- .fitting_mixture_breaks(masses, nbins, equal, minwindow = minwindow, maxwindow = maxwindow)
  .adjust_breaks(mixed, digits = digits)
}

# Stops unless `value` is a single number > 0 (if `positive`) or >= 0, and a whole number if `whole`; `name` is the
# argument name for the message
.check_number <- function(value, name, positive, whole) {
  scalar <- is.numeric(value) && length(value) == 1 && !is.na(value)
  valid <- scalar && all(value >= 0, value > 0 || !positive, value == round(value) || !whole)
  if (!valid) {
    kind <- paste0(if (positive) "positive" else "non-negative", if (whole) " whole" else "")
    stop("`", name, "` must be a single ", kind, " number, not ", deparse1(value), ".", call. = FALSE)
  }
}

# Stops unless `masses` is a non-empty numeric vector of finite values
.check_masses <- function(masses) {
  if (!is.numeric(masses) || length(masses) == 0) {
    stop(
      "`masses` must be a non-empty numeric vector, not ",
      class(masses)[1],
      " of length ",
      length(masses),
      ".",
      call. = FALSE
    )
  }
  if (!all(is.finite(masses))) {
    stop("`masses` must not contain NA, NaN or infinite values; found ", sum(!is.finite(masses)), ".", call. = FALSE)
  }
}

# The `method` argument: the first of `choices` for the default (all choices), otherwise exactly one of them
.match_method <- function(method, choices) {
  if (identical(method, choices)) {
    return(choices[1])
  }
  if (!(length(method) == 1 && method %in% choices)) {
    stop(
      "`method` must be one of ",
      paste0("\"", choices, "\"", collapse = ", "),
      "; not ",
      deparse1(method),
      ".",
      call. = FALSE
    )
  }
  method
}

# Stops for an assignment to the read-only field `name` of an object of the R6 class `class`
.stop_read_only <- function(name, class) {
  stop("`", name, "` of a ", class, " object is read-only; create a new object with ", class, "$new().", call. = FALSE)
}


#' Compute dynamic SWATH (DIA) windows
#'
#' @description
#' An R6 class that splits the MS1 mass range into `nbins` windows with about equal numbers of precursors. The window
#' boundaries (breaks) are computed once, when the object is created, with one of three methods:
#' * `"sampling"` (default): the quantile breaks of the masses mixed with uniformly spread data, with the smallest share
#'   of uniform data that keeps all window widths within `[minwindow, maxwindow]`; the quantile breaks themselves if
#'   they fit. Stops if even windows of equal width violate the limits.
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
#' @family DIA windows
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw$new(masses, maxwindow = 100)
#' cdsw
#' cdsw$plot_breaks()
#' cdsw$plot()
#' head(cdsw$as_table())
#' head(cdsw$optimize_windows())
#' cdsw$show_cycle()
#' # one object per method
#' quantile_windows <- Cdsw$new(masses, method = "quantile")
#' c(sampling = cdsw$error()$score1, quantile = quantile_windows$error()$score1)
Cdsw <- R6::R6Class(
  "Cdsw",
  public = list(
    #' @description
    #' Compute the breaks of `nbins` windows over the masses with `method`.
    #' @param masses MS1 masses (m/z): a numeric vector of finite values.
    #' @param nbins Number of windows.
    #' @param method How to place the breaks: `"sampling"` (default), `"quantile"` or `"constant"`; see Description.
    #' @param maxwindow,minwindow Largest and smallest window width in m/z; used by `method = "sampling"` only.
    #' @param digits Number of decimal places of the breaks.
    #' @return A new `Cdsw` object.
    initialize = function(
      masses,
      nbins = 25,
      method = c("sampling", "quantile", "constant"),
      maxwindow = 150,
      minwindow = 5,
      digits = 2
    ) {
      method <- .match_method(method, c("sampling", "quantile", "constant"))
      .check_masses(masses)
      .check_number(nbins, "nbins", positive = TRUE, whole = TRUE)
      .check_number(maxwindow, "maxwindow", positive = TRUE, whole = FALSE)
      .check_number(minwindow, "minwindow", positive = FALSE, whole = FALSE)
      .check_number(digits, "digits", positive = FALSE, whole = TRUE)
      private$.masses <- masses
      private$.nbins <- nbins
      private$.method <- method
      private$.digits <- digits
      private$.breaks <- switch(
        method,
        sampling = .sampling_breaks(masses, nbins, maxwindow = maxwindow, minwindow = minwindow, digits = digits),
        quantile = .quantile_breaks(masses, nbins, digits = digits),
        constant = .constant_breaks(masses, nbins, digits = digits)
      )
    },
    #' @description
    #' Print the number of masses and windows, the method and the range of window widths.
    #' @param ... Not used.
    #' @return The object, invisibly.
    print = function(...) {
      widths <- formatC(range(diff(private$.breaks)), format = "f", digits = private$.digits)
      cat(
        "<Cdsw> ",
        format(length(private$.masses), big.mark = ","),
        " masses in ",
        private$.nbins,
        " windows, method \"",
        private$.method,
        "\"\n",
        "window widths ",
        widths[1],
        " to ",
        widths[2],
        " m/z\n",
        sep = ""
      )
      invisible(self)
    },
    #' @description
    #' The windows as a data frame.
    #' @param overlap Overlap of neighboring windows in m/z: each window is widened by `overlap / 2` on both sides.
    #' @return A data frame with one row per window: start `from`, end `to`, center `mid`, `width` and the number of
    #'   precursors with `from < m/z <= to` (`counts`).
    as_table = function(overlap = 1) {
      from <- utils::head(private$.breaks, -1) - overlap * 0.5
      to <- utils::tail(private$.breaks, -1) + overlap * 0.5
      width <- to - from
      counts <- .window_counts(private$.masses, from, to)
      data.frame(from = from, to = to, mid = from + width * 0.5, width = width, counts = counts)
    },
    #' @description
    #' Deviation of the windows from equal numbers of precursors.
    #' @return A list with `score1` (Manhattan distance), `score2` (Euclidean distance), `counts` (precursors per
    #'   window) and `optimum_n` (equal counts), both counts normalized to length 1.
    error = function() {
      objective_ms1(private$.breaks, private$.masses)
    },
    #' @description
    #' Move window starts and ends to regions with few precursors, with [readjust_windows()].
    #' @param bin_width Width of the bins (m/z) into which the masses are binned.
    #' @param maxbin Bins with fewer than `maxbin` precursors count as empty.
    #' @param overlap Overlap of neighboring windows in m/z, as in `as_table()`.
    #' @return A data frame with the columns of `as_table()`.
    optimize_windows = function(bin_width = 0.1, maxbin = 15, overlap = 0) {
      readjust_windows(self$as_table(overlap = overlap), private$.masses, bin_width = bin_width, maxbin = maxbin)
    },
    #' @description
    #' Plot the precursors per window (`from < m/z <= to`); the dashed line marks equal counts.
    #' @return A ggplot object.
    plot = function() {
      from <- utils::head(private$.breaks, -1)
      to <- utils::tail(private$.breaks, -1)
      windows <- data.frame(from = from, to = to, count = .window_counts(private$.masses, from, to))
      ggplot2::ggplot(windows) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin = .data$from, xmax = .data$to, ymin = 0, ymax = .data$count),
          fill = "#56B4E9",
          colour = "white"
        ) +
        ggplot2::geom_hline(yintercept = length(private$.masses) / private$.nbins, linetype = "dashed") +
        ggplot2::labs(x = "precursor m/z", y = "precursors per window")
    },
    #' @description
    #' Plot the breaks of this object against the constant and quantile breaks of the same masses.
    #' @return A ggplot object.
    plot_breaks = function() {
      masses <- private$.masses
      nbins <- private$.nbins
      layouts <- c("constant", "quantile", "this object")
      breaks_long <- data.frame(
        mz = c(
          .constant_breaks(masses, nbins, digits = private$.digits),
          .quantile_breaks(masses, nbins, digits = private$.digits),
          private$.breaks
        ),
        index = rep(seq(0, nbins), 3),
        layout = factor(rep(layouts, each = nbins + 1), levels = layouts)
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
    show_cycle = function(overlap = 1) {
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
    masses = function(value) {
      if (!missing(value)) {
        .stop_read_only("masses", "Cdsw")
      }
      private$.masses
    },
    #' @field nbins Number of windows (read-only).
    nbins = function(value) {
      if (!missing(value)) {
        .stop_read_only("nbins", "Cdsw")
      }
      private$.nbins
    },
    #' @field method Method that placed the breaks (read-only).
    method = function(value) {
      if (!missing(value)) {
        .stop_read_only("method", "Cdsw")
      }
      private$.method
    },
    #' @field breaks Window boundaries, `nbins + 1` values (read-only).
    breaks = function(value) {
      if (!missing(value)) {
        .stop_read_only("breaks", "Cdsw")
      }
      private$.breaks
    }
  ),
  private = list(
    .masses = NULL,
    .nbins = NULL,
    .method = NULL,
    .digits = NULL,
    .breaks = NULL
  )
)
