.adjustBreaks <- function(breaks, digits) {
  breaks[1] <- breaks[1] - 1 / (10 ^ digits)
  breaks[length(breaks)] <-
    breaks[length(breaks)] + 1 / (10 ^ digits)
  round(breaks, digits = digits)
}

#' MS masses
#' A dataset containing approx 150000 MS1 precursor masses
#'
"masses"

#' compute the deviation from optimum: equal number of MS1 per bin
#' @param splits the new window boundaries
#' @param data the data
#' @return list with score1 - manhattan distance, score2 - euclidean distance, counts - observed counts, optimumN - optimum counts
#'
objectiveMS1Function <- function(splits, data) {
  counts <- graphics::hist(data, breaks = splits, plot = FALSE)$counts
  nbins <- length(splits) - 1
  optimumN <- length(data) / (length(splits) - 1)
  optimumN <- rep(optimumN, nbins)
  counts <- counts / sqrt(sum((counts) ^ 2))
  optimumN <- optimumN / sqrt(sum(optimumN ^ 2))
  score2 <- sqrt(sum((counts - optimumN) ^ 2)) / length(counts)
  score1 <-
    sum(abs(counts / max(abs(counts)) - optimumN / max(abs(optimumN)))) / length(counts)
  return(list(
    score1 = score1,
    score2 = score2,
    counts = counts,
    optimumN = (optimumN)
  ))
}

#' tests hard constraints
#' @param splits the proposed splits
#' @param minwindow min window size
#' @param maxwindow max window size
#' @export
hardconstrain <- function(splits,
                          minwindow = 5,
                          maxwindow = 50) {
  difsp <- diff(splits)
  return(sum(difsp >= minwindow) == length(difsp) &
           sum(difsp <= maxwindow) == length(difsp))
}
#' Compute dynamic swath windows
#' @field masses MS1 masses
#' @field breaks the breaks
#' @field nbins number of bins
#' @field digits mass accuracy in result
#' @import methods
#' @importFrom methods new
#' @include readjustWindow.R
#' @export Cdsw
#' @exportClass Cdsw
#' @examples
#' data(masses)
#' cdsw <- Cdsw(masses)
#' cdsw$sampling_breaks(maxwindow=100,plot=TRUE)
#' cdsw$plot()
#' cdsw$asTable()
#' cdsw$breaks
Cdsw <- setRefClass(
  "Cdsw",
  fields = list(
    masses = "numeric",
    breaks = "numeric",
    nbins = "numeric",
    digits = "numeric"
  )
  ,
  methods = list(
    initialize = function(masses,
                          nbins = 25,
                          digits = 1) {
      .self$masses = masses
      .self$nbins = nbins
      .self$digits = digits
      constant_breaks()
    },
    constant_breaks = function(digits = 2) {
      minmass <- round(min(.self$masses) - 1 / 10 ^ digits, digits = 2)
      maxmass <-
        round(max(.self$masses) + 1 / 10 ^ digits, digits = 2)
      .self$breaks <-
        round(seq(minmass, maxmass, length.out = nbins + 1), digits = digits)
      invisible(.self$breaks)
    }
    ,
    quantile_breaks = function(digits = 2) {
      "same number of MS1 in each window but might violate hard constraints"
      qqs <-
        quantile(.self$masses, seq(0, 1, length = nbins + 1))
      .self$breaks <-
        .adjustBreaks(qqs, digits = digits)
      invisible(.self$breaks)
    },
    sampling_breaks = function(maxwindow = 150,
                               minwindow = 5,
                               digits = 2,
                               plot = FALSE)
    {
      "starts with quantile breaks but mixes with uniform data to satisfy had constraints"
      nrbreaks <- nbins + 1
      qqs <-
        quantile(masses, probs = seq(0, 1, by = 1 / (nbins)))
      unif <-
        seq(min(masses), max(masses), length = (nrbreaks))

      if (plot) {
        graphics::plot(qqs, 1:nrbreaks, type = "b")
        graphics::legend("topleft", legend = c(paste("maxwindow = ", maxwindow),
                                               paste("nbins = ", nrbreaks)))
        # equidistant spaced bins
        graphics::lines(unif, 1:nrbreaks, col = 2, type =
                          "b")

      }

      if (!hardconstrain(unif, minwindow = minwindow, maxwindow)) {
        warning(
          "there is no way to generate bins given this number of bins " ,
          nbins,
          "and minwindow size :",
          minwindow,
          " , maxwindow size ",
          maxwindow,
          "\n"
        )
      }

      mixeddata <- masses

      while (!hardconstrain(qqs, minwindow, maxwindow)) {
        uniformdata <-
          runif(round(length(masses) / 20),
                min = min(masses),
                max = max(masses))
        mixeddata <- c(mixeddata, uniformdata)
        qqs <-
          quantile(mixeddata, probs = seq(0, 1, by = 1 / (nbins)))
        if (plot) {
          graphics::lines(qqs, 1:nrbreaks, type = "b", col = "#00DD00AA")
        }
      }
      .self$breaks <-
        .adjustBreaks(qqs, digits = digits)
      invisible(.self$breaks)

    },
    plot = function() {
      tmp <- graphics::hist(x = .self$masses,
                            breaks = .self$breaks,
                            plot = FALSE)
      names(tmp$counts) <-
        round(tmp$mids, digits = 2)
      barplot(tmp$counts, las = 2)
    },
    asTable = function(overlap = 1) {
      "make windows"
      q <- .self$breaks
      n <- length(q) - 1
      idx <- 1:n
      from <- q[idx] - overlap * 0.5
      to <- q[idx + 1] + overlap * 0.5
      width <-  (to - from)
      mid <- from + width * 0.5

      tmp <- data.frame(from, to, mid, width)
      counts <- vector("integer")
      for (i in 1:nrow(tmp)) {
        counts[i] <-
          sum(.self$masses > tmp$from[i] & .self$masses < tmp$to[i])
      }
      tmp$counts <- counts
      return(tmp)
    },
    show = function() {
      "get a table with windows start and end"
      print(asTable())
    },
    error = function() {
      "show error"
      objectiveMS1Function(.self$breaks, .self$masses)
    },
    optimizeWindows = function(digits = 1,
                               maxbin = 15,
                               plot = FALSE,
                               overlap = 0.) {
      "optimizes the windows"
      data.frame(readjustWindows(
        asTable(overlap = overlap) ,
        masses,
        digits = digits,
        maxbin = maxbin,
        plot = plot
      ))
    },
    showCycle = function(overlap = 1) {
      "shows the generated DIA cycle"
      tmp <- .self$asTable(overlap = overlap)
      graphics::plot(
        c(1, nrow(tmp) + 1),
        c(min(tmp$from), max(tmp$to)),
        type = "n",
        xlab = "window #",
        ylab = "mz",
        main = "DIA Cycle"
      )
      rect(1:nrow(tmp), tmp$from, 1:nrow(tmp) + 1, tmp$to)
    }
  )
)
