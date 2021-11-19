# moves the windows start and end to regions where no peaks are observed
.makenewfromto <- function(windfrom, empty , isfrom = TRUE) {
  newfrom <- NULL
  for (from in windfrom) {
    idx <- which.min(abs(from - empty))
    startmass <- 0
    if (isfrom) {
      if (idx > 1) {
        if (empty[idx] < from) {
          startmass <- empty[idx]
        } else {
          startmass <- empty[idx - 1]
        }
      } else{
        startmass <- from
      }
    } else{
      if (idx < length(empty)) {
        if (empty[idx] > from) {
          startmass <- empty[idx]
        } else {
          startmass <- empty[idx + 1]
        }
      } else{
        startmass <- from
      }
    }
    newfrom <- c(newfrom, startmass)
  }
  return(newfrom)
}

#' Readjust windows so that boundaries in regions of few peaks.
#'
#' @param wind a data frame with columns from and to
#' @param ms1data masses
#' @param digits mass accuracy
#' @param maxbin maximum number of bins
#' @param plot diagnostic plots (default FALSE)
#' @export
#' @examples
#' data(masses)
#' cdsw <- Cdsw(masses)
#' breaks <- cdsw$sampling_breaks(maxwindow=100,plot=TRUE)
#' table <- cdsw$asTable()
#' dim(table)
#' head(table)
#'
#' tmp <- readjustWindows(table, masses,maxbin=10)
#' data.frame(tmp)
#'
readjustWindows <-
  function(wind ,
           ms1data,
           digits = 1,
           maxbin = 15,
           plot = FALSE) {
    breaks <-
      seq(
        round(min(ms1data) - 1 / 10 ^ digits, digits = 1),
        round(max(ms1data) + 1 / 10 ^ digits, digits = 1),
        by = 1 / 10 ^ digits
      )

    reshist <- graphics::hist(ms1data, breaks = breaks, plot = plot)
    if (plot) {
      graphics::abline(v = wind$from,
                       col = 2,
                       lty = 2)
      graphics::abline(v = wind$to,
                       col = 3,
                       lty = 2)
    }
    empty <- reshist$mids[which(reshist$counts < maxbin)]
    newfrom <- .makenewfromto(wind$from , empty)
    newto <- .makenewfromto(wind$to , empty , isfrom = FALSE)
    if (plot) {
      graphics::plot(reshist, xlim = c(newfrom[round(length(newfrom) / 2)] - 3, newfrom[round(length(newfrom) /
                                                                                                2)] + 3))
      graphics::abline(v = empty, col = "gray")

      graphics::abline(v = newfrom, lwd = 0.5, col = "red")
      graphics::abline(v = newto , lwd = 0.5, col = "green")
      graphics::plot(reshist, xlim = c(newfrom[round(length(newfrom) / 4)] -
                                         3, newfrom[round(length(newfrom) / 4)] + 3))
      graphics::abline(v = newfrom, lwd = 0.5, col = "red")
      graphics::abline(v = newto , lwd = 0.5, col = "green")
    }
    width <- (newto - newfrom)
    mid <- (newfrom + newto) * 0.5
    newCounts <- NULL

    for (i in 1:length(newfrom))
    {
      newCounts <-
        c(newCounts, sum(ms1data >= newfrom[i] & ms1data <= newto[i]))
    }
    list(
      from = newfrom,
      to = newto,
      mid = mid,
      width = width,
      counts = newCounts
    )
  }
