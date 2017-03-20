#' compute FDR given a score
#' @param score a vector with scores
#' @param proteinID - list with protein id's
#' @param decoy decoy pattern, default "REV_"
#' @param decreasing If the scores sorted from small (better),
#'  than decreasing=FALSE, or if from large (better),
#'  than decreasing TRUE (default).
#' @return list with proteinID, decoy_hit (indicates if decoy), score the search engine score,
#' FDR1 false discovery rate estimated using the method of Gyggi, FDR2 - estimated using the method of Kell.
#'
computeFDR <-function(score, proteinID, decoy = "REV_", decreasing =TRUE)
{
    idx <- order(score, decreasing = TRUE)
    score <- score[idx]
    decoy_hit <- grepl(decoy, proteinID[idx])

    FP <- cumsum(decoy_hit)
    TP <- 1:length(idx) - FP

    FDR1 <- (2 * FP) / (TP + FP)
    FDR2 <- FP / TP
    return(list(porteinID = proteinID[idx],
                decoy_hit = decoy_hit,
                score = score,
                FDR1 = FDR1,
                FDR2 = FDR2))
}
