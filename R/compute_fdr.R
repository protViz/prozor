#' Compute the FDR of target-decoy search results from protein IDs
#'
#' Same as [compute_fdr()], but takes protein IDs and finds the decoys with the `decoy` pattern. See
#' `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param score Search engine scores.
#' @param id Protein IDs, one per score.
#' @param decoy Regular expression that matches decoy IDs; the default matches IDs that start with `REV_`.
#' @param larger_better `TRUE` if a larger score is better, `FALSE` if a smaller score is better.
#' @return The list returned by [compute_fdr()], plus the element `ID`: the protein IDs in the sorted order.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr_with_id(fdr_sample$score2, fdr_sample$proteinID, larger_better = TRUE)
#' head(as.data.frame(fdr))
compute_fdr_with_id <- function(score, id, decoy = "^REV_", larger_better = TRUE) {
  res <- compute_fdr(score, grepl(decoy, id), larger_better = larger_better)
  res$ID <- id[res$order]
  res
}


#' Compute the FDR of target-decoy search results
#'
#' Sorts the hits from best to worst score and estimates the false discovery rate at each score. Within tied
#' scores, decoys are ranked before targets, so the result does not depend on the input order and tied targets
#' get the conservative FDR. See `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param score Search engine scores.
#' @param decoy_hit Logical, `TRUE` for decoy hits, one per score.
#' @param larger_better `TRUE` if a larger score is better, `FALSE` if a smaller score is better.
#' @return A list of vectors in the sorted order:
#' * `larger_better`: the argument.
#' * `order`: the positions of the sorted hits in the input.
#' * `decoy_hit`, `score`: the sorted input.
#' * `FPR`: the FDR estimated as 2 * decoys / all hits (Elias and Gygi).
#' * `SimpleFDR`: the FDR estimated as decoys / targets (Kaell).
#' * `qValue_FPR`, `qValue_SimpleFDR`: the q-values of both estimates.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr(fdr_sample$score, grepl("^REV_", fdr_sample$proteinID), larger_better = FALSE)
#' head(as.data.frame(fdr))
compute_fdr <- function(score, decoy_hit, larger_better = TRUE) {
  ord <- order(score, decoy_hit, decreasing = c(larger_better, TRUE), method = "radix")
  score <- score[ord]
  decoy_hit <- decoy_hit[ord]

  fp <- cumsum(decoy_hit)
  tp <- seq_along(ord) - fp
  fpr <- (2 * fp) / (tp + fp)
  simple_fdr <- fp / tp
  list(
    larger_better = larger_better,
    order = ord,
    decoy_hit = decoy_hit,
    score = score,
    FPR = fpr,
    SimpleFDR = simple_fdr,
    qValue_FPR = rev(cummin(rev(fpr))),
    qValue_SimpleFDR = rev(cummin(rev(simple_fdr)))
  )
}


#' Find the score cut-off for a given FDR
#'
#' See `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param fdr A result of [compute_fdr()] or [compute_fdr_with_id()].
#' @param q_value The FDR in percent.
#' @param method The FDR estimate to use: `"SimpleFDR"` or `"FPR"`.
#' @return The least strict score whose q-value is at most `q_value` percent. An error if no score reaches it.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr_with_id(fdr_sample$score, fdr_sample$proteinID, larger_better = FALSE)
#' predict_score_fdr(fdr, q_value = 5)
#' fdr2 <- compute_fdr_with_id(fdr_sample$score2, fdr_sample$proteinID, larger_better = TRUE)
#' predict_score_fdr(fdr2, q_value = 5)
predict_score_fdr <- function(fdr, q_value = 1, method = "SimpleFDR") {
  q <- switch(
    method,
    FPR = fdr$qValue_FPR,
    SimpleFDR = fdr$qValue_SimpleFDR,
    stop("`method` must be \"FPR\" or \"SimpleFDR\", not \"", method, "\".")
  )
  valid_scores <- fdr$score[q <= q_value / 100]
  if (length(valid_scores) == 0) {
    stop("No score reaches a q-value of ", q_value, " % (method \"", method, "\").")
  }
  if (fdr$larger_better) min(valid_scores) else max(valid_scores)
}
