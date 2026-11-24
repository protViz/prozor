#' Compute the FDR of target-decoy search results from protein IDs
#'
#' Same as [compute_fdr()], but takes protein IDs and finds the decoys with the `decoy` pattern. See
#' `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param score Search engine scores, without `NA`.
#' @param id Protein IDs, one per score.
#' @param decoy Regular expression that matches decoy IDs; the default matches IDs that start with `REV_`.
#' @param larger_better `TRUE` if a larger score is better, `FALSE` if a smaller score is better.
#' @return The list returned by [compute_fdr()], plus the element `id`: the protein IDs in the sorted order.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr_with_id(fdr_sample$score2, fdr_sample$protein_id, larger_better = TRUE)
#' head(as.data.frame(fdr))
compute_fdr_with_id <- function(score, id, decoy = "^REV_", larger_better = TRUE) {
  res <- compute_fdr(score, grepl(decoy, id), larger_better = larger_better)
  res$id <- id[res$order]
  res
}


#' Compute the FDR of target-decoy search results
#'
#' Sorts the hits from best to worst score and estimates the false discovery rate at each score. Within tied
#' scores, decoys are ranked before targets, so the result does not depend on the input order and tied targets
#' get the conservative FDR. See `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param score Search engine scores, without `NA`.
#' @param decoy_hit Logical, `TRUE` for decoy hits, one per score, without `NA`.
#' @param larger_better `TRUE` if a larger score is better, `FALSE` if a smaller score is better.
#' @return A list of vectors in the sorted order:
#' * `larger_better`: the argument.
#' * `order`: the positions of the sorted hits in the input.
#' * `decoy_hit`, `score`: the sorted input.
#' * `fdr_concatenated`: the FDR of a concatenated target-decoy search, estimated as 2 * decoys / all hits (Elias
#'   and Gygi).
#' * `simple_fdr`: the FDR estimated as decoys / targets (Kaell).
#' * `q_value_fdr_concatenated`, `q_value_simple_fdr`: the q-values of both estimates, as fractions.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr(fdr_sample$score, grepl("^REV_", fdr_sample$protein_id), larger_better = FALSE)
#' head(as.data.frame(fdr))
compute_fdr <- function(score, decoy_hit, larger_better = TRUE) {
  if (!is.logical(decoy_hit)) {
    stop("`decoy_hit` must be logical, not ", class(decoy_hit)[1], ".")
  }
  if (length(score) != length(decoy_hit)) {
    stop("`score` and `decoy_hit` must have the same length, not ", length(score), " and ", length(decoy_hit), ".")
  }
  if (anyNA(score)) {
    stop("`score` must not contain NA; found ", sum(is.na(score)), ".")
  }
  if (anyNA(decoy_hit)) {
    stop("`decoy_hit` must not contain NA; found ", sum(is.na(decoy_hit)), ".")
  }
  ord <- order(score, decoy_hit, decreasing = c(larger_better, TRUE), method = "radix")
  score <- score[ord]
  decoy_hit <- decoy_hit[ord]

  fp <- cumsum(decoy_hit)
  tp <- seq_along(ord) - fp
  fdr_concatenated <- (2 * fp) / (tp + fp)
  simple_fdr <- fp / tp
  list(
    larger_better = larger_better,
    order = ord,
    decoy_hit = decoy_hit,
    score = score,
    fdr_concatenated = fdr_concatenated,
    simple_fdr = simple_fdr,
    q_value_fdr_concatenated = rev(cummin(rev(fdr_concatenated))),
    q_value_simple_fdr = rev(cummin(rev(simple_fdr)))
  )
}


#' Find the score cut-off for a given FDR
#'
#' See `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param fdr A result of [compute_fdr()] or [compute_fdr_with_id()].
#' @param q_value The FDR as a fraction, like the q-values in `fdr`: `0.01` for 1 %.
#' @param method The FDR estimate to use: `"simple_fdr"` or `"fdr_concatenated"`.
#' @return The least strict score whose q-value is at most `q_value`. An error if no score reaches it.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr_with_id(fdr_sample$score, fdr_sample$protein_id, larger_better = FALSE)
#' predict_score_fdr(fdr, q_value = 0.01)
#' fdr2 <- compute_fdr_with_id(fdr_sample$score2, fdr_sample$protein_id, larger_better = TRUE)
#' predict_score_fdr(fdr2, q_value = 0.01, method = "fdr_concatenated")
predict_score_fdr <- function(fdr, q_value = 0.01, method = "simple_fdr") {
  if (!(q_value >= 0 && q_value <= 1)) {
    stop("`q_value` must be a fraction between 0 and 1 (0.01 for 1 %), not ", q_value, ".")
  }
  q <- switch(
    method,
    fdr_concatenated = fdr$q_value_fdr_concatenated,
    simple_fdr = fdr$q_value_simple_fdr,
    stop("`method` must be \"simple_fdr\" or \"fdr_concatenated\", not \"", method, "\".")
  )
  valid_scores <- fdr$score[q <= q_value]
  if (length(valid_scores) == 0) {
    stop("No score reaches a q-value of ", q_value, " (method \"", method, "\").")
  }
  if (fdr$larger_better) min(valid_scores) else max(valid_scores)
}
