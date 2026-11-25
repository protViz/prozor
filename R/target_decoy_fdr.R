# Stops unless `score` is numeric, `decoy_hit` is logical, `score`, `decoy_hit` and `id` (if given) have the same
# length, `score` and `decoy_hit` have no NA, and `larger_better` is TRUE or FALSE
.check_target_decoy <- function(score, decoy_hit, larger_better, id) {
  if (!is.numeric(score)) {
    stop("`score` must be numeric, not ", class(score)[1], ".", call. = FALSE)
  }
  if (!is.logical(decoy_hit)) {
    stop("`decoy_hit` must be logical, not ", class(decoy_hit)[1], ".", call. = FALSE)
  }
  if (length(score) != length(decoy_hit)) {
    stop(
      "`score` and `decoy_hit` must have the same length, not ",
      length(score),
      " and ",
      length(decoy_hit),
      ".",
      call. = FALSE
    )
  }
  if (!is.null(id) && length(id) != length(score)) {
    stop("`score` and `id` must have the same length, not ", length(score), " and ", length(id), ".", call. = FALSE)
  }
  if (anyNA(score)) {
    stop("`score` must not contain NA; found ", sum(is.na(score)), ".", call. = FALSE)
  }
  if (anyNA(decoy_hit)) {
    stop("`decoy_hit` must not contain NA; found ", sum(is.na(decoy_hit)), ".", call. = FALSE)
  }
  if (!(isTRUE(larger_better) || isFALSE(larger_better))) {
    stop("`larger_better` must be TRUE or FALSE, not ", deparse1(larger_better), ".", call. = FALSE)
  }
}

# Stops unless `q_value` is a single number between 0 and 1
.check_q_value <- function(q_value) {
  if (!(is.numeric(q_value) && length(q_value) == 1 && isTRUE(q_value >= 0 & q_value <= 1))) {
    stop("`q_value` must be a fraction between 0 and 1 (0.01 for 1 %), not ", deparse1(q_value), ".", call. = FALSE)
  }
}


#' Target-decoy false discovery rate
#'
#' @description
#' An R6 class that estimates the false discovery rate (FDR) of target-decoy search results. When the object is
#' created, the hits are sorted from best to worst score and two FDR estimates are computed at each score:
#' * `fdr_concatenated`: the FDR of a concatenated target-decoy search, estimated as 2 * decoys / all hits (Elias and
#'   Gygi 2007).
#' * `fdr_simple`: the FDR estimated as decoys / targets (Käll et al. 2008).
#'
#' Both estimates are fractions in \[0, 1\]: values above 1, e.g. where the best hits are decoys, are set to 1.
#'
#' The q-value of a hit is the lowest FDR at which it is accepted: the smallest estimate at its own or any worse
#' score. Within tied scores, decoys are ranked before targets, so the result does not depend on the input order and
#' tied targets get the conservative FDR.
#'
#' The decoys are given as logical flags. For protein IDs whose decoys start with `REV_`, pass
#' `decoy_hit = grepl("^REV_", id)`, and `id = id` to keep the IDs in `table`.
#'
#' The object cannot be changed: its fields are read-only and its methods only read them. `clone()` copies an object.
#' See `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @return `TargetDecoyFDR$new()` returns a `TargetDecoyFDR` R6 object; see Methods.
#' @references
#' Elias JE, Gygi SP (2007). Target-decoy search strategy for increased confidence in large-scale protein
#' identifications by mass spectrometry. *Nature Methods*, 4(3), 207-214. \doi{10.1038/nmeth1019}
#'
#' Käll L, Storey JD, MacCoss MJ, Noble WS (2008). Assigning significance to peptides identified by tandem mass
#' spectrometry using decoy databases. *Journal of Proteome Research*, 7(1), 29-34. \doi{10.1021/pr700600n}
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' # smaller scores are better; decoy IDs start with REV_
#' fdr <- TargetDecoyFDR$new(
#'   fdr_sample$score,
#'   decoy_hit = grepl("^REV_", fdr_sample$protein_id),
#'   larger_better = FALSE,
#'   id = fdr_sample$protein_id
#' )
#' fdr
#' head(fdr$table)
#' fdr$score_at(q_value = 0.01)
#' fdr$score_at(q_value = 0.01, method = "fdr_concatenated")
#' fdr$plot()
#' # larger scores are better
#' fdr2 <- TargetDecoyFDR$new(
#'   fdr_sample$neg_log10_score,
#'   decoy_hit = grepl("^REV_", fdr_sample$protein_id)
#' )
#' fdr2$score_at(q_value = 0.01)
#' fdr2$plot()
TargetDecoyFDR <- R6::R6Class(
  "TargetDecoyFDR",
  public = list(
    #' @description
    #' Sort the hits from best to worst score and compute both FDR estimates and their q-values.
    #' @param score Search engine scores, without `NA`.
    #' @param decoy_hit Logical, `TRUE` for decoy hits, one per score, without `NA`.
    #' @param larger_better `TRUE` if a larger score is better, `FALSE` if a smaller score is better.
    #' @param id IDs of the hits, e.g. protein IDs, one per score; `NULL` (default) for none.
    #' @return A new `TargetDecoyFDR` object.
    initialize = function(score, decoy_hit, larger_better = TRUE, id = NULL) {
      .check_target_decoy(score, decoy_hit, larger_better, id)
      ord <- order(score, decoy_hit, decreasing = c(larger_better, TRUE), method = "radix")
      decoy_hit <- decoy_hit[ord]
      fp <- cumsum(decoy_hit)
      tp <- seq_along(ord) - fp
      # capped at 1: fp / tp is Inf while all hits are decoys
      fdr_concatenated <- pmin(1, (2 * fp) / (tp + fp))
      fdr_simple <- pmin(1, fp / tp)
      columns <- c(
        list(order = ord),
        if (!is.null(id)) list(id = id[ord]),
        list(
          score = score[ord],
          decoy_hit = decoy_hit,
          fdr_concatenated = fdr_concatenated,
          fdr_simple = fdr_simple,
          q_value_fdr_concatenated = rev(cummin(rev(fdr_concatenated))),
          q_value_fdr_simple = rev(cummin(rev(fdr_simple)))
        )
      )
      private$.larger_better <- larger_better
      private$.table <- data.frame(columns, row.names = NULL)
    },
    #' @description
    #' Print the numbers of hits and decoys, `larger_better`, and the scores at `fdr_simple` q-values of 0.01 and 0.05
    #' (`none` if no score reaches them).
    #' @param ... Not used.
    #' @return The object, invisibly.
    print = function(...) {
      cutoffs <- vapply(
        c(0.01, 0.05),
        function(q_value) {
          if (any(private$.table$q_value_fdr_simple <= q_value)) format(self$score_at(q_value), digits = 4) else "none"
        },
        character(1)
      )
      n_hits <- format(nrow(private$.table), big.mark = ",")
      n_decoys <- format(sum(private$.table$decoy_hit), big.mark = ",")
      cat(sprintf(
        "<TargetDecoyFDR> hits: %s, decoys: %s, larger_better: %s\n",
        n_hits,
        n_decoys,
        private$.larger_better
      ))
      cat(sprintf("score at fdr_simple q-value 0.01: %s, at 0.05: %s\n", cutoffs[1], cutoffs[2]))
      invisible(self)
    },
    #' @description
    #' Find the score cut-off for a given FDR.
    #' @param q_value The FDR as a fraction, like the q-values in `table`: `0.01` for 1 %.
    #' @param method The FDR estimate to use: `"fdr_simple"` (default) or `"fdr_concatenated"`.
    #' @return The least strict score whose q-value is at most `q_value`. An error if no score reaches it.
    score_at = function(q_value = 0.01, method = c("fdr_simple", "fdr_concatenated")) {
      .check_q_value(q_value)
      method <- .match_method(method, c("fdr_simple", "fdr_concatenated"))
      q <- private$.table[[paste0("q_value_", method)]]
      valid_scores <- private$.table$score[q <= q_value]
      if (length(valid_scores) == 0) {
        stop("No score reaches a q-value of ", q_value, " (method \"", method, "\").", call. = FALSE)
      }
      if (private$.larger_better) min(valid_scores) else max(valid_scores)
    },
    #' @description
    #' Plot the score distributions and q-values: two panels on a shared score axis, the histogram of target and decoy
    #' scores, and the q-values of both FDR estimates in percent. When smaller scores are better and all scores are
    #' positive, the score axis is logarithmic.
    #' @param bins Number of histogram bins.
    #' @return A ggplot object.
    plot = function(bins = 100) {
      fdr <- private$.table
      panels <- c("hits", "q-value (%)")
      scores <- data.frame(
        score = fdr$score,
        # decoys drawn last, on top of the far more numerous targets
        hit = factor(ifelse(fdr$decoy_hit, "decoy", "target"), levels = c("target", "decoy")),
        panel = factor(panels[1], levels = panels)
      )
      q_values <- data.frame(
        score = rep(fdr$score, 2),
        q_value = 100 * c(fdr$q_value_fdr_concatenated, fdr$q_value_fdr_simple),
        estimate = rep(c("fdr_concatenated", "fdr_simple"), each = nrow(fdr)),
        panel = factor(panels[2], levels = panels)
      )
      p <- ggplot2::ggplot() +
        ggplot2::geom_histogram(
          data = scores,
          ggplot2::aes(x = .data$score, fill = .data$hit),
          bins = bins,
          position = "identity",
          alpha = 0.7
        ) +
        ggplot2::geom_line(
          data = q_values,
          ggplot2::aes(x = .data$score, y = .data$q_value, colour = .data$estimate)
        ) +
        ggplot2::facet_grid(panel ~ ., scales = "free_y") +
        ggplot2::scale_fill_manual(values = c(target = "#0072B2", decoy = "#D55E00"), name = NULL) +
        ggplot2::scale_colour_manual(
          values = c(fdr_concatenated = "#000000", fdr_simple = "#E69F00"),
          name = "FDR estimate"
        ) +
        ggplot2::labs(
          x = if (private$.larger_better) "score (larger is better)" else "score (smaller is better)",
          y = NULL
        )
      if (!private$.larger_better && all(fdr$score > 0)) {
        p <- p + ggplot2::scale_x_log10()
      }
      p
    }
  ),
  active = list(
    #' @field larger_better `TRUE` if a larger score is better (read-only).
    larger_better = function(value) {
      if (!missing(value)) {
        .stop_read_only("larger_better", "TargetDecoyFDR")
      }
      private$.larger_better
    },
    #' @field n Number of hits (read-only).
    n = function(value) {
      if (!missing(value)) {
        .stop_read_only("n", "TargetDecoyFDR")
      }
      nrow(private$.table)
    },
    #' @field table The hits from best to worst score, as a data frame (read-only) with the columns:
    #' * `order`: the position of the hit in the input.
    #' * `id`: the ID of the hit; only when `id` was given.
    #' * `score`, `decoy_hit`: the sorted input.
    #' * `fdr_concatenated`, `fdr_simple`: the FDR estimates at the score, as fractions in \[0, 1\].
    #' * `q_value_fdr_concatenated`, `q_value_fdr_simple`: the q-values of both estimates, as fractions in \[0, 1\].
    table = function(value) {
      if (!missing(value)) {
        .stop_read_only("table", "TargetDecoyFDR")
      }
      private$.table
    }
  ),
  private = list(
    .larger_better = NULL,
    .table = NULL
  )
)
