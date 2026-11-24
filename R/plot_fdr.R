#' Plot score distributions and q-values
#'
#' Two panels on a shared score axis: the histogram of target and decoy scores, and the q-values of both FDR
#' estimates in percent. When smaller scores are better, the score axis is logarithmic. See
#' `vignette("TargetDecoyFDR_Example", package = "prozor")`.
#'
#' @param fdr A result of [compute_fdr()] or [compute_fdr_with_id()].
#' @param bins Number of histogram bins.
#' @return A ggplot object.
#' @family FDR
#' @export
#' @examples
#' data(fdr_sample)
#' fdr <- compute_fdr_with_id(fdr_sample$score, fdr_sample$protein_id, larger_better = FALSE)
#' plot_fdr(fdr)
#' fdr2 <- compute_fdr_with_id(fdr_sample$score2, fdr_sample$protein_id, larger_better = TRUE)
#' plot_fdr(fdr2)
plot_fdr <- function(fdr, bins = 100) {
  panels <- c("hits", "q-value (%)")
  scores <- data.frame(
    score = fdr$score,
    # decoys drawn last, on top of the far more numerous targets
    hit = factor(ifelse(fdr$decoy_hit, "decoy", "target"), levels = c("target", "decoy")),
    panel = factor(panels[1], levels = panels)
  )
  q_values <- data.frame(
    score = rep(fdr$score, 2),
    q_value = 100 * c(fdr$q_value_fpr, fdr$q_value_simple_fdr),
    estimate = rep(c("FPR", "SimpleFDR"), each = length(fdr$score)),
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
    ggplot2::geom_line(data = q_values, ggplot2::aes(x = .data$score, y = .data$q_value, colour = .data$estimate)) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(target = "#0072B2", decoy = "#D55E00"), name = NULL) +
    ggplot2::scale_colour_manual(values = c(FPR = "#000000", SimpleFDR = "#E69F00"), name = "FDR estimate") +
    ggplot2::labs(
      x = if (fdr$larger_better) "score (larger is better)" else "score (smaller is better)",
      y = NULL
    )
  if (!fdr$larger_better) {
    p <- p + ggplot2::scale_x_log10()
  }
  p
}
