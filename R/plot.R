# ============================================================================
# plot.R — Visualization functions for AAPA
# ============================================================================

#' Plot score distribution
#'
#' Visualize the distribution of composite scores across all test
#' individuals and candidate families, highlighting assigned vs.
#' rejected individuals.
#'
#' @param result An \code{aapa_result} object.
#' @param type One of \code{"histogram"} or \code{"density"}.
#'   Default: \code{"histogram"}.
#' @return A ggplot2 object (if ggplot2 is available), otherwise a base
#'   R plot is produced and NULL is returned invisibly.
#' @family visualization
#' @export
plot_score_distribution <- function(result, type = "histogram") {
  checkmate::assert_class(result, "aapa_result")
  checkmate::assert_choice(type, c("histogram", "density"))
  asgn <- result$assignment

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    score_col <- "score"
    status_col <- "status"
    p <- ggplot2::ggplot(asgn, ggplot2::aes(x = .data[[score_col]],
                                             fill = .data[[status_col]])) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.7,
                              position = "identity") +
      ggplot2::labs(title = "AAPA Score Distribution",
                    x = "Composite Score (best family)",
                    y = "Count", fill = "Status") +
      ggplot2::theme_minimal()
    if (type == "density") {
      p <- ggplot2::ggplot(asgn, ggplot2::aes(x = .data[[score_col]],
                                               fill = .data[[status_col]])) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::labs(title = "AAPA Score Distribution",
                      x = "Composite Score (best family)",
                      y = "Density", fill = "Status") +
        ggplot2::theme_minimal()
    }
    return(p)
  }

  # Base R fallback
  assigned_scores <- asgn$score[asgn$status == "ASSIGNED"]
  rejected_scores <- asgn$score[asgn$status == "REJECT"]

  xlim <- range(c(assigned_scores, rejected_scores), na.rm = TRUE)
  hist(assigned_scores, col = grDevices::rgb(0.2, 0.6, 0.2, 0.5),
       xlim = xlim, main = "AAPA Score Distribution",
       xlab = "Composite Score", ylab = "Count")
  if (length(rejected_scores) > 0) {
    hist(rejected_scores, col = grDevices::rgb(0.8, 0.2, 0.2, 0.5),
         add = TRUE)
  }
  legend("topright", legend = c("Assigned", "Rejected"),
         fill = c(grDevices::rgb(0.2, 0.6, 0.2, 0.5),
                  grDevices::rgb(0.8, 0.2, 0.2, 0.5)))
  invisible(NULL)
}

#' Plot top-k candidates for an individual
#'
#' Show the top-k candidate families and their scores for a specific
#' test individual, useful for diagnosing ambiguous or rejected cases.
#'
#' @param result An \code{aapa_result} object.
#' @param individual_id Character; the ID of the individual to plot.
#' @return A ggplot2 object (if available), otherwise base R plot.
#' @family visualization
#' @export
plot_topk <- function(result, individual_id) {
  checkmate::assert_class(result, "aapa_result")
  checkmate::assert_string(individual_id)
  if (!individual_id %in% names(result$topk)) {
    cli::cli_abort(
      "Individual {.val {individual_id}} not found in results."
    )
  }

  topk <- result$topk[[individual_id]]
  topk$family_id <- factor(topk$family_id,
                           levels = rev(topk$family_id))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(topk, ggplot2::aes(x = .data[["score"]],
                                             y = .data[["family_id"]])) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::geom_text(ggplot2::aes(
        label = sprintf("C=%.3f", .data[["conflict"]])),
        hjust = -0.1, size = 3) +
      ggplot2::labs(title = paste("Top-k Candidates:", individual_id),
                    x = "Composite Score", y = "Family") +
      ggplot2::theme_minimal()
    return(p)
  }

  # Base R fallback
  barplot(rev(topk$score), names.arg = rev(topk$family_id),
          horiz = TRUE, col = "steelblue",
          main = paste("Top-k Candidates:", individual_id),
          xlab = "Composite Score")
  invisible(NULL)
}

#' Plot rejection diagnostics
#'
#' Visualize the relationship between score, confidence (gap), and
#' conflict rate, with rejection threshold lines.
#'
#' @param result An \code{aapa_result} object.
#' @return A ggplot2 object (if available), otherwise base R plot.
#' @family visualization
#' @export
plot_rejection_diagnostics <- function(result) {
  checkmate::assert_class(result, "aapa_result")
  asgn <- result$assignment
  params <- result$params

  # Get conflict for best family per individual
  best_conflict <- vapply(seq_len(nrow(asgn)), function(i) {
    tid <- asgn$individual_id[i]
    result$topk[[tid]]$conflict[1]
  }, numeric(1))

  plot_data <- data.frame(
    individual_id = asgn$individual_id,
    score = asgn$score,
    confidence = pmin(asgn$confidence, 2), # cap for visualization
    conflict = best_conflict,
    status = asgn$status,
    stringsAsFactors = FALSE
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = .data[["score"]],
                                      y = .data[["confidence"]],
                                      color = .data[["status"]],
                                      size = .data[["conflict"]])) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_vline(xintercept = params$tau_conf,
                          linetype = "dashed", color = "red") +
      ggplot2::geom_hline(yintercept = params$tau_rej,
                          linetype = "dashed", color = "orange") +
      ggplot2::labs(title = "AAPA Rejection Diagnostics",
                    x = "Best Score",
                    y = "Confidence (top1 - top2 gap)",
                    color = "Status",
                    size = "Conflict Rate") +
      ggplot2::theme_minimal()
    return(p)
  }

  # Base R fallback
  cols <- ifelse(plot_data$status == "ASSIGNED", "blue", "red")
  plot(plot_data$score, plot_data$confidence,
       col = cols, pch = 19,
       cex = 1 + plot_data$conflict * 3,
       main = "AAPA Rejection Diagnostics",
       xlab = "Best Score", ylab = "Confidence (gap)")
  abline(v = params$tau_conf, lty = 2, col = "red")
  abline(h = params$tau_rej, lty = 2, col = "orange")
  legend("topright", legend = c("Assigned", "Rejected"),
         col = c("blue", "red"), pch = 19)
  invisible(NULL)
}
