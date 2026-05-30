# ============================================================================
# assign.R — Main assignment pipeline for AAPA
# ============================================================================

#' AAPA family assignment
#'
#' Perform anchor-assisted pedigree assignment for test individuals.
#' This is the main user-facing function that orchestrates the full
#' pipeline: Mendelian conflict scoring, anchor kinship scoring,
#' composite scoring, top-k pruning, and rejection filtering.
#'
#' @param genotype Numeric matrix (individuals x SNPs) with dosage
#'   encoding 0/1/2/NA.
#' @param parents An \code{aapa_parents} object or a named list of
#'   families (see \code{read_parents()}).
#' @param anchors An \code{aapa_anchors} object (see \code{read_anchors()}).
#' @param test_ids Character vector of individual IDs to assign. If NULL,
#'   all individuals not in parents or anchors are used.
#' @param alpha Weight for Mendelian conflict penalty (default: 1.0).
#' @param beta Weight for anchor kinship reward (default: 1.0).
#' @param top_k Number of top candidate families to retain (default: 5).
#' @param tau_conf Minimum score threshold for confident assignment
#'   (default: -0.05). Individuals whose best score falls below this
#'   are rejected.
#' @param tau_rej Minimum score gap between top-1 and top-2 families
#'   (default: 0.1). Individuals with smaller gaps are rejected.
#' @param max_conflict Maximum allowable Mendelian conflict rate
#'   (default: 0.1). Individuals whose best family exceeds this are
#'   rejected.
#' @return An object of class \code{aapa_result}, a list containing:
#'   \describe{
#'     \item{assignment}{Data frame with columns: individual_id,
#'       assigned_family, score, confidence, status (ASSIGNED/REJECT).}
#'     \item{topk}{List of per-individual top-k candidate data frames.}
#'     \item{conflict_matrix}{Full Mendelian conflict rate matrix.}
#'     \item{kinship_matrix}{Full anchor kinship score matrix.}
#'     \item{score_matrix}{Full composite score matrix.}
#'     \item{params}{List of parameters used.}
#'   }
#' @family assignment
#' @examples
#' # Simulate data and run the full pipeline
#' sim <- simulate_aapa_data(n_families = 3, n_snps = 100, seed = 1)
#'
#' # Write temp files for reading functions
#' tmp_parents <- tempfile(fileext = ".csv")
#' write.csv(sim$parents, tmp_parents, row.names = FALSE)
#' parents <- read_parents(tmp_parents, sim$genotype)
#'
#' tmp_anchors <- tempfile(fileext = ".csv")
#' write.csv(sim$anchors, tmp_anchors, row.names = FALSE)
#' anchors <- read_anchors(tmp_anchors, sim$genotype)
#'
#' result <- aapa_assign(sim$genotype, parents, anchors)
#' print(result)
#'
#' unlink(c(tmp_parents, tmp_anchors))
#' @export
aapa_assign <- function(genotype, parents, anchors,
                        test_ids = NULL,
                        alpha = 1.0, beta = 1.0,
                        top_k = 5,
                        tau_conf = -0.05,
                        tau_rej = 0.1,
                        max_conflict = 0.1) {
  # --- Parameter validation ---
  checkmate::assert_matrix(genotype,
    mode = "numeric", min.rows = 1,
    min.cols = 1
  )
  checkmate::assert_list(parents, min.len = 1)
  checkmate::assert_character(test_ids, null.ok = TRUE)
  checkmate::assert_number(alpha, lower = 0)
  checkmate::assert_number(beta, lower = 0)
  checkmate::assert_count(top_k, positive = TRUE)
  checkmate::assert_number(tau_conf)
  checkmate::assert_number(tau_rej)
  checkmate::assert_number(max_conflict, lower = 0, upper = 1)

  .validate_genotype_matrix(genotype)
  .validate_parents_object(parents)
  parents <- .align_parents_to_markers(parents, colnames(genotype))
  .validate_anchors_object(anchors, family_ids = names(parents))
  anchors <- .align_anchors_to_markers(anchors, colnames(genotype))

  # Determine test individuals
  if (is.null(test_ids)) {
    parent_ids <- unique(unlist(lapply(parents, function(x) {
      c(x$sire_id, x$dam_id)
    })))
    anchor_ids <- if (inherits(anchors, "aapa_anchors")) {
      anchors$individual_id
    } else {
      character(0)
    }
    test_ids <- setdiff(
      rownames(genotype),
      union(parent_ids, anchor_ids)
    )
  }
  if (length(test_ids) == 0) {
    cli::cli_abort("No test individuals found.")
  }

  # Verify test_ids exist in genotype matrix
  missing <- setdiff(test_ids, rownames(genotype))
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Test individuals not found in genotype matrix.",
      "x" = "Missing: {.val {head(missing, 5)}}"
    ))
  }

  # --- Step 1: Mendelian conflict ---
  cli::cli_alert_info("Step 1/4: Computing Mendelian conflict rates...")
  conflict_mat <- mendelian_conflict(genotype, parents, test_ids)

  # --- Step 2: Anchor kinship ---
  cli::cli_alert_info("Step 2/4: Computing anchor kinship scores...")
  kinship_mat <- anchor_kinship(genotype, anchors, test_ids)

  # --- Step 3: Composite score ---
  cli::cli_alert_info("Step 3/4: Computing composite scores...")
  score_mat <- composite_score(conflict_mat, kinship_mat, alpha, beta)

  # --- Step 4: Top-k pruning and assignment ---
  cli::cli_alert_info("Step 4/4: Top-k pruning and rejection filtering...")
  n_fam <- ncol(score_mat)
  effective_k <- min(top_k, n_fam)

  assignment <- data.frame(
    individual_id = character(0),
    assigned_family = character(0),
    score = numeric(0),
    confidence = numeric(0),
    status = character(0),
    reject_reason = character(0),
    stringsAsFactors = FALSE
  )

  topk_list <- list()

  for (ii in seq_along(test_ids)) {
    tid <- test_ids[ii]
    scores <- score_mat[tid, ]
    conflicts <- conflict_mat[tid, ]

    # Sort by score (descending)
    ord <- order(scores, decreasing = TRUE)
    topk_idx <- ord[seq_len(effective_k)]

    topk_df <- data.frame(
      rank = seq_len(effective_k),
      family_id = colnames(score_mat)[topk_idx],
      score = scores[topk_idx],
      conflict = conflicts[topk_idx],
      stringsAsFactors = FALSE
    )
    topk_list[[tid]] <- topk_df

    # Best candidate
    best_fam <- topk_df$family_id[1]
    best_score <- topk_df$score[1]
    best_conflict <- topk_df$conflict[1]

    # Confidence: gap between top-1 and top-2
    if (effective_k >= 2) {
      gap <- topk_df$score[1] - topk_df$score[2]
    } else {
      gap <- Inf
    }

    # --- Step 5: Rejection rules ---
    reject <- FALSE
    reason <- ""

    # Rule 1: Absolute score threshold
    if (!is.na(best_score) && best_score < tau_conf) {
      reject <- TRUE
      reason <- paste0(reason, "score_below_threshold;")
    }

    # Rule 2: Relative gap threshold
    if (!is.na(gap) && gap < tau_rej) {
      reject <- TRUE
      reason <- paste0(reason, "insufficient_gap;")
    }

    # Rule 3: Conflict upper limit
    if (!is.na(best_conflict) && best_conflict > max_conflict) {
      reject <- TRUE
      reason <- paste0(reason, "conflict_too_high;")
    }

    # Rule 4: All candidates low confidence (unknown family suspect)
    if (all(!is.na(conflicts) & conflicts > max_conflict)) {
      reject <- TRUE
      reason <- paste0(reason, "all_candidates_high_conflict;")
    }

    status <- if (reject) "REJECT" else "ASSIGNED"
    assigned <- if (reject) "REJECT" else best_fam

    assignment <- rbind(assignment, data.frame(
      individual_id = tid,
      assigned_family = assigned,
      score = best_score,
      confidence = gap,
      status = status,
      reject_reason = reason,
      stringsAsFactors = FALSE
    ))
  }

  cli::cli_alert_success(
    sprintf(
      "Assignment complete: %d assigned, %d rejected",
      sum(assignment$status == "ASSIGNED"),
      sum(assignment$status == "REJECT")
    )
  )

  result <- list(
    assignment = assignment,
    topk = topk_list,
    conflict_matrix = conflict_mat,
    kinship_matrix = kinship_mat,
    score_matrix = score_mat,
    params = list(
      alpha = alpha,
      beta = beta,
      top_k = top_k,
      tau_conf = tau_conf,
      tau_rej = tau_rej,
      max_conflict = max_conflict
    )
  )
  class(result) <- "aapa_result"
  result
}

#' Print method for aapa_result
#'
#' @param x An \code{aapa_result} object.
#' @param ... Additional arguments (ignored).
#' @family assignment
#' @export
print.aapa_result <- function(x, ...) {
  asgn <- x$assignment
  n_total <- nrow(asgn)
  n_assigned <- sum(asgn$status == "ASSIGNED")
  n_rejected <- sum(asgn$status == "REJECT")
  pct_assigned <- round(100 * n_assigned / n_total, 1)
  pct_rejected <- round(100 * n_rejected / n_total, 1)

  cli::cli_h1("AAPA Assignment Result")
  cli::cli_text(sprintf("Total individuals: %d", n_total))
  cli::cli_text(sprintf("Assigned: %d (%.1f%%)", n_assigned, pct_assigned))
  cli::cli_text(sprintf("Rejected: %d (%.1f%%)", n_rejected, pct_rejected))
  cli::cli_text(sprintf("Families: %d", ncol(x$score_matrix)))
  cli::cli_text(
    sprintf(
      "Parameters: alpha=%s, beta=%s, top_k=%s",
      x$params$alpha,
      x$params$beta,
      x$params$top_k
    )
  )
  invisible(x)
}

#' Summary method for aapa_result
#'
#' @param object An \code{aapa_result} object.
#' @param ... Additional arguments (ignored).
#' @family assignment
#' @export
summary.aapa_result <- function(object, ...) {
  asgn <- object$assignment
  cli::cli_h1("AAPA Assignment Summary")

  # Per-family counts
  assigned <- asgn[asgn$status == "ASSIGNED", ]
  if (nrow(assigned) > 0) {
    fam_counts <- table(assigned$assigned_family)
    cli::cli_h2("Assignments per family")
    print(fam_counts)
  }

  # Rejection reasons
  rejected <- asgn[asgn$status == "REJECT", ]
  if (nrow(rejected) > 0) {
    cli::cli_h2("Rejection summary")
    reasons <- unlist(strsplit(rejected$reject_reason, ";"))
    reasons <- reasons[reasons != ""]
    print(table(reasons))
  }

  # Score distribution
  cli::cli_h2("Score statistics (best candidate)")
  print(summary(asgn$score))
  cli::cli_h2("Confidence (top1-top2 gap) statistics")
  finite_conf <- asgn$confidence[is.finite(asgn$confidence)]
  if (length(finite_conf) > 0) print(summary(finite_conf))
  invisible(object)
}
