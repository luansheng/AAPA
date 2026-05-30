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
#'   families (see [read_parents()]).
#' @param anchors An \code{aapa_anchors} object (see [read_anchors()]).
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
  stopifnot(is.matrix(genotype), is.numeric(genotype))
  stopifnot(is.list(parents), length(parents) > 0)
  stopifnot(alpha >= 0, beta >= 0)
  stopifnot(top_k >= 1)

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
    test_ids <- setdiff(rownames(genotype),
                        union(parent_ids, anchor_ids))
  }
  if (length(test_ids) == 0) {
    stop("No test individuals found.")
  }

  # Verify test_ids exist in genotype matrix
  missing <- setdiff(test_ids, rownames(genotype))
  if (length(missing) > 0) {
    stop("Test individuals not found in genotype matrix: ",
         paste(head(missing, 5), collapse = ", "))
  }

  # --- Step 1: Mendelian conflict ---
  message("Computing Mendelian conflict rates...")
  conflict_mat <- mendelian_conflict(genotype, parents, test_ids)

  # --- Step 2: Anchor kinship ---
  message("Computing anchor kinship scores...")
  kinship_mat <- anchor_kinship(genotype, anchors, test_ids)

  # --- Step 3: Composite score ---
  message("Computing composite scores...")
  score_mat <- composite_score(conflict_mat, kinship_mat, alpha, beta)

  # --- Step 4: Top-k pruning and assignment ---
  message("Performing top-k pruning and rejection filtering...")
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

  cat("AAPA Assignment Result\n")
  cat("======================\n")
  cat("Total individuals: ", n_total, "\n")
  cat("Assigned:          ", n_assigned,
      sprintf(" (%.1f%%)\n", 100 * n_assigned / n_total))
  cat("Rejected:          ", n_rejected,
      sprintf(" (%.1f%%)\n", 100 * n_rejected / n_total))
  cat("Families:          ", ncol(x$score_matrix), "\n")
  cat("Parameters:        ",
      sprintf("alpha=%.2f, beta=%.2f, top_k=%d\n",
              x$params$alpha, x$params$beta, x$params$top_k))
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
  cat("AAPA Assignment Summary\n")
  cat("=======================\n\n")

  # Per-family counts
  assigned <- asgn[asgn$status == "ASSIGNED", ]
  if (nrow(assigned) > 0) {
    fam_counts <- table(assigned$assigned_family)
    cat("Assignments per family:\n")
    print(fam_counts)
    cat("\n")
  }

  # Rejection reasons
  rejected <- asgn[asgn$status == "REJECT", ]
  if (nrow(rejected) > 0) {
    cat("Rejection summary:\n")
    reasons <- unlist(strsplit(rejected$reject_reason, ";"))
    reasons <- reasons[reasons != ""]
    print(table(reasons))
    cat("\n")
  }

  # Score distribution
  cat("Score statistics (best candidate):\n")
  print(summary(asgn$score))
  cat("\nConfidence (top1-top2 gap) statistics:\n")
  finite_conf <- asgn$confidence[is.finite(asgn$confidence)]
  if (length(finite_conf) > 0) print(summary(finite_conf))
  invisible(object)
}
