# ============================================================================
# scoring.R — Core scoring engine for AAPA
# ============================================================================

#' Compute Mendelian conflict rate
#'
#' For each test individual and each candidate family, calculate the
#' proportion of loci where the individual's genotype is incompatible
#' with the expected offspring genotypes given the parental cross.
#'
#' @param genotype Numeric matrix (individuals x SNPs), dosage 0/1/2/NA.
#' @param parents An \code{aapa_parents} object (from [read_parents()] or
#'   a named list of families with \code{sire_geno} and \code{dam_geno}).
#' @param test_ids Character vector of individual IDs to test. If NULL,
#'   all non-parent individuals in genotype are used.
#' @return A numeric matrix of dimension (n_test x n_families) with
#'   Mendelian conflict rates in [0, 1].
#' @family scoring
#' @examples
#' sim <- simulate_aapa_data(n_families = 3, n_snps = 100)
#' parents <- read_parents(
#'   file = textConnection(paste(
#'     apply(sim$parents, 1, paste, collapse = ","), collapse = "\n"
#'   )),
#'   genotype_matrix = sim$genotype
#' )
#' cm <- mendelian_conflict(sim$genotype, parents)
#' @export
mendelian_conflict <- function(genotype, parents, test_ids = NULL) {
  if (is.null(test_ids)) {
    parent_ids <- unique(unlist(lapply(parents, function(x) {
      c(x$sire_id, x$dam_id)
    })))
    test_ids <- setdiff(rownames(genotype), parent_ids)
  }

  test_geno <- genotype[test_ids, , drop = FALSE]
  n_test <- nrow(test_geno)
  n_fam  <- length(parents)
  n_snp  <- ncol(test_geno)

  # Pre-compute allowed offspring genotypes for each family and locus
  # For dosage coding: sire_alleles x dam_alleles -> possible offspring dosages
  conflict_mat <- matrix(NA_real_, nrow = n_test, ncol = n_fam)
  rownames(conflict_mat) <- test_ids
  colnames(conflict_mat) <- names(parents)

  for (fi in seq_along(parents)) {
    fam <- parents[[fi]]
    sg <- fam$sire_geno  # length M
    dg <- fam$dam_geno   # length M

    # Compute allowed genotype set per locus
    # Sire dosage s -> alleles (0: {0,0}, 1: {0,1}, 2: {1,1})
    # Dam dosage d  -> alleles (0: {0,0}, 1: {0,1}, 2: {1,1})
    # Possible offspring dosage = sire_allele + dam_allele
    allowed <- .compute_allowed_genotypes(sg, dg)

    for (ii in seq_len(n_test)) {
      ig <- test_geno[ii, ]
      # Valid loci: both individual and parents are non-missing
      valid <- !is.na(ig) & !is.na(sg) & !is.na(dg)
      n_valid <- sum(valid)
      if (n_valid == 0) {
        conflict_mat[ii, fi] <- NA_real_
        next
      }
      # Count conflicts
      n_conflict <- sum(!.is_compatible(ig[valid], allowed[valid, ,
                                                           drop = FALSE]))
      conflict_mat[ii, fi] <- n_conflict / n_valid
    }
  }
  conflict_mat
}

#' Compute allowed offspring genotypes per locus
#'
#' @param sire_geno Integer vector of sire dosages (0/1/2/NA).
#' @param dam_geno Integer vector of dam dosages (0/1/2/NA).
#' @return A logical matrix (n_loci x 3) where columns correspond to
#'   dosage 0, 1, 2 indicating whether each is an allowed offspring genotype.
#' @keywords internal
.compute_allowed_genotypes <- function(sire_geno, dam_geno) {
  n <- length(sire_geno)
  allowed <- matrix(FALSE, nrow = n, ncol = 3)
  colnames(allowed) <- c("0", "1", "2")

  for (j in seq_len(n)) {
    s <- sire_geno[j]
    d <- dam_geno[j]
    if (is.na(s) || is.na(d)) {
      allowed[j, ] <- TRUE  # unknown parent -> all allowed
      next
    }
    # Sire alleles
    s_alleles <- switch(as.character(s),
                        "0" = c(0, 0),
                        "1" = c(0, 1),
                        "2" = c(1, 1),
                        c(0, 1))  # fallback
    # Dam alleles
    d_alleles <- switch(as.character(d),
                        "0" = c(0, 0),
                        "1" = c(0, 1),
                        "2" = c(1, 1),
                        c(0, 1))
    # All possible offspring dosages
    possible <- unique(outer(s_alleles, d_alleles, "+"))
    allowed[j, as.character(possible)] <- TRUE
  }
  allowed
}

#' Check compatibility of individual genotypes with allowed set
#'
#' @param ind_geno Integer vector of individual dosages.
#' @param allowed Logical matrix (n_loci x 3).
#' @return Logical vector; TRUE if compatible.
#' @keywords internal
.is_compatible <- function(ind_geno, allowed) {
  vapply(seq_along(ind_geno), function(j) {
    g <- ind_geno[j]
    if (is.na(g) || g < 0 || g > 2) return(TRUE)
    allowed[j, as.character(g)]
  }, logical(1))
}

#' Compute anchor kinship scores
#'
#' For each test individual and each candidate family, compute the mean
#' IBS (Identity-By-State) similarity with the family's anchor individuals.
#'
#' @param genotype Numeric matrix (individuals x SNPs).
#' @param anchors An \code{aapa_anchors} object (from [read_anchors()]).
#' @param test_ids Character vector of test individual IDs.
#' @param method Kinship estimation method. Currently only \code{"ibs"}
#'   (proportion of IBS matches) is supported.
#' @return A numeric matrix (n_test x n_families) with kinship scores.
#' @family scoring
#' @export
anchor_kinship <- function(genotype, anchors, test_ids,
                           method = "ibs") {
  anchor_geno <- attr(anchors, "geno")
  families <- unique(anchors$family_id)

  kin_mat <- matrix(NA_real_, nrow = length(test_ids), ncol = length(families))
  rownames(kin_mat) <- test_ids
  colnames(kin_mat) <- families

  for (fi in seq_along(families)) {
    fid <- families[fi]
    anc_ids <- anchors$individual_id[anchors$family_id == fid]
    anc_weights <- anchors$weight[anchors$family_id == fid]

    if (length(anc_ids) == 0) next

    for (ii in seq_along(test_ids)) {
      tid <- test_ids[ii]
      ig <- genotype[tid, ]

      # Compute weighted mean IBS with anchors
      ibs_scores <- numeric(length(anc_ids))
      for (ai in seq_along(anc_ids)) {
        ag <- anchor_geno[anc_ids[ai], ]
        valid <- !is.na(ig) & !is.na(ag)
        if (sum(valid) == 0) {
          ibs_scores[ai] <- NA_real_
          next
        }
        # IBS proportion: fraction of matching allele dosages
        # More nuanced: IBS0/1/2 scoring
        diff <- abs(ig[valid] - ag[valid])
        # IBS score: 1 - mean(diff) / 2  (normalized to [0,1])
        ibs_scores[ai] <- 1 - mean(diff) / 2
      }
      # Weighted mean
      valid_scores <- !is.na(ibs_scores)
      if (any(valid_scores)) {
        w <- anc_weights[valid_scores]
        kin_mat[ii, fi] <- stats::weighted.mean(ibs_scores[valid_scores], w)
      }
    }
  }
  kin_mat
}

#' Compute composite assignment score
#'
#' Combine Mendelian conflict rate and anchor kinship score into a single
#' composite score for family assignment.
#'
#' \deqn{S_{i,f} = -\alpha \cdot C_{i,f} + \beta \cdot K_{i,f}}
#'
#' @param conflict_mat Numeric matrix of Mendelian conflict rates
#'   (from [mendelian_conflict()]).
#' @param kinship_mat Numeric matrix of anchor kinship scores
#'   (from [anchor_kinship()]).
#' @param alpha Weight for conflict penalty (default: 1.0).
#' @param beta Weight for kinship reward (default: 1.0).
#' @return A numeric matrix (n_test x n_families) of composite scores.
#'   Higher scores indicate better match.
#' @family scoring
#' @export
composite_score <- function(conflict_mat, kinship_mat,
                            alpha = 1.0, beta = 1.0) {
  # Ensure dimensions match
  if (!identical(dim(conflict_mat), dim(kinship_mat))) {
    # Align by shared row/column names
    shared_rows <- intersect(rownames(conflict_mat), rownames(kinship_mat))
    shared_cols <- intersect(colnames(conflict_mat), colnames(kinship_mat))
    if (length(shared_rows) == 0 || length(shared_cols) == 0) {
      stop("No shared individuals or families between conflict and kinship matrices.")
    }
    conflict_mat <- conflict_mat[shared_rows, shared_cols, drop = FALSE]
    kinship_mat  <- kinship_mat[shared_rows, shared_cols, drop = FALSE]
  }

  score_mat <- -alpha * conflict_mat + beta * kinship_mat
  score_mat
}
