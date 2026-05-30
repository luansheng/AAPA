# ============================================================================
# qc.R — Quality control for AAPA
# ============================================================================

#' Quality control filter for genotype data
#'
#' Filter SNP markers and individuals based on missing data rates,
#' minor allele frequency, and other quality metrics.
#'
#' @param genotype Numeric matrix (individuals x SNPs), dosage 0/1/2/NA.
#' @param max_snp_missing Maximum allowable missing rate per SNP
#'   (default: 0.1). SNPs exceeding this are removed.
#' @param max_ind_missing Maximum allowable missing rate per individual
#'   (default: 0.2). Individuals exceeding this are removed.
#' @param min_maf Minimum minor allele frequency (default: 0.01).
#'   SNPs below this are removed.
#' @param verbose Logical; print QC summary? Default TRUE.
#' @return A list with components:
#'   \describe{
#'     \item{genotype}{Filtered genotype matrix.}
#'     \item{removed_snps}{Names of removed SNPs.}
#'     \item{removed_inds}{Names of removed individuals.}
#'     \item{summary}{QC summary statistics.}
#'   }
#' @family quality-control
#' @examples
#' sim <- simulate_aapa_data(n_families = 3, n_snps = 100)
#' qc_result <- qc_filter(sim$genotype, verbose = FALSE)
#' dim(qc_result$genotype)
#' @export
qc_filter <- function(genotype,
                      max_snp_missing = 0.1,
                      max_ind_missing = 0.2,
                      min_maf = 0.01,
                      verbose = TRUE) {
  stopifnot(is.matrix(genotype), is.numeric(genotype))
  n_ind_orig <- nrow(genotype)
  n_snp_orig <- ncol(genotype)

  # --- SNP missing rate filter ---
  snp_missing_rate <- colMeans(is.na(genotype))
  keep_snp_miss <- snp_missing_rate <= max_snp_missing
  removed_snp_miss <- colnames(genotype)[!keep_snp_miss]

  # --- Minor allele frequency filter ---
  # Compute MAF on non-missing data
  maf <- colMeans(genotype, na.rm = TRUE) / 2
  maf <- pmin(maf, 1 - maf)  # fold to minor allele
  keep_maf <- !is.na(maf) & maf >= min_maf
  removed_snp_maf <- colnames(genotype)[!keep_maf & keep_snp_miss]

  # Combined SNP filter
  keep_snp <- keep_snp_miss & keep_maf
  genotype <- genotype[, keep_snp, drop = FALSE]

  # --- Individual missing rate filter ---
  ind_missing_rate <- rowMeans(is.na(genotype))
  keep_ind <- ind_missing_rate <= max_ind_missing
  removed_inds <- rownames(genotype)[!keep_ind]
  genotype <- genotype[keep_ind, , drop = FALSE]

  removed_snps <- c(removed_snp_miss, removed_snp_maf)
  summary_stats <- list(
    n_ind_before = n_ind_orig,
    n_ind_after = nrow(genotype),
    n_ind_removed = length(removed_inds),
    n_snp_before = n_snp_orig,
    n_snp_after = ncol(genotype),
    n_snp_removed_miss = length(removed_snp_miss),
    n_snp_removed_maf = length(removed_snp_maf)
  )

  if (verbose) {
    cat("AAPA Quality Control Summary\n")
    cat("============================\n")
    cat(sprintf("Individuals: %d -> %d (removed %d)\n",
                n_ind_orig, nrow(genotype), length(removed_inds)))
    cat(sprintf("SNPs: %d -> %d (removed %d missing, %d low MAF)\n",
                n_snp_orig, ncol(genotype),
                length(removed_snp_miss), length(removed_snp_maf)))
  }

  list(
    genotype = genotype,
    removed_snps = removed_snps,
    removed_inds = removed_inds,
    summary = summary_stats
  )
}
