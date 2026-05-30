# ============================================================================
# test-qc.R — Unit tests for AAPA quality control
# ============================================================================

test_that("qc_filter removes high-missing SNPs", {
  geno <- matrix(c(
    0, 1, NA, 2,
    1, NA, NA, 0,
    2, 1, NA, 1,
    0, 0, NA, 2,
    1, 1, NA, 0
  ), nrow = 5, byrow = TRUE)
  rownames(geno) <- paste0("IND", 1:5)
  colnames(geno) <- paste0("SNP", 1:4)

  # SNP3 has 100% missing
  result <- qc_filter(geno, max_snp_missing = 0.5,
                      min_maf = 0, verbose = FALSE)
  expect_false("SNP3" %in% colnames(result$genotype))
  expect_true("SNP3" %in% result$removed_snps)
})

test_that("qc_filter removes low-MAF SNPs", {
  # Create a SNP that's nearly monomorphic
  geno <- matrix(c(
    0, 0, 1,
    0, 0, 2,
    0, 0, 0,
    0, 0, 1,
    0, 0, 2
  ), nrow = 5, byrow = TRUE)
  rownames(geno) <- paste0("IND", 1:5)
  colnames(geno) <- paste0("SNP", 1:3)

  # SNP1 has MAF=0 (monomorphic)
  result <- qc_filter(geno, max_snp_missing = 1.0,
                      min_maf = 0.05, verbose = FALSE)
  expect_false("SNP1" %in% colnames(result$genotype))
})

test_that("qc_filter removes high-missing individuals", {
  geno <- matrix(c(
    0,  1,  2,  1,
    NA, NA, NA, NA,  # 100% missing
    1,  0,  2,  1,
    2,  NA, 1,  0
  ), nrow = 4, byrow = TRUE)
  rownames(geno) <- paste0("IND", 1:4)
  colnames(geno) <- paste0("SNP", 1:4)

  result <- qc_filter(geno, max_snp_missing = 1.0,
                      max_ind_missing = 0.5,
                      min_maf = 0, verbose = FALSE)
  expect_false("IND2" %in% rownames(result$genotype))
  expect_true("IND2" %in% result$removed_inds)
})
