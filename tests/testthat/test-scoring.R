# ============================================================================
# test-scoring.R — Unit tests for AAPA scoring functions
# ============================================================================

test_that("mendelian_conflict returns correct dimensions", {
  sim <- simulate_aapa_data(
    n_families = 3, n_snps = 100,
    n_offspring_per_family = 5,
    n_anchors_per_family = 2,
    n_unknown = 2,
    missing_rate = 0, error_rate = 0
  )

  geno <- sim$genotype
  parents_df <- sim$parents
  anchors_df <- sim$anchors

  # Build parents object manually
  parents <- lapply(seq_len(nrow(parents_df)), function(i) {
    list(
      family_id = parents_df$family_id[i],
      sire_id   = parents_df$sire_id[i],
      dam_id    = parents_df$dam_id[i],
      sire_geno = geno[parents_df$sire_id[i], ],
      dam_geno  = geno[parents_df$dam_id[i], ]
    )
  })
  names(parents) <- parents_df$family_id
  class(parents) <- "aapa_parents"

  # Get test individuals (exclude parents and anchors)
  parent_ids <- c(parents_df$sire_id, parents_df$dam_id)
  anchor_ids <- anchors_df$individual_id
  test_ids <- setdiff(rownames(geno), c(parent_ids, anchor_ids))

  conflict_mat <- mendelian_conflict(geno, parents, test_ids)

  expect_equal(nrow(conflict_mat), length(test_ids))
  expect_equal(ncol(conflict_mat), 3)  # 3 families
  expect_true(all(conflict_mat >= 0 & conflict_mat <= 1, na.rm = TRUE))
})

test_that("true offspring have low conflict with own family", {
  sim <- simulate_aapa_data(
    n_families = 3, n_snps = 200,
    n_offspring_per_family = 5,
    n_anchors_per_family = 2,
    n_unknown = 0,
    missing_rate = 0, error_rate = 0
  )

  geno <- sim$genotype
  parents_df <- sim$parents

  parents <- lapply(seq_len(nrow(parents_df)), function(i) {
    list(
      family_id = parents_df$family_id[i],
      sire_id   = parents_df$sire_id[i],
      dam_id    = parents_df$dam_id[i],
      sire_geno = geno[parents_df$sire_id[i], ],
      dam_geno  = geno[parents_df$dam_id[i], ]
    )
  })
  names(parents) <- parents_df$family_id
  class(parents) <- "aapa_parents"

  # Test offspring of FAM001
  offspring_ids <- grep("^FAM001_OFF", rownames(geno), value = TRUE)
  conflict_mat <- mendelian_conflict(geno, parents, offspring_ids)

  # With no errors, offspring should have 0 conflict with own family
  expect_true(all(conflict_mat[, "FAM001"] == 0))
})

test_that("composite_score correctly combines conflict and kinship", {
  conflict <- matrix(c(0.0, 0.5, 0.1, 0.3), nrow = 2, ncol = 2)
  kinship  <- matrix(c(0.8, 0.3, 0.6, 0.4), nrow = 2, ncol = 2)
  rownames(conflict) <- rownames(kinship) <- c("IND1", "IND2")
  colnames(conflict) <- colnames(kinship) <- c("FAM1", "FAM2")

  score <- composite_score(conflict, kinship, alpha = 1, beta = 1)
  expected <- -1 * conflict + 1 * kinship
  expect_equal(score, expected)

  # With different weights
  score2 <- composite_score(conflict, kinship, alpha = 2, beta = 0.5)
  expected2 <- -2 * conflict + 0.5 * kinship
  expect_equal(score2, expected2)
})

test_that(".compute_allowed_genotypes handles all parental combinations", {
  # Both parents homozygous 0 -> offspring must be 0
  allowed <- aapa:::.compute_allowed_genotypes(0L, 0L)
  expect_true(allowed[1, "0"])
  expect_false(allowed[1, "1"])
  expect_false(allowed[1, "2"])

  # Both parents homozygous 2 -> offspring must be 2
  allowed <- aapa:::.compute_allowed_genotypes(2L, 2L)
  expect_false(allowed[1, "0"])
  expect_false(allowed[1, "1"])
  expect_true(allowed[1, "2"])

  # One parent 0, other 2 -> offspring must be 1
  allowed <- aapa:::.compute_allowed_genotypes(0L, 2L)
  expect_false(allowed[1, "0"])
  expect_true(allowed[1, "1"])
  expect_false(allowed[1, "2"])

  # Both heterozygous -> offspring can be 0, 1, or 2
  allowed <- aapa:::.compute_allowed_genotypes(1L, 1L)
  expect_true(all(allowed[1, ]))
})
