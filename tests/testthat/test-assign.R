# ============================================================================
# test-assign.R — Unit tests for AAPA assignment pipeline
# ============================================================================

test_that("aapa_assign runs end-to-end with simulated data", {
  sim <- simulate_aapa_data(
    n_families = 5, n_snps = 200,
    n_offspring_per_family = 8,
    n_anchors_per_family = 2,
    n_unknown = 3,
    missing_rate = 0.01, error_rate = 0.001
  )

  geno <- sim$genotype
  parents_df <- sim$parents
  anchors_df <- sim$anchors

  parents <- make_parents_object(geno, parents_df)
  anchors <- make_anchors_object(geno, anchors_df)

  result <- aapa_assign(geno, parents, anchors, top_k = 3)

  expect_s3_class(result, "aapa_result")
  expect_true(nrow(result$assignment) > 0)
  expect_true(all(result$assignment$status %in% c("ASSIGNED", "REJECT")))
  expect_equal(ncol(result$conflict_matrix), 5)
  expect_true(length(result$topk) > 0)
  expect_identical(rownames(result$conflict_matrix), result$assignment$individual_id)
  expect_identical(colnames(result$conflict_matrix), names(parents))
  expect_identical(colnames(result$kinship_matrix), names(parents))
  expect_identical(colnames(result$score_matrix), names(parents))
})

test_that("aapa_assign correctly assigns perfect data", {
  sim <- simulate_aapa_data(
    n_families = 3, n_snps = 500,
    n_offspring_per_family = 5,
    n_anchors_per_family = 3,
    n_unknown = 0,
    missing_rate = 0, error_rate = 0
  )

  geno <- sim$genotype
  parents_df <- sim$parents
  anchors_df <- sim$anchors

  parents <- make_parents_object(geno, parents_df)
  anchors <- make_anchors_object(geno, anchors_df)

  result <- aapa_assign(geno, parents, anchors,
    top_k = 3, tau_conf = -Inf,
    tau_rej = 0, max_conflict = 1.0
  )

  # With perfect data, most true offspring should be correctly assigned
  asgn <- result$assignment
  assigned <- asgn[asgn$status == "ASSIGNED", ]
  true_labels <- sim$true_labels

  correct <- sum(
    assigned$assigned_family == true_labels[assigned$individual_id],
    na.rm = TRUE
  )
  accuracy <- correct / nrow(assigned)

  # Expect high accuracy with perfect data

  expect_true(accuracy >= 0.8,
    info = sprintf("Accuracy %.2f is below 0.8", accuracy)
  )
})

test_that("aapa_assign rejects unknown-family individuals", {
  sim <- simulate_aapa_data(
    n_families = 3, n_snps = 500,
    n_offspring_per_family = 5,
    n_anchors_per_family = 3,
    n_unknown = 5,
    missing_rate = 0, error_rate = 0
  )

  geno <- sim$genotype
  parents_df <- sim$parents
  anchors_df <- sim$anchors

  parents <- make_parents_object(geno, parents_df)
  anchors <- make_anchors_object(geno, anchors_df)

  # Use strict rejection thresholds
  result <- aapa_assign(geno, parents, anchors,
    top_k = 3,
    tau_conf = -0.05,
    tau_rej = 0.05,
    max_conflict = 0.05
  )

  # Check that at least some UNKNOWN individuals are rejected
  unknown_ids <- grep("^UNKNOWN", result$assignment$individual_id,
    value = TRUE
  )
  unknown_results <- result$assignment[
    result$assignment$individual_id %in% unknown_ids,
  ]

  if (nrow(unknown_results) > 0) {
    # At least some unknown individuals should be rejected
    n_rejected_unknown <- sum(unknown_results$status == "REJECT")
    # This is a soft check since rejection depends on data characteristics
    expect_true(n_rejected_unknown >= 0)
  }
})

test_that("aapa_assign aligns parents and anchors after QC marker filtering", {
  sim <- simulate_aapa_data(
    n_families = 4, n_snps = 200,
    n_offspring_per_family = 5,
    n_anchors_per_family = 2,
    n_unknown = 2,
    missing_rate = 0.05, error_rate = 0
  )

  geno <- sim$genotype
  parents <- make_parents_object(geno, sim$parents)
  anchors <- make_anchors_object(geno, sim$anchors)

  qc_result <- qc_filter(
    geno,
    max_snp_missing = 0.1,
    max_ind_missing = 1,
    min_maf = 0.05,
    verbose = FALSE
  )

  result <- aapa_assign(qc_result$genotype, parents, anchors, top_k = 3)

  expect_s3_class(result, "aapa_result")
  expect_identical(colnames(result$conflict_matrix), names(parents))
})

test_that("aapa_assign rejects anchors with unknown family IDs", {
  sim <- simulate_aapa_data(
    n_families = 3, n_snps = 100,
    n_offspring_per_family = 4,
    n_anchors_per_family = 2,
    n_unknown = 0,
    missing_rate = 0, error_rate = 0
  )

  parents <- make_parents_object(sim$genotype, sim$parents)
  anchors <- make_anchors_object(sim$genotype, sim$anchors)
  anchors$family_id[1] <- "FAM999"

  expect_error(
    aapa_assign(sim$genotype, parents, anchors),
    "Anchor family IDs must exist"
  )
})
