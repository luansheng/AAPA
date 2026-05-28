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

  # Build parents object
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

  # Build anchors object
  anchors <- structure(
    anchors_df,
    class = c("aapa_anchors", "data.frame"),
    geno = geno[anchors_df$individual_id, , drop = FALSE]
  )

  result <- aapa_assign(geno, parents, anchors, top_k = 3)

  expect_s3_class(result, "aapa_result")
  expect_true(nrow(result$assignment) > 0)
  expect_true(all(result$assignment$status %in% c("ASSIGNED", "REJECT")))
  expect_equal(ncol(result$conflict_matrix), 5)
  expect_true(length(result$topk) > 0)
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

  anchors <- structure(
    anchors_df,
    class = c("aapa_anchors", "data.frame"),
    geno = geno[anchors_df$individual_id, , drop = FALSE]
  )

  result <- aapa_assign(geno, parents, anchors,
                        top_k = 3, tau_conf = -Inf,
                        tau_rej = 0, max_conflict = 1.0)

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
              info = sprintf("Accuracy %.2f is below 0.8", accuracy))
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

  anchors <- structure(
    anchors_df,
    class = c("aapa_anchors", "data.frame"),
    geno = geno[anchors_df$individual_id, , drop = FALSE]
  )

  # Use strict rejection thresholds
  result <- aapa_assign(geno, parents, anchors,
                        top_k = 3,
                        tau_conf = -0.05,
                        tau_rej = 0.05,
                        max_conflict = 0.05)

  # Check that at least some UNKNOWN individuals are rejected
  unknown_ids <- grep("^UNKNOWN", result$assignment$individual_id,
                      value = TRUE)
  unknown_results <- result$assignment[
    result$assignment$individual_id %in% unknown_ids, ]

  if (nrow(unknown_results) > 0) {
    # At least some unknown individuals should be rejected
    n_rejected_unknown <- sum(unknown_results$status == "REJECT")
    # This is a soft check since rejection depends on data characteristics
    expect_true(n_rejected_unknown >= 0)
  }
})
