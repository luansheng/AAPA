# ============================================================================
# test-data-io.R — Unit tests for AAPA data input helpers
# ============================================================================

test_that("read_genotype returns a valid genotype matrix", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  write.csv(
    data.frame(
      id = c("IND1", "IND2"),
      SNP1 = c(0, 1),
      SNP2 = c(2, NA)
    ),
    tmp,
    row.names = FALSE
  )

  geno <- read_genotype(tmp)

  expect_true(is.matrix(geno))
  expect_identical(storage.mode(geno), "integer")
  expect_identical(rownames(geno), c("IND1", "IND2"))
  expect_identical(colnames(geno), c("SNP1", "SNP2"))
})

test_that("read_genotype rejects duplicate individual IDs", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  write.csv(
    data.frame(
      id = c("IND1", "IND1"),
      SNP1 = c(0, 1),
      SNP2 = c(2, 1)
    ),
    tmp,
    row.names = FALSE
  )

  expect_error(
    read_genotype(tmp),
    "unique, non-missing individual IDs"
  )
})

test_that("read_parents returns a valid aapa_parents object", {
  sim <- simulate_aapa_data(
    n_families = 3,
    n_snps = 50,
    n_offspring_per_family = 3,
    n_anchors_per_family = 1,
    n_unknown = 0,
    missing_rate = 0,
    error_rate = 0
  )

  con <- textConnection(
    utils::capture.output(write.csv(sim$parents, row.names = FALSE))
  )
  on.exit(close(con), add = TRUE)

  parents <- read_parents(con, sim$genotype)

  expect_s3_class(parents, "aapa_parents")
  expect_identical(names(parents), sim$parents$family_id)
  expect_true(all(vapply(parents, `[[`, character(1), "family_id") == names(parents)))
})

test_that("read_parents rejects duplicate family IDs", {
  sim <- simulate_aapa_data(
    n_families = 2,
    n_snps = 20,
    n_offspring_per_family = 2,
    n_anchors_per_family = 1,
    n_unknown = 0,
    missing_rate = 0,
    error_rate = 0
  )

  parents_df <- sim$parents
  parents_df$family_id[2] <- parents_df$family_id[1]
  con <- textConnection(
    utils::capture.output(write.csv(parents_df, row.names = FALSE))
  )
  on.exit(close(con), add = TRUE)

  expect_error(
    read_parents(con, sim$genotype),
    "unique, non-missing `family_id`"
  )
})

test_that("read_anchors returns a valid aapa_anchors object", {
  sim <- simulate_aapa_data(
    n_families = 3,
    n_snps = 50,
    n_offspring_per_family = 3,
    n_anchors_per_family = 2,
    n_unknown = 0,
    missing_rate = 0,
    error_rate = 0
  )

  con <- textConnection(
    utils::capture.output(write.csv(sim$anchors, row.names = FALSE))
  )
  on.exit(close(con), add = TRUE)

  anchors <- read_anchors(con, sim$genotype)

  expect_s3_class(anchors, "aapa_anchors")
  expect_identical(rownames(attr(anchors, "geno")), anchors$individual_id)
  expect_identical(colnames(attr(anchors, "geno")), colnames(sim$genotype))
})

test_that("read_anchors rejects duplicate anchor IDs", {
  sim <- simulate_aapa_data(
    n_families = 2,
    n_snps = 20,
    n_offspring_per_family = 2,
    n_anchors_per_family = 2,
    n_unknown = 0,
    missing_rate = 0,
    error_rate = 0
  )

  anchors_df <- sim$anchors
  anchors_df$individual_id[2] <- anchors_df$individual_id[1]
  con <- textConnection(
    utils::capture.output(write.csv(anchors_df, row.names = FALSE))
  )
  on.exit(close(con), add = TRUE)

  expect_error(
    read_anchors(con, sim$genotype),
    "unique, non-missing `individual_id`"
  )
})
