# Small example dataset for AAPA
#
# This script generates a minimal example dataset that can be
# used to test the AAPA package functions. Run this script to
# regenerate the CSV files in this directory.
#
# Generated files:
#   genotype.csv    - Genotype dosage matrix (0/1/2/NA)
#   parents.csv     - Candidate family parent table
#   anchors.csv     - Anchor individual table
#   true_labels.csv - True family assignments (for evaluation)

library(aapa)

sim <- simulate_aapa_data(
  n_families = 10,
  n_snps = 500,
  n_offspring_per_family = 10,
  n_anchors_per_family = 2,
  n_unknown = 5,
  missing_rate = 0.01,
  error_rate = 0.001,
  seed = 42
)

# Save genotype matrix
geno_df <- data.frame(individual_id = rownames(sim$genotype),
                      sim$genotype, check.names = FALSE)
write.csv(geno_df, "genotype.csv", row.names = FALSE)

# Save parents table
write.csv(sim$parents, "parents.csv", row.names = FALSE)

# Save anchors table
write.csv(sim$anchors, "anchors.csv", row.names = FALSE)

# Save true labels
labels_df <- data.frame(
  individual_id = names(sim$true_labels),
  true_family = sim$true_labels,
  stringsAsFactors = FALSE
)
write.csv(labels_df, "true_labels.csv", row.names = FALSE)

cat("Example data generated successfully.\n")
cat(sprintf("  Individuals: %d\n", nrow(sim$genotype)))
cat(sprintf("  SNPs: %d\n", ncol(sim$genotype)))
cat(sprintf("  Families: %d\n", nrow(sim$parents)))
cat(sprintf("  Anchors: %d\n", nrow(sim$anchors)))
