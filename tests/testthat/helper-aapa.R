# Load the aapa package for testing
library(aapa)

make_parents_object <- function(genotype, parents_df) {
  parents <- lapply(seq_len(nrow(parents_df)), function(i) {
    list(
      family_id = parents_df$family_id[i],
      sire_id = parents_df$sire_id[i],
      dam_id = parents_df$dam_id[i],
      sire_geno = genotype[parents_df$sire_id[i], ],
      dam_geno = genotype[parents_df$dam_id[i], ]
    )
  })

  names(parents) <- parents_df$family_id
  class(parents) <- "aapa_parents"
  parents
}

make_anchors_object <- function(genotype, anchors_df) {
  structure(
    anchors_df,
    class = c("aapa_anchors", "data.frame"),
    geno = genotype[anchors_df$individual_id, , drop = FALSE]
  )
}

compute_allowed_genotypes <- getFromNamespace(".compute_allowed_genotypes", "aapa")
