# ============================================================================
# data_io.R — Data input helpers for AAPA
# ============================================================================

#' @import cli
#' @import checkmate
#' @importFrom data.table fread
NULL

#' Validate genotype matrix contract
#'
#' Check that a genotype matrix satisfies the AAPA object contract:
#' unique sample IDs, unique marker IDs, and dosage encoding 0/1/2/NA.
#'
#' @param genotype A genotype matrix.
#' @param arg Argument name used in error messages.
#' @return The input matrix, invisibly.
#' @keywords internal
.validate_genotype_matrix <- function(genotype, arg = "genotype") {
  checkmate::assert_matrix(genotype,
    mode = "numeric", min.rows = 1,
    min.cols = 1
  )

  sample_ids <- rownames(genotype)
  marker_ids <- colnames(genotype)

  if (is.null(sample_ids) || anyNA(sample_ids) || any(sample_ids == "") || anyDuplicated(sample_ids)) {
    cli::cli_abort(c(
      "{.arg {arg}} must have unique, non-missing row names.",
      "x" = "Sample IDs are required for all downstream alignment steps."
    ))
  }

  if (is.null(marker_ids) || anyNA(marker_ids) || any(marker_ids == "") || anyDuplicated(marker_ids)) {
    cli::cli_abort(c(
      "{.arg {arg}} must have unique, non-missing column names.",
      "x" = "Marker IDs are required for all downstream alignment steps."
    ))
  }

  valid_vals <- is.na(genotype) | genotype %in% c(0L, 1L, 2L)
  if (!all(valid_vals)) {
    cli::cli_abort(c(
      "{.arg {arg}} must use dosage encoding 0/1/2/NA.",
      "x" = "Found genotype values outside the supported encoding."
    ))
  }

  invisible(genotype)
}

#' Validate parent identifiers and family naming
#'
#' @param family A single parent entry.
#' @param family_name Expected family name from `names(parents)`.
#' @return The input family entry, invisibly.
#' @keywords internal
.validate_parent_ids <- function(family, family_name) {
  required_fields <- c("family_id", "sire_id", "dam_id", "sire_geno", "dam_geno")

  if (!all(required_fields %in% names(family))) {
    cli::cli_abort(c(
      "Each `parents` entry must contain the required fields.",
      "x" = sprintf("Family %s is missing one or more required fields.", family_name)
    ))
  }

  if (!identical(as.character(family$family_id), family_name)) {
    cli::cli_abort(c(
      "`names(parents)` must match each entry's `family_id`.",
      "x" = sprintf("Found mismatch for family %s.", family_name)
    ))
  }

  if (!checkmate::test_string(family$sire_id) || !checkmate::test_string(family$dam_id)) {
    cli::cli_abort("Each `parents` entry must contain scalar `sire_id` and `dam_id` strings.")
  }

  invisible(family)
}

#' Validate parent genotype vector types and lengths
#'
#' @param family A single parent entry.
#' @param family_name Expected family name from `names(parents)`.
#' @return The input family entry, invisibly.
#' @keywords internal
.validate_parent_geno_vectors <- function(family, family_name) {
  invalid_sire_geno <- !is.atomic(family$sire_geno) || is.matrix(family$sire_geno)
  invalid_dam_geno <- !is.atomic(family$dam_geno) || is.matrix(family$dam_geno)

  if (invalid_sire_geno || invalid_dam_geno) {
    cli::cli_abort("`sire_geno` and `dam_geno` must be named atomic vectors.")
  }

  if (!identical(length(family$sire_geno), length(family$dam_geno))) {
    cli::cli_abort(c(
      "Parent genotype lengths must match within each family.",
      "x" = sprintf("Family %s has incompatible sire/dam genotype lengths.", family_name)
    ))
  }

  invisible(family)
}

#' Validate parent marker names
#'
#' @param family A single parent entry.
#' @param family_name Expected family name from `names(parents)`.
#' @return The input family entry, invisibly.
#' @keywords internal
.validate_parent_marker_names <- function(family, family_name) {
  sire_markers <- names(family$sire_geno)
  dam_markers <- names(family$dam_geno)

  invalid_sire_markers <- is.null(sire_markers) || anyDuplicated(sire_markers) || anyNA(sire_markers)
  invalid_dam_markers <- is.null(dam_markers) || anyDuplicated(dam_markers) || anyNA(dam_markers)

  if (invalid_sire_markers || invalid_dam_markers) {
    cli::cli_abort(c(
      "Parent genotype vectors must carry unique marker names.",
      "x" = sprintf(
        "Family %s is missing marker names on `sire_geno` or `dam_geno`.",
        family_name
      )
    ))
  }

  invisible(family)
}

#' Validate parent genotype vectors
#'
#' @param family A single parent entry.
#' @param family_name Expected family name from `names(parents)`.
#' @return The input family entry, invisibly.
#' @keywords internal
.validate_parent_genotypes <- function(family, family_name) {
  .validate_parent_geno_vectors(family, family_name)
  .validate_parent_marker_names(family, family_name)

  invisible(family)
}

#' Validate a single parent entry
#'
#' Check that one family entry in an `aapa_parents` object satisfies the
#' required field, naming, and marker constraints.
#'
#' @param family A single parent entry.
#' @param family_name Expected family name from `names(parents)`.
#' @return The input family entry, invisibly.
#' @keywords internal
.validate_parent_entry <- function(family, family_name) {
  .validate_parent_ids(family, family_name)
  .validate_parent_genotypes(family, family_name)

  invisible(family)
}

#' Validate parents object contract
#'
#' Check that a parent list contains the required fields, stable family names,
#' and marker-named parent genotype vectors.
#'
#' @param parents An `aapa_parents` object or compatible named list.
#' @return The input parents object, invisibly.
#' @keywords internal
.validate_parents_object <- function(parents) {
  checkmate::assert_list(parents, min.len = 1)

  parent_names <- names(parents)
  if (is.null(parent_names) || anyNA(parent_names) || any(parent_names == "") || anyDuplicated(parent_names)) {
    cli::cli_abort(c(
      "`parents` must be a named list with unique family IDs.",
      "x" = "Use `names(parents) <- family_id`."
    ))
  }

  for (ii in seq_along(parents)) {
    .validate_parent_entry(parents[[ii]], parent_names[[ii]])
  }

  invisible(parents)
}

#' Align parent genotype vectors to genotype markers
#'
#' Reorder each family's `sire_geno` and `dam_geno` vectors to the supplied
#' marker order after checking that all markers are present.
#'
#' @param parents An `aapa_parents` object or compatible named list.
#' @param marker_ids Character vector of marker IDs defining the target order.
#' @return A parents object aligned to `marker_ids`.
#' @keywords internal
.align_parents_to_markers <- function(parents, marker_ids) {
  .validate_parents_object(parents)

  aligned <- lapply(seq_along(parents), function(ii) {
    family <- parents[[ii]]

    missing_sire <- setdiff(marker_ids, names(family$sire_geno))
    missing_dam <- setdiff(marker_ids, names(family$dam_geno))

    if (length(missing_sire) > 0 || length(missing_dam) > 0) {
      cli::cli_abort(c(
        "Parent genotypes must contain all genotype markers.",
        "x" = sprintf(
          "Family %s is missing marker names required for alignment.",
          names(parents)[[ii]]
        )
      ))
    }

    family$sire_geno <- family$sire_geno[marker_ids]
    family$dam_geno <- family$dam_geno[marker_ids]
    family
  })

  names(aligned) <- names(parents)
  structure(aligned, class = class(parents))
}

#' Validate anchors object contract
#'
#' Check that an anchors object has required columns, stable IDs, and an
#' aligned genotype matrix in `attr(anchors, "geno")`.
#'
#' @param anchors An `aapa_anchors` object.
#' @param family_ids Optional candidate family IDs used to validate anchor
#'   family membership.
#' @return The input anchors object, invisibly.
#' @keywords internal
.validate_anchors_object <- function(anchors, family_ids = NULL) {
  checkmate::assert_class(anchors, "aapa_anchors")

  required_cols <- c("individual_id", "family_id")
  if (!all(required_cols %in% names(anchors))) {
    cli::cli_abort(c(
      "`anchors` must contain required columns.",
      "x" = "Required columns: {.field {required_cols}}"
    ))
  }

  if (anyNA(anchors$individual_id) || any(anchors$individual_id == "") || anyDuplicated(anchors$individual_id)) {
    cli::cli_abort("`anchors$individual_id` must be unique and non-missing.")
  }

  if (anyNA(anchors$family_id) || any(anchors$family_id == "")) {
    cli::cli_abort("`anchors$family_id` must be non-missing.")
  }

  if (!is.null(family_ids) && !all(anchors$family_id %in% family_ids)) {
    missing_families <- unique(setdiff(anchors$family_id, family_ids))
    cli::cli_abort(c(
      "Anchor family IDs must exist in the candidate family set.",
      "x" = sprintf("Unknown anchor families: %s", paste(missing_families, collapse = ", "))
    ))
  }

  anchor_geno <- attr(anchors, "geno")
  if (!is.matrix(anchor_geno)) {
    cli::cli_abort("`anchors` must carry a genotype matrix in `attr(anchors, \"geno\")`.")
  }

  if (is.null(rownames(anchor_geno)) || anyDuplicated(rownames(anchor_geno)) || anyNA(rownames(anchor_geno))) {
    cli::cli_abort("`attr(anchors, \"geno\")` must have unique row names matching `individual_id`.")
  }

  invisible(anchors)
}

#' Align anchor genotype matrix to genotype markers
#'
#' Reorder the anchor genotype submatrix to the supplied marker order after
#' checking that all required anchor IDs and marker IDs are present.
#'
#' @param anchors An `aapa_anchors` object.
#' @param marker_ids Character vector of marker IDs defining the target order.
#' @return An anchors object aligned to `marker_ids`.
#' @keywords internal
.align_anchors_to_markers <- function(anchors, marker_ids) {
  anchor_geno <- attr(anchors, "geno")

  missing_anchor_ids <- setdiff(anchors$individual_id, rownames(anchor_geno))
  if (length(missing_anchor_ids) > 0) {
    cli::cli_abort(c(
      "Anchor genotype matrix must contain all anchor individuals.",
      "x" = sprintf("Missing anchor IDs: %s", paste(missing_anchor_ids, collapse = ", "))
    ))
  }

  if (is.null(colnames(anchor_geno)) || anyDuplicated(colnames(anchor_geno)) || anyNA(colnames(anchor_geno))) {
    cli::cli_abort("`attr(anchors, \"geno\")` must have unique marker names.")
  }

  missing_markers <- setdiff(marker_ids, colnames(anchor_geno))
  if (length(missing_markers) > 0) {
    cli::cli_abort(c(
      "Anchor genotype matrix must contain all genotype markers.",
      "x" = sprintf("Missing markers: %s", paste(missing_markers, collapse = ", "))
    ))
  }

  attr(anchors, "geno") <- anchor_geno[anchors$individual_id, marker_ids, drop = FALSE]
  anchors
}

#' Read genotype dosage matrix
#'
#' Read a genotype matrix in dosage format (0/1/2, NA for missing).
#' Rows are individuals, columns are SNP markers.
#'
#' @param file Path to a CSV/TSV file. First column is individual ID,
#'   remaining columns are SNP dosage values (0, 1, 2, or NA).
#' @param sep Field separator (default: comma).
#' @param header Logical; does the file contain a header row? Default TRUE.
#' @return A numeric matrix with rownames = individual IDs,
#'   colnames = marker names.
#' @family data-io
#' @examples
#' # Create a temporary genotype file
#' tmp <- tempfile(fileext = ".csv")
#' geno_data <- data.frame(
#'   id = c("IND1", "IND2"),
#'   SNP1 = c(0, 1),
#'   SNP2 = c(2, 1)
#' )
#' write.csv(geno_data, tmp, row.names = FALSE)
#' geno <- read_genotype(tmp)
#' unlink(tmp)
#' @export
read_genotype <- function(file, sep = ",", header = TRUE) {
  checkmate::assert_string(file)
  checkmate::assert_file_exists(file)
  checkmate::assert_string(sep)
  checkmate::assert_flag(header)

  cli::cli_alert_info("Reading genotype file: {.file {file}}")
  dat <- data.table::fread(file,
    sep = sep, header = header,
    stringsAsFactors = FALSE, check.names = FALSE,
    data.table = FALSE
  )
  ids <- as.character(dat[[1]])
  geno <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(geno) <- "integer"
  rownames(geno) <- ids

  if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
    cli::cli_abort("Genotype file must contain unique, non-missing individual IDs in the first column.")
  }

  if (is.null(colnames(geno)) || anyNA(colnames(geno)) || any(colnames(geno) == "") || anyDuplicated(colnames(geno))) {
    cli::cli_abort("Genotype file must contain unique, non-missing marker names.")
  }

  # Validate dosage values
  valid_vals <- c(0L, 1L, 2L, NA_integer_)
  if (!all(geno %in% valid_vals)) {
    cli::cli_warn(c(
      "Genotype matrix contains values other than 0, 1, 2, or NA.",
      "i" = "Non-standard values will be treated as missing."
    ))
    geno[!geno %in% c(0L, 1L, 2L)] <- NA_integer_
  }

  .validate_genotype_matrix(geno)

  cli::cli_alert_success(
    "Read {nrow(geno)} individual{?s} x {ncol(geno)} marker{?s}"
  )
  geno
}

#' Read parents table
#'
#' Read a table of candidate families with parent genotypes.
#'
#' @param file Path to a CSV/TSV file with columns: family_id, sire_id,
#'   dam_id.
#' @param genotype_matrix A genotype matrix (from \code{read_genotype()}) that
#'   includes sire and dam genotypes.
#' @param sep Field separator (default: comma).
#' @return A list of class \code{aapa_parents}, where each element is a
#'   list with \code{family_id}, \code{sire_id}, \code{dam_id},
#'   \code{sire_geno}, \code{dam_geno}.
#' @family data-io
#' @export
read_parents <- function(file, genotype_matrix, sep = ",") {
  checkmate::assert(
    checkmate::check_string(file),
    checkmate::check_class(file, "connection")
  )
  if (is.character(file)) checkmate::assert_file_exists(file)
  .validate_genotype_matrix(genotype_matrix, arg = "genotype_matrix")
  checkmate::assert_string(sep)

  if (is.character(file)) {
    cli::cli_alert_info("Reading parents file: {.file {file}}")
  }
  dat <- if (is.character(file)) {
    data.table::fread(file,
      sep = sep, header = TRUE,
      stringsAsFactors = FALSE, data.table = FALSE
    )
  } else {
    utils::read.csv(file,
      sep = sep, header = TRUE,
      stringsAsFactors = FALSE
    )
  }
  required_cols <- c("family_id", "sire_id", "dam_id")
  if (!all(required_cols %in% names(dat))) {
    cli::cli_abort(c(
      "Parents file is missing required columns.",
      "x" = "Required: {.field {required_cols}}",
      "i" = "Found: {.field {names(dat)}}"
    ))
  }

  dat$family_id <- as.character(dat$family_id)
  dat$sire_id <- as.character(dat$sire_id)
  dat$dam_id <- as.character(dat$dam_id)

  if (anyNA(dat$family_id) || any(dat$family_id == "") || anyDuplicated(dat$family_id)) {
    cli::cli_abort("Parents file must contain unique, non-missing `family_id` values.")
  }

  families <- lapply(seq_len(nrow(dat)), function(i) {
    fid <- as.character(dat$family_id[i])
    sid <- as.character(dat$sire_id[i])
    did <- as.character(dat$dam_id[i])

    if (!sid %in% rownames(genotype_matrix)) {
      cli::cli_abort("Sire {.val {sid}} not found in genotype matrix.")
    }
    if (!did %in% rownames(genotype_matrix)) {
      cli::cli_abort("Dam {.val {did}} not found in genotype matrix.")
    }

    list(
      family_id = fid,
      sire_id = sid,
      dam_id = did,
      sire_geno = genotype_matrix[sid, ],
      dam_geno = genotype_matrix[did, ]
    )
  })
  names(families) <- vapply(families, `[[`, character(1), "family_id")
  .validate_parents_object(families)
  cli::cli_alert_success("Read {length(families)} candidate famil{?y/ies}")
  structure(families, class = "aapa_parents")
}

#' Read anchors table
#'
#' Read a table of anchor individuals with known family assignments.
#'
#' @param file Path to a CSV/TSV file with columns: individual_id,
#'   family_id, and optionally weight.
#' @param genotype_matrix A genotype matrix that includes anchor genotypes.
#' @param sep Field separator (default: comma).
#' @return A data.frame of class \code{aapa_anchors} with columns:
#'   individual_id, family_id, weight, and a \code{geno} attribute
#'   containing the anchor genotype sub-matrix.
#' @family data-io
#' @export
read_anchors <- function(file, genotype_matrix, sep = ",") {
  checkmate::assert(
    checkmate::check_string(file),
    checkmate::check_class(file, "connection")
  )
  if (is.character(file)) checkmate::assert_file_exists(file)
  .validate_genotype_matrix(genotype_matrix, arg = "genotype_matrix")
  checkmate::assert_string(sep)

  if (is.character(file)) {
    cli::cli_alert_info("Reading anchors file: {.file {file}}")
  }
  dat <- if (is.character(file)) {
    data.table::fread(file,
      sep = sep, header = TRUE,
      stringsAsFactors = FALSE, data.table = FALSE
    )
  } else {
    utils::read.csv(file,
      sep = sep, header = TRUE,
      stringsAsFactors = FALSE
    )
  }
  required_cols <- c("individual_id", "family_id")
  if (!all(required_cols %in% names(dat))) {
    cli::cli_abort(c(
      "Anchors file is missing required columns.",
      "x" = "Required: {.field {required_cols}}",
      "i" = "Found: {.field {names(dat)}}"
    ))
  }

  dat$individual_id <- as.character(dat$individual_id)
  dat$family_id <- as.character(dat$family_id)

  if (anyNA(dat$individual_id) || any(dat$individual_id == "") || anyDuplicated(dat$individual_id)) {
    cli::cli_abort("Anchors file must contain unique, non-missing `individual_id` values.")
  }

  if (!"weight" %in% names(dat)) {
    dat$weight <- 1.0
  }

  missing_ids <- setdiff(dat$individual_id, rownames(genotype_matrix))
  if (length(missing_ids) > 0) {
    cli::cli_abort(c(
      "Anchor individuals not found in genotype matrix.",
      "x" = "Missing: {.val {missing_ids}}"
    ))
  }

  anchor_geno <- genotype_matrix[dat$individual_id, , drop = FALSE]
  cli::cli_alert_success(
    sprintf(
      "Read %d anchors from %d families",
      nrow(dat),
      length(unique(dat$family_id))
    )
  )
  anchors <- structure(dat,
    class = c("aapa_anchors", "data.frame"),
    geno = anchor_geno
  )
  .validate_anchors_object(anchors)
  anchors
}

#' Simulate AAPA example data
#'
#' Generate a small simulated dataset for testing and demonstration.
#'
#' @param n_families Number of candidate families (default: 10).
#' @param n_snps Number of SNP markers (default: 500).
#' @param n_offspring_per_family Number of offspring per family (default: 10).
#' @param n_anchors_per_family Number of anchors per family (default: 2).
#' @param n_unknown Number of unknown-family test individuals (default: 5).
#' @param missing_rate Proportion of missing genotypes (default: 0.01).
#' @param error_rate Genotyping error rate (default: 0.001).
#' @param seed Random seed for reproducibility.
#' @return A list with components: \code{genotype} (matrix), \code{parents}
#'   (data.frame), \code{anchors} (data.frame), \code{true_labels}
#'   (named character vector).
#' @family data-io
#' @examples
#' sim <- simulate_aapa_data(n_families = 3, n_snps = 100)
#' str(sim, max.level = 1)
#' @export
simulate_aapa_data <- function(n_families = 10, n_snps = 500,
                               n_offspring_per_family = 10,
                               n_anchors_per_family = 2,
                               n_unknown = 5,
                               missing_rate = 0.01,
                               error_rate = 0.001,
                               seed = 42) {
  checkmate::assert_count(n_families, positive = TRUE)
  checkmate::assert_count(n_snps, positive = TRUE)
  checkmate::assert_count(n_offspring_per_family, positive = TRUE)
  checkmate::assert_count(n_anchors_per_family)
  checkmate::assert_count(n_unknown)
  checkmate::assert_number(missing_rate, lower = 0, upper = 1)
  checkmate::assert_number(error_rate, lower = 0, upper = 1)
  checkmate::assert_int(seed)

  cli::cli_alert_info(
    "Simulating data: {n_families} famil{?y/ies}, {n_snps} SNP{?s}"
  )
  set.seed(seed)

  # Allele frequencies
  p <- stats::runif(n_snps, 0.1, 0.9)

  # Helper: simulate diploid genotype from allele freq
  sim_geno <- function(n, p) {
    allele1 <- matrix(stats::rbinom(n * length(p), 1, rep(p, each = n)),
      nrow = n, ncol = length(p)
    )
    allele2 <- matrix(stats::rbinom(n * length(p), 1, rep(p, each = n)),
      nrow = n, ncol = length(p)
    )
    allele1 + allele2
  }

  # Helper: simulate offspring from sire and dam genotypes
  sim_offspring <- function(sire_geno, dam_geno, n) {
    n_snps_local <- length(sire_geno)
    offspring <- matrix(0L, nrow = n, ncol = n_snps_local)
    for (j in seq_len(n_snps_local)) {
      # Sire allele
      s <- sire_geno[j]
      s_allele <- if (s == 0) {
        rep(0L, n)
      } else if (s == 2) {
        rep(1L, n)
      } else {
        stats::rbinom(n, 1, 0.5)
      }
      # Dam allele
      d <- dam_geno[j]
      d_allele <- if (d == 0) {
        rep(0L, n)
      } else if (d == 2) {
        rep(1L, n)
      } else {
        stats::rbinom(n, 1, 0.5)
      }
      offspring[, j] <- s_allele + d_allele
    }
    offspring
  }

  snp_names <- paste0("SNP", seq_len(n_snps))

  # Generate parents
  parents_df <- data.frame(
    family_id = character(0),
    sire_id = character(0),
    dam_id = character(0),
    stringsAsFactors = FALSE
  )
  all_geno <- list()
  anchor_df <- data.frame(
    individual_id = character(0),
    family_id = character(0),
    weight = numeric(0),
    stringsAsFactors = FALSE
  )
  true_labels <- character(0)

  for (f in seq_len(n_families)) {
    fid <- paste0("FAM", sprintf("%03d", f))
    sid <- paste0("SIRE", sprintf("%03d", f))
    did <- paste0("DAM", sprintf("%03d", f))

    # Parent genotypes
    sire_g <- sim_geno(1, p)[1, ]
    dam_g <- sim_geno(1, p)[1, ]

    all_geno[[sid]] <- sire_g
    all_geno[[did]] <- dam_g

    parents_df <- rbind(parents_df, data.frame(
      family_id = fid, sire_id = sid, dam_id = did,
      stringsAsFactors = FALSE
    ))

    # Offspring (includes anchors + test individuals)
    n_total <- n_anchors_per_family + n_offspring_per_family
    offspring_g <- sim_offspring(sire_g, dam_g, n_total)

    for (k in seq_len(n_total)) {
      if (k <= n_anchors_per_family) {
        oid <- paste0(fid, "_ANC", k)
        anchor_df <- rbind(anchor_df, data.frame(
          individual_id = oid, family_id = fid, weight = 1.0,
          stringsAsFactors = FALSE
        ))
      } else {
        oid <- paste0(fid, "_OFF", k - n_anchors_per_family)
      }
      all_geno[[oid]] <- offspring_g[k, ]
      true_labels[oid] <- fid
    }
  }

  # Unknown-family individuals (generated from random allele freqs)
  if (n_unknown > 0) {
    unknown_g <- sim_geno(n_unknown, p)
    for (u in seq_len(n_unknown)) {
      uid <- paste0("UNKNOWN", sprintf("%03d", u))
      all_geno[[uid]] <- unknown_g[u, ]
      true_labels[uid] <- "UNKNOWN"
    }
  }

  # Assemble genotype matrix
  geno_mat <- do.call(rbind, all_geno)
  colnames(geno_mat) <- snp_names
  storage.mode(geno_mat) <- "integer"

  # Introduce missing values
  if (missing_rate > 0) {
    n_missing <- round(length(geno_mat) * missing_rate)
    miss_idx <- sample(length(geno_mat), n_missing)
    geno_mat[miss_idx] <- NA_integer_
  }

  # Introduce genotyping errors
  if (error_rate > 0) {
    n_errors <- round(length(geno_mat) * error_rate)
    err_idx <- sample(length(geno_mat), n_errors)
    # Flip to a different dosage value
    for (idx in err_idx) {
      if (!is.na(geno_mat[idx])) {
        current <- geno_mat[idx]
        geno_mat[idx] <- sample(setdiff(c(0L, 1L, 2L), current), 1)
      }
    }
  }

  cli::cli_alert_success(
    sprintf("Simulated %d individuals x %d markers", nrow(geno_mat), n_snps)
  )

  list(
    genotype = geno_mat,
    parents = parents_df,
    anchors = anchor_df,
    true_labels = true_labels
  )
}
