# ============================================================================
# data_io.R — Data input helpers for AAPA
# ============================================================================

#' @import cli
#' @import checkmate
NULL

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
  dat <- utils::read.csv(file, sep = sep, header = header,
                         stringsAsFactors = FALSE, check.names = FALSE)
  ids <- as.character(dat[[1]])
  geno <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(geno) <- "integer"
  rownames(geno) <- ids

  # Validate dosage values
  valid_vals <- c(0L, 1L, 2L, NA_integer_)
  if (!all(geno %in% valid_vals)) {
    cli::cli_warn(c(
      "Genotype matrix contains values other than 0, 1, 2, or NA.",
      "i" = "Non-standard values will be treated as missing."
    ))
    geno[!geno %in% c(0L, 1L, 2L)] <- NA_integer_
  }
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
#' @param genotype_matrix A genotype matrix (from [read_genotype()]) that
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
  checkmate::assert_matrix(genotype_matrix, mode = "numeric")
  checkmate::assert_string(sep)

  if (is.character(file)) {
    cli::cli_alert_info("Reading parents file: {.file {file}}")
  }
  dat <- utils::read.csv(file, sep = sep, header = TRUE,
                         stringsAsFactors = FALSE)
  required_cols <- c("family_id", "sire_id", "dam_id")
  if (!all(required_cols %in% names(dat))) {
    cli::cli_abort(c(
      "Parents file is missing required columns.",
      "x" = "Required: {.field {required_cols}}",
      "i" = "Found: {.field {names(dat)}}"
    ))
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
  checkmate::assert_matrix(genotype_matrix, mode = "numeric")
  checkmate::assert_string(sep)

  if (is.character(file)) {
    cli::cli_alert_info("Reading anchors file: {.file {file}}")
  }
  dat <- utils::read.csv(file, sep = sep, header = TRUE,
                         stringsAsFactors = FALSE)
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
  n_fam <- length(unique(dat$family_id))
  cli::cli_alert_success(
    "Read {nrow(dat)} anchor{?s} from {n_fam} famil{?y/ies}"
  )
  structure(dat, class = c("aapa_anchors", "data.frame"),
            geno = anchor_geno)
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
                      nrow = n, ncol = length(p))
    allele2 <- matrix(stats::rbinom(n * length(p), 1, rep(p, each = n)),
                      nrow = n, ncol = length(p))
    allele1 + allele2
  }

  # Helper: simulate offspring from sire and dam genotypes
  sim_offspring <- function(sire_geno, dam_geno, n) {
    n_snps_local <- length(sire_geno)
    offspring <- matrix(0L, nrow = n, ncol = n_snps_local)
    for (j in seq_len(n_snps_local)) {
      # Sire allele
      s <- sire_geno[j]
      s_allele <- if (s == 0) rep(0L, n)
                  else if (s == 2) rep(1L, n)
                  else stats::rbinom(n, 1, 0.5)
      # Dam allele
      d <- dam_geno[j]
      d_allele <- if (d == 0) rep(0L, n)
                  else if (d == 2) rep(1L, n)
                  else stats::rbinom(n, 1, 0.5)
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
  all_ids <- character(0)
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

  n_total_ind <- nrow(geno_mat)
  cli::cli_alert_success(
    "Simulated {n_total_ind} individual{?s} x {n_snps} marker{?s}"
  )

  list(
    genotype = geno_mat,
    parents = parents_df,
    anchors = anchor_df,
    true_labels = true_labels
  )
}
