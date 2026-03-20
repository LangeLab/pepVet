.required_digest_columns <- c(
  "protein_id",
  "peptide",
  "start",
  "end",
  "length",
  "missed_cleavages"
)

.validate_digest_result <- function(digest_result, arg_name = "digest_result") {
  if (!inherits(digest_result, "data.frame")) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        "} must be a digest tibble or data frame with pepVet digest columns."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  missing_columns <- setdiff(.required_digest_columns, names(digest_result))

  if (length(missing_columns) > 0L) {
    cli::cli_abort(
      c(
        paste0(
          "{.arg ",
          arg_name,
          "} is missing required digest columns."
        ),
        "i" = paste("Missing columns:", paste(missing_columns, collapse = ", "))
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  if (nrow(digest_result) == 0L) {
    cli::cli_abort(
      paste0("{.arg ", arg_name, "} must contain at least one peptide row."),
      class = "pepvet_error_invalid_digest"
    )
  }

  digest_result <- digest_result[, .required_digest_columns, drop = FALSE]

  if (
    !is.character(digest_result$protein_id) ||
      !is.character(digest_result$peptide)
  ) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        paste(
          "} must store {.field protein_id} and {.field peptide}",
          "as character columns."
        )
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  numeric_columns <- c("start", "end", "length", "missed_cleavages")

  if (!all(vapply(digest_result[numeric_columns], is.numeric, logical(1)))) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        "} must store coordinate and count columns as numeric values."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  if (anyNA(digest_result)) {
    cli::cli_abort(
      paste0("{.arg ", arg_name, "} must not contain missing values."),
      class = "pepvet_error_invalid_digest"
    )
  }

  peptide_lengths <- nchar(digest_result$peptide, type = "chars")
  row_lengths <- as.integer(digest_result$length)
  row_starts <- as.integer(digest_result$start)
  row_ends <- as.integer(digest_result$end)
  row_missed <- as.integer(digest_result$missed_cleavages)

  if (any(row_starts < 1L) || any(row_ends < row_starts)) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        "} contains invalid start/end coordinates."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  if (
    any(row_lengths != peptide_lengths) ||
      any(row_lengths != (row_ends - row_starts + 1L))
  ) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        "} must have internally consistent peptide lengths."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  if (any(row_missed < 0L)) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        "} must not contain negative missed-cleavage counts."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  digest_result$start <- row_starts
  digest_result$end <- row_ends
  digest_result$length <- row_lengths
  digest_result$missed_cleavages <- row_missed
  tibble::as_tibble(digest_result)
}

.valid_length_mask <- function(protein_digest) {
  protein_digest$length >= 7L & protein_digest$length <= 25L
}

.extract_valid_digest <- function(protein_digest) {
  protein_digest[.valid_length_mask(protein_digest), , drop = FALSE]
}

.score_length <- function(protein_digest) {
  sum(.valid_length_mask(protein_digest)) / nrow(protein_digest)
}

.score_coverage <- function(protein_digest) {
  valid_digest <- .extract_valid_digest(protein_digest)

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  protein_length <- max(protein_digest$end)
  covered_ranges <- IRanges::reduce(
    IRanges::IRanges(valid_digest$start, valid_digest$end)
  )

  sum(IRanges::width(covered_ranges)) / protein_length
}

.score_count <- function(protein_digest) {
  valid_count <- nrow(.extract_valid_digest(protein_digest))
  protein_length <- max(protein_digest$end)
  expected_count <- protein_length / 12

  min(valid_count / expected_count, 1)
}

.score_hydro <- function(protein_digest) {
  valid_digest <- .extract_valid_digest(protein_digest)

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  gravy_values <- vapply(valid_digest$peptide, .calculate_gravy, numeric(1))
  tolerance <- sqrt(.Machine$double.eps)

  mean(
    gravy_values >= (-1.0 - tolerance) &
      gravy_values <= (0.6 + tolerance)
  )
}

.score_charge <- function(protein_digest) {
  valid_digest <- .extract_valid_digest(protein_digest)

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  aa_properties <- .get_aa_properties()
  basic_residues <- aa_properties$amino_acid[aa_properties$is_basic]
  chargeable <- vapply(
    valid_digest$peptide,
    function(peptide) {
      internal_residues <- substr(peptide, 1L, nchar(peptide) - 1L)
      any(
        strsplit(internal_residues, split = "", fixed = TRUE)[[1]] %in%
          basic_residues
      )
    },
    logical(1)
  )

  mean(chargeable)
}

.score_unique <- function(protein_digest, proteome_index) {
  valid_digest <- .extract_valid_digest(protein_digest)

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  protein_id <- unique(protein_digest$protein_id)

  if (length(protein_id) != 1L) {
    cli::cli_abort(
      "{.arg protein_digest} must contain rows from exactly one protein.",
      class = "pepvet_error_invalid_digest"
    )
  }

  is_unique <- vapply(
    valid_digest$peptide,
    function(peptide) {
      indexed_proteins <- get0(
        peptide,
        envir = proteome_index,
        inherits = FALSE
      )

      if (is.null(indexed_proteins)) {
        return(TRUE)
      }

      length(indexed_proteins) == 1L && identical(indexed_proteins, protein_id)
    },
    logical(1)
  )

  mean(is_unique)
}

.classify_verdict <- function(composite_score) {
  ifelse(
    composite_score >= 0.7,
    "Good",
    ifelse(composite_score >= 0.4, "Moderate", "Poor")
  )
}

.score_components <- function(protein_digest, proteome_index = NULL) {
  component_scores <- c(
    S_length = .score_length(protein_digest),
    S_coverage = .score_coverage(protein_digest),
    S_count = .score_count(protein_digest),
    S_hydro = .score_hydro(protein_digest),
    S_charge = .score_charge(protein_digest)
  )

  if (!is.null(proteome_index)) {
    component_scores <- c(
      component_scores,
      S_unique = .score_unique(protein_digest, proteome_index)
    )
  }

  component_scores
}

#' Score Digested Peptides for Proteomics Suitability
#'
#' `score_peptides()` summarizes a pepVet digest tibble into per-protein scoring
#' components and a weighted composite suitability score.
#'
#' @param digest_result A digest tibble produced by [digest_protein()] or an
#'   equivalent table containing the columns `protein_id`, `peptide`, `start`,
#'   `end`, `length`, and `missed_cleavages`.
#' @param proteome Optional digest tibble representing the comparison proteome
#'   used for peptide uniqueness scoring. When omitted, `S_unique` is excluded
#'   and protein-only default weights are used.
#' @param weights Optional numeric weight vector. In protein-only mode the
#'   default weights are `c(S_length = 0.25, S_coverage = 0.25, S_count = 0.20,
#'   S_hydro = 0.15, S_charge = 0.15)`. In proteome-aware mode the default
#'   weights are `c(S_length = 0.20, S_coverage = 0.20, S_count = 0.15,
#'   S_hydro = 0.15, S_charge = 0.10, S_unique = 0.20)`.
#'
#' @return A tibble with one row per `protein_id` and the component score
#'   columns `S_length`, `S_coverage`, `S_count`, `S_hydro`, `S_charge`,
#'   optional `S_unique`, plus `composite_score` and `verdict`.
#'
#' @details Valid peptides are defined as peptides with lengths between 7 and
#'   25 residues inclusive. Coverage is computed from valid peptide coordinates
#'   with overlapping intervals reduced before coverage is summed. Composite
#'   verdicts are classified as `Good` for scores >= 0.7, `Moderate` for
#'   scores >= 0.4, and `Poor` otherwise.
#'
#' @examples
#' digest_result <- digest_protein("MKWVTFISLLFLFSSAYSR")
#' score_peptides(digest_result)
#'
#' target_and_background <- digest_protein(
#'   c(target = "AAAAAAARAAAAAAAK", background = "AAAAAAARGGGGGGGK")
#' )
#' score_peptides(
#'   target_and_background[target_and_background$protein_id == "target", ],
#'   proteome = target_and_background
#' )
#' @export
score_peptides <- function(digest_result, proteome = NULL, weights = NULL) {
  validated_digest <- .validate_digest_result(digest_result)
  has_proteome <- !is.null(proteome)
  proteome_index <- NULL

  if (has_proteome) {
    validated_proteome <- .validate_digest_result(
      proteome,
      arg_name = "proteome"
    )
    proteome_index <- .build_proteome_index(validated_proteome)
  }

  normalized_weights <- .validate_weights(weights, has_proteome)
  protein_levels <- unique(validated_digest$protein_id)
  protein_groups <- split(
    seq_len(nrow(validated_digest)),
    factor(validated_digest$protein_id, levels = protein_levels)
  )

  scored_rows <- lapply(
    protein_groups,
    function(row_indices) {
      protein_digest <- validated_digest[row_indices, , drop = FALSE]
      component_scores <- .score_components(protein_digest, proteome_index)
      weighted_components <- component_scores[names(normalized_weights)]
      composite_score <- sum(weighted_components * normalized_weights)

      tibble::as_tibble(
        c(
          list(protein_id = protein_digest$protein_id[[1]]),
          as.list(component_scores),
          list(
            composite_score = composite_score,
            verdict = .classify_verdict(composite_score)
          )
        )
      )
    }
  )

  tibble::as_tibble(do.call(rbind, scored_rows))
}
