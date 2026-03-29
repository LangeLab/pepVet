.cleavage_ranges <- get("cleavageRanges", envir = asNamespace("cleaver"))

.cleavage_annotation_trypsin_enzymes <- c(
  "trypsin",
  "trypsin-high",
  "trypsin-low",
  "trypsin-simple"
)

.supports_cleavage_efficiency_annotations <- function(enzyme) {
  enzyme %in% .cleavage_annotation_trypsin_enzymes
}

.annotate_trypsin_cleavage_sites <- function(sequence) {
  residues <- strsplit(sequence, split = "", fixed = TRUE)[[1]]
  residue_count <- length(residues)

  candidate_positions <- which(
    residues %in% c("K", "R") & seq_len(residue_count) < residue_count
  )

  if (length(candidate_positions) == 0L) {
    return(tibble::tibble(
      position = integer(0),
      residue = character(0),
      flanking_context = character(0),
      efficiency = character(0),
      rule_applied = character(0)
    ))
  }

  next_residues <- residues[candidate_positions + 1L]
  efficiency <- ifelse(
    next_residues == "P",
    "low",
    ifelse(
      next_residues %in% c("D", "E"),
      "low",
      ifelse(next_residues %in% c("K", "R"), "medium", "high")
    )
  )
  rule_applied <- ifelse(
    next_residues == "P",
    "proline_block",
    ifelse(
      next_residues %in% c("D", "E"),
      "acidic_p1_prime",
      ifelse(next_residues %in% c("K", "R"), "adjacent_basic_residues", "default_trypsin_site")
    )
  )
  flanking_context <- vapply(
    candidate_positions,
    function(position) {
      window_start <- max(1L, position - 2L)
      window_end <- min(residue_count, position + 2L)
      substr(sequence, window_start, window_end)
    },
    character(1)
  )

  tibble::tibble(
    position = as.integer(candidate_positions),
    residue = residues[candidate_positions],
    flanking_context = flanking_context,
    efficiency = efficiency,
    rule_applied = rule_applied
  )
}

.least_efficient_class <- function(efficiency) {
  if (length(efficiency) == 0L || all(is.na(efficiency))) {
    return(NA_character_)
  }

  severity <- c(low = 1L, medium = 2L, high = 3L)
  efficiency[[which.min(severity[efficiency])]]
}

.map_peptide_cleavage_efficiency <- function(start, end, site_annotations) {
  relevant_efficiency <- site_annotations$efficiency[
    site_annotations$position >= max(1L, start - 1L) &
      site_annotations$position <= end
  ]

  .least_efficient_class(relevant_efficiency)
}

.cleavage_efficiency_summary <- function(sequence, enzyme) {
  if (!.supports_cleavage_efficiency_annotations(enzyme)) {
    return(list(
      n_high_efficiency_sites = NA_integer_,
      n_low_efficiency_sites = NA_integer_
    ))
  }

  site_annotations <- .annotate_trypsin_cleavage_sites(sequence)

  list(
    n_high_efficiency_sites = sum(site_annotations$efficiency == "high"),
    n_low_efficiency_sites = sum(site_annotations$efficiency == "low")
  )
}

#' Annotate Cleavage-Site Efficiency
#'
#' `annotate_cleavage_sites()` classifies sequence-local cleavage-site
#' efficiency for trypsin-family digests. It is intended as a companion to
#' [digest_protein()] when you want to inspect cleavage hotspots before or
#' alongside peptide generation.
#'
#' @param sequence Protein input supplied as a single sequence, a named
#'   character vector of length 1, a `Biostrings::AAString`, or a FASTA file
#'   path resolving to exactly one protein.
#' @param enzyme Cleavage rule name. Cleavage-efficiency annotations are
#'   currently implemented for the trypsin family: `trypsin`, `trypsin-high`,
#'   `trypsin-low`, and `trypsin-simple`.
#'
#' @return A tibble with one row per candidate cleavage site and the columns
#'   `position`, `residue`, `flanking_context`, `efficiency`, and
#'   `rule_applied`.
#'
#' @details The current annotations are sequence-local and based on P1-P1'
#'   context only. They do not model higher-order structural accessibility,
#'   extended subsite preferences beyond P1', or PTMs that block cleavage.
#'   Unsupported enzymes currently raise an error rather than returning a
#'   partially annotated table.
#'
#' @examples
#' annotate_cleavage_sites("AKRTPK", enzyme = "trypsin")
#' @export
annotate_cleavage_sites <- function(sequence, enzyme = "trypsin") {
  normalized_input <- .read_input(sequence)

  if (length(normalized_input) != 1L) {
    cli::cli_abort(
      paste(
        "{.arg sequence} must resolve to exactly one protein for cleavage-site",
        "annotation."
      ),
      class = "pepvet_error_invalid_input"
    )
  }

  normalized_enzyme <- .normalize_enzyme(enzyme)

  if (!.supports_cleavage_efficiency_annotations(normalized_enzyme)) {
    cli::cli_abort(
      c(
        paste(
          "Cleavage-efficiency annotations are currently implemented only for",
          "the trypsin family."
        ),
        "i" = paste(
          "Supported annotation enzymes:",
          paste(.cleavage_annotation_trypsin_enzymes, collapse = ", ")
        )
      ),
      class = "pepvet_error_unsupported_cleavage_annotation"
    )
  }

  .annotate_trypsin_cleavage_sites(as.character(normalized_input[[1L]]))
}

#' Simulate a Proteolytic Digest
#'
#' `digest_protein()` normalizes protein-like input, validates the amino-acid
#' alphabet, applies a cleaver-compatible enzyme rule, and returns peptide
#' coordinates in a tidy tibble. It is the entry point for all higher-level
#' pepVet workflows, including scoring, enzyme comparison, and batch
#' evaluation.
#'
#' @param sequence Protein input supplied as a single sequence, a named
#'   character vector of sequences, a `Biostrings::AAString`, a
#'   `Biostrings::AAStringSet`, or a path to a FASTA file. File extension is
#'   not used to detect FASTA input; any existing file that
#'   `Biostrings::readAAStringSet()` can parse as FASTA is accepted.
#' @param enzyme Cleavage rule name. pepVet validates this against its
#'   hard-coded registry of cleaver-compatible enzyme names, including
#'   `trypsin`, `lysc`, `glutamyl endopeptidase`, `asp-n endopeptidase`,
#'   `chymotrypsin-high`, and `thermolysin`.
#' @param missed_cleavages Maximum number of missed cleavages to include.
#' @param include_cleavage_efficiency Logical flag indicating whether to append
#'   a `cleavage_efficiency` column to the peptide output. Trypsin-family
#'   digests receive sequence-local high/medium/low annotations; unsupported
#'   enzymes currently return `NA` in this optional column.
#' @import cleaver
#'
#' @return A tibble with one row per peptide and the columns `protein_id`,
#'   `peptide`, `start`, `end`, `length`, and `missed_cleavages`. Each row
#'   represents one observed cleavage product for one protein under the
#'   selected enzyme rule and missed-cleavage allowance. When
#'   `include_cleavage_efficiency = TRUE`, the output also includes a
#'   `cleavage_efficiency` column.
#'
#' @details FASTA record names are preserved as `protein_id` values when they
#'   are present, including irregular headers that do not use UniProt pipe
#'   formatting. Unnamed input sequences receive generated `sequence_<n>` IDs.
#'   pepVet uses cleaver-compatible cleavage rules for the strict cut sites and
#'   expands missed cleavages itself so repeated peptides and overlapping ranges
#'   retain exact start and end coordinates. Peptides are returned whether or
#'   not they later count as valid in the pepVet scoring model. Validity is a
#'   separate decision controlled by `score_peptides()` through `length_range`.
#'
#'   When cleavage-efficiency annotations are requested, pepVet records the
#'   weakest annotated efficiency class across cleavage sites that delimit or
#'   fall within a peptide. These annotations reflect local P1-P1' sequence
#'   context only. They do not model extended subsite preferences, structural
#'   protection, or PTMs that alter cleavage behavior.
#' @examples
#' digest_protein("MKWVTFISLLFLFSSAYSR")
#' digest_protein(Biostrings::AAString("MKWVTFISLLFLFSSAYSR"))
#' digest_protein("AKRTPK", include_cleavage_efficiency = TRUE)
#' @export
# nolint start: object_usage_linter.
digest_protein <- function(sequence,
                           enzyme = "trypsin",
                           missed_cleavages = 0L,
                           include_cleavage_efficiency = FALSE) {
  normalized_input <- .read_input(sequence)
  normalized_enzyme <- .normalize_enzyme(enzyme)
  max_missed_cleavages <- .validate_missed_cleavages(missed_cleavages)
  include_efficiency <- .validate_include_cleavage_efficiency(include_cleavage_efficiency)
  sequence_strings <- as.character(normalized_input)
  protein_ids <- names(normalized_input)

  digest_tables <- lapply(
    seq_along(sequence_strings),
    function(index) {
      protein_sequence <- sequence_strings[[index]]
      protein_id <- protein_ids[[index]]
      strict_ranges <- .cleavage_ranges(
        Biostrings::AAString(protein_sequence),
        enzym = normalized_enzyme
      )
      digest_ranges <- .build_digest_ranges(
        strict_ranges,
        max_missed_cleavages
      )

      digest_table <- tibble::tibble(
        protein_id = rep(protein_id, length(digest_ranges)),
        peptide = vapply(
          digest_ranges,
          function(range_row) {
            substr(protein_sequence, range_row$start, range_row$end)
          },
          character(1)
        ),
        start = as.integer(vapply(digest_ranges, `[[`, integer(1), "start")),
        end = as.integer(vapply(digest_ranges, `[[`, integer(1), "end")),
        length = as.integer(vapply(
          digest_ranges,
          function(range_row) {
            range_row$end - range_row$start + 1L
          },
          integer(1)
        )),
        missed_cleavages = as.integer(vapply(
          digest_ranges,
          `[[`,
          integer(1),
          "missed_cleavages"
        ))
      )

      if (!isTRUE(include_efficiency)) {
        return(digest_table)
      }

      if (!.supports_cleavage_efficiency_annotations(normalized_enzyme)) {
        digest_table$cleavage_efficiency <- rep(NA_character_, nrow(digest_table))
        return(digest_table)
      }

      site_annotations <- .annotate_trypsin_cleavage_sites(protein_sequence)
      digest_table$cleavage_efficiency <- vapply(
        seq_len(nrow(digest_table)),
        function(row_index) {
          .map_peptide_cleavage_efficiency(
            digest_table$start[[row_index]],
            digest_table$end[[row_index]],
            site_annotations
          )
        },
        character(1)
      )

      digest_table
    }
  )

  do.call(rbind, digest_tables)
}
# nolint end
