.cleavage_ranges <- get("cleavageRanges", envir = asNamespace("cleaver"))

#' Simulate a Proteolytic Digest
#'
#' `digest_protein()` normalizes protein-like input, validates the amino-acid
#' alphabet, applies a cleaver-compatible enzyme rule, and returns peptide
#' coordinates in a tidy tibble.
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
#' @import cleaver
#'
#' @return A tibble with one row per peptide and the columns `protein_id`,
#'   `peptide`, `start`, `end`, `length`, and `missed_cleavages`.
#'
#' @details FASTA record names are preserved as `protein_id` values when they
#'   are present, including irregular headers that do not use UniProt pipe
#'   formatting. Unnamed input sequences receive generated `sequence_<n>` IDs.
#'   pepVet uses cleaver-compatible cleavage rules for the strict cut sites and
#'   expands missed cleavages itself so repeated peptides and overlapping ranges
#'   retain exact start and end coordinates.
#' @examples
#' digest_protein("MKWVTFISLLFLFSSAYSR")
#' digest_protein(Biostrings::AAString("MKWVTFISLLFLFSSAYSR"))
#' @export
# nolint start: object_usage_linter.
digest_protein <- function(sequence,
                           enzyme = "trypsin",
                           missed_cleavages = 0L) {
  normalized_input <- .read_input(sequence)
  normalized_enzyme <- .normalize_enzyme(enzyme)
  max_missed_cleavages <- .validate_missed_cleavages(missed_cleavages)
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

      tibble::tibble(
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
    }
  )

  do.call(rbind, digest_tables)
}
# nolint end
