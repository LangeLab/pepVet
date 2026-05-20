#' Export a Peptide List for Downstream Tools
#'
#' `export_peptide_list()` filters a pepVet peptide tibble to valid peptides
#' and returns or writes the result in a format compatible with downstream
#' proteomics tools. Supported formats are `"skyline"`, `"generic"`, and
#' `"fasta"`.
#'
#' @param peptides A peptide tibble produced by [digest_protein()] or
#'   accessible via [evaluate_digest()]`$peptides`. Must contain at minimum
#'   the columns `protein_id`, `peptide`, and `length`.
#' @param format Export format. One of `"skyline"`, `"generic"`, or `"fasta"`.
#'   \describe{
#'     \item{`"skyline"`}{A tibble (or CSV when `file` is specified) with
#'       columns `Protein`, `Peptide Sequence`, `Precursor Charge`, and
#'       `Precursor Mz`. One row per peptide per charge state. M/z values are
#'       computed via [calculate_peptide_mass()].}
#'     \item{`"generic"`}{A tibble (or CSV when `file` is specified) with all
#'       pepVet peptide columns plus a computed `gravy` column and a `valid`
#'       logical column marking peptides that pass `length_range`.}
#'     \item{`"fasta"`}{A character vector (or file when `file` is specified)
#'       of FASTA-formatted records for valid peptides only. Each record has a
#'       header of the form `>protein_id|peptide_start-end`.}
#'   }
#' @param charges Integer vector of precursor charge states for the
#'   `"skyline"` format. Defaults to `2:3`. Ignored for other formats.
#' @param length_range Integer vector of length 2 defining the inclusive valid
#'   peptide length window. Defaults to `c(7L, 25L)`.
#' @param file Optional file path. When `NULL` (default), returns the result.
#'   When specified, writes to that path and returns `file` invisibly.
#'   `"skyline"` and `"generic"` are written as CSV via [utils::write.csv()];
#'   `"fasta"` is written with [base::writeLines()].
#'
#' @return When `file = NULL`: a tibble for `"skyline"` and `"generic"`, or a
#'   character vector for `"fasta"`. When `file` is specified: `file`,
#'   invisibly.
#'
#' @details Precursor m/z for Skyline export is computed as
#'   \eqn{(M + z \times 1.007276) / z} where \eqn{M} is the neutral
#'   monoisotopic peptide mass and \eqn{z} is the charge state. Skyline
#'   accepts this format via File > Import > Transition List.
#'
#' @seealso [digest_protein()], [evaluate_digest()], [calculate_peptide_mass()]
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' peps <- digest_protein(bsa_path, enzyme = "trypsin")
#' export_peptide_list(peps, format = "skyline")
#' export_peptide_list(peps, format = "generic")
#' export_peptide_list(peps, format = "fasta")
#' @export
# nolint start: object_usage_linter.
export_peptide_list <- function(peptides,
                                format = "skyline",
                                charges = 2:3,
                                length_range = c(7L, 25L),
                                file = NULL) {
  .validate_export_peptides(peptides)
  normalized_format <- .validate_export_format(format)
  normalized_length_range <- .validate_length_range(length_range)

  if (identical(normalized_format, "skyline")) {
    .validate_export_charges(charges)
  }

  if (
    !is.null(file) &&
      (!is.character(file) || length(file) != 1L || is.na(file))
  ) {
    .abort(
      "{.arg file} must be a single file path string or {.val NULL}.",
      class = "pepvet_error_invalid_file"
    )
  }

  valid_mask <- peptides$length >= normalized_length_range[[1]] &
    peptides$length <= normalized_length_range[[2]]
  valid_peps <- peptides[valid_mask, , drop = FALSE]

  result <- switch(normalized_format,
    skyline = .export_skyline(valid_peps, as.integer(charges)),
    generic = .export_generic(peptides, valid_mask),
    fasta   = .export_fasta(valid_peps)
  )

  if (is.null(file)) {
    return(result)
  }

  if (identical(normalized_format, "fasta")) {
    writeLines(result, file)
  } else {
    utils::write.csv(result, file, row.names = FALSE)
  }

  invisible(file)
}
# nolint end

# ---- Private export helpers ----

.validate_export_peptides <- function(peptides) {
  if (!inherits(peptides, "data.frame")) {
    .abort(
      "{.arg peptides} must be a peptide tibble from {.fn digest_protein} or {.fn evaluate_digest}{.code $peptides}.",
      class = "pepvet_error_invalid_export_input"
    )
  }

  required <- c("protein_id", "peptide", "length")
  missing_cols <- setdiff(required, names(peptides))

  if (length(missing_cols) > 0L) {
    .abort(
      c(
        "{.arg peptides} is missing required columns.",
        "i" = "Missing: {.val {missing_cols}}"
      ),
      class = "pepvet_error_invalid_export_input"
    )
  }

  if (nrow(peptides) == 0L) {
    .abort(
      "{.arg peptides} must contain at least one peptide row.",
      class = "pepvet_error_invalid_export_input"
    )
  }
}

.validate_export_format <- function(format) {
  if (!is.character(format) || length(format) != 1L || is.na(format)) {
    .abort(
      "{.arg format} must be a single string.",
      class = "pepvet_error_invalid_export_format"
    )
  }

  normalized <- tolower(trimws(format))
  supported <- c("skyline", "generic", "fasta")

  if (!normalized %in% supported) {
    .abort(
      c(
        "{.arg format} {.val {format}} is not supported.",
        "i" = "Supported formats: {.val {supported}}"
      ),
      class = "pepvet_error_invalid_export_format"
    )
  }

  normalized
}

.validate_export_charges <- function(charges) {
  if (
    !is.numeric(charges) ||
      length(charges) == 0L ||
      anyNA(charges) ||
      any(charges < 1L) ||
      any(charges != as.integer(charges))
  ) {
    .abort(
      "{.arg charges} must be a non-empty integer vector of positive charge states (e.g., {.code 2:3}).",
      class = "pepvet_error_invalid_charges"
    )
  }
}

.export_skyline <- function(valid_peps, charges) {
  if (nrow(valid_peps) == 0L) {
    return(tibble::tibble(
      `Protein`          = character(0),
      `Peptide Sequence` = character(0),
      `Precursor Charge` = integer(0),
      `Precursor Mz`     = numeric(0)
    ))
  }

  rows <- lapply(seq_len(nrow(valid_peps)), function(i) {
    pep <- valid_peps$peptide[[i]]
    pid <- valid_peps$protein_id[[i]]

    lapply(charges, function(z) {
      tibble::tibble(
        `Protein`          = pid,
        `Peptide Sequence` = pep,
        `Precursor Charge` = as.integer(z),
        `Precursor Mz`     = as.numeric(calculate_peptide_mass(pep, charge = z))
      )
    })
  })

  .bind_rows(do.call(c, rows))
}

.export_generic <- function(peptides, valid_mask) {
  result <- peptides
  result$gravy <- .calculate_gravy_vec(peptides$peptide)
  result$pI <- as.numeric(calculate_pI(peptides$peptide))
  result$valid <- valid_mask
  result
}

.export_fasta <- function(valid_peps) {
  if (nrow(valid_peps) == 0L) {
    return(character(0))
  }

  records <- lapply(seq_len(nrow(valid_peps)), function(i) {
    c(
      paste0(
        ">", valid_peps$protein_id[[i]],
        "|peptide_", valid_peps$start[[i]], "-", valid_peps$end[[i]]
      ),
      valid_peps$peptide[[i]]
    )
  })

  unlist(records)
}
