#' Export a peptide list for downstream tools
#'
#' `export_peptide_list()` filters a pepVet peptide tibble to valid peptides
#' and returns or writes the result in a format compatible with downstream
#' proteomics tools. Supported formats are `"skyline"`, `"generic"`, and
#' `"fasta"`.
#'
#' @param peptides A peptide tibble produced by [digest_protein()] or
#'   accessible via [evaluate_digest()]`$peptides`. Must contain at minimum
#'   the columns `protein_id`, `peptide`, and `length`. If `NULL` or `NA`,
#'   raises an error.
#' @param format Export format. Defaults to `"skyline"`.
#'   \describe{
#'     \item{\code{"skyline"}}{A tibble (or CSV when `file` is specified) with
#'       columns \code{Protein}, \code{Peptide Sequence},
#'       \code{Precursor Charge}, and \code{Precursor Mz}. One row per peptide
#'       per charge state. M/z values are computed via
#'       [calculate_peptide_mass()].}
#'     \item{\code{"generic"}}{A tibble (or CSV when `file` is specified) with
#'       all pepVet peptide columns plus a computed `gravy` column and a
#'       `valid` logical column marking peptides that pass `length_range`.}
#'     \item{\code{"fasta"}}{A character vector (or file when `file` is
#'       specified) of FASTA-formatted records for valid peptides only. Each
#'       record has a header of the form `>protein_id|peptide_start-end`.
#'       Requires `start` and `end` columns in `peptides`.}
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
#' @details Precursor m/z for Skyline export is computed as
#'   \eqn{(M + z \times 1.007276) / z} where \eqn{M} is the neutral
#'   monoisotopic peptide mass and \eqn{z} is the charge state. Skyline
#'   accepts this format via File > Import > Transition List.
#'
#'   When no peptides pass the `length_range` filter, skyline format returns
#'   an empty tibble and fasta format returns an empty character vector.
#'
#' @family export
#' @section Limitations:
#'   Only peptides passing the specified `length_range` are exported.
#'
#' @return When `file = NULL`:
#'   \itemize{
#'     \item For `"skyline"`: a tibble with columns \code{Protein},
#'       \code{Peptide Sequence}, \code{Precursor Charge}, and
#'       \code{Precursor Mz}.
#'     \item For `"generic"`: a tibble with all original columns plus
#'       \code{gravy}, \code{pI}, and \code{valid}.
#'     \item For `"fasta"`: a character vector of FASTA records.
#'   }
#'   When `file` is specified: `file`, invisibly.
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
  normalized_format <- .validate_export_format(format)
  .validate_export_peptides(peptides, normalized_format)
  normalized_length_range <- .validate_length_range(length_range)
  normalized_file <- .validate_export_file(file)

  normalized_charges <- if (identical(normalized_format, "skyline")) {
    .validate_export_charges(charges)
  } else {
    NULL
  }

  valid_mask <- peptides$length >= normalized_length_range[[1]] &
    peptides$length <= normalized_length_range[[2]]
  valid_peps <- peptides[valid_mask, , drop = FALSE]

  result <- switch(normalized_format,
    skyline = .export_skyline(valid_peps, normalized_charges),
    generic = .export_generic(peptides, valid_mask),
    fasta   = .export_fasta(valid_peps)
  )

  if (is.null(normalized_file)) {
    return(result)
  }

  write_failure <- function(condition) {
    .abort(
      c(
        "Could not write the export file.",
        "i" = conditionMessage(condition)
      ),
      class = "pepvet_error_invalid_file"
    )
  }

  tryCatch(
    if (identical(normalized_format, "fasta")) {
      writeLines(result, normalized_file)
    } else {
      utils::write.csv(result, normalized_file, row.names = FALSE)
    },
    warning = write_failure,
    error = write_failure
  )

  invisible(normalized_file)
}
# nolint end

## Private export helpers

.validate_export_peptides <- function(peptides, format = NULL) {
  if (!inherits(peptides, "data.frame")) {
    .abort(
      paste0(
        "{.arg peptides} must be a peptide tibble from ",
        "{.fn digest_protein} or {.fn evaluate_digest}{.code $peptides}."
      ),
      class = "pepvet_error_invalid_export_input"
    )
  }

  if (anyDuplicated(names(peptides)) > 0L) {
    .abort(
      "{.arg peptides} must have unique column names.",
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

  if (
    !is.character(peptides$protein_id) ||
      !is.character(peptides$peptide) ||
    !is.numeric(peptides$length)
  ) {
    .abort(
      paste0(
        "{.arg peptides} must contain character identifiers and sequences ",
        "plus a numeric length column."
      ),
      class = "pepvet_error_invalid_export_input"
    )
  }

  if (
    anyNA(peptides$protein_id) ||
      anyNA(peptides$peptide) ||
      anyNA(peptides$length) ||
      any(!is.finite(peptides$length)) ||
      any(peptides$length < 1) ||
      any(peptides$length > .Machine$integer.max) ||
      any(peptides$length != floor(peptides$length))
  ) {
    .abort(
      "{.arg peptides} contains missing, non-finite, or invalid peptide lengths.",
      class = "pepvet_error_invalid_export_input"
    )
  }

  if (
    any(!nzchar(peptides$protein_id)) ||
      any(!nzchar(peptides$peptide))
  ) {
    .abort(
      "{.arg peptides} must not contain empty protein identifiers or peptide sequences.",
      class = "pepvet_error_invalid_export_input"
    )
  }

  if (nrow(peptides) > 0L) {
    tryCatch(
      .normalize_peptide_sequences(peptides$peptide, arg_name = "peptides"),
      error = function(error) {
        .abort(
          "{.arg peptides} contains an invalid amino-acid sequence.",
          class = "pepvet_error_invalid_export_input"
        )
      }
    )

    if (any(as.integer(peptides$length) != nchar(peptides$peptide))) {
      .abort(
        "{.arg peptides} has lengths that do not match its peptide sequences.",
        class = "pepvet_error_invalid_export_input"
      )
    }
  }

  if (identical(format, "fasta")) {
    coordinate_columns <- c("start", "end")
    missing_coordinates <- setdiff(coordinate_columns, names(peptides))

    if (length(missing_coordinates) > 0L) {
      .abort(
        c(
          "{.arg peptides} is missing FASTA coordinate columns.",
          "i" = "Missing: {.val {missing_coordinates}}"
        ),
        class = "pepvet_error_invalid_export_input"
      )
    }

    if (
      !is.numeric(peptides$start) ||
        !is.numeric(peptides$end) ||
        anyNA(peptides$start) ||
        anyNA(peptides$end) ||
        any(!is.finite(peptides$start)) ||
        any(!is.finite(peptides$end)) ||
        any(peptides$start < 1) ||
        any(peptides$end < peptides$start) ||
        any(peptides$start != floor(peptides$start)) ||
        any(peptides$end != floor(peptides$end))
    ) {
      .abort(
        "{.arg peptides} contains invalid FASTA coordinates.",
        class = "pepvet_error_invalid_export_input"
      )
    }

    if (any(peptides$length != peptides$end - peptides$start + 1L)) {
      .abort(
        "{.arg peptides} has coordinates that do not match its lengths.",
        class = "pepvet_error_invalid_export_input"
      )
    }
  }

  invisible(peptides)
}

.validate_export_file <- function(file) {
  if (is.null(file)) {
    return(NULL)
  }

  if (
    !is.character(file) ||
      length(file) != 1L ||
      is.na(file) ||
      !nzchar(trimws(file)) ||
      any(as.integer(charToRaw(file)) == 0L)
  ) {
    .abort(
      "{.arg file} must be a non-empty file path string or {.val NULL}.",
      class = "pepvet_error_invalid_file"
    )
  }

  file
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
      any(!is.finite(charges)) ||
      any(charges < 1L) ||
      any(charges > .Machine$integer.max) ||
      any(charges != floor(charges))
  ) {
    .abort(
      paste0(
        "{.arg charges} must be a non-empty integer vector ",
        "of positive charge states (e.g., {.code 2:3})."
      ),
      class = "pepvet_error_invalid_charges"
    )
  }

  as.integer(charges)
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
  result$gravy <- .calculate_gravy(peptides$peptide)
  result$pI <- if (nrow(peptides) == 0L) {
    numeric(0)
  } else {
    as.numeric(calculate_pI(peptides$peptide))
  }
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
