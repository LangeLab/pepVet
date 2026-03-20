.get_aa_properties <- local({
  aa_properties_cache <- NULL

  function() {
    if (is.null(aa_properties_cache)) {
      data_env <- new.env(parent = emptyenv())
      utils::data("aa_properties", package = "pepVet", envir = data_env)
      aa_properties_cache <<- data_env$aa_properties
    }

    aa_properties_cache
  }
})

.calculate_gravy <- function(peptide_sequence) {
  if (!is.character(peptide_sequence) || length(peptide_sequence) != 1L) {
    cli::cli_abort("{.arg peptide_sequence} must be a single character string.")
  }

  if (is.na(peptide_sequence)) {
    cli::cli_abort("{.arg peptide_sequence} must not be missing.")
  }

  if (!nzchar(peptide_sequence)) {
    cli::cli_abort("{.arg peptide_sequence} must not be empty.")
  }

  peptide_sequence <- toupper(peptide_sequence)
  aa_properties <- .get_aa_properties()

  residues <- strsplit(peptide_sequence, split = "", fixed = TRUE)[[1]]
  residue_index <- match(residues, aa_properties$amino_acid)

  if (anyNA(residue_index)) {
    cli::cli_abort(
      paste0(
        "Unknown amino acid code(s): ",
        paste(unique(residues[is.na(residue_index)]), collapse = ", "),
        "."
      )
    )
  }

  mean(aa_properties$hydrophobicity[residue_index])
}
