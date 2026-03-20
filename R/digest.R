#' Simulate a Proteolytic Digest
#'
#' `digest_protein()` is the first public entry point for pepVet. The current
#' implementation is a scaffold stub that reserves the package API while the
#' digestion pipeline is built out.
#'
#' @param sequences Protein sequence input. This will eventually accept
#'   character vectors and `Biostrings::AAStringSet` objects.
#' @param enzyme Enzyme name or cleavage rule identifier.
#' @param min_length Minimum peptide length to retain.
#' @param max_length Maximum peptide length to retain.
#'
#' @return A tibble-like object describing candidate peptides.
#' @export
digest_protein <- function(sequences,
                           enzyme = "trypsin",
                           min_length = 7L,
                           max_length = 35L) {
  pepvet_touch_runtime_imports()
  pepvet_not_implemented("digest_protein()")
}
