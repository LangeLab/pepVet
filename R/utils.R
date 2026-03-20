pepvet_touch_runtime_imports <- function() {
  invisible(list(
    cleaver = cleaver::cleave,
    Biostrings = Biostrings::AAStringSet,
    IRanges = IRanges::IRanges,
    S4Vectors = S4Vectors::DataFrame,
    tibble = tibble::tibble,
    dplyr = dplyr::mutate,
    rlang = rlang::check_installed
  ))
}

pepvet_not_implemented <- function(feature) {
  cli::cli_abort(
    c(
      "{.pkg pepVet} does not implement {.val {feature}} yet.",
      "i" = "This repository currently contains the package skeleton only."
    ),
    class = "pepvet_not_implemented"
  )
}
