pepvet_not_implemented <- function(feature) {
  cli::cli_abort(
    c(
      "{.pkg pepVet} does not implement {.val {feature}} yet.",
      "i" = "This repository currently contains the package skeleton only."
    ),
    class = "pepvet_not_implemented"
  )
}
