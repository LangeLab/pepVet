.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    sprintf("%s %s", pkgname, utils::packageVersion(pkgname))
  )
}

utils::globalVariables(c(
  "badge", "category", "charge_state", "comp_id", "component",
  "composite_score", "count", "display_id", "enzyme", "flag",
  "gravy", "is_best", "label", "length_class", "mz", "n",
  "n_wins", "pct", "pct_wins", "protein_label", "protein_length",
  "score", "tier", "value", "verdict", "x_idx", "x_val",
  "xmax", "xmin", "y_val"
))
