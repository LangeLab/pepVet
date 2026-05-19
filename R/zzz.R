.onLoad <- function(libname, pkgname) {
  # Replace .pepvet_pal with an active binding that reads from the
  # mutable config environment, so pepvet_plot_config() changes take effect.
  ns <- topenv()
  if (exists(".pepvet_config_env", envir = ns, inherits = FALSE)) {
    env <- get(".pepvet_config_env", envir = ns)
    # Copy the static palette/params into the mutable environment
    if (exists(".pepvet_pal", envir = ns, inherits = FALSE)) {
      env$pal <- get(".pepvet_pal", envir = ns)
      env$pal_default <- env$pal
    }
    if (exists(".pepvet_params", envir = ns, inherits = FALSE)) {
      env$params <- get(".pepvet_params", envir = ns)
      env$params_default <- env$params
    }
    if (is.null(env$theme_overrides)) env$theme_overrides <- list()

    # Unlock bindings and replace with active bindings
    if (exists(".pepvet_pal", envir = ns, inherits = FALSE)) {
      unlockBinding(".pepvet_pal", ns)
      rm(".pepvet_pal", envir = ns)
      makeActiveBinding(".pepvet_pal", function() env$pal, ns)
      lockBinding(".pepvet_pal", ns)
    }
    if (exists(".pepvet_params", envir = ns, inherits = FALSE)) {
      unlockBinding(".pepvet_params", ns)
      rm(".pepvet_params", envir = ns)
      makeActiveBinding(".pepvet_params", function() env$params, ns)
      lockBinding(".pepvet_params", ns)
    }
  }
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
