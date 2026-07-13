package_source_root <- function() {
  candidates <- unique(c(
    testthat::test_path("../.."),
    getwd(),
    find.package("pepVet", quiet = TRUE)
  ))
  candidates <- candidates[nzchar(candidates)]
  has_description <- vapply(
    candidates,
    function(path) file.exists(file.path(path, "DESCRIPTION")),
    logical(1)
  )

  if (!any(has_description)) {
    stop("Could not locate the pepVet package source root.", call. = FALSE)
  }

  normalizePath(candidates[[which(has_description)[[1L]]]], mustWork = TRUE)
}

quote_process_arg <- function(value) {
  shell_type <- if (.Platform$OS.type == "windows") "cmd" else "sh"
  shQuote(value, type = shell_type)
}

process_path <- function(path, must_work = TRUE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

is_pkgload_package <- function() {
  "pkgload" %in% loadedNamespaces() && isTRUE(
    get("is_dev_package", asNamespace("pkgload"))("pepVet")
  )
}

test_that("an installed package loads cleanly and initializes state", {
  skip_if(
    identical(Sys.getenv("R_COVR"), "true"),
    "nested package installation is checked outside coverage runs"
  )

  temp_root <- withr::local_tempdir("pepVet-package-load-")
  rscript_command <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )

  if (is_pkgload_package()) {
    source_root <- process_path(package_source_root())
    install_library <- file.path(temp_root, "library")
    dir.create(install_library)
    install_library <- process_path(install_library)
    r_command <- file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "R.exe" else "R"
    )

    install_output <- system2(
      r_command,
      args = c(
        "CMD",
        "INSTALL",
        "--no-multiarch",
        quote_process_arg(paste0("--library=", install_library)),
        quote_process_arg(source_root)
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    install_status <- attr(install_output, "status")
    if (is.null(install_status)) install_status <- 0L

    expect_identical(
      install_status,
      0L,
      info = paste(install_output, collapse = "\n")
    )

    if (install_status != 0L) return(invisible(NULL))
  } else {
    install_library <- dirname(
      process_path(find.package("pepVet", quiet = FALSE))
    )
  }

  child_script <- file.path(temp_root, "check-load.R")
  writeLines(
    c(
      "lib <- commandArgs(trailingOnly = TRUE)[[1L]]",
      "library(pepVet, lib.loc = lib)",
      "ns <- asNamespace(\"pepVet\")",
      "config <- get(\".pepvet_config_env\", envir = ns)",
      "stopifnot(bindingIsActive(\".pepvet_pal\", ns))",
      "stopifnot(bindingIsActive(\".pepvet_params\", ns))",
      "stopifnot(identical(config$pal, config$pal_default))",
      "stopifnot(identical(config$params, config$params_default))",
      "pepVet::pepvet_plot_config(",
      "  palette = list(brand = \"#004488\"),",
      "  params = list(verdict_good = 0.7)",
      ")",
      "custom_config <- pepVet::pepvet_plot_config()",
      "get(\".onLoad\", envir = ns)(NULL, \"pepVet\")",
      "reloaded_config <- pepVet::pepvet_plot_config()",
      "stopifnot(identical(custom_config$palette, reloaded_config$palette))",
      "stopifnot(identical(custom_config$params, reloaded_config$params))",
      "data_path <- system.file(\"data\", \"aa_properties.rda\", package = \"pepVet\")",
      "stopifnot(nzchar(data_path), file.exists(data_path))",
      "extdata_path <- system.file(\"extdata\", \"P02769.fasta\", package = \"pepVet\")",
      "stopifnot(nzchar(extdata_path), file.exists(extdata_path))",
      "unloadNamespace(\"pepVet\")",
      "library(pepVet, lib.loc = lib)",
      "ns <- asNamespace(\"pepVet\")",
      "config <- get(\".pepvet_config_env\", envir = ns)",
      "stopifnot(bindingIsActive(\".pepvet_pal\", ns))",
      "stopifnot(bindingIsActive(\".pepvet_params\", ns))",
      "stopifnot(identical(config$pal, config$pal_default))",
      "stopifnot(identical(config$params, config$params_default))"
    ),
    child_script
  )

  child_output <- system2(
    rscript_command,
    args = c(
      "--vanilla",
      quote_process_arg(process_path(child_script)),
      quote_process_arg(install_library)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  child_status <- attr(child_output, "status")
  if (is.null(child_status)) child_status <- 0L

  expect_identical(
    child_status,
    0L,
    info = paste(child_output, collapse = "\n")
  )
  expect_length(child_output, 0L)
})
