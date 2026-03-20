.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    sprintf("%s %s", pkgname, utils::packageVersion(pkgname))
  )
}
