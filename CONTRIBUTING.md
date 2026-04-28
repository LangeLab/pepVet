# Contributing to pepVet

Thanks for considering a contribution to pepVet.

pepVet is an R package for proteolytic digest evaluation in proteomics workflows. The project values small, reviewable changes, explicit tests, and documentation that states scope and limitations clearly.

## Before you open an issue or pull request

- Search existing issues and pull requests first.
- Use the issue templates when they fit.
- For user questions or setup help, start with [SUPPORT.md](SUPPORT.md).
- For security-sensitive reports, follow [SECURITY.md](SECURITY.md) instead of opening a public issue.

## Development setup

pepVet depends on Bioconductor packages and uses standard R package tooling.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("Biostrings", "IRanges", "S4Vectors", "cleaver"))

install.packages(c("devtools", "roxygen2", "testthat", "pkgdown", "lintr"))
```

Then clone the repository and install the package locally.

```r
devtools::load_all()
```

## What to include in a pull request

- A focused change with a clear reason for the change.
- Tests for behavior changes when the code path is testable.
- Documentation updates when the public API, scoring interpretation, or workflow guidance changes.
- Regenerated roxygen output when function signatures or docs change.

## Recommended local checks

Run the narrowest relevant checks first, then the broader package checks before opening a PR.

```r
devtools::document()
devtools::test()
devtools::check(document = FALSE, manual = FALSE)
pkgdown::build_site()
```

`R CMD check` may emit a warning if the external `qpdf` binary is missing. That is an environment issue, not a package-code failure.

## Style expectations

- Keep changes minimal and scoped to the problem.
- Preserve existing naming and table schemas unless the change requires an API update.
- Document scientific limitations directly when adding heuristics or annotations.
- Avoid speculative claims in README, roxygen, or vignettes.

## Review notes

- Explain user-visible behavior changes in the pull request description.
- Call out any fixture updates, snapshot updates, or changes to pinned values.
- If a workflow or scoring assumption changes, state whether it is literature-backed, heuristic, or an expert prior.

## Community standards

By participating in this project, you agree to follow [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
