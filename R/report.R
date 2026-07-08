.verdict_bullet <- function(verdict) {
  switch(verdict,
    Good     = cli::col_green(cli::symbol$tick),
    Moderate = cli::col_yellow(cli::symbol$warning),
    Poor     = cli::col_red(cli::symbol$cross),
    cli::symbol$info
  )
}

.verdict_colour <- function(verdict, text) {
  switch(verdict,
    Good     = cli::col_green(text),
    Moderate = cli::col_yellow(text),
    Poor     = cli::col_red(text),
    text
  )
}

.format_score <- function(x) {
  formatC(x, format = "f", digits = 3)
}

.format_component_bar <- function(value, width = 10L) {
  filled <- max(0L, min(width, round(value * width)))
  paste0(
    cli::col_blue(strrep("\u2588", filled)),
    cli::col_silver(strrep("\u2591", width - filled)),
    sprintf(" %s", .format_score(value))
  )
}

.print_single_result <- function(scores_row, protein_id = NULL) {
  pid <- if (!is.null(protein_id)) protein_id else scores_row$protein_id
  verdict <- scores_row$verdict

  cli::cat_line(
    cli::style_bold(pid),
    "  ",
    .verdict_bullet(verdict),
    " ",
    .verdict_colour(verdict, cli::style_bold(verdict)),
    "  (composite: ",
    cli::style_bold(.format_score(scores_row$composite_score)),
    ")"
  )

  cli::cat_line(
    "  S_length   ", .format_component_bar(scores_row$S_length)
  )
  cli::cat_line(
    "  S_coverage ", .format_component_bar(scores_row$S_coverage)
  )
  cli::cat_line(
    "  S_count    ", .format_component_bar(scores_row$S_count)
  )
  cli::cat_line(
    "  S_hydro    ", .format_component_bar(scores_row$S_hydro)
  )
  cli::cat_line(
    "  S_charge   ", .format_component_bar(scores_row$S_charge)
  )

  if ("S_unique" %in% names(scores_row)) {
    cli::cat_line(
      "  S_unique   ", .format_component_bar(scores_row$S_unique)
    )
  }
}

.print_comparison_table <- function(comparison, protein_id = NULL) {
  pid <- if (!is.null(protein_id)) protein_id else comparison$protein_id[[1]]
  best_enzyme <- comparison$enzyme[[1]]

  cli::cat_rule(
    left = cli::style_bold(pid),
    right = paste0("best: ", cli::style_bold(best_enzyme))
  )

  header <- sprintf(
    "  %-32s  %s  %s  %s  %s  %s  %s  %s",
    "enzyme",
    "S_len", "S_cov", "S_cnt", "S_hyd", "S_chg",
    "composite", "verdict"
  )
  cli::cat_line(cli::style_italic(header))
  cli::cat_line(cli::col_silver(strrep("-", nchar(header))))

  for (i in seq_len(nrow(comparison))) {
    row <- comparison[i, , drop = FALSE]
    mark <- if (i == 1L) cli::col_green(">") else " "
    verb <- sprintf("%-8s", row$verdict)

    cli::cat_line(sprintf(
      "%s %-32s  %s  %s  %s  %s  %s  %s  %s",
      mark,
      row$enzyme,
      .format_score(row$S_length),
      .format_score(row$S_coverage),
      .format_score(row$S_count),
      .format_score(row$S_hydro),
      .format_score(row$S_charge),
      .format_score(row$composite_score),
      .verdict_colour(row$verdict, verb)
    ))
  }
}

#' Print a styled console report for a proteolytic digest
#'
#' `digest_report()` formats the output of [evaluate_digest()] or
#' [compare_digests()] as a human-readable styled console summary. The
#' function returns its input invisibly so it can be used in pipelines. Use it
#' when you want a compact review of component scores during interactive enzyme
#' selection or package-level demonstrations.
#'
#' @param x The object to report on. Accepts:
#'   \describe{
#'     \item{Named list from [evaluate_digest()]}{Prints a single-protein
#'       component-bar summary for the evaluated enzyme.}
#'     \item{Tibble from [compare_digests()]}{Prints a multi-enzyme ranking
#'       table with the best enzyme highlighted.}
#'   }
#'   If `NULL` or an unrecognised type, raises an error.
#' @param title Optional character string printed as a section header above
#'   the report. When `NULL` (default), the protein ID is used as the header.
#'
#' @family report
#' @section Limitations:
#'   Output is printed to the console only. File output is not supported.
#'
#' @return `x`, invisibly.
#'
#' @seealso [evaluate_digest()], [compare_digests()]
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' ev <- evaluate_digest(bsa_path)
#' digest_report(ev)
#'
#' comp <- compare_digests(bsa_path, enzymes = c("trypsin", "lysc"))
#' digest_report(comp)
#' @export
# nolint start: object_usage_linter.
digest_report <- function(x, title = NULL) {
  is_evaluate_result <- is.list(x) && !is.data.frame(x) &&
    all(c("scores", "peptides", "params") %in% names(x))
  is_comparison_result <- is.data.frame(x) &&
    "enzyme" %in% names(x) && "composite_score" %in% names(x)

  if (is_evaluate_result) {
    .report_evaluate(x, title)
  } else if (is_comparison_result) {
    .report_comparison(x, title)
  } else {
    .abort(
      paste0(
        "{.arg x} must be a list from {.fn evaluate_digest} ",
        "or a tibble from {.fn compare_digests}."
      ),
      class = "pepvet_error_invalid_report_input"
    )
  }

  invisible(x)
}
# nolint end

.report_evaluate <- function(ev, title) {
  scores <- ev$scores
  params <- ev$params

  header_text <- if (!is.null(title)) {
    title
  } else {
    paste0(params$protein_ids[[1]], "  /  ", params$enzyme,
           "  (mc=", params$missed_cleavages, ")")
  }

  cli::cat_rule(left = cli::style_bold(header_text))

  for (i in seq_len(nrow(scores))) {
    row <- scores[i, , drop = FALSE]

    .print_single_result(row, protein_id = row$protein_id)

    if (i < nrow(scores)) {
      cli::cat_line()
    }
  }

  cli::cat_rule()
}

.report_comparison <- function(comparison, title) {
  protein_id <- comparison$protein_id[[1]]
  header_text <- if (!is.null(title)) title else protein_id

  .print_comparison_table(comparison, protein_id = header_text)
  cli::cat_rule()
}

#' Quick digest check for a single protein
#'
#' `pepvet_check()` is a convenience wrapper that evaluates a protein digest
#' and immediately prints a styled console report. It is intended for
#' interactive use and first-time exploration where a single call is more
#' useful than manually wiring [evaluate_digest()] and [digest_report()].
#'
#' @param sequence Protein input. Accepts the same forms as [evaluate_digest()].
#'   If `NULL`, raises an error.
#' @param enzyme Enzyme name. Defaults to `"trypsin"`.  If `NULL`, raises an
#'   error.
#' @param ... Additional arguments passed to [evaluate_digest()], such as
#'   `missed_cleavages`, `include_cleavage_efficiency`, `weights`,
#'   `gravy_range`, and `length_range`.
#'
#' @family report
#'
#' @return The [evaluate_digest()] result list, invisibly.
#'
#' @seealso [evaluate_digest()], [digest_report()]
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' pepvet_check(bsa_path, enzyme = "trypsin")
#' @export
pepvet_check <- function(sequence, enzyme = "trypsin", ...) {
  result <- evaluate_digest(sequence, enzyme = enzyme, ...)
  digest_report(result)
  invisible(result)
}
