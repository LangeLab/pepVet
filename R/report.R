.report_component_columns <- c(
  "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
)

.verdict_bullet <- function(verdict) {
  if (length(verdict) != 1L || is.na(verdict)) {
    return(cli::symbol$info)
  }

  switch(as.character(verdict),
    Good     = cli::col_green(cli::symbol$tick),
    Moderate = cli::col_yellow(cli::symbol$warning),
    Poor     = cli::col_red(cli::symbol$cross),
    cli::symbol$info
  )
}

.verdict_colour <- function(verdict, text) {
  if (length(verdict) != 1L || is.na(verdict)) {
    return(text)
  }

  switch(as.character(verdict),
    Good     = cli::col_green(text),
    Moderate = cli::col_yellow(text),
    Poor     = cli::col_red(text),
    text
  )
}

.format_score <- function(x) {
  result <- formatC(x, format = "f", digits = 3)
  result[is.na(x)] <- "NA"
  result
}

.format_component_bar <- function(value, width = 10L) {
  if (!is.numeric(value) || length(value) != 1L || is.nan(value)) {
    .abort(
      "Component bar values must be a single numeric value.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (
    !is.numeric(width) ||
      length(width) != 1L ||
      is.na(width) ||
      !is.finite(width) ||
      width < 0 ||
      width != floor(width)
  ) {
    .abort(
      "Component bar width must be a non-negative integer.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  width <- as.integer(width)

  if (is.na(value)) {
    return(paste0(
      cli::col_silver(strrep("\u2591", width)),
      " NA"
    ))
  }

  filled <- if (is.infinite(value)) {
    if (value > 0) width else 0L
  } else {
    max(0L, min(width, round(value * width)))
  }

  paste0(
    cli::col_blue(strrep("\u2588", filled)),
    cli::col_silver(strrep("\u2591", width - filled)),
    sprintf(" %s", .format_score(value))
  )
}

.validate_report_title <- function(title) {
  if (
    !is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))
  ) {
    .abort(
      "{.arg title} must be a single character string or {.val NULL}.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  invisible(title)
}

.validate_report_score_table <- function(scores) {
  if (!inherits(scores, "data.frame") || nrow(scores) == 0L) {
    .abort(
      "A report score table must be a non-empty data frame.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (anyDuplicated(names(scores)) > 0L) {
    .abort(
      "The report score table must have unique column names.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  required <- c(
    "protein_id", .report_component_columns,
    "composite_score", "verdict"
  )
  missing_columns <- setdiff(required, names(scores))

  if (length(missing_columns) > 0L) {
    .abort(
      c(
        "The report score table is missing required columns.",
        "i" = "Missing: {.val {missing_columns}}"
      ),
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (
    !is.character(scores$protein_id) ||
      anyNA(scores$protein_id) ||
      any(!nzchar(trimws(scores$protein_id)))
  ) {
    .abort(
      "The report score table must contain non-empty character protein IDs.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  numeric_columns <- c(
    .report_component_columns,
    "composite_score",
    if ("S_unique" %in% names(scores)) "S_unique"
  )

  if (!all(vapply(scores[numeric_columns], is.numeric, logical(1)))) {
    .abort(
      "The report score table must contain numeric score columns.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (any(vapply(
    scores[numeric_columns],
    function(values) anyNA(values) || any(!is.finite(values)),
    logical(1)
  ))) {
    .abort(
      "The report score table must contain finite score values.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (any(vapply(
    scores[numeric_columns],
    function(values) any(values < 0) || any(values > 1),
    logical(1)
  ))) {
    .abort(
      "The report score table must contain scores between 0 and 1.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (
    !is.character(scores$verdict) ||
      anyNA(scores$verdict) ||
      any(!scores$verdict %in% c("Good", "Moderate", "Poor"))
  ) {
    .abort(
      "The report score table must contain valid verdict labels.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  invisible(scores)
}

.validate_report_evaluate <- function(result) {
  if (
    !is.data.frame(result$scores) ||
      !is.data.frame(result$peptides) ||
      !is.list(result$params)
  ) {
    .abort(
      "The evaluation result does not contain valid scores, peptides, and parameters.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  .validate_report_score_table(result$scores)

  required_params <- c(
    "protein_ids", "enzyme", "missed_cleavages"
  )
  missing_params <- setdiff(required_params, names(result$params))

  if (length(missing_params) > 0L) {
    .abort(
      c(
        "The evaluation result is missing report parameters.",
        "i" = "Missing: {.val {missing_params}}"
      ),
      class = "pepvet_error_invalid_report_input"
    )
  }

  protein_ids <- result$params$protein_ids
  missed_cleavages <- result$params$missed_cleavages

  if (
    !is.character(protein_ids) ||
      length(protein_ids) == 0L ||
      anyNA(protein_ids) ||
      any(!nzchar(trimws(protein_ids))) ||
      !is.character(result$params$enzyme) ||
      length(result$params$enzyme) != 1L ||
      is.na(result$params$enzyme) ||
      !nzchar(trimws(result$params$enzyme)) ||
      !is.numeric(missed_cleavages) ||
      length(missed_cleavages) != 1L ||
      is.na(missed_cleavages) ||
      !is.finite(missed_cleavages) ||
      missed_cleavages < 0 ||
      missed_cleavages != floor(missed_cleavages)
  ) {
    .abort(
      "The evaluation result contains invalid report parameters.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  if (any(!result$scores$protein_id %in% protein_ids)) {
    .abort(
      "The evaluation scores contain an unknown protein identifier.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  invisible(result)
}

.validate_report_comparison <- function(comparison) {
  required <- c(
    "enzyme", "protein_id", .report_component_columns,
    "composite_score", "verdict"
  )
  missing_columns <- setdiff(required, names(comparison))

  if (length(missing_columns) > 0L) {
    .abort(
      c(
        "The comparison result is missing report columns.",
        "i" = "Missing: {.val {missing_columns}}"
      ),
      class = "pepvet_error_invalid_report_input"
    )
  }

  .validate_report_score_table(comparison)

  if (
    !is.character(comparison$enzyme) ||
      anyNA(comparison$enzyme) ||
      any(!nzchar(trimws(comparison$enzyme))) ||
      length(unique(comparison$protein_id)) != 1L
  ) {
    .abort(
      "The comparison result must contain one protein and non-empty enzyme names.",
      class = "pepvet_error_invalid_report_input"
    )
  }

  invisible(comparison)
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
  best_index <- which.max(comparison$composite_score)
  best_enzyme <- comparison$enzyme[[best_index]]

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
    mark <- if (i == best_index) cli::col_green(">") else " "
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
#'   A non-`NULL` value must be a single character string.
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
  .validate_report_title(title)

  is_evaluate_result <- is.list(x) && !is.data.frame(x) &&
    any(c("scores", "peptides", "params") %in% names(x))
  is_comparison_result <- is.data.frame(x) &&
    any(c("enzyme", "composite_score") %in% names(x))

  if (is_evaluate_result) {
    .validate_report_evaluate(x)
    .report_evaluate(x, title)
  } else if (is_comparison_result) {
    .validate_report_comparison(x)
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
