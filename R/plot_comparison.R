# ── pepVet Comparison plots ────────────────────────────────────────────────
# ── plot_enzyme_comparison ────────────────────────────────────────────────────

#' Validate comparison tibble for plotting (internal helper)
#'
#' @noRd
.validate_comparison_for_plot <- function(comparison) {
  if (!is.data.frame(comparison)) {
    .abort(
      c(
        "{.arg comparison} must be a data.frame / tibble returned by {.fn compare_digests}.",
        "x" = "Got {.cls {class(comparison)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  required <- c("enzyme", "composite_score")
  missing <- setdiff(required, names(comparison))
  if (length(missing) > 0L) {
    .abort(
      c(
        "{.arg comparison} is missing required columns.",
        "x" = "Missing: {.field {missing}}.",
        "i" = "Pass the tibble returned by {.fn compare_digests} directly."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  if (nrow(comparison) < 2L) {
    .abort(
      c(
        "{.arg comparison} must contain at least 2 enzymes to compare.",
        "x" = "Only {nrow(comparison)} row{?s} found."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  invisible(comparison)
}


#' Enzyme Comparison Chart
#'
#' `plot_enzyme_comparison()` visualises the output of [compare_digests()] as
#' a dual-panel chart: component score bars (Panel A) and a composite score
#' lollipop with verdict badge (Panel B).
#'
#' Panel A shows a horizontal grouped bar chart where each enzyme occupies one
#' row and each component score is a separate colored bar, dodged side-by-side.
#' Reference lines at 0.50 and 0.75 divide the axis into poor / moderate / good
#' regions. Enzymes are sorted by composite score with the best at the top.
#'
#' Panel B shows the composite score as a lollipop, color-coded by verdict tier
#' (green >= 0.65, amber 0.40-0.64, red < 0.40). When `recommend = TRUE` a
#' gold "* Recommended" badge is appended next to the top-ranked enzyme.
#'
#' @param comparison A tibble returned by [compare_digests()].  Must contain
#'   at least the columns `enzyme` and `composite_score`, plus whichever
#'   component-score columns are requested in `scores`.
#' @param scores Character vector of component-score column names to display
#'   in Panel A.  Any column absent from `comparison` is silently dropped.
#'   Defaults to all five standard component scores.
#' @param recommend Logical.  When `TRUE` (default) a "Recommended" badge
#'   marks the enzyme with the highest composite score in Panel B.
#' @param title Optional character string for the overall plot title.
#'   Auto-generated from the protein accession when `NULL` (default).
#'
#' @return A `patchwork` object (two ggplot panels side by side) that can be
#'   printed, further customised, or saved with [ggplot2::ggsave()].
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'   requireNamespace("patchwork", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   comp <- compare_digests(bsa_path,
#'     enzymes = c(
#'       "trypsin", "lysc",
#'       "glutamyl endopeptidase"
#'     )
#'   )
#'   p <- plot_enzyme_comparison(comp)
#'   print(p)
#' }
#'
#' @seealso [compare_digests()], [recommend_enzyme()],
#'   [plot_digest_profile()], [plot_coverage_map()], [plot_batch_comparison()]
#' @family plot-comparison
#' @export
plot_enzyme_comparison <- function(
  comparison,
  scores = c("S_coverage", "S_length", "S_count", "S_hydro", "S_charge"),
  recommend = TRUE,
  title = NULL
) {
  rlang::check_installed("ggplot2", reason = "to use plot_enzyme_comparison()")
  rlang::check_installed("patchwork", reason = "to use plot_enzyme_comparison()")

  .validate_comparison_for_plot(comparison)

  # ── Restrict to requested scores that actually exist ──────────────────────
  scores <- intersect(scores, names(comparison))
  if (length(scores) == 0L) {
    .abort(
      c(
        "None of the requested {.arg scores} columns were found in
         {.arg comparison}.",
        "i" = "Available columns: {.field {names(comparison)}}."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }

  # ── Display names for score components ────────────────────────────────────
  score_labels <- c(
    S_coverage = "Coverage",
    S_length   = "Length",
    S_count    = "Peptide count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge"
  )
  display_names <- ifelse(
    scores %in% names(score_labels),
    score_labels[scores],
    scores
  )

  # ── Distinct JCO-inspired colors, one per component ───────────────────────
  component_colors <- .pepvet_pal$component
  col_map <- component_colors[scores]
  names(col_map) <- display_names

  # ── Enzyme factor: sorted worst → best composite (top of chart = best) ────
  ordered_enzymes <- comparison$enzyme[order(comparison$composite_score)]
  comparison$enzyme <- factor(comparison$enzyme, levels = ordered_enzymes)

  # ── Reshape to long for Panel A ───────────────────────────────────────────
  long <- .bind_rows(lapply(seq_along(scores), function(i) {
    data.frame(
      enzyme = comparison$enzyme,
      score_name = display_names[[i]],
      value = comparison[[scores[[i]]]],
      stringsAsFactors = FALSE
    )
  }))
  long$score_name <- factor(long$score_name, levels = rev(display_names))
  long$enzyme <- factor(long$enzyme, levels = levels(comparison$enzyme))

  # ── Verdict tier color for composite lollipop heads ───────────────────────
  tier_color <- function(x) {
    ifelse(x >= .get_param("verdict_good"), .pepvet_pal$good,
      ifelse(x >= .get_param("verdict_moderate"), .pepvet_pal$moderate, .pepvet_pal$poor)
    )
  }
  comparison$comp_color <- tier_color(comparison$composite_score)
  comparison$comp_label <- sprintf("%.2f", comparison$composite_score)

  # ── Auto title ────────────────────────────────────────────────────────────
  auto_title <- if (!is.null(title)) {
    title
  } else if ("protein_id" %in% names(comparison)) {
    pid <- comparison$protein_id[[1L]]
    paste0(.tidy_protein_id(pid), "  \u00b7  Enzyme comparison")
  } else {
    "Enzyme comparison"
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # Panel A: component scores grouped bar chart
  # ═══════════════════════════════════════════════════════════════════════════
  pa <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$value, y = .data$enzyme, fill = .data$score_name)
  ) +
    # Good-region shading
    ggplot2::annotate("rect",
      xmin = .get_param("verdict_good"), xmax = 1.0, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.55
    ) +
    # Threshold reference lines
    ggplot2::geom_vline(
      xintercept = c(.get_param("verdict_moderate"), .get_param("verdict_good")),
      color = .pepvet_pal$separator, linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.68, alpha = .get_param("scatter_alpha")
    ) +
    # Composite score: bold dark vertical tick spanning the full enzyme row
    ggplot2::geom_segment(
      data = comparison,
      ggplot2::aes(
        x    = .data$composite_score,
        xend = .data$composite_score,
        y    = as.integer(.data$enzyme) - 0.43,
        yend = as.integer(.data$enzyme) + 0.43
      ),
      color = .pepvet_pal$brand_dark, linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    # Threshold region labels
    ggplot2::annotate("text",
      x = c(
        .get_param("verdict_moderate") + 0.01,
        .get_param("verdict_good") + 0.06
      ), y = 0.52,
      label = c("Moderate", "Good"),
      hjust = 0, size = 2.4,
      color = .pepvet_pal$text_secondary, fontface = "italic"
    ) +
    ggplot2::scale_fill_manual(
      values = col_map,
      name = "Component",
      guide = ggplot2::guide_legend(
        reverse = TRUE,
        override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.0),
      labels = c("0", ".25", ".50", ".75", "1"),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      tag      = "A",
      x        = "Score  [0 \u2013 1]",
      y        = NULL,
      subtitle = "Component scores per enzyme  \u00b7  | = composite"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position    = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )

  # ═══════════════════════════════════════════════════════════════════════════
  # Panel B: composite score lollipop + verdict badge
  # ═══════════════════════════════════════════════════════════════════════════
  best_enzyme <- as.character(
    comparison$enzyme[which.max(comparison$composite_score)]
  )
  badge_df <- comparison[
    as.character(comparison$enzyme) == best_enzyme, ,
    drop = FALSE
  ]

  pb <- ggplot2::ggplot(
    comparison,
    ggplot2::aes(x = .data$composite_score, y = .data$enzyme)
  ) +
    ggplot2::annotate("rect",
      xmin = .get_param("verdict_good"), xmax = 1.1, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.55
    ) +
    ggplot2::geom_vline(
      xintercept = c(.get_param("verdict_moderate"), .get_param("verdict_good")),
      color = .pepvet_pal$separator, linetype = "dashed", linewidth = 0.4
    ) +
    # Lollipop stems
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0, xend = .data$composite_score,
        y = .data$enzyme, yend = .data$enzyme
      ),
      color = .pepvet_pal$separator, linewidth = 0.6
    ) +
    # Lollipop heads colored by verdict tier
    ggplot2::geom_point(
      ggplot2::aes(color = I(.data$comp_color)),
      size = 4.5
    ) +
    # Score value label
    ggplot2::geom_text(
      ggplot2::aes(label = .data$comp_label),
      hjust = -0.35, size = 3.0,
      color = .pepvet_pal$text_axis_title, fontface = "bold"
    ) +
    # Recommended badge on the best enzyme
    {
      if (recommend) {
        ggplot2::annotate(
          "label",
          x = badge_df$composite_score[[1L]] + 0.01,
          y = best_enzyme,
          label = "\u2605 Recommended",
          hjust = -0.05,
          vjust = -0.55,
          size = 2.6,
          color = .pepvet_pal$badge_gold_text,
          fill = .pepvet_pal$badge_gold_fill,
          linewidth = 0,
          fontface = "bold"
        )
      }
    } +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.0),
      labels = c("0", ".25", ".50", ".75", "1"),
      expand = ggplot2::expansion(mult = c(0, 0.28))
    ) +
    ggplot2::labs(
      tag      = "B",
      x        = "Composite score",
      y        = NULL,
      subtitle = "Overall ranking"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  # ═══════════════════════════════════════════════════════════════════════════
  # Assemble with patchwork
  # ═══════════════════════════════════════════════════════════════════════════
  (pa | pb) +
    patchwork::plot_layout(widths = c(2.8, 1)) +
    patchwork::plot_annotation(
      title = auto_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 15, face = "bold",
          color = .pepvet_pal$brand_dark,
          margin = ggplot2::margin(b = 6)
        )
      )
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        size = 14, face = "bold", color = .pepvet_pal$brand
      )
    )
}
