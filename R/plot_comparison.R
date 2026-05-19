# ── pepVet Comparison plots ────────────────────────────────────────────────
# ── plot_enzyme_comparison ────────────────────────────────────────────────────

#' Validate comparison tibble for plotting (internal helper)
#'
#' @noRd
.validate_comparison_for_plot <- function(comparison) {
  if (!is.data.frame(comparison)) {
    .abort(
      c(
        "{.arg comparison} must be a data.frame / tibble returned by
         {.fn compare_digests}.",
        "x" = "Got {.cls {class(comparison)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  required <- c("enzyme", "composite_score")
  missing  <- setdiff(required, names(comparison))
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
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   comp <- compare_digests(bsa_path,
#'                           enzymes = c("trypsin", "lysc",
#'                                       "glutamyl endopeptidase"))
#'   p <- plot_enzyme_comparison(comp)
#'   print(p)
#' }
#'
#' @seealso [compare_digests()], [recommend_enzyme()],
#'   [plot_digest_profile()], [plot_coverage_map()]
#' @export
plot_enzyme_comparison <- function(
    comparison,
    scores    = c("S_coverage", "S_length", "S_count", "S_hydro", "S_charge"),
    recommend = TRUE,
    title     = NULL
) {
  rlang::check_installed("ggplot2",   reason = "to use plot_enzyme_comparison()")
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
  component_colors <- c(
    S_coverage = "#2C5F8A",   # brand blue
    S_length   = "#27AE60",   # good green
    S_count    = "#E8A838",   # amber
    S_hydro    = "#8B5E99",   # purple
    S_charge   = "#4AAFB0"    # teal
  )
  col_map        <- component_colors[scores]
  names(col_map) <- display_names

  # ── Enzyme factor: sorted worst → best composite (top of chart = best) ────
  ordered_enzymes <- comparison$enzyme[order(comparison$composite_score)]
  comparison$enzyme <- factor(comparison$enzyme, levels = ordered_enzymes)

  # ── Reshape to long for Panel A ───────────────────────────────────────────
  long <- do.call(rbind, lapply(seq_along(scores), function(i) {
    data.frame(
      enzyme     = comparison$enzyme,
      score_name = display_names[[i]],
      value      = comparison[[scores[[i]]]],
      stringsAsFactors = FALSE
    )
  }))
  long$score_name <- factor(long$score_name, levels = rev(display_names))
  long$enzyme     <- factor(long$enzyme,     levels = levels(comparison$enzyme))

  # ── Verdict tier color for composite lollipop heads ───────────────────────
  tier_color <- function(x) {
    ifelse(x >= .get_param("verdict_good"), .pepvet_pal$good,
    ifelse(x >= .get_param("verdict_moderate"), .pepvet_pal$moderate, .pepvet_pal$poor))
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
      x = c(0.41, 0.71), y = 0.52,
      label  = c("Moderate", "Good"),
      hjust  = 0, size = 2.4,
      color  = "#999999", fontface = "italic"
    ) +
    ggplot2::scale_fill_manual(
      values = col_map,
      name   = "Component",
      guide  = ggplot2::guide_legend(
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
    as.character(comparison$enzyme) == best_enzyme, , drop = FALSE
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
      ggplot2::aes(x = 0, xend = .data$composite_score,
                   y = .data$enzyme, yend = .data$enzyme),
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
      color = "#444444", fontface = "bold"
    ) +
    # Recommended badge on the best enzyme
    {
      if (recommend) {
        ggplot2::annotate("label",
          x        = badge_df$composite_score[[1L]] + 0.01,
          y        = best_enzyme,
          label    = "\u2605 Recommended",
          hjust    = -0.05, vjust = -0.55,
          size     = 2.6,
          color    = "#7A5A00",
          fill     = "#FFF5CC",
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
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
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
          size   = 13, face = "bold",
          color  = .pepvet_pal$brand_dark,
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


# ── plot_protein_comparison ───────────────────────────────────────────────────

#' Protein Comparison: Component Scores Across Multiple Proteins
#'
#' `plot_protein_comparison()` draws a grouped bar chart of component scores
#' for multiple proteins digested with the same enzyme, allowing direct
#' side-by-side comparison.  It is the multi-protein mirror of
#' [plot_enzyme_comparison()].
#'
#' @param results A named list of [evaluate_digest()] results (each for a
#'   different protein, same enzyme), **or** a `batch_evaluate()` tibble.
#'   List names are used as protein labels; when unnamed the function falls
#'   back to protein IDs extracted from `result$params`.
#' @param components Character vector of component score column names to show.
#'   Defaults to all standard components (`S_length`, `S_coverage`, `S_count`,
#'   `S_hydro`, `S_charge`) plus `composite_score`.
#' @param show_verdict Logical.  When `TRUE` (default) a verdict badge is
#'   drawn above each protein's group.
#' @param title Optional character string for the plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
#'   results  <- list(
#'     BSA = evaluate_digest(bsa_path, enzyme = "trypsin"),
#'     H3  = evaluate_digest(h3_path,  enzyme = "trypsin")
#'   )
#'   p <- plot_protein_comparison(results)
#'   print(p)
#' }
#'
#' @seealso [plot_enzyme_comparison()], [evaluate_digest()], [batch_evaluate()]
#' @export
plot_protein_comparison <- function(
    results,
    components    = c("S_length", "S_coverage", "S_count",
                      "S_hydro", "S_charge", "composite_score"),
    show_verdict  = TRUE,
    title         = NULL
) {
  rlang::check_installed("ggplot2", reason = "to use plot_protein_comparison()")

  # ── Parse input ───────────────────────────────────────────────────────────
  score_df <- if (is.data.frame(results)) {
    # batch_evaluate() tibble
    .validate_batch_result(results)
    results
  } else if (.is_named_results_list(results)) {
    labels <- if (!is.null(names(results))) names(results) else
      vapply(results, .result_label, character(1L))
    rows <- lapply(seq_along(results), function(i) {
      s <- results[[i]]$scores[1L, , drop = FALSE]
      s$protein_label <- labels[[i]]
      s
    })
    do.call(rbind, rows)
  } else {
    .abort(
      c(
        "{.arg results} must be a named list of {.fn evaluate_digest} results or a {.fn batch_evaluate} tibble.",
        "x" = "Got {.cls {class(results)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }

  if ("protein_label" %in% names(score_df)) {
    # already set (from named list path above)
  } else {
    score_df$protein_label <- vapply(score_df$protein_id, .tidy_protein_id, character(1L))
  }

  # ── Validate & restrict components ────────────────────────────────────────
  avail_comps <- intersect(components, names(score_df))
  if (length(avail_comps) == 0L) {
    .abort("None of the requested component columns were found in the score data.")
  }
  components <- avail_comps

  # ── Order proteins by composite_score (descending = best on left) ─────────
  if ("composite_score" %in% names(score_df)) {
    ord <- order(score_df$composite_score, decreasing = TRUE)
    score_df <- score_df[ord, , drop = FALSE]
  }
  protein_levels <- score_df$protein_label

  # ── Tidy long format ──────────────────────────────────────────────────────
  long <- do.call(rbind, lapply(components, function(comp) {
    data.frame(
      protein   = factor(score_df$protein_label, levels = protein_levels),
      component = comp,
      value     = as.numeric(score_df[[comp]]),
      stringsAsFactors = FALSE
    )
  }))
  long$component <- factor(long$component, levels = components)

  # ── Verdict annotation data ───────────────────────────────────────────────
  verdict_df <- NULL
  if (show_verdict && "verdict" %in% names(score_df)) {
    verdict_df <- data.frame(
      protein  = factor(score_df$protein_label, levels = protein_levels),
      verdict  = score_df$verdict,
      score    = if ("composite_score" %in% names(score_df)) score_df$composite_score else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # ── Colors: components get brand palette; composite gets brand_dark ───────
  n_comps    <- length(components)
  comp_cols  <- grDevices::hcl.colors(n_comps, palette = "viridis", alpha = 0.85)
  names(comp_cols) <- components
  if ("composite_score" %in% names(comp_cols))
    comp_cols[["composite_score"]] <- .pepvet_pal$brand_dark

  # Verdict label colors
  verdict_colors <- c(
    "Good"     = .pepvet_pal$good,
    "Moderate" = .pepvet_pal$moderate,
    "Poor"     = .pepvet_pal$poor
  )

  # ── Threshold lines ───────────────────────────────────────────────────────
  thresholds <- data.frame(
    y     = c(.get_param("verdict_good"), .get_param("verdict_moderate")),
    color = c(.pepvet_pal$good, .pepvet_pal$moderate),
    label = c("Good", "Moderate")
  )

  auto_title <- if (!is.null(title)) {
    title
  } else {
    enzymes <- unique(unlist(lapply(
      if (is.data.frame(results)) list() else results,
      function(r) r$params$enzyme
    )))
    if (length(enzymes) == 1L)
      paste0(enzymes, "  \u00b7  Protein comparison")
    else
      "Protein comparison"
  }

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    long, ggplot2::aes(x = .data$protein, y = .data$value,
                       fill = .data$component)
  ) +
    # Threshold reference lines
    ggplot2::geom_hline(
      data = thresholds,
      ggplot2::aes(yintercept = .data$y),
      color     = thresholds$color,
      linewidth = 0.5,
      linetype  = "dashed",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_col(
      position  = ggplot2::position_dodge(width = 0.85),
      width     = 0.80,
      color     = "white",
      linewidth = 0.15
    ) +
    ggplot2::scale_fill_manual(
      values = comp_cols,
      name   = NULL,
      labels = c(
        S_length        = "Length",
        S_coverage      = "Coverage",
        S_count         = "Peptide count",
        S_hydro         = "Hydrophobicity",
        S_charge        = "Charge",
        composite_score = "Composite"
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.12),
      breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(mult = c(0, 0))
    )

  # Verdict badges
  if (!is.null(verdict_df)) {
    p <- p + ggplot2::geom_text(
      data = verdict_df,
      ggplot2::aes(x = .data$protein, y = 1.06,
                   label = .data$verdict,
                   color = .data$verdict),
      size        = 2.6,
      fontface    = "bold",
      inherit.aes = FALSE
    ) +
      ggplot2::scale_color_manual(
        values = verdict_colors, name = NULL, guide = "none"
      )
  }

  p +
    ggplot2::labs(
      title    = auto_title,
      subtitle = sprintf("%d proteins  \u00b7  sorted by composite score",
                         nrow(score_df)),
      x        = NULL,
      y        = "Score"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position  = "bottom",
      axis.text.x      = ggplot2::element_text(
        angle = 35, hjust = 1, size = 9),
      plot.title       = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark),
      plot.subtitle    = ggplot2::element_text(size = 9, color = "#666666")
    )
}


# ── plot_enzyme_protein_heatmap ───────────────────────────────────────────────

#' Enzyme × Protein Score Heatmap
#'
#' `plot_enzyme_protein_heatmap()` draws a 2D tile matrix: rows = proteins
#' (y-axis), columns = enzymes (x-axis), fill = score.  It is the only
#' pepVet function that shows the full multi-protein × multi-enzyme comparison
#' space simultaneously.
#'
#' @param results A named list of named lists.  Outer names = protein labels,
#'   inner names = enzyme labels, each leaf = an [evaluate_digest()] result.
#'   **Or** a long-format data.frame / tibble with columns `protein_label`,
#'   `enzyme`, and all component score columns (as returned by iterating
#'   [batch_evaluate()] across multiple enzymes and binding rows with an
#'   `enzyme` column).
#' @param component Character.  Which score to use as the fill.  One of
#'   `c("composite_score","S_length","S_coverage","S_count","S_hydro","S_charge")`.
#'   Defaults to `"composite_score"`.
#' @param show_verdict Logical.  When `TRUE` (default) prints a one-letter
#'   verdict badge (G/M/P) inside each tile.
#' @param title Optional character title.
#'
#' @return A `ggplot` object.
#' @export
plot_enzyme_protein_heatmap <- function(
    results,
    component    = "composite_score",
    show_verdict = TRUE,
    title        = NULL) {

  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots")

  valid_components <- c("composite_score", "S_length", "S_coverage",
                        "S_count", "S_hydro", "S_charge")
  if (!component %in% valid_components) {
    .abort(
      c("!" = "{.arg component} must be one of {.val {valid_components}}.",
        "i" = "Received: {.val {component}}."),
      class = "pepvet_error_invalid_component"
    )
  }

  # ── Build long data.frame ───────────────────────────────────────────────
  if (is.data.frame(results)) {
    # Pre-built long tibble path
    required_cols <- c("protein_label", "enzyme", component)
    missing_cols  <- setdiff(required_cols, names(results))
    if (length(missing_cols) > 0L) {
      .abort(
        c("!" = "When {.arg results} is a data.frame it must contain columns: {.val {required_cols}}.",
          "i" = "Missing: {.val {missing_cols}}."),
        class = "pepvet_error_invalid_digest_result"
      )
    }
    df <- data.frame(
      protein_label = as.character(results$protein_label),
      enzyme        = as.character(results$enzyme),
      score         = as.numeric(results[[component]]),
      stringsAsFactors = FALSE
    )
    # Add verdict column if not present
    if ("verdict" %in% names(results)) {
      df$verdict <- as.character(results$verdict)
    } else if ("composite_score" %in% names(results)) {
      cs <- as.numeric(results$composite_score)
      good_thresh <- .get_param("verdict_good")
      mod_thresh  <- .get_param("verdict_moderate")
      df$verdict <- ifelse(cs >= good_thresh, "Good",
                      ifelse(cs >= mod_thresh, "Moderate", "Poor"))
    } else {
      df$verdict <- NA_character_
    }
  } else if (is.list(results) && !is.null(names(results)) &&
             is.list(results[[1L]]) && !is.null(names(results[[1L]]))) {
    # Nested named list: outer=protein, inner=enzyme
    rows <- list()
    for (prot_name in names(results)) {
      prot_block <- results[[prot_name]]
      for (enz_name in names(prot_block)) {
        leaf <- prot_block[[enz_name]]
        if (!is.list(leaf) || !"scores" %in% names(leaf)) next
        sc  <- as.numeric(leaf$scores[[component]][[1L]])
        vrd <- as.character(leaf$scores$verdict[[1L]])
        rows[[length(rows) + 1L]] <- data.frame(
          protein_label = prot_name,
          enzyme        = enz_name,
          score         = sc,
          verdict       = vrd,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) == 0L) {
      .abort(
        "!" = "No valid leaves found in the nested results list.",
        class = "pepvet_error_invalid_digest_result"
      )
    }
    df <- do.call(rbind, rows)
  } else {
    .abort(
      c("!" = paste0("{.arg results} must be a nested named list (outer=protein,",
                     " inner=enzyme) or a data.frame with columns",
                     " {.code protein_label}, {.code enzyme}, and score columns."),
        "i" = "See {.fn plot_enzyme_protein_heatmap} documentation."),
      class = "pepvet_error_invalid_digest_result"
    )
  }

  # ── Order axes by mean score ────────────────────────────────────────────
  prot_order <- names(sort(
    tapply(df$score, df$protein_label, mean, na.rm = TRUE),
    decreasing = TRUE
  ))
  enz_order <- names(sort(
    tapply(df$score, df$enzyme, mean, na.rm = TRUE),
    decreasing = FALSE  # best enzyme on right
  ))
  df$protein_label <- factor(df$protein_label, levels = rev(prot_order))
  df$enzyme        <- factor(df$enzyme, levels = enz_order)

  # ── Best enzyme per protein (bold border) ───────────────────────────────
  best_df <- do.call(rbind, lapply(split(df, df$protein_label), function(g) {
    g[which.max(g$score), , drop = FALSE]
  }))

  # ── Verdict badge label ──────────────────────────────────────────────────
  df$badge <- ""
  if (show_verdict && "verdict" %in% names(df)) {
    df$badge <- ifelse(df$verdict == "Good", "G",
                  ifelse(df$verdict == "Moderate", "M",
                    ifelse(df$verdict == "Poor", "P", "")))
  }

  # ── Component label for subtitle ────────────────────────────────────────
  comp_label <- switch(component,
    composite_score = "Composite score",
    S_length        = "Length score",
    S_coverage      = "Coverage score",
    S_count         = "Count score",
    S_hydro         = "Hydrophobicity score",
    S_charge        = "Charge score",
    component
  )
  auto_title <- title %||% paste0("Enzyme \u00d7 Protein Score Matrix  \u00b7  ", comp_label)

  p <- ggplot2::ggplot(df,
    ggplot2::aes(x = enzyme, y = protein_label, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    # Best-enzyme border highlight
    ggplot2::geom_tile(
      data = best_df,
      ggplot2::aes(x = enzyme, y = protein_label),
      fill  = NA,
      color = .pepvet_pal$brand_dark,
      linewidth = 1.6
    ) +
    ggplot2::scale_fill_gradient2(
      low      = .pepvet_pal$poor,
      mid      = .pepvet_pal$moderate,
      high     = .pepvet_pal$good,
      midpoint = 0.5,
      limits   = c(0, 1),
      name     = comp_label,
      breaks = c(0, 0.25, 0.50, 0.75, 1.0),
      labels = c("0", "0.25", "0.50", "0.75", "1.00"),
      guide    = ggplot2::guide_colorbar(
        barwidth  = ggplot2::unit(160, "pt"),
        barheight = ggplot2::unit(7, "pt"),
        title.position = "left",
        title.vjust = 0.9
      )
    ) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = "Dark border = best enzyme per protein  \u00b7  Badges: G=Good, M=Moderate, P=Poor",
      x = "Enzyme",
      y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid      = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (show_verdict && any(nchar(df$badge) > 0L)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = badge),
      size     = 3.5,
      fontface = "bold",
      color    = "white"
    )
  }

  p
}
