## pepVet Batch and proteome-scale plots
##
## Functions in this file operate on multi-protein batch results from
## batch_evaluate() or equivalent structures.
##

## plot_proteome_overview

#' Proteome digest overview
#'
#' `plot_proteome_overview()` produces a three-panel portrait of a
#' single-enzyme proteome digest from [batch_evaluate()]:
#'
#' - **(A) Score distribution:** histogram of composite scores, verdict-colored,
#'   with background zone shading for Good (>= 0.65), Moderate (0.40-0.65),
#'   and Poor (< 0.40) regions. A percentage badge is anchored in the Good zone.
#' - **(B) Component profile:** horizontal lollipop chart showing the median of
#'   each component score across all proteins.  Each component uses its
#'   designated color from the pepVet component palette.  Immediately reveals
#'   which digest quality dimension is limiting the proteome.
#' - **(C) Difficulty flags:** 100% stacked horizontal bars showing the
#'   proportion of proteins carrying each difficulty flag (short protein, no
#'   valid peptides, hydrophobic, low-complexity).  Flags are ordered by
#'   prevalence (most common at top).
#'
#' @param batch A tibble returned by [batch_evaluate()], with columns
#'   `protein_id`, `composite_score`, `verdict`, the five component score
#'   columns, and the four difficulty flag columns. If `NULL` or empty, raises
#'   an error.
#' @param title Optional character title for the combined figure. When `NULL`
#'   (default), generates an auto-title with protein count.
#'
#' @details When the input batch lacks difficulty flag columns, panel C shows
#'   an empty placeholder message instead of stacked bars.
#' @section Limitations:
#' At proteome sizes exceeding 2000 proteins the figure panels become dense
#' and text labels may overlap. Consider filtering the batch to a
#' representative subset or increasing the output dimensions.
#' @return A `patchwork` object with three panels: score distribution
#'   histogram (A), component median lollipop chart (B), and difficulty-flag
#'   prevalence bars (C).
#' @seealso [batch_evaluate()] for the upstream evaluation step.
#' @family plot-batch
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'   requireNamespace("patchwork", quietly = TRUE)) {
#'   small <- system.file(
#'     "extdata", "small_proteome_50_proteins.fasta",
#'     package = "pepVet"
#'   )
#'   batch <- batch_evaluate(small, enzyme = "trypsin")
#'   plot_proteome_overview(batch)
#' }
#' @export
plot_proteome_overview <- function(batch, title = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots"
  )
  rlang::check_installed("patchwork",
    reason = "to assemble plot_proteome_overview panels"
  )

  required_cols <- c(
    "protein_id", "composite_score", "verdict",
    "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
  )
  missing_cols <- setdiff(required_cols, names(batch))
  if (length(missing_cols) > 0L) {
    .abort(
      c(
        "!" = "{.arg batch} must be a tibble from {.fn batch_evaluate}.",
        "i" = "Missing columns: {.val {missing_cols}}."
      ),
      class = "pepvet_error_invalid_batch"
    )
  }

  .validate_nonempty(batch, "batch", class = "pepvet_error_invalid_batch")

  batch$composite_score <- as.numeric(batch$composite_score)
  batch$verdict <- factor(
    as.character(batch$verdict),
    levels = c("Good", "Moderate", "Poor")
  )

  verdict_colors <- .pepvet_pal$verdict

  n_total <- nrow(batch)
  n_good <- sum(batch$verdict == "Good", na.rm = TRUE)
  n_moderate <- sum(batch$verdict == "Moderate", na.rm = TRUE)
  n_poor <- sum(batch$verdict == "Poor", na.rm = TRUE)
  pct_good <- round(100 * n_good / n_total)

  ## Panel A: Score distribution
  pa <- ggplot2::ggplot(
    batch,
    ggplot2::aes(x = composite_score, fill = verdict)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = .get_param("verdict_good"), xmax = 1.01, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.50
    ) +
    ggplot2::annotate(
      "rect",
      xmin = .get_param("verdict_moderate"),
      xmax = .get_param("verdict_good"),
      ymin = -Inf,
      ymax = Inf,
      fill = .pepvet_pal$zone_moderate, alpha = 0.45
    ) +
    ggplot2::annotate(
      "rect",
      xmin = -0.01, xmax = .get_param("verdict_moderate"),
      ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$zone_poor, alpha = 0.35
    ) +
    ggplot2::geom_histogram(
      binwidth = 0.05,
      color = "white",
      linewidth = 0.2,
      alpha = .get_param("scatter_alpha"),
      position = "stack"
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_moderate"),
      linetype   = "dashed",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.7,
      alpha      = 0.9
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_good"),
      linetype   = "dashed",
      color      = .pepvet_pal$good,
      linewidth  = 0.7,
      alpha      = 0.9
    ) +
    ggplot2::annotate(
      "text",
      x = (.get_param("verdict_good") + 1.0) / 2,
      y = Inf, hjust = 0.5, vjust = 1.7,
      label = sprintf("Good\n%d%%", pct_good),
      size = 3.2,
      fontface = "bold",
      color = .pepvet_pal$good
    ) +
    ggplot2::scale_fill_manual(
      values = verdict_colors,
      name   = "Verdict",
      guide  = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(add = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::coord_cartesian(xlim = c(-0.01, 1.01)) +
    ggplot2::labs(
      title = "Score Distribution",
      subtitle = sprintf(
        "%d proteins  \u00b7  Good: %d  \u00b7  Moderate: %d  \u00b7  Poor: %d",
        n_total, n_good, n_moderate, n_poor
      ),
      x = "Composite score",
      y = "Proteins"
    ) +
    .pepvet_theme()

  ## Panel B: Component profile (lollipop)
  ## Component colors follow the pepVet component-score color map
  comp_cols <- c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge")
  comp_names <- c(
    S_length   = "Length",
    S_coverage = "Coverage",
    S_count    = "Count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge"
  )
  comp_colors <- .pepvet_pal$component

  comp_medians <- vapply(
    comp_cols,
    function(cc) stats::median(as.numeric(batch[[cc]]), na.rm = TRUE),
    numeric(1L)
  )
  comp_df <- data.frame(
    comp_id = comp_cols,
    label = factor(
      comp_names[comp_cols],
      levels = rev(comp_names[comp_cols]) # bottom-to-top order
    ),
    score = comp_medians,
    stringsAsFactors = FALSE
  )

  pb <- ggplot2::ggplot(comp_df, ggplot2::aes(y = label, x = score)) +
    ggplot2::annotate(
      "rect",
      xmin = .get_param("verdict_good"), xmax = 1.02, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.55
    ) +
    ggplot2::annotate(
      "rect",
      xmin = .get_param("verdict_moderate"),
      xmax = .get_param("verdict_good"),
      ymin = -Inf,
      ymax = Inf,
      fill = .pepvet_pal$zone_moderate, alpha = 0.40
    ) +
    ggplot2::annotate(
      "rect",
      xmin = 0, xmax = .get_param("verdict_moderate"), ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$zone_poor, alpha = 0.30
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_moderate"),
      linetype   = "dashed",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.5,
      alpha      = 0.75
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_good"),
      linetype   = "dashed",
      color      = .pepvet_pal$good,
      linewidth  = 0.5,
      alpha      = 0.75
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0, xend = score,
        y = label, yend = label,
        color = comp_id
      ),
      linewidth = 2.0,
      alpha = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = comp_id),
      size = 9,
      alpha = 0.92
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", score)),
      size = 2.8,
      fontface = "bold",
      color = "white"
    ) +
    ggplot2::scale_color_manual(values = comp_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.0),
      breaks = c(
        0, .get_param("verdict_moderate"),
        .get_param("verdict_good"), 1.0
      ),
      labels = c(
        "0",
        format(.get_param("verdict_moderate")),
        format(.get_param("verdict_good")),
        "1.0"
      ),
      expand = ggplot2::expansion(add = c(0.02, 0.05))
    ) +
    ggplot2::labs(
      title = "Component Profile",
      subtitle = paste(
        "Median component score across all proteins",
        "Color = one digest quality dimension",
        sep = "\n"
      ),
      x = "Median score",
      y = NULL
    ) +
    .pepvet_theme()

  ## Panel C: Difficulty flags
  flag_cols <- c(
    "flag_short_protein", "flag_no_valid_peptides",
    "flag_hydrophobic",   "flag_low_complexity"
  )
  flag_labels <- c(
    flag_short_protein     = "Short protein  (<100 aa)",
    flag_no_valid_peptides = "No valid peptides",
    flag_hydrophobic       = "Hydrophobic proteome",
    flag_low_complexity    = "Low complexity sequence"
  )
  present_flags <- flag_cols[flag_cols %in% names(batch)]

  if (length(present_flags) > 0L) {
    flag_pcts <- vapply(present_flags, function(fc) {
      vals <- as.logical(batch[[fc]])
      vals[is.na(vals)] <- FALSE
      sum(vals) / n_total
    }, numeric(1L))

    ordered_flags <- present_flags[order(flag_pcts)]

    flag_long <- .bind_rows(lapply(ordered_flags, function(fc) {
      vals <- as.logical(batch[[fc]])
      vals[is.na(vals)] <- FALSE
      n_flagged <- sum(vals)
      pct_flag <- 100 * n_flagged / n_total
      rbind(
        data.frame(
          flag = flag_labels[[fc]], category = "Flagged",
          pct = pct_flag, n = n_flagged,
          stringsAsFactors = FALSE
        ),
        data.frame(
          flag = flag_labels[[fc]], category = "Not flagged",
          pct = 100 - pct_flag, n = n_total - n_flagged,
          stringsAsFactors = FALSE
        )
      )
    }))
    flag_long$flag <- factor(
      flag_long$flag,
      levels = flag_labels[ordered_flags]
    )
    flag_long$category <- factor(
      flag_long$category,
      levels = c("Not flagged", "Flagged")
    )

    pc <- ggplot2::ggplot(
      flag_long,
      ggplot2::aes(y = flag, x = pct, fill = category)
    ) +
      ggplot2::geom_col(alpha = .get_param("scatter_alpha"), width = 0.55) +
      ggplot2::geom_text(
        data = flag_long[flag_long$pct >= 5, ],
        ggplot2::aes(
          label = sprintf("%d\u00a0%%  (%d)", round(pct), n)
        ),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 2.8,
        fontface = "bold",
        color = "white"
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          Flagged       = .pepvet_pal$poor,
          "Not flagged" = .pepvet_pal$good
        ),
        name = NULL
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100),
        labels = c("0%", "25%", "50%", "75%", "100%"),
        expand = ggplot2::expansion(add = c(0, 0))
      ) +
      ggplot2::labs(
        title = "Difficulty Flags",
        subtitle = paste0(
          "Proportion of proteins carrying each digest difficulty flag  ",
          "\u00b7  Ordered by prevalence"
        ),
        x = "Proportion of proteome",
        y = NULL
      ) +
      .pepvet_theme()
  } else {
    # Fallback: empty placeholder when no flag columns present
    pc <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.5, y = 0.5,
        label = "Difficulty flag columns not available",
        size = 3.5,
        color = .pepvet_pal$brand_dark
      ) +
      ggplot2::theme_void()
  }

  ## Compose
  auto_title <- if (is.null(title)) {
    sprintf("Proteome Digest Overview  \u00b7  %d proteins", n_total)
  } else {
    title
  }

  (pa | pb) / pc +
    patchwork::plot_layout(heights = c(1.35, 1)) +
    patchwork::plot_annotation(
      title = auto_title,
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face   = "bold",
          size   = .get_param("patchwork_title_size"),
          color  = .pepvet_pal$brand_dark,
          margin = ggplot2::margin(b = 8)
        )
      )
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face  = "bold",
        size  = .get_param("patchwork_tag_size"),
        color = .pepvet_pal$brand
      )
    )
}

## plot_batch_comparison

#' Multi-enzyme proteome comparison
#'
#' `plot_batch_comparison()` produces a four-panel side-by-side comparison of
#' enzyme performance across a full proteome, using output from
#' [batch_compare_enzymes()]:
#'
#' - **(A) Verdict summary:** 100% stacked horizontal bars showing the
#'   Good/Moderate/Poor verdict breakdown per enzyme.  Enzymes are ordered
#'   by descending Good%, and the top-performing enzyme is marked with a star.
#' - **(B) Score distributions:** horizontal violin plots for each enzyme,
#'   showing the full distribution of composite scores.  An IQR boxplot is
#'   overlaid on each violin.  Violin fill color reflects the enzyme's median
#'   verdict.
#' - **(C) Component heatmap:** median component scores in an
#'   enzyme-by-component
#'   grid, filled by the verdict gradient (red to amber to green). Reveals which
#'   digest quality dimension differentiates the enzymes.
#' - **(D) Per-protein win rate:** bar chart showing the proportion of proteins
#'   for which each enzyme achieves the highest composite score. The starred
#'   enzyme
#'   matches the recommendation in panel A.
#'
#' @param comparison A `pepvet_batch_comparison` tibble returned by
#'   [batch_compare_enzymes()], with columns `protein_id`, `enzyme`,
#'   `composite_score`, `verdict`, and the five component score columns. If
#'   `NULL` or empty, raises an error.
#' @param title Optional character title for the combined figure. When `NULL`
#'   (default), generates an auto-title with protein and enzyme counts.
#'
#' @section Limitations:
#' At proteome sizes exceeding 2000 proteins the figure panels become dense
#' and text labels may overlap. Consider filtering the comparison to a
#' representative subset or increasing the output dimensions.
#' @return A `patchwork` object with four panels: verdict summary bars (A),
#'   score distribution violins (B), component-score heatmap (C), and
#'   per-protein win-rate bars (D).
#' @seealso [batch_compare_enzymes()] for the upstream comparison step.
#' @family plot-batch
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'   requireNamespace("patchwork", quietly = TRUE)) {
#'   small <- system.file(
#'     "extdata", "small_proteome_50_proteins.fasta",
#'     package = "pepVet"
#'   )
#'   comp <- batch_compare_enzymes(small, enzymes = c("trypsin", "lysc"))
#'   plot_batch_comparison(comp)
#' }
#' @export
plot_batch_comparison <- function(comparison, title = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots"
  )
  rlang::check_installed("patchwork",
    reason = "to assemble plot_batch_comparison panels"
  )

  required_cols <- c(
    "protein_id", "enzyme", "composite_score", "verdict",
    "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
  )
  missing_cols <- setdiff(required_cols, names(comparison))
  if (length(missing_cols) > 0L) {
    .abort(
      c(
        "!" = paste0(
          "{.arg comparison} must be a tibble from {.fn batch_compare_enzymes}."
        ),
        "i" = "Missing columns: {.val {missing_cols}}."
      ),
      class = "pepvet_error_invalid_batch"
    )
  }

  .validate_nonempty(
    comparison, "comparison",
    class = "pepvet_error_invalid_batch"
  )

  comparison$composite_score <- as.numeric(comparison$composite_score)
  comparison$enzyme <- as.character(comparison$enzyme)
  comparison$verdict <- factor(
    as.character(comparison$verdict),
    levels = c("Good", "Moderate", "Poor")
  )

  enz_levels <- if (inherits(comparison, "pepvet_batch_comparison")) {
    attr(comparison, "enzymes")
  } else {
    unique(comparison$enzyme)
  }
  n_enz <- length(enz_levels)
  n_proteins <- length(unique(comparison$protein_id))

  verdict_colors <- .pepvet_pal$verdict

  ## Per-enzyme verdict summary table
  enz_stats <- .bind_rows(lapply(enz_levels, function(enz) {
    sub <- comparison[comparison$enzyme == enz, ]
    n <- nrow(sub)
    n_g <- sum(sub$verdict == "Good", na.rm = TRUE)
    n_m <- sum(sub$verdict == "Moderate", na.rm = TRUE)
    n_p <- sum(sub$verdict == "Poor", na.rm = TRUE)
    data.frame(
      enzyme = enz,
      n_total = n,
      n_good = n_g,
      n_moderate = n_m,
      n_poor = n_p,
      pct_good = if (n > 0L) 100 * n_g / n else NA_real_,
      pct_moderate = if (n > 0L) 100 * n_m / n else NA_real_,
      pct_poor = if (n > 0L) 100 * n_p / n else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  enz_order <- enz_stats$enzyme[order(-enz_stats$pct_good, na.last = TRUE)]
  best_enzyme <- enz_order[[1L]]

  ## Panel A: Verdict summary (stacked 100% horizontal bars)
  enz_long <- .bind_rows(lapply(enz_levels, function(enz) {
    row <- enz_stats[enz_stats$enzyme == enz, ]
    data.frame(
      enzyme = enz,
      verdict = c("Poor", "Moderate", "Good"),
      pct = c(row$pct_poor, row$pct_moderate, row$pct_good),
      n = c(row$n_poor, row$n_moderate, row$n_good),
      stringsAsFactors = FALSE
    )
  }))
  enz_long$verdict <- factor(
    enz_long$verdict,
    levels = c("Poor", "Moderate", "Good")
  )
  enz_long$enzyme <- factor(enz_long$enzyme, levels = rev(enz_order))

  y_labels <- setNames(enz_order, enz_order)
  y_labels[[best_enzyme]] <- paste0(best_enzyme, "  \u2605")

  pa <- ggplot2::ggplot(
    enz_long,
    ggplot2::aes(y = enzyme, x = pct, fill = verdict)
  ) +
    ggplot2::geom_col(
      alpha = .get_param("scatter_alpha"),
      width = 0.65,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      data = enz_long[!is.na(enz_long$pct) & enz_long$pct >= 7, ],
      ggplot2::aes(label = sprintf("%.0f%%", pct)),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.8,
      fontface = "bold",
      color = "white"
    ) +
    ggplot2::scale_fill_manual(
      values = verdict_colors,
      name = "Verdict",
      breaks = c("Good", "Moderate", "Poor"),
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1),
        reverse      = FALSE
      )
    ) +
    ggplot2::scale_y_discrete(labels = y_labels) +
    ggplot2::scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Verdict Summary",
      subtitle = sprintf(
        "%d proteins  \u00b7  %d enzymes  \u00b7  \u2605 = recommended",
        n_proteins, n_enz
      ),
      x = "Proportion of proteome",
      y = NULL
    ) +
    .pepvet_theme()

  ## Panel B: Score distributions (horizontal violins)
  median_verdict_per_enz <- vapply(enz_levels, function(enz) {
    ms <- stats::median(
      comparison$composite_score[comparison$enzyme == enz],
      na.rm = TRUE
    )
    if (is.na(ms) || ms >= .get_param("verdict_good")) {
      return("Good")
    }
    if (ms >= .get_param("verdict_moderate")) {
      return("Moderate")
    }
    "Poor"
  }, character(1L))
  violin_fills <- setNames(
    verdict_colors[median_verdict_per_enz], enz_levels
  )

  comp_viol <- comparison
  comp_viol$enzyme <- factor(comp_viol$enzyme, levels = rev(enz_order))

  pb <- ggplot2::ggplot(
    comp_viol,
    ggplot2::aes(x = composite_score, y = enzyme, fill = enzyme)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = .get_param("verdict_good"), xmax = 1.01, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.40
    ) +
    ggplot2::annotate(
      "rect",
      xmin = 0, xmax = .get_param("verdict_moderate"), ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$zone_poor, alpha = 0.25
    ) +
    ggplot2::geom_violin(
      alpha = 0.78,
      trim  = TRUE,
      scale = "width",
      color = NA
    ) +
    ggplot2::geom_boxplot(
      width         = 0.20,
      outlier.size  = 0.8,
      outlier.alpha = 0.30,
      fill          = "white",
      color         = .pepvet_pal$brand_dark,
      linewidth     = 0.45,
      alpha         = 0.75
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_moderate"),
      linetype   = "dashed",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.55,
      alpha      = 0.85
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_good"),
      linetype   = "dashed",
      color      = .pepvet_pal$good,
      linewidth  = 0.55,
      alpha      = 0.85
    ) +
    ggplot2::scale_fill_manual(values = violin_fills, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = c(
        0, .get_param("verdict_moderate"),
        .get_param("verdict_good"), 1.0
      ),
      labels = c(
        "0",
        format(.get_param("verdict_moderate")),
        format(.get_param("verdict_good")),
        "1.0"
      ),
      expand = ggplot2::expansion(add = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = "Score Distributions",
      subtitle = sprintf(
        paste0(
          "Violin + IQR box  \u00b7  Fill = median verdict  ",
          "\u00b7  Dashed lines = %s / %s"
        ),
        .get_param("verdict_moderate"), .get_param("verdict_good")
      ),
      x = "Composite score",
      y = NULL
    ) +
    .pepvet_theme()

  ## Panel C: Component heatmap (enzyme \u00d7 component)
  comp_cols <- c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge")
  comp_labels <- c(
    S_length   = "Length",
    S_coverage = "Coverage",
    S_count    = "Count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge"
  )

  tile_rows <- .bind_rows(lapply(enz_levels, function(enz) {
    sub <- comparison[comparison$enzyme == enz, ]
    meds <- vapply(
      comp_cols,
      function(cc) stats::median(as.numeric(sub[[cc]]), na.rm = TRUE),
      numeric(1L)
    )
    data.frame(
      enzyme = enz,
      component = comp_labels[comp_cols],
      score = meds,
      stringsAsFactors = FALSE
    )
  }))
  tile_rows$enzyme <- factor(tile_rows$enzyme, levels = rev(enz_order))
  tile_rows$component <- factor(
    tile_rows$component,
    levels = comp_labels[comp_cols]
  )

  verdict_grad <- grDevices::colorRampPalette(
    c(.pepvet_pal$poor, .pepvet_pal$heatmap_mid, .pepvet_pal$good)
  )(100)

  pc <- ggplot2::ggplot(
    tile_rows,
    ggplot2::aes(x = component, y = enzyme, fill = score)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1.0) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", score)),
      size = 3.0,
      fontface = "bold",
      color = .pepvet_pal$brand_dark
    ) +
    ggplot2::scale_fill_gradientn(
      colors = verdict_grad,
      limits = c(0, 1),
      breaks = c(
        0, .get_param("verdict_moderate"),
        .get_param("verdict_good"), 1.0
      ),
      labels = c(
        "0",
        format(.get_param("verdict_moderate")),
        format(.get_param("verdict_good")),
        "1.0"
      ),
      name = "Median\nscore",
      guide = ggplot2::guide_colorbar(
        barwidth       = ggplot2::unit(10, "pt"),
        barheight      = ggplot2::unit(60, "pt"),
        title.position = "top"
      )
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::labs(
      title = "Component Score Profile",
      subtitle = paste0(
        "Median component score per enzyme  ",
        "\u00b7  Color = verdict gradient"
      ),
      x = NULL,
      y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "right"
    )

  ## Panel D: Per-protein win rate
  # For each protein, find the enzyme with the highest composite score
  cs_df <- data.frame(
    protein_id = as.character(comparison$protein_id),
    enzyme = as.character(comparison$enzyme),
    composite_score = comparison$composite_score,
    stringsAsFactors = FALSE
  )
  pids <- unique(cs_df$protein_id)
  best_for <- vapply(pids, function(pid) {
    sub <- cs_df[cs_df$protein_id == pid, ]
    sub$enzyme[which.max(sub$composite_score)]
  }, character(1L))

  win_counts <- tabulate(
    factor(best_for, levels = enz_order),
    nbins = length(enz_order)
  )
  win_df <- data.frame(
    enzyme = enz_order,
    n_wins = win_counts,
    pct_wins = 100 * win_counts / length(pids),
    is_best = enz_order == best_enzyme,
    stringsAsFactors = FALSE
  )
  win_df$enzyme <- factor(win_df$enzyme, levels = rev(enz_order))

  pd <- ggplot2::ggplot(
    win_df,
    ggplot2::aes(y = enzyme, x = pct_wins, fill = is_best)
  ) +
    ggplot2::geom_col(
      alpha = .get_param("scatter_alpha"),
      width = 0.65, color = NA
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x     = pct_wins + 0.8,
        label = sprintf("%.0f%%  (%d)", pct_wins, n_wins)
      ),
      hjust = 0,
      size = 2.8,
      color = .pepvet_pal$brand_dark,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = .pepvet_pal$good, `FALSE` = .pepvet_pal$brand_light),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 115),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Per-Protein Win Rate",
      subtitle = paste0(
        "Enzyme with highest score per protein  ",
        "\u00b7  Ties = first-listed enzyme"
      ),
      x = "Proteins (%)",
      y = NULL
    ) +
    .pepvet_theme()

  ## Compose
  auto_title <- if (is.null(title)) {
    sprintf(
      "Proteome Enzyme Comparison  \u00b7  %d proteins  \u00b7  %d enzymes",
      n_proteins, n_enz
    )
  } else {
    title
  }

  (pa | pb) / (pc | pd) +
    patchwork::plot_layout(
      widths  = c(1, 1),
      heights = c(1.1, 1)
    ) +
    patchwork::plot_annotation(
      title = auto_title,
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face   = "bold",
          size   = .get_param("patchwork_title_size"),
          color  = .pepvet_pal$brand_dark,
          margin = ggplot2::margin(b = 10)
        )
      )
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face  = "bold",
        size  = .get_param("patchwork_tag_size"),
        color = .pepvet_pal$brand
      )
    )
}
