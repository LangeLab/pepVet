# ── pepVet Comparison plots ────────────────────────────────────────────────
# ── plot_enzyme_comparison ────────────────────────────────────────────────────

#' Validate comparison tibble for plotting (internal helper)
#'
#' @noRd
.validate_comparison_for_plot <- function(comparison) {
  if (!is.data.frame(comparison)) {
    cli::cli_abort(
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
    cli::cli_abort(
      c(
        "{.arg comparison} is missing required columns.",
        "x" = "Missing: {.field {missing}}.",
        "i" = "Pass the tibble returned by {.fn compare_digests} directly."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  if (nrow(comparison) < 2L) {
    cli::cli_abort(
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
#' (green >= 0.70, amber 0.40-0.69, red < 0.40). When `recommend = TRUE` a
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
    cli::cli_abort(
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
    ifelse(x >= 0.70, .pepvet_pal$good,
    ifelse(x >= 0.40, .pepvet_pal$moderate, .pepvet_pal$poor))
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
  # Panel A — component scores grouped bar chart
  # ═══════════════════════════════════════════════════════════════════════════
  pa <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$value, y = .data$enzyme, fill = .data$score_name)
  ) +
    # Good-region shading
    ggplot2::annotate("rect",
      xmin = 0.70, xmax = 1.0, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.55
    ) +
    # Threshold reference lines
    ggplot2::geom_vline(
      xintercept = c(0.40, 0.70),
      color = .pepvet_pal$separator, linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.68, alpha = 0.88
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
  # Panel B — composite score lollipop + verdict badge
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
      xmin = 0.70, xmax = 1.1, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.55
    ) +
    ggplot2::geom_vline(
      xintercept = c(0.40, 0.70),
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


# ── plot_score_radar ──────────────────────────────────────────────────────────

#' Build a closed polygon data.frame for one radar layer (internal helper)
#'
#' Converts a named numeric vector of scores (values in [0, 1]) into
#' Cartesian (x, y) coordinates for a regular n-gon centred at the origin.
#' Axes are placed at equal angular intervals starting from the top (pi/2)
#' and proceeding clockwise.
#'
#' @param values  Named numeric vector; names must match `axis_names`.
#' @param axis_names Character vector giving axis order.
#' @param group   Label identifying this polygon (enzyme name, etc.).
#' @return data.frame with columns `x`, `y`, `group`.
#' @noRd
.radar_polygon <- function(values, axis_names, group) {
  n     <- length(axis_names)
  angs  <- pi / 2 - 2 * pi * (seq_len(n) - 1L) / n
  vals  <- as.numeric(values[axis_names])
  vals[is.na(vals)] <- 0
  # Close the polygon
  data.frame(
    x     = c(vals * cos(angs), vals[[1L]] * cos(angs[[1L]])),
    y     = c(vals * sin(angs), vals[[1L]] * sin(angs[[1L]])),
    group = group,
    stringsAsFactors = FALSE
  )
}

#' Build a regular-polygon grid ring (internal helper)
#'
#' @param r        Radius of the ring (0 < r <= 1).
#' @param n_axes   Number of axes.
#' @return data.frame with columns `x`, `y`, `r`.
#' @noRd
.radar_ring <- function(r, n_axes) {
  angs <- pi / 2 - 2 * pi * (seq_len(n_axes) - 1L) / n_axes
  data.frame(
    x = c(r * cos(angs), r * cos(angs[[1L]])),
    y = c(r * sin(angs), r * sin(angs[[1L]])),
    r = r,
    stringsAsFactors = FALSE
  )
}

#' Validate input for plot_score_radar (internal helper)
#'
#' Accepts either a compare_digests tibble (data.frame with `enzyme` column)
#' or a single evaluate_digest list.  Returns a normalised one-row-per-enzyme
#' data.frame with at minimum columns `enzyme` and `composite_score`.
#' @noRd
.prepare_radar_data <- function(comparison) {
  # Single evaluate_digest() list
  if (is.list(comparison) && !is.data.frame(comparison) &&
      all(c("scores", "peptides", "params") %in% names(comparison))) {
    df         <- comparison$scores
    df$enzyme  <- comparison$params$enzyme
    return(df)
  }
  # compare_digests tibble
  if (!is.data.frame(comparison)) {
    cli::cli_abort(
      c(
        "{.arg comparison} must be a tibble from {.fn compare_digests} or a
         list from {.fn evaluate_digest}.",
        "x" = "Got {.cls {class(comparison)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  required <- c("enzyme", "composite_score")
  missing  <- setdiff(required, names(comparison))
  if (length(missing) > 0L) {
    cli::cli_abort(
      c(
        "{.arg comparison} is missing required columns.",
        "x" = "Missing: {.field {missing}}."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }
  comparison
}


#' Score Radar (Spider) Chart
#'
#' `plot_score_radar()` draws a radar / spider chart of component scores for
#' one or more enzymes, overlaid as filled polygons on a regular grid.  The
#' polygon *shape* reveals tradeoffs that bar charts hide: an enzyme with high
#' coverage but low peptide count produces a very different silhouette from one
#' that scores evenly across all components.
#'
#' The chart is built from pure ggplot2 Cartesian geometry — no additional
#' packages beyond those already in `Suggests` are required.  Axes radiate
#' from the centre at equal angular intervals, starting from the top and
#' proceeding clockwise.  Concentric reference polygons at 0.25, 0.50, 0.75,
#' and 1.00 act as a grid.
#'
#' @param comparison A tibble returned by [compare_digests()], **or** a named
#'   list returned by [evaluate_digest()] (single-enzyme radar).  When a
#'   single `evaluate_digest()` result is supplied the polygon is labeled with
#'   the enzyme name stored in `result$params$enzyme`.
#' @param scores Character vector of component-score column names to use as
#'   radar axes.  Any column absent from `comparison` is silently dropped.
#'   Defaults to all five standard component scores.
#' @param title Optional character string for the plot title.  Auto-generated
#'   when `NULL` (default).
#' @param legend_ncol Integer.  Number of columns in the enzyme legend.
#'   Defaults to `min(n_enzymes, 4)` so up to four enzymes sit in one row and
#'   larger sets wrap onto additional rows automatically.
#' @param legend_nrow Integer.  Number of rows in the enzyme legend.  Leave as
#'   `NULL` (default) to let `legend_ncol` control wrapping.
#'
#' @return A `ggplot` object that can be printed, further customised, or saved
#'   with [ggplot2::ggsave()].
#'
#' @details Axis labels use short human-readable names (e.g. "Coverage",
#'   "Length").  The label for each axis is placed just outside the outer ring
#'   (radius 1.18) and auto-justified toward the centre so it does not clip
#'   against the panel border.
#'
#'   When more than five enzymes are compared, enzyme colors cycle through the
#'   eight-color JCO-inspired palette.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   comp <- compare_digests(bsa_path,
#'                           enzymes = c("trypsin", "lysc",
#'                                       "glutamyl endopeptidase"))
#'   p <- plot_score_radar(comp)
#'   print(p)
#' }
#'
#' @seealso [compare_digests()], [plot_enzyme_comparison()],
#'   [plot_digest_profile()]
#' @export
plot_score_radar <- function(
    comparison,
    scores       = c("S_coverage", "S_length", "S_count", "S_hydro", "S_charge"),
    title        = NULL,
    legend_ncol  = NULL,
    legend_nrow  = NULL
) {
  rlang::check_installed("ggplot2", reason = "to use plot_score_radar()")

  df     <- .prepare_radar_data(comparison)
  scores <- intersect(scores, names(df))
  if (length(scores) == 0L) {
    cli::cli_abort(
      c(
        "None of the requested {.arg scores} columns found in {.arg comparison}.",
        "i" = "Available: {.field {names(df)}}."
      ),
      class = "pepvet_error_invalid_comparison"
    )
  }

  # ── Axis display names ───────────────────────────────────────────────────
  score_labels <- c(
    S_coverage = "Coverage",
    S_length   = "Length",
    S_count    = "Count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge"
  )
  axis_labels <- ifelse(scores %in% names(score_labels),
                        score_labels[scores], scores)

  n_axes <- length(scores)
  # Axes start at top (pi/2) and run clockwise
  angs <- pi / 2 - 2 * pi * (seq_len(n_axes) - 1L) / n_axes

  # ── Enzyme color palette (8 colours, cycling for > 8 enzymes) ────────────
  enzyme_palette <- c(
    "#2C5F8A", "#27AE60", "#E8A838", "#8B5E99",
    "#4AAFB0", "#C94040", "#1A3D5C", "#7BAED4"
  )
  enzymes    <- as.character(df$enzyme)
  n_enz      <- length(enzymes)
  enz_colors <- enzyme_palette[((seq_len(n_enz) - 1L) %% length(enzyme_palette)) + 1L]
  names(enz_colors) <- enzymes

  # ── Legend layout: default to ≤ 4 per row so long lists wrap cleanly ─────
  if (is.null(legend_ncol) && is.null(legend_nrow)) {
    legend_ncol <- min(n_enz, 4L)
  }

  # ── Grid rings (closed n-gon polygons at r = 0.25 / 0.50 / 0.75 / 1.00) ─
  grid_r    <- c(0.25, 0.50, 0.75, 1.00)
  grid_data <- do.call(rbind, lapply(grid_r, .radar_ring, n_axes = n_axes))

  # ── Axis spokes: from centre (0,0) outward to the unit ring ─────────────
  spoke_df <- data.frame(
    x    = cos(angs),
    y    = sin(angs),
    xend = rep(0, n_axes),
    yend = rep(0, n_axes)
  )

  # ── Axis labels at r = 1.20; align text toward the outer edge ───────────
  label_r   <- 1.20
  hjust_map <- ifelse(cos(angs) > 0.1, 0, ifelse(cos(angs) < -0.1, 1, 0.5))
  vjust_map <- ifelse(sin(angs) > 0.1, 0, ifelse(sin(angs) < -0.1, 1, 0.5))

  axis_label_df <- data.frame(
    x      = label_r * cos(angs),
    y      = label_r * sin(angs),
    label  = axis_labels,
    hjust  = hjust_map,
    vjust  = vjust_map,
    stringsAsFactors = FALSE
  )

  # ── Grid level labels: placed along the midpoint angle between spokes 1&2
  # so they never sit on top of a spoke or the polygon edge ─────────────────
  mid_ang       <- (angs[[1L]] + angs[[min(2L, n_axes)]]) / 2
  grid_label_df <- data.frame(
    x     = grid_r * cos(mid_ang) * 0.92,
    y     = grid_r * sin(mid_ang) * 0.92,
    label = as.character(grid_r),
    stringsAsFactors = FALSE
  )

  # ── Enzyme polygons ──────────────────────────────────────────────────────
  poly_data <- do.call(rbind, lapply(seq_len(n_enz), function(i) {
    row  <- df[i, , drop = FALSE]
    vals <- unlist(row[scores])
    grp  <- enzymes[[i]]
    .radar_polygon(vals, scores, grp)
  }))
  poly_data$enzyme <- poly_data$group

  # ── Vertex points for each enzyme ────────────────────────────────────────
  vertex_data <- do.call(rbind, lapply(seq_len(n_enz), function(i) {
    row  <- df[i, , drop = FALSE]
    vals <- as.numeric(row[scores])
    data.frame(
      x      = vals * cos(angs),
      y      = vals * sin(angs),
      enzyme = enzymes[[i]],
      stringsAsFactors = FALSE
    )
  }))

  # ── Auto title ───────────────────────────────────────────────────────────
  auto_title <- if (!is.null(title)) {
    title
  } else if ("protein_id" %in% names(df)) {
    pid <- df$protein_id[[1L]]
    paste0(.tidy_protein_id(pid), "  \u00b7  Score radar")
  } else {
    "Score radar"
  }

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +

    # Background: shade the "Good" zone (r 0.70 – 1.00) in the same green
    # as the rest of the package; mask below 0.70 back to white
    ggplot2::geom_polygon(
      data  = .radar_ring(1.00, n_axes),
      ggplot2::aes(x = .data$x, y = .data$y),
      fill  = .pepvet_pal$shade, alpha = 0.40, color = NA
    ) +
    ggplot2::geom_polygon(
      data  = .radar_ring(0.70, n_axes),
      ggplot2::aes(x = .data$x, y = .data$y),
      fill  = "white", alpha = 1, color = NA
    ) +

    # Grid rings (drawn on top of shading so they show as ring borders)
    ggplot2::geom_polygon(
      data      = grid_data,
      ggplot2::aes(x = .data$x, y = .data$y, group = factor(.data$r)),
      fill      = NA,
      color     = .pepvet_pal$separator,
      linewidth = 0.30
    ) +

    # Axis spokes (dotted, same separator colour)
    ggplot2::geom_segment(
      data = spoke_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      color     = .pepvet_pal$separator,
      linewidth = 0.30,
      linetype  = "dotted"
    ) +

    # Grid-level reference numbers along the inter-spoke gap
    ggplot2::geom_text(
      data = grid_label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size  = 2.0,
      color = "#BBBBBB",
      hjust = 0.5,
      vjust = 0.5
    ) +

    # Enzyme score polygons (filled + outlined)
    ggplot2::geom_polygon(
      data = poly_data,
      ggplot2::aes(x     = .data$x,
                   y     = .data$y,
                   group = .data$enzyme,
                   fill  = .data$enzyme,
                   color = .data$enzyme),
      alpha     = 0.18,
      linewidth = 0.9
    ) +

    # Vertex dots
    ggplot2::geom_point(
      data = vertex_data,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$enzyme),
      size = 2.4
    ) +

    # Axis names: bold brand-dark text; hjust/vjust push each label away from
    # the chart centre toward its own edge for clean separation
    ggplot2::geom_text(
      data = axis_label_df,
      ggplot2::aes(x     = .data$x,
                   y     = .data$y,
                   label = .data$label,
                   hjust = .data$hjust,
                   vjust = .data$vjust),
      size     = 3.1,
      fontface = "bold",
      color    = .pepvet_pal$brand_dark
    ) +

    # ── Scales ──────────────────────────────────────────────────────────────
    ggplot2::scale_fill_manual(values = enz_colors, name = "Enzyme") +
    ggplot2::scale_color_manual(values = enz_colors, name = "Enzyme") +

    # Merge fill + color into one legend; respect caller-supplied ncol/nrow
    ggplot2::guides(
      fill  = ggplot2::guide_legend(ncol = legend_ncol, nrow = legend_nrow,
                                    byrow = TRUE),
      color = ggplot2::guide_legend(ncol = legend_ncol, nrow = legend_nrow,
                                    byrow = TRUE)
    ) +

    ggplot2::coord_equal(
      xlim = c(-1.50, 1.50),
      ylim = c(-1.15, 1.50)
    ) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = sprintf(
        "%d component%s  \u00b7  %d enzyme%s  \u00b7  shaded region \u2265 0.70 (Good)",
        n_axes, if (n_axes == 1L) "" else "s",
        n_enz,  if (n_enz  == 1L) "" else "s"
      ),
      x = NULL, y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.text       = ggplot2::element_blank(),
      axis.ticks      = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(face = "bold",
                                              color = .pepvet_pal$brand_dark),
      plot.title      = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark
      )
    )

  p
}



# ── plot_protein_comparison ───────────────────────────────────────────────────

#' Protein Comparison — Component Scores Across Multiple Proteins
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
    cli::cli_abort(
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
    cli::cli_abort("None of the requested component columns were found in the score data.")
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
    y     = c(0.7, 0.4),
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
        S_length        = "Length (S\u2097)",
        S_coverage      = "Coverage (S\u2099)",
        S_count         = "Count (S\u2095)",
        S_hydro         = "Hydrophobicity (S\u02b0)",
        S_charge        = "Charge (S\u2c)",
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
    cli::cli_abort(
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
      cli::cli_abort(
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
      df$verdict <- ifelse(cs >= 0.7, "Good",
                      ifelse(cs >= 0.4, "Moderate", "Poor"))
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
      cli::cli_abort(
        "!" = "No valid leaves found in the nested results list.",
        class = "pepvet_error_invalid_digest_result"
      )
    }
    df <- do.call(rbind, rows)
  } else {
    cli::cli_abort(
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
      guide    = ggplot2::guide_colorbar(
        barwidth  = ggplot2::unit(80, "pt"),
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
