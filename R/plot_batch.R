# ── pepVet Batch & proteome-scale plots ─────────────────────────────────────
#
# Functions in this file operate on multi-protein batch results from
# batch_evaluate() or equivalent structures.
# ─────────────────────────────────────────────────────────────────────────────

# ── plot_batch_summary ────────────────────────────────────────────────────────

#' Batch Digest Score Summary
#'
#' `plot_batch_summary()` produces a two-panel overview for proteome-scale
#' results from [batch_evaluate()]:
#'
#' - **(A) Score distribution:** histogram of composite scores across all
#'   proteins, colored by verdict (Good=green, Moderate=amber, Poor=red).
#'   Vertical dashed lines at the 0.40 and 0.70 thresholds.
#' - **(B) Score vs. length scatter:** composite score on y, protein length on
#'   x.  Points colored by verdict.  Poor-verdict proteins labeled with their
#'   tidy protein ID.  A LOESS trend line reveals whether longer proteins
#'   systematically score differently.
#'
#' @param batch A tibble returned by [batch_evaluate()], with columns
#'   `protein_id`, `protein_length`, `composite_score`, and `verdict`.
#' @param label_poor Logical.  When `TRUE` (default), Poor-verdict proteins
#'   are labeled with their protein ID in panel B.
#' @param title Optional character title for the combined figure.
#'
#' @return A `patchwork` object.
#' @export
plot_batch_summary <- function(batch,
                               label_poor = TRUE,
                               title      = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots")
  rlang::check_installed("patchwork",
    reason = "to assemble plot_batch_summary panels")

  required_cols <- c("protein_id", "protein_length", "composite_score", "verdict")
  missing_cols  <- setdiff(required_cols, names(batch))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      c("!" = "{.arg batch} must be a tibble from {.fn batch_evaluate}.",
        "i" = "Missing columns: {.val {missing_cols}}."),
      class = "pepvet_error_invalid_batch"
    )
  }

  batch$composite_score <- as.numeric(batch$composite_score)
  batch$protein_length  <- as.numeric(batch$protein_length)
  batch$verdict         <- as.character(batch$verdict)

  verdict_colors <- c(
    Good     = .pepvet_pal$good,
    Moderate = .pepvet_pal$moderate,
    Poor     = .pepvet_pal$poor
  )
  verdict_levels <- c("Good", "Moderate", "Poor")
  batch$verdict  <- factor(batch$verdict, levels = verdict_levels)

  n_total    <- nrow(batch)
  n_good     <- sum(batch$verdict == "Good",     na.rm = TRUE)
  n_moderate <- sum(batch$verdict == "Moderate", na.rm = TRUE)
  n_poor     <- sum(batch$verdict == "Poor",     na.rm = TRUE)

  # ── Panel A: Score histogram ─────────────────────────────────────────────
  pa <- ggplot2::ggplot(
    batch,
    ggplot2::aes(x = composite_score, fill = verdict)
  ) +
    ggplot2::geom_histogram(
      binwidth = 0.05,
      color    = "white",
      linewidth = 0.2,
      alpha    = 0.88,
      position = "stack"
    ) +
    ggplot2::geom_vline(
      xintercept = 0.40,
      linetype   = "dashed",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = 0.70,
      linetype   = "dashed",
      color      = .pepvet_pal$good,
      linewidth  = 0.8
    ) +
    ggplot2::scale_fill_manual(
      values = verdict_colors,
      name   = "Verdict",
      guide  = ggplot2::guide_legend(
        override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(add = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(
      title    = "Score Distribution",
      subtitle = sprintf(
        "%d proteins  \u00b7  Good: %d  \u00b7  Moderate: %d  \u00b7  Poor: %d",
        n_total, n_good, n_moderate, n_poor
      ),
      x = "Composite score",
      y = "Number of proteins"
    ) +
    .pepvet_theme()

  # ── Panel B: Score vs. length scatter ────────────────────────────────────
  poor_df <- batch[!is.na(batch$verdict) & batch$verdict == "Poor", , drop = FALSE]
  if (nrow(poor_df) > 0L) {
    poor_df$display_id <- vapply(
      as.character(poor_df$protein_id), .tidy_protein_id, character(1L))
  }

  pb <- ggplot2::ggplot(
    batch,
    ggplot2::aes(x = protein_length, y = composite_score, color = verdict)
  ) +
    ggplot2::geom_hline(
      yintercept = 0.40,
      linetype   = "dashed",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.5,
      alpha      = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = 0.70,
      linetype   = "dashed",
      color      = .pepvet_pal$good,
      linewidth  = 0.5,
      alpha      = 0.7
    ) +
    ggplot2::geom_smooth(
      method    = "loess",
      formula   = y ~ x,
      se        = TRUE,
      color     = .pepvet_pal$brand_dark,
      fill      = .pepvet_pal$brand_light,
      linewidth = 0.7,
      alpha     = 0.15
    ) +
    ggplot2::geom_point(
      shape = 21,
      size  = 2.4,
      stroke = 0.3,
      ggplot2::aes(fill = verdict),
      color = "white",
      alpha = 0.85
    ) +
    ggplot2::scale_color_manual(values = verdict_colors, guide = "none") +
    ggplot2::scale_fill_manual(
      values = verdict_colors,
      name   = "Verdict",
      guide  = "none"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title    = "Composite Score vs. Protein Length",
      subtitle = "LOESS trend  \u00b7  Poor-verdict proteins labeled",
      x = "Protein length (aa)",
      y = "Composite score"
    ) +
    .pepvet_theme()

  # Add Poor protein labels
  if (label_poor && nrow(poor_df) > 0L) {
    pb <- pb + ggplot2::geom_text(
      data = poor_df,
      ggplot2::aes(
        x     = protein_length,
        y     = composite_score,
        label = display_id
      ),
      size          = 2.4,
      hjust         = -0.12,
      color         = .pepvet_pal$poor,
      check_overlap = TRUE
    )
  }

  auto_title <- title %||% "Batch Digest Summary"

  (pa | pb) +
    patchwork::plot_layout(widths = c(1, 1.4)) +
    patchwork::plot_annotation(
      title     = auto_title,
      tag_levels = "A",
      theme      = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face   = "bold",
          size   = 14,
          color  = .pepvet_pal$brand_dark,
          margin = ggplot2::margin(b = 8)
        )
      )
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face  = "bold",
        size  = 13,
        color = .pepvet_pal$brand
      )
    )
}


# ── plot_proteome_heatmap ─────────────────────────────────────────────────────

#' Proteome Component Score Heatmap
#'
#' `plot_proteome_heatmap()` draws a clustered heatmap of all component scores
#' across all proteins in a [batch_evaluate()] result.  Rows (proteins) are
#' hierarchically clustered by score profile similarity so proteins that fail
#' for the same reason appear together.  A verdict color sidebar is drawn on
#' the left.
#'
#' Requires the `pheatmap` package (suggested dependency).
#'
#' @param batch A tibble from [batch_evaluate()], with columns `protein_id`,
#'   `verdict`, and the component score columns.
#' @param components Character vector of component columns to include.
#'   Defaults to `c("S_length","S_coverage","S_count","S_hydro","S_charge")`.
#' @param cluster_rows Logical.  Cluster rows (proteins) by score similarity.
#'   Default `TRUE`.
#' @param cluster_cols Logical.  Cluster columns (components).  Default `FALSE`
#'   (keeps the natural component order).
#' @param title Optional character title passed to `pheatmap::pheatmap()`.
#'
#' @return A `pheatmap` object (not a ggplot).
#' @export
plot_proteome_heatmap <- function(
    batch,
    components   = c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge"),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    title        = NULL) {

  rlang::check_installed("pheatmap",
    reason = "to draw clustered heatmaps in plot_proteome_heatmap()")

  required_cols <- c("protein_id", "verdict")
  missing_cols  <- setdiff(required_cols, names(batch))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      c("!" = "{.arg batch} must be a tibble from {.fn batch_evaluate}.",
        "i" = "Missing columns: {.val {missing_cols}}."),
      class = "pepvet_error_invalid_batch"
    )
  }

  valid_comp <- components[components %in% names(batch)]
  if (length(valid_comp) == 0L) {
    cli::cli_abort(
      c("!" = "None of the requested {.arg components} are present in {.arg batch}.",
        "i" = "Available columns: {.val {names(batch)}}."),
      class = "pepvet_error_invalid_batch"
    )
  }

  # ── Build matrix ─────────────────────────────────────────────────────────
  mat <- as.matrix(batch[, valid_comp, drop = FALSE])
  # Tidy row names
  row.names(mat) <- vapply(
    as.character(batch$protein_id), .tidy_protein_id, character(1L))
  # Clean column names
  col_labels <- c(
    S_length   = "Length",
    S_coverage = "Coverage",
    S_count    = "Count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge"
  )
  colnames(mat) <- ifelse(
    colnames(mat) %in% names(col_labels),
    col_labels[colnames(mat)],
    colnames(mat)
  )

  # ── Sidebar annotation ────────────────────────────────────────────────────
  verdict_char <- as.character(batch$verdict)
  ann_row <- data.frame(
    Verdict = factor(verdict_char, levels = c("Good", "Moderate", "Poor")),
    row.names = row.names(mat)
  )
  ann_colors <- list(
    Verdict = c(
      Good     = .pepvet_pal$good,
      Moderate = .pepvet_pal$moderate,
      Poor     = .pepvet_pal$poor
    )
  )

  auto_title <- title %||% "Proteome Component Score Heatmap"

  pheatmap::pheatmap(
    mat,
    color            = grDevices::colorRampPalette(
      c(.pepvet_pal$poor, "#FFFAEC", .pepvet_pal$good)
    )(100),
    breaks           = seq(0, 1, length.out = 101),
    cluster_rows     = cluster_rows,
    cluster_cols     = cluster_cols,
    annotation_row   = ann_row,
    annotation_colors = ann_colors,
    fontsize         = 9,
    fontsize_row     = 7,
    fontsize_col     = 9,
    border_color     = "white",
    cellheight       = max(8, floor(400 / nrow(mat))),
    cellwidth        = 50,
    main             = auto_title
  )
}


# ── plot_component_scatter ────────────────────────────────────────────────────

#' Proteome Component Score 2D Scatter
#'
#' `plot_component_scatter()` plots any two component scores against each
#' other for all proteins in a [batch_evaluate()] tibble.  Points are colored
#' by verdict.  Poor-verdict proteins are labeled; the top-N proteins by
#' composite score can also be labeled via `label_top_n`.  Rug marks on both
#' axes show the marginal distributions.
#'
#' @param batch A tibble from [batch_evaluate()].
#' @param x_component Character.  Column for the x-axis.  Default
#'   `"S_hydro"`.
#' @param y_component Character.  Column for the y-axis.  Default
#'   `"S_charge"`.
#' @param size_by_length Logical.  When `TRUE` point size is proportional to
#'   `protein_length`.  Default `FALSE`.
#' @param label_top_n Integer.  Additionally label the top-N proteins by
#'   composite score.  Default `0` (no extra labels).
#' @param title Optional character title.
#'
#' @return A `ggplot` object.
#' @export
plot_component_scatter <- function(
    batch,
    x_component    = "S_hydro",
    y_component    = "S_charge",
    size_by_length = FALSE,
    label_top_n    = 0L,
    title          = NULL) {

  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots")

  valid_components <- c("composite_score", "S_length", "S_coverage",
                        "S_count", "S_hydro", "S_charge")
  for (comp in c(x_component, y_component)) {
    if (!comp %in% valid_components) {
      cli::cli_abort(
        c("!" = "{.arg {comp}} must be one of {.val {valid_components}}.",
          "i" = "Received: {.val {comp}}."),
        class = "pepvet_error_invalid_component"
      )
    }
  }

  required_cols <- c("protein_id", "verdict", x_component, y_component)
  missing_cols  <- setdiff(required_cols, names(batch))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      c("!" = "{.arg batch} must be a tibble from {.fn batch_evaluate}.",
        "i" = "Missing columns: {.val {missing_cols}}."),
      class = "pepvet_error_invalid_batch"
    )
  }

  batch$verdict <- factor(
    as.character(batch$verdict), levels = c("Good", "Moderate", "Poor"))
  batch$x_val <- as.numeric(batch[[x_component]])
  batch$y_val <- as.numeric(batch[[y_component]])

  verdict_colors <- c(
    Good     = .pepvet_pal$good,
    Moderate = .pepvet_pal$moderate,
    Poor     = .pepvet_pal$poor
  )

  # ── Label set: Poor + optional top-N ────────────────────────────────────
  poor_mask <- !is.na(batch$verdict) & batch$verdict == "Poor"
  label_mask <- poor_mask

  if (label_top_n > 0L && "composite_score" %in% names(batch)) {
    cs_rank <- rank(-as.numeric(batch$composite_score), ties.method = "first")
    label_mask <- label_mask | cs_rank <= label_top_n
  }
  label_df <- batch[label_mask, , drop = FALSE]
  if (nrow(label_df) > 0L) {
    label_df$display_id <- vapply(
      as.character(label_df$protein_id), .tidy_protein_id, character(1L))
  }

  # ── Human-readable axis labels ─────────────────────────────────────────
  comp_labels <- c(
    composite_score = "Composite score",
    S_length        = "Length score",
    S_coverage      = "Coverage score",
    S_count         = "Count score",
    S_hydro         = "Hydrophobicity score",
    S_charge        = "Charge score"
  )
  x_label <- comp_labels[x_component]
  y_label <- comp_labels[y_component]

  auto_title <- title %||% paste0(
    "Proteome Component Scatter  \u00b7  ",
    x_label, " vs. ", y_label
  )

  p <- ggplot2::ggplot(
    batch,
    ggplot2::aes(x = x_val, y = y_val, color = verdict, fill = verdict)
  ) +
    # Threshold reference lines
    ggplot2::geom_vline(
      xintercept = c(0.40, 0.70),
      linetype   = c("dashed", "dashed"),
      color      = c(.pepvet_pal$moderate, .pepvet_pal$good),
      linewidth  = 0.5,
      alpha      = 0.65
    ) +
    ggplot2::geom_hline(
      yintercept = c(0.40, 0.70),
      linetype   = c("dashed", "dashed"),
      color      = c(.pepvet_pal$moderate, .pepvet_pal$good),
      linewidth  = 0.5,
      alpha      = 0.65
    ) +
    ggplot2::geom_rug(
      alpha    = 0.35,
      linewidth = 0.4,
      sides    = "bl"
    )

  # Point layer: optionally size by protein length
  if (size_by_length && "protein_length" %in% names(batch)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(size = protein_length),
        shape  = 21,
        stroke = 0.3,
        color  = "white",
        alpha  = 0.82
      ) +
      ggplot2::scale_size_continuous(
        range  = c(2, 6),
        name   = "Protein length (aa)",
        guide  = ggplot2::guide_legend(
          override.aes = list(color = NA, fill = .pepvet_pal$brand),
          title.position = "top"
        )
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        shape  = 21,
        size   = 2.8,
        stroke = 0.3,
        color  = "white",
        alpha  = 0.82
      )
  }

  p <- p +
    ggplot2::scale_color_manual(values = verdict_colors, name = "Verdict") +
    ggplot2::scale_fill_manual( values = verdict_colors, name = "Verdict") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = "Dashed lines at 0.40 / 0.70  \u00b7  Poor proteins labeled  \u00b7  Rug = marginal density",
      x = x_label,
      y = y_label
    ) +
    .pepvet_theme()

  if (nrow(label_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x_val, y = y_val, label = display_id),
      size          = 2.5,
      hjust         = -0.1,
      color         = "#333333",
      check_overlap = TRUE
    )
  }

  p
}
