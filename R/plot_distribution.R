# ── pepVet Distribution & landscape plots ──────────────────────────────────
# ── plot_length_distribution ──────────────────────────────────────────────────

#' Standalone Peptide Length Distribution
#'
#' `plot_length_distribution()` draws a histogram of peptide lengths
#' colour-coded by validity class (Valid / Too short / Too long), with the
#' valid range shaded in the package's green, per-category percentage
#' annotations, and an optional density-curve overlay.
#'
#' @param result A named list returned by [evaluate_digest()], **or** a
#'   data.frame / tibble with at least a `length` column (e.g. the
#'   `$peptides` slot of such a result).
#' @param length_range Integer vector of length 2 giving the valid length
#'   window `c(lo, hi)`.  Defaults to `c(7L, 25L)`.  Ignored (and read from
#'   `result$params`) when a full `evaluate_digest()` result is supplied.
#' @param show_density Logical.  When `TRUE` (default) a scaled kernel-density
#'   curve is overlaid on the histogram.
#' @param title Optional character string for the plot title.  Auto-generated
#'   when `NULL` (default).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   p   <- plot_length_distribution(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [plot_digest_profile()],
#'   [plot_gravy_landscape()]
#' @export
plot_length_distribution <- function(
    result,
    length_range  = c(7L, 25L),
    show_density  = TRUE,
    title         = NULL
) {
  rlang::check_installed("ggplot2", reason = "to use plot_length_distribution()")

  # ── Multi-input mode: named list of evaluate_digest() results ────────────
  if (.is_named_results_list(result)) {
    return(.plot_length_distribution_multi(result, length_range = length_range,
                                           show_density = show_density,
                                           title = title))
  }

  # ── Accept evaluate_digest() list or a bare peptide data.frame ───────────
  if (is.list(result) && !is.data.frame(result) &&
      all(c("peptides", "params") %in% names(result))) {
    peps         <- result$peptides
    length_range <- result$params$length_range %||% length_range
  } else if (is.data.frame(result)) {
    peps <- result
  } else {
    .abort(
      c(
        "{.arg result} must be an {.fn evaluate_digest} list or a peptide data.frame.",
        "x" = "Got {.cls {class(result)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }

  if (!"length" %in% names(peps)) {
    .abort(
      "{.arg result} must contain a {.field length} column.",
      class = "pepvet_error_invalid_digest_result"
    )
  }

  length_lo <- as.integer(length_range[[1L]])
  length_hi <- as.integer(length_range[[2L]])

  # ── Classify peptides ─────────────────────────────────────────────────────
  peps$length_class <- factor(
    ifelse(peps$length < length_lo, "Too short",
      ifelse(peps$length > length_hi, "Too long", "Valid")),
    levels = c("Valid", "Too short", "Too long")
  )

  n_total    <- nrow(peps)
  tbl        <- table(peps$length_class)
  pct_valid  <- round(100 * sum(peps$length_class == "Valid")    / n_total, 1)
  pct_short  <- round(100 * sum(peps$length_class == "Too short") / n_total, 1)
  pct_long   <- round(100 * sum(peps$length_class == "Too long")  / n_total, 1)

  class_colors <- c(
    "Valid"     = .pepvet_pal$valid,
    "Too short" = .pepvet_pal$too_short,
    "Too long"  = .pepvet_pal$too_long
  )

  x_max <- max(peps$length) + 1L
  x_min <- max(0L, min(peps$length) - 1L)

  # ── Category annotation positions: each label sits at the centre of its
  #    x-range, just below the top of the panel ─────────────────────────────
  cat_labels <- data.frame(
    x     = c((x_min + length_lo - 1) / 2,
              (length_lo + length_hi) / 2,
              (length_hi + x_max + 1) / 2),
    label = c(
      sprintf("Too short\n< %d aa\n(%.0f%%)", length_lo, pct_short),
      sprintf("Valid\n[%d\u2013%d aa]\n(%.0f%%)", length_lo, length_hi, pct_valid),
      sprintf("Too long\n> %d aa\n(%.0f%%)", length_hi, pct_long)
    ),
    color = unname(class_colors[c("Too short", "Valid", "Too long")]),
    stringsAsFactors = FALSE
  )
  # Drop labels for empty categories
  cat_labels <- cat_labels[
    c(pct_short > 0, TRUE, pct_long > 0), , drop = FALSE
  ]

  # ── Auto title ────────────────────────────────────────────────────────────
  auto_title <- if (!is.null(title)) {
    title
  } else if (is.list(result) && !is.data.frame(result) &&
             "params" %in% names(result)) {
    pid    <- result$params$protein_ids[[1L]]
    enzyme <- result$params$enzyme
    paste0(.tidy_protein_id(pid), "  \u00b7  ", enzyme,
           "  \u00b7  Length distribution")
  } else {
    "Peptide length distribution"
  }

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(peps, ggplot2::aes(x = .data$length,
                                           fill = .data$length_class)) +
    # Valid-range background shading
    ggplot2::annotate(
      "rect",
      xmin = length_lo - 0.5, xmax = length_hi + 0.5,
      ymin = -Inf,             ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.55
    ) +
    # Boundary lines at valid range edges
    ggplot2::geom_vline(
      xintercept = length_lo - 0.5,
      color = .pepvet_pal$good, linewidth = 0.6, linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = length_hi + 0.5,
      color = .pepvet_pal$poor, linewidth = 0.6, linetype = "dashed"
    ) +
    # Histogram bars
    ggplot2::geom_histogram(
      binwidth  = 1L,
      color     = "white",
      linewidth = 0.15,
      alpha = .get_param("scatter_alpha")
    ) +
    # Optional density overlay: computed over ALL peptides (one curve,
    # independent of fill grouping) to avoid group-size warnings
    {
      if (show_density) {
        ggplot2::stat_density(
          ggplot2::aes(x = .data$length,
                       y = ggplot2::after_stat(count)),
          data        = peps,
          geom        = "line",
          color       = .pepvet_pal$brand_dark,
          linewidth   = 0.8,
          linetype    = "solid",
          adjust      = 1.2,
          inherit.aes = FALSE
        )
      }
    } +
    # Per-category annotation labels
    ggplot2::geom_text(
      data = cat_labels,
      ggplot2::aes(x = .data$x, label = .data$label),
      y        = Inf,
      vjust    = 1.2,
      size     = 2.8,
      fontface = "bold",
      color    = cat_labels$color,
      inherit.aes = FALSE
    ) +
    # Scales
    ggplot2::scale_fill_manual(
      values = class_colors,
      name   = NULL,
      guide  = ggplot2::guide_legend(
        override.aes = list(alpha = 1, color = NA)
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, x_max + 4L, by = 5L),
      expand = ggplot2::expansion(add = c(0.5, 1))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max + 1), clip = "off") +
    ggplot2::labs(
      title    = auto_title,
      subtitle = sprintf(
        "%d peptides total  \u00b7  %d valid (%.0f%%)  \u00b7  range [%d\u2013%d aa]",
        n_total,
        as.integer(tbl[["Valid"]]),
        pct_valid,
        length_lo, length_hi
      ),
      x = "Peptide length (aa)",
      y = "Count"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark
      )
    )

  p
}

# Private helper: detect a named list of evaluate_digest() results
.is_named_results_list <- function(x) {
  is.list(x) && !is.data.frame(x) &&
    length(x) >= 1L &&
    !all(c("peptides", "params") %in% names(x)) &&
    all(vapply(x, function(r)
      is.list(r) && !is.data.frame(r) &&
        all(c("peptides", "params") %in% names(r)),
      logical(1L)))
}

# Private helper: extract auto-label for a single evaluate_digest result
.result_label <- function(r) {
  paste0(.tidy_protein_id(r$params$protein_ids[[1L]]),
         " / ", r$params$enzyme)
}

# Private: multi-input length distribution (faceted)
.plot_length_distribution_multi <- function(results, length_range, show_density, title) {
  rlang::check_installed("ggplot2", reason = "to use plot_length_distribution()")

  labels <- if (!is.null(names(results))) names(results) else
    vapply(results, .result_label, character(1L))

  peps_list <- lapply(seq_along(results), function(i) {
    r   <- results[[i]]
    lr  <- r$params$length_range %||% length_range
    df  <- r$peptides
    df$length_range_lo <- as.integer(lr[[1L]])
    df$length_range_hi <- as.integer(lr[[2L]])
    df$.label <- factor(labels[[i]], levels = labels)
    df
  })
  peps <- do.call(rbind, peps_list)

  # Use the first result's range as global default for shading
  g_lo <- peps_list[[1L]]$length_range_lo[[1L]]
  g_hi <- peps_list[[1L]]$length_range_hi[[1L]]

  peps$length_class <- factor(
    ifelse(peps$length < peps$length_range_lo, "Too short",
      ifelse(peps$length > peps$length_range_hi, "Too long", "Valid")),
    levels = c("Valid", "Too short", "Too long")
  )
  class_colors <- c("Valid" = .pepvet_pal$valid,
                    "Too short" = .pepvet_pal$too_short,
                    "Too long"  = .pepvet_pal$too_long)

  x_lo <- max(0L, min(peps$length) - 1L)
  x_hi <- max(peps$length) + 1L

  auto_title <- title %||% "Peptide length distribution: comparison"

  ggplot2::ggplot(peps, ggplot2::aes(x = .data$length, fill = .data$length_class)) +
    ggplot2::annotate("rect",
      xmin = g_lo - 0.5, xmax = g_hi + 0.5,
      ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.45) +
    ggplot2::geom_vline(xintercept = g_lo - 0.5,
      color = .pepvet_pal$good, linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = g_hi + 0.5,
      color = .pepvet_pal$poor, linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_histogram(binwidth = 1L, color = "white",
      linewidth = 0.12, alpha = .get_param("scatter_alpha")) +
    {
      if (show_density) {
        ggplot2::stat_density(
          ggplot2::aes(x = .data$length, y = ggplot2::after_stat(count)),
          data = peps, geom = "line",
          color = .pepvet_pal$brand_dark, linewidth = 0.7,
          adjust = 1.2, inherit.aes = FALSE
        )
      }
    } +
    ggplot2::facet_wrap(ggplot2::vars(.data$.label), scales = "free_y") +
    ggplot2::scale_fill_manual(values = class_colors, name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(alpha = 1, color = NA))) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, x_hi + 4L, by = 5L),
      expand = ggplot2::expansion(add = c(0.5, 1))) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi + 1L), clip = "off") +
    ggplot2::labs(title = auto_title,
                  x = "Peptide length (aa)", y = "Count") +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(
        size = 9, face = "bold", color = .pepvet_pal$brand_dark),
      plot.title = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark)
    )
}


# ── plot_gravy_landscape ──────────────────────────────────────────────────────

#' GRAVY Landscape: 2D Scatter of Peptide Length vs. Hydrophobicity
#'
#' `plot_gravy_landscape()` plots each peptide as a point in the
#' length \eqn{\times} GRAVY physicochemical space.  The LC-friendly valid
#' region is shaded, points are colour-coded by validity class, and marginal
#' density panels show the 1D distributions on each axis.  Outlier peptides
#' (outside the valid region) are labelled with their sequences when there are
#' fewer than `label_outliers_n` of them.
#'
#' @param result A named list returned by [evaluate_digest()], **or** a
#'   data.frame / tibble with at least `length` and `gravy` columns.  When a
#'   bare data.frame lacks a `gravy` column but has a `peptide` column, GRAVY
#'   scores are computed automatically.
#' @param length_range Integer vector of length 2.  Valid length window.
#'   Defaults to `c(7L, 25L)`.  Read from `result$params` when a full
#'   [evaluate_digest()] result is supplied.
#' @param gravy_range Numeric vector of length 2.  LC-friendly GRAVY window.
#'   Defaults to `c(-1.0, 0.6)`.  Read from `result$params` when available.
#' @param label_outliers_n Integer.  Maximum number of outlier points to label
#'   with their peptide sequences.  Labels are suppressed when the count
#'   exceeds this threshold.  Defaults to `15L`.
#' @param title Optional character string for the plot title.
#'
#' @return A `patchwork` composite `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   p   <- plot_gravy_landscape(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [plot_length_distribution()],
#'   [plot_digest_profile()]
#' @export
plot_gravy_landscape <- function(
    result,
    length_range     = c(7L, 25L),
    gravy_range      = c(-1.0, 0.6),
    label_outliers_n = 15L,
    title            = NULL
) {
  rlang::check_installed("ggplot2",   reason = "to use plot_gravy_landscape()")
  rlang::check_installed("patchwork", reason = "to assemble panels in plot_gravy_landscape()")

  # ── Multi-input mode: named list of evaluate_digest() results ────────────
  if (.is_named_results_list(result)) {
    return(.plot_gravy_landscape_multi(result, length_range = length_range,
                                       gravy_range = gravy_range, title = title))
  }

  # ── Parse input ───────────────────────────────────────────────────────────
  if (is.list(result) && !is.data.frame(result) &&
      all(c("peptides", "params") %in% names(result))) {
    peps         <- result$peptides
    length_range <- result$params$length_range %||% length_range
    gravy_range  <- result$params$gravy_range  %||% gravy_range
  } else if (is.data.frame(result)) {
    peps <- result
  } else {
    .abort(
      c(
        "{.arg result} must be an {.fn evaluate_digest} list or a peptide data.frame.",
        "x" = "Got {.cls {class(result)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }

  if (!"length" %in% names(peps)) {
    .abort(
      "{.arg result} must contain a {.field length} column.",
      class = "pepvet_error_invalid_digest_result"
    )
  }

  # Compute GRAVY if not already present
  if (!"gravy" %in% names(peps)) {
    if (!"peptide" %in% names(peps)) {
      .abort(
        "Cannot compute GRAVY: {.arg result} needs a {.field gravy} or {.field peptide} column.",
        class = "pepvet_error_invalid_digest_result"
      )
    }
    peps$gravy <- .calculate_gravy_vec(peps$peptide)
  }

  length_lo <- as.integer(length_range[[1L]])
  length_hi <- as.integer(length_range[[2L]])
  gravy_lo  <- gravy_range[[1L]]
  gravy_hi  <- gravy_range[[2L]]

  # ── Classify peptides ─────────────────────────────────────────────────────
  peps$valid_length <- peps$length >= length_lo & peps$length <= length_hi
  peps$valid_gravy  <- peps$gravy  >= gravy_lo  & peps$gravy  <= gravy_hi
  peps$class <- factor(
    ifelse(
       peps$valid_length &  peps$valid_gravy,  "Valid",
      ifelse(
        !peps$valid_length &  peps$valid_gravy, "Invalid length",
        ifelse(
           peps$valid_length & !peps$valid_gravy, "Invalid GRAVY",
          "Invalid both"
        )
      )
    ),
    levels = c("Valid", "Invalid length", "Invalid GRAVY", "Invalid both")
  )

  class_colors <- c(
    "Valid"          = .pepvet_pal$valid,
    "Invalid length" = .pepvet_pal$too_short,
    "Invalid GRAVY"  = .pepvet_pal$moderate,
    "Invalid both"   = .pepvet_pal$poor
  )

  n_total   <- nrow(peps)
  n_valid   <- sum(peps$class == "Valid")
  pct_valid <- round(100 * n_valid / n_total, 1)

  # ── Outlier labels ────────────────────────────────────────────────────────
  outliers  <- peps[peps$class != "Valid", , drop = FALSE]
  do_label  <- "peptide" %in% names(peps) &&
    nrow(outliers) > 0L && nrow(outliers) <= as.integer(label_outliers_n)

  # ── Axis limits with padding ──────────────────────────────────────────────
  x_pad  <- 1.5
  y_pad  <- 0.15
  x_lims <- c(max(0, min(peps$length) - x_pad), max(peps$length) + x_pad)
  y_lims <- c(min(peps$gravy, na.rm = TRUE) - y_pad,
              max(peps$gravy, na.rm = TRUE) + y_pad)

  # ── Auto title ────────────────────────────────────────────────────────────
  auto_title <- if (!is.null(title)) {
    title
  } else if (is.list(result) && !is.data.frame(result) &&
             "params" %in% names(result)) {
    pid    <- result$params$protein_ids[[1L]]
    enzyme <- result$params$enzyme
    paste0(.tidy_protein_id(pid), "  \u00b7  ", enzyme,
           "  \u00b7  GRAVY landscape")
  } else {
    "GRAVY landscape"
  }

  # ── Central scatter ───────────────────────────────────────────────────────
  p_scatter <- ggplot2::ggplot(
    peps, ggplot2::aes(x = .data$length, y = .data$gravy,
                       color = .data$class, fill = .data$class)
  ) +
    # Valid region background
    ggplot2::annotate(
      "rect",
      xmin = length_lo - 0.5, xmax = length_hi + 0.5,
      ymin = gravy_lo,         ymax = gravy_hi,
      fill = .pepvet_pal$shade, alpha = 0.55, color = NA
    ) +
    # Valid region dashed border
    ggplot2::annotate(
      "rect",
      xmin = length_lo - 0.5, xmax = length_hi + 0.5,
      ymin = gravy_lo,         ymax = gravy_hi,
      fill = NA, color = .pepvet_pal$good,
      linewidth = 0.45, linetype = "dashed"
    ) +
    # Reference lines at valid boundaries
    ggplot2::geom_hline(
      yintercept = c(gravy_lo, gravy_hi),
      color = .pepvet_pal$separator, linewidth = 0.3, linetype = "dotted"
    ) +
    ggplot2::geom_vline(
      xintercept = c(length_lo - 0.5, length_hi + 0.5),
      color = .pepvet_pal$separator, linewidth = 0.3, linetype = "dotted"
    ) +
    # Points (jittered horizontally to avoid over-plotting at integer lengths)
    ggplot2::geom_jitter(
      shape  = 21,
      size   = 2.2,
      stroke = 0.35,
      width  = 0.25,
      height = 0,
      alpha  = 0.80
    ) +
    # Outlier sequence labels when count <= label_outliers_n
    {
      if (do_label) {
        ggplot2::geom_text(
          data    = outliers,
          ggplot2::aes(label = .data$peptide),
          size          = 1.9,
          fontface      = "italic",
          nudge_y       = 0.06,
          check_overlap = TRUE,
          inherit.aes   = TRUE
        )
      }
    } +
    ggplot2::scale_color_manual(values = class_colors, name = NULL,
      guide = "none") +
    ggplot2::scale_fill_manual(values = class_colors, name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(shape = 21, size = 3, alpha = 1,
                            stroke = 0.5, color = "white")
      )) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, ceiling(x_lims[[2L]] / 5) * 5L, by = 5L)
    ) +
    ggplot2::coord_cartesian(xlim = x_lims, ylim = y_lims) +
    ggplot2::labs(
      x        = "Peptide length (aa)",
      y        = "GRAVY score",
      subtitle = sprintf(
        "%d / %d peptides fully valid (%.0f%%)  \u00b7  valid region [%d\u2013%d aa, %.1f\u2013%.1f GRAVY]",
        n_valid, n_total, pct_valid, length_lo, length_hi, gravy_lo, gravy_hi
      )
    ) +
    .pepvet_theme() +
    ggplot2::theme(legend.position = "bottom")

  # ── Top marginal: length density by class ────────────────────────────────
  p_top <- ggplot2::ggplot(
    peps, ggplot2::aes(x = .data$length, fill = .data$class)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = length_lo - 0.5, xmax = length_hi + 0.5,
      ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.30
    ) +
    ggplot2::geom_density(alpha = 0.50, linewidth = 0.25, adjust = 1.2,
                          color = NA) +
    ggplot2::scale_fill_manual(values = class_colors, guide = "none") +
    ggplot2::coord_cartesian(xlim = x_lims) +
    ggplot2::labs(x = NULL, y = NULL) +
    .pepvet_theme() +
    ggplot2::theme(
      panel.border      = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      axis.text         = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(0, 0, 2, 0)
    )

  # ── Right marginal: GRAVY density by class (coord_flip) ───────────────────
  p_right <- ggplot2::ggplot(
    peps, ggplot2::aes(x = .data$gravy, fill = .data$class)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = gravy_lo, xmax = gravy_hi,
      ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.30
    ) +
    ggplot2::geom_density(alpha = 0.50, linewidth = 0.25, adjust = 1.2,
                          color = NA) +
    ggplot2::scale_fill_manual(values = class_colors, guide = "none") +
    ggplot2::coord_flip(xlim = y_lims) +
    ggplot2::labs(x = NULL, y = NULL) +
    .pepvet_theme() +
    ggplot2::theme(
      panel.border      = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      axis.text         = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(0, 0, 0, 2)
    )

  # ── Assemble: top-marginal / (scatter + right-marginal) ───────────────────
  (p_top | patchwork::plot_spacer()) /
    (p_scatter | p_right) +
    patchwork::plot_layout(heights = c(1, 4), widths = c(4, 1)) +
    patchwork::plot_annotation(
      title = auto_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size   = 13, face = "bold",
          color  = .pepvet_pal$brand_dark,
          margin = ggplot2::margin(b = 4)
        )
      )
    )
}

# Private: multi-input GRAVY landscape (faceted, no marginals)
.plot_gravy_landscape_multi <- function(results, length_range, gravy_range, title) {
  rlang::check_installed("ggplot2", reason = "to use plot_gravy_landscape()")

  labels <- if (!is.null(names(results))) names(results) else
    vapply(results, .result_label, character(1L))

  peps_list <- lapply(seq_along(results), function(i) {
    r   <- results[[i]]
    lr  <- r$params$length_range %||% length_range
    gr  <- r$params$gravy_range  %||% gravy_range
    df  <- r$peptides
    if (!"gravy" %in% names(df) && "peptide" %in% names(df))
      df$gravy <- .calculate_gravy_vec(df$peptide)
    df$length_lo <- as.integer(lr[[1L]])
    df$length_hi <- as.integer(lr[[2L]])
    df$gravy_lo  <- gr[[1L]]
    df$gravy_hi  <- gr[[2L]]
    df$.label    <- factor(labels[[i]], levels = labels)
    df
  })
  peps <- do.call(rbind, peps_list)

  g_len_lo <- peps_list[[1L]]$length_lo[[1L]]
  g_len_hi <- peps_list[[1L]]$length_hi[[1L]]
  g_grav_lo <- peps_list[[1L]]$gravy_lo[[1L]]
  g_grav_hi <- peps_list[[1L]]$gravy_hi[[1L]]

  peps$valid_length <- peps$length >= peps$length_lo & peps$length <= peps$length_hi
  peps$valid_gravy  <- peps$gravy  >= peps$gravy_lo  & peps$gravy  <= peps$gravy_hi
  peps$class <- factor(
    ifelse( peps$valid_length &  peps$valid_gravy, "Valid",
      ifelse(!peps$valid_length &  peps$valid_gravy, "Invalid length",
        ifelse( peps$valid_length & !peps$valid_gravy, "Invalid GRAVY",
          "Invalid both"))),
    levels = c("Valid", "Invalid length", "Invalid GRAVY", "Invalid both")
  )
  class_colors <- c("Valid" = .pepvet_pal$valid,
                    "Invalid length" = .pepvet_pal$too_short,
                    "Invalid GRAVY"  = .pepvet_pal$moderate,
                    "Invalid both"   = .pepvet_pal$poor)

  auto_title <- title %||% "GRAVY landscape: comparison"

  ggplot2::ggplot(
    peps, ggplot2::aes(x = .data$length, y = .data$gravy,
                       fill = .data$class, color = .data$class)
  ) +
    ggplot2::annotate("rect",
      xmin = g_len_lo - 0.5, xmax = g_len_hi + 0.5,
      ymin = g_grav_lo, ymax = g_grav_hi,
      fill = .pepvet_pal$shade, alpha = 0.50, color = NA) +
    ggplot2::annotate("rect",
      xmin = g_len_lo - 0.5, xmax = g_len_hi + 0.5,
      ymin = g_grav_lo, ymax = g_grav_hi,
      fill = NA, color = .pepvet_pal$good,
      linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_jitter(shape = 21, size = 1.8, stroke = 0.3,
                         width = 0.2, height = 0, alpha = 0.75) +
    ggplot2::facet_wrap(ggplot2::vars(.data$.label)) +
    ggplot2::scale_color_manual(values = class_colors, name = NULL,
                                guide = "none") +
    ggplot2::scale_fill_manual(values = class_colors, name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(shape = 21, size = 3, alpha = 1,
                            stroke = 0.5, color = "white"))) +
    ggplot2::labs(title = auto_title,
                  x = "Peptide length (aa)", y = "GRAVY score") +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(
        size = 9, face = "bold", color = .pepvet_pal$brand_dark),
      plot.title = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark)
    )
}


# ── plot_pI_distribution ──────────────────────────────────────────────────────

#' pI Distribution: Histogram of Peptide Isoelectric Points
#'
#' `plot_pI_distribution()` draws a histogram of peptide isoelectric points
#' coloured by SCX fraction bin (e.g., pH 3–4, 4–5, …) to preview
#' fractionation outcomes.  Vertical boundary lines and per-fraction count
#' annotations are optionally overlaid.
#'
#' @param result Accepted inputs:
#'   * A named list returned by [evaluate_digest()].  pI values are computed
#'     automatically from the valid-peptide sequences.
#'   * A tibble returned by [score_peptides()] with `include_pI = TRUE`
#'     (contains a `pI` list column).
#'   * A plain data.frame / tibble with a numeric `pI` column.
#'   * A bare numeric vector of pI values.
#' @param fraction_breaks Numeric vector of pH boundary values defining the
#'   fraction bins.  Defaults to `c(3, 4, 5, 6, 7, 8, 9, 10)`, which produces
#'   eight bins: `<3`, `3–4`, …, `9–10`, `>10`.
#' @param show_fraction_lines Logical.  When `TRUE` (default) vertical dashed
#'   lines are drawn at each interior fraction boundary.
#' @param title Optional character string for the plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   p   <- plot_pI_distribution(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [score_peptides()], [plot_length_distribution()],
#'   [plot_gravy_landscape()]
#' @export
plot_pI_distribution <- function(
    result,
    fraction_breaks     = c(3, 4, 5, 6, 7, 8, 9, 10),
    show_fraction_lines = TRUE,
    title               = NULL
) {
  rlang::check_installed("ggplot2", reason = "to use plot_pI_distribution()")

  # ── Multi-input mode: named list of evaluate_digest() results ────────────
  if (.is_named_results_list(result)) {
    return(.plot_pI_distribution_multi(result,
                                       fraction_breaks = fraction_breaks,
                                       show_fraction_lines = show_fraction_lines,
                                       title = title))
  }

  # ── Extract pI values ─────────────────────────────────────────────────────
  pI_vals <- if (is.numeric(result)) {
    result[!is.na(result)]

  } else if (is.list(result) && !is.data.frame(result) &&
             all(c("peptides", "params") %in% names(result))) {
    # evaluate_digest() list: compute pI for valid peptides
    peps   <- result$peptides
    lr     <- result$params$length_range %||% c(7L, 25L)
    valid  <- peps[peps$length >= lr[[1L]] & peps$length <= lr[[2L]], , drop = FALSE]
    if (nrow(valid) == 0L || !"peptide" %in% names(valid)) {
      .abort(
        "No valid peptides found in {.arg result} to compute pI values.",
        class = "pepvet_error_invalid_digest_result"
      )
    }
    as.numeric(calculate_pI(valid$peptide))

  } else if (is.data.frame(result)) {
    if ("pI" %in% names(result) && is.list(result$pI)) {
      # score_peptides(include_pI = TRUE): unlist the list column
      vals <- unlist(result$pI, use.names = FALSE)
      as.numeric(vals[!is.na(vals)])
    } else if ("pI" %in% names(result)) {
      as.numeric(result$pI[!is.na(result$pI)])
    } else {
      .abort(
        "{.arg result} data.frame must contain a {.field pI} column.",
        class = "pepvet_error_invalid_digest_result"
      )
    }

  } else {
    .abort(
      c(
        "{.arg result} must be an {.fn evaluate_digest} list, a {.fn score_peptides} tibble, a data.frame with a {.field pI} column, or a numeric vector.",
        "x" = "Got {.cls {class(result)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }

  if (length(pI_vals) == 0L) {
    .abort("No pI values to plot.", class = "pepvet_error_invalid_digest_result")
  }

  # ── Fraction bins ─────────────────────────────────────────────────────────
  breaks   <- sort(unique(as.numeric(fraction_breaks)))
  lo_break <- breaks[[1L]]
  hi_break <- breaks[[length(breaks)]]

  # Build label vector: "<lo", "lo–b2", "b2–b3", …, ">hi"
  bin_labels <- character(length(breaks) + 1L)
  bin_labels[[1L]] <- sprintf("< %.0f", lo_break)
  for (i in seq_len(length(breaks) - 1L)) {
    bin_labels[[i + 1L]] <- sprintf("%.0f\u2013%.0f", breaks[[i]], breaks[[i + 1L]])
  }
  bin_labels[[length(bin_labels)]] <- sprintf("> %.0f", hi_break)

  all_breaks <- c(-Inf, breaks, Inf)
  pI_bin <- cut(pI_vals, breaks = all_breaks, labels = bin_labels,
                right = TRUE, include.lowest = FALSE)

  df <- data.frame(pI = pI_vals, bin = pI_bin)

  # ── Fraction count annotation data ────────────────────────────────────────
  bin_counts <- as.integer(table(pI_bin))
  n_bins     <- length(bin_labels)

  # Mid-x for each label: for <lo and >hi, place 1 unit outside; for intervals, mid-point
  bin_mids <- numeric(n_bins)
  bin_mids[[1L]] <- lo_break - 1
  for (i in seq_len(length(breaks) - 1L)) {
    bin_mids[[i + 1L]] <- (breaks[[i]] + breaks[[i + 1L]]) / 2
  }
  bin_mids[[n_bins]] <- hi_break + 1

  ann_df <- data.frame(
    pI    = bin_mids,
    count = bin_counts,
    bin   = factor(bin_labels, levels = bin_labels),
    label = ifelse(bin_counts == 0L, "", as.character(bin_counts))
  )

  # ── Colours: viridis-based sequential across bins ─────────────────────────
  n_colours <- n_bins
  # Use a hand-picked sequential palette that works with the pepVet brand
  # Acidic (low pI) → cool blues; basic (high pI) → warm orange-reds
  viridis_cols <- grDevices::hcl.colors(n_colours, palette = "viridis", alpha = 0.85)

  # ── Auto title ────────────────────────────────────────────────────────────
  auto_title <- if (!is.null(title)) {
    title
  } else if (is.list(result) && !is.data.frame(result) &&
             "params" %in% names(result)) {
    pid    <- result$params$protein_ids[[1L]]
    enzyme <- result$params$enzyme
    paste0(.tidy_protein_id(pid), "  \u00b7  ", enzyme, "  \u00b7  pI distribution")
  } else {
    "pI distribution"
  }

  x_lo <- min(c(pI_vals, lo_break)) - 0.5
  x_hi <- max(c(pI_vals, hi_break)) + 0.5

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pI, fill = .data$bin)) +
    ggplot2::geom_histogram(
      binwidth = 0.25,
      color    = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(
      values = viridis_cols,
      name   = "SCX fraction",
      drop   = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = breaks,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    )

  # Optional fraction boundary lines
  if (isTRUE(show_fraction_lines)) {
    interior_breaks <- breaks[-c(1L, length(breaks))]
    if (length(interior_breaks) > 0L) {
      p <- p + ggplot2::geom_vline(
        xintercept = interior_breaks,
        color      = .pepvet_pal$separator,
        linewidth  = 0.5,
        linetype   = "dashed"
      )
    }
  }

  # Per-fraction count annotations (above bars, suppressed for zero-count bins)
  ann_nonzero <- ann_df[ann_df$count > 0L, , drop = FALSE]
  if (nrow(ann_nonzero) > 0L) {
    p <- p + ggplot2::geom_text(
      data        = ann_nonzero,
      ggplot2::aes(x = .data$pI, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      vjust       = 1.4,
      size        = 2.7,
      fontface    = "bold",
      color       = "#555555"
    )
  }

  p +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    ggplot2::labs(
      x     = "Isoelectric point (pI)",
      y     = "Peptide count",
      title = auto_title,
      subtitle = sprintf(
        "%d peptides  \u00b7  median pI %.2f  \u00b7  range [%.1f, %.1f]",
        length(pI_vals),
        stats::median(pI_vals),
        min(pI_vals),
        max(pI_vals)
      )
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "right",
      plot.subtitle   = ggplot2::element_text(size = 9, color = "#666666"),
      plot.title      = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark
      )
    )
}

# Private: extract pI values from a single evaluate_digest result
.pI_from_result <- function(r, length_range = c(7L, 25L)) {
  lr    <- r$params$length_range %||% length_range
  peps  <- r$peptides
  valid <- peps[peps$length >= lr[[1L]] & peps$length <= lr[[2L]], , drop = FALSE]
  if (nrow(valid) == 0L || !"peptide" %in% names(valid))
    return(numeric(0L))
  as.numeric(calculate_pI(valid$peptide))
}

# Private: multi-input pI distribution (overlaid density curves)
.plot_pI_distribution_multi <- function(results, fraction_breaks,
                                        show_fraction_lines, title) {
  rlang::check_installed("ggplot2", reason = "to use plot_pI_distribution()")

  labels <- if (!is.null(names(results))) names(results) else
    vapply(results, .result_label, character(1L))

  pI_list <- lapply(seq_along(results), function(i) {
    vals <- .pI_from_result(results[[i]])
    if (length(vals) == 0L) return(NULL)
    data.frame(pI = vals, .label = factor(labels[[i]], levels = labels),
               stringsAsFactors = FALSE)
  })
  pI_list <- Filter(Negate(is.null), pI_list)

  if (length(pI_list) == 0L) {
    .abort("No pI values found in any of the supplied results.",
                   class = "pepvet_error_invalid_digest_result")
  }
  df <- do.call(rbind, pI_list)

  breaks   <- sort(unique(as.numeric(fraction_breaks)))
  lo_break <- breaks[[1L]]
  hi_break <- breaks[[length(breaks)]]
  x_lo     <- min(c(df$pI, lo_break)) - 0.5
  x_hi     <- max(c(df$pI, hi_break)) + 0.5

  # Colour palette: one line per protein/enzyme
  n     <- length(labels)
  cols  <- grDevices::hcl.colors(n, palette = "Dark 2")
  names(cols) <- labels

  auto_title <- title %||% "pI distribution: comparison"

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$pI, color = .data$.label, fill = .data$.label
  )) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1.2)

  if (isTRUE(show_fraction_lines)) {
    interior <- breaks[-c(1L, length(breaks))]
    if (length(interior) > 0L) {
      p <- p + ggplot2::geom_vline(
        xintercept = interior,
        color = .pepvet_pal$separator, linewidth = 0.45, linetype = "dashed"
      )
    }
  }

  p +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual( values = cols, name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi)) +
    ggplot2::labs(
      title = auto_title,
      x     = "Isoelectric point (pI)",
      y     = "Density"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark)
    )
}



# ── plot_missed_cleavage_impact ───────────────────────────────────────────────

#' Missed Cleavage Impact Plot
#'
#' `plot_missed_cleavage_impact()` visualizes how allowing more missed
#' cleavages changes each component score and the composite.  The user runs
#' [evaluate_digest()] at MC = 0, 1, 2 (or more) and passes the results as a
#' named list.  Each component score is drawn as a connected line; the
#' composite score is drawn as a bold line.  An annotation marks the MC count
#' that maximizes the composite.
#'
#' @param results A named list of [evaluate_digest()] results.  Names should
#'   be the MC level (e.g., `list("MC=0" = r0, "MC=1" = r1, "MC=2" = r2)`).
#'   **Or** an unnamed list of length 2–4, in which case names are auto-assigned
#'   as `"MC=0"`, `"MC=1"`, etc.
#' @param components Character vector of component score columns to show.
#'   Defaults to `c("S_length","S_coverage","S_count","S_hydro","S_charge")`.
#' @param title Optional character title.  Auto-generated when `NULL`.
#'
#' @return A `ggplot` object.
#' @export
plot_missed_cleavage_impact <- function(
    results,
    components = c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge"),
    title      = NULL) {

  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots")

  if (!is.list(results) || length(results) < 2L) {
    .abort(
      c("!" = "{.arg results} must be a list of at least 2 {.fn evaluate_digest} results.",
        "i" = "Create one result per missed-cleavage level, e.g. MC=0, MC=1, MC=2."),
      class = "pepvet_error_invalid_digest_result"
    )
  }
  # Auto-name if unnamed
  if (is.null(names(results))) {
    names(results) <- paste0("MC=", seq_along(results) - 1L)
  }

  # Validate each leaf
  for (nm in names(results)) {
    r <- results[[nm]]
    if (!is.list(r) || !"scores" %in% names(r)) {
      .abort(
        c("!" = "Element {.val {nm}} is not a valid {.fn evaluate_digest} result.",
          "i" = "Each element must be a list with a {.code scores} component."),
        class = "pepvet_error_invalid_digest_result"
      )
    }
  }

  # ── Extract scores ───────────────────────────────────────────────────────
  valid_comp <- components[components %in% names(results[[1L]]$scores)]
  all_comp   <- c(valid_comp, "composite_score")

  rows <- lapply(names(results), function(nm) {
    sc <- results[[nm]]$scores
    as.data.frame(
      c(list(mc_label = nm),
        lapply(sc[all_comp], function(col) unname(as.numeric(col)))),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df$mc_label <- factor(df$mc_label, levels = names(results))
  df$x_idx    <- as.integer(df$mc_label)

  # ── Reshape to long ──────────────────────────────────────────────────────
  component_labels <- c(
    S_length   = "Length",
    S_coverage = "Coverage",
    S_count    = "Count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge"
  )
  long_rows <- lapply(valid_comp, function(cn) {
    data.frame(
      mc_label  = df$mc_label,
      x_idx     = df$x_idx,
      score     = df[[cn]],
      component = unname(component_labels[cn] %||% cn),
      is_composite = FALSE,
      stringsAsFactors = FALSE
    )
  })
  # Add composite as its own group
  long_rows[[length(long_rows) + 1L]] <- data.frame(
    mc_label  = df$mc_label,
    x_idx     = df$x_idx,
    score     = df$composite_score,
    component = "Composite",
    is_composite = TRUE,
    stringsAsFactors = FALSE
  )
  long_df <- do.call(rbind, long_rows)
  rownames(long_df) <- NULL

  # ── Best MC for composite ────────────────────────────────────────────────
  comp_vals     <- df$composite_score
  best_idx      <- which.max(comp_vals)
  best_mc_label <- as.character(df$mc_label[[best_idx]])
  best_score    <- comp_vals[[best_idx]]

  # ── Colors: component lines pale, composite bold brand_dark ─────────────
  n_comp  <- length(valid_comp)
  comp_colors <- stats::setNames(
    grDevices::hcl.colors(n_comp, palette = "Dark 2"),
    component_labels[valid_comp]
  )
  comp_colors["Composite"] <- .pepvet_pal$brand_dark

  # ── Line widths ──────────────────────────────────────────────────────────
  line_widths <- stats::setNames(
    rep(0.8, n_comp + 1L),
    c(component_labels[valid_comp], "Composite")
  )
  line_widths["Composite"] <- 2.2

  # ── Line types ──────────────────────────────────────────────────────────
  line_types <- stats::setNames(
    rep("solid", n_comp + 1L),
    c(component_labels[valid_comp], "Composite")
  )
  line_types["Composite"] <- "solid"

  # ── Auto title ───────────────────────────────────────────────────────────
  protein_id  <- results[[1L]]$params$protein_ids[[1L]]
  enzyme_name <- results[[1L]]$params$enzyme
  display_id  <- .tidy_protein_id(protein_id)
  auto_title  <- title %||% paste0(display_id, "    \u00b7    ",
                                    enzyme_name, "    \u00b7    Missed Cleavage Impact")

  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x     = x_idx,
      y     = score,
      color = component,
      group = component
    )
  ) +
    # Threshold reference lines
    ggplot2::geom_hline(yintercept = 0.40, linetype = "dotted",
      color = .pepvet_pal$moderate, linewidth = 0.5, alpha = 0.7) +
    ggplot2::geom_hline(yintercept = .get_param("verdict_good"), linetype = "dotted",
      color = .pepvet_pal$good, linewidth = 0.5, alpha = 0.7) +
    # Component lines
    ggplot2::geom_line(
      data = long_df[!long_df$is_composite, ],
      ggplot2::aes(linewidth = component),
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = long_df[!long_df$is_composite, ],
      size = 2.0, alpha = 0.80
    ) +
    # Composite bold line
    ggplot2::geom_line(
      data = long_df[long_df$is_composite, ],
      ggplot2::aes(linewidth = component)
    ) +
    ggplot2::geom_point(
      data = long_df[long_df$is_composite, ],
      size = 3.5, shape = 18
    ) +
    # Best-MC annotation
    ggplot2::annotate(
      "label",
      x     = best_idx,
      y     = best_score + 0.06,
      label = sprintf("Best: %s\nComposite = %.3f", best_mc_label, best_score),
      size  = 2.8,
      hjust = 0.5,
      fill  = .pepvet_pal$shade,
      color = .pepvet_pal$brand_dark,
      fontface = "bold",
      linewidth = 0.3
    ) +
    ggplot2::scale_color_manual(values = comp_colors, name = "Component") +
    ggplot2::scale_linewidth_manual(values = line_widths, guide = "none") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(levels(df$mc_label)),
      labels = levels(df$mc_label),
      expand = ggplot2::expansion(mult = c(0.12, 0.18))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = paste0("Bold line = composite  \u00b7  Dotted thresholds at 0.40 / 0.65",
                        "  \u00b7  Best MC highlighted"),
      x = "Missed cleavages allowed",
      y = "Score (0 \u2013 1)"
    ) +
    .pepvet_theme() +
    ggplot2::theme(legend.position = "right")
}


# ── plot_mz_distribution ──────────────────────────────────────────────────────

#' Precursor m/z Distribution
#'
#' `plot_mz_distribution()` draws overlapping density fills for precursor
#' m/z values at charge states \eqn{z = +2} and \eqn{z = +3}, overlaid on a
#' shaded instrument scan window.  This closes the key coverage gap in the
#' pepVet distribution suite: length, GRAVY, and pI characterise peptide
#' properties; m/z characterises instrument detectability.
#'
#' An enzyme can produce ideal-length, well-behaved peptides that still fall
#' outside the instrument's MS1 scan window at the predominant charge states.
#' This function makes that problem immediately visible.
#'
#' @param result Accepted inputs:
#'   * A named list returned by [evaluate_digest()].  Valid peptides are
#'     extracted automatically using `length_range`, and m/z values are
#'     computed via [calculate_peptide_mass()].
#'   * A named list of [evaluate_digest()] results (multi-input mode).
#'     Produces a faceted plot with one panel per result.
#'   * A data.frame / tibble with a `peptide` column.  m/z values are
#'     computed from sequences.
#'   * A data.frame / tibble with `mz` and `charge_state` columns (pre-
#'     computed m/z values; `charge_states` argument is ignored).
#' @param scan_range Numeric vector of length 2 giving the instrument's MS1
#'   scan window boundaries in m/z.  Defaults to `c(350, 1500)` (typical
#'   DDA on Orbitrap / Q-TOF instruments).  Use `c(400, 1000)` for targeted
#'   methods.
#' @param charge_states Integer vector of charge states to compute.  Defaults
#'   to `2:3`.  Ignored when `result` already contains an `mz` column.
#' @param length_range Integer vector of length 2.  Valid peptide length
#'   window used to filter input peptides.  Defaults to `c(7L, 25L)`.  Read
#'   from `result$params` when a full [evaluate_digest()] result is supplied.
#' @param show_rug Logical.  When `TRUE` (default) a rug of individual
#'   peptide m/z values is added below the density fills at `alpha = 0.30`.
#' @param title Optional character string for the plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   p   <- plot_mz_distribution(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [calculate_peptide_mass()],
#'   [plot_length_distribution()], [plot_gravy_landscape()],
#'   [plot_pI_distribution()]
#' @export
plot_mz_distribution <- function(
    result,
    scan_range    = c(350, 1500),
    charge_states = 2:3,
    length_range  = c(7L, 25L),
    show_rug      = TRUE,
    title         = NULL
) {
  rlang::check_installed("ggplot2",
    reason = "to use plot_mz_distribution()")

  # ── Validate scan_range ───────────────────────────────────────────────────
  if (!is.numeric(scan_range) || length(scan_range) != 2L ||
      anyNA(scan_range) || scan_range[[1L]] >= scan_range[[2L]]) {
    .abort(
      "{.arg scan_range} must be a numeric vector of length 2 in ascending order.",
      class = "pepvet_error_invalid_input"
    )
  }
  scan_lo <- as.numeric(scan_range[[1L]])
  scan_hi <- as.numeric(scan_range[[2L]])

  # ── Multi-input mode: named list of evaluate_digest() results ────────────
  if (.is_named_results_list(result)) {
    return(.plot_mz_distribution_multi(
      result,
      scan_range    = scan_range,
      charge_states = charge_states,
      length_range  = length_range,
      show_rug      = show_rug,
      title         = title
    ))
  }

  # ── Extract peptide data ──────────────────────────────────────────────────
  if (is.list(result) && !is.data.frame(result) &&
      all(c("peptides", "params") %in% names(result))) {
    peps         <- result$peptides
    length_range <- result$params$length_range %||% length_range
    auto_label   <- if (!is.null(result$params$protein_ids)) {
      pid    <- result$params$protein_ids[[1L]]
      enzyme <- result$params$enzyme
      paste0(.tidy_protein_id(pid), "  \u00b7  ", enzyme)
    } else { NULL }

  } else if (is.data.frame(result)) {
    peps       <- result
    auto_label <- NULL

  } else {
    .abort(
      c(
        paste0("{.arg result} must be an {.fn evaluate_digest} list, a",
          " named list of such results, or a data.frame with a",
          " {.field peptide} column."),
        "x" = "Got {.cls {class(result)[[1L]]}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }

  # ── If data already has mz + charge_state, use directly ──────────────────
  if (all(c("mz", "charge_state") %in% names(peps))) {
    mz_long <- peps[, c("mz", "charge_state"), drop = FALSE]
    mz_long$charge_state <- as.character(mz_long$charge_state)
    n_total <- nrow(mz_long)

  } else {
    # Filter to valid length range, compute m/z per charge state
    if (!"peptide" %in% names(peps)) {
      .abort(
        paste0("{.arg result} must contain a {.field peptide} column",
          " to compute m/z values."),
        class = "pepvet_error_invalid_digest_result"
      )
    }

    length_lo <- as.integer(length_range[[1L]])
    length_hi <- as.integer(length_range[[2L]])

    if ("length" %in% names(peps)) {
      valid_peps <- peps[peps$length >= length_lo & peps$length <= length_hi,
                         , drop = FALSE]
    } else {
      valid_peps <- peps
    }

    if (nrow(valid_peps) == 0L) {
      .abort(
        "No valid peptides found in {.arg result} to compute m/z values.",
        class = "pepvet_error_invalid_digest_result"
      )
    }

    charge_states_int <- sort(unique(as.integer(charge_states)))
    n_total           <- nrow(valid_peps)

    mz_long <- do.call(rbind, lapply(charge_states_int, function(z) {
      mz_vals <- as.numeric(
        calculate_peptide_mass(valid_peps$peptide, charge = z))
      data.frame(
        mz           = mz_vals,
        charge_state = paste0("z = +", z),
        stringsAsFactors = FALSE
      )
    }))
  }

  mz_long$charge_state <- factor(mz_long$charge_state,
    levels = unique(mz_long$charge_state))

  # ── Assign brand colors to charge states ─────────────────────────────────
  # z=+2 = primary brand color; z=+3 = brand_light; z=+4+ = moderate, etc.
  z_color_palette <- c(
    .pepvet_pal$brand,
    .pepvet_pal$brand_light,
    .pepvet_pal$moderate,
    .pepvet_pal$poor
  )
  z_levels  <- levels(mz_long$charge_state)
  z_colors  <- setNames(
    z_color_palette[seq_along(z_levels)],
    z_levels
  )

  # ── Per-charge-state % inside scan window ─────────────────────────────────
  window_stats <- do.call(rbind, lapply(z_levels, function(z) {
    sub    <- mz_long$mz[mz_long$charge_state == z]
    n_in   <- sum(sub >= scan_lo & sub <= scan_hi, na.rm = TRUE)
    n_all  <- sum(!is.na(sub))
    pct    <- if (n_all > 0L) round(100 * n_in / n_all, 1) else NA_real_
    data.frame(charge_state = z, n_in = n_in, n_all = n_all,
               pct = pct, stringsAsFactors = FALSE)
  }))

  # ── x-axis range: pad to nearest 100 m/z beyond data, min at 0 ───────────
  mz_min <- max(0,   floor(min(mz_long$mz, na.rm = TRUE)  / 100) * 100 - 50)
  mz_max <- ceiling(max(mz_long$mz, na.rm = TRUE) / 100) * 100 + 50
  x_lo   <- min(mz_min, scan_lo - 50)
  x_hi   <- max(mz_max, scan_hi + 50)

  # ── Annotation label: stacked text per charge state ───────────────────────
  ann_text <- paste(
    vapply(seq_len(nrow(window_stats)), function(i) {
      sprintf("%s: %.0f%% within window (%d / %d peptides)",
        window_stats$charge_state[[i]],
        window_stats$pct[[i]],
        window_stats$n_in[[i]],
        window_stats$n_all[[i]])
    }, character(1L)),
    collapse = "\n"
  )

  # ── Auto title ────────────────────────────────────────────────────────────
  auto_title <- if (!is.null(title)) {
    title
  } else if (!is.null(auto_label)) {
    paste0(auto_label, "  \u00b7  Precursor m/z distribution")
  } else {
    "Precursor m/z distribution"
  }

  n_peptides_label <- if (!all(c("mz", "charge_state") %in% names(
      if (is.data.frame(result)) result else result$peptides))) {
    n_total
  } else {
    nrow(mz_long) / length(z_levels)
  }

  subtitle_text <- sprintf(
    "%d valid peptides  \u00b7  scan window %.0f\u2013%.0f m/z  \u00b7  z = %s",
    n_total,
    scan_lo, scan_hi,
    paste(gsub("z = \\+", "+", z_levels), collapse = " & ")
  )

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    mz_long,
    ggplot2::aes(x = mz, color = charge_state, fill = charge_state)
  ) +
    # ── Background: outside-window zone shading (neutral) ─────────────────
    ggplot2::annotate(
      "rect", xmin = x_lo, xmax = scan_lo,
      ymin = -Inf, ymax = Inf,
      fill  = "#FFF3E0", alpha = 0.45
    ) +
    ggplot2::annotate(
      "rect", xmin = scan_hi, xmax = x_hi,
      ymin = -Inf, ymax = Inf,
      fill  = "#FFF3E0", alpha = 0.45
    ) +
    # ── Valid window shading (green) ──────────────────────────────────────
    ggplot2::annotate(
      "rect", xmin = scan_lo, xmax = scan_hi,
      ymin = -Inf, ymax = Inf,
      fill  = .pepvet_pal$shade, alpha = 0.50
    ) +
    # ── Window boundary lines ─────────────────────────────────────────────
    ggplot2::geom_vline(
      xintercept = scan_lo,
      color      = .pepvet_pal$good,
      linewidth  = 0.65,
      linetype   = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = scan_hi,
      color      = .pepvet_pal$poor,
      linewidth  = 0.65,
      linetype   = "dashed"
    ) +
    # ── Density fills (overlapping, semi-transparent) ─────────────────────
    ggplot2::geom_density(
      alpha   = 0.40,
      linewidth = 0.90,
      adjust  = 0.9
    ) +
    # ── % inside window annotation block ──────────────────────────────────
    ggplot2::annotate(
      "text",
      x        = scan_lo + (scan_hi - scan_lo) * 0.97,
      y        = Inf,
      hjust    = 1,
      vjust    = 1.5,
      label    = ann_text,
      size     = 2.7,
      fontface = "bold",
      color    = .pepvet_pal$brand_dark,
      lineheight = 1.45
    ) +
    # ── Zone boundary labels ──────────────────────────────────────────────
    ggplot2::annotate(
      "text",
      x = scan_lo + (scan_hi - scan_lo) * 0.015,
      y = Inf, hjust = 0, vjust = 1.6,
      label    = sprintf("%.0f m/z", scan_lo),
      size     = 2.5,
      fontface = "italic",
      color    = .pepvet_pal$good
    ) +
    ggplot2::annotate(
      "text",
      x = scan_hi - (scan_hi - scan_lo) * 0.015,
      y = Inf, hjust = 1, vjust = 1.6,
      label    = sprintf("%.0f m/z", scan_hi),
      size     = 2.5,
      fontface = "italic",
      color    = .pepvet_pal$poor
    ) +
    # ── Rug marks ─────────────────────────────────────────────────────────
    {
      if (isTRUE(show_rug)) {
        ggplot2::geom_rug(
          sides  = "b",
          alpha  = 0.30,
          length = ggplot2::unit(0.02, "npc")
        )
      }
    } +
    # ── Scales ────────────────────────────────────────────────────────────
    ggplot2::scale_color_manual(
      values = z_colors,
      name   = "Charge state",
      guide  = ggplot2::guide_legend(
        override.aes = list(fill = z_colors, alpha = 1, linewidth = 1.2)
      )
    ) +
    ggplot2::scale_fill_manual(
      values = z_colors,
      name   = "Charge state"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, x_hi + 200, by = 200),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.22))
    ) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
    ggplot2::labs(
      title    = auto_title,
      subtitle = subtitle_text,
      x        = "Precursor m/z",
      y        = "Density"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title      = ggplot2::element_text(
        size  = 13,
        face  = "bold",
        color = .pepvet_pal$brand_dark
      )
    )

  p
}

# Private: multi-input mz distribution (faceted)
.plot_mz_distribution_multi <- function(
    results,
    scan_range,
    charge_states,
    length_range,
    show_rug,
    title
) {
  rlang::check_installed("ggplot2",
    reason = "to use plot_mz_distribution()")

  labels <- if (!is.null(names(results))) names(results) else
    vapply(results, .result_label, character(1L))

  scan_lo <- as.numeric(scan_range[[1L]])
  scan_hi <- as.numeric(scan_range[[2L]])

  charge_states_int <- sort(unique(as.integer(charge_states)))

  mz_all <- do.call(rbind, lapply(seq_along(results), function(i) {
    r    <- results[[i]]
    lr   <- r$params$length_range %||% length_range
    peps <- r$peptides
    peps <- peps[peps$length >= lr[[1L]] & peps$length <= lr[[2L]], ,
                 drop = FALSE]
    if (nrow(peps) == 0L || !"peptide" %in% names(peps))
      return(NULL)

    rows <- do.call(rbind, lapply(charge_states_int, function(z) {
      mz_vals <- as.numeric(
        calculate_peptide_mass(peps$peptide, charge = z))
      data.frame(
        mz           = mz_vals,
        charge_state = paste0("z = +", z),
        .label       = factor(labels[[i]], levels = labels),
        stringsAsFactors = FALSE
      )
    }))
    rows
  }))

  if (is.null(mz_all) || nrow(mz_all) == 0L) {
    .abort("No valid peptide m/z values could be computed.",
      class = "pepvet_error_invalid_digest_result")
  }

  z_levels <- paste0("z = +", charge_states_int)
  mz_all$charge_state <- factor(mz_all$charge_state, levels = z_levels)

  z_color_palette <- c(
    .pepvet_pal$brand,
    .pepvet_pal$brand_light,
    .pepvet_pal$moderate,
    .pepvet_pal$poor
  )
  z_colors <- setNames(
    z_color_palette[seq_along(z_levels)],
    z_levels
  )

  x_lo <- max(0, floor(min(mz_all$mz, na.rm = TRUE) / 100) * 100 - 50)
  x_hi <- ceiling(max(mz_all$mz, na.rm = TRUE) / 100) * 100 + 50
  x_lo <- min(x_lo, scan_lo - 50)
  x_hi <- max(x_hi, scan_hi + 50)

  auto_title <- if (!is.null(title)) title else
    "Precursor m/z distribution: comparison"

  ggplot2::ggplot(
    mz_all,
    ggplot2::aes(x = mz, color = charge_state, fill = charge_state)
  ) +
    ggplot2::annotate("rect",
      xmin = x_lo, xmax = scan_lo, ymin = -Inf, ymax = Inf,
      fill = "#FFF3E0", alpha = 0.40) +
    ggplot2::annotate("rect",
      xmin = scan_hi, xmax = x_hi, ymin = -Inf, ymax = Inf,
      fill = "#FFF3E0", alpha = 0.40) +
    ggplot2::annotate("rect",
      xmin = scan_lo, xmax = scan_hi, ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.45) +
    ggplot2::geom_vline(xintercept = scan_lo,
      color = .pepvet_pal$good, linewidth = 0.55, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = scan_hi,
      color = .pepvet_pal$poor, linewidth = 0.55, linetype = "dashed") +
    ggplot2::geom_density(alpha = 0.38, linewidth = 0.75, adjust = 0.9) +
    {
      if (isTRUE(show_rug)) {
        ggplot2::geom_rug(sides = "b", alpha = 0.20,
          length = ggplot2::unit(0.025, "npc"))
      }
    } +
    ggplot2::facet_wrap(ggplot2::vars(.data$.label), scales = "free_y") +
    ggplot2::scale_color_manual(values = z_colors, name = "Charge state",
      guide = ggplot2::guide_legend(
        override.aes = list(fill = z_colors, alpha = 1, linewidth = 1.2))) +
    ggplot2::scale_fill_manual(values = z_colors, name = "Charge state") +
    ggplot2::scale_x_continuous(
      breaks = seq(0, x_hi + 200, by = 200),
      expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.20))) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi)) +
    ggplot2::labs(
      title    = auto_title,
      subtitle = sprintf("Scan window %.0f\u2013%.0f m/z  \u00b7  z = %s",
        scan_lo, scan_hi,
        paste(gsub("z = \\+", "+", z_levels), collapse = " & ")),
      x = "Precursor m/z",
      y = "Density"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(
        size = 9, face = "bold", color = .pepvet_pal$brand_dark),
      plot.title = ggplot2::element_text(
        size = 13, face = "bold", color = .pepvet_pal$brand_dark)
    )
}
