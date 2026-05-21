#' Four-Panel Digest Diagnostic Plot
#'
#' `plot_digest_profile()` assembles a four-panel figure for a single
#' protein–enzyme pair from an [evaluate_digest()] result.  The panels are:
#'
#' * **(A) Length distribution:** histogram of peptide lengths with the valid
#'   window shaded.  Bars are colored by length class: valid (blue), too short
#'   (amber), too long (rose).
#' * **(B) GRAVY distribution:** histogram of GRAVY hydrophobicity scores.
#'   The LC-friendly range is shaded and bounded by dashed lines.
#' * **(C) Coverage map:** protein drawn as horizontal lanes (one per
#'   missed-cleavage level) with valid-length peptides stacked via a greedy
#'   interval-packing algorithm.  Uncovered regions are highlighted in red.
#'   Peptide length labels appear inside segments of 8 aa or longer.
#' * **(D) Component scores:** horizontal bar chart for each scoring
#'   component, colored by tier (green \eqn{\geq} 0.65, amber 0.40–0.64, red
#'   < 0.40).  The composite score is marked with a dashed vertical line.
#'
#' @param result A named list returned by [evaluate_digest()].  Must describe
#'   a single protein (one unique `protein_id` in `result$peptides`).
#' @param length_range Integer vector of length 2.  Defines the valid peptide
#'   length window, passed to the length and coverage panels.  Defaults to
#'   `c(7L, 25L)`.
#' @param gravy_range Numeric vector of length 2.  Defines the LC-friendly
#'   GRAVY range shaded in panel B.  Defaults to `c(-1.0, 0.6)`.
#' @param title Optional character string for the figure title.  When `NULL`
#'   (default) a title is auto-generated from the protein accession, enzyme,
#'   and missed-cleavage count.
#'
#' @return A `patchwork` object combining all four ggplot panels.  The object
#'   can be printed directly, saved with [ggplot2::ggsave()], or composed
#'   further with other patchwork operators.
#'
#' @details GRAVY scores are computed internally from the peptide sequences in
#'   `result$peptides` using the Kyte–Doolittle scale.  No external columns are
#'   required beyond the standard [evaluate_digest()] output.
#'
#'   Panel C labels peptide lengths inside segments of 8 aa or longer.  For
#'   heavily digested proteins this keeps the map readable without overlap.
#'   When multiple missed-cleavage levels are present, each level occupies its
#'   own horizontal lane with peptides stacked using the same greedy
#'   interval-packing algorithm as [plot_coverage_map()].
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'   requireNamespace("patchwork", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin", missed_cleavages = 1L)
#'   p <- plot_digest_profile(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [plot_coverage_map()],
#'   [plot_enzyme_comparison()]
#' @family plot-single
#' @export

plot_digest_profile <- function(result,
                                length_range = c(7L, 25L),
                                gravy_range = c(-1.0, 0.6),
                                title = NULL) {
  rlang::check_installed(
    "ggplot2",
    reason = "to produce pepVet visualization plots"
  )
  rlang::check_installed(
    "patchwork",
    reason = "to assemble multi-panel figures in plot_digest_profile()"
  )

  .validate_digest_result_for_plot(result)

  length_range <- .validate_length_range(length_range)
  gravy_range  <- .validate_gravy_range(gravy_range)

  peps   <- result$peptides
  scores <- result$scores
  params <- result$params

  protein_id  <- params$protein_ids[[1L]]
  enzyme_name <- params$enzyme
  display_id <- .tidy_protein_id(protein_id)

  # Add GRAVY scores to the peptide table (computed from sequences)
  peps$gravy <- .calculate_gravy_vec(peps$peptide)

  protein_length <- max(peps$end, na.rm = TRUE)

  # ── Build panels ──────────────────────────────────────────────────────────
  pa <- .panel_length(peps, length_range)
  pb <- .panel_gravy(peps, gravy_range)
  mc_val <- params$missed_cleavages %||% 1L
  pc <- .panel_coverage(peps, protein_length, length_range,
    missed_cleavages = mc_val)
  pd <- .panel_scores(scores)

  # ── Assemble with patchwork ───────────────────────────────────────────────
  # Layout: top row = A | B (two equal columns)
  #         middle  = C     (full width)
  #         bottom  = D     (full width)
  figure <- (pa | pb) /
    pc /
    pd +
    patchwork::plot_layout(heights = c(3, 1.8, 2.2))

  auto_title <- if (is.null(title)) {
    paste0(display_id, "    \u00b7    ", enzyme_name,
      "    \u00b7    MC = ", mc_val)
  } else {
    title
  }

  figure +
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
        face   = "bold",
        size   = .get_param("patchwork_tag_size"),
        color  = .pepvet_pal$brand,
        margin = ggplot2::margin(t = 1, r = 6, b = 1, l = 2)
      )
    )
}


#' Protein Coverage Map
#'
#' `plot_coverage_map()` draws a genome-browser-style view of how proteolytic
#' peptides map onto the full protein sequence.  When the input contains
#' multiple missed-cleavage levels (e.g. MC = 0, 1, 2), each level occupies a
#' separate horizontal lane so the cumulative coverage gain from allowing missed
#' cleavages is immediately visible.  Cleavage-site efficiency ticks and
#' optional domain annotations can be overlaid to add mechanistic context.
#'
#' @param result A named list returned by [evaluate_digest()].  Must describe
#'   a single protein.
#' @param color_by Character string selecting the peptide color scheme.
#'   \describe{
#'     \item{`"validity"` (default)}{Blue = valid-length peptide; gray =
#'       invalid-length peptide.}
#'     \item{`"length_class"`}{Three-way coloring: valid (blue), too short
#'       (amber), too long (rose-red).}
#'     \item{`"hydrophobicity"`}{Continuous GRAVY gradient for every peptide:
#'       brand blue (very hydrophilic) → green (LC-optimal) → amber
#'       (borderline) → rose-red (very hydrophobic).}
#'   }
#' @param length_range Integer vector of length 2.  Valid peptide window.
#'   Defaults to `c(7L, 25L)`.
#' @param cleavage_sites Optional data.frame from [annotate_cleavage_sites()].
#'   When provided, vertical ticks below the lanes mark each cleavage site,
#'   colored by efficiency: high = green, medium = amber, low = rose-red.
#' @param domains Optional data.frame with columns `name`, `start`, `end`
#'   (integer, in residue positions).  Each domain is drawn as a translucent
#'   background rectangle spanning all lanes with its name labelled at the top.
#' @param title Optional character string for the plot title.  Auto-generated
#'   from the protein accession and enzyme when `NULL` (default).
#'
#' @return A `ggplot` object that can be printed, further customised, or saved
#'   with [ggplot2::ggsave()].
#'
#' @details Valid peptides are drawn taller than invalid peptides within each
#'   lane, maintaining a visual hierarchy that keeps validity legible regardless
#'   of `color_by`.  When missed-cleavage levels overlap (MC ≥ 1), peptides
#'   within each lane are distributed into non-overlapping tracks using a greedy
#'   interval-packing algorithm, mirroring genome-browser stacking.  Gap regions
#'   (residues not covered by any valid MC = 0 peptide) are highlighted with a
#'   translucent red overlay.  Peptide lengths are labelled inside bars whose
#'   span is at least 2.5 % of the protein length (ensuring legibility at
#'   typical export widths); shorter bars are left unlabelled.
#'
#'   The `color_by = "hydrophobicity"` mode calls the internal Kyte–Doolittle
#'   GRAVY calculator and requires no additional input columns.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin", missed_cleavages = 1)
#'   cs <- annotate_cleavage_sites(bsa_path, enzyme = "trypsin")
#'   p <- plot_coverage_map(res, cleavage_sites = cs)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [annotate_cleavage_sites()],
#'   [plot_digest_profile()]
#' @family plot-single
#' @export
plot_coverage_map <- function(result,
                              color_by = c(
                                "validity",
                                "length_class",
                                "hydrophobicity"
                              ),
                              length_range = c(7L, 25L),
                              cleavage_sites = NULL,
                              domains = NULL,
                              title = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots"
  )

  color_by <- match.arg(color_by)

  .validate_digest_result_for_plot(result)

  length_range <- .validate_length_range(length_range)

  # ── Input validation for optional args ─────────────────────────────────────
  if (!is.null(cleavage_sites)) {
    if (
      !is.data.frame(cleavage_sites) ||
        !all(c("position", "efficiency") %in% names(cleavage_sites))
    ) {
      .abort(
        c(
          "!" = "{.arg cleavage_sites} must be a data.frame from {.fn annotate_cleavage_sites}.",
          "i" = "Required columns: {.code position}, {.code efficiency}."
        ),
        class = "pepvet_error_invalid_cleavage_sites"
      )
    }
  }
  if (!is.null(domains)) {
    if (
      !is.data.frame(domains) ||
        !all(c("name", "start", "end") %in% names(domains))
    ) {
      .abort(
        c(
          "!" = "{.arg domains} must be a data.frame with columns {.code name}, {.code start}, {.code end}.",
          "i" = "Each row describes one annotated protein domain."
        ),
        class = "pepvet_error_invalid_domains"
      )
    }
  }

  # ── Extract data ────────────────────────────────────────────────────────────
  peps <- result$peptides
  params <- result$params
  protein_id <- params$protein_ids[[1L]]
  enzyme_name <- params$enzyme
  display_id <- .tidy_protein_id(protein_id)
  protein_length <- max(peps$end, na.rm = TRUE)

  # ── GRAVY (computed once, used for coloring or silently ignored) ────────────
  peps$gravy <- .calculate_gravy_vec(peps$peptide)

  # ── Determine MC levels present and build lane coordinates ─────────────────
  has_mc <- "missed_cleavages" %in% names(peps)
  mc_levels <- if (has_mc) sort(unique(peps$missed_cleavages[!is.na(peps$missed_cleavages)])) else 0L
  tick_h <- if (!is.null(cleavage_sites)) 0.11 else 0.0
  lanes <- .lane_y_coords(mc_levels, tick_height = tick_h)

  # ── Coverage stats on MC = 0 peptides (for subtitle) ───────────────────────
  cs0 <- .compute_coverage_stats(peps, protein_length, length_range,
    mc_filter = if (has_mc) 0L else NULL
  )
  pct_cov <- cs0$pct_cov
  gap_df <- cs0$gap_df

  # ── Build fill aesthetic and scale ─────────────────────────────────────────
  #   For "validity"/"length_class": factor column `fill_cat` → scale_fill_manual
  #   For "hydrophobicity":          numeric column `fill_val` (GRAVY)
  #                                  → scale_fill_gradientn
  if (color_by == "hydrophobicity") {
    peps$fill_val <- peps$gravy
  } else {
    # Validity categories (used by both "validity" and "length_class")
    lo <- length_range[[1L]]
    hi <- length_range[[2L]]
    peps$fill_cat <- factor(
      ifelse(peps$length < lo, "Too short",
        ifelse(peps$length > hi, "Too long", "Valid")
      ),
      levels = c("Valid", "Too short", "Too long", "Invalid")
    )
    if (color_by == "validity") {
      levels(peps$fill_cat)[levels(peps$fill_cat) %in% c("Too short", "Too long")] <-
        "Invalid"
    }
  }

  # x-axis step
  x_step <- max(50L, as.integer(round(protein_length / 10.0 / 50.0) * 50L))

  # ── Initialize base plot ───────────────────────────────────────────────────
  p <- ggplot2::ggplot()

  # ── Domain backgrounds (behind everything else) ────────────────────────────
  # Pre-defined palette of 8 soft distinguishable fills for domains
  domain_fills <- c(
    "#D9EAF7", "#FFF0C2", "#D4EDD4", "#F7D9D9",
    "#EDD4F7", "#D4F7F0", "#F7ECD4", "#D4D9F7"
  )
  if (!is.null(domains)) {
    n_dom <- nrow(domains)
    for (k in seq_len(n_dom)) {
      dom_fill <- domain_fills[((k - 1L) %% length(domain_fills)) + 1L]
      p <- p +
        ggplot2::annotate("rect",
          xmin = domains$start[[k]], xmax = domains$end[[k]],
          ymin = tick_h, ymax = 1.0,
          fill = dom_fill, color = NA, alpha = 0.55
        ) +
        ggplot2::annotate("text",
          x = (domains$start[[k]] + domains$end[[k]]) / 2.0,
          y = 0.995,
          label = domains$name[[k]],
          size = 2.8, hjust = 0.5, vjust = 1,
          color = .pepvet_pal$text_axis_title, fontface = "italic"
        )
    }
  }

  # ── Draw each lane ─────────────────────────────────────────────────────────
  # IMPORTANT: all y-positions used inside aes() must be stored as data-frame
  # columns accessed via .data$, not as loop-scoped variables.  ggplot2 layers
  # are lazy: aes() expressions are evaluated at render time, after the loop
  # completes, so bare loop variables always capture the final iteration value.
  for (k in seq_len(nrow(lanes))) {
    mc_val <- lanes$mc[[k]]
    y_lo <- lanes$y_lo[[k]]
    y_hi <- lanes$y_hi[[k]]
    y_mid <- lanes$y_mid[[k]]
    lane_h <- y_hi - y_lo

    # Peptides for this lane
    lane_peps <- if (has_mc) {
      peps[peps$missed_cleavages == mc_val, ,
        drop = FALSE
      ]
    } else {
      peps
    }

    lane_valid <- lane_peps[
      lane_peps$length >= length_range[[1L]] &
        lane_peps$length <= length_range[[2L]], ,
      drop = FALSE
    ]
    lane_invalid <- lane_peps[
      lane_peps$length < length_range[[1L]] |
        lane_peps$length > length_range[[2L]], ,
      drop = FALSE
    ]

    # ── Greedy packing: stack overlapping peptides into sub-rows ─────────────
    # Each peptide gets a `track` integer (1 = bottom, 2 = above, ...).
    # For MC=0 tryptic digests there are no overlaps so n_tracks == 1.
    # For MC>=1 merged peptides share residues; packing yields 2-4 tracks.
    # y-bounds are pre-computed into data-frame columns to avoid the closure trap.
    if (nrow(lane_invalid) > 0L) {
      lane_invalid <- .pack_peptides(lane_invalid)
      n_tracks_i <- max(lane_invalid$track)
      inv_inner_lo <- y_mid - 0.28 * lane_h
      inv_inner_hi <- y_mid + 0.28 * lane_h
      track_h_i <- (inv_inner_hi - inv_inner_lo) / n_tracks_i
      margin_i <- 0.05 * track_h_i
      lane_invalid$.y_lo <- inv_inner_lo +
        (lane_invalid$track - 1L) * track_h_i + margin_i
      lane_invalid$.y_hi <- inv_inner_lo +
        lane_invalid$track * track_h_i - margin_i
    }
    if (nrow(lane_valid) > 0L) {
      lane_valid <- .pack_peptides(lane_valid)
      n_tracks_v <- max(lane_valid$track)
      inner_lo <- y_lo + 0.04 * lane_h
      inner_hi <- y_hi - 0.04 * lane_h
      track_h_v <- (inner_hi - inner_lo) / n_tracks_v
      margin_v <- 0.05 * track_h_v
      lane_valid$.y_lo <- inner_lo +
        (lane_valid$track - 1L) * track_h_v + margin_v
      lane_valid$.y_hi <- inner_lo +
        lane_valid$track * track_h_v - margin_v
      lane_valid$.y_mid <- (lane_valid$.y_lo + lane_valid$.y_hi) / 2.0
    }

    # Protein backbone (thin bar at lane centre) – annotate() is eager, safe
    p <- p + ggplot2::annotate("rect",
      xmin = 0.5, xmax = protein_length + 0.5,
      ymin = y_mid - 0.07 * lane_h, ymax = y_mid + 0.07 * lane_h,
      fill = .pepvet_pal$backbone_fill, color = .pepvet_pal$backbone_brd, linewidth = 0.35
    )

    # ── Gap overlays (annotate() is eager, safe with loop vars) ──────────────
    if (nrow(gap_df) > 0L) {
      p <- p + ggplot2::annotate("rect",
        xmin = gap_df$xmin, xmax = gap_df$xmax,
        ymin = y_lo + 0.03 * lane_h, ymax = y_hi - 0.03 * lane_h,
        fill = .pepvet_pal$gap, alpha = 0.12
      )
    }

    # ── Invalid peptides (thinner bars, behind valid) ─────────────────────────
    if (nrow(lane_invalid) > 0L) {
      if (color_by == "hydrophobicity") {
        p <- p + ggplot2::geom_rect(
          data = lane_invalid,
          ggplot2::aes(
            xmin = .data$start - 0.3, xmax = .data$end + 0.3,
            ymin = .data$.y_lo,       ymax = .data$.y_hi,
            fill = .data$fill_val
          ),
          alpha = 0.65, color = "white", linewidth = 0.1
        )
      } else {
        p <- p + ggplot2::geom_rect(
          data = lane_invalid,
          ggplot2::aes(
            xmin = .data$start - 0.3, xmax = .data$end + 0.3,
            ymin = .data$.y_lo,       ymax = .data$.y_hi,
            fill = .data$fill_cat
          ),
          alpha = 0.65, color = "white", linewidth = 0.1
        )
      }
    }

    # ── Valid peptides (taller bars, on top) ──────────────────────────────────
    if (nrow(lane_valid) > 0L) {
      if (color_by == "hydrophobicity") {
        p <- p + ggplot2::geom_rect(
          data = lane_valid,
          ggplot2::aes(
            xmin = .data$start - 0.4, xmax = .data$end + 0.4,
            ymin = .data$.y_lo,       ymax = .data$.y_hi,
            fill = .data$fill_val
          ),
          alpha = 0.90, color = "white", linewidth = 0.15
        )
      } else {
        p <- p + ggplot2::geom_rect(
          data = lane_valid,
          ggplot2::aes(
            xmin = .data$start - 0.4, xmax = .data$end + 0.4,
            ymin = .data$.y_lo,       ymax = .data$.y_hi,
            fill = .data$fill_cat
          ),
          alpha = 0.90, color = "white", linewidth = 0.15
        )
      }

      # Label peptides whose bar is wide enough to hold the digit(s).
      # A fixed aa threshold (e.g. >= 8) is misleading: for a long protein every
      # short peptide bar is too narrow.  Use a proportional minimum instead:
      # at least 2.5% of the protein length guarantees ~5-6 pixels per digit at
      # typical export widths (14 in / 150 dpi).  Also suppress labels when
      # packing produces many tracks and bars become too short vertically.
      min_label_span <- max(5L, as.integer(protein_length * 0.025))
      label_v <- lane_valid[
        (lane_valid$end - lane_valid$start + 1L) >= min_label_span, ,
        drop = FALSE
      ]
      label_v$label_x <- (label_v$start + label_v$end) / 2.0
      if (nrow(label_v) > 0L && n_tracks_v <= 4L) {
        p <- p + ggplot2::geom_text(
          data = label_v,
          ggplot2::aes(x = .data$label_x, y = .data$.y_mid, label = .data$length),
          size = 2.3, color = "white", fontface = "bold"
        )
      }
    }

    # Lane label on the right margin
    mc_label <- if (has_mc) sprintf("MC = %d", mc_val) else "Peptides"
    p <- p + ggplot2::annotate("text",
      x = protein_length + 3L,
      y = y_mid,
      label = mc_label,
      hjust = 0, vjust = 0.5,
      size = 3.0, color = .pepvet_pal$text_axis_title, fontface = "bold"
    )
  }

  # ── Lane separator lines (between lanes only) ─────────────────────────────
  if (nrow(lanes) > 1L) {
    for (k in seq_len(nrow(lanes) - 1L)) {
      sep_y <- (lanes$y_hi[[k]] + lanes$y_lo[[k + 1L]]) / 2.0
      p <- p + ggplot2::annotate("segment",
        x = 0, xend = protein_length,
        y = sep_y, yend = sep_y,
        color = .pepvet_pal$separator, linewidth = 0.4, linetype = "dotted"
      )
    }
  }

  # ── Cleavage-site ticks ────────────────────────────────────────────────────
  # efficiency is used as a proper color aesthetic (not I()) so ggplot2 adds it
  # to the legend panel alongside the fill legend.
  has_cs <- !is.null(cleavage_sites) && tick_h > 0
  if (has_cs) {
    cs_df <- cleavage_sites
    cs_df$eff_level <- factor(
      cs_df$efficiency,
      levels = c("high", "medium", "low"),
      labels = c("High", "Medium", "Low")
    )
    tick_top <- tick_h * 0.85

    p <- p +
      ggplot2::geom_segment(
        data = cs_df,
        ggplot2::aes(
          x = .data$position, xend = .data$position,
          y = 0, yend = tick_top,
          color = .data$eff_level
        ),
        linewidth = 0.6, alpha = 0.85
      ) +
      # Thin line separating tick zone from lane area
      ggplot2::annotate("segment",
        x = 0, xend = protein_length,
        y = tick_h, yend = tick_h,
        color = .pepvet_pal$separator, linewidth = 0.35
      )
  }

  # ── Fill scale ─────────────────────────────────────────────────────────────
  if (color_by == "hydrophobicity") {
    # Gradient: brand blue → green → amber → poor red
    # 4 evenly-spaced color stops spanning the data GRAVY range
    rescaled <- seq(0, 1, length.out = 4)
    p <- p + ggplot2::scale_fill_gradientn(
      colors = c(
        .pepvet_pal$brand, .pepvet_pal$good,
        .pepvet_pal$moderate, .pepvet_pal$poor
      ),
      values = rescaled,
      name = "GRAVY",
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(160, "pt"),
        barheight = ggplot2::unit(7, "pt"),
        direction = "horizontal",
        title.position = "top"
      )
    )
  } else {
    color_vals <- c(
      "Valid"     = .pepvet_pal$valid,
      "Too short" = .pepvet_pal$too_short,
      "Too long"  = .pepvet_pal$too_long,
      "Invalid"   = .pepvet_pal$na_gray
    )
    if (color_by == "validity") {
      color_vals <- color_vals[c("Valid", "Invalid")]
    }
    p <- p + ggplot2::scale_fill_manual(
      values = color_vals,
      name = NULL,
      na.value = .pepvet_pal$na_gray,
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, color = NA, size = 5)
      )
    )
  }

  # ── Cleavage efficiency color scale (only when ticks are shown) ────────────
  if (has_cs) {
    p <- p + ggplot2::scale_color_manual(
      values = c(
        "High"   = .pepvet_pal$good,
        "Medium" = .pepvet_pal$moderate,
        "Low"    = .pepvet_pal$poor
      ),
      name = "Cleavage efficiency",
      guide = ggplot2::guide_legend(
        override.aes = list(linewidth = 2, alpha = 1)
      )
    )
  } else {
    p <- p + ggplot2::guides(color = "none")
  }

  # ── Scales, axes, theme ────────────────────────────────────────────────────
  n_sites_str <- if (!is.null(cleavage_sites)) {
    eff_tbl <- table(cleavage_sites$efficiency)
    parts <- paste0(as.integer(eff_tbl), " ", tolower(names(eff_tbl)))
    paste0("  \u00b7  ", paste(parts, collapse = " / "), " efficiency sites")
  } else {
    ""
  }

  auto_title <- if (is.null(title)) {
    paste0(display_id, "    \u00b7    ", enzyme_name)
  } else {
    title
  }

  p +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, protein_length + x_step, by = x_step),
      limits = c(0, protein_length + protein_length * 0.10),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::labs(
      title = auto_title,
      subtitle = sprintf(
        paste(
          "%.0f%% sequence coverage (MC=0 valid peptides).",
          "%d gap region(s). Protein length %d aa%s"
        ),
        pct_cov, nrow(gap_df), protein_length, n_sites_str
      ),
      x = "Residue position",
      y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position    = "bottom"
    )
}


#' Amino-Acid Peptide Overlap Map
#'
#' `plot_peptide_overlap_map()` renders a wrapped residue-level view of a
#' single protein and colors each amino-acid tile by how many peptides
#' cover that residue (MC=0 only by default). Uncovered residues are white;
#' increasing overlap depth gets progressively stronger blue fills.
#' The subtitle reports the maximum overlap, which helps gauge whether the
#' 6-color scale (0, 1, 2-3, 4-7, 8-15, 16+) is sufficient for your data.
#'
#' @param result A named list returned by [evaluate_digest()]. Must describe a
#'   single protein.
#' @param length_range Optional integer vector of length 2 defining which
#'   peptide lengths count toward overlap. Defaults to `c(7L, 25L)`. Peptides
#'   are counted from the MC level specified by `missed_cleavages`. Use
#'   `length_range = NULL` to count all digested peptides at all
#'   missed-cleavage levels.
#' @param missed_cleavages Integer. Which missed-cleavage level to use for
#'   the overlap count. Defaults to `1L` (standard bottom-up practice).
#'   Pass `0L` for strict non-overlapping fragments, or `NULL` to count
#'   all MC levels (automatically set when `length_range = NULL`).
#' @param residues_per_line Positive integer. Number of residues to show before
#'   wrapping to the next display row. Defaults to `50L`.
#' @param title Optional character string for the plot title. Auto-generated
#'   from the protein accession and enzyme when `NULL` (default).
#'
#' @return A `ggplot` object showing one tile per residue.
#'
#' @details The plotted letters are reconstructed from the MC=0 peptide set in
#' `result$peptides`, so the plot can be generated directly from an
#' [evaluate_digest()] result without re-reading the original FASTA input.
#' Overlap counts are computed from valid-length MC=1 peptides by default
#' (configurable via `missed_cleavages`). Pass `missed_cleavages = 0L` for
#' strict non-overlapping fragments, or `length_range = NULL` to count all
#' digested peptides at all missed-cleavage levels.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin", missed_cleavages = 1)
#'   p <- plot_peptide_overlap_map(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [plot_coverage_map()], [plot_cleavage_map()]
#' @family plot-single
#' @export
plot_peptide_overlap_map <- function(result,
                                     length_range = c(7L, 25L),
                                     missed_cleavages = 1L,
                                     residues_per_line = 50L,
                                     title = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots"
  )

  .validate_digest_result_for_plot(result)

  normalized_length_range <- if (is.null(length_range)) {
    NULL
  } else {
    .validate_length_range(length_range)
  }

  if (
    !is.numeric(residues_per_line) ||
      length(residues_per_line) != 1L ||
      is.na(residues_per_line)
  ) {
    .abort(
      "{.arg residues_per_line} must be a single positive integer.",
      class = "pepvet_error_invalid_input"
    )
  }

  normalized_wrap <- as.integer(residues_per_line)
  if (
    normalized_wrap < 1L ||
      !isTRUE(all.equal(as.numeric(residues_per_line), as.numeric(normalized_wrap)))
  ) {
    .abort(
      "{.arg residues_per_line} must be a single positive integer.",
      class = "pepvet_error_invalid_input"
    )
  }

  peps <- result$peptides
  params <- result$params
  protein_id <- params$protein_ids[[1L]]
  enzyme_name <- params$enzyme
  display_id <- .tidy_protein_id(protein_id)
  protein_length <- max(peps$end, na.rm = TRUE)

  tile_df <- .build_peptide_overlap_df(
    peps,
    protein_length = protein_length,
    length_range = normalized_length_range,
    missed_cleavages = if (is.null(normalized_length_range)) NULL else missed_cleavages,
    residues_per_line = normalized_wrap
  )

  n_detected <- sum(tile_df$overlap_count > 0L)
  pct_detected <- round(100 * n_detected / nrow(tile_df), 1)
  n_overlapped <- sum(tile_df$overlap_count >= 2L)
  max_overlap <- max(tile_df$overlap_count)
  overlap_source <- if (is.null(normalized_length_range)) {
    "all digested peptides"
  } else {
    paste0(
      "valid peptides [",
      normalized_length_range[[1L]],
      "-",
      normalized_length_range[[2L]],
      " aa]"
    )
  }

  auto_title <- if (is.null(title)) {
    paste0(display_id, "    \u00b7    ", enzyme_name)
  } else {
    title
  }

  letter_size <- if (normalized_wrap <= 40L) {
    3.4
  } else if (normalized_wrap <= 60L) {
    3.0
  } else if (normalized_wrap <= 80L) {
    2.4
  } else {
    1.9
  }

  ggplot2::ggplot(
    tile_df,
    ggplot2::aes(x = .data$column_index, y = 1, fill = .data$overlap_class)
  ) +
    ggplot2::geom_tile(
      width = 0.98,
      height = 0.92,
      color = .pepvet_pal$separator,
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$residue, color = .data$letter_color),
      size = letter_size,
      fontface = "bold"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$line_label),
      switch = "y"
    ) +
    ggplot2::scale_fill_manual(
      values = .pepvet_pal$overlap,
      drop = FALSE,
      name = "Peptide overlaps"
    ) +
    ggplot2::scale_color_manual(
      values = c(dark = .pepvet_pal$brand_dark, light = "white"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(1L, normalized_wrap, by = 10L),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = NULL,
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = auto_title,
      subtitle = sprintf(
        paste(
          "%.0f%% residues detected by %s.",
          "%d residues have overlapping support; max overlap = %d."
        ),
        pct_detected,
        overlap_source,
        n_overlapped,
        max_overlap
      ),
      x = "Residue position within wrapped line",
      y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(
        fill = .pepvet_pal$neutral,
        color = .pepvet_pal$separator
      ),
      strip.text.y.left = ggplot2::element_text(
        angle = 0,
        hjust = 0,
        face = "bold",
        color = .pepvet_pal$brand_dark
      ),
      legend.position = "bottom"
    )
}


# ── plot_cleavage_map ─────────────────────────────────────────────────────────

#' Cleavage Site Map
#'
#' `plot_cleavage_map()` draws the full protein as a horizontal bar and marks
#' every cleavage site as a vertical tick, colored by efficiency (high=green,
#' medium=amber, low=red).  Peptide fragments between consecutive cleavage
#' sites are drawn as colored blocks, with invalid peptides dimmed.  When
#' `cleavage_sites` data is not available, sites are inferred from the peptide
#' boundaries and all rendered as the same default color.
#'
#' @param result A named list returned by [evaluate_digest()].
#' @param cleavage_sites Optional data.frame from [annotate_cleavage_sites()]
#'   with columns `position`, `efficiency` (character: `"high"`, `"medium"`,
#'   `"low"`), and optionally `rule`.  When `NULL` (default) sites are inferred
#'   from peptide boundaries.
#' @param length_range Integer vector of length 2.  Valid peptide window.
#'   Defaults to `c(7L, 25L)`.
#' @param title Optional character title.  Auto-generated when `NULL`.
#'
#' @return A `ggplot` object.
#' @seealso [evaluate_digest()], [annotate_cleavage_sites()],
#'   [plot_coverage_map()]
#' @family plot-single
#' @export
plot_cleavage_map <- function(result,
                              cleavage_sites = NULL,
                              length_range = c(7L, 25L),
                              title = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots"
  )

  .validate_digest_result_for_plot(result)

  length_range <- .validate_length_range(length_range)

  peps <- result$peptides
  params <- result$params
  protein_id <- params$protein_ids[[1L]]
  enzyme_name <- params$enzyme
  display_id <- .tidy_protein_id(protein_id)
  protein_length <- max(peps$end, na.rm = TRUE)

  length_lo <- length_range[[1L]]
  length_hi <- length_range[[2L]]
  peps$valid <- peps$length >= length_lo & peps$length <= length_hi

  # Build cleavage site table (infer if not provided)
  if (
    !is.null(cleavage_sites) &&
      all(c("position", "efficiency") %in% names(cleavage_sites))
  ) {
    site_df <- cleavage_sites
    # Map efficiency to colors
    eff_levels <- c(
      "high" = .pepvet_pal$good,
      "medium" = .pepvet_pal$moderate,
      "low" = .pepvet_pal$poor
    )
    site_df$site_color <- ifelse(
      site_df$efficiency %in% names(eff_levels),
      eff_levels[site_df$efficiency],
      .pepvet_pal$brand
    )
    site_df$efficiency <- factor(
      site_df$efficiency,
      levels = c("high", "medium", "low")
    )
    has_efficiency <- TRUE
  } else {
    # Infer sites from peptide end positions (not including C-terminus)
    site_positions <- unique(peps$end[peps$end < protein_length])
    site_df <- data.frame(
      position = sort(site_positions),
      efficiency = "unknown",
      site_color = .pepvet_pal$brand,
      stringsAsFactors = FALSE
    )
    has_efficiency <- FALSE
  }

  # Fragment blocks between cleavage sites (using MC=0 peptides only)
  mc0 <- if ("missed_cleavages" %in% names(peps)) {
    peps[peps$missed_cleavages == 0L, , drop = FALSE]
  } else {
    peps
  }

  # Peptide label: length if valid and wide enough (> 2.5% protein width)
  min_label_width <- ceiling(protein_length * 0.025)
  mc0$label <- ifelse(mc0$valid & mc0$length >= min_label_width,
    as.character(mc0$length), ""
  )

  x_step <- max(50L, as.integer(round(protein_length / 10.0 / 50.0) * 50L))

  p <- ggplot2::ggplot() +
    # Full protein background bar
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = protein_length + 0.5,
      ymin = 0.28, ymax = 0.72,
      fill = .pepvet_pal$cleavage_bg, color = .pepvet_pal$protein_brd, linewidth = 0.4
    ) +
    # Fragment blocks
    ggplot2::geom_rect(
      data = mc0,
      ggplot2::aes(
        xmin = .data$start - 0.3,
        xmax = .data$end + 0.3,
        ymin = 0.32, ymax = 0.68,
        fill = .data$valid,
        alpha = .data$valid
      ),
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = .pepvet_pal$brand, `FALSE` = .pepvet_pal$invalid_pep),
      labels = c(`TRUE` = "Valid", `FALSE` = "Invalid"),
      name = NULL
    ) +
    ggplot2::scale_alpha_manual(
      values = c(`TRUE` = 0.75, `FALSE` = 0.35),
      guide = "none"
    )

  # Fragment length labels
  label_mc0 <- mc0[nchar(mc0$label) > 0L, , drop = FALSE]
  if (nrow(label_mc0) > 0L) {
    label_mc0$mid_x <- (label_mc0$start + label_mc0$end) / 2.0
    p <- p + ggplot2::geom_text(
      data = label_mc0,
      ggplot2::aes(x = .data$mid_x, y = 0.50, label = .data$label),
      size = 2.5, color = "white", fontface = "bold"
    )
  }

  # Cleavage site ticks
  if (nrow(site_df) > 0L) {
    if (has_efficiency) {
      p <- p + ggplot2::geom_segment(
        data = site_df,
        ggplot2::aes(
          x = .data$position + 0.5, xend = .data$position + 0.5,
          y = 0.10, yend = 0.88,
          color = .data$efficiency
        ),
        linewidth = 0.8,
        alpha = 0.85
      ) +
        ggplot2::scale_color_manual(
          values = c(
            high = .pepvet_pal$good,
            medium = .pepvet_pal$moderate,
            low = .pepvet_pal$poor
          ),
          name = "Cleavage efficiency",
          guide = ggplot2::guide_legend(
            override.aes = list(linewidth = 2)
          )
        )
    } else {
      p <- p + ggplot2::geom_segment(
        data = site_df,
        ggplot2::aes(
          x = .data$position + 0.5, xend = .data$position + 0.5,
          y = 0.10, yend = 0.88
        ),
        color = .pepvet_pal$brand_dark,
        linewidth = 1.0,
        alpha = 1.0
      )
    }
  }

  n_valid <- sum(mc0$valid)
  n_total <- nrow(mc0)
  n_sites <- nrow(site_df)
  auto_title <- title %||% paste0(display_id, "    \u00b7    ", enzyme_name)

  p +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, protein_length + x_step, by = x_step),
      limits = c(0, protein_length + 1),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = auto_title,
      subtitle = sprintf(
        "%d cleavage sites  \u00b7  %d / %d valid fragments  \u00b7  protein length %d aa",
        n_sites, n_valid, n_total, protein_length
      ),
      caption = if (!has_efficiency) {
        "Tip: pass annotate_cleavage_sites() output to see efficiency coloring"
      } else {
        NULL
      },
      x = "Residue position",
      y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position    = if (has_efficiency) "bottom" else "none"
    )
}


#' Weight Sensitivity Plot
#'
#' `plot_weight_sensitivity()` visualises the result of
#' [sensitivity_analysis()].  For single-protein input it draws a density plot
#' of composite scores with verdict shading, threshold lines, and a rug mark
#' at the default composite.  For batch input it draws a faceted histogram of
#' per-protein verdict instability with per-enzyme mean and median annotations.
#'
#' @param x A list returned by [sensitivity_analysis()].
#' @param title Optional character string for the plot title.  When `NULL`
#'   (default) a title is auto-generated.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   sens <- sensitivity_analysis(res, n_iter = 1000L)
#'   p <- plot_weight_sensitivity(sens)
#'   print(p)
#' }
#'
#' @seealso [sensitivity_analysis()]
#' @family plot-single
#' @export
plot_weight_sensitivity <- function(x, title = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots"
  )

  if ("per_protein" %in% names(x)) {
    .plot_sensitivity_batch(x, title)
  } else {
    .plot_sensitivity_single(x, title)
  }
}


.plot_sensitivity_single <- function(x, title) {
  iter <- x$iterations
  summ <- x$summary
  default_comp <- summ$composite_mean
  n_iter <- nrow(iter)

  auto_title <- title %||% "Weight Sensitivity Analysis"

  iter_df <- iter
  n_good <- sum(iter_df$verdict == "Good")
  n_mod  <- sum(iter_df$verdict == "Moderate")
  n_poor <- sum(iter_df$verdict == "Poor")

  verdict_pct_str <- paste(
    sprintf("%s (%.0f%%)", names(summ$verdict_pct), summ$verdict_pct * 100),
    collapse = " / "
  )

  fill_scale <- c(
    "Good"     = .pepvet_pal$good,
    "Moderate" = .pepvet_pal$moderate,
    "Poor"     = .pepvet_pal$poor
  )

  p <- ggplot2::ggplot(iter_df, ggplot2::aes(x = .data$composite_score))

  if (n_mod > 0L && n_poor > 0L) {
    p <- p + ggplot2::geom_density(
      ggplot2::aes(fill = .data$verdict, color = .data$verdict),
      alpha = 0.30, bw = 0.025, linewidth = 0.3
    )
  } else if (n_mod > 0L) {
    p <- p +
      ggplot2::geom_density(fill = .pepvet_pal$good, colour = .pepvet_pal$good,
        alpha = 0.30, bw = 0.025, linewidth = 0.3) +
      ggplot2::geom_density(
        data = iter_df[iter_df$verdict == "Moderate", , drop = FALSE],
        fill = .pepvet_pal$moderate, colour = .pepvet_pal$moderate,
        alpha = 0.30, bw = 0.025, linewidth = 0.3
      )
  } else {
    p <- p + ggplot2::geom_density(
      fill = .pepvet_pal$good, color = .pepvet_pal$brand_dark,
      alpha = 0.20, bw = 0.025, linewidth = 0.4
    )
  }

  ci_lo <- summ$composite_ci[[1L]]
  ci_hi <- summ$composite_ci[[2L]]

  p <- p +
    # Verdict zone annotations
    ggplot2::annotate("text", x = 0.20, y = Inf,
      label = "Poor", hjust = 0.5, vjust = 1.8, size = 3,
      colour = .pepvet_pal$poor, fontface = "italic", alpha = 0.6) +
    ggplot2::annotate("text", x = 0.525, y = Inf,
      label = "Moderate", hjust = 0.5, vjust = 1.8, size = 3,
      colour = .pepvet_pal$moderate, fontface = "italic", alpha = 0.6) +
    ggplot2::annotate("text", x = 0.825, y = Inf,
      label = "Good", hjust = 0.5, vjust = 1.8, size = 3,
      colour = .pepvet_pal$good, fontface = "italic", alpha = 0.6) +
    # Threshold lines
    ggplot2::geom_vline(
      xintercept = c(.get_param("verdict_moderate"), .get_param("verdict_good")),
      linetype = "dashed", colour = .pepvet_pal$text_axis_title,
      linewidth = 0.5) +
    # CI ribbon
    ggplot2::annotate("rect",
      xmin = ci_lo, xmax = ci_hi, ymin = 0, ymax = Inf,
      fill = .pepvet_pal$brand, alpha = 0.06) +
    # Default composite rug
    ggplot2::geom_rug(
      data = data.frame(x = default_comp),
      ggplot2::aes(x = .data$x),
      sides = "b", colour = .pepvet_pal$brand_dark, linewidth = 0.8) +
    ggplot2::annotate("text",
      x = default_comp, y = 0,
      label = sprintf("default = %.3f", default_comp),
      hjust = 0, vjust = -0.5, size = 2.8,
      colour = .pepvet_pal$brand_dark, fontface = "bold") +
    ggplot2::scale_fill_manual(values = fill_scale, guide = "none") +
    ggplot2::scale_colour_manual(values = fill_scale, guide = "none") +
    ggplot2::labs(
      title = auto_title,
      subtitle = sprintf(
        "%s  |  95%% CI [%.3f, %.3f]  |  %d iterations",
        verdict_pct_str, ci_lo, ci_hi, n_iter
      ),
      x = "Composite score",
      y = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = 0.05)
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank()
    )
  p
}


# Batch sensitivity histogram
.plot_sensitivity_batch <- function(x, title) {
  pp <- x$per_protein
  total_inst <- x$summary$total_instability
  n_prot <- length(unique(pp$protein_id))
  n_iter_est <- NULL  # not stored in summary, infer from data if possible

  auto_title <- title %||% "Proteome Verdict Stability"

  has_enzyme <- "enzyme" %in% names(pp)

  if (!has_enzyme) {
    median_inst <- stats::median(pp$verdict_instability, na.rm = TRUE)
    pct_stable <- mean(pp$verdict_instability < 0.05, na.rm = TRUE) * 100

    ggplot2::ggplot(pp, ggplot2::aes(x = .data$verdict_instability)) +
      ggplot2::geom_histogram(bins = 30, fill = .pepvet_pal$brand,
        colour = "white", alpha = 0.85, linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = total_inst,
        colour = .pepvet_pal$poor, linewidth = 0.7, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = median_inst,
        colour = .pepvet_pal$moderate, linewidth = 0.5, linetype = "dotted") +
      ggplot2::annotate("text",
        x = total_inst, y = Inf,
        label = sprintf("Mean = %.1f%%", total_inst * 100),
        hjust = if (total_inst > 0.5) 1 else -0.1,
        vjust = 1.5, size = 2.8, colour = .pepvet_pal$poor,
        fontface = "italic") +
      ggplot2::labs(
        title = auto_title,
        subtitle = sprintf(
          "%.0f%% of proteins have <5%% verdict instability  |  Median = %.1f%%  |  %d proteins",
          pct_stable, median_inst * 100, n_prot
        ),
        x = "Verdict instability (fraction of iterations where verdict changed)",
        y = "Count"
      ) +
      .pepvet_theme()
  } else {
    enz_means <- tapply(pp$verdict_instability, pp$enzyme, mean, na.rm = TRUE)
    enz_medians <- tapply(pp$verdict_instability, pp$enzyme, stats::median, na.rm = TRUE)
    enz_n <- tapply(pp$verdict_instability, pp$enzyme, function(v) sum(!is.na(v)))

    enz_order <- names(sort(enz_means))
    pp$enzyme <- factor(pp$enzyme, levels = enz_order)

    ann_df <- data.frame(
      enzyme = names(enz_means),
      mean_inst = as.numeric(enz_means),
      median_inst = as.numeric(enz_medians),
      n_prots = as.integer(enz_n),
      stringsAsFactors = FALSE
    )
    ann_df$enzyme <- factor(ann_df$enzyme, levels = enz_order)
    ann_df$label <- sprintf(
      "mean = %.1f%%\nmedian = %.1f%%\nn = %d",
      ann_df$mean_inst * 100, ann_df$median_inst * 100, ann_df$n_prots
    )

    pct_stable <- mean(pp$verdict_instability < 0.05, na.rm = TRUE) * 100

    # Wrap long enzyme names for facet labels
    enz_levels <- levels(pp$enzyme)
    wrapped <- vapply(enz_levels, function(e) {
      paste(strwrap(e, width = 20), collapse = "\n")
    }, character(1))
    names(wrapped) <- enz_levels
    pp$enzyme_wrapped <- factor(wrapped[as.character(pp$enzyme)],
      levels = wrapped[enz_order])

    ann_df$enzyme_wrapped <- factor(wrapped[as.character(ann_df$enzyme)],
      levels = wrapped[enz_order])

    ggplot2::ggplot(pp, ggplot2::aes(x = .data$verdict_instability)) +
      ggplot2::geom_histogram(bins = 25, fill = .pepvet_pal$brand,
        colour = "white", alpha = 0.75, linewidth = 0.15) +
      ggplot2::geom_vline(data = ann_df,
        ggplot2::aes(xintercept = .data$mean_inst),
        colour = .pepvet_pal$poor, linewidth = 0.5, linetype = "dashed") +
      ggplot2::geom_text(data = ann_df,
        ggplot2::aes(x = Inf, y = Inf, label = .data$label),
        hjust = 1.05, vjust = 1.2, size = 2.5,
        colour = .pepvet_pal$text_subtitle, lineheight = 0.9) +
      ggplot2::facet_wrap(~ .data$enzyme_wrapped, ncol = 5, scales = "free_y") +
      ggplot2::labs(
        title = auto_title,
        subtitle = sprintf(
          "%.0f%% of protein-enzyme pairs have <5%% verdict instability  |  %d proteins  |  %d enzymes",
          pct_stable, n_prot, length(enz_means)
        ),
        x = "Verdict instability (fraction of iterations where verdict changed)",
        y = "Count"
      ) +
      .pepvet_theme() +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 7, face = "bold",
          colour = .pepvet_pal$brand_dark, lineheight = 0.85)
      )
  }
}
