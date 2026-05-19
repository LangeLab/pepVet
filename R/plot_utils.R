# -- pepVet Visualization Suite ----------------------------------------------
#
# All public functions in this file guard their ggplot2 dependency at runtime
# via rlang::check_installed().  This means ggplot2 (and patchwork for
# multi-panel functions) live in Suggests, not Imports, keeping the package
# installable on servers where graphical output is not needed.
#
# Internal helpers (.pepvet_pal, .pepvet_theme, etc.) may be called freely by
# any plotting function in this file.
# ---------------------------------------------------------------------------

#' @importFrom rlang .data check_installed
#' @importFrom stats setNames
#' @importFrom utils head
NULL


# -- Internal color palette --------------------------------------------------

#' pepVet internal color palette
#'
#' A named list of hex color strings used consistently across all pepVet plots.
#' Inspired by the JCO (Journal of Clinical Oncology) palette but adapted for
#' scientific clarity: the brand blue anchors valid peptides and primary data
#' marks, while traffic-light green/amber/red encode scoring tiers.
#'
#' @noRd
.pepvet_pal <- list(
  # Brand / primary
  brand = "#2C5F8A", # pepVet blue - valid peptides, primary bars
  brand_dark = "#1A3D5C", # darker blue - borders, contrast text
  brand_light = "#7BAED4", # lighter blue - fills, highlights

  # Length-class categories
  valid = "#2C5F8A", # valid-length peptides
  too_short = "#E8A838", # amber - too-short peptides
  too_long = "#C94040", # rose-red - too-long peptides

  # Scoring tiers
  good = "#27AE60", # green  - score >= 0.65
  moderate = "#E8A838", # amber  - score 0.40-0.69
  poor = "#C94040", # red    - score < 0.40

  # Coverage
  gap = "#C94040", # red - uncovered protein regions
  covered = "#2C5F8A", # blue - covered by valid peptides

  # Background shading for valid ranges
  shade = "#EDF6F0", # very light green - valid-range backdrop
  neutral = "#F4F6F9", # off-white - panel backgrounds
  separator = "#DDDDDD" # light gray - borders, grid lines
)


# -- Internal ggplot2 theme --------------------------------------------------

#' pepVet base ggplot2 theme
#'
#' Returns a clean, publication-ready ggplot2 theme used by all pepVet
#' visualization functions.  Extends `theme_minimal` with tighter grid lines,
#' branded color accents, and consistent typography.
#'
#' @param base_size Numeric base font size.  Defaults to `11`.
#' @noRd
.pepvet_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # Title / subtitle
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1,
        color = .pepvet_pal$brand_dark,
        margin = ggplot2::margin(b = 3)
      ),
      plot.subtitle = ggplot2::element_text(
        size   = base_size - 1.5,
        color  = "#666666",
        margin = ggplot2::margin(b = 6)
      ),

      # Axes
      axis.title = ggplot2::element_text(size = base_size - 1, color = "#444444"),
      axis.text = ggplot2::element_text(size = base_size - 2, color = "#555555"),
      axis.ticks = ggplot2::element_line(color = "#CCCCCC", linewidth = 0.3),

      # Grid
      panel.grid.major = ggplot2::element_line(color = "#EBEBEB", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill      = NA,
        color     = "#DDDDDD",
        linewidth = 0.5
      ),

      # Strip (for faceted plots)
      strip.background = ggplot2::element_rect(fill = "#F0F4F8", color = "#DDDDDD"),
      strip.text = ggplot2::element_text(
        face  = "bold",
        size  = base_size - 0.5,
        color = "#2C5F8A"
      ),

      # Legend
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.55, "lines"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(
        size = base_size - 1.5,
        face = "bold",
        color = "#444444"
      ),
      legend.background = ggplot2::element_rect(fill = NA, color = NA),

      # Plot margins
      plot.margin = ggplot2::margin(8, 10, 6, 10)
    )
}


# -- Internal score-to-color helper ------------------------------------------

#' Map numeric scores to tier colors
#'
#' @param x Numeric vector of score values between 0 and 1.
#' @return Character vector of hex color strings.
#' @noRd
.score_color <- function(x) {
  colors <- character(length(x))
  colors[x >= .get_param("verdict_good")] <- .pepvet_pal$good
  colors[
    x >= .get_param("verdict_moderate") &
      x < .get_param("verdict_good")
  ] <- .pepvet_pal$moderate
  colors[x < .get_param("verdict_moderate")] <- .pepvet_pal$poor
  colors
}


# -- Internal: tidy protein display ID ---------------------------------------

#' Shorten a FASTA header to an accession + gene label
#'
#' Strips the `sp|ACC|GENE` prefix and returns `"ACC (GENE)"`.  For
#' non-standard headers, falls back to the first 40 characters.
#'
#' @param protein_id Character string - the raw protein ID from the pepVet
#'   peptide table.
#' @noRd
.tidy_protein_id <- function(protein_id) {
  m <- regmatches(
    protein_id,
    regexpr("^[a-z]+\\|([A-Z0-9]+)\\|([A-Z0-9_]+)", protein_id)
  )
  if (length(m) == 1L && nchar(m) > 0L) {
    parts <- strsplit(m, "\\|")[[1L]]
    return(paste0(parts[[2L]], "  (", parts[[3L]], ")"))
  }
  if (nchar(protein_id) > 42L) {
    paste0(substr(protein_id, 1L, 39L), "...")
  } else {
    protein_id
  }
}


# -- Internal: input validation ----------------------------------------------

#' Validate an evaluate_digest result for plotting
#'
#' @param result Object to validate.
#' @noRd
.validate_digest_result_for_plot <- function(result) {
  if (
    !is.list(result) ||
      !all(c("scores", "peptides", "params") %in% names(result))
  ) {
    .abort(
      c(
        "!" = "{.arg result} must be a named list returned by {.fn evaluate_digest}.",
        "i" = "Expected elements {.code scores}, {.code peptides}, and {.code params}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }
  peps <- result$peptides
  required_cols <- c("protein_id", "peptide", "start", "end", "length")
  if (
    !is.data.frame(peps) ||
      !all(required_cols %in% names(peps))
  ) {
    .abort(
      c(
        "!" = "{.code result$peptides} must be a tibble from {.fn digest_protein}.",
        "i" = "Missing required columns: {.val {setdiff(required_cols, names(peps))}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }
  n_proteins <- length(unique(peps$protein_id))
  if (n_proteins > 1L) {
    .abort(
      c(
        "!" = "{.fn plot_digest_profile} is designed for single-protein results.",
        "i" = "Found {.val {n_proteins}} distinct protein IDs in {.code result$peptides}.",
        "i" = "Run {.fn evaluate_digest} with a single-entry FASTA or a bare sequence."
      ),
      class = "pepvet_error_multi_protein"
    )
  }
  invisible(NULL)
}


# -- Panel builders (internal) -----------------------------------------------

#' Panel A - Peptide length distribution
#' @noRd
.panel_length <- function(peps, length_range) {
  length_lo <- length_range[[1L]]
  length_hi <- length_range[[2L]]

  peps$length_class <- factor(
    ifelse(peps$length < length_lo, "Too short",
      ifelse(peps$length > length_hi, "Too long", "Valid")
    ),
    levels = c("Valid", "Too short", "Too long")
  )

  n_valid <- sum(peps$length_class == "Valid")
  n_total <- nrow(peps)
  pct <- round(100 * n_valid / n_total, 1)

  class_colors <- c(
    "Valid"     = .pepvet_pal$valid,
    "Too short" = .pepvet_pal$too_short,
    "Too long"  = .pepvet_pal$too_long
  )

  # Sensible x-axis upper limit
  x_max <- max(peps$length) + 1L

  ggplot2::ggplot(peps, ggplot2::aes(x = length, fill = length_class)) +
    # Valid-range shading
    ggplot2::annotate(
      "rect",
      xmin = length_lo - 0.5, xmax = length_hi + 0.5,
      ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.70
    ) +
    ggplot2::geom_histogram(
      binwidth = 1L,
      color = "white",
      linewidth = 0.2,
      alpha = 0.90
    ) +
    ggplot2::scale_fill_manual(
      values = class_colors,
      name = NULL,
      guide = ggplot2::guide_legend(
        override.aes = list(alpha = 1, color = NA)
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, x_max + 4L, by = 5L)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, x_max + 1)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = "Peptide Length Distribution",
      subtitle = sprintf(
        "%d / %d peptides in valid range [%d\u2013%d aa]  \u00b7  %.0f%% valid",
        n_valid, n_total, length_lo, length_hi, pct
      ),
      x = "Peptide length (aa)",
      y = "Count"
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(
        angle  = 90,
        hjust  = 0.5,
        margin = ggplot2::margin(r = 4)
      )
    )
}


#' Panel B - GRAVY hydrophobicity distribution
#' @noRd
.panel_gravy <- function(peps, gravy_range) {
  gravy_lo <- gravy_range[[1L]]
  gravy_hi <- gravy_range[[2L]]

  n_inside <- sum(peps$gravy >= gravy_lo & peps$gravy <= gravy_hi,
    na.rm = TRUE
  )
  n_total <- sum(!is.na(peps$gravy))

  # Suitable bin count for the data range
  n_bins <- max(15L, min(40L, as.integer(diff(range(peps$gravy, na.rm = TRUE)) / 0.08)))

  ggplot2::ggplot(peps, ggplot2::aes(x = gravy)) +
    # LC-friendly range shading
    ggplot2::annotate(
      "rect",
      xmin = gravy_lo, xmax = gravy_hi,
      ymin = -Inf, ymax = Inf,
      fill = .pepvet_pal$shade, alpha = 0.70
    ) +
    ggplot2::geom_histogram(
      bins = n_bins,
      fill = .pepvet_pal$brand,
      color = "white",
      linewidth = 0.2,
      alpha = .get_param("scatter_alpha")
    ) +
    # Range boundary lines
    ggplot2::geom_vline(
      xintercept = gravy_lo,
      linetype   = "dashed",
      color      = .pepvet_pal$good,
      linewidth  = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = gravy_hi,
      linetype   = "dashed",
      color      = .pepvet_pal$poor,
      linewidth  = 0.7
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title = "GRAVY Hydrophobicity",
      subtitle = sprintf(
        "LC-friendly range [%.1f, %.1f]  \u00b7  %d / %d peptides inside",
        gravy_lo, gravy_hi, n_inside, n_total
      ),
      x = "GRAVY score",
      y = NULL
    ) +
    .pepvet_theme()
}


#' Compute coverage statistics from a peptide table (shared helper)
#'
#' Splits peptides into valid/invalid sets, computes per-residue coverage,
#' summary percentage, and gap runs.  Used by both `.panel_coverage()` and
#' `plot_coverage_map()` to avoid duplicating this logic.
#'
#' @param peps      Peptide tibble from [evaluate_digest()].
#' @param protein_length Integer. Length of the protein in amino acids.
#' @param length_range   Integer vector of length 2.
#' @param mc_filter  Optional integer.  When not `NULL`, coverage is computed
#'   only on peptides with `missed_cleavages == mc_filter`.
#' @return A named list: `valid_peps`, `invalid_peps`, `covered` (logical
#'   vector of length `protein_length`), `pct_cov` (numeric), `gap_df`
#'   (data.frame with `xmin`/`xmax` for each uncovered run).
#' @noRd
.compute_coverage_stats <- function(peps, protein_length,
                                    length_range, mc_filter = NULL) {
  length_lo <- length_range[[1L]]
  length_hi <- length_range[[2L]]
  valid_peps <- peps[peps$length >= length_lo & peps$length <= length_hi, ,
    drop = FALSE
  ]
  invalid_peps <- peps[peps$length < length_lo | peps$length > length_hi, ,
    drop = FALSE
  ]

  # Restrict coverage calculation to a single MC level when requested
  cov_peps <- if (!is.null(mc_filter) && "missed_cleavages" %in% names(peps)) {
    valid_peps[valid_peps$missed_cleavages == mc_filter, , drop = FALSE]
  } else {
    valid_peps
  }

  covered <- rep(FALSE, protein_length)
  for (i in seq_len(nrow(cov_peps))) {
    s <- cov_peps$start[[i]]
    e <- min(cov_peps$end[[i]], protein_length)
    covered[seq.int(s, e)] <- TRUE
  }
  pct_cov <- round(100 * sum(covered) / protein_length, 1L)

  rl <- rle(covered)
  rl_ends <- cumsum(rl$lengths)
  rl_starts <- c(1L, head(rl_ends, -1L) + 1L)
  gap_mask <- !rl$values

  list(
    valid_peps = valid_peps,
    invalid_peps = invalid_peps,
    covered = covered,
    pct_cov = pct_cov,
    gap_df = data.frame(
      xmin = rl_starts[gap_mask],
      xmax = rl_ends[gap_mask]
    )
  )
}


#' Map GRAVY values to hex colors via a 4-stop gradient (shared helper)
#'
#' Color stops follow the pepVet physicochemical story, reusing `.pepvet_pal`:
#' very hydrophilic (GRAVY << -1) maps to brand blue; neutral (about 0)
#' maps to good green; borderline (about 0.6) maps to amber; very
#' hydrophobic (>> 0.6) maps to poor red.
#'
#' @param gravy_values Numeric vector.
#' @param lo,hi       Clamp limits. Values outside the lo-to-hi interval are clamped.
#' @return Character vector of hex color strings, same length as input.
#' @noRd
.gravy_to_color <- function(gravy_values, lo = -2.0, hi = 2.0) {
  stops <- c(
    .pepvet_pal$brand, .pepvet_pal$good,
    .pepvet_pal$moderate, .pepvet_pal$poor
  )
  ramp <- grDevices::colorRamp(stops, interpolate = "spline")
  scaled <- pmax(0, pmin(1, (gravy_values - lo) / (hi - lo)))
  m <- ramp(scaled)
  grDevices::rgb(m[, 1L], m[, 2L], m[, 3L], maxColorValue = 255)
}


#' Greedy interval packing for non-overlapping peptide display (shared helper)
#'
#' Assigns each peptide to the lowest-numbered "track" (sub-row) where it does
#' not overlap any previously placed peptide.  Peptides are processed in order
#' of start position.  The result is a new column `track` (1-based integer).
#' For MC=0 tryptic peptides the result is always track=1 (no overlaps exist).
#' For MC>=1 peptides, consecutive merged peptides share residues, so the
#' packer distributes them across 2-3 tracks.
#'
#' @param peps Data.frame with columns `start` and `end` (integer).
#' @return `peps` with an added integer column `track`.
#' @noRd
.pack_peptides <- function(peps) {
  if (nrow(peps) == 0L) {
    peps$track <- integer(0L)
    return(peps)
  }
  o <- order(peps$start, peps$end)
  peps <- peps[o, , drop = FALSE]
  n <- nrow(peps)

  track_ends <- integer(0L) # last 'end' position in each open track
  tracks <- integer(n)

  for (i in seq_len(n)) {
    # Find tracks where the last end is strictly before this peptide's start
    fit <- which(track_ends < peps$start[[i]])
    if (length(fit) == 0L) {
      track_ends <- c(track_ends, peps$end[[i]])
      tracks[[i]] <- length(track_ends)
    } else {
      t <- min(fit)
      track_ends[t] <- peps$end[[i]]
      tracks[[i]] <- t
    }
  }
  peps$track <- tracks
  peps
}


#' Compute y-band coordinates for multi-lane coverage plots (shared helper)
#'
#' Divides the plotting area (y between 0 and 1) into equal horizontal lanes, one
#' per missed-cleavage level, optionally reserving space at the bottom for
#' cleavage-site tick marks.
#'
#' @param mc_levels    Integer vector of MC levels present (e.g. `0:2`).
#' @param tick_height  Fraction of total height reserved for ticks.  0 = none.
#' @return A data.frame with columns `mc`, `y_lo`, `y_mid`, `y_hi`.
#' @noRd
.lane_y_coords <- function(mc_levels, tick_height = 0.0) {
  n <- length(mc_levels)
  gap <- 0.03
  avail <- 1.0 - tick_height - gap * max(0L, n - 1L)
  lane_h <- avail / n

  rows <- lapply(seq_along(mc_levels), function(i) {
    y_lo <- tick_height + (i - 1L) * (lane_h + gap)
    y_hi <- y_lo + lane_h
    data.frame(
      mc = mc_levels[[i]],
      y_lo = y_lo,
      y_mid = (y_lo + y_hi) / 2.0,
      y_hi = y_hi,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' Panel C - Sequence coverage map
#' @noRd
.panel_coverage <- function(peps, protein_length, length_range) {
  # Delegate statistics to the shared helper (MC=0 only for coverage %)
  cs <- .compute_coverage_stats(peps, protein_length, length_range,
    mc_filter = 0L
  )
  valid_peps <- cs$valid_peps
  invalid_peps <- cs$invalid_peps
  pct_cov <- cs$pct_cov
  gap_df <- cs$gap_df

  # Label peptides >= 8 aa (enough room for a 2-digit number)
  label_peps <- valid_peps[valid_peps$length >= 8L, , drop = FALSE]
  label_peps$label_x <- (label_peps$start + label_peps$end) / 2.0

  # x-axis break step
  x_step <- max(50L, as.integer(round(protein_length / 10.0 / 50.0) * 50L))

  p <- ggplot2::ggplot() +
    # Full protein background bar
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = protein_length + 0.5,
      ymin = 0.30, ymax = 0.70,
      fill = "#D8DDE6",
      color = "#AAAAAA",
      linewidth = 0.4
    )

  # Invalid peptides (dimmed, behind valid ones)
  if (nrow(invalid_peps) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = invalid_peps,
      ggplot2::aes(
        xmin = .data$start - 0.3,
        xmax = .data$end + 0.3,
        ymin = 0.35, ymax = 0.65
      ),
      fill = "#C5CDD8",
      color = "white",
      linewidth = 0.15,
      alpha = 0.60
    )
  }

  # Valid peptide segments
  if (nrow(valid_peps) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = valid_peps,
      ggplot2::aes(
        xmin = .data$start - 0.4,
        xmax = .data$end + 0.4,
        ymin = 0.22, ymax = 0.78
      ),
      fill = .pepvet_pal$covered,
      color = "white",
      linewidth = 0.15,
      alpha = .get_param("scatter_alpha")
    )
  }

  # Gap overlays (red tint on uncovered regions)
  if (nrow(gap_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = gap_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0.30, ymax = 0.70),
      fill = .pepvet_pal$gap,
      alpha = 0.22
    )
  }

  # Peptide length labels inside valid segments
  if (nrow(label_peps) > 0L) {
    p <- p + ggplot2::geom_text(
      data = label_peps,
      ggplot2::aes(x = .data$label_x, y = 0.50, label = .data$length),
      size = 2.3,
      color = "white",
      fontface = "bold"
    )
  }

  p +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, protein_length + x_step, by = x_step),
      limits = c(0, protein_length + 1),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.05, 0.05))
    ) +
    ggplot2::labs(
      title = "Sequence Coverage",
      subtitle = sprintf(
        paste(
          "%.0f%% covered by valid-length peptides.",
          "%d uncovered region(s). Protein length %d aa"
        ),
        pct_cov, nrow(gap_df), protein_length
      ),
      x = "Residue position",
      y = NULL
    ) +
    .pepvet_theme() +
    ggplot2::theme(
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )
}


#' Panel D - Component score bar chart
#' @noRd
.panel_scores <- function(scores) {
  score_cols <- c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge")
  if ("S_unique" %in% names(scores)) {
    score_cols <- c(score_cols, "S_unique")
  }
  score_cols <- score_cols[score_cols %in% names(scores)]

  score_labels <- c(
    S_length   = "Length",
    S_coverage = "Coverage",
    S_count    = "Count",
    S_hydro    = "Hydrophobicity",
    S_charge   = "Charge richness",
    S_unique   = "Uniqueness"
  )

  vals <- as.numeric(scores[1L, score_cols])

  df <- data.frame(
    score = score_cols,
    label = score_labels[score_cols],
    value = vals,
    stringsAsFactors = FALSE
  )
  good_thresh <- .get_param("verdict_good")
  mod_thresh <- .get_param("verdict_moderate")
  df$tier <- ifelse(df$value >= good_thresh, "Good",
    ifelse(df$value >= mod_thresh, "Moderate", "Poor")
  )
  # Ordered factor so highest-priority score is at top of chart
  df$label <- factor(df$label, levels = rev(df$label))

  tier_colors <- c(
    Good     = .pepvet_pal$good,
    Moderate = .pepvet_pal$moderate,
    Poor     = .pepvet_pal$poor
  )

  composite <- as.numeric(scores$composite_score[[1L]])
  verdict <- as.character(scores$verdict[[1L]])

  # Verdict badge background color
  badge_fill <- switch(verdict,
    Good     = .pepvet_pal$good,
    Moderate = .pepvet_pal$moderate,
    Poor     = .pepvet_pal$poor,
    .pepvet_pal$neutral
  )

  n_scores <- nrow(df)

  ggplot2::ggplot(df, ggplot2::aes(x = value, y = label, fill = tier)) +
    # Tier boundary guide lines
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_moderate"),
      linetype   = "dotted",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.55,
      alpha      = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = .get_param("verdict_good"),
      linetype   = "dotted",
      color      = .pepvet_pal$good,
      linewidth  = 0.55,
      alpha      = 0.8
    ) +
    ggplot2::geom_col(width = 0.62, alpha = 0.90) +
    # Score value labels (outside bars)
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.3f", value),
        x     = value + 0.015
      ),
      hjust = 0,
      size = 3.1,
      color = "#333333",
      fontface = "bold"
    ) +
    # Composite score reference line
    ggplot2::geom_vline(
      xintercept = composite,
      linetype   = "dashed",
      color      = .pepvet_pal$brand_dark,
      linewidth  = 1.0
    ) +
    # Composite label (annotated at the top of the chart)
    ggplot2::annotate(
      "label",
      x         = composite,
      y         = n_scores + 0.60,
      label     = sprintf("Composite  %.3f\n%s", composite, verdict),
      size      = 3.0,
      hjust     = 0.5,
      fontface  = "bold",
      fill      = badge_fill,
      color     = "white",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(values = tier_colors, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.18),
      breaks = seq(0, 1, by = 0.2),
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Component Scores",
      subtitle = paste0(
        "Dotted thresholds at ", .get_param("verdict_moderate"),
        " (Moderate) and ", .get_param("verdict_good"),
        " (Good). Dashed line = composite"
      ),
      x = "Score (0 \u2013 1)",
      y = NULL
    ) +
    .pepvet_theme()
}


# ===========================================================================
# Public API
# ===========================================================================
