# ── pepVet Visualization Suite ───────────────────────────────────────────────
#
# All public functions in this file guard their ggplot2 dependency at runtime
# via rlang::check_installed().  This means ggplot2 (and patchwork for
# multi-panel functions) live in Suggests, not Imports, keeping the package
# installable on servers where graphical output is not needed.
#
# Internal helpers (.pepvet_pal, .pepvet_theme, etc.) may be called freely by
# any plotting function in this file.
# ─────────────────────────────────────────────────────────────────────────────


# ── Internal color palette ───────────────────────────────────────────────────

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
  brand      = "#2C5F8A",   # pepVet blue – valid peptides, primary bars
  brand_dark = "#1A3D5C",   # darker blue – borders, contrast text
  brand_light = "#7BAED4",  # lighter blue – fills, highlights

  # Length-class categories
  valid      = "#2C5F8A",   # valid-length peptides
  too_short  = "#E8A838",   # amber – too-short peptides
  too_long   = "#C94040",   # rose-red – too-long peptides

  # Scoring tiers
  good       = "#27AE60",   # green  – score >= 0.70
  moderate   = "#E8A838",   # amber  – score  0.40–0.69
  poor       = "#C94040",   # red    – score <  0.40

  # Coverage
  gap        = "#C94040",   # red – uncovered protein regions
  covered    = "#2C5F8A",   # blue – covered by valid peptides

  # Background shading for valid ranges
  shade      = "#EDF6F0",   # very light green – valid-range backdrop
  neutral    = "#F4F6F9",   # off-white – panel backgrounds
  separator  = "#DDDDDD"    # light gray – borders, grid lines
)


# ── Internal ggplot2 theme ────────────────────────────────────────────────────

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
        face  = "bold",
        size  = base_size + 1,
        color = "#1A1A2E",
        margin = ggplot2::margin(b = 3)
      ),
      plot.subtitle = ggplot2::element_text(
        size   = base_size - 1.5,
        color  = "#666666",
        margin = ggplot2::margin(b = 6)
      ),

      # Axes
      axis.title   = ggplot2::element_text(size = base_size - 1, color = "#444444"),
      axis.text    = ggplot2::element_text(size = base_size - 2, color = "#555555"),
      axis.ticks   = ggplot2::element_line(color = "#CCCCCC", linewidth = 0.3),

      # Grid
      panel.grid.major = ggplot2::element_line(color = "#EBEBEB", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(
        fill      = NA,
        color     = "#DDDDDD",
        linewidth = 0.5
      ),

      # Strip (for faceted plots)
      strip.background = ggplot2::element_rect(fill = "#F0F4F8", color = "#DDDDDD"),
      strip.text       = ggplot2::element_text(
        face  = "bold",
        size  = base_size - 0.5,
        color = "#2C5F8A"
      ),

      # Legend
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.55, "lines"),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.title     = ggplot2::element_text(
        size = base_size - 1.5,
        face = "bold",
        color = "#444444"
      ),
      legend.background = ggplot2::element_rect(fill = NA, color = NA),

      # Plot margins
      plot.margin = ggplot2::margin(8, 10, 6, 10)
    )
}


# ── Internal score-to-color helper ───────────────────────────────────────────

#' Map numeric scores to tier colors
#'
#' @param x Numeric vector of score values in [0, 1].
#' @return Character vector of hex color strings.
#' @noRd
.score_color <- function(x) {
  colors <- character(length(x))
  colors[x >= 0.7]              <- .pepvet_pal$good
  colors[x >= 0.4 & x < 0.7]   <- .pepvet_pal$moderate
  colors[x < 0.4]               <- .pepvet_pal$poor
  colors
}


# ── Internal: tidy protein display ID ────────────────────────────────────────

#' Shorten a FASTA header to an accession + gene label
#'
#' Strips the `sp|ACC|GENE` prefix and returns `"ACC (GENE)"`.  For
#' non-standard headers, falls back to the first 40 characters.
#'
#' @param protein_id Character string – the raw protein ID from the pepVet
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
    paste0(substr(protein_id, 1L, 39L), "\u2026")
  } else {
    protein_id
  }
}


# ── Internal: input validation ────────────────────────────────────────────────

#' Validate an evaluate_digest result for plotting
#'
#' @param result Object to validate.
#' @noRd
.validate_digest_result_for_plot <- function(result) {
  if (!is.list(result) ||
      !all(c("scores", "peptides", "params") %in% names(result))) {
    cli::cli_abort(
      c(
        "!" = "{.arg result} must be a named list returned by {.fn evaluate_digest}.",
        "i" = "Expected elements {.code scores}, {.code peptides}, and {.code params}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }
  peps <- result$peptides
  required_cols <- c("protein_id", "peptide", "start", "end", "length")
  if (!is.data.frame(peps) ||
      !all(required_cols %in% names(peps))) {
    cli::cli_abort(
      c(
        "!" = "{.code result$peptides} must be a tibble from {.fn digest_protein}.",
        "i" = "Missing required columns: {.val {setdiff(required_cols, names(peps))}}."
      ),
      class = "pepvet_error_invalid_digest_result"
    )
  }
  n_proteins <- length(unique(peps$protein_id))
  if (n_proteins > 1L) {
    cli::cli_abort(
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


# ── Panel builders (internal) ─────────────────────────────────────────────────

#' Panel A – Peptide length distribution
#' @noRd
.panel_length <- function(peps, length_range) {
  length_lo <- length_range[[1L]]
  length_hi <- length_range[[2L]]

  peps$length_class <- factor(
    ifelse(peps$length < length_lo, "Too short",
      ifelse(peps$length > length_hi, "Too long", "Valid")),
    levels = c("Valid", "Too short", "Too long")
  )

  n_valid <- sum(peps$length_class == "Valid")
  n_total <- nrow(peps)
  pct     <- round(100 * n_valid / n_total, 1)

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
      color    = "white",
      linewidth = 0.2,
      alpha    = 0.90
    ) +
    ggplot2::scale_fill_manual(
      values = class_colors,
      name   = NULL,
      guide  = ggplot2::guide_legend(
        override.aes = list(alpha = 1, color = NA)
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0L, x_max + 4L, by = 5L)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, x_max + 1)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(
      title    = "Peptide Length Distribution",
      subtitle = sprintf(
        "%d / %d peptides in valid range [%d\u2013%d aa]  \u00b7  %.0f%% valid",
        n_valid, n_total, length_lo, length_hi, pct
      ),
      x = "Peptide length (aa)",
      y = "Count"
    ) +
    .pepvet_theme()
}


#' Panel B – GRAVY hydrophobicity distribution
#' @noRd
.panel_gravy <- function(peps, gravy_range) {
  gravy_lo <- gravy_range[[1L]]
  gravy_hi <- gravy_range[[2L]]

  n_inside <- sum(peps$gravy >= gravy_lo & peps$gravy <= gravy_hi,
                  na.rm = TRUE)
  n_total  <- sum(!is.na(peps$gravy))

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
      bins      = n_bins,
      fill      = .pepvet_pal$brand,
      color     = "white",
      linewidth = 0.2,
      alpha     = 0.88
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
      title    = "GRAVY Hydrophobicity",
      subtitle = sprintf(
        "LC-friendly range [%.1f, %.1f]  \u00b7  %d / %d peptides inside",
        gravy_lo, gravy_hi, n_inside, n_total
      ),
      x = "GRAVY score",
      y = "Count"
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
  length_lo    <- length_range[[1L]]
  length_hi    <- length_range[[2L]]
  valid_peps   <- peps[peps$length >= length_lo & peps$length <= length_hi, ,
                       drop = FALSE]
  invalid_peps <- peps[peps$length  < length_lo | peps$length  > length_hi, ,
                       drop = FALSE]

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

  rl        <- rle(covered)
  rl_ends   <- cumsum(rl$lengths)
  rl_starts <- c(1L, head(rl_ends, -1L) + 1L)
  gap_mask  <- !rl$values

  list(
    valid_peps   = valid_peps,
    invalid_peps = invalid_peps,
    covered      = covered,
    pct_cov      = pct_cov,
    gap_df       = data.frame(
      xmin = rl_starts[gap_mask],
      xmax = rl_ends[gap_mask]
    )
  )
}


#' Map GRAVY values to hex colors via a 4-stop gradient (shared helper)
#'
#' Color stops follow the pepVet physicochemical story, reusing `.pepvet_pal`:
#' very hydrophilic (GRAVY << -1) → brand blue; neutral (≈ 0) → good green;
#' borderline (≈ 0.6) → amber; very hydrophobic (>> 0.6) → poor red.
#'
#' @param gravy_values Numeric vector.
#' @param lo,hi       Clamp limits.  Values outside [lo, hi] are clamped.
#' @return Character vector of hex color strings, same length as input.
#' @noRd
.gravy_to_color <- function(gravy_values, lo = -2.0, hi = 2.0) {
  stops  <- c(.pepvet_pal$brand, .pepvet_pal$good,
               .pepvet_pal$moderate, .pepvet_pal$poor)
  ramp   <- grDevices::colorRamp(stops, interpolate = "spline")
  scaled <- pmax(0, pmin(1, (gravy_values - lo) / (hi - lo)))
  m      <- ramp(scaled)
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
  o    <- order(peps$start, peps$end)
  peps <- peps[o, , drop = FALSE]
  n    <- nrow(peps)

  track_ends <- integer(0L)   # last 'end' position in each open track
  tracks     <- integer(n)

  for (i in seq_len(n)) {
    # Find tracks where the last end is strictly before this peptide's start
    fit <- which(track_ends < peps$start[[i]])
    if (length(fit) == 0L) {
      track_ends    <- c(track_ends, peps$end[[i]])
      tracks[[i]]   <- length(track_ends)
    } else {
      t             <- min(fit)
      track_ends[t] <- peps$end[[i]]
      tracks[[i]]   <- t
    }
  }
  peps$track <- tracks
  peps
}


#' Compute y-band coordinates for multi-lane coverage plots (shared helper)
#'
#' Divides the plotting area (y in [0, 1]) into equal horizontal lanes, one
#' per missed-cleavage level, optionally reserving space at the bottom for
#' cleavage-site tick marks.
#'
#' @param mc_levels    Integer vector of MC levels present (e.g. `0:2`).
#' @param tick_height  Fraction of total height reserved for ticks.  0 = none.
#' @return A data.frame with columns `mc`, `y_lo`, `y_mid`, `y_hi`.
#' @noRd
.lane_y_coords <- function(mc_levels, tick_height = 0.0) {
  n      <- length(mc_levels)
  gap    <- 0.03
  avail  <- 1.0 - tick_height - gap * max(0L, n - 1L)
  lane_h <- avail / n

  rows <- lapply(seq_along(mc_levels), function(i) {
    y_lo <- tick_height + (i - 1L) * (lane_h + gap)
    y_hi <- y_lo + lane_h
    data.frame(mc    = mc_levels[[i]],
               y_lo  = y_lo,
               y_mid = (y_lo + y_hi) / 2.0,
               y_hi  = y_hi,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


#' Panel C – Sequence coverage map
#' @noRd
.panel_coverage <- function(peps, protein_length, length_range) {
  # Delegate statistics to the shared helper (MC=0 only for coverage %)
  cs <- .compute_coverage_stats(peps, protein_length, length_range,
                                 mc_filter = 0L)
  valid_peps   <- cs$valid_peps
  invalid_peps <- cs$invalid_peps
  pct_cov      <- cs$pct_cov
  gap_df       <- cs$gap_df

  # Label peptides >= 8 aa (enough room for a 2-digit number)
  label_peps         <- valid_peps[valid_peps$length >= 8L, , drop = FALSE]
  label_peps$label_x <- (label_peps$start + label_peps$end) / 2.0

  # x-axis break step
  x_step <- max(50L, as.integer(round(protein_length / 10.0 / 50.0) * 50L))

  p <- ggplot2::ggplot() +
    # Full protein background bar
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = protein_length + 0.5,
      ymin = 0.30, ymax = 0.70,
      fill  = "#D8DDE6",
      color = "#AAAAAA",
      linewidth = 0.4
    )

  # Invalid peptides (dimmed, behind valid ones)
  if (nrow(invalid_peps) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = invalid_peps,
      ggplot2::aes(
        xmin = .data$start - 0.3,
        xmax = .data$end   + 0.3,
        ymin = 0.35, ymax = 0.65
      ),
      fill      = "#C5CDD8",
      color     = "white",
      linewidth = 0.15,
      alpha     = 0.60
    )
  }

  # Valid peptide segments
  if (nrow(valid_peps) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = valid_peps,
      ggplot2::aes(
        xmin = .data$start - 0.4,
        xmax = .data$end   + 0.4,
        ymin = 0.22, ymax = 0.78
      ),
      fill      = .pepvet_pal$covered,
      color     = "white",
      linewidth = 0.15,
      alpha     = 0.88
    )
  }

  # Gap overlays (red tint on uncovered regions)
  if (nrow(gap_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = gap_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0.30, ymax = 0.70),
      fill  = .pepvet_pal$gap,
      alpha = 0.22
    )
  }

  # Peptide length labels inside valid segments
  if (nrow(label_peps) > 0L) {
    p <- p + ggplot2::geom_text(
      data = label_peps,
      ggplot2::aes(x = .data$label_x, y = 0.50, label = .data$length),
      size      = 2.3,
      color     = "white",
      fontface  = "bold"
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
      title    = "Sequence Coverage",
      subtitle = sprintf(
        "%.0f%% covered by valid-length peptides  \u00b7  %d uncovered region(s)  \u00b7  protein length %d aa",
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


#' Panel D – Component score bar chart
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
  df$tier  <- ifelse(df$value >= 0.7, "Good",
                ifelse(df$value >= 0.4, "Moderate", "Poor"))
  # Ordered factor so highest-priority score is at top of chart
  df$label <- factor(df$label, levels = rev(df$label))

  tier_colors <- c(
    Good     = .pepvet_pal$good,
    Moderate = .pepvet_pal$moderate,
    Poor     = .pepvet_pal$poor
  )

  composite <- as.numeric(scores$composite_score[[1L]])
  verdict   <- as.character(scores$verdict[[1L]])

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
      xintercept = 0.4,
      linetype   = "dotted",
      color      = .pepvet_pal$moderate,
      linewidth  = 0.55,
      alpha      = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = 0.7,
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
      hjust    = 0,
      size     = 3.1,
      color    = "#333333",
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
      title    = "Component Scores",
      subtitle = "Dotted thresholds at 0.40 (Moderate) and 0.70 (Good)  \u00b7  Dashed line = composite",
      x        = "Score (0 \u2013 1)",
      y        = NULL
    ) +
    .pepvet_theme()
}


# ═══════════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════════

#' Four-Panel Digest Diagnostic Plot
#'
#' `plot_digest_profile()` assembles a four-panel figure for a single
#' protein–enzyme pair from an [evaluate_digest()] result.  The panels are:
#'
#' * **(A) Length distribution** — histogram of peptide lengths with the valid
#'   window shaded.  Bars are colored by length class: valid (blue), too short
#'   (amber), too long (rose).
#' * **(B) GRAVY distribution** — histogram of GRAVY hydrophobicity scores.
#'   The LC-friendly range is shaded and bounded by dashed lines.
#' * **(C) Coverage map** — protein drawn as a horizontal track with
#'   valid-length peptides overlaid as colored segments.  Uncovered regions
#'   are highlighted in red.  Peptide length labels appear inside segments of
#'   8 aa or longer.
#' * **(D) Component scores** — horizontal bar chart for each scoring
#'   component, colored by tier (green \eqn{\geq} 0.70, amber 0.40–0.69, red
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
#'   (default) a title is auto-generated from the protein accession and enzyme.
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
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   p   <- plot_digest_profile(res)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [plot_coverage_map()],
#'   [plot_enzyme_comparison()]
#' @export
plot_digest_profile <- function(result,
                                length_range = c(7L, 25L),
                                gravy_range  = c(-1.0, 0.6),
                                title        = NULL) {
  rlang::check_installed(
    "ggplot2",
    reason = "to produce pepVet visualization plots"
  )
  rlang::check_installed(
    "patchwork",
    reason = "to assemble multi-panel figures in plot_digest_profile()"
  )

  .validate_digest_result_for_plot(result)

  peps   <- result$peptides
  scores <- result$scores
  params <- result$params

  protein_id  <- params$protein_ids[[1L]]
  enzyme_name <- params$enzyme
  display_id  <- .tidy_protein_id(protein_id)

  # Add GRAVY scores to the peptide table (computed from sequences)
  peps$gravy <- vapply(peps$peptide, .calculate_gravy, numeric(1L))

  protein_length <- max(peps$end)

  # ── Build panels ──────────────────────────────────────────────────────────
  pa <- .panel_length(peps,   length_range)
  pb <- .panel_gravy(peps,    gravy_range)
  pc <- .panel_coverage(peps, protein_length, length_range)
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
    paste0(display_id, "    \u00b7    ", enzyme_name)
  } else {
    title
  }

  figure +
    patchwork::plot_annotation(
      title      = auto_title,
      tag_levels = "A",
      theme      = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face   = "bold",
          size   = 15,
          color  = "#1A1A2E",
          margin = ggplot2::margin(b = 10)
        ),
        plot.tag = ggplot2::element_text(
          face   = "bold",
          size   = 14,
          color  = "#2C5F8A",
          margin = ggplot2::margin(t = 1, r = 6, b = 1, l = 2)
        )
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
#'     \item{`"hydrophobicity"`}{Continuous GRAVY gradient for every peptide —
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
#'   cs  <- annotate_cleavage_sites(bsa_path, enzyme = "trypsin")
#'   p   <- plot_coverage_map(res, cleavage_sites = cs)
#'   print(p)
#' }
#'
#' @seealso [evaluate_digest()], [annotate_cleavage_sites()],
#'   [plot_digest_profile()]
#' @export
plot_coverage_map <- function(result,
                              color_by       = c("validity",
                                                 "length_class",
                                                 "hydrophobicity"),
                              length_range   = c(7L, 25L),
                              cleavage_sites = NULL,
                              domains        = NULL,
                              title          = NULL) {
  rlang::check_installed("ggplot2",
    reason = "to produce pepVet visualization plots")

  color_by <- match.arg(color_by)

  .validate_digest_result_for_plot(result)

  # ── Input validation for optional args ─────────────────────────────────────
  if (!is.null(cleavage_sites)) {
    if (!is.data.frame(cleavage_sites) ||
        !all(c("position", "efficiency") %in% names(cleavage_sites))) {
      cli::cli_abort(
        c("!" = paste0("{.arg cleavage_sites} must be a data.frame from ",
                       "{.fn annotate_cleavage_sites}."),
          "i" = "Required columns: {.code position}, {.code efficiency}."),
        class = "pepvet_error_invalid_cleavage_sites"
      )
    }
  }
  if (!is.null(domains)) {
    if (!is.data.frame(domains) ||
        !all(c("name", "start", "end") %in% names(domains))) {
      cli::cli_abort(
        c("!" = paste0("{.arg domains} must be a data.frame with columns ",
                       "{.code name}, {.code start}, {.code end}."),
          "i" = "Each row describes one annotated protein domain."),
        class = "pepvet_error_invalid_domains"
      )
    }
  }

  # ── Extract data ────────────────────────────────────────────────────────────
  peps           <- result$peptides
  params         <- result$params
  protein_id     <- params$protein_ids[[1L]]
  enzyme_name    <- params$enzyme
  display_id     <- .tidy_protein_id(protein_id)
  protein_length <- max(peps$end)

  # ── GRAVY (computed once, used for coloring or silently ignored) ────────────
  peps$gravy <- vapply(peps$peptide, .calculate_gravy, numeric(1L))

  # ── Determine MC levels present and build lane coordinates ─────────────────
  has_mc     <- "missed_cleavages" %in% names(peps)
  mc_levels  <- if (has_mc) sort(unique(peps$missed_cleavages)) else 0L
  tick_h     <- if (!is.null(cleavage_sites)) 0.11 else 0.0
  lanes      <- .lane_y_coords(mc_levels, tick_height = tick_h)

  # ── Coverage stats on MC = 0 peptides (for subtitle) ───────────────────────
  cs0 <- .compute_coverage_stats(peps, protein_length, length_range,
                                  mc_filter = if (has_mc) 0L else NULL)
  pct_cov <- cs0$pct_cov
  gap_df  <- cs0$gap_df

  # ── Build fill aesthetic and scale ─────────────────────────────────────────
  #   For "validity"/"length_class": factor column `fill_cat` → scale_fill_manual
  #   For "hydrophobicity":          numeric column `fill_val` (GRAVY)
  #                                  → scale_fill_gradientn
  if (color_by == "hydrophobicity") {
    peps$fill_val <- peps$gravy
  } else {
    # Validity categories (used by both "validity" and "length_class")
    lo <- length_range[[1L]];  hi <- length_range[[2L]]
    peps$fill_cat <- factor(
      ifelse(peps$length < lo,  "Too short",
        ifelse(peps$length > hi, "Too long", "Valid")),
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
  domain_fills <- c("#D9EAF7", "#FFF0C2", "#D4EDD4", "#F7D9D9",
                    "#EDD4F7", "#D4F7F0", "#F7ECD4", "#D4D9F7")
  if (!is.null(domains)) {
    n_dom <- nrow(domains)
    for (k in seq_len(n_dom)) {
      dom_fill <- domain_fills[((k - 1L) %% length(domain_fills)) + 1L]
      p <- p +
        ggplot2::annotate("rect",
          xmin = domains$start[[k]], xmax = domains$end[[k]],
          ymin = tick_h, ymax = 1.0,
          fill  = dom_fill, color = NA, alpha = 0.55
        ) +
        ggplot2::annotate("text",
          x        = (domains$start[[k]] + domains$end[[k]]) / 2.0,
          y        = 0.995,
          label    = domains$name[[k]],
          size     = 2.8, hjust = 0.5, vjust = 1,
          color    = "#444444", fontface = "italic"
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
    y_lo   <- lanes$y_lo[[k]]
    y_hi   <- lanes$y_hi[[k]]
    y_mid  <- lanes$y_mid[[k]]
    lane_h <- y_hi - y_lo

    # Peptides for this lane
    lane_peps <- if (has_mc) peps[peps$missed_cleavages == mc_val, ,
                                   drop = FALSE] else peps

    lane_valid   <- lane_peps[lane_peps$length >= length_range[[1L]] &
                                lane_peps$length <= length_range[[2L]], ,
                               drop = FALSE]
    lane_invalid <- lane_peps[lane_peps$length  < length_range[[1L]] |
                                lane_peps$length  > length_range[[2L]], ,
                               drop = FALSE]

    # ── Greedy packing: stack overlapping peptides into sub-rows ─────────────
    # Each peptide gets a `track` integer (1 = bottom, 2 = above, ...).
    # For MC=0 tryptic digests there are no overlaps so n_tracks == 1.
    # For MC>=1 merged peptides share residues; packing yields 2-4 tracks.
    # y-bounds are pre-computed into data-frame columns to avoid the closure trap.
    if (nrow(lane_invalid) > 0L) {
      lane_invalid  <- .pack_peptides(lane_invalid)
      n_tracks_i    <- max(lane_invalid$track)
      inv_inner_lo  <- y_mid - 0.28 * lane_h
      inv_inner_hi  <- y_mid + 0.28 * lane_h
      track_h_i     <- (inv_inner_hi - inv_inner_lo) / n_tracks_i
      margin_i      <- 0.05 * track_h_i
      lane_invalid$.y_lo <- inv_inner_lo +
        (lane_invalid$track - 1L) * track_h_i + margin_i
      lane_invalid$.y_hi <- inv_inner_lo +
        lane_invalid$track       * track_h_i - margin_i
    }
    if (nrow(lane_valid) > 0L) {
      lane_valid   <- .pack_peptides(lane_valid)
      n_tracks_v   <- max(lane_valid$track)
      inner_lo     <- y_lo + 0.04 * lane_h
      inner_hi     <- y_hi - 0.04 * lane_h
      track_h_v    <- (inner_hi - inner_lo) / n_tracks_v
      margin_v     <- 0.05 * track_h_v
      lane_valid$.y_lo  <- inner_lo +
        (lane_valid$track - 1L) * track_h_v + margin_v
      lane_valid$.y_hi  <- inner_lo +
        lane_valid$track       * track_h_v - margin_v
      lane_valid$.y_mid <- (lane_valid$.y_lo + lane_valid$.y_hi) / 2.0
    }

    # Protein backbone (thin bar at lane centre) – annotate() is eager, safe
    p <- p + ggplot2::annotate("rect",
      xmin = 0.5, xmax = protein_length + 0.5,
      ymin = y_mid - 0.07 * lane_h, ymax = y_mid + 0.07 * lane_h,
      fill = "#D0D6E0", color = "#A8AEBA", linewidth = 0.35
    )

    # ── Gap overlays (annotate() is eager, safe with loop vars) ──────────────
    if (nrow(gap_df) > 0L) {
      p <- p + ggplot2::annotate("rect",
        xmin = gap_df$xmin, xmax = gap_df$xmax,
        ymin = y_lo + 0.03 * lane_h, ymax = y_hi - 0.03 * lane_h,
        fill  = .pepvet_pal$gap, alpha = 0.12
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
      label_v        <- lane_valid[
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
      x        = protein_length + 3L,
      y        = y_mid,
      label    = mc_label,
      hjust    = 0, vjust = 0.5,
      size     = 3.0, color = "#444444", fontface = "bold"
    )
  }

  # ── Lane separator lines (between lanes only) ─────────────────────────────
  if (nrow(lanes) > 1L) {
    for (k in seq_len(nrow(lanes) - 1L)) {
      sep_y <- (lanes$y_hi[[k]] + lanes$y_lo[[k + 1L]]) / 2.0
      p <- p + ggplot2::annotate("segment",
        x    = 0, xend = protein_length,
        y    = sep_y, yend = sep_y,
        color = .pepvet_pal$separator, linewidth = 0.4, linetype = "dotted"
      )
    }
  }

  # ── Cleavage-site ticks ────────────────────────────────────────────────────
  # efficiency is used as a proper color aesthetic (not I()) so ggplot2 adds it
  # to the legend panel alongside the fill legend.
  has_cs <- !is.null(cleavage_sites) && tick_h > 0
  if (has_cs) {
    cs_df              <- cleavage_sites
    cs_df$eff_level    <- factor(
      cs_df$efficiency,
      levels = c("high", "medium", "low"),
      labels = c("High", "Medium", "Low")
    )
    tick_top <- tick_h * 0.85

    p <- p +
      ggplot2::geom_segment(
        data = cs_df,
        ggplot2::aes(x    = .data$position, xend = .data$position,
                     y    = 0,              yend = tick_top,
                     color = .data$eff_level),
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
    # values rescaled over the display range [-2, 2]
    gravy_stops <- c(-2.0, -0.5, 0.6, 2.0)
    rescaled    <- (gravy_stops - (-2.0)) / 4.0
    p <- p + ggplot2::scale_fill_gradientn(
      colors = c(.pepvet_pal$brand, .pepvet_pal$good,
                  .pepvet_pal$moderate, .pepvet_pal$poor),
      values = rescaled,
      name   = "GRAVY",
      guide  = ggplot2::guide_colorbar(
        barwidth  = ggplot2::unit(8, "lines"),
        barheight = ggplot2::unit(0.55, "lines"),
        direction = "horizontal",
        title.position = "top"
      )
    )
  } else {
    color_vals <- c(
      "Valid"     = .pepvet_pal$valid,
      "Too short" = .pepvet_pal$too_short,
      "Too long"  = .pepvet_pal$too_long,
      "Invalid"   = "#B8C2CC"
    )
    if (color_by == "validity") {
      color_vals <- color_vals[c("Valid", "Invalid")]
    }
    p <- p + ggplot2::scale_fill_manual(
      values = color_vals,
      name   = NULL,
      na.value = "#B8C2CC",
      guide  = ggplot2::guide_legend(
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
      name  = "Cleavage efficiency",
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
    parts   <- paste0(as.integer(eff_tbl), " ", tolower(names(eff_tbl)))
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
      title    = auto_title,
      subtitle = sprintf(
        "%.0f%% sequence coverage (MC=0 valid peptides)  \u00b7  %d gap region(s)  \u00b7  protein length %d aa%s",
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
