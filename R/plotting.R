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
        color = .pepvet_pal$brand_dark,
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
    .pepvet_theme() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(
        angle  = 90,
        hjust  = 0.5,
        margin = ggplot2::margin(r = 4)
      )
    )
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
          color  = .pepvet_pal$brand_dark,
          margin = ggplot2::margin(b = 10)
        )
      )
    ) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face   = "bold",
        size   = 14,
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

