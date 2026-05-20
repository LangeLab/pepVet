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

#' @importFrom rlang .data check_installed %||%
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom tools file_ext
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
  brand       = "#2C5F8A",
  brand_dark  = "#1A3D5C",
  brand_light = "#7BAED4",

  # Length-class categories
  valid      = "#2C5F8A",
  too_short  = "#E8A838",
  too_long   = "#C94040",

  # Scoring tiers
  good       = "#27AE60",
  moderate   = "#E8A838",
  poor       = "#C94040",

  # Verdict colors (named, matches the tier values above)
  verdict = c(
    Good     = "#27AE60",
    Moderate = "#E8A838",
    Poor     = "#C94040"
  ),

  # Component score colors (5 + S_unique)
  component = c(
    S_coverage = "#2C5F8A",
    S_length   = "#27AE60",
    S_count    = "#E8A838",
    S_hydro    = "#8B5E99",
    S_charge   = "#4AAFB0",
    S_unique   = "#B8C2CC"
  ),

  # Coverage / cleavage map structural elements
  protein_bg  = "#D8DDE6",
  protein_brd = "#AAAAAA",
  backbone_fill = "#D0D6E0",
  backbone_brd  = "#A8AEBA",
  invalid_pep   = "#C5CDD8",
  covered       = "#2C5F8A",
  gap           = "#C94040",
  cleavage_bg   = "#E8EDF3",
  na_gray       = "#B8C2CC",
  overlap = c(
    "Not detected"      = "#FFFFFF",
    "Detected once"     = "#DCEAF5",
    "Detected twice"    = "#7BAED4",
    "Detected 3+ times" = "#2C5F8A"
  ),

  # Background shading for valid ranges
  shade      = "#EDF6F0",
  neutral    = "#F4F6F9",
  separator  = "#DDDDDD",

  # Theme element colors
  text_subtitle   = "#666666",
  text_axis_title = "#444444",
  text_axis_tick  = "#555555",
  axis_tick       = "#CCCCCC",
  grid_major      = "#EBEBEB",
  strip_bg        = "#F0F4F8",
  text_secondary  = "#999999",
  text_dark       = "#333333",

  # Zone shading (verdict background bands)
  zone_moderate = "#FFF3E0",
  zone_poor     = "#FFEBEE",

  # Badge colors
  badge_gold_text = "#7A5A00",
  badge_gold_fill = "#FFF5CC",

  # Heatmap gradient midpoint
  heatmap_mid = "#FFFAEC"
)

# Snapshot defaults at build time for pepvet_plot_config_reset()
# Store in an environment (not locked like namespace bindings)
.pepvet_config_env <- new.env(parent = emptyenv())
.pepvet_config_env$pal_default <- .pepvet_pal
.pepvet_config_env$params_default <- NULL  # set after .pepvet_params loads
.pepvet_config_env$theme_overrides <- list()

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
  out <- ggplot2::theme_minimal(base_size = base_size) +
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
        color  = .pepvet_pal$text_subtitle,
        margin = ggplot2::margin(b = 6)
      ),

      # Axes
      axis.title = ggplot2::element_text(size = base_size - 1, color = .pepvet_pal$text_axis_title),
      axis.text = ggplot2::element_text(size = base_size - 2, color = .pepvet_pal$text_axis_tick),
      axis.ticks = ggplot2::element_line(color = .pepvet_pal$axis_tick, linewidth = 0.3),

      # Grid
      panel.grid.major = ggplot2::element_line(color = .pepvet_pal$grid_major, linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill      = NA,
        color     = .pepvet_pal$separator,
        linewidth = 0.5
      ),

      # Strip (for faceted plots)
      strip.background = ggplot2::element_rect(fill = .pepvet_pal$strip_bg, color = .pepvet_pal$separator),
      strip.text = ggplot2::element_text(
        face  = "bold",
        size  = base_size - 0.5,
        color = .pepvet_pal$brand
      ),

      # Legend
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.55, "lines"),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(
        size = base_size - 1.5,
        face = "bold",
        color = .pepvet_pal$text_axis_title
      ),
      legend.background = ggplot2::element_rect(fill = NA, color = NA),

      # Plot margins
      plot.margin = ggplot2::margin(8, 10, 6, 10)
    )

  # Apply user theme overrides if set
  overrides <- .pepvet_config_env$theme_overrides
  if (length(overrides) > 0L) {
    out <- out + do.call(ggplot2::theme, overrides)
  }

  out
}


# -- Shared classification helpers -------------------------------------------

#' Classify peptide lengths into validity categories
#' @noRd
.classify_length <- function(lengths, length_range) {
  lo <- length_range[[1L]]
  hi <- length_range[[2L]]
  factor(
    ifelse(lengths < lo, "Too short",
      ifelse(lengths > hi, "Too long", "Valid")),
    levels = c("Valid", "Too short", "Too long")
  )
}

#' Length-class color map
#' @noRd
.length_class_colors <- function() {
  c(
    "Valid"     = .pepvet_pal$valid,
    "Too short" = .pepvet_pal$too_short,
    "Too long"  = .pepvet_pal$too_long
  )
}

#' Classify numeric scores into verdict tiers
#' @noRd
.classify_verdict <- function(x) {
  ifelse(x >= .get_param("verdict_good"), "Good",
    ifelse(x >= .get_param("verdict_moderate"), "Moderate", "Poor"))
}

#' Nice x-axis step for protein-length scales
#' @noRd
.nice_x_step <- function(protein_length) {
  max(50L, as.integer(round(protein_length / 10.0 / 50.0) * 50L))
}

#' Guard against empty peptide tables
#' @noRd
.validate_nonempty <- function(data, name = "data", class = "pepvet_error_invalid_input") {
  if (is.data.frame(data) && nrow(data) == 0L) {
    .abort(
      "{.arg {name}} must have at least one row.",
      class = class
    )
  }
  invisible(data)
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
  # UniProt format: sp|ACC|GENE or tr|ACC|GENE
  m <- regmatches(
    protein_id,
    regexpr("^[a-z]+\\|([A-Z0-9]+)\\|([A-Z0-9_./-]+)", protein_id)
  )
  if (length(m) == 1L && nchar(m) > 0L) {
    parts <- strsplit(m, "\\|")[[1L]]
    return(paste0(parts[[2L]], "  (", parts[[3L]], ")"))
  }
  # NCBI RefSeq: NP_001234.1 or XP_...
  m2 <- regmatches(protein_id, regexpr("^([A-Z]{2}_[0-9]+(\\.[0-9]+)?)", protein_id))
  if (length(m2) == 1L && nchar(m2) > 0L) {
    return(m2)
  }
  # Generic: take first space-delimited token (FASTA >header convention)
  m3 <- regmatches(protein_id, regexpr("^[^ ]+", protein_id))
  if (length(m3) == 1L && nchar(m3) > 0L && nchar(m3) <= 42L) {
    return(m3)
  }
  # Truncate as last resort
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
  if (nrow(peps) == 0L) {
    .abort(
      "No peptides found in {.arg result}. The digest produced zero peptides.",
      class = "pepvet_error_invalid_digest_result"
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

  peps$length_class <- .classify_length(peps$length, length_range)

  n_valid <- sum(peps$length_class == "Valid")
  n_total <- nrow(peps)
  pct <- round(100 * n_valid / n_total, 1)

  class_colors <- .length_class_colors()

  # Sensible x-axis upper limit
  x_max <- max(peps$length, na.rm = TRUE) + 1L

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

  # Suitable bin count for the data range (Freedman-Diaconis rule, clamped)
  gravy_iqr <- diff(stats::quantile(peps$gravy, c(0.25, 0.75), na.rm = TRUE, names = FALSE))
  fd_bw <- 2 * gravy_iqr / length(peps$gravy)^(1/3)
  n_bins <- if (is.finite(fd_bw) && fd_bw > 0) {
    max(15L, min(40L, as.integer(diff(range(peps$gravy, na.rm = TRUE)) / fd_bw)))
  } else {
    30L
  }

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


#' Reconstruct a full single-protein sequence from MC=0 peptides
#'
#' Uses the non-overlapping MC=0 digest fragments to rebuild the original
#' sequence so residue-level plots do not need the raw FASTA input again.
#'
#' @param peps Peptide tibble from [evaluate_digest()].
#' @return Length-1 character string containing the full protein sequence.
#' @noRd
.reconstruct_sequence_from_peptides <- function(peps) {
  protein_length <- max(peps$end, na.rm = TRUE)
  mc0_peps <- if ("missed_cleavages" %in% names(peps)) {
    peps[peps$missed_cleavages == 0L, , drop = FALSE]
  } else {
    peps
  }

  if (nrow(mc0_peps) == 0L) {
    .abort(
      "Could not reconstruct the full protein sequence because no MC=0 peptides were available.",
      class = "pepvet_error_invalid_digest_result"
    )
  }

  mc0_peps <- mc0_peps[order(mc0_peps$start, mc0_peps$end), , drop = FALSE]
  sequence_chars <- rep(NA_character_, protein_length)

  for (index in seq_len(nrow(mc0_peps))) {
    peptide_chars <- strsplit(mc0_peps$peptide[[index]], "", fixed = TRUE)[[1L]]
    residue_positions <- seq.int(mc0_peps$start[[index]], mc0_peps$end[[index]])

    if (length(peptide_chars) != length(residue_positions)) {
      .abort(
        "Peptide coordinates do not match peptide sequence length during sequence reconstruction.",
        class = "pepvet_error_invalid_digest_result"
      )
    }

    existing_chars <- sequence_chars[residue_positions]
    has_mismatch <- !is.na(existing_chars) & existing_chars != peptide_chars
    if (any(has_mismatch)) {
      .abort(
        "Peptide table contains inconsistent residue assignments and cannot be rendered as a sequence map.",
        class = "pepvet_error_invalid_digest_result"
      )
    }

    sequence_chars[residue_positions] <- peptide_chars
  }

  if (anyNA(sequence_chars)) {
    .abort(
      "Could not reconstruct a complete protein sequence from the peptide table.",
      class = "pepvet_error_invalid_digest_result"
    )
  }

  paste(sequence_chars, collapse = "")
}


#' Build wrapped residue tiles with peptide overlap counts
#'
#' Counts how many peptides cover each residue, then formats the sequence into
#' wrapped display rows for residue-letter heatmap style plots.
#'
#' @param peps Peptide tibble from [evaluate_digest()].
#' @param protein_length Integer protein length.
#' @param length_range Optional integer vector of length 2. When `NULL`, all
#'   peptides at all MC levels are counted.
#' @param missed_cleavages Integer MC level to filter to when `length_range`
#'   is not `NULL`. Default `1L`.
#' @param residues_per_line Positive integer wrap width.
#' @return Data frame with one row per residue and columns for plotting.
#' @noRd
.build_peptide_overlap_df <- function(peps, protein_length,
                                      length_range = NULL,
                                      missed_cleavages = 1L,
                                      residues_per_line = 50L) {
  protein_sequence <- .reconstruct_sequence_from_peptides(peps)
  overlap_counts <- integer(protein_length)

  # When length_range is NULL, count all MC levels (user opted in).
  # Otherwise filter to the specified MC level (default 1L).
  has_mc <- "missed_cleavages" %in% names(peps)
  count_all_mc <- is.null(length_range) && has_mc

  overlap_peps <- if (count_all_mc) {
    peps
  } else if (has_mc) {
    peps[peps$missed_cleavages == missed_cleavages, , drop = FALSE]
  } else {
    peps
  }

  if (!is.null(length_range)) {
    normalized_length_range <- .validate_length_range(length_range)
    overlap_peps <- overlap_peps[
      overlap_peps$length >= normalized_length_range[[1L]] &
        overlap_peps$length <= normalized_length_range[[2L]],
      ,
      drop = FALSE
    ]
  }

  if (nrow(overlap_peps) > 0L) {
    for (index in seq_len(nrow(overlap_peps))) {
      residue_positions <- seq.int(overlap_peps$start[[index]], overlap_peps$end[[index]])
      overlap_counts[residue_positions] <- overlap_counts[residue_positions] + 1L
    }
  }

  residue_positions <- seq_len(protein_length)
  residues <- strsplit(protein_sequence, "", fixed = TRUE)[[1L]]
  line_index <- ((residue_positions - 1L) %/% residues_per_line) + 1L
  column_index <- ((residue_positions - 1L) %% residues_per_line) + 1L

  line_start <- ((line_index - 1L) * residues_per_line) + 1L
  line_end <- pmin(line_index * residues_per_line, protein_length)
  line_label <- paste0("Residues ", line_start, "-", line_end)

  overlap_class <- ifelse(
    overlap_counts >= 3L, "Detected 3+ times",
    ifelse(
      overlap_counts == 2L, "Detected twice",
      ifelse(overlap_counts == 1L, "Detected once", "Not detected")
    )
  )

  data.frame(
    position = residue_positions,
    residue = residues,
    overlap_count = overlap_counts,
    overlap_class = factor(overlap_class, levels = names(.pepvet_pal$overlap)),
    letter_color = ifelse(overlap_counts >= 2L, "light", "dark"),
    line_label = factor(line_label, levels = unique(line_label)),
    column_index = column_index,
    stringsAsFactors = FALSE
  )
}

# -- Greedy interval packing for non-overlapping peptide display (shared helper)

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
  x_step <- .nice_x_step(protein_length)

  p <- ggplot2::ggplot() +
    # Full protein background bar
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = protein_length + 0.5,
      ymin = 0.30, ymax = 0.70,
      fill = .pepvet_pal$protein_bg,
      color = .pepvet_pal$protein_brd,
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
      fill = .pepvet_pal$invalid_pep,
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
      breaks = seq(0L, protein_length, by = x_step),
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
      color = .pepvet_pal$text_dark,
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
      limits = c(0, max(max(vals) + 0.08, 1.05)),
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

#' Configure pepVet plot appearance
#'
#' Override default colors, numeric parameters, and theme elements used by
#' all pepVet plotting functions. Changes persist for the session.
#' Call with no arguments to view the current configuration.
#'
#' @param palette Named list of color overrides. Names must match existing
#'   `.pepvet_pal` entries (e.g. `list(brand = "#004488", good = "#2ECC71")`).
#'   Sub-lists like `verdict`, `component`, and `overlap` are replaced entirely.
#' @param params Named list of parameter overrides. Names must match existing
#'   `.pepvet_params` entries (e.g. `list(verdict_good = 0.70, scatter_alpha = 0.90)`).
#' @param theme Named list of ggplot2 theme element overrides. Passed directly
#'   to [ggplot2::theme()] and applied on top of the base `.pepvet_theme()`.
#'   (e.g. `list(legend.position = "right", plot.title = element_text(size = 14))`).
#'
#' @return Invisibly returns a list with current `palette`, `params`, and `theme`.
#' @family plot-utils
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # View current config
#'   pepvet_plot_config()
#'
#'   # Customize brand color and Good threshold
#'   pepvet_plot_config(
#'     palette = list(brand = "#004488", good = "#2ECC71"),
#'     params  = list(verdict_good = 0.70)
#'   )
#'
#'   # Reset to defaults
#'   pepvet_plot_config_reset()
#' }
pepvet_plot_config <- function(palette = NULL, params = NULL, theme = NULL) {
  env <- .pepvet_config_env

  if (is.null(palette) && is.null(params) && is.null(theme)) {
    return(invisible(list(
      palette = env$pal,
      params  = env$params,
      theme   = env$theme_overrides
    )))
  }

  if (!is.null(palette)) {
    if (!is.list(palette) || is.null(names(palette))) {
      .abort(
        "{.arg palette} must be a named list.",
        class = "pepvet_error_invalid_config"
      )
    }
    unknown <- setdiff(names(palette), names(env$pal))
    if (length(unknown) > 0L) {
      .abort(
        "Unknown palette {.field {unknown}}.",
        "i" = "Valid names: {.val {names(env$pal)}}.",
        class = "pepvet_error_invalid_config"
      )
    }
    env$pal[names(palette)] <- palette
  }

  if (!is.null(params)) {
    if (!is.list(params) || is.null(names(params))) {
      .abort(
        "{.arg params} must be a named list.",
        class = "pepvet_error_invalid_config"
      )
    }
    unknown <- setdiff(names(params), names(env$params))
    if (length(unknown) > 0L) {
      .abort(
        "Unknown param {.field {unknown}}.",
        "i" = "Valid names: {.val {names(env$params)}}.",
        class = "pepvet_error_invalid_config"
      )
    }
    env$params[names(params)] <- params
  }

  if (!is.null(theme)) {
    if (!is.list(theme) || is.null(names(theme))) {
      .abort(
        "{.arg theme} must be a named list.",
        class = "pepvet_error_invalid_config"
      )
    }
    env$theme_overrides <- theme
  }

  invisible(list(
    palette = env$pal,
    params  = env$params,
    theme   = env$theme_overrides
  ))
}


#' Reset pepVet plot configuration to defaults
#'
#' Restores all colors, parameters, and theme overrides to the package defaults.
#'
#' @return Invisibly returns a list with current `palette`, `params`, and `theme`.
#' @family plot-utils
#' @export
pepvet_plot_config_reset <- function() {
  env <- .pepvet_config_env
  env$pal[] <- env$pal_default
  if (!is.null(env$params_default)) env$params[] <- env$params_default
  env$theme_overrides <- list()
  invisible(list(
    palette = env$pal,
    params  = env$params,
    theme   = env$theme_overrides
  ))
}


#' Save a pepVet figure with publication-ready defaults
#'
#' Wraps [ggplot2::ggsave()] with pepVet's recommended defaults: auto-sizing
#' based on whether the plot is a multi-panel patchwork or a single panel,
#' 300 DPI, anti-aliased PNG via ragg when available, and white background.
#' All arguments in `...` are passed to [ggplot2::ggsave()] and can override
#' the defaults.
#'
#' @param plot A ggplot or patchwork object produced by any pepVet plot function.
#' @param filename Character path for the output file. Extensions `.png`, `.pdf`,
#'   `.svg`, etc. are handled by [ggplot2::ggsave()]. Defaults to `"pepvet_plot.png"`
#'   in the working directory.
#' @param width,height Numeric. Plot dimensions in inches. When `NULL` (default),
#'   auto-sized: single-panel = 10x7, multi-panel patchwork = 14x10.
#' @param dpi Numeric. Resolution in dots per inch. Defaults to `300` (publication).
#' @param bg Character. Background color. Defaults to `"white"`.
#' @param device Device to use. When `NULL` (default) and the filename extension
#'   is `.png`, tries [ragg::agg_png()] for anti-aliased output, falling back
#'   to `"png"`. For other extensions, the device is inferred from the filename.
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return The saved file path invisibly.
#' @family plot-utils
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   p <- plot_digest_profile(res)
#'   tmp <- tempfile(fileext = ".png")
#'   pepvet_save_figure(p, tmp)
#' }
pepvet_save_figure <- function(plot,
                               filename = "pepvet_plot.png",
                               width    = NULL,
                               height   = NULL,
                               dpi      = 300,
                               bg       = "white",
                               device   = NULL,
                               ...) {
  # Auto-size: patchwork gets larger default canvas
  if (is.null(width) || is.null(height)) {
    is_patchwork <- inherits(plot, "patchwork")
    if (is_patchwork) {
      if (is.null(width))  width  <- 14
      if (is.null(height)) height <- 10
    } else {
      if (is.null(width))  width  <- 10
      if (is.null(height)) height <- 7
    }
  }

  # Auto-device: prefer ragg for anti-aliased PNG
  if (is.null(device)) {
    ext <- tolower(tools::file_ext(filename))
    if (identical(ext, "png") &&
        requireNamespace("ragg", quietly = TRUE)) {
      device <- ragg::agg_png
    }
  }

  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = bg,
    device   = device,
    ...
  )

  invisible(normalizePath(filename, mustWork = FALSE))
}


#' pepVet manuscript theme
#'
#' A compact theme for journal figures. Smaller base text (9pt), thinner grid
#' lines, tighter margins. Use by adding to any pepVet plot:
#' `plot + pepvet_theme_manuscript()`.
#'
#' @return A ggplot2 theme object.
#' @family plot-utils
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   plot_length_distribution(res) + pepvet_theme_manuscript()
#' }
pepvet_theme_manuscript <- function() {
  .pepvet_theme(base_size = 9) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.25),
      plot.margin = ggplot2::margin(4, 6, 4, 6),
      legend.key.size = ggplot2::unit(0.4, "lines")
    )
}


#' pepVet presentation theme
#'
#' A bold theme for talks and posters. Larger base text (14pt), thicker grid
#' lines, wider margins. Use by adding to any pepVet plot:
#' `plot + pepvet_theme_presentation()`.
#'
#' @return A ggplot2 theme object.
#' @family plot-utils
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#'   res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#'   plot_length_distribution(res) + pepvet_theme_presentation()
#' }
pepvet_theme_presentation <- function() {
  .pepvet_theme(base_size = 14) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.5),
      plot.margin = ggplot2::margin(12, 14, 10, 14),
      legend.key.size = ggplot2::unit(0.7, "lines"),
      axis.ticks = ggplot2::element_line(linewidth = 0.5),
      panel.border = ggplot2::element_rect(
        fill = NA, color = "#DDDDDD", linewidth = 0.8
      )
    )
}
