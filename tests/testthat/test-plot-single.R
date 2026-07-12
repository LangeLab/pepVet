# Tests aligned with the corresponding plotting source domain.
test_that("plot_digest_profile returns a patchwork object for BSA / trypsin", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  p <- plot_digest_profile(res)

  expect_s3_class(p, "patchwork")
  expect_gte(length(p$patches$plots), 2L)
  expect_match(p$patches$annotation$title, "trypsin")
  expect_equal(as.numeric(p$patches$layout$heights), c(3, 1.8, 2.2))
})


test_that("plot_digest_profile: custom title and length/gravy ranges", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  p <- plot_digest_profile(
    res,
    length_range = c(6L, 30L),
    gravy_range  = c(-1.5, 1.0),
    title        = "Custom ranges"
  )
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "Custom ranges")
  panel_a <- p$patches$plots[[1L]]$patches$plots[[1L]]
  expect_true(any(vapply(
    panel_a$layers,
    function(layer) {
      is.data.frame(layer$data) &&
        any(layer$data$xmin == 5.5 & layer$data$xmax == 30.5)
    },
    logical(1L)
  )))
})

test_that("plot_digest_profile rejects malformed ranges without warnings", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  expect_no_warning(
    expect_error(
      plot_digest_profile(res, length_range = "bad"),
      class = "pepvet_error_invalid_length_range"
    )
  )
  expect_no_warning(
    expect_error(
      plot_digest_profile(res, gravy_range = "bad"),
      class = "pepvet_error_invalid_gravy_range"
    )
  )
})

test_that("single-result plots inherit resolved range metadata", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  metadata_result <- .bsa_result(mc = 1L)
  metadata_result$params$length_range <- c(1L, 100L)
  metadata_result$params$gravy_range <- c(-1.5, 1.0)

  profile <- plot_digest_profile(metadata_result)
  profile_length <- profile$patches$plots[[1L]]$patches$plots[[1L]]
  expect_true(any(vapply(profile_length$layers, function(layer) {
    is.data.frame(layer$data) && "xmin" %in% names(layer$data) &&
      any(layer$data$xmin == 0.5 & layer$data$xmax == 100.5)
  }, logical(1L))))

  coverage <- plot_coverage_map(metadata_result)
  fill_values <- unlist(lapply(coverage$layers, function(layer) {
    if (is.data.frame(layer$data) && "fill_cat" %in% names(layer$data)) {
      return(as.character(layer$data$fill_cat))
    }
    character(0L)
  }))
  expect_true(length(fill_values) > 0L && all(fill_values == "Valid"))

  overlap_default <- plot_peptide_overlap_map(.bsa_result(mc = 1L))
  overlap_metadata <- plot_peptide_overlap_map(metadata_result)
  expect_gt(
    sum(overlap_metadata$data$overlap_count > 0L),
    sum(overlap_default$data$overlap_count > 0L)
  )

  cleavage <- plot_cleavage_map(metadata_result)
  expect_match(cleavage$labels$subtitle, "79 / 79 valid fragments")
})

test_that("plot_digest_profile: proteome-aware result includes S_unique panel", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  proteome_path <- system.file(
    "extdata", "small_proteome_50_proteins.fasta",
    package = "pepVet"
  )
  proteome_digest <- digest_protein(proteome_path, enzyme = "trypsin")
  res <- evaluate_digest(.bsa_path, enzyme = "trypsin", proteome = proteome_digest)

  p <- plot_digest_profile(res)
  expect_s3_class(p, "patchwork")
})

test_that("plot_digest_profile: chymotrypsin-high renders without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  # chymotrypsin-high produces many short aromatic peptides: exercises
  # the length/GRAVY panels with a very different peptide population
  res <- .fix_bsa_chymotryp
  p <- plot_digest_profile(res)
  expect_s3_class(p, "patchwork")
})

# Input validation.

test_that("plot_digest_profile errors on non-list input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  expect_error(
    plot_digest_profile("not a list"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_digest_profile errors on multi-protein result", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  small_path <- system.file(
    "extdata", "small_proteome_50_proteins.fasta",
    package = "pepVet"
  )
  # compare_digests returns a flat tibble, not an evaluate_digest list.
  # we need a proper evaluate_digest result with multiple protein IDs.
  # Fake one by duplicating BSA rows with a different protein_id.
  res <- .bsa_result()
  extra <- res$peptides
  extra$protein_id <- "fake|P99999|FAKE_PROTEIN"
  res$peptides <- rbind(res$peptides, extra)

  expect_error(
    plot_digest_profile(res),
    class = "pepvet_error_multi_protein"
  )
})

test_that("plot_digest_profile errors on missing required columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  bad_peps <- res$peptides[, c("protein_id", "peptide"), drop = FALSE]
  bad_result <- list(
    scores   = res$scores,
    peptides = bad_peps,
    params   = res$params
  )

  expect_error(
    plot_digest_profile(bad_result),
    class = "pepvet_error_invalid_digest_result"
  )
})

# plot_coverage_map.

test_that("plot_coverage_map returns a ggplot for BSA / trypsin MC=0", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_coverage_map(res)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "Residue position")
  expect_true(any(vapply(
    p$layers,
    function(layer) {
      inherits(layer$geom, "GeomRect") &&
        is.data.frame(layer$data) && "start" %in% names(layer$data)
    },
    logical(1L)
  )))
})

test_that("plot_coverage_map supports peptide tables without a missed-cleavage column", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  res$peptides$missed_cleavages <- NULL
  p <- plot_coverage_map(res)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "MC=0 valid peptides")
})

test_that("plot_coverage_map returns a ggplot for MC=1 (multi-lane + packing)", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  p <- plot_coverage_map(res)
  expect_s3_class(p, "ggplot")
  lane_labels <- unlist(lapply(p$layers, function(layer) {
    if (inherits(layer$geom, "GeomText") &&
        "label" %in% names(layer$aes_params)) {
      return(as.character(layer$aes_params$label))
    }
    character(0L)
  }))
  expect_true(all(c("MC = 0", "MC = 1") %in% lane_labels))
})

test_that("plot_coverage_map: color_by = 'length_class' renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 2L)
  p <- plot_coverage_map(res, color_by = "length_class")
  expect_s3_class(p, "ggplot")
  fill_scale <- p$scales$get_scales("fill")
  expect_true(all(c("Valid", "Too short", "Too long") %in%
    names(fill_scale$palette(3L))))
})

test_that("plot_coverage_map: color_by = 'hydrophobicity' renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  p <- plot_coverage_map(res, color_by = "hydrophobicity")
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(
    p$layers,
    function(layer) {
      is.data.frame(layer$data) && "fill_val" %in% names(layer$data)
    },
    logical(1L)
  )))
})

test_that("plot_coverage_map: cleavage_sites overlay renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  cs <- .bsa_cs()
  p <- plot_coverage_map(res, cleavage_sites = cs)
  expect_s3_class(p, "ggplot")
  tick_layer <- p$layers[[which(vapply(
    p$layers,
    function(layer) {
      inherits(layer$geom, "GeomSegment") &&
        is.data.frame(layer$data) && "eff_level" %in% names(layer$data)
    },
    logical(1L)
  ))[[1L]]]]
  expect_setequal(as.character(tick_layer$data$eff_level),
    c("High", "Medium", "Low"))
})

test_that("plot_coverage_map: domains overlay renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  domains <- data.frame(
    name = c("Domain A", "Domain B"),
    start = c(1L, 200L),
    end = c(150L, 400L),
    stringsAsFactors = FALSE
  )
  p <- plot_coverage_map(res, domains = domains)
  expect_s3_class(p, "ggplot")
  domain_labels <- unlist(lapply(p$layers, function(layer) {
    if (inherits(layer$geom, "GeomText") &&
        "label" %in% names(layer$aes_params)) {
      return(as.character(layer$aes_params$label))
    }
    character(0L)
  }))
  expect_true(all(c("Domain A", "Domain B") %in% domain_labels))
})

test_that("plot_coverage_map: custom title is accepted", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_coverage_map(res, title = "My custom title")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "My custom title")
})

test_that("plot_coverage_map: H3.1 short protein renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result()
  cs <- .h3_cs()
  p <- expect_no_warning(plot_coverage_map(res, cleavage_sites = cs))
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "sequence coverage")
})

test_that("plot_coverage_map: chymotrypsin-high produces different lane layout", {
  skip_if_not_installed("ggplot2")

  # Chymotrypsin-high yields shorter, more numerous peptides than trypsin.
  # exercises the lane-packing algorithm with a denser peptide map
  res <- .fix_bsa_chymotryp
  p <- plot_coverage_map(res)
  expect_s3_class(p, "gg")
  expect_true(sum(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomRect") && is.data.frame(layer$data) &&
      "track" %in% names(layer$data)
  }, logical(1L))) >= 2L)
})

test_that("plot_coverage_map errors on non-list input", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_coverage_map("not a list"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_coverage_map errors on multi-protein result", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  extra <- res$peptides
  extra$protein_id <- "fake|P99999|FAKE_PROTEIN"
  res$peptides <- rbind(res$peptides, extra)

  expect_error(
    plot_coverage_map(res),
    class = "pepvet_error_multi_protein"
  )
})

test_that("plot_coverage_map validates ranges and overlay schemas", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  expect_error(
    plot_coverage_map(res, length_range = c(25L, 7L)),
    class = "pepvet_error_invalid_length_range"
  )
  expect_error(
    plot_coverage_map(
      res,
      cleavage_sites = data.frame(position = 2L)
    ),
    class = "pepvet_error_invalid_cleavage_sites"
  )
  expect_error(
    plot_coverage_map(
      res,
      cleavage_sites = data.frame(position = 2L, efficiency = "unknown")
    ),
    class = "pepvet_error_invalid_cleavage_sites"
  )
  expect_no_warning(
    expect_error(
      plot_coverage_map(
        res,
        cleavage_sites = data.frame(
          position = 2L,
          efficiency = I(list("high"))
        )
      ),
      class = "pepvet_error_invalid_cleavage_sites"
    )
  )
  expect_error(
    plot_coverage_map(
      res,
      domains = data.frame(name = "bad", start = 0L, end = 10L)
    ),
    class = "pepvet_error_invalid_domains"
  )
  expect_error(
    plot_coverage_map(
      res,
      domains = data.frame(name = "bad", start = 1L, end = 9999L)
    ),
    class = "pepvet_error_invalid_domains"
  )
  expect_error(
    plot_coverage_map(
      res,
      domains = data.frame(name = "bad", start = 1L)
    ),
    class = "pepvet_error_invalid_domains"
  )
  expect_error(
    plot_coverage_map(res, color_by = "not-a-coloring"),
    class = "pepvet_error_invalid_input"
  )
})

# Restored section from the original plotting test surface.

# Restored section from the original plotting test surface.
test_that("plot_peptide_overlap_map returns a ggplot for BSA / trypsin MC=1", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  p <- plot_peptide_overlap_map(res)

  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), max(res$peptides$end))
  expect_true(all(c("overlap_count", "overlap_class", "line_label") %in%
    names(p$data)))
})

test_that("plot_peptide_overlap_map supports wrapped rows and all-peptide mode", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result(mc = 1L)

  p <- plot_peptide_overlap_map(
    res,
    length_range = NULL,
    residues_per_line = 20L
  )
  expect_s3_class(p, "ggplot")
  expect_equal(max(p$data$column_index), 20L)
  expect_gt(nlevels(p$data$line_label), 1L)
})

test_that("plot_peptide_overlap_map rejects invalid filtering inputs", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  expect_error(
    plot_peptide_overlap_map(res, length_range = c(25L, 7L)),
    class = "pepvet_error_invalid_length_range"
  )
  expect_error(
    plot_peptide_overlap_map(res, residues_per_line = 0L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    plot_peptide_overlap_map(res, residues_per_line = 2.5),
    class = "pepvet_error_invalid_input"
  )
  expect_no_warning(
    expect_error(
      plot_peptide_overlap_map(res, residues_per_line = 1e20),
      class = "pepvet_error_invalid_input"
    )
  )
  expect_error(
    plot_peptide_overlap_map(res, missed_cleavages = -1L),
    class = "pepvet_error_invalid_input"
  )
  expect_no_warning(
    expect_error(
      plot_peptide_overlap_map(res, missed_cleavages = 1e20),
      class = "pepvet_error_invalid_input"
    )
  )
  expect_error(
    plot_peptide_overlap_map(res, residues_per_line = "20"),
    class = "pepvet_error_invalid_input"
  )
  p <- plot_peptide_overlap_map(
    res, residues_per_line = 100L, title = "Wrapped overlap"
  )
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Wrapped overlap")
  expect_s3_class(
    plot_peptide_overlap_map(res, residues_per_line = 70L),
    "ggplot"
  )
})

test_that("plot_cleavage_map returns ggplot from evaluate_digest result", {
  res <- .fix_bsa_trypsin
  p <- plot_cleavage_map(res)
  expect_s3_class(p, "gg")
  expect_match(p$labels$subtitle, "cleavage sites")
  expect_true(any(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomSegment") && is.data.frame(layer$data)
  }, logical(1L))))
})

test_that("plot_cleavage_map renders without warnings (no cleavage_sites)", {
  res <- .fix_bsa_trypsin
  expect_no_warning(plot_cleavage_map(res))
})

test_that("plot_cleavage_map: custom title accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_cleavage_map(res, title = "BSA cleavage")
  expect_s3_class(p, "gg")
  expect_identical(p$labels$title, "BSA cleavage")
})

test_that("plot_cleavage_map uses supplied efficiency annotations", {
  cs <- .bsa_cs()
  p <- plot_cleavage_map(.fix_bsa_trypsin, cleavage_sites = cs)
  expect_s3_class(p, "ggplot")
  efficiency_layer <- p$layers[[which(vapply(
    p$layers,
    function(layer) {
      inherits(layer$geom, "GeomSegment") && is.data.frame(layer$data) &&
        "efficiency" %in% names(layer$data)
    },
    logical(1L)
  ))[[1L]]]]
  expect_setequal(as.character(efficiency_layer$data$efficiency),
    c("high", "medium", "low"))
  expect_identical(p$labels$caption, NULL)
})

test_that("plot_cleavage_map: asp-n endopeptidase (N-terminal cutter) renders", {
  # Asp-N cuts N-terminal to D. This exercises site inference from the opposite
  # end of peptides compared to trypsin/lysc
  res <- .fix_bsa_aspn
  p <- plot_cleavage_map(res)
  expect_s3_class(p, "gg")
  expect_gt(sum(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomSegment") && is.data.frame(layer$data) &&
      nrow(layer$data) > 0L
  }, logical(1L))), 0L)
})

test_that("plot_cleavage_map handles a digest with no internal sites", {
  res <- NULL
  expect_warning(
    res <- evaluate_digest("AAAA", enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  p <- plot_cleavage_map(res)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "0 cleavage sites")
  expect_false(any(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomSegment") && is.data.frame(layer$data) &&
      nrow(layer$data) > 0L
  }, logical(1L))))
})

test_that("plot_cleavage_map supports peptide tables without missed-cleavage values", {
  res <- .fix_bsa_trypsin
  res$peptides$missed_cleavages <- NULL
  p <- plot_cleavage_map(res)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "cleavage sites")
})

test_that("plot_cleavage_map validates annotation data", {
  expect_error(
    plot_cleavage_map(
      .fix_bsa_trypsin,
      cleavage_sites = data.frame(position = 2L)
    ),
    class = "pepvet_error_invalid_cleavage_sites"
  )
  expect_error(
    plot_cleavage_map(
      .fix_bsa_trypsin,
      cleavage_sites = data.frame(position = 2L, efficiency = "unknown")
    ),
    class = "pepvet_error_invalid_cleavage_sites"
  )
})

test_that("plot_cleavage_map errors on invalid result", {
  expect_error(
    plot_cleavage_map("not valid"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_weight_sensitivity renders batch instability semantics", {
  skip_if_not_installed("ggplot2")
  sensitivity <- sensitivity_analysis(
    .fix_batch_small,
    n_iter = 20L,
    chunk_size = 17L
  )
  plot <- plot_weight_sensitivity(sensitivity)

  expect_s3_class(plot, "ggplot")
  expect_true("verdict_instability" %in% names(plot$data))
  expect_identical(plot$labels$x,
    "Verdict instability (fraction of iterations where verdict changed)")
  expect_true(any(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomBar"),
    logical(1)
  )))
})

test_that("plot_weight_sensitivity renders single-result density semantics", {
  skip_if_not_installed("ggplot2")

  sensitivity <- list(
    iterations = data.frame(
      composite_score = c(0.20, 0.45, 0.82, 0.68),
      verdict = c("Poor", "Moderate", "Good", "Good"),
      stringsAsFactors = FALSE
    ),
    summary = list(
      verdict_pct = c(Good = 0.5, Moderate = 0.25, Poor = 0.25),
      composite_ci = c(0.20, 0.82),
      reference_composite = 0.58
    )
  )
  p <- plot_weight_sensitivity(sensitivity, title = "Single sensitivity")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 4L)
  expect_identical(p$labels$title, "Single sensitivity")
  expect_true(any(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomRug")
  }, logical(1L))))
  expect_true(any(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomVline")
  }, logical(1L))))
})

test_that("plot_weight_sensitivity handles a moderate-only iteration set", {
  skip_if_not_installed("ggplot2")

  sensitivity <- list(
    iterations = data.frame(
      composite_score = c(0.42, 0.48, 0.55),
      verdict = c("Moderate", "Moderate", "Good"),
      stringsAsFactors = FALSE
    ),
    summary = list(
      verdict_pct = c(Good = 1 / 3, Moderate = 2 / 3, Poor = 0),
      composite_ci = c(0.42, 0.55),
      reference_composite = 0.48
    )
  )
  p <- plot_weight_sensitivity(sensitivity)
  expect_s3_class(p, "ggplot")
  expect_equal(sum(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomDensity")
  }, logical(1L))), 2L)
})

test_that("plot_weight_sensitivity rejects malformed sensitivity results", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_weight_sensitivity("not a sensitivity result"),
    class = "pepvet_error_invalid_input"
  )

  expect_error(
    plot_weight_sensitivity(
      list(iterations = tibble::tibble(), summary = list())
    ),
    class = "pepvet_error_invalid_input"
  )

  malformed_batch <- list(
    per_protein = tibble::tibble(
      protein_id = "protein_1",
      verdict_instability = 1.1
    ),
    summary = list(total_instability = 0.5)
  )
  expect_error(
    plot_weight_sensitivity(malformed_batch),
    class = "pepvet_error_invalid_input"
  )

  missing_batch_column <- list(
    per_protein = tibble::tibble(protein_id = "protein_1"),
    summary = list(total_instability = 0.5)
  )
  expect_error(
    plot_weight_sensitivity(missing_batch_column),
    class = "pepvet_error_invalid_input"
  )

  invalid_batch_enzyme <- malformed_batch
  invalid_batch_enzyme$per_protein$verdict_instability <- 0.5
  invalid_batch_enzyme$per_protein$enzyme <- NA_character_
  expect_error(
    plot_weight_sensitivity(invalid_batch_enzyme),
    class = "pepvet_error_invalid_input"
  )

  invalid_batch_summary <- malformed_batch
  invalid_batch_summary$per_protein$verdict_instability <- 0.5
  invalid_batch_summary$summary$total_instability <- Inf
  expect_error(
    plot_weight_sensitivity(invalid_batch_summary),
    class = "pepvet_error_invalid_input"
  )

  valid_single <- sensitivity_analysis(.fix_bsa_trypsin, n_iter = 5L)
  invalid_single_values <- valid_single
  invalid_single_values$iterations$composite_score[[1L]] <- Inf
  expect_error(
    plot_weight_sensitivity(invalid_single_values),
    class = "pepvet_error_invalid_input"
  )

  invalid_single_summary <- valid_single
  invalid_single_summary$summary$composite_ci <- c(0.9, 0.1)
  expect_error(
    plot_weight_sensitivity(invalid_single_summary),
    class = "pepvet_error_invalid_input"
  )
})

test_that("plot_weight_sensitivity facets batch enzyme results", {
  skip_if_not_installed("ggplot2")
  sequences <- Biostrings::AAStringSet(c(
    alpha = "AKAAAAAAK",
    beta = "AKRTPK"
  ))
  comparison <- suppressMessages(batch_compare_enzymes(
    sequences,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 0L
  ))
  sensitivity <- sensitivity_analysis(
    comparison,
    n_iter = 20L,
    chunk_size = 2L
  )
  plot <- plot_weight_sensitivity(sensitivity)

  expect_s3_class(plot, "ggplot")
  expect_true("enzyme" %in% names(plot$data))
  expect_true("enzyme_wrapped" %in% names(plot$data))
  expect_true(any(vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomText") && is.data.frame(layer$data) &&
      "label" %in% names(layer$data)
  }, logical(1L))))
})
