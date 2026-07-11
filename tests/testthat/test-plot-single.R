# Tests aligned with the corresponding plotting source domain.
test_that("plot_digest_profile returns a patchwork object for BSA / trypsin", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  p <- plot_digest_profile(res)

  expect_s3_class(p, "patchwork")
})


test_that("plot_digest_profile: custom title and length/gravy ranges", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  expect_no_error(
    plot_digest_profile(
      res,
      length_range = c(6L, 30L),
      gravy_range  = c(-1.5, 1.0),
      title        = "Custom ranges"
    )
  )
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
})

test_that("plot_coverage_map returns a ggplot for MC=1 (multi-lane + packing)", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  p <- plot_coverage_map(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_coverage_map: color_by = 'length_class' renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 2L)
  expect_no_error(plot_coverage_map(res, color_by = "length_class"))
})

test_that("plot_coverage_map: color_by = 'hydrophobicity' renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  expect_no_error(plot_coverage_map(res, color_by = "hydrophobicity"))
})

test_that("plot_coverage_map: cleavage_sites overlay renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  cs <- .bsa_cs()
  expect_no_error(plot_coverage_map(res, cleavage_sites = cs))
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
  expect_no_error(plot_coverage_map(res, domains = domains))
})

test_that("plot_coverage_map: custom title is accepted", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_coverage_map(res, title = "My custom title")
  expect_s3_class(p, "ggplot")
})

test_that("plot_coverage_map: H3.1 short protein renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result()
  cs <- .h3_cs()
  expect_no_error(plot_coverage_map(res, cleavage_sites = cs))
})

test_that("plot_coverage_map: chymotrypsin-high produces different lane layout", {
  skip_if_not_installed("ggplot2")

  # Chymotrypsin-high yields shorter, more numerous peptides than trypsin.
  # exercises the lane-packing algorithm with a denser peptide map
  res <- .fix_bsa_chymotryp
  p <- plot_coverage_map(res)
  expect_s3_class(p, "gg")
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

# Restored section from the original plotting test surface.

# Restored section from the original plotting test surface.
test_that("plot_peptide_overlap_map returns a ggplot for BSA / trypsin MC=1", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  p <- plot_peptide_overlap_map(res)

  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_peptide_overlap_map supports wrapped rows and all-peptide mode", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result(mc = 1L)

  expect_no_error(
    plot_peptide_overlap_map(
      res,
      length_range = NULL,
      residues_per_line = 20L
    )
  )
})

test_that("plot_cleavage_map returns ggplot from evaluate_digest result", {
  res <- .fix_bsa_trypsin
  p <- plot_cleavage_map(res)
  expect_s3_class(p, "gg")
})

test_that("plot_cleavage_map renders without warnings (no cleavage_sites)", {
  res <- .fix_bsa_trypsin
  expect_no_warning(plot_cleavage_map(res))
})

test_that("plot_cleavage_map: custom title accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_cleavage_map(res, title = "BSA cleavage")
  expect_s3_class(p, "gg")
})

test_that("plot_cleavage_map: asp-n endopeptidase (N-terminal cutter) renders", {
  # Asp-N cuts N-terminal to D. This exercises site inference from the opposite
  # end of peptides compared to trypsin/lysc
  res <- .fix_bsa_aspn
  p <- plot_cleavage_map(res)
  expect_s3_class(p, "gg")
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
  expect_true(any(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomBar"),
    logical(1)
  )))
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
})
