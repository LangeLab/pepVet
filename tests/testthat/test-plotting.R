# ── test-plotting.R ───────────────────────────────────────────────────────────
# vdiffr snapshot tests for pepVet visualization functions.
# Run `vdiffr::manage_cases()` to review new / changed snapshots interactively.
# ─────────────────────────────────────────────────────────────────────────────

# Helper shared across tests in this file
.bsa_result <- function(mc = 0L) {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  evaluate_digest(bsa_path, enzyme = "trypsin", missed_cleavages = mc)
}

.h3_result <- function(mc = 0L) {
  h3_path <- system.file("extdata", "P68431.fasta", package = "pepVet")
  evaluate_digest(h3_path, enzyme = "trypsin", missed_cleavages = mc)
}

# ── plot_digest_profile ───────────────────────────────────────────────────────

test_that("plot_digest_profile returns a patchwork object for BSA / trypsin", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  res <- .bsa_result()
  p   <- plot_digest_profile(res)

  expect_s3_class(p, "patchwork")
})

test_that("plot_digest_profile snapshot: BSA / trypsin (standard)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("vdiffr")

  res <- .bsa_result()
  vdiffr::expect_doppelganger(
    "digest_profile_bsa_trypsin",
    plot_digest_profile(res, title = "BSA – trypsin")
  )
})

test_that("plot_digest_profile snapshot: Histone H3.1 / trypsin (difficult)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("vdiffr")

  res <- .h3_result()
  vdiffr::expect_doppelganger(
    "digest_profile_h3_trypsin",
    plot_digest_profile(res, title = "Histone H3.1 – trypsin")
  )
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

  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  proteome_path <- system.file(
    "extdata", "small_proteome_50_proteins.fasta", package = "pepVet"
  )
  proteome_digest <- digest_protein(proteome_path, enzyme = "trypsin")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin", proteome = proteome_digest)

  p <- plot_digest_profile(res)
  expect_s3_class(p, "patchwork")
})

# ── Input validation ──────────────────────────────────────────────────────────

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
    "extdata", "small_proteome_50_proteins.fasta", package = "pepVet"
  )
  # compare_digests returns a flat tibble, not an evaluate_digest list;
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

  res          <- .bsa_result()
  bad_peps     <- res$peptides[, c("protein_id", "peptide"), drop = FALSE]
  bad_result   <- list(
    scores   = res$scores,
    peptides = bad_peps,
    params   = res$params
  )

  expect_error(
    plot_digest_profile(bad_result),
    class = "pepvet_error_invalid_digest_result"
  )
})

# ── plot_coverage_map ─────────────────────────────────────────────────────────

# Shared fixture helpers
.bsa_cs <- function() {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  annotate_cleavage_sites(bsa_path, enzyme = "trypsin")
}

.h3_cs <- function() {
  h3_path <- system.file("extdata", "P68431.fasta", package = "pepVet")
  annotate_cleavage_sites(h3_path, enzyme = "trypsin")
}

test_that("plot_coverage_map returns a ggplot for BSA / trypsin MC=0", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p   <- plot_coverage_map(res)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_coverage_map returns a ggplot for MC=1 (multi-lane + packing)", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result(mc = 1L)
  p   <- plot_coverage_map(res)
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
  cs  <- .bsa_cs()
  expect_no_error(plot_coverage_map(res, cleavage_sites = cs))
})

test_that("plot_coverage_map: domains overlay renders without error", {
  skip_if_not_installed("ggplot2")

  res     <- .bsa_result()
  domains <- data.frame(
    name  = c("Domain A", "Domain B"),
    start = c(1L, 200L),
    end   = c(150L, 400L),
    stringsAsFactors = FALSE
  )
  expect_no_error(plot_coverage_map(res, domains = domains))
})

test_that("plot_coverage_map: custom title is accepted", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p   <- plot_coverage_map(res, title = "My custom title")
  expect_s3_class(p, "ggplot")
})

test_that("plot_coverage_map: H3.1 short protein renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result()
  cs  <- .h3_cs()
  expect_no_error(plot_coverage_map(res, cleavage_sites = cs))
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

  res   <- .bsa_result()
  extra <- res$peptides
  extra$protein_id <- "fake|P99999|FAKE_PROTEIN"
  res$peptides <- rbind(res$peptides, extra)

  expect_error(
    plot_coverage_map(res),
    class = "pepvet_error_multi_protein"
  )
})

test_that("plot_coverage_map snapshot: BSA / trypsin MC=1 with cleavage sites", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("vdiffr")

  res <- .bsa_result(mc = 1L)
  cs  <- .bsa_cs()
  vdiffr::expect_doppelganger(
    "coverage_map_bsa_trypsin_mc1",
    plot_coverage_map(res, cleavage_sites = cs,
                      title = "BSA – trypsin – MC=1")
  )
})

# ── plot_enzyme_comparison ────────────────────────────────────────────────────

# Shared fixture helper
.bsa_comparison <- function(enzymes = c("trypsin", "lysc",
                                        "glutamyl endopeptidase")) {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  compare_digests(bsa_path, enzymes = enzymes)
}

test_that("plot_enzyme_comparison returns a patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p    <- plot_enzyme_comparison(comp)
  expect_s3_class(p, "patchwork")
})

test_that("plot_enzyme_comparison: recommend = FALSE suppresses badge", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  expect_no_error(plot_enzyme_comparison(comp, recommend = FALSE))
})

test_that("plot_enzyme_comparison: subset of scores renders without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  expect_no_error(
    plot_enzyme_comparison(comp, scores = c("S_coverage", "S_count"))
  )
})

test_that("plot_enzyme_comparison: custom title is accepted", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p    <- plot_enzyme_comparison(comp, title = "My comparison")
  expect_s3_class(p, "patchwork")
})

test_that("plot_enzyme_comparison: two-enzyme comparison works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison(enzymes = c("trypsin", "lysc"))
  expect_no_error(plot_enzyme_comparison(comp))
})

test_that("plot_enzyme_comparison errors on non-data.frame input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  expect_error(
    plot_enzyme_comparison("not a tibble"),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_enzyme_comparison errors on missing required columns", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  bad <- data.frame(enzyme = c("a", "b"), foo = 1:2)
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_enzyme_comparison errors when fewer than 2 enzymes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  bad <- data.frame(enzyme = "trypsin", composite_score = 0.6,
                    S_coverage = 0.8, S_length = 0.7,
                    S_count = 0.5, S_hydro = 0.6, S_charge = 0.4)
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_enzyme_comparison snapshot: BSA 3-enzyme comparison", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("vdiffr")

  comp <- .bsa_comparison()
  vdiffr::expect_doppelganger(
    "enzyme_comparison_bsa_3",
    plot_enzyme_comparison(comp, title = "BSA – enzyme comparison")
  )
})

# ── plot_score_radar ──────────────────────────────────────────────────────────

test_that("plot_score_radar returns a ggplot from compare_digests tibble", {
  skip_if_not_installed("ggplot2")

  comp <- .bsa_comparison()
  p    <- plot_score_radar(comp)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_score_radar: single evaluate_digest() result renders", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  expect_no_error(plot_score_radar(res))
})

test_that("plot_score_radar: four enzymes renders without error", {
  skip_if_not_installed("ggplot2")

  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  comp4 <- compare_digests(bsa_path,
    enzymes = c("trypsin", "lysc", "glutamyl endopeptidase",
                "chymotrypsin-high"))
  expect_no_error(plot_score_radar(comp4))
})

test_that("plot_score_radar: subset of scores axes renders without error", {
  skip_if_not_installed("ggplot2")

  comp <- .bsa_comparison()
  expect_no_error(
    plot_score_radar(comp, scores = c("S_coverage", "S_count", "S_charge"))
  )
})

test_that("plot_score_radar: custom title is accepted", {
  skip_if_not_installed("ggplot2")

  comp <- .bsa_comparison()
  p    <- plot_score_radar(comp, title = "My radar")
  expect_s3_class(p, "ggplot")
})

test_that("plot_score_radar errors on non-data.frame / non-list input", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_score_radar("not valid"),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_score_radar errors on missing required columns", {
  skip_if_not_installed("ggplot2")

  bad <- data.frame(enzyme = c("a", "b"), foo = 1:2)
  expect_error(
    plot_score_radar(bad),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_score_radar snapshot: BSA 3-enzyme radar", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("vdiffr")

  comp <- .bsa_comparison()
  vdiffr::expect_doppelganger(
    "score_radar_bsa_3",
    plot_score_radar(comp, title = "BSA – score radar")
  )
})

# ── plot_length_distribution ──────────────────────────────────────────────────

test_that("plot_length_distribution returns a ggplot from evaluate_digest result", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p   <- plot_length_distribution(res)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_length_distribution accepts a bare peptide data.frame", {
  skip_if_not_installed("ggplot2")

  res  <- .bsa_result()
  peps <- res$peptides
  expect_no_error(plot_length_distribution(peps))
})

test_that("plot_length_distribution: show_density = FALSE renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  expect_no_error(plot_length_distribution(res, show_density = FALSE))
})

test_that("plot_length_distribution: custom length_range is respected", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  expect_no_error(plot_length_distribution(res, length_range = c(5L, 30L)))
})

test_that("plot_length_distribution: custom title is accepted", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p   <- plot_length_distribution(res, title = "My length plot")
  expect_s3_class(p, "ggplot")
})

test_that("plot_length_distribution: H3.1 (short protein) renders without warnings", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result()
  expect_no_warning(plot_length_distribution(res))
})

test_that("plot_length_distribution errors on invalid input", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_length_distribution("not valid"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_length_distribution errors when length column missing", {
  skip_if_not_installed("ggplot2")

  bad <- data.frame(x = 1:5)
  expect_error(
    plot_length_distribution(bad),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_length_distribution snapshot: BSA trypsin", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("vdiffr")

  res <- .bsa_result()
  vdiffr::expect_doppelganger(
    "length_dist_bsa_trypsin",
    plot_length_distribution(res, title = "BSA – trypsin length distribution")
  )
})


# ── plot_gravy_landscape ──────────────────────────────────────────────────────

test_that("plot_gravy_landscape returns a patchwork/ggplot from evaluate_digest result", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  p   <- plot_gravy_landscape(res)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("plot_gravy_landscape accepts a bare data.frame with length + gravy", {
  peps <- data.frame(
    length = c(8L, 12L, 30L, 5L),
    gravy  = c(0.2, -0.5, 0.8, -1.2)
  )
  expect_no_error(plot_gravy_landscape(peps))
})

test_that("plot_gravy_landscape auto-computes GRAVY from peptide column", {
  peps <- data.frame(
    length  = c(8L, 12L),
    peptide = c("ACDEFGHIK", "LMNPQRSTVW")
  )
  expect_no_error(plot_gravy_landscape(peps))
})

test_that("plot_gravy_landscape: custom length_range and gravy_range accepted", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  expect_no_error(
    plot_gravy_landscape(res, length_range = c(5L, 30L), gravy_range = c(-1.5, 1.0))
  )
})

test_that("plot_gravy_landscape: label_outliers_n = 0L suppresses labels silently", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  expect_no_error(plot_gravy_landscape(res, label_outliers_n = 0L))
})

test_that("plot_gravy_landscape: custom title is accepted", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  p   <- plot_gravy_landscape(res, title = "My GRAVY plot")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("plot_gravy_landscape: H3.1 (short protein) renders without warnings", {
  h3_path <- system.file("extdata", "P68431.fasta", package = "pepVet")
  res <- evaluate_digest(h3_path, enzyme = "trypsin")
  expect_no_warning(plot_gravy_landscape(res))
})

test_that("plot_gravy_landscape errors on invalid input", {
  expect_error(
    plot_gravy_landscape("not valid"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_gravy_landscape errors when length column is missing", {
  bad <- data.frame(gravy = c(0.1, 0.2))
  expect_error(
    plot_gravy_landscape(bad),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_gravy_landscape snapshot: BSA trypsin", {
  skip_if_not_installed("vdiffr")
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  vdiffr::expect_doppelganger(
    "gravy_landscape_bsa_trypsin",
    plot_gravy_landscape(res, title = "BSA – trypsin GRAVY landscape")
  )
})


# ── plot_pI_distribution ──────────────────────────────────────────────────────

test_that("plot_pI_distribution returns a ggplot from evaluate_digest result", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  p   <- plot_pI_distribution(res)
  expect_s3_class(p, "gg")
})

test_that("plot_pI_distribution accepts score_peptides(include_pI=TRUE) tibble", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res    <- evaluate_digest(bsa_path, enzyme = "trypsin")
  scored <- score_peptides(res$peptides, enzyme = "trypsin", include_pI = TRUE)
  expect_no_error(plot_pI_distribution(scored))
})

test_that("plot_pI_distribution accepts data.frame with numeric pI column", {
  df <- data.frame(pI = c(4.2, 5.7, 6.1, 8.3, 9.5))
  expect_no_error(plot_pI_distribution(df))
})

test_that("plot_pI_distribution accepts a bare numeric vector", {
  expect_no_error(plot_pI_distribution(c(4.2, 5.7, 6.1, 8.3, 9.5)))
})

test_that("plot_pI_distribution: show_fraction_lines = FALSE renders without error", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  expect_no_error(plot_pI_distribution(res, show_fraction_lines = FALSE))
})

test_that("plot_pI_distribution: custom fraction_breaks accepted", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  expect_no_error(plot_pI_distribution(res, fraction_breaks = c(4, 6, 8, 10)))
})

test_that("plot_pI_distribution: custom title is accepted", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  p   <- plot_pI_distribution(res, title = "My pI plot")
  expect_s3_class(p, "gg")
})

test_that("plot_pI_distribution: H3.1 renders without warnings", {
  h3_path <- system.file("extdata", "P68431.fasta", package = "pepVet")
  res <- evaluate_digest(h3_path, enzyme = "trypsin")
  expect_no_warning(plot_pI_distribution(res))
})

test_that("plot_pI_distribution errors on invalid input", {
  expect_error(
    plot_pI_distribution("not valid"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_pI_distribution errors when data.frame lacks pI column", {
  bad <- data.frame(length = c(8L, 12L))
  expect_error(
    plot_pI_distribution(bad),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_pI_distribution snapshot: BSA trypsin", {
  skip_if_not_installed("vdiffr")
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  res <- evaluate_digest(bsa_path, enzyme = "trypsin")
  vdiffr::expect_doppelganger(
    "pi_dist_bsa_trypsin",
    plot_pI_distribution(res, title = "BSA – trypsin pI distribution")
  )
})


# ── multi-input: plot_length_distribution ────────────────────────────────────

test_that("plot_length_distribution: named list produces faceted ggplot", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  p <- plot_length_distribution(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_length_distribution: multi-input renders without warnings", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  expect_no_warning(
    plot_length_distribution(list(BSA = bsa_res, H3 = h3_res))
  )
})

# ── multi-input: plot_gravy_landscape ────────────────────────────────────────

test_that("plot_gravy_landscape: named list produces faceted ggplot", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  p <- plot_gravy_landscape(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_gravy_landscape: multi-input renders without warnings", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  expect_no_warning(
    plot_gravy_landscape(list(BSA = bsa_res, H3 = h3_res))
  )
})

# ── multi-input: plot_pI_distribution ────────────────────────────────────────

test_that("plot_pI_distribution: named list produces overlaid density ggplot", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  p <- plot_pI_distribution(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_pI_distribution: multi-input renders without warnings", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  expect_no_warning(
    plot_pI_distribution(list(BSA = bsa_res, H3 = h3_res))
  )
})

# ── plot_protein_comparison ───────────────────────────────────────────────────

test_that("plot_protein_comparison returns ggplot from named list", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  p <- plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_protein_comparison accepts batch_evaluate tibble", {
  skip_if_not_installed("Biostrings")
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_seq  <- Biostrings::readAAStringSet(bsa_path)
  h3_seq   <- Biostrings::readAAStringSet(h3_path)
  batch    <- batch_evaluate(c(bsa_seq, h3_seq), enzyme = "trypsin")
  p <- plot_protein_comparison(batch)
  expect_s3_class(p, "gg")
})

test_that("plot_protein_comparison: show_verdict = FALSE renders without error", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  expect_no_error(
    plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res),
                            show_verdict = FALSE)
  )
})

test_that("plot_protein_comparison: custom title accepted", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  p <- plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res),
                                title = "My comparison")
  expect_s3_class(p, "gg")
})

test_that("plot_protein_comparison: renders without warnings", {
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  expect_no_warning(
    plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res))
  )
})

test_that("plot_protein_comparison errors on invalid input", {
  expect_error(
    plot_protein_comparison("not valid"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_protein_comparison snapshot: BSA vs H3 trypsin", {
  skip_if_not_installed("vdiffr")
  bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
  h3_path  <- system.file("extdata", "P68431.fasta", package = "pepVet")
  bsa_res  <- evaluate_digest(bsa_path, enzyme = "trypsin")
  h3_res   <- evaluate_digest(h3_path,  enzyme = "trypsin")
  vdiffr::expect_doppelganger(
    "protein_comparison_bsa_h3_trypsin",
    plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res),
                            title = "BSA vs H3.1 — trypsin")
  )
})
