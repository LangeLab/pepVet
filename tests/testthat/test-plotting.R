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
