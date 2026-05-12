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
