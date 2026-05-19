# ── test-plotting.R ───────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────

# Shared fixtures are pre-computed in helper-fixtures.R (sourced once).
# .bsa_result() / .h3_result() are thin wrappers that return the cached
# fixture, falling back to live computation for non-default MC levels.
.bsa_result <- function(mc = 0L) {
  if (mc == 0L) {
    return(.fix_bsa_mc0)
  }
  if (mc == 1L) {
    return(.fix_bsa_mc1)
  }
  if (mc == 2L) {
    return(.fix_bsa_mc2)
  }
  evaluate_digest(.bsa_path, enzyme = "trypsin", missed_cleavages = mc)
}

.h3_result <- function(mc = 0L) {
  if (mc == 0L) {
    return(.fix_h3_trypsin)
  }
  evaluate_digest(.h3_path, enzyme = "trypsin", missed_cleavages = mc)
}

# ── plot_digest_profile ───────────────────────────────────────────────────────

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
    "extdata", "small_proteome_50_proteins.fasta",
    package = "pepVet"
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

# ── plot_coverage_map ─────────────────────────────────────────────────────────

# Shared fixture helpers
.bsa_cs <- function() {
  annotate_cleavage_sites(.bsa_path, enzyme = "trypsin")
}

.h3_cs <- function() {
  annotate_cleavage_sites(.h3_path, enzyme = "trypsin")
}

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

  # chymotrypsin-high yields shorter, more numerous peptides than trypsin;
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


# ── plot_enzyme_comparison ────────────────────────────────────────────────────

# Shared fixture helper
.bsa_comparison <- function(
  enzymes = c(
    "trypsin", "lysc",
    "glutamyl endopeptidase"
  )
) {
  compare_digests(.bsa_path, enzymes = enzymes)
}

test_that("plot_enzyme_comparison returns a patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp)
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
  p <- plot_enzyme_comparison(comp, title = "My comparison")
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

  bad <- data.frame(
    enzyme = "trypsin", composite_score = 0.6,
    S_coverage = 0.8, S_length = 0.7,
    S_count = 0.5, S_hydro = 0.6, S_charge = 0.4
  )
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
})


# ── plot_length_distribution ──────────────────────────────────────────────────

test_that("plot_length_distribution returns a ggplot from evaluate_digest result", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_length_distribution(res)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_length_distribution accepts a bare peptide data.frame", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
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
  p <- plot_length_distribution(res, title = "My length plot")
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


# ── plot_gravy_landscape ──────────────────────────────────────────────────────

test_that("plot_gravy_landscape returns a patchwork/ggplot from evaluate_digest result", {
  res <- .fix_bsa_trypsin
  p <- plot_gravy_landscape(res)
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
  res <- .fix_bsa_trypsin
  expect_no_error(
    plot_gravy_landscape(res, length_range = c(5L, 30L), gravy_range = c(-1.5, 1.0))
  )
})

test_that("plot_gravy_landscape: label_outliers_n = 0L suppresses labels silently", {
  res <- .fix_bsa_trypsin
  expect_no_error(plot_gravy_landscape(res, label_outliers_n = 0L))
})

test_that("plot_gravy_landscape: custom title is accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_gravy_landscape(res, title = "My GRAVY plot")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("plot_gravy_landscape: H3.1 (short protein) renders without warnings", {
  res <- .fix_h3_trypsin
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


# ── plot_pI_distribution ──────────────────────────────────────────────────────

test_that("plot_pI_distribution returns a ggplot from evaluate_digest result", {
  res <- .fix_bsa_trypsin
  p <- plot_pI_distribution(res)
  expect_s3_class(p, "gg")
})

test_that("plot_pI_distribution accepts score_peptides(include_pI=TRUE) tibble", {
  res <- .fix_bsa_trypsin
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
  res <- .fix_bsa_trypsin
  expect_no_error(plot_pI_distribution(res, show_fraction_lines = FALSE))
})

test_that("plot_pI_distribution: custom fraction_breaks accepted", {
  res <- .fix_bsa_trypsin
  expect_no_error(plot_pI_distribution(res, fraction_breaks = c(4, 6, 8, 10)))
})

test_that("plot_pI_distribution: custom title is accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_pI_distribution(res, title = "My pI plot")
  expect_s3_class(p, "gg")
})

test_that("plot_pI_distribution: H3.1 renders without warnings", {
  res <- .fix_h3_trypsin
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


# ── multi-input: plot_length_distribution ────────────────────────────────────

test_that("plot_length_distribution: named list produces faceted ggplot", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_length_distribution(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_length_distribution: multi-input renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  expect_no_warning(
    plot_length_distribution(list(BSA = bsa_res, H3 = h3_res))
  )
})

# ── multi-input: plot_gravy_landscape ────────────────────────────────────────

test_that("plot_gravy_landscape: named list produces faceted ggplot", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_gravy_landscape(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_gravy_landscape: multi-input renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  expect_no_warning(
    plot_gravy_landscape(list(BSA = bsa_res, H3 = h3_res))
  )
})

# ── multi-input: plot_pI_distribution ────────────────────────────────────────

test_that("plot_pI_distribution: named list produces overlaid density ggplot", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_pI_distribution(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_pI_distribution: multi-input renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  expect_no_warning(
    plot_pI_distribution(list(BSA = bsa_res, H3 = h3_res))
  )
})

# ── plot_protein_comparison ───────────────────────────────────────────────────

test_that("plot_protein_comparison returns ggplot from named list", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
})

test_that("plot_protein_comparison accepts batch_evaluate tibble", {
  skip_if_not_installed("Biostrings")
  batch <- .fix_batch_trypsin
  p <- plot_protein_comparison(batch)
  expect_s3_class(p, "gg")
})

test_that("plot_protein_comparison: show_verdict = FALSE renders without error", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  expect_no_error(
    plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res),
      show_verdict = FALSE
    )
  )
})

test_that("plot_protein_comparison: custom title accepted", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_protein_comparison(list(BSA = bsa_res, H3 = h3_res),
    title = "My comparison"
  )
  expect_s3_class(p, "gg")
})

test_that("plot_protein_comparison: renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
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


# ── plot_cleavage_map ─────────────────────────────────────────────────────────

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
  # asp-n cuts N-terminal to D; exercises site-inference from the opposite
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

# ── plot_enzyme_protein_heatmap ───────────────────────────────────────────────

test_that("plot_enzyme_protein_heatmap returns ggplot from nested list", {
  nested <- list(
    BSA = list(
      trypsin = .fix_bsa_trypsin,
      lysc    = .fix_bsa_lysc
    ),
    H3 = list(
      trypsin = .fix_h3_trypsin,
      lysc    = .fix_h3_lysc
    )
  )
  p <- plot_enzyme_protein_heatmap(nested)
  expect_s3_class(p, "gg")
})

test_that("plot_enzyme_protein_heatmap: component argument works", {
  nested <- list(
    BSA = list(trypsin = .fix_bsa_trypsin),
    H3  = list(trypsin = .fix_h3_trypsin)
  )
  p <- plot_enzyme_protein_heatmap(nested, component = "S_length")
  expect_s3_class(p, "gg")
})

test_that("plot_enzyme_protein_heatmap renders without warnings", {
  nested <- list(
    BSA = list(
      trypsin = .fix_bsa_trypsin,
      lysc = .fix_bsa_lysc
    ),
    H3 = list(
      trypsin = .fix_h3_trypsin,
      lysc = .fix_h3_lysc
    )
  )
  expect_no_warning(plot_enzyme_protein_heatmap(nested))
})

test_that("plot_enzyme_protein_heatmap errors on bad component", {
  nested <- list(BSA = list(trypsin = .fix_bsa_trypsin))
  expect_error(
    plot_enzyme_protein_heatmap(nested, component = "not_a_score"),
    class = "pepvet_error_invalid_component"
  )
})

# ── plot_missed_cleavage_impact ───────────────────────────────────────────────

test_that("plot_missed_cleavage_impact returns ggplot from MC list", {
  r0 <- .fix_bsa_mc0
  r1 <- .fix_bsa_mc1
  r2 <- .fix_bsa_mc2
  p <- plot_missed_cleavage_impact(list("MC=0" = r0, "MC=1" = r1, "MC=2" = r2))
  expect_s3_class(p, "gg")
})

test_that("plot_missed_cleavage_impact: unnamed list auto-named", {
  r0 <- .fix_bsa_mc0
  r1 <- .fix_bsa_mc1
  p <- plot_missed_cleavage_impact(list(r0, r1))
  expect_s3_class(p, "gg")
})

test_that("plot_missed_cleavage_impact renders without warnings", {
  r0 <- .fix_bsa_mc0
  r1 <- .fix_bsa_mc1
  r2 <- .fix_bsa_mc2
  expect_no_warning(
    plot_missed_cleavage_impact(list("MC=0" = r0, "MC=1" = r1, "MC=2" = r2))
  )
})

test_that("plot_missed_cleavage_impact: chymotrypsin-high MC=0/1 renders", {
  # chymotrypsin produces many more peptides than trypsin; tests the line-plot
  # rendering when S_count and S_length change more dramatically across MC levels
  r0 <- .fix_bsa_chymotryp
  r1 <- .fix_bsa_chymotryp_mc1
  p <- plot_missed_cleavage_impact(list("MC=0" = r0, "MC=1" = r1))
  expect_s3_class(p, "gg")
})

test_that("plot_missed_cleavage_impact errors on single-element list", {
  r0 <- .fix_bsa_mc0
  expect_error(
    plot_missed_cleavage_impact(list(r0)),
    class = "pepvet_error_invalid_digest_result"
  )
})

# ── plot_batch_summary ────────────────────────────────────────────────────────

test_that("plot_batch_summary returns patchwork from batch_evaluate", {
  skip_if_not_installed("Biostrings")
  batch <- .fix_batch_trypsin
  p <- plot_batch_summary(batch)
  expect_s3_class(p, "patchwork")
})

test_that("plot_batch_summary renders without warnings", {
  skip_if_not_installed("Biostrings")
  batch <- .fix_batch_trypsin
  expect_no_warning(plot_batch_summary(batch))
})

test_that("plot_batch_summary: chymotrypsin-high batch renders", {
  skip_if_not_installed("Biostrings")
  # chymotrypsin yields higher peptide counts and lower scores than trypsin;
  # exercises the histogram color-coding and LOESS scatter with a different
  # verdict distribution
  batch <- .fix_batch_chymotryp
  p <- plot_batch_summary(batch)
  expect_s3_class(p, "patchwork")
})

test_that("plot_batch_summary errors on missing required columns", {
  bad_batch <- data.frame(x = 1:3)
  expect_error(
    plot_batch_summary(bad_batch),
    class = "pepvet_error_invalid_batch"
  )
})

# ── plot_component_scatter ────────────────────────────────────────────────────

test_that("plot_component_scatter returns ggplot from batch_evaluate", {
  skip_if_not_installed("Biostrings")
  batch <- .fix_batch_trypsin
  p <- plot_component_scatter(batch)
  expect_s3_class(p, "gg")
})

test_that("plot_component_scatter: custom component axes work", {
  skip_if_not_installed("Biostrings")
  batch <- .fix_batch_trypsin
  p <- plot_component_scatter(batch,
    x_component = "S_length",
    y_component = "S_coverage"
  )
  expect_s3_class(p, "gg")
})

test_that("plot_component_scatter renders without warnings", {
  skip_if_not_installed("Biostrings")
  batch <- .fix_batch_trypsin
  expect_no_warning(plot_component_scatter(batch))
})

test_that("plot_component_scatter: chymotrypsin-high batch renders", {
  skip_if_not_installed("Biostrings")
  # chymotrypsin batch has a distinctive score cluster; exercises label
  # rendering for Poor proteins and rug marks with a different distribution
  batch <- .fix_batch_chymotryp
  p <- plot_component_scatter(batch)
  expect_s3_class(p, "gg")
})

test_that("plot_component_scatter errors on invalid component", {
  skip_if_not_installed("Biostrings")
  batch <- batch_evaluate(
    Biostrings::readAAStringSet(.bsa_path),
    enzyme = "trypsin"
  )
  expect_error(
    plot_component_scatter(batch, x_component = "not_a_col"),
    class = "pepvet_error_invalid_component"
  )
})
