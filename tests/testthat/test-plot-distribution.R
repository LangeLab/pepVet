# Tests aligned with the corresponding plotting source domain.
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


# plot_gravy_landscape.

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


# plot_pI_distribution.

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


# multi-input: plot_length_distribution.

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

# multi-input: plot_gravy_landscape.

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

# multi-input: plot_pI_distribution.

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



# Restored section from the original plotting test surface.

# Restored section from the original plotting test surface.
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
  # Chymotrypsin produces many more peptides than trypsin. Test the line plot
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

test_that("plot_mz_distribution returns ggplot from evaluate_digest", {
  p <- plot_mz_distribution(.fix_bsa_trypsin)
  expect_s3_class(p, "gg")
})

test_that("plot_mz_distribution: show_rug = TRUE and FALSE both work", {
  p1 <- plot_mz_distribution(.fix_bsa_trypsin, show_rug = TRUE)
  p2 <- plot_mz_distribution(.fix_bsa_trypsin, show_rug = FALSE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
})

test_that("plot_mz_distribution: custom scan_range", {
  p <- plot_mz_distribution(.fix_bsa_trypsin, scan_range = c(400, 1000))
  expect_s3_class(p, "gg")
})

test_that("plot_mz_distribution: multi-input mode", {
  p <- plot_mz_distribution(
    list(BSA = .fix_bsa_trypsin, H3 = .fix_h3_trypsin)
  )
  expect_s3_class(p, "gg")
})
