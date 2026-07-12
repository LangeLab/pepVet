# Tests aligned with the corresponding plotting source domain.
test_that("plot utility classifiers preserve inclusive boundaries", {
  classes <- .classify_length(c(6, 7, 25, 26), c(7L, 25L))
  expect_identical(as.character(classes),
    c("Too short", "Valid", "Valid", "Too long"))
  expect_identical(levels(classes), c("Valid", "Too short", "Too long"))
  expect_equal(.nice_x_step(597L), 50L)
})

test_that("tidy protein labels handle supported and fallback headers", {
  expect_identical(
    .tidy_protein_id("sp|P02769|ALBU_BOVIN Albumin"),
    "P02769  (ALBU_BOVIN)"
  )
  expect_identical(.tidy_protein_id("NP_001234.1 protein"), "NP_001234.1")
  expect_identical(.tidy_protein_id("generic_id description"), "generic_id")
  expect_identical(.tidy_protein_id(""), "")
  long_id <- paste(rep("X", 50L), collapse = "")
  expect_identical(.tidy_protein_id(long_id), paste0(strrep("X", 39L), "..."))
})

test_that("digest plot validation protects score and coordinate consumers", {
  res <- .fix_bsa_trypsin
  expect_null(.validate_digest_result_for_plot(res))

  bad <- res
  bad$scores$composite_score[[1L]] <- Inf
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$scores$S_unique <- Inf
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$scores$protein_id[[1L]] <- "other"
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$params$protein_ids[[1L]] <- "other"
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  empty_scores <- res
  empty_scores$scores <- res$scores[0, , drop = FALSE]
  expect_error(
    .validate_digest_result_for_plot(empty_scores),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$peptides$end[[1L]] <- bad$peptides$end[[1L]] + 1L
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$peptides$end[[1L]] <- .Machine$integer.max + 1
  expect_no_warning(
    expect_error(
      .validate_digest_result_for_plot(bad),
      class = "pepvet_error_invalid_digest_result"
    )
  )
  bad <- res
  bad$params$enzyme <- ""
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$params$missed_cleavages <- 1e20
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$scores$verdict[[1L]] <- "Unknown"
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$peptides$protein_id[[1L]] <- ""
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$peptides$peptide[[1L]] <- paste0(
    "Z",
    substr(
      bad$peptides$peptide[[1L]],
      2L,
      nchar(bad$peptides$peptide[[1L]])
    )
  )
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  bad <- res
  bad$peptides$missed_cleavages[[1L]] <- 0.5
  expect_error(
    .validate_digest_result_for_plot(bad),
    class = "pepvet_error_invalid_digest_result"
  )
  empty <- res
  empty$peptides <- empty$peptides[0, , drop = FALSE]
  expect_error(
    .validate_digest_result_for_plot(empty),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("coverage and sequence helpers use coordinate semantics", {
  peps <- data.frame(
    protein_id = rep("synthetic", 3L),
    peptide = c("ABC", "CDE", "FG"),
    start = c(1L, 3L, 6L),
    end = c(3L, 5L, 7L),
    length = c(3L, 3L, 2L),
    missed_cleavages = c(0L, 1L, 0L),
    stringsAsFactors = FALSE
  )
  stats <- .compute_coverage_stats(peps, 7L, c(2L, 3L), mc_filter = 0L)
  expect_identical(stats$covered, c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(stats$pct_cov, 71.4)
  expect_equal(stats$gap_df$xmin, 4L)
  expect_equal(stats$gap_df$xmax, 5L)
  no_mc_stats <- .compute_coverage_stats(
    peps[, setdiff(names(peps), "missed_cleavages"), drop = FALSE],
    7L, c(2L, 3L)
  )
  expect_identical(no_mc_stats$pct_cov, 100)
  reconstruct_peps <- peps
  reconstruct_peps$missed_cleavages <- 0L
  expect_identical(
    .reconstruct_sequence_from_peptides(reconstruct_peps),
    "ABCDEFG"
  )

  mismatch <- reconstruct_peps
  mismatch$peptide[[2L]] <- "XXE"
  expect_error(
    .reconstruct_sequence_from_peptides(mismatch),
    class = "pepvet_error_invalid_digest_result"
  )

  no_mc <- peps
  no_mc$missed_cleavages <- NULL
  expect_identical(.reconstruct_sequence_from_peptides(no_mc), "ABCDEFG")
  no_mc_empty <- peps
  no_mc_empty$missed_cleavages <- 1L
  expect_error(
    .reconstruct_sequence_from_peptides(no_mc_empty),
    class = "pepvet_error_invalid_digest_result"
  )
  incomplete <- peps[2:3, , drop = FALSE]
  expect_error(
    .reconstruct_sequence_from_peptides(incomplete),
    class = "pepvet_error_invalid_digest_result"
  )
  coordinate_mismatch <- peps[1L, , drop = FALSE]
  coordinate_mismatch$peptide <- "AB"
  expect_error(
    .reconstruct_sequence_from_peptides(coordinate_mismatch),
    class = "pepvet_error_invalid_digest_result"
  )

  sparse <- .fix_bsa_trypsin$peptides[c(1L, nrow(.fix_bsa_trypsin$peptides)), ,
    drop = FALSE
  ]
  coverage_panel <- .panel_coverage(
    sparse, protein_length = max(sparse$end), length_range = c(7L, 25L),
    missed_cleavages = 0L
  )
  expect_s3_class(coverage_panel, "ggplot")
  expect_match(coverage_panel$labels$caption, "uncovered")
  complete_panel <- .panel_coverage(
    peps, protein_length = 7L, length_range = c(2L, 3L)
  )
  expect_identical(complete_panel$labels$title, "Sequence Coverage")
  expect_null(complete_panel$labels$caption)
  gravy_panel <- .panel_gravy(
    data.frame(length = c(8L, 12L, 20L), gravy = c(0.2, 0.2, 0.2)),
    c(-1, 0.6)
  )
  expect_s3_class(gravy_panel, "ggplot")
})

test_that("layout helpers pack intervals and reserve tick space", {
  peps <- data.frame(start = c(1L, 2L, 5L), end = c(4L, 3L, 6L))
  packed <- .pack_peptides(peps)
  expect_identical(packed$track, c(1L, 2L, 1L))
  expect_true(all(packed$track >= 1L))
  empty <- .pack_peptides(peps[0, , drop = FALSE])
  expect_type(empty$track, "integer")
  expect_length(empty$track, 0L)

  lanes <- .lane_y_coords(c(0L, 1L, 2L), tick_height = 0.1)
  expect_equal(nrow(lanes), 3L)
  expect_true(all(lanes$y_lo < lanes$y_mid & lanes$y_mid < lanes$y_hi))
  expect_true(min(lanes$y_lo) >= 0.1)
  expect_true(max(lanes$y_hi) <= 1)
})

test_that("plot configuration validators enforce shape and relationships", {
  current <- .pepvet_pal
  updated <- .validate_plot_palette(
    list(verdict = c(Good = "#00AA00", Moderate = "#AAAA00", Poor = "#AA0000")),
    current
  )
  expect_identical(updated$verdict[["Good"]], "#00AA00")
  expect_error(
    .validate_plot_palette(list(verdict = c("#000000")), current),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    .validate_plot_palette(c(brand = "#000000"), current),
    class = "pepvet_error_invalid_config"
  )

  params <- .validate_plot_params(
    list(length_lo = 8, length_hi = 24, verdict_good = 0.7),
    .pepvet_params
  )
  expect_identical(params$length_lo, 8L)
  expect_identical(params$length_hi, 24L)
  expect_equal(params$verdict_good, 0.7)
  expect_error(
    .validate_plot_params(list(gravy_lo = 1, gravy_hi = -1), .pepvet_params),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    .validate_plot_params(c(verdict_good = 0.7), .pepvet_params),
    class = "pepvet_error_invalid_config"
  )
})

test_that("plot_peptide_overlap_map overlap helper counts residue support", {
  peps <- data.frame(
    protein_id = rep("synthetic", 4L),
    peptide = c("AB", "CD", "EF", "ABCD"),
    start = c(1L, 3L, 5L, 1L),
    end = c(2L, 4L, 6L, 4L),
    length = c(2L, 2L, 2L, 4L),
    missed_cleavages = c(0L, 0L, 0L, 1L),
    stringsAsFactors = FALSE
  )

  tile_df <- .build_peptide_overlap_df(
    peps,
    protein_length = 6L,
    length_range = NULL,
    residues_per_line = 4L
  )

  expect_equal(tile_df$residue, c("A", "B", "C", "D", "E", "F"))
  expect_equal(tile_df$overlap_count, c(2L, 2L, 2L, 2L, 1L, 1L))
  expect_equal(
    as.character(tile_df$overlap_class),
    c(
      rep("Detected twice", 4L),
      rep("Detected once", 2L)
    )
  )

  no_mc <- peps
  no_mc$missed_cleavages <- NULL
  no_mc_tiles <- .build_peptide_overlap_df(
    no_mc, protein_length = 6L, length_range = c(2L, 4L),
    residues_per_line = 4L
  )
  expect_equal(no_mc_tiles$overlap_count, c(2L, 2L, 2L, 2L, 1L, 1L))
})


# Restored section from the original plotting test surface.

# Restored section from the original plotting test surface.
test_that("pepvet_plot_config returns current config when called with no args", {
  cfg <- pepvet_plot_config()
  expect_type(cfg, "list")
  expect_named(cfg, c("palette", "params", "theme"))
  expect_type(cfg$palette, "list")
  expect_type(cfg$params, "list")
})

test_that("pepvet_plot_config validates palette names", {
  expect_error(
    pepvet_plot_config(palette = list(nonexistent = "#000000")),
    class = "pepvet_error_invalid_config"
  )
})

test_that("pepvet_plot_config validates params names", {
  expect_error(
    pepvet_plot_config(params = list(nonexistent = 99)),
    class = "pepvet_error_invalid_config"
  )
})

test_that("pepvet_plot_config validates palette and parameter values", {
  expect_error(
    pepvet_plot_config(palette = list(brand = 1)),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    pepvet_plot_config(palette = list(brand = "not-a-color")),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    pepvet_plot_config(params = list(verdict_good = 1.1)),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    pepvet_plot_config(params = list(length_lo = 0L)),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    pepvet_plot_config(
      params = list(verdict_good = 0.3, verdict_moderate = 0.4)
    ),
    class = "pepvet_error_invalid_config"
  )
})

test_that("pepvet_plot_config commits vector overrides atomically", {
  pepvet_plot_config(
    palette = list(verdict = c(
      Good = "#008800", Moderate = "#AA8800", Poor = "#880000"
    )),
    params = list(length_lo = 8L, length_hi = 24L)
  )
  on.exit(pepvet_plot_config_reset(), add = TRUE)
  cfg <- pepvet_plot_config()
  expect_identical(cfg$palette$verdict[["Good"]], "#008800")
  expect_identical(cfg$params$length_lo, 8L)
  expect_identical(cfg$params$length_hi, 24L)

  before <- cfg
  expect_error(
    pepvet_plot_config(
      palette = list(brand = "#004488"),
      params = list(length_lo = 30L)
    ),
    class = "pepvet_error_invalid_config"
  )
  expect_identical(pepvet_plot_config(), before)
})

test_that("pepvet_plot_config validates and applies theme overrides", {
  pepvet_plot_config(theme = list(legend.position = "right"))
  on.exit(pepvet_plot_config_reset(), add = TRUE)
  p <- plot_length_distribution(.fix_bsa_trypsin)
  expect_s3_class(p, "ggplot")
  expect_identical(p$theme$legend.position, "right")

  before <- pepvet_plot_config()
  expect_error(
    pepvet_plot_config(theme = list(not_a_theme_element = 1)),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    pepvet_plot_config(theme = 1),
    class = "pepvet_error_invalid_config"
  )
  expect_error(
    pepvet_plot_config(theme = list(complete = 1)),
    class = "pepvet_error_invalid_config"
  )
  expect_identical(pepvet_plot_config(), before)
})

test_that("pepvet_plot_config palette changes propagate to plots", {
  pepvet_plot_config(palette = list(valid = "#FF0000"))
  on.exit(pepvet_plot_config_reset(), add = TRUE)
  res <- .fix_bsa_trypsin
  p <- plot_length_distribution(res)
  expect_s3_class(p, "gg")
  expect_identical(p$scales$get_scales("fill")$palette(3L)[["Valid"]],
    "#FF0000")
})

test_that("plot palette overrides reach diagnostics plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  pepvet_plot_config(palette = list(good = "#008800"))
  on.exit(pepvet_plot_config_reset(), add = TRUE)
  diagnostics <- score_diagnostics(.fix_batch_small)
  p <- plot_score_diagnostics(diagnostics)
  expect_s3_class(p, "patchwork")
  vif_panel <- p$patches$plots[[1L]]
  fill_scale <- vif_panel$scales$get_scales("fill")
  expect_identical(unname(fill_scale$palette(4L)[[1L]]), "#008800")
})

test_that("pepvet_plot_config params changes propagate to plots", {
  pepvet_plot_config(params = list(verdict_good = 0.80))
  on.exit(pepvet_plot_config_reset())
  expect_equal(.pepvet_params$verdict_good, 0.80)
})

test_that("pepvet_plot_config_reset restores defaults", {
  pepvet_plot_config(params = list(verdict_good = 0.99))
  pepvet_plot_config(palette = list(valid = "#FF0000"))
  pepvet_plot_config_reset()
  expect_equal(.pepvet_params$verdict_good, 0.65)
  expect_identical(.pepvet_pal$valid, "#2C5F8A")
})

# pepvet_save_figure.

test_that("pepvet_save_figure saves single ggplot", {
  p <- plot_length_distribution(.fix_bsa_trypsin)
  f <- tempfile(fileext = ".png")
  result <- pepvet_save_figure(p, f)
  expect_identical(result, normalizePath(f, mustWork = FALSE))
  expect_true(file.exists(result))
  expect_identical(
    readBin(f, what = "raw", n = 8L),
    as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  )
  unlink(f)
})

test_that("pepvet_save_figure saves patchwork", {
  p <- plot_digest_profile(.fix_bsa_trypsin)
  f <- tempfile(fileext = ".png")
  result <- pepvet_save_figure(p, f)
  expect_identical(result, normalizePath(f, mustWork = FALSE))
  expect_true(file.exists(result))
  unlink(f)
})

test_that("pepvet_save_figure respects custom dimensions", {
  p <- plot_length_distribution(.fix_bsa_trypsin)
  f <- tempfile(fileext = ".pdf")
  suppressWarnings(
    pepvet_save_figure(p, f, width = 5, height = 4, device = "pdf")
  )
  pdf_lines <- suppressWarnings(readLines(f, warn = FALSE))
  expect_match(pdf_lines[[1L]], "%PDF")
  expect_true(any(suppressWarnings(grepl(
    "/MediaBox \\[0 0 360 288\\]", pdf_lines
  ))))
  unlink(f)
})

test_that("pepvet_save_figure errors on invalid plot object", {
  f <- tempfile(fileext = ".png")
  expect_error(
    pepvet_save_figure("not_a_plot", f),
    class = "pepvet_error_invalid_plot"
  )
})

# pepvet_theme_manuscript and pepvet_theme_presentation.

test_that("pepvet_theme_manuscript returns a theme object", {
  t <- pepvet_theme_manuscript()
  expect_s3_class(t, "theme")
  expect_equal(ggplot2::calc_element("plot.title", t)$size, 10)
})

test_that("pepvet_theme_presentation returns a theme object", {
  t <- pepvet_theme_presentation()
  expect_s3_class(t, "theme")
  expect_equal(ggplot2::calc_element("plot.title", t)$size, 15)
})

test_that("theme presets can be added to plots", {
  p <- plot_length_distribution(.fix_bsa_trypsin)
  p1 <- p + pepvet_theme_manuscript()
  p2 <- p + pepvet_theme_presentation()
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
  expect_equal(ggplot2::calc_element("plot.title", p1$theme)$size, 10)
  expect_equal(ggplot2::calc_element("plot.title", p2$theme)$size, 15)
})
