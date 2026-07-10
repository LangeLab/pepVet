# Tests aligned with the corresponding plotting source domain.
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

test_that("pepvet_plot_config palette changes propagate to plots", {
  pepvet_plot_config(palette = list(brand = "#FF0000"))
  on.exit(pepvet_plot_config_reset())
  res <- .fix_bsa_trypsin
  p <- plot_length_distribution(res)
  expect_s3_class(p, "gg")
})

test_that("pepvet_plot_config params changes propagate to plots", {
  pepvet_plot_config(params = list(verdict_good = 0.80))
  on.exit(pepvet_plot_config_reset())
  expect_equal(.pepvet_params$verdict_good, 0.80)
})

test_that("pepvet_plot_config_reset restores defaults", {
  pepvet_plot_config(params = list(verdict_good = 0.99))
  pepvet_plot_config_reset()
  expect_equal(.pepvet_params$verdict_good, 0.65)
})

# pepvet_save_figure.

test_that("pepvet_save_figure saves single ggplot", {
  p <- plot_length_distribution(.fix_bsa_trypsin)
  f <- tempfile(fileext = ".png")
  result <- pepvet_save_figure(p, f)
  expect_true(file.exists(f))
  expect_true(file.size(f) > 1000)
  unlink(f)
})

test_that("pepvet_save_figure saves patchwork", {
  p <- plot_digest_profile(.fix_bsa_trypsin)
  f <- tempfile(fileext = ".png")
  result <- pepvet_save_figure(p, f)
  expect_true(file.exists(f))
  unlink(f)
})

test_that("pepvet_save_figure respects custom dimensions", {
  p <- plot_length_distribution(.fix_bsa_trypsin)
  f <- tempfile(fileext = ".png")
  pepvet_save_figure(p, f, width = 5, height = 4, dpi = 72)
  expect_true(file.exists(f))
  unlink(f)
})

test_that("pepvet_save_figure errors on invalid plot object", {
  f <- tempfile(fileext = ".png")
  expect_error(pepvet_save_figure("not_a_plot", f))
})

# pepvet_theme_manuscript and pepvet_theme_presentation.

test_that("pepvet_theme_manuscript returns a theme object", {
  t <- pepvet_theme_manuscript()
  expect_s3_class(t, "theme")
})

test_that("pepvet_theme_presentation returns a theme object", {
  t <- pepvet_theme_presentation()
  expect_s3_class(t, "theme")
})

test_that("theme presets can be added to plots", {
  p <- plot_length_distribution(.fix_bsa_trypsin)
  p1 <- p + pepvet_theme_manuscript()
  p2 <- p + pepvet_theme_presentation()
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
})
