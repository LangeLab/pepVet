# Tests for score_diagnostics() and plot_score_diagnostics()

# ---- Happy path ----

test_that("score_diagnostics returns correct structure", {
  d <- score_diagnostics(.fix_batch_small)

  expect_type(d, "list")
  expect_named(d, c("vif", "pca", "ablation", "n_proteins",
    "n_components", "weights"))
  expect_named(d$vif, c("S_length", "S_coverage", "S_count",
    "S_hydro", "S_charge"))
  expect_named(d$pca, c("var_explained", "loadings", "sdev", "x"))
  expect_true(is.data.frame(d$ablation))
  expect_named(d$ablation, c("component", "weight", "mean_drop",
    "sd_drop", "max_drop", "n_verdict_flipped"))
  expect_equal(d$n_proteins, 50L)
  expect_equal(d$n_components, 5L)
})

test_that("VIF values are finite and positive", {
  d <- score_diagnostics(.fix_batch_small)

  expect_true(all(is.finite(d$vif)))
  expect_true(all(d$vif > 0))
})

test_that("PCA variance sums to 1", {
  d <- score_diagnostics(.fix_batch_small)

  expect_equal(sum(d$pca$var_explained), 1, tolerance = 1e-10)
})

test_that("PCA loadings are a square matrix", {
  d <- score_diagnostics(.fix_batch_small)

  expect_equal(dim(d$pca$loadings), c(5, 5))
})

test_that("Ablation drops are non-negative", {
  d <- score_diagnostics(.fix_batch_small)

  expect_true(all(d$ablation$mean_drop >= 0))
  expect_true(all(d$ablation$max_drop >= 0))
})

test_that("Ablation rows match component count", {
  d <- score_diagnostics(.fix_batch_small)

  expect_equal(nrow(d$ablation), 5L)
})

# ---- Edge cases ----

test_that("too-few-proteins produces NA VIF with warning", {
  batch3 <- .fix_batch_small[1:3, ]

  expect_warning(
    d <- score_diagnostics(batch3),
    class = "pepvet_warning_diagnostics_vif"
  )
  expect_true(all(is.na(d$vif)))
  expect_true(all(d$ablation$mean_drop >= 0))
})

test_that("single component column raises error", {
  bad <- .fix_batch_small[, c("protein_id", "S_length"), drop = FALSE]

  expect_error(
    score_diagnostics(bad),
    class = "pepvet_error_invalid_diagnostics_input"
  )
})

test_that("non-data-frame input raises error", {
  expect_error(
    score_diagnostics(NULL),
    class = "pepvet_error_invalid_diagnostics_input"
  )
  expect_error(
    score_diagnostics("not_a_tibble"),
    class = "pepvet_error_invalid_diagnostics_input"
  )
})

# ---- Custom weights ----

test_that("custom weights produce proportional ablation drops", {
  custom_weights <- c(
    S_length = 0.3, S_coverage = 0.3, S_count = 0.2,
    S_hydro = 0.1, S_charge = 0.1)
  d <- score_diagnostics(.fix_batch_small, weights = custom_weights)

  expect_equal(d$weights, custom_weights / sum(custom_weights),
    tolerance = 1e-10)
  expect_true(all(d$ablation$mean_drop >= 0))
})

test_that("score_diagnostics rejects invalid custom weights", {
  expect_error(
    score_diagnostics(.fix_batch_small, weights = rep(0.2, 5L)),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    score_diagnostics(
      .fix_batch_small,
      weights = c(
        S_length = 0.3, S_coverage = 0.3, S_count = 0.2,
        S_hydro = -0.1, S_charge = 0.3
      )
    ),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    score_diagnostics(
      .fix_batch_small,
      weights = c(
        S_length = 0.2, S_coverage = 0.2, S_count = 0.2,
        S_hydro = 0.2, S_other = 0.2
      )
    ),
    class = "pepvet_error_invalid_weights"
  )
})

# ---- Proteome-aware mode ----

test_that("proteome-aware mode detects S_unique", {
  prot_digest <- digest_protein(.small_path, enzyme = "trypsin")
  batch_pa <- batch_evaluate(.small_path, enzyme = "trypsin",
    proteome = prot_digest)
  d <- score_diagnostics(batch_pa)

  expect_equal(d$n_components, 6L)
  expect_true("S_unique" %in% names(d$vif))
  expect_true("S_unique" %in% d$ablation$component)
})

# ---- Different enzyme ----

test_that("different enzyme produces valid diagnostics", {
  expect_warning(
    batch_lysc <- batch_evaluate(.small_path, enzyme = "lysc"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  d <- score_diagnostics(batch_lysc)

  expect_true(all(is.finite(d$vif)))
  expect_true(all(d$ablation$mean_drop >= 0))
})

# ---- Plot function ----

test_that("plot_score_diagnostics returns a patchwork", {
  d <- score_diagnostics(.fix_batch_small)

  expect_error(
    plot_score_diagnostics(NULL),
    class = "pepvet_error_invalid_diagnostics_input"
  )
  expect_error(
    plot_score_diagnostics(list(a = 1)),
    class = "pepvet_error_invalid_diagnostics_input"
  )
})

test_that("plot_score_diagnostics works with ggplot2 and patchwork", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  d <- score_diagnostics(.fix_batch_small)
  p <- plot_score_diagnostics(d)

  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "patchwork")
})

# ---- Determinism ----

test_that("score_diagnostics is deterministic", {
  d1 <- score_diagnostics(.fix_batch_small)
  d2 <- score_diagnostics(.fix_batch_small)

  expect_equal(d1$vif, d2$vif)
  expect_equal(d1$ablation$mean_drop, d2$ablation$mean_drop)
})
