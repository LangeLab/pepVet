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
  expect_equal(
    d$weights,
    pepVet:::.default_scoring_weights$protein_only,
    tolerance = 1e-12
  )
})

test_that("VIF values are finite and positive", {
  d <- score_diagnostics(.fix_batch_small)

  expect_true(all(is.finite(d$vif)))
  expect_true(all(d$vif > 0))
})

test_that("orthogonal component fixtures have unit VIF", {
  batch <- data.frame(
    S_length = c(0, 0, 0, 1, 1, 1),
    S_coverage = c(0, 1, 0, 0, 1, 0)
  )

  result <- score_diagnostics(batch)

  expect_equal(result$vif, c(S_length = 1, S_coverage = 1), tolerance = 1e-12)
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

test_that("one-row diagnostic input fails before PCA or regression", {
  one_row <- data.frame(
    S_length = 0.2,
    S_coverage = 0.4
  )

  expect_error(
    score_diagnostics(one_row),
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

test_that("ablation values match an independent weighted-drop oracle", {
  # Keep composites away from verdict thresholds. Boundary classification is
  # tested separately; this fixture isolates the ablation arithmetic contract.
  batch <- data.frame(
    S_length = c(0.2, 0.76, 0.5),
    S_coverage = c(0.4, 0.6, 0.2),
    stringsAsFactors = FALSE
  )
  weights <- c(S_length = 0.25, S_coverage = 0.75)
  true_composite <- weights[["S_length"]] * batch$S_length +
    weights[["S_coverage"]] * batch$S_coverage

  expect_warning(
    result <- score_diagnostics(batch, weights = weights),
    class = "pepvet_warning_diagnostics_vif"
  )

  expected_length_drop <- weights[["S_length"]] * batch$S_length
  expected_coverage_drop <- weights[["S_coverage"]] * batch$S_coverage
  expected_means <- c(
    mean(expected_length_drop),
    mean(expected_coverage_drop)
  )
  expected_sds <- c(
    stats::sd(expected_length_drop),
    stats::sd(expected_coverage_drop)
  )
  expected_maxima <- c(
    max(expected_length_drop),
    max(expected_coverage_drop)
  )
  true_verdict <- ifelse(
    true_composite >= 0.65, "Good",
    ifelse(true_composite >= 0.40, "Moderate", "Poor")
  )
  perturbed_length <- true_composite - expected_length_drop
  perturbed_coverage <- true_composite - expected_coverage_drop
  expected_flips <- c(
    sum(ifelse(perturbed_length >= 0.65, "Good",
      ifelse(perturbed_length >= 0.40, "Moderate", "Poor")) != true_verdict),
    sum(ifelse(perturbed_coverage >= 0.65, "Good",
      ifelse(perturbed_coverage >= 0.40, "Moderate", "Poor")) != true_verdict)
  )

  expect_equal(result$ablation$mean_drop, expected_means, tolerance = 1e-12)
  expect_equal(result$ablation$sd_drop, expected_sds, tolerance = 1e-12)
  expect_equal(result$ablation$max_drop, expected_maxima, tolerance = 1e-12)
  expect_identical(result$ablation$n_verdict_flipped, expected_flips)
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
  for (invalid in list(
    NA_real_,
    numeric(0),
    c(S_length = 0.2, S_coverage = 0.2, S_count = 0.2,
      S_hydro = 0.2, S_charge = Inf),
    c(S_length = 0.2, S_coverage = 0.2, S_count = 0.2,
      S_hydro = 0.2, S_charge = 0.2, S_other = 0),
    "not numeric",
    c(S_length = 0.2, S_length = 0.3, S_coverage = 0.5,
      S_count = 0, S_hydro = 0, S_charge = 0),
    c(S_length = 0, S_coverage = 0, S_count = 0,
      S_hydro = 0, S_charge = 0)
  )) {
    expect_error(
      score_diagnostics(.fix_batch_small, weights = invalid),
      class = "pepvet_error_invalid_weights"
    )
  }
})

test_that("diagnostics rejects empty, nonnumeric, missing, non-finite, and constant matrices", {
  base <- data.frame(
    S_length = c(0.2, 0.4, 0.6),
    S_coverage = c(0.3, 0.5, 0.7),
    stringsAsFactors = FALSE
  )
  invalid <- list(
    empty = base[0, , drop = FALSE],
    nonnumeric = data.frame(
      S_length = c("0.2", "0.4", "0.6"),
      S_coverage = base$S_coverage
    ),
    missing = within(base, S_length <- c(0.2, NA_real_, 0.6)),
    infinite = within(base, S_coverage <- c(0.3, Inf, 0.7)),
    outside = within(base, S_length <- c(0.2, 1.2, 0.6)),
    constant = within(base, S_length <- rep(0.2, 3L))
  )
  invalid$unknown_component <- within(base, S_unknown <- S_length)
  invalid$bad_composite <- within(
    base, composite_score <- c(0.2, Inf, 0.6)
  )
  invalid$bad_verdict <- within(
    base, verdict <- c("Good", "invalid", "Poor")
  )

  for (case in invalid) {
    expect_error(
      score_diagnostics(case),
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  duplicate <- data.frame(
    S_length = c(0.2, 0.4, 0.6),
    S_coverage = c(0.3, 0.5, 0.7),
    S_length = c(0.3, 0.5, 0.7),
    check.names = FALSE
  )
  expect_error(
    score_diagnostics(duplicate),
    class = "pepvet_error_invalid_diagnostics_input"
  )

  duplicate_columns <- as.data.frame(base)
  names(duplicate_columns)[[2L]] <- names(duplicate_columns)[[1L]]
  expect_error(
    score_diagnostics(duplicate_columns),
    class = "pepvet_error_invalid_diagnostics_input"
  )
})

test_that("perfect collinearity is reported as infinite VIF without a warning", {
  batch <- data.frame(
    S_length = seq(0.1, 0.6, by = 0.1),
    S_coverage = seq(0.1, 0.6, by = 0.1)
  )

  expect_silent(result <- score_diagnostics(batch))
  expect_true(all(is.infinite(result$vif)))
  expect_equal(result$pca$var_explained[[1L]], 1, tolerance = 1e-12)
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
  expect_equal(
    d$weights,
    pepVet:::.default_scoring_weights$proteome_aware,
    tolerance = 1e-12
  )
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

test_that("plot_score_diagnostics rejects malformed input before plotting", {
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
  expect_identical(p$patches$layout$ncol, 3)
  vif_plot <- p$patches$plots[[1L]]
  built <- ggplot2::ggplot_build(vif_plot)
  fills <- unlist(lapply(built$data, function(data) {
    if ("fill" %in% names(data)) data$fill else character()
  }))
  expect_true(all(fills %in% unname(c(
    pepVet:::.pepvet_pal$good,
    pepVet:::.pepvet_pal$moderate,
    pepVet:::.pepvet_pal$poor,
    pepVet:::.pepvet_pal$na_gray
  ))))
})

test_that("plot_score_diagnostics validates the complete result structure", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  d <- score_diagnostics(.fix_batch_small)
  invalid_results <- list(
    missing_vif = within(d, vif <- NULL),
    wrong_vif_length = within(d, vif <- vif[-1L]),
    invalid_vif = within(d, vif[[1L]] <- NaN),
    unknown_vif = within(d, names(vif)[[1L]] <- "S_unknown"),
    invalid_weights = within(d, weights[[1L]] <- -0.1),
    bad_pca_structure = within(d, pca <- list()),
    large_dimensions = within(
      d, n_proteins <- .Machine$integer.max + 1
    ),
    bad_loadings_type = within(
      d, pca$loadings <- matrix("bad", nrow = nrow(pca$loadings),
        ncol = ncol(pca$loadings))
    ),
    bad_scores_type = within(
      d, pca$x <- matrix("bad", nrow = nrow(pca$x), ncol = ncol(pca$x))
    ),
    negative_sdev = within(d, pca$sdev[[1L]] <- -1),
    bad_pca = within(d, pca$var_explained[[1L]] <- NA_real_),
    bad_ablation = within(d, ablation$mean_drop[[1L]] <- NA_real_),
    bad_dimensions = within(d, n_components <- 1L)
  )

  for (invalid in invalid_results) {
    expect_error(
      plot_score_diagnostics(invalid),
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  expect_error(
    plot_score_diagnostics(d, title = integer(0)),
    class = "pepvet_error_invalid_diagnostics_input"
  )
})

test_that("plot_score_diagnostics labels NA and infinite VIF values", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  expect_warning(
    too_few <- score_diagnostics(.fix_batch_small[1:3, ]),
    class = "pepvet_warning_diagnostics_vif"
  )
  na_plot <- plot_score_diagnostics(too_few)
  na_data <- ggplot2::ggplot_build(na_plot$patches$plots[[1L]])$data
  na_labels <- unlist(lapply(na_data, function(data) {
    if ("label" %in% names(data)) data$label else character()
  }))
  expect_true("N/A" %in% na_labels)

  perfect <- data.frame(
    S_length = seq(0.1, 0.6, by = 0.1),
    S_coverage = seq(0.1, 0.6, by = 0.1)
  )
  infinite_plot <- plot_score_diagnostics(score_diagnostics(perfect))
  infinite_data <- ggplot2::ggplot_build(
    infinite_plot$patches$plots[[1L]]
  )$data
  infinite_labels <- unlist(lapply(infinite_data, function(data) {
    if ("label" %in% names(data)) data$label else character()
  }))
  expect_true("Inf" %in% infinite_labels)
})

test_that("plot_score_diagnostics supports six-component results and titles", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  prot_digest <- digest_protein(.small_path, enzyme = "trypsin")
  batch <- batch_evaluate(.small_path, enzyme = "trypsin",
    proteome = prot_digest
  )
  result <- plot_score_diagnostics(
    score_diagnostics(batch),
    title = "Proteome diagnostic fixture"
  )

  expect_s3_class(result, "patchwork")
  expect_identical(result$patches$layout$ncol, 3)
})

# ---- Determinism ----

test_that("score_diagnostics is deterministic", {
  d1 <- score_diagnostics(.fix_batch_small)
  d2 <- score_diagnostics(.fix_batch_small)

  expect_equal(d1$vif, d2$vif)
  expect_equal(d1$ablation$mean_drop, d2$ablation$mean_drop)
})
