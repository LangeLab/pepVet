# Tests aligned with the corresponding plotting source domain.
test_that("plot_length_distribution returns a ggplot from evaluate_digest result", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_length_distribution(res)
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  expect_identical(levels(p$data$length_class),
    c("Valid", "Too short", "Too long"))
  expect_identical(p$labels$x, "Peptide length (aa)")
})

test_that("plot_length_distribution accepts a bare peptide data.frame", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  peps <- res$peptides
  p <- plot_length_distribution(peps)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(peps))
})

test_that("plot_length_distribution: show_density = FALSE renders without error", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_length_distribution(res, show_density = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false(any(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomLine")
  }, logical(1L))))
})

test_that("plot_length_distribution: custom length_range is respected", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_length_distribution(
    res$peptides, length_range = c(5L, 30L)
  )
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(layer) {
    is.data.frame(layer$data) && "xmin" %in% names(layer$data) &&
      any(layer$data$xmin == 4.5 & layer$data$xmax == 30.5)
  }, logical(1L))))
})

test_that("plot_length_distribution: custom title is accepted", {
  skip_if_not_installed("ggplot2")

  res <- .bsa_result()
  p <- plot_length_distribution(res, title = "My length plot")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "My length plot")
})

test_that("distribution plots use resolved result ranges", {
  metadata_result <- .fix_bsa_trypsin
  metadata_result$params$length_range <- c(1L, 100L)
  metadata_result$params$gravy_range <- c(-5, 5)

  length_plot <- plot_length_distribution(metadata_result)
  expect_match(length_plot$labels$subtitle, "100")

  gravy_plot <- plot_gravy_landscape(metadata_result)
  gravy_data <- gravy_plot$patches$plots[[2L]]$patches$plots[[1L]]$data
  expect_true(all(gravy_data$valid_length))
  expect_true(all(gravy_data$valid_gravy))

  pI_default <- plot_pI_distribution(.fix_bsa_trypsin)
  pI_metadata <- plot_pI_distribution(metadata_result)
  expect_gt(nrow(pI_metadata$data), nrow(pI_default$data))

  mz_default <- plot_mz_distribution(.fix_bsa_trypsin)
  mz_metadata <- plot_mz_distribution(metadata_result)
  expect_gt(nrow(mz_metadata$data), nrow(mz_default$data))
})

test_that("plot_length_distribution: H3.1 (short protein) renders without warnings", {
  skip_if_not_installed("ggplot2")

  res <- .h3_result()
  p <- expect_no_warning(plot_length_distribution(res))
  expect_s3_class(p, "ggplot")
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

test_that("plot_length_distribution rejects empty data and invalid ranges", {
  skip_if_not_installed("ggplot2")

  expect_error(
    plot_length_distribution(data.frame(length = integer(0))),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_length_distribution(.bsa_result(), length_range = c(25L, 7L)),
    class = "pepvet_error_invalid_length_range"
  )
  expect_no_warning(
    expect_error(
      plot_length_distribution(data.frame(length = 1e20)),
      class = "pepvet_error_invalid_digest_result"
    )
  )
  expect_error(
    plot_length_distribution(.bsa_result(), show_density = "yes"),
    class = "pepvet_error_invalid_input"
  )
  p <- plot_length_distribution(.bsa_result(), show_density = NULL)
  expect_s3_class(p, "ggplot")
})


# plot_gravy_landscape.

test_that("plot_gravy_landscape returns a patchwork/ggplot from evaluate_digest result", {
  res <- .fix_bsa_trypsin
  p <- plot_gravy_landscape(res)
  expect_s3_class(p, "patchwork")
  scatter <- p$patches$plots[[2L]]$patches$plots[[1L]]
  expect_true(all(c("length", "gravy", "class") %in% names(scatter$data)))
})

test_that("plot_gravy_landscape accepts a bare data.frame with length + gravy", {
  peps <- data.frame(
    length = c(8L, 12L, 30L, 5L),
    gravy  = c(0.2, -0.5, 0.8, -1.2)
  )
  p <- plot_gravy_landscape(peps)
  expect_s3_class(p, "patchwork")
  expect_true(any(p$patches$plots[[2L]]$patches$plots[[1L]]$data$class ==
    "Outside both"))
})

test_that("plot_gravy_landscape auto-computes GRAVY from peptide column", {
  peps <- data.frame(
    length  = c(8L, 12L),
    peptide = c("ACDEFGHIK", "LMNPQRSTVW")
  )
  p <- plot_gravy_landscape(peps)
  expect_s3_class(p, "patchwork")
  expect_true("gravy" %in% names(p$patches$plots[[2L]]$patches$plots[[1L]]$data))
})

test_that("plot_gravy_landscape labels a bounded set of outliers", {
  peps <- data.frame(
    length = c(10L, 12L, 14L),
    gravy = c(0.1, 2.0, -0.5),
    peptide = c("AAAAAAAAAA", "CCCCCCCCCC", "DDDDDDDDDD")
  )
  p <- plot_gravy_landscape(peps, label_outliers_n = 1L)
  expect_s3_class(p, "patchwork")
  scatter <- p$patches$plots[[2L]]$patches$plots[[1L]]
  expect_true(any(vapply(scatter$layers, function(layer) {
    inherits(layer$geom, "GeomText") && is.data.frame(layer$data)
  }, logical(1L))))
})

test_that("plot_gravy_landscape: custom length_range and gravy_range accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_gravy_landscape(
    res, length_range = c(5L, 30L), gravy_range = c(-1.5, 1.0)
  )
  expect_s3_class(p, "patchwork")
  expect_match(p$patches$annotation$title, "GRAVY landscape")
})

test_that("plot_gravy_landscape: label_outliers_n = 0L suppresses labels silently", {
  res <- .fix_bsa_trypsin
  p <- plot_gravy_landscape(res, label_outliers_n = 0L)
  expect_s3_class(p, "patchwork")
  scatter <- p$patches$plots[[2L]]$patches$plots[[1L]]
  expect_false(any(vapply(scatter$layers, function(layer) {
    inherits(layer$geom, "GeomText")
  }, logical(1L))))
})

test_that("plot_gravy_landscape handles constant GRAVY values", {
  peps <- data.frame(
    length = c(8L, 12L, 20L),
    gravy = c(0.2, 0.2, 0.2)
  )
  p <- expect_no_warning(plot_gravy_landscape(peps))
  expect_s3_class(p, "patchwork")
  expect_true(all(p$patches$plots[[2L]]$patches$plots[[1L]]$data$valid_gravy))
})

test_that("plot_gravy_landscape: custom title is accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_gravy_landscape(res, title = "My GRAVY plot")
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "My GRAVY plot")
})

test_that("plot_gravy_landscape: H3.1 (short protein) renders without warnings", {
  res <- .fix_h3_trypsin
  p <- expect_no_warning(plot_gravy_landscape(res))
  expect_s3_class(p, "patchwork")
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

test_that("plot_gravy_landscape rejects empty data and invalid ranges", {
  expect_error(
    plot_gravy_landscape(data.frame(length = integer(0), gravy = numeric(0))),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_gravy_landscape(.fix_bsa_trypsin, gravy_range = c(1, -1)),
    class = "pepvet_error_invalid_gravy_range"
  )
  expect_error(
    plot_gravy_landscape(.fix_bsa_trypsin, label_outliers_n = -1L),
    class = "pepvet_error_invalid_input"
  )
  expect_no_warning(
    expect_error(
      plot_gravy_landscape(.fix_bsa_trypsin, label_outliers_n = 1e20),
      class = "pepvet_error_invalid_input"
    )
  )
  expect_error(
    plot_gravy_landscape(data.frame(length = 10L)),
    class = "pepvet_error_invalid_digest_result"
  )
})


# plot_pI_distribution.

test_that("plot_pI_distribution returns a ggplot from evaluate_digest result", {
  res <- .fix_bsa_trypsin
  p <- plot_pI_distribution(res)
  expect_s3_class(p, "gg")
  expect_true(all(c("pI", "bin") %in% names(p$data)))
  expect_identical(p$labels$x, "Isoelectric point (pI)")
})

test_that("plot_pI_distribution accepts score_peptides(include_pI=TRUE) tibble", {
  res <- .fix_bsa_trypsin
  scored <- score_peptides(res$peptides, enzyme = "trypsin", include_pI = TRUE)
  p <- plot_pI_distribution(scored)
  expect_s3_class(p, "ggplot")
  expect_gt(nrow(p$data), 0L)
})

test_that("plot_pI_distribution accepts data.frame with numeric pI column", {
  df <- data.frame(pI = c(4.2, 5.7, 6.1, 8.3, 9.5))
  p <- plot_pI_distribution(df)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(df))
})

test_that("plot_pI_distribution accepts a bare numeric vector", {
  p <- plot_pI_distribution(c(4.2, 5.7, 6.1, 8.3, 9.5))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 5L)
})

test_that("plot_pI_distribution: show_fraction_lines = FALSE renders without error", {
  res <- .fix_bsa_trypsin
  p <- plot_pI_distribution(res, show_fraction_lines = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false(any(vapply(p$layers, function(layer) {
    inherits(layer$geom, "GeomVline")
  }, logical(1L))))
})

test_that("plot_pI_distribution: custom fraction_breaks accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_pI_distribution(res, fraction_breaks = c(4, 6, 8, 10))
  expect_s3_class(p, "ggplot")
  expect_equal(p$scales$get_scales("fill")$name, "SCX fraction")
})

test_that("plot_pI_distribution: custom title is accepted", {
  res <- .fix_bsa_trypsin
  p <- plot_pI_distribution(res, title = "My pI plot")
  expect_s3_class(p, "gg")
  expect_identical(p$labels$title, "My pI plot")
})

test_that("plot_pI_distribution: H3.1 renders without warnings", {
  res <- .fix_h3_trypsin
  p <- expect_no_warning(plot_pI_distribution(res))
  expect_s3_class(p, "ggplot")
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

test_that("plot_pI_distribution rejects empty, non-finite, and invalid breaks", {
  expect_error(
    plot_pI_distribution(numeric(0)),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_pI_distribution(c(4, Inf)),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_pI_distribution(c(4, 15)),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_pI_distribution(c(4, 5), fraction_breaks = c(4, 4)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    plot_pI_distribution(c(4, 5), show_fraction_lines = "yes"),
    class = "pepvet_error_invalid_input"
  )
})

test_that("plot_pI_distribution validates evaluate and list-column inputs", {
  no_valid <- NULL
  expect_warning(
    no_valid <- evaluate_digest("AAAA", enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  expect_error(
    plot_pI_distribution(no_valid),
    class = "pepvet_error_invalid_digest_result"
  )
  bad_list <- data.frame(pI = I(list("not numeric")))
  expect_error(
    plot_pI_distribution(bad_list),
    class = "pepvet_error_invalid_digest_result"
  )
  bad_vector <- data.frame(pI = "4.2")
  expect_error(
    plot_pI_distribution(bad_vector),
    class = "pepvet_error_invalid_digest_result"
  )
})


# multi-input: plot_length_distribution.

test_that("plot_length_distribution: named list produces faceted ggplot", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_length_distribution(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
  expect_true(".label" %in% names(p$data))
  expect_equal(nlevels(p$data$.label), 2L)
})

test_that("plot distribution auto-labels unnamed results and rejects malformed labels", {
  p <- plot_length_distribution(list(.fix_bsa_trypsin, .fix_h3_trypsin))
  expect_s3_class(p, "ggplot")
  expect_true(all(grepl(" / trypsin$", levels(p$data$.label))))
  expect_match(.result_label(.fix_bsa_trypsin), "P02769")
  expect_error(
    .result_label(list(params = list(protein_ids = "P"))),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_length_distribution: multi-input renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- expect_no_warning(
    plot_length_distribution(list(BSA = bsa_res, H3 = h3_res))
  )
  expect_s3_class(p, "ggplot")
})

test_that("multi-input distribution plots reject malformed leaves and labels", {
  multi_protein <- .fix_bsa_trypsin
  extra_peptides <- multi_protein$peptides
  extra_peptides$protein_id <- "second"
  multi_protein$peptides <- rbind(multi_protein$peptides, extra_peptides)
  extra_score <- multi_protein$scores
  extra_score$protein_id <- "second"
  multi_protein$scores <- rbind(multi_protein$scores, extra_score)
  multi_protein$params$protein_ids <- c(
    multi_protein$params$protein_ids, "second"
  )
  multi_plot <- plot_length_distribution(list(combined = multi_protein))
  expect_s3_class(multi_plot, "ggplot")
  expect_equal(nrow(multi_plot$data), nrow(multi_protein$peptides))

  bad_scores <- .fix_bsa_trypsin
  bad_scores$scores$composite_score[[1L]] <- Inf
  expect_error(
    plot_length_distribution(list(BSA = bad_scores, H3 = .fix_h3_trypsin)),
    class = "pepvet_error_invalid_digest_result"
  )

  bad_gravy <- .fix_bsa_trypsin
  bad_gravy$peptides$gravy <- Inf
  expect_error(
    plot_gravy_landscape(list(BSA = bad_gravy, H3 = .fix_h3_trypsin)),
    class = "pepvet_error_invalid_digest_result"
  )

  expect_error(
    plot_length_distribution(list(A = .fix_bsa_trypsin, A = .fix_h3_trypsin)),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_length_distribution(list(.fix_bsa_trypsin, .fix_bsa_trypsin)),
    class = "pepvet_error_invalid_digest_result"
  )
})

# multi-input: plot_gravy_landscape.

test_that("plot_gravy_landscape: named list produces faceted ggplot", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_gravy_landscape(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
  expect_true(".label" %in% names(p$data))
  expect_equal(nlevels(p$data$.label), 2L)
})

test_that("plot_gravy_landscape: multi-input renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- expect_no_warning(
    plot_gravy_landscape(list(BSA = bsa_res, H3 = h3_res))
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_gravy_landscape auto-labels unnamed result collections", {
  p <- plot_gravy_landscape(list(.fix_bsa_trypsin, .fix_h3_trypsin))
  expect_s3_class(p, "ggplot")
  expect_true(all(grepl(" / trypsin$", levels(p$data$.label))))
})

# multi-input: plot_pI_distribution.

test_that("plot_pI_distribution: named list produces overlaid density ggplot", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- plot_pI_distribution(list(BSA = bsa_res, H3 = h3_res))
  expect_s3_class(p, "gg")
  expect_true(".label" %in% names(p$data))
})

test_that("plot_pI_distribution: multi-input renders without warnings", {
  bsa_res <- .fix_bsa_trypsin
  h3_res <- .fix_h3_trypsin
  p <- expect_no_warning(
    plot_pI_distribution(list(BSA = bsa_res, H3 = h3_res))
  )
  expect_s3_class(p, "ggplot")
})



# Restored section from the original plotting test surface.

# Restored section from the original plotting test surface.
test_that("plot_missed_cleavage_impact returns ggplot from MC list", {
  r0 <- .fix_bsa_mc0
  r1 <- .fix_bsa_mc1
  r2 <- .fix_bsa_mc2
  p <- plot_missed_cleavage_impact(list("MC=0" = r0, "MC=1" = r1, "MC=2" = r2))
  expect_s3_class(p, "gg")
  expect_true(all(c("x_idx", "score", "component") %in% names(p$data)))
  expect_true("Composite" %in% as.character(p$data$component))
})

test_that("plot_missed_cleavage_impact: unnamed list auto-named", {
  r0 <- .fix_bsa_mc0
  r1 <- .fix_bsa_mc1
  p <- plot_missed_cleavage_impact(list(r0, r1))
  expect_s3_class(p, "gg")
  expect_equal(levels(p$data$mc_label), c("MC=0", "MC=1"))
})

test_that("plot_missed_cleavage_impact renders without warnings", {
  r0 <- .fix_bsa_mc0
  r1 <- .fix_bsa_mc1
  r2 <- .fix_bsa_mc2
  p <- expect_no_warning(
    plot_missed_cleavage_impact(list("MC=0" = r0, "MC=1" = r1, "MC=2" = r2))
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_missed_cleavage_impact: chymotrypsin-high MC=0/1 renders", {
  # Chymotrypsin produces many more peptides than trypsin. Test the line plot
  # rendering when S_count and S_length change more dramatically across MC levels
  r0 <- .fix_bsa_chymotryp
  r1 <- .fix_bsa_chymotryp_mc1
  p <- plot_missed_cleavage_impact(list("MC=0" = r0, "MC=1" = r1))
  expect_s3_class(p, "gg")
  expect_match(p$labels$title, "chymotrypsin")
})

test_that("plot_missed_cleavage_impact errors on single-element list", {
  r0 <- .fix_bsa_mc0
  expect_error(
    plot_missed_cleavage_impact(list(r0)),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_missed_cleavage_impact(list(
      .fix_bsa_mc0, .fix_bsa_mc1, .fix_bsa_mc2,
      .fix_bsa_mc0, .fix_bsa_mc1
    )),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_missed_cleavage_impact validates selected components and levels", {
  expect_error(
    plot_missed_cleavage_impact(list(.fix_bsa_mc0, .fix_bsa_mc1), components = NULL),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_missed_cleavage_impact(
      list("MC=0" = .fix_bsa_mc0, "MC=1" = .fix_bsa_mc1),
      components = "missing_score"
    ),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    {
      bad_names <- list(.fix_bsa_mc0, .fix_bsa_mc1)
      names(bad_names) <- c("", "MC=1")
      plot_missed_cleavage_impact(bad_names)
    },
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_missed_cleavage_impact(list(.fix_bsa_mc0, list(scores = data.frame()))),
    class = "pepvet_error_invalid_digest_result"
  )
  mixed_enzyme <- .fix_bsa_mc1
  mixed_enzyme$params$enzyme <- "lysc"
  expect_error(
    plot_missed_cleavage_impact(
      list("MC=0" = .fix_bsa_mc0, "MC=1" = mixed_enzyme)
    ),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_missed_cleavage_impact(
      list("MC=0" = .fix_bsa_mc0, "MC=1" = .fix_h3_trypsin)
    ),
    class = "pepvet_error_invalid_digest_result"
  )
  unique_score <- .fix_bsa_mc0
  unique_score$scores$S_unique <- 0.5
  p <- plot_missed_cleavage_impact(
    list("MC=0" = unique_score, "MC=1" = unique_score),
    components = "S_unique"
  )
  expect_s3_class(p, "ggplot")
  expect_true("Uniqueness" %in% as.character(p$data$component))
})

test_that("plot_mz_distribution returns ggplot from evaluate_digest", {
  p <- plot_mz_distribution(.fix_bsa_trypsin)
  expect_s3_class(p, "gg")
  expect_true(all(c("mz", "charge_state") %in% names(p$data)))
  expect_equal(levels(p$data$charge_state), c("z = +2", "z = +3"))
})

test_that("plot_mz_distribution: show_rug = TRUE and FALSE both work", {
  p1 <- plot_mz_distribution(.fix_bsa_trypsin, show_rug = TRUE)
  p2 <- plot_mz_distribution(.fix_bsa_trypsin, show_rug = FALSE)
  expect_s3_class(p1, "gg")
  expect_s3_class(p2, "gg")
  expect_true(sum(vapply(p1$layers, function(layer) {
    inherits(layer$geom, "GeomRug")
  }, logical(1L))) > sum(vapply(p2$layers, function(layer) {
    inherits(layer$geom, "GeomRug")
  }, logical(1L))))
})

test_that("plot_mz_distribution: custom scan_range", {
  p <- plot_mz_distribution(.fix_bsa_trypsin, scan_range = c(400, 1000))
  expect_s3_class(p, "gg")
  expect_match(p$labels$subtitle, "400")
  expect_error(
    plot_mz_distribution(.fix_bsa_trypsin, scan_range = c(-1, 1000)),
    class = "pepvet_error_invalid_input"
  )
  expect_no_warning(
    expect_error(
      plot_mz_distribution(.fix_bsa_trypsin, charge_states = 1e20),
      class = "pepvet_error_invalid_input"
    )
  )
  precomputed <- data.frame(
    mz = c(400, -1),
    charge_state = c("z = +2", "z = +2")
  )
  expect_error(
    plot_mz_distribution(precomputed),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_mz_distribution: multi-input mode", {
  p <- plot_mz_distribution(
    list(BSA = .fix_bsa_trypsin, H3 = .fix_h3_trypsin)
  )
  expect_s3_class(p, "gg")
  expect_true(".label" %in% names(p$data))
  expect_equal(nlevels(p$data$.label), 2L)
})

test_that("plot_mz_distribution: multi-input palette covers extra charge states", {
  p <- plot_mz_distribution(
    list(BSA = .fix_bsa_trypsin, H3 = .fix_h3_trypsin),
    charge_states = 2:6
  )
  fill_scale <- p$scales$get_scales("fill")
  expect_equal(length(fill_scale$palette(5L)), 5L)
  expect_false(anyNA(fill_scale$palette(5L)))
  expect_equal(levels(p$data$charge_state), paste0("z = +", 2:6))
})

test_that("plot_mz_distribution: unnamed and empty multi-input branches are classed", {
  p <- plot_mz_distribution(list(.fix_bsa_trypsin, .fix_h3_trypsin))
  expect_s3_class(p, "ggplot")
  expect_true(all(grepl(" / trypsin$", levels(p$data$.label))))

  no_valid <- NULL
  expect_warning(
    no_valid <- evaluate_digest("AAAA", enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  p <- plot_mz_distribution(
    list(no_valid, .fix_bsa_trypsin), title = "m/z comparison"
  )
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "m/z comparison")
  expect_error(
    plot_mz_distribution(list(no_valid, no_valid)),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_pI_distribution ignores a result with no valid peptides in a comparison", {
  no_valid <- NULL
  expect_warning(
    no_valid <- evaluate_digest("AAAA", enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  p <- plot_pI_distribution(list(no_valid, .fix_bsa_trypsin))
  expect_s3_class(p, "ggplot")
  expect_equal(nlevels(p$data$.label), 2L)
  expect_error(
    plot_pI_distribution(list(no_valid, no_valid)),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_mz_distribution accepts precomputed values and validates edges", {
  precomputed <- data.frame(
    mz = c(400, 500, 700, 800),
    charge_state = c(2L, 2L, 3L, 3L)
  )
  p <- plot_mz_distribution(
    precomputed, charge_states = NULL, scan_range = c(450, 750)
  )
  expect_s3_class(p, "ggplot")
  expect_setequal(levels(p$data$charge_state), c("2", "3"))
  expect_error(
    plot_mz_distribution(precomputed, scan_range = c(750, 450)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    plot_mz_distribution(.fix_bsa_trypsin, charge_states = c(0L, 2L)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    plot_mz_distribution(.fix_bsa_trypsin, scan_range = c(NA_real_, 1)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    plot_mz_distribution(
      data.frame(peptide = character(0)),
      charge_states = 2L
    ),
    class = "pepvet_error_invalid_digest_result"
  )
  p <- plot_mz_distribution(
    data.frame(peptide = c("ACDEFGHIK", "LMNPQRSTV")),
    charge_states = c(2L, 3L, 4L),
    show_rug = NULL
  )
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$charge_state),
    c("z = +2", "z = +3", "z = +4"))
  expect_error(
    plot_mz_distribution("not a result"),
    class = "pepvet_error_invalid_digest_result"
  )
})

test_that("plot_mz_distribution validates precomputed labels and computed ranges", {
  precomputed <- data.frame(
    mz = c(400, 500),
    charge_state = c("z = +2", "z = +3"),
    stringsAsFactors = FALSE
  )
  p <- plot_mz_distribution(precomputed, title = "Precomputed m/z")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Precomputed m/z")

  bad_labels <- precomputed
  bad_labels$charge_state[[1L]] <- ""
  expect_error(
    plot_mz_distribution(bad_labels),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_mz_distribution(
      data.frame(mz = numeric(0), charge_state = character(0))
    ),
    class = "pepvet_error_invalid_digest_result"
  )
  bad_charge <- precomputed
  bad_charge$charge_state <- c(.Machine$integer.max + 1, 2)
  expect_error(
    plot_mz_distribution(bad_charge),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_mz_distribution(
      data.frame(length = c(4L, 5L), peptide = c("AAAA", "CCCCC")),
      length_range = c(7L, 25L)
    ),
    class = "pepvet_error_invalid_digest_result"
  )
  expect_error(
    plot_mz_distribution(data.frame(length = 8L), charge_states = 2L),
    class = "pepvet_error_invalid_digest_result"
  )
  no_id <- .fix_bsa_trypsin
  no_id$params$protein_ids <- NULL
  p <- plot_mz_distribution(no_id)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Precursor m/z distribution")
})

test_that("distribution plots reject non-scalar titles", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  mc_results <- list(mc0 = .fix_bsa_mc0, mc1 = .fix_bsa_mc1)
  runners <- list(
    length = function(title) {
      plot_length_distribution(.fix_bsa_trypsin, title = title)
    },
    gravy = function(title) {
      plot_gravy_landscape(.fix_bsa_trypsin, title = title)
    },
    pI = function(title) {
      plot_pI_distribution(.fix_bsa_trypsin, title = title)
    },
    missed = function(title) {
      plot_missed_cleavage_impact(mc_results, title = title)
    },
    mz = function(title) {
      plot_mz_distribution(.fix_bsa_trypsin, title = title)
    }
  )

  for (runner_name in names(runners)) {
    expect_error(
      runners[[runner_name]](c("one", "two")),
      class = "pepvet_error_invalid_input",
      info = runner_name
    )
  }
})
