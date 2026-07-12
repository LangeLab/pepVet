# Tests aligned with the corresponding plotting source domain.
grob_text_labels <- function(grob) {
  labels <- character(0L)
  if (is.character(grob$label)) {
    labels <- c(labels, grob$label)
  }
  if (!is.null(grob$grobs)) {
    labels <- c(labels, unlist(lapply(grob$grobs, grob_text_labels)))
  }
  if (!is.null(grob$children)) {
    labels <- c(labels, unlist(lapply(grob$children, grob_text_labels)))
  }
  labels
}

test_that("plot_enzyme_comparison returns a patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp)
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "P02769  (ALBU_BOVIN)  ·  Enzyme comparison")
  expect_true(any(c("Coverage", "Length") %in%
    as.character(p$patches$plots[[1L]]$data$score_name)))
  expect_true("★ Recommended" %in%
    grob_text_labels(patchwork::patchworkGrob(p)))
})

test_that("plot_enzyme_comparison: recommend = FALSE suppresses badge", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp, recommend = FALSE)
  expect_s3_class(p, "patchwork")
  expect_equal(as.numeric(p$patches$layout$widths), c(2.8, 1))
  expect_false("★ Recommended" %in%
    grob_text_labels(patchwork::patchworkGrob(p)))
})

test_that("plot_enzyme_comparison: subset of scores renders without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp, scores = c("S_coverage", "S_count"))
  expect_s3_class(p, "patchwork")
  expect_setequal(
    unique(as.character(p$patches$plots[[1L]]$data$score_name)),
    c("Coverage", "Peptide count")
  )
  comp$S_unique <- seq(0.2, 0.4, length.out = nrow(comp))
  uniqueness <- plot_enzyme_comparison(comp, scores = "S_unique")
  expect_s3_class(uniqueness, "patchwork")
  expect_identical(
    unique(as.character(uniqueness$patches$plots[[1L]]$data$score_name)),
    "Uniqueness"
  )
})

test_that("plot_enzyme_comparison supports score-free fallback titles and rejects unknown scores", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison(enzymes = c("trypsin", "lysc"))
  comp$protein_id <- NULL
  p <- plot_enzyme_comparison(comp, scores = "S_coverage")
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "Enzyme comparison")
  expect_error(
    plot_enzyme_comparison(comp, scores = "not_a_score"),
    class = "pepvet_error_invalid_comparison"
  )
})

test_that("plot_enzyme_comparison: custom title is accepted", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison()
  p <- plot_enzyme_comparison(comp, title = "My comparison")
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "My comparison")
})

test_that("plot_enzyme_comparison: two-enzyme comparison works", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison(enzymes = c("trypsin", "lysc"))
  p <- plot_enzyme_comparison(comp)
  expect_s3_class(p, "patchwork")
  expect_equal(nlevels(p$patches$plots[[1L]]$data$enzyme), 2L)
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

test_that("plot_enzyme_comparison rejects invalid score and recommendation inputs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  comp <- .bsa_comparison(enzymes = c("trypsin", "lysc"))
  expect_error(
    plot_enzyme_comparison(comp, scores = NULL),
    class = "pepvet_error_invalid_comparison"
  )
  expect_error(
    plot_enzyme_comparison(comp, recommend = NULL),
    class = "pepvet_error_invalid_comparison"
  )
  bad <- comp
  bad$composite_score[[1L]] <- Inf
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
  bad <- comp
  bad$S_coverage[[1L]] <- -0.1
  expect_error(
    plot_enzyme_comparison(bad),
    class = "pepvet_error_invalid_comparison"
  )
  bad <- comp
  bad$protein_id[[1L]] <- "different-protein"
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
