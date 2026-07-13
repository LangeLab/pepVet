# Tests aligned with the corresponding plotting source domain.
batch_grob_text_labels <- function(grob) {
  labels <- if (is.character(grob$label)) grob$label else character(0L)
  if (!is.null(grob$grobs)) {
    labels <- c(labels, unlist(lapply(grob$grobs, batch_grob_text_labels)))
  }
  if (!is.null(grob$children)) {
    labels <- c(labels, unlist(lapply(grob$children, batch_grob_text_labels)))
  }
  labels
}

test_that("plot_proteome_overview returns patchwork from batch_evaluate", {
  batch <- batch_evaluate(
    Biostrings::readAAStringSet(.bsa_path),
    enzyme = "trypsin"
  )
  p <- plot_proteome_overview(batch)
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title,
    "Proteome Digest Overview  ·  1 proteins")
  expect_identical(as.numeric(p$patches$layout$heights), c(1.35, 1))
})

test_that("plot_proteome_overview errors on empty batch", {
  empty <- data.frame(
    protein_id = character(0), composite_score = numeric(0),
    verdict = character(0), stringsAsFactors = FALSE
  )
  expect_error(
    plot_proteome_overview(empty),
    class = "pepvet_error_invalid_batch"
  )
})

test_that("plot_proteome_overview preserves verdict and flag semantics", {
  batch <- .fix_batch_trypsin
  p <- plot_proteome_overview(batch, title = "Overview")
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "Overview")
  score_panel <- p$patches$plots[[1L]]
  expect_true(inherits(score_panel, "patchwork"))
  expect_true(any(c("Good", "Moderate", "Poor") %in%
    as.character(score_panel$patches$plots[[1L]]$data$verdict)))

  with_unique <- batch
  with_unique$S_unique <- seq(0.2, 0.8, length.out = nrow(with_unique))
  unique_plot <- plot_proteome_overview(with_unique)
  expect_true("Uniqueness" %in%
    batch_grob_text_labels(patchwork::patchworkGrob(unique_plot)))
})

test_that("plot_proteome_overview rejects malformed scores, verdicts, and flags", {
  batch <- .fix_batch_bsa
  bad <- batch
  bad$composite_score[[1L]] <- Inf
  expect_error(
    plot_proteome_overview(bad),
    class = "pepvet_error_invalid_batch"
  )
  bad <- batch
  bad$verdict[[1L]] <- "Unknown"
  expect_error(
    plot_proteome_overview(bad),
    class = "pepvet_error_invalid_batch"
  )
  bad <- batch
  bad$flag_hydrophobic[[1L]] <- NA
  expect_error(
    plot_proteome_overview(bad),
    class = "pepvet_error_invalid_batch"
  )
  bad <- batch
  bad$protein_id[[1L]] <- ""
  expect_error(
    plot_proteome_overview(bad),
    class = "pepvet_error_invalid_batch"
  )
  expect_error(
    plot_proteome_overview(data.frame(composite_score = 0.5)),
    class = "pepvet_error_invalid_batch"
  )
  expect_error(
    plot_proteome_overview("not a batch"),
    class = "pepvet_error_invalid_batch"
  )
  duplicate <- rbind(batch, batch)
  expect_error(
    plot_proteome_overview(duplicate),
    class = "pepvet_error_invalid_batch"
  )
})

test_that("plot_proteome_overview has an explicit fallback without flags", {
  batch <- .fix_batch_bsa
  batch <- batch[, !names(batch) %in% c(
    "flag_short_protein", "flag_no_valid_peptides",
    "flag_hydrophobic", "flag_low_complexity"
  ), drop = FALSE]
  p <- plot_proteome_overview(batch, title = "No flags")
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "No flags")
})

# plot_batch_comparison.

test_that("plot_batch_comparison returns patchwork from batch_compare_enzymes", {
  comp <- batch_compare_enzymes(
    Biostrings::readAAStringSet(.bsa_path),
    enzymes = c("trypsin", "lysc")
  )
  p <- plot_batch_comparison(comp)
  expect_s3_class(p, "patchwork")
  expect_match(p$patches$annotation$title, "Proteome Enzyme Comparison")
  expect_equal(as.numeric(p$patches$layout$heights), c(1.1, 1))
})

test_that("plot_batch_comparison errors on empty comparison", {
  empty <- data.frame(
    protein_id = character(0), enzyme = character(0),
    composite_score = numeric(0), verdict = character(0),
    S_length = numeric(0), S_coverage = numeric(0),
    S_count = numeric(0), S_hydro = numeric(0),
    S_charge = numeric(0), stringsAsFactors = FALSE
  )
  expect_error(
    plot_batch_comparison(empty),
    class = "pepvet_error_invalid_batch"
  )
  expect_error(
    plot_batch_comparison("not a comparison"),
    class = "pepvet_error_invalid_batch"
  )
})

test_that("plot_batch_comparison validates enzyme grids and score semantics", {
  comp <- batch_compare_enzymes(
    Biostrings::readAAStringSet(.bsa_path),
    enzymes = c("trypsin", "lysc")
  )
  bad <- comp[-1L, , drop = FALSE]
  expect_error(
    plot_batch_comparison(bad),
    class = "pepvet_error_invalid_batch"
  )
  bad <- comp
  bad$verdict[[1L]] <- "Unknown"
  expect_error(
    plot_batch_comparison(bad),
    class = "pepvet_error_invalid_batch"
  )
  bad <- comp
  bad$enzyme[[1L]] <- NA_character_
  expect_error(
    plot_batch_comparison(bad),
    class = "pepvet_error_invalid_batch"
  )
  duplicate <- as.data.frame(comp)
  duplicate <- rbind(duplicate, duplicate[1L, , drop = FALSE])
  expect_error(
    plot_batch_comparison(duplicate),
    class = "pepvet_error_invalid_batch"
  )
})

test_that("plot_batch_comparison supports raw grids, custom titles, and verdict tiers", {
  comp <- batch_compare_enzymes(
    Biostrings::readAAStringSet(.bsa_path),
    enzymes = c("trypsin", "lysc")
  )
  raw <- as.data.frame(comp)
  raw$composite_score <- c(0.50, 0.30)
  raw$verdict <- c("Moderate", "Poor")
  p <- plot_batch_comparison(raw, title = "Raw comparison")
  expect_s3_class(p, "patchwork")
  expect_identical(p$patches$annotation$title, "Raw comparison")
  expect_true(all(c("Moderate", "Poor") %in%
    as.character(p$patches$plots[[1L]]$patches$plots[[1L]]$data$verdict)))

  raw_unique <- rbind(raw, raw)
  second_rows <- seq.int(nrow(raw) + 1L, nrow(raw_unique))
  raw_unique$protein_id[second_rows] <- paste0(
    raw_unique$protein_id[second_rows], "_2"
  )
  raw_unique$S_unique <- seq(0.2, 0.8, length.out = nrow(raw_unique))
  unique_plot <- plot_batch_comparison(raw_unique)
  expect_true("Uniqueness" %in%
    batch_grob_text_labels(patchwork::patchworkGrob(unique_plot)))
})

test_that("batch plots reject duplicate columns and invalid titles", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  duplicate <- as.data.frame(.fix_batch_trypsin)
  names(duplicate)[[2L]] <- names(duplicate)[[1L]]
  expect_error(
    plot_proteome_overview(duplicate),
    class = "pepvet_error_invalid_batch"
  )

  comparison <- batch_compare_enzymes(
    Biostrings::readAAStringSet(.bsa_path),
    enzymes = c("trypsin", "lysc")
  )
  expect_error(
    plot_proteome_overview(.fix_batch_trypsin, title = c("one", "two")),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    plot_batch_comparison(comparison, title = c("one", "two")),
    class = "pepvet_error_invalid_input"
  )
})
