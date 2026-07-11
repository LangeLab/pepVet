# Snapshot tests capture the exact console output.
# Run testthat::snapshot_review("report") after intentional layout changes.

# Single-protein bar reports.

test_that("BSA trypsin report matches snapshot (Good/Moderate verdict)", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"),
    enzyme = "trypsin",
    missed_cleavages = 1L
  )
  expect_snapshot(digest_report(ev))
})

test_that("Histone H3 trypsin report matches snapshot (Moderate verdict)", {
  ev <- evaluate_digest(reference_fasta("P68431.fasta"), enzyme = "trypsin")
  expect_snapshot(digest_report(ev))
})

test_that("BSA lysc report matches snapshot", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"), enzyme = "lysc")
  expect_snapshot(digest_report(ev))
})

test_that("proteome-aware report includes S_unique bar", {
  multi <- reference_fasta("P37840_isoforms.fasta")
  isoforms <- Biostrings::readAAStringSet(multi)
  proteome_digest <- digest_protein(multi)
  ev <- evaluate_digest(isoforms[1], proteome = proteome_digest)
  # S_unique must appear in output
  output <- capture.output(digest_report(ev))
  expect_true(any(grepl("S_unique", output)))
})

test_that("report output uses the short S_charge label with clarified meaning", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  output <- capture.output(digest_report(ev))

  expect_true(any(grepl("^  S_charge\\s", output)))
})

test_that("report format helpers handle boundaries and missing display values", {
  expect_identical(
    .format_score(c(0, 1, NA_real_)),
    c("0.000", "1.000", "NA")
  )
  expect_match(
    cli::ansi_strip(.format_component_bar(0)),
    "0.000",
    fixed = TRUE
  )
  expect_match(
    cli::ansi_strip(.format_component_bar(1)),
    "1.000",
    fixed = TRUE
  )
  expect_match(
    cli::ansi_strip(.format_component_bar(NA_real_)),
    "NA",
    fixed = TRUE
  )
  expect_match(
    cli::ansi_strip(.format_component_bar(2, width = 0L)),
    "2.000",
    fixed = TRUE
  )
  expect_identical(.verdict_bullet("unrecognised"), cli::symbol$info)
  expect_identical(.verdict_bullet(NA_character_), cli::symbol$info)
  expect_identical(
    cli::ansi_strip(.verdict_bullet("Good")),
    cli::ansi_strip(cli::symbol$tick)
  )
  expect_identical(
    cli::ansi_strip(.verdict_bullet("Moderate")),
    cli::ansi_strip(cli::symbol$warning)
  )
  expect_identical(
    cli::ansi_strip(.verdict_bullet("Poor")),
    cli::ansi_strip(cli::symbol$cross)
  )
  expect_identical(cli::ansi_strip(.verdict_colour("Good", "score")), "score")
  expect_identical(
    cli::ansi_strip(.verdict_colour("Moderate", "score")), "score"
  )
  expect_identical(cli::ansi_strip(.verdict_colour("Poor", "score")), "score")
  expect_identical(
    .verdict_colour("unrecognised", "score"), "score"
  )
  expect_identical(.verdict_colour(NA_character_, "score"), "score")
  expect_error(
    .format_component_bar(0.5, width = -1L),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    .format_component_bar("0.5"),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    .format_component_bar(NaN),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    .format_component_bar(0.5, width = "10"),
    class = "pepvet_error_invalid_report_input"
  )
  expect_match(
    cli::ansi_strip(.format_component_bar(Inf)),
    "Inf",
    fixed = TRUE
  )
  expect_match(
    cli::ansi_strip(.format_component_bar(-Inf)),
    "-Inf",
    fixed = TRUE
  )
})

# Multi-enzyme comparison reports.

test_that("multi-enzyme comparison report matches snapshot", {
  comp <- compare_digests(
    reference_fasta("P02769.fasta"),
    enzymes = c(
      "trypsin", "lysc",
      "glutamyl endopeptidase", "asp-n endopeptidase"
    )
  )
  expect_snapshot(digest_report(comp))
})

test_that("Histone H3 multi-enzyme comparison report matches snapshot", {
  comp <- compare_digests(
    reference_fasta("P68431.fasta"),
    enzymes = c("trypsin", "lysc", "asp-n endopeptidase")
  )
  expect_snapshot(digest_report(comp))
})

# Behavioural invariants.

test_that("digest_report returns its input invisibly", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"))
  result <- withVisible(digest_report(ev))
  expect_false(result$visible)
  expect_identical(result$value, ev)
})

test_that("score values appear in the evaluate_digest report output", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"))
  output <- capture.output(digest_report(ev))
  composite_text <- formatC(ev$scores$composite_score, format = "f", digits = 3)
  expect_true(any(grepl(composite_text, output, fixed = TRUE)))
})

test_that("custom title overrides the protein ID header", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"))
  output <- capture.output(digest_report(ev, title = "My Custom Title"))
  expect_true(any(grepl("My Custom Title", output, fixed = TRUE)))
})

test_that("comparison reports identify the actual best enzyme after reordering", {
  comparison <- .bsa_comparison(c("trypsin", "lysc", "chymotrypsin-high"))
  comparison <- comparison[rev(seq_len(nrow(comparison))), , drop = FALSE]
  expected_best <- comparison$enzyme[[which.max(comparison$composite_score)]]
  output <- cli::ansi_strip(capture.output(
    digest_report(comparison, title = "reordered comparison")
  ))

  expect_true(any(grepl(
    paste0("best: ", expected_best), output, fixed = TRUE
  )))
  marked <- grep("^>", output, value = TRUE)
  expect_length(marked, 1L)
  expect_match(marked, expected_best, fixed = TRUE)
})

# Edge cases.

test_that("digest_report handles a single-peptide protein without error", {
  # MKWVTFISLLFLFSSAYSR has one tryptic peptide (no internal K/R)
  ev <- evaluate_digest("MKWVTFISLLFLFSSAYSR", enzyme = "trypsin")
  output <- capture.output(result <- digest_report(ev))
  expect_identical(result, ev)
  expect_true(any(grepl("S_charge", output, fixed = TRUE)))
})

test_that("digest_report separates multiple evaluated score rows", {
  ev <- .fix_bsa_trypsin
  ev$scores <- rbind(ev$scores, ev$scores)

  output <- capture.output(result <- digest_report(ev))

  expect_identical(result, ev)
  expect_equal(sum(grepl("S_length", output, fixed = TRUE)), 2L)
})

test_that("digest_report rejects invalid input with a classed error", {
  expect_error(
    digest_report(NULL),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    digest_report(list(a = 1, b = 2)),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    digest_report("not a result"),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores$S_length[[1L]] <- NA_real_
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores <- malformed$scores[, setdiff(
    names(malformed$scores), "S_charge"
  ), drop = FALSE]
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores$S_length <- as.character(malformed$scores$S_length)
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores$S_length[[1L]] <- 1.1
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$params$protein_ids <- NULL
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$params$enzyme <- NA_character_
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$params$enzyme <- 1L
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$params$missed_cleavages <- 1.5
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores$protein_id <- "unknown"
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores$protein_id <- ""
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  names(malformed$scores)[[2L]] <- names(malformed$scores)[[1L]]
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$peptides <- NULL
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  malformed <- .fix_bsa_trypsin
  malformed$scores <- NULL
  expect_error(
    digest_report(malformed),
    class = "pepvet_error_invalid_report_input"
  )

  comparison <- .bsa_comparison()
  comparison$verdict[[1L]] <- "Unknown"
  expect_error(
    digest_report(comparison),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    digest_report(comparison[0, , drop = FALSE]),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    digest_report(comparison[, setdiff(names(comparison), "S_charge"),
      drop = FALSE]),
    class = "pepvet_error_invalid_report_input"
  )

  comparison <- .bsa_comparison()
  comparison$enzyme[[1L]] <- ""
  expect_error(
    digest_report(comparison),
    class = "pepvet_error_invalid_report_input"
  )

  comparison <- .bsa_comparison()
  comparison$protein_id[[2L]] <- "another-protein"
  expect_error(
    digest_report(comparison),
    class = "pepvet_error_invalid_report_input"
  )
})

test_that("report titles must be scalar character values", {
  ev <- .fix_bsa_trypsin

  for (invalid in list(integer(0), NA_character_, c("a", "b"), factor("a"))) {
    expect_error(
      digest_report(ev, title = invalid),
      class = "pepvet_error_invalid_report_input"
    )
  }
})

test_that("pepvet_check forwards scoring and returns a complete result", {
  result <- pepvet_check(
    .bsa_path,
    enzyme = "trypsin",
    gravy_range = c(-2, 2),
    length_range = c(6L, 30L)
  )

  expect_s3_class(result$scores, "tbl_df")
  expect_equal(result$params$gravy_range, c(-2, 2))
  expect_identical(result$params$length_range, c(6L, 30L))
})

test_that("pepvet_check prints the same report as evaluate plus digest_report", {
  expected <- evaluate_digest(.bsa_path, enzyme = "trypsin")
  expected_output <- cli::ansi_strip(capture.output(digest_report(expected)))
  actual_output <- cli::ansi_strip(capture.output(
    pepvet_check(.bsa_path, enzyme = "trypsin")
  ))

  expect_identical(actual_output, expected_output)
})

# pepvet_check.

test_that("pepvet_check returns evaluate_digest result invisibly and unchanged", {
  bsa_path <- reference_fasta("P02769.fasta")
  ev <- evaluate_digest(bsa_path, enzyme = "trypsin")
  result <- withVisible(pepvet_check(bsa_path, enzyme = "trypsin"))

  expect_false(result$visible)
  expect_identical(result$value$scores, ev$scores)
  expect_identical(result$value$peptides, ev$peptides)
})

test_that("pepvet_check passes ... arguments through to evaluate_digest", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- pepvet_check(bsa_path, enzyme = "trypsin", missed_cleavages = 1L)

  expect_equal(result$params$missed_cleavages, 1L)
})
