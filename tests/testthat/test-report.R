reference_fasta <- function(file_name) {
  system.file("extdata", file_name, package = "pepVet")
}

# Snapshot tests capture the exact console output.
# Run testthat::snapshot_review("report") after intentional layout changes.

# ── Single-protein bar reports ────────────────────────────────────────────

test_that("BSA trypsin report matches snapshot (Good/Moderate verdict)", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"),
    enzyme = "trypsin",
    missed_cleavages = 1L
  )
  expect_snapshot(digest_report(ev))
})

test_that("Histone H3 trypsin report matches snapshot (Poor verdict)", {
  ev <- evaluate_digest(reference_fasta("P68431.fasta"), enzyme = "trypsin")
  expect_snapshot(digest_report(ev))
})

test_that("BSA lysc report matches snapshot", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"), enzyme = "lysc")
  expect_snapshot(digest_report(ev))
})

test_that("proteome-aware report includes S_unique bar", {
  multi <- reference_fasta("P37840_isoforms.fasta")
  proteome_digest <- digest_protein(multi)
  ev <- evaluate_digest(multi[1], proteome = proteome_digest)
  # S_unique must appear in output
  output <- capture.output(digest_report(ev))
  expect_true(any(grepl("S_unique", output)))
})

test_that("report output uses the short S_charge label with clarified meaning", {
  ev <- evaluate_digest(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  output <- capture.output(digest_report(ev))

  expect_true(any(grepl("^  S_charge\\s", output)))
})

# ── Multi-enzyme comparison reports ───────────────────────────────────────────

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

# ── Behavioural invariants ───────────────────────────────────────────────────

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

# ── Edge cases ───────────────────────────────────────────────────────────────

test_that("digest_report handles a single-peptide protein without error", {
  # MKWVTFISLLFLFSSAYSR has one tryptic peptide (no internal K/R)
  ev <- evaluate_digest("MKWVTFISLLFLFSSAYSR", enzyme = "trypsin")
  expect_no_error(digest_report(ev))
})

test_that("digest_report rejects invalid input with a classed error", {
  expect_error(
    digest_report(list(a = 1, b = 2)),
    class = "pepvet_error_invalid_report_input"
  )
  expect_error(
    digest_report("not a result"),
    class = "pepvet_error_invalid_report_input"
  )
})

# ── pepvet_check ──────────────────────────────────────────────────────────────

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
