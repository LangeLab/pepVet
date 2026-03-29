reference_fasta <- function(file_name) {
  system.file("extdata", file_name, package = "pepVet")
}

# ── Cross-function consistency ─────────────────────────────────────────────

test_that("evaluate_digest gives the same result as manual pipeline", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- evaluate_digest(bsa_path, enzyme = "trypsin")

  manual_peptides <- digest_protein(bsa_path, enzyme = "trypsin")
  manual_scores <- score_peptides(manual_peptides)

  expect_identical(result$peptides, manual_peptides)
  expect_identical(result$scores, manual_scores)
})

test_that("batch_evaluate matches individual evaluate_digest calls exactly", {
  small_path <- reference_fasta("small_proteome_50_proteins.fasta")
  batch <- batch_evaluate(small_path, enzyme = "trypsin")

  sequences <- Biostrings::readAAStringSet(small_path)

  for (protein_id in names(sequences)[1:3]) {
    individual <- evaluate_digest(sequences[protein_id], enzyme = "trypsin")
    expect_identical(batch[[protein_id]], individual)
  }
})

test_that("batch_evaluate works for a single protein", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- batch_evaluate(bsa_path, enzyme = "trypsin")

  expect_length(result, 1L)
  expect_type(result, "list")
  expect_s3_class(result[[1]]$scores, "tbl_df")
})

test_that("batch_evaluate threads the proteome argument to every evaluation", {
  multi_path <- reference_fasta("P37840_isoforms.fasta")
  proteome_digest <- digest_protein(multi_path)

  batch_with_proteome <- batch_evaluate(multi_path, proteome = proteome_digest)
  batch_no_proteome <- batch_evaluate(multi_path)

  expect_true("S_unique" %in% names(batch_with_proteome[[1]]$scores))
  expect_false("S_unique" %in% names(batch_no_proteome[[1]]$scores))
})

test_that("evaluate_digest passes include_pI through to score output", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- evaluate_digest(bsa_path, enzyme = "trypsin", include_pI = TRUE)

  expect_true("pI" %in% names(result$scores))
  expect_type(result$scores$pI, "list")
})

# ── Comparison & recommendation ────────────────────────────────────────────

test_that("compare_digests output is sorted by composite_score descending", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- compare_digests(bsa_path, enzymes = c("trypsin", "lysc"))

  expect_s3_class(result, "tbl_df")
  expect_true(
    all(diff(result$composite_score) <= 0),
    info = "composite_score must be non-increasing across rows"
  )
})

test_that("compare_digests output has enzyme column plus all score columns", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- compare_digests(bsa_path, enzymes = c("trypsin", "lysc"))

  expect_identical(names(result)[[1L]], "enzyme")
  expected_score_cols <- c(
    "protein_id", "S_length", "S_coverage", "S_count",
    "S_hydro", "S_charge", "composite_score", "verdict",
    "median_peptide_length"
  )
  expect_true(all(expected_score_cols %in% names(result)))
  expect_identical(nrow(result), 2L)
})

test_that("compare_digests rejects multi-protein input", {
  multi_path <- reference_fasta("P37840_isoforms.fasta")
  expect_error(
    compare_digests(multi_path, enzymes = c("trypsin", "lysc")),
    class = "pepvet_error_invalid_input"
  )
})

test_that("recommend_enzyme selects trypsin for BSA at one missed cleavage", {
  # At missed_cleavages = 0, trypsin over-digests BSA (many short sub-7 AA K/R
  # peptides), and lysc wins on S_length. With missed_cleavages = 1, the merged
  # spans are in the valid range and trypsin's higher peptide yield wins out.
  # mc = 1 also reflects how BSA is used as a calibration standard in practice.
  bsa_path <- reference_fasta("P02769.fasta")
  result <- recommend_enzyme(
    bsa_path,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 1L
  )

  expect_identical(result, "trypsin")
})

test_that("recommend_enzyme does not select trypsin for Histone H3", {
  h3_path <- reference_fasta("P68431.fasta")
  result <- recommend_enzyme(h3_path, enzymes = c("trypsin", "lysc"))

  expect_false("trypsin" %in% result)
})

test_that("recommend_enzyme returns all tied enzymes in alphabetical order", {
  # A poly-alanine sequence has no trypsin or lysc cut sites; both return an
  # identical single-peptide digest and receive the same composite score.
  poly_ala <- strrep("A", 20L)
  result <- suppressWarnings(
    recommend_enzyme(poly_ala, enzymes = c("trypsin", "lysc"))
  )

  expect_type(result, "character")
  expect_identical(result, c("lysc", "trypsin"))
})

# ── Return structure ───────────────────────────────────────────────────────

test_that("evaluate_digest returns named list with scores, peptides, params", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- evaluate_digest(bsa_path)

  expect_type(result, "list")
  expect_identical(names(result), c("scores", "peptides", "params"))
  expect_s3_class(result$scores, "tbl_df")
  expect_s3_class(result$peptides, "tbl_df")
  expect_type(result$params, "list")
  expect_identical(
    names(result$params),
    c("enzyme", "missed_cleavages", "protein_ids")
  )
})

test_that("params reflects the resolved enzyme name and missed_cleavages", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- evaluate_digest(bsa_path, enzyme = "Trypsin", missed_cleavages = 1L)

  expect_identical(result$params$enzyme, "trypsin")
  expect_identical(result$params$missed_cleavages, 1L)
  expect_type(result$params$protein_ids, "character")
})

test_that("evaluate_digest peptides matches direct digest_protein output", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- evaluate_digest(bsa_path, enzyme = "lysc", missed_cleavages = 1L)
  direct <- digest_protein(bsa_path, enzyme = "lysc", missed_cleavages = 1L)

  expect_identical(result$peptides, direct)
})

test_that("protein_id is preserved across scores, peptides, and params", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- evaluate_digest(bsa_path)

  expect_identical(result$scores$protein_id, result$params$protein_ids)
  expect_true(all(result$peptides$protein_id %in% result$params$protein_ids))
})

# ── Error handling ─────────────────────────────────────────────────────────

test_that("invalid sequence in batch_evaluate propagates a classed error", {
  expect_error(
    batch_evaluate("MXBZ123"),
    class = "pepvet_error_invalid_sequence"
  )
})

test_that("empty AAStringSet in batch_evaluate raises a classed error", {
  expect_error(
    batch_evaluate(Biostrings::AAStringSet()),
    class = "pepvet_error_invalid_input"
  )
})
