# Tests aligned with the corresponding plotting source domain.
test_that("plot_proteome_overview returns patchwork from batch_evaluate", {
  batch <- batch_evaluate(
    Biostrings::readAAStringSet(.bsa_path),
    enzyme = "trypsin"
  )
  p <- plot_proteome_overview(batch)
  expect_s3_class(p, "patchwork")
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

# plot_batch_comparison.

test_that("plot_batch_comparison returns patchwork from batch_compare_enzymes", {
  comp <- batch_compare_enzymes(
    Biostrings::readAAStringSet(.bsa_path),
    enzymes = c("trypsin", "lysc")
  )
  p <- plot_batch_comparison(comp)
  expect_s3_class(p, "patchwork")
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
})
