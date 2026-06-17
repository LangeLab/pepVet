# ── Cross-function consistency ─────────────────────────────────────────────

test_that("evaluate_digest gives the same result as manual pipeline", {
  result <- evaluate_digest(.bsa_path,
    enzyme = "trypsin",
    missed_cleavages = 0L
  )

  manual_peptides <- digest_protein(.bsa_path,
    enzyme = "trypsin",
    missed_cleavages = 0L
  )
  manual_scores <- score_peptides(manual_peptides)
  manual_annotations <- annotate_cleavage_sites(.bsa_path, enzyme = "trypsin")
  manual_scores <- tibble::add_column(
    manual_scores,
    n_high_efficiency_sites = sum(manual_annotations$efficiency == "high"),
    n_low_efficiency_sites = sum(manual_annotations$efficiency == "low"),
    .after = "preset_used"
  )

  expect_equal(result$peptides, manual_peptides, tolerance = 1e-15)
  expect_equal(result$scores, manual_scores, tolerance = 1e-15)
})

test_that("batch_evaluate returns a tibble with one row per protein and required columns", {
  batch <- .fix_batch_small

  sequences <- Biostrings::readAAStringSet(.small_path)
  expect_s3_class(batch, "tbl_df")
  expect_equal(nrow(batch), length(sequences))
  expect_true(
    all(
      c(
        "protein_id", "protein_length", "n_peptides", "n_valid_peptides",
        "composite_score", "verdict", "median_peptide_length",
        "flag_short_protein", "flag_hydrophobic",
        "flag_low_complexity", "flag_no_valid_peptides"
      ) %in% names(batch)
    )
  )
})

test_that("batch_evaluate composite_score and verdict match evaluate_digest for the same protein", {
  batch <- .fix_batch_bsa
  individual <- .fix_bsa_trypsin

  expect_equal(batch$composite_score[[1L]], individual$scores$composite_score)
  expect_equal(batch$verdict[[1L]], individual$scores$verdict)
})

test_that("batch_evaluate includes S_unique column when proteome is provided", {
  proteome_digest <- digest_protein(.multi_path)

  batch_with <- batch_evaluate(.multi_path, proteome = proteome_digest)
  batch_without <- .fix_batch_multi

  expect_true("S_unique" %in% names(batch_with))
  expect_false("S_unique" %in% names(batch_without))
})

test_that("evaluate_digest passes include_pI through to score output", {
  result <- evaluate_digest(.bsa_path, enzyme = "trypsin", include_pI = TRUE)

  expect_true("pI" %in% names(result$scores))
  expect_type(result$scores$pI, "list")
})

test_that("evaluate_digest can append peptide-level cleavage efficiency", {
  result <- evaluate_digest(
    "AKRTPK",
    enzyme = "trypsin",
    missed_cleavages = 0L,
    include_cleavage_efficiency = TRUE
  )

  expect_true("cleavage_efficiency" %in% names(result$peptides))
  expect_identical(result$peptides$cleavage_efficiency, c("medium", "medium", "high"))
})

# ── Comparison & recommendation ────────────────────────────────────────────

test_that("compare_digests output is sorted by composite_score descending", {
  result <- compare_digests(.bsa_path, enzymes = c("trypsin", "lysc"))

  expect_s3_class(result, "tbl_df")
  expect_true(
    all(diff(result$composite_score) <= 0),
    info = "composite_score must be non-increasing across rows"
  )
})

test_that("compare_digests output has enzyme column plus all score columns", {
  result <- compare_digests(.bsa_path, enzymes = c("trypsin", "lysc"))

  expect_identical(names(result)[[1L]], "enzyme")
  expected_score_cols <- c(
    "protein_id", "S_length", "S_coverage", "S_count",
    "S_hydro", "S_charge", "composite_score", "verdict",
    "median_peptide_length", "preset_used", "n_high_efficiency_sites",
    "n_low_efficiency_sites"
  )
  expect_true(all(expected_score_cols %in% names(result)))
  expect_identical(nrow(result), 2L)
})

test_that("subsetted batch_compare_enzymes objects print without error", {
  result <- suppressWarnings(
    batch_compare_enzymes(.small_path, enzymes = c("trypsin", "lysc"))
  )
  subsetted <- result[result$enzyme == "trypsin", c("protein_id", "composite_score")]

  expect_no_error(print(subsetted))
})

test_that("compare_digests rejects multi-protein input", {
  expect_error(
    compare_digests(.multi_path, enzymes = c("trypsin", "lysc")),
    class = "pepvet_error_invalid_input"
  )
})

test_that("recommend_enzyme selects trypsin for BSA at one missed cleavage", {
  # At missed_cleavages = 0, trypsin over-digests BSA (many short sub-7 AA K/R
  # peptides), and lysc wins on S_length. With missed_cleavages = 1, the merged
  # spans are in the valid range and trypsin's higher peptide yield wins out.
  # mc = 1 also reflects how BSA is used as a calibration standard in practice.
  result <- recommend_enzyme(
    .bsa_path,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 1L
  )

  expect_identical(result, "trypsin")
})

test_that("recommend_enzyme does not select trypsin for Histone H3", {
  result <- recommend_enzyme(.h3_path, enzymes = c("trypsin", "lysc"))

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
  result <- .fix_bsa_trypsin

  expect_type(result, "list")
  expect_identical(names(result), c("scores", "peptides", "params"))
  expect_s3_class(result$scores, "tbl_df")
  expect_s3_class(result$peptides, "tbl_df")
  expect_type(result$params, "list")
  expect_identical(
    names(result$params),
    c("enzyme", "missed_cleavages", "protein_ids", "preset_used")
  )
  expect_identical(result$params$preset_used, "standard")
})

test_that("evaluate_digest records cleavage-site counts for trypsin-family digests", {
  result <- .fix_bsa_trypsin
  annotations <- annotate_cleavage_sites(.bsa_path, enzyme = "trypsin")

  expect_identical(
    result$scores$n_high_efficiency_sites,
    sum(annotations$efficiency == "high")
  )
  expect_identical(
    result$scores$n_low_efficiency_sites,
    sum(annotations$efficiency == "low")
  )
})

test_that("unsupported enzymes get NA cleavage-site counts", {
  result <- evaluate_digest("AKRTPK", enzyme = "lysc")

  expect_true(all(is.na(result$scores$n_high_efficiency_sites)))
  expect_true(all(is.na(result$scores$n_low_efficiency_sites)))
})

test_that("params reflects the resolved enzyme name and missed_cleavages", {
  result <- evaluate_digest(.bsa_path, enzyme = "Trypsin", missed_cleavages = 1L)

  expect_identical(result$params$enzyme, "trypsin")
  expect_identical(result$params$missed_cleavages, 1L)
  expect_type(result$params$protein_ids, "character")
  expect_identical(result$params$preset_used, "standard")
})

test_that("evaluate_digest records preset_used in params for named presets", {
  result <- do.call(
    evaluate_digest,
    c(list(sequence = .bsa_path, enzyme = "trypsin"), pepvet_preset("fractionated"))
  )

  expect_identical(result$params$preset_used, "fractionated")
})

test_that("evaluate_digest peptides matches direct digest_protein output", {
  result <- .fix_bsa_lysc_mc1
  direct <- digest_protein(.bsa_path, enzyme = "lysc", missed_cleavages = 1L)

  expect_identical(result$peptides, direct)
})

test_that("protein_id is preserved across scores, peptides, and params", {
  result <- .fix_bsa_trypsin

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

# ── summarize_batch ─────────────────────────────────────────────────────────

test_that("summarize_batch returns a list with expected element names", {
  batch <- .fix_batch_small
  summary <- summarize_batch(batch)

  expect_type(summary, "list")
  expect_setequal(
    names(summary),
    c(
      "verdict_counts", "score_distribution", "component_summary",
      "problem_proteins", "enzyme_switch_candidates"
    )
  )
})

test_that("summarize_batch verdict_counts n sums to number of proteins", {
  batch <- .fix_batch_small
  summary <- summarize_batch(batch)

  expect_equal(sum(summary$verdict_counts$n), nrow(batch))
})

test_that("summarize_batch verdict_counts has the three verdict levels", {
  batch <- .fix_batch_bsa
  summary <- summarize_batch(batch)

  expect_equal(summary$verdict_counts$verdict, c("Good", "Moderate", "Poor"))
})

test_that("summarize_batch score_distribution has expected statistic names", {
  batch <- .fix_batch_small
  summary <- summarize_batch(batch)

  expect_named(
    summary$score_distribution,
    c("mean", "median", "sd", "q25", "q75", "min", "max")
  )
  expect_true(all(is.finite(summary$score_distribution)))
})

test_that("summarize_batch component_summary contains the five core components", {
  batch <- .fix_batch_bsa
  summary <- summarize_batch(batch)

  expect_true(
    all(
      c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge") %in%
        names(summary$component_summary)
    )
  )
})

test_that("summarize_batch problem_proteins is a tibble ordered by score", {
  batch <- .fix_batch_small
  summary <- summarize_batch(batch)

  expect_s3_class(summary$problem_proteins, "tbl_df")
  scores <- summary$problem_proteins$composite_score
  expect_true(all(diff(scores) >= 0))
})

test_that("summarize_batch rejects a non-tibble input with a classed error", {
  expect_error(
    summarize_batch("not a tibble"),
    class = "pepvet_error_invalid_batch_result"
  )
})

test_that("summarize_batch rejects an empty tibble with a classed error", {
  expect_error(
    summarize_batch(tibble::tibble()),
    class = "pepvet_error_invalid_batch_result"
  )
})

# ── triage_proteins ─────────────────────────────────────────────────────────

test_that("triage_proteins returns a tibble with an action column", {
  batch <- .fix_batch_small
  triaged <- triage_proteins(batch)

  expect_s3_class(triaged, "tbl_df")
  expect_true("action" %in% names(triaged))
})

test_that("triage_proteins action values are from the expected set", {
  batch <- .fix_batch_small
  triaged <- triage_proteins(batch)

  valid_actions <- c(
    "proceed", "consider_alternative",
    "try_other_enzyme", "skip"
  )
  expect_true(all(triaged$action %in% valid_actions))
})

test_that("triage_proteins returns one row per protein", {
  batch <- .fix_batch_small
  triaged <- triage_proteins(batch)

  expect_equal(nrow(triaged), nrow(batch))
})

test_that("triage_proteins categorizes BSA trypsin (mc=1) as proceed", {
  batch <- .fix_batch_bsa_mc1
  triaged <- triage_proteins(batch)

  expect_equal(triaged$action[[1]], "proceed")
})

test_that("triage_proteins categorizes Histone H3.1 trypsin as try_other_enzyme", {
  batch <- batch_evaluate(system.file("extdata", "P68431.fasta", package = "pepVet"),
    enzyme = "trypsin", missed_cleavages = 0L
  )
  triaged <- triage_proteins(batch)

  expect_equal(triaged$action[[1]], "try_other_enzyme")
})

test_that("triage_proteins flat tibble contains expected score columns", {
  batch <- .fix_batch_bsa
  triaged <- triage_proteins(batch)

  expect_true(
    all(
      c(
        "protein_id", "protein_length", "n_peptides", "n_valid_peptides",
        "composite_score", "verdict",
        "flag_short_protein", "flag_hydrophobic",
        "flag_low_complexity", "flag_no_valid_peptides"
      ) %in% names(triaged)
    )
  )
})

# ── Weight sensitivity analysis ──────────────────────────────────────────────

test_that("sensitivity_analysis returns correct structure for evaluate_digest input", {
  res <- .fix_bsa_trypsin
  withr::local_seed(42)
  sens <- sensitivity_analysis(res, n_iter = 500L)
  expect_type(sens, "list")
  expect_named(sens, c("iterations", "convergence", "summary"))
  expect_s3_class(sens$iterations, "tbl_df")
  expect_true("composite_score" %in% names(sens$iterations))
  expect_true("verdict" %in% names(sens$iterations))
  expect_equal(nrow(sens$iterations), 500L)
  expect_true(sens$summary$verdict_pct["Good"] > 0.95)
  expect_length(sens$summary$composite_ci, 2L)
})

test_that("sensitivity_analysis importance returns R2 values", {
  res <- .fix_bsa_trypsin
  withr::local_seed(42)
  sens <- sensitivity_analysis(res, n_iter = 500L, importance = TRUE)
  expect_named(sens$summary$weight_importance,
    c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge"))
  expect_true(all(sens$summary$weight_importance >= 0))
  expect_true(all(sens$summary$weight_importance <= 1))
})

test_that("sensitivity_analysis corner_cases returns table", {
  res <- .fix_bsa_trypsin
  withr::local_seed(42)
  sens <- sensitivity_analysis(res, n_iter = 500L, corner_cases = TRUE)
  expect_s3_class(sens$summary$corner_cases, "tbl_df")
  expect_true("composite_at_lo" %in% names(sens$summary$corner_cases))
  expect_true("composite_at_hi" %in% names(sens$summary$corner_cases))
})

test_that("sensitivity_analysis verdict matches expectation for zero-cleavage protein", {
  no_cleave <- "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  res <- evaluate_digest(no_cleave, enzyme = "trypsin", missed_cleavages = 1L)
  expect_identical(res$scores$verdict, "Poor")
  expect_equal(res$scores$composite_score, 0)
})

test_that("sensitivity_analysis on multi-enzyme input returns rank stability", {
  bsa_path <- .bsa_path
  comp <- compare_digests(bsa_path,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 1L
  )
  withr::local_seed(42)
  sens <- sensitivity_analysis(comp, n_iter = 500L)
  expect_named(sens$summary, c("top1_stability", "kendall_mean"))
  expect_true("trypsin" %in% names(sens$summary$top1_stability))
  expect_true(sens$summary$kendall_mean > 0)
})

test_that("sensitivity_analysis on batch input returns per-protein instability", {
  batch <- .fix_batch_small
  withr::local_seed(42)
  sens <- sensitivity_analysis(batch, n_iter = 500L, chunk_size = 50L)
  expect_named(sens, c("per_protein", "summary"))
  expect_s3_class(sens$per_protein, "tbl_df")
  expect_true("verdict_instability" %in% names(sens$per_protein))
  expect_true("protein_id" %in% names(sens$per_protein))
  expect_equal(nrow(sens$per_protein), nrow(batch))
  expect_named(sens$summary, c("total_instability", "mean_ci_width", "ci_width_quantiles"))
})

test_that("Dirichlet sampler produces correct mean and variance", {
  w0 <- c(a = 0.2, b = 0.3, c = 0.5)
  nu <- 100
  set.seed(42)
  samples <- .rdirichlet(10000, nu * w0)
  expect_equal(nrow(samples), 10000)
  expect_equal(ncol(samples), 3L)
  expect_equal(rowSums(samples), rep(1, 10000), tolerance = 1e-10)
  col_means <- colMeans(samples)
  expect_equal(as.numeric(col_means), as.numeric(w0), tolerance = 0.02)
  col_vars <- apply(samples, 2, var)
  expected_var <- as.numeric(w0 * (1 - w0) / (nu + 1))
  expect_equal(as.numeric(col_vars), expected_var, tolerance = 0.005)
})

test_that("plot_weight_sensitivity returns a ggplot for single-protein input", {
  skip_if_not_installed("ggplot2")
  res <- .fix_bsa_trypsin
  sens <- sensitivity_analysis(res, n_iter = 500L)
  p <- plot_weight_sensitivity(sens)
  expect_s3_class(p, "ggplot")
})

# ── ... argument passthrough ──────────────────────────────────────────────

test_that("evaluate_digest passes ... to score_peptides", {
  bsa <- .bsa_path
  default <- evaluate_digest(bsa, enzyme = "trypsin")
  wider <- evaluate_digest(bsa, enzyme = "trypsin", gravy_range = c(-2.0, 2.0))
  expect_true(wider$scores$S_hydro >= default$scores$S_hydro)
})

test_that("compare_digests passes ... to evaluate_digest", {
  bsa <- .bsa_path
  default <- compare_digests(bsa, enzymes = c("trypsin", "lysc"))
  wider <- compare_digests(bsa, enzymes = c("trypsin", "lysc"),
    gravy_range = c(-2.0, 2.0))
  expect_true(wider$S_hydro[[1]] >= default$S_hydro[[1]])
})

test_that("recommend_enzyme passes ... to compare_digests", {
  bsa <- .bsa_path
  default <- recommend_enzyme(bsa, enzymes = c("trypsin", "lysc"))
  wider <- recommend_enzyme(bsa, enzymes = c("trypsin", "lysc"),
    gravy_range = c(-2.0, 2.0))
  expect_type(default, "character")
  expect_type(wider, "character")
})

test_that("batch_evaluate passes ... to score_peptides", {
  small <- .small_path
  default <- batch_evaluate(small, enzyme = "trypsin")
  wider <- batch_evaluate(small, enzyme = "trypsin", gravy_range = c(-2.0, 2.0))
  expect_true(wider$S_hydro[[1]] >= default$S_hydro[[1]])
})

test_that("batch_compare_enzymes passes ... to batch_evaluate", {
  small <- .small_path
  default <- batch_compare_enzymes(small, enzymes = c("trypsin", "lysc"))
  wider <- batch_compare_enzymes(small, enzymes = c("trypsin", "lysc"),
    gravy_range = c(-2.0, 2.0))
  expect_true(wider$S_hydro[[1]] >= default$S_hydro[[1]])
})

test_that("pepvet_check passes ... to evaluate_digest", {
  bsa <- .bsa_path
  default <- pepvet_check(bsa, enzyme = "trypsin")
  wider <- pepvet_check(bsa, enzyme = "trypsin", gravy_range = c(-2.0, 2.0))
  expect_true(wider$scores$S_hydro >= default$scores$S_hydro)
})

# ── Error class tests ─────────────────────────────────────────────────────

test_that("batch_evaluate rejects invalid cores", {
  expect_error(
    batch_evaluate(.bsa_path, cores = 0L),
    class = "pepvet_error_invalid_cores"
  )
  expect_error(
    batch_evaluate(.bsa_path, cores = -1L),
    class = "pepvet_error_invalid_cores"
  )
})

test_that("compare_digests rejects invalid enzymes", {
  expect_error(
    compare_digests(.bsa_path, enzymes = NULL),
    class = "pepvet_error_invalid_enzymes"
  )
  expect_error(
    compare_digests(.bsa_path, enzymes = NA_character_),
    class = "pepvet_error_invalid_enzymes"
  )
})
