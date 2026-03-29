reference_fasta <- function(file_name) {
  system.file("extdata", file_name, package = "pepVet")
}

make_digest_result <- function(peptides,
                               starts = NULL,
                               protein_id = "sequence_1",
                               missed_cleavages = 0L) {
  lengths <- nchar(peptides, type = "chars")

  if (is.null(starts)) {
    starts <- c(1L, head(cumsum(as.integer(lengths)) + 1L, -1L))
  }

  if (length(protein_id) == 1L) {
    protein_id <- rep(protein_id, length(peptides))
  }

  if (length(missed_cleavages) == 1L) {
    missed_cleavages <- rep(as.integer(missed_cleavages), length(peptides))
  }

  tibble::tibble(
    protein_id = protein_id,
    peptide = peptides,
    start = as.integer(starts),
    end = as.integer(starts + lengths - 1L),
    length = as.integer(lengths),
    missed_cleavages = as.integer(missed_cleavages)
  )
}

combine_digest_results <- function(...) {
  tibble::as_tibble(do.call(rbind, list(...)))
}

expected_protein_only_weights <- c(
  S_length = 0.25,
  S_coverage = 0.25,
  S_count = 0.20,
  S_hydro = 0.15,
  S_charge = 0.15
)

expected_proteome_weights <- c(
  S_length = 0.20,
  S_coverage = 0.20,
  S_count = 0.15,
  S_hydro = 0.15,
  S_charge = 0.10,
  S_unique = 0.20
)

reference_files <- c(
  "P02769.fasta",
  "P00698.fasta",
  "P56817.fasta",
  "Q8WZ42.fasta",
  "P0CG48.fasta",
  "P37840_isoforms.fasta",
  "P68431.fasta"
)

reference_enzymes <- c(
  "trypsin",
  "lysc",
  "glutamyl endopeptidase",
  "asp-n endopeptidase",
  "chymotrypsin-high"
)

test_that("weight validation returns the documented defaults", {
  expect_equal(
    pepVet:::.validate_weights(NULL, has_proteome = FALSE),
    expected_protein_only_weights
  )
  expect_equal(
    pepVet:::.validate_weights(NULL, has_proteome = TRUE),
    expected_proteome_weights
  )
})

test_that("weight validation accepts named weights and reorders them", {
  weights <- c(
    S_charge = 0.15,
    S_hydro = 0.15,
    S_count = 0.20,
    S_coverage = 0.25,
    S_length = 0.25
  )

  expect_equal(
    pepVet:::.validate_weights(weights, has_proteome = FALSE),
    expected_protein_only_weights
  )

  preset_weights <- pepvet_preset("standard")$weights
  expect_equal(
    pepVet:::.validate_weights(preset_weights, has_proteome = FALSE),
    expected_protein_only_weights
  )
})

test_that("weight validation rejects bad sums, negatives, and wrong lengths", {
  expect_error(
    pepVet:::.validate_weights(c(0.5, 0.5, 0.5, 0.5, 0.5), FALSE),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    pepVet:::.validate_weights(c(-0.1, 0.3, 0.3, 0.3, 0.2), FALSE),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    pepVet:::.validate_weights(rep(0.2, 6), FALSE),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    pepVet:::.validate_weights(rep(1 / 5, 5), TRUE),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    pepVet:::.validate_weights(pepvet_preset("targeted")$weights, FALSE),
    class = "pepvet_error_invalid_weights"
  )
  expect_error(
    pepVet:::.validate_weights(
      c(foo = 1, bar = 0, baz = 0, qux = 0, quux = 0),
      FALSE
    ),
    class = "pepvet_error_invalid_weights"
  )
})

test_that("score_peptides returns the documented schema without proteome", {
  digest_result <- digest_protein(
    reference_fasta("P02769.fasta"),
    enzyme = "trypsin"
  )
  result <- score_peptides(digest_result)

  expect_s3_class(result, "tbl_df")
  expect_identical(
    names(result),
    c(
      "protein_id",
      "S_length",
      "S_coverage",
      "S_count",
      "S_hydro",
      "S_charge",
      "composite_score",
      "verdict",
      "median_peptide_length",
      "preset_used"
    )
  )
  expect_type(result$protein_id, "character")
  expect_type(result$S_length, "double")
  expect_type(result$S_coverage, "double")
  expect_type(result$S_count, "double")
  expect_type(result$S_hydro, "double")
  expect_type(result$S_charge, "double")
  expect_type(result$composite_score, "double")
  expect_type(result$verdict, "character")
  expect_type(result$median_peptide_length, "double")
  expect_type(result$preset_used, "character")
  expect_false(anyNA(result))
  expect_identical(result$preset_used, "standard")
})

test_that("score_peptides includes S_unique in proteome-aware mode", {
  proteome_digest <- digest_protein(
    c(target = "AAAAAAARAAAAAAAK", background = "AAAAAAARGGGGGGGK")
  )
  target_digest <- proteome_digest[
    proteome_digest$protein_id == "target",
    ,
    drop = FALSE
  ]
  result <- score_peptides(target_digest, proteome = proteome_digest)

  expect_identical(
    names(result),
    c(
      "protein_id",
      "S_length",
      "S_coverage",
      "S_count",
      "S_hydro",
      "S_charge",
      "S_unique",
      "composite_score",
      "verdict",
      "median_peptide_length",
      "preset_used"
    )
  )
  expect_identical(result$S_unique, 0.5)
  expect_identical(result$preset_used, "custom")
})

test_that("score_peptides returns the digest-derived median peptide length", {
  digest_result <- make_digest_result(
    c(rep("AAAAAAAAAA", 5), strrep("A", 70)),
    starts = c(1L, 11L, 21L, 31L, 41L, 51L)
  )

  result <- score_peptides(digest_result, enzyme = "trypsin")

  expect_identical(result$median_peptide_length, 10)
})

test_that("preset helper returns modifiable scoring configuration lists", {
  preset <- pepvet_preset("standard")

  expect_identical(
    names(preset),
    c("gravy_range", "length_range", "weights", "include_pI")
  )
  expect_equal(preset$gravy_range, c(-1.0, 0.6))
  expect_equal(preset$length_range, c(7L, 25L))
  expect_true(all(c("S_length", "S_unique") %in% names(preset$weights)))
  expect_false(preset$include_pI)
})

test_that("score_peptides respects configurable length and GRAVY ranges", {
  digest_result <- make_digest_result(
    c("AAAAAA", "AAAAAAAA", "AIVVWWWK"),
    starts = c(1L, 7L, 14L)
  )

  default_length_result <- score_peptides(digest_result)
  widened_length_result <- score_peptides(
    digest_result,
    length_range = c(6L, 25L)
  )
  widened_hydro_result <- score_peptides(
    digest_result,
    gravy_range = c(-1.0, 1.5)
  )

  expect_lt(default_length_result$S_length, widened_length_result$S_length)
  expect_lt(default_length_result$S_hydro, widened_hydro_result$S_hydro)
})

test_that("presets can be applied directly to evaluate_digest", {
  bsa_path <- reference_fasta("P02769.fasta")
  standard_result <- do.call(
    evaluate_digest,
    c(list(sequence = bsa_path, enzyme = "trypsin"), pepvet_preset("standard"))
  )

  expect_s3_class(standard_result$scores, "tbl_df")
  expect_identical(standard_result$scores$median_peptide_length, 7)
  expect_identical(standard_result$scores$preset_used, "standard")
  expect_false("pI" %in% names(standard_result$scores))
})

test_that("fractionated preset enables peptide pI annotation", {
  bsa_path <- reference_fasta("P02769.fasta")
  result <- do.call(
    evaluate_digest,
    c(list(sequence = bsa_path, enzyme = "trypsin"), pepvet_preset("fractionated"))
  )

  expect_true("pI" %in% names(result$scores))
  expect_type(result$scores$pI, "list")
  expect_true(length(result$scores$pI[[1]]) > 0L)
  expect_type(unname(result$scores$pI[[1]]), "double")
  expect_identical(result$scores$preset_used, "fractionated")
})

test_that("score_peptides can append peptide pI values without affecting scores", {
  digest_result <- make_digest_result(c("PEPTIDEK", "AAAAAAAR"))
  baseline <- score_peptides(digest_result)
  with_pi <- score_peptides(digest_result, include_pI = TRUE)

  expect_equal(with_pi$composite_score, baseline$composite_score)
  expect_true("pI" %in% names(with_pi))
  expect_type(with_pi$pI, "list")
  expect_identical(names(with_pi$pI[[1]]), c("PEPTIDEK", "AAAAAAAR"))
  expect_identical(with_pi$preset_used, "fractionated")
})

test_that("preset_used falls back to custom when weights do not match a named preset", {
  digest_result <- make_digest_result(c("PEPTIDEK", "AAAAAAAR"))
  custom_result <- score_peptides(
    digest_result,
    weights = c(S_length = 0.30, S_coverage = 0.20, S_count = 0.20, S_hydro = 0.15, S_charge = 0.15)
  )

  expect_identical(custom_result$preset_used, "custom")
})

test_that("score_peptides falls back to enzyme-class expected lengths", {
  digest_result <- make_digest_result(
    c(strrep("A", 20), strrep("A", 24)),
    starts = c(1L, 21L)
  )

  lysc_result <- score_peptides(digest_result, enzyme = "lysc")
  trypsin_result <- score_peptides(digest_result, enzyme = "trypsin")

  expect_identical(lysc_result$median_peptide_length, 24)
  expect_identical(trypsin_result$median_peptide_length, 12)
})

test_that("score_peptides warns and zeros S_count for no-cleavage digests", {
  digest_result <- digest_protein(strrep("A", 20L), enzyme = "trypsin")

  expect_warning(
    result <- score_peptides(digest_result, enzyme = "trypsin"),
    "has no cleavage sites"
  )

  expect_identical(result$S_count, 0)
  expect_identical(result$median_peptide_length, 12)
})

test_that("score_peptides rejects invalid digest and proteome inputs", {
  valid_digest <- make_digest_result("AAAAAAAAAAAK")

  expect_error(
    score_peptides(tibble::tibble(foo = 1)),
    class = "pepvet_error_invalid_digest"
  )
  expect_error(
    score_peptides(valid_digest, proteome = tibble::tibble(foo = 1)),
    class = "pepvet_error_invalid_digest"
  )
})

test_that("component scores stay bounded on the reference grid", {
  for (file_name in reference_files) {
    fasta_path <- reference_fasta(file_name)

    for (enzyme in reference_enzymes) {
      digest_result <- digest_protein(fasta_path, enzyme = enzyme)
      score_result <- score_peptides(digest_result, enzyme = enzyme)
      numeric_columns <- setdiff(
        names(score_result),
        c("protein_id", "verdict", "median_peptide_length", "preset_used")
      )
      numeric_scores <- unlist(score_result[numeric_columns])

      expect_true(
        all(numeric_scores >= 0),
        info = paste(file_name, enzyme, "lower bound")
      )
      expect_true(
        all(numeric_scores <= 1),
        info = paste(file_name, enzyme, "upper bound")
      )
      expect_false(anyNA(numeric_scores), info = paste(file_name, enzyme, "NA"))
      expect_false(
        any(is.nan(numeric_scores)),
        info = paste(file_name, enzyme, "NaN")
      )
      expect_false(
        any(is.infinite(numeric_scores)),
        info = paste(file_name, enzyme, "Inf")
      )
    }
  }
})

test_that("composite equals the weighted sum and is deterministic", {
  digest_result <- digest_protein(
    reference_fasta("P02769.fasta"),
    enzyme = "trypsin"
  )
  result_one <- score_peptides(digest_result)
  result_two <- score_peptides(digest_result)
  expect_identical(result_one, result_two)
  expect_equal(
    result_one$composite_score,
    with(
      result_one,
      S_length * expected_protein_only_weights[["S_length"]] +
        S_coverage * expected_protein_only_weights[["S_coverage"]] +
        S_count * expected_protein_only_weights[["S_count"]] +
        S_hydro * expected_protein_only_weights[["S_hydro"]] +
        S_charge * expected_protein_only_weights[["S_charge"]]
    ),
    tolerance = 1e-10
  )
})

test_that("weight isolation maps the composite score to each component", {
  digest_result <- make_digest_result(
    c("AAAAAAAA", "AAAAAAAA", "AAAAAAAA", "AAA", "AAA"),
    starts = c(1L, 9L, 17L, 25L, 28L)
  )
  baseline <- score_peptides(digest_result)

  expect_equal(
    score_peptides(digest_result, weights = c(1, 0, 0, 0, 0))$composite_score,
    baseline$S_length
  )
  expect_equal(
    score_peptides(digest_result, weights = c(0, 1, 0, 0, 0))$composite_score,
    baseline$S_coverage
  )
  expect_equal(
    score_peptides(digest_result, weights = c(0, 0, 1, 0, 0))$composite_score,
    baseline$S_count
  )
  expect_equal(
    score_peptides(digest_result, weights = c(0, 0, 0, 1, 0))$composite_score,
    baseline$S_hydro
  )
  expect_equal(
    score_peptides(digest_result, weights = c(0, 0, 0, 0, 1))$composite_score,
    baseline$S_charge
  )
})

test_that("length scoring handles valid, invalid, ratio, and boundaries", {
  all_valid <- make_digest_result("AAAAAAAAAAAK")
  all_invalid <- make_digest_result(c("AAR", "AAR", "AA"))
  exact_ratio <- make_digest_result(
    c("AAAAAAAK", "AAAAAAAR", "AAAAAAAE", "AAK", "AAR")
  )
  boundary <- make_digest_result(
    c(
      "AAAAAAA",
      "AAAAAA",
      "AAAAAAAAAAAAAAAAAAAAAAAAA",
      "AAAAAAAAAAAAAAAAAAAAAAAAAA"
    )
  )

  expect_identical(pepVet:::.score_length(all_valid), 1)
  expect_identical(pepVet:::.score_length(all_invalid), 0)
  expect_identical(pepVet:::.score_length(exact_ratio), 0.6)
  expect_equal(
    pepVet:::.valid_length_mask(boundary),
    c(TRUE, FALSE, TRUE, FALSE)
  )
})

test_that("coverage scoring handles full, zero, partial, and overlap", {
  full_coverage <- make_digest_result(
    c(rep("AAAAAAAAAA", 10)),
    starts = seq(1L, 91L, by = 10L)
  )
  zero_coverage <- make_digest_result(c("AAA", "AAA", "AAA"))
  partial_coverage <- make_digest_result(
    c(strrep("A", 20), strrep("A", 20), strrep("A", 20), strrep("A", 40)),
    starts = c(1L, 21L, 41L, 61L)
  )
  overlapping <- make_digest_result(
    c("AAAAAAAAAA", "AAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA"),
    starts = c(1L, 11L, 1L),
    missed_cleavages = c(0L, 0L, 1L)
  )

  expect_identical(pepVet:::.score_coverage(full_coverage), 1)
  expect_identical(pepVet:::.score_coverage(zero_coverage), 0)
  expect_identical(pepVet:::.score_coverage(partial_coverage), 0.6)
  expect_identical(pepVet:::.score_coverage(overlapping), 1)
})

test_that("count scoring handles caps, ratios, and edge protein sizes", {
  capped <- make_digest_result(
    rep("AAAAAAAA", 15),
    starts = seq(1L, 113L, by = 8L)
  )
  exact_ratio <- make_digest_result(
    c(rep("AAAAAAAAAA", 5), strrep("A", 70)),
    starts = c(1L, 11L, 21L, 31L, 41L, 51L)
  )
  zero_valid <- make_digest_result(c("AAA", "AAA", "AAA"))
  short_protein <- make_digest_result("AAAAAAA")

  expect_identical(pepVet:::.score_count(capped), 1)
  expect_equal(pepVet:::.score_count(exact_ratio), 5 / 12)
  expect_identical(pepVet:::.score_count(zero_valid), 0)
  expect_warning(
    expect_identical(pepVet:::.score_count(short_protein), 0),
    "has no cleavage sites"
  )
})

test_that("hydrophobicity scoring respects thresholds and zero-valid cases", {
  in_range <- make_digest_result(c("ALIVDEK", "ALIVDEK"))
  too_hydrophobic <- make_digest_result("IIIIIIIIK")
  too_hydrophilic <- make_digest_result("RRRRRRRRK")
  boundary_low <- make_digest_result("SWWWWYY")
  boundary_high <- make_digest_result("STVVWWW")
  zero_valid <- make_digest_result(c("AAA", "AAA"))

  expect_equal(pepVet:::.calculate_gravy("SWWWWYY"), -1)
  expect_equal(pepVet:::.calculate_gravy("STVVWWW"), 0.6)
  expect_identical(pepVet:::.score_hydro(in_range), 1)
  expect_identical(pepVet:::.score_hydro(too_hydrophobic), 0)
  expect_identical(pepVet:::.score_hydro(too_hydrophilic), 0)
  expect_identical(pepVet:::.score_hydro(boundary_low), 1)
  expect_identical(pepVet:::.score_hydro(boundary_high), 1)
  expect_identical(pepVet:::.score_hydro(zero_valid), 0)
})

test_that("charge scoring checks internal basic residues only", {
  internal_basic <- make_digest_result("AAHAAAAAAK")
  no_internal_basic <- make_digest_result("AAAAAAAAAK")
  non_tryptic_pass <- make_digest_result("AKAAAAAAAAE")
  non_tryptic_fail <- make_digest_result("AAAAAAAAAAE")
  zero_valid <- make_digest_result(c("AAA", "AAA"))
  all_basic <- make_digest_result(c("AAHAAAAAAK", "AARAAAAAAK"))

  expect_identical(pepVet:::.score_charge(internal_basic), 1)
  expect_identical(pepVet:::.score_charge(no_internal_basic), 0)
  expect_identical(pepVet:::.score_charge(non_tryptic_pass), 1)
  expect_identical(pepVet:::.score_charge(non_tryptic_fail), 0)
  expect_identical(pepVet:::.score_charge(zero_valid), 0)
  expect_identical(pepVet:::.score_charge(all_basic), 1)
})

test_that("uniqueness scoring handles shared and unique peptides", {
  target <- make_digest_result(
    c("AAAAAAAAR", "LLLLLLLLLK"),
    protein_id = "target"
  )
  other <- make_digest_result(
    c("AAAAAAAAR", "GGGGGGGGGR"),
    protein_id = "other"
  )
  duplicate <- make_digest_result(
    c("AAAAAAAAR", "LLLLLLLLLK"),
    protein_id = "duplicate"
  )
  target_index <- pepVet:::.build_proteome_index(target)
  mixed_index <- pepVet:::.build_proteome_index(
    combine_digest_results(target, other)
  )
  shared_index <- pepVet:::.build_proteome_index(
    combine_digest_results(target, duplicate)
  )

  expect_identical(pepVet:::.score_unique(target, mixed_index), 0.5)
  expect_identical(pepVet:::.score_unique(target, target_index), 1)
  expect_identical(pepVet:::.score_unique(target, shared_index), 0)
  expect_false("S_unique" %in% names(score_peptides(target)))
})

test_that("verdict classification respects both decision boundaries", {
  expect_identical(
    pepVet:::.classify_verdict(c(0.7, 0.6999, 0.4, 0.3999, 0, 1)),
    c("Good", "Moderate", "Moderate", "Poor", "Poor", "Good")
  )
})

test_that("score_peptides reproduces boundary verdicts via weight isolation", {
  good_boundary <- make_digest_result(c(rep("AAAAAAAA", 7), rep("AAA", 3)))
  moderate_boundary <- make_digest_result(c(rep("AAAAAAAA", 2), rep("AAA", 3)))

  expect_identical(
    score_peptides(good_boundary, weights = c(1, 0, 0, 0, 0))$verdict,
    "Good"
  )
  expect_identical(
    score_peptides(moderate_boundary, weights = c(1, 0, 0, 0, 0))$verdict,
    "Moderate"
  )
})

test_that("fixture-backed scoring captures expected biological separation", {
  bsa <- score_peptides(
    digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  )
  bace1 <- score_peptides(
    digest_protein(reference_fasta("P56817.fasta"), enzyme = "trypsin")
  )
  histone_h3 <- score_peptides(
    digest_protein(reference_fasta("P68431.fasta"), enzyme = "trypsin")
  )

  expect_gt(bsa$S_length, histone_h3$S_length)
  expect_gt(bsa$composite_score, histone_h3$composite_score)
  expect_lt(bace1$composite_score, bsa$composite_score)
})
