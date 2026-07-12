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
  pepVet:::.bind_rows(list(...))
}

expected_protein_only_weights <- c(
  S_length   = 0.200,
  S_coverage = 0.348,
  S_count    = 0.226,
  S_hydro    = 0.138,
  S_charge   = 0.088
)

expected_proteome_weights <- c(
  S_length   = 0.160,
  S_coverage = 0.279,
  S_count    = 0.181,
  S_hydro    = 0.110,
  S_charge   = 0.070,
  S_unique   = 0.200
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
    S_charge   = 0.088,
    S_hydro    = 0.138,
    S_count    = 0.226,
    S_coverage = 0.348,
    S_length   = 0.200
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
    pepVet:::.validate_weights(rep(0, 5L), FALSE),
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
  expect_identical(
    pepVet:::.validate_weights(
      c(expected_protein_only_weights, S_unique = 0),
      has_proteome = FALSE
    ),
    expected_protein_only_weights
  )
  expect_identical(
    pepVet:::.validate_weights(
      unname(c(expected_protein_only_weights, S_unique = 0)),
      has_proteome = FALSE
    ),
    expected_protein_only_weights
  )
  expect_identical(
    pepVet:::.validate_weights(
      expected_proteome_weights,
      has_proteome = TRUE
    ),
    expected_proteome_weights
  )
  expect_identical(
    pepVet:::.validate_weights(
      unname(expected_proteome_weights),
      has_proteome = TRUE
    ),
    expected_proteome_weights
  )
})

test_that("weight normalization rejects incomplete and ambiguous names", {
  empty_name <- rep(0.2, 5L)
  names(empty_name) <- c(
    "S_length", "S_coverage", "S_count", "S_hydro", ""
  )
  invalid_weights <- list(
    wrong_length = rep(0.2, 4L),
    missing_name = c(
      S_length = 0.2, S_coverage = 0.2, S_count = 0.2, S_hydro = 0.2
    ),
    duplicate_name = c(
      S_length = 0.1, S_length = 0.1, S_coverage = 0.2,
      S_count = 0.3, S_hydro = 0.3
    ),
    empty_name = empty_name,
    unknown_name = c(
      foo = 0.2, S_coverage = 0.2, S_count = 0.2,
      S_hydro = 0.2, S_charge = 0.2
    ),
    missing_value = c(0.2, 0.2, 0.2, 0.2, NA_real_),
    non_finite = c(0.2, 0.2, 0.2, 0.2, Inf),
    wrong_type = rep("0.2", 5L)
  )

  for (input_name in names(invalid_weights)) {
    expect_error(
      pepVet:::.normalize_weights(
        invalid_weights[[input_name]],
        expected_protein_only_weights
      ),
      class = "pepvet_error_invalid_weights",
      info = input_name
    )
  }
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
    c(target = "AAAAAAARAAAAAAAK", background = "AAAAAAARGGGGGGGK"),
    missed_cleavages = 0L
  )
  target_digest <- proteome_digest[
    proteome_digest$protein_id == "target", ,
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

test_that("preset helper returns every documented scoring configuration", {
  expected_ranges <- list(
    standard = list(
      gravy = c(-1.0, 0.6), length = c(7L, 25L), pI = FALSE,
      weights = c(
        S_length = 0.200, S_coverage = 0.348, S_count = 0.226,
        S_hydro = 0.138, S_charge = 0.088, S_unique = 0.000
      )
    ),
    dia = list(
      gravy = c(-1.0, 0.8), length = c(7L, 30L), pI = FALSE,
      weights = c(
        S_length = 0.20, S_coverage = 0.30, S_count = 0.20,
        S_hydro = 0.10, S_charge = 0.10, S_unique = 0.10
      )
    ),
    targeted = list(
      gravy = c(-0.8, 0.4), length = c(8L, 20L), pI = FALSE,
      weights = c(
        S_length = 0.15, S_coverage = 0.10, S_count = 0.15,
        S_hydro = 0.15, S_charge = 0.15, S_unique = 0.30
      )
    ),
    membrane = list(
      gravy = c(-1.0, 2.0), length = c(7L, 30L), pI = FALSE,
      weights = c(
        S_length = 0.25, S_coverage = 0.25, S_count = 0.20,
        S_hydro = 0.05, S_charge = 0.15, S_unique = 0.10
      )
    ),
    ffpe_degraded = list(
      gravy = c(-1.0, 0.8), length = c(6L, 30L), pI = FALSE,
      weights = c(
        S_length = 0.20, S_coverage = 0.20, S_count = 0.30,
        S_hydro = 0.10, S_charge = 0.10, S_unique = 0.10
      )
    ),
    fractionated = list(
      gravy = c(-1.0, 0.6), length = c(7L, 25L), pI = TRUE,
      weights = c(
        S_length = 0.200, S_coverage = 0.348, S_count = 0.226,
        S_hydro = 0.138, S_charge = 0.088, S_unique = 0.000
      )
    )
  )

  for (preset_name in names(expected_ranges)) {
    preset <- pepvet_preset(preset_name)
    expected <- expected_ranges[[preset_name]]

    expect_identical(
      names(preset),
      c("gravy_range", "length_range", "weights", "include_pI"),
      info = preset_name
    )
    expect_equal(preset$gravy_range, expected$gravy, info = preset_name)
    expect_identical(preset$length_range, expected$length, info = preset_name)
    expect_identical(preset$include_pI, expected$pI, info = preset_name)
    expect_identical(
      names(preset$weights),
      names(expected_proteome_weights),
      info = preset_name
    )
    expect_equal(preset$weights, expected$weights, tolerance = 1e-12)
    expect_true(all(is.finite(preset$weights)), info = preset_name)
    expect_true(all(preset$weights >= 0), info = preset_name)
    expect_equal(sum(preset$weights), 1, tolerance = 1e-12)
  }

  expect_identical(
    pepvet_preset(" DIA ")$length_range,
    c(7L, 30L)
  )

  standard <- pepvet_preset("standard")
  standard$weights[["S_length"]] <- 0
  expect_equal(
    pepvet_preset("standard")$weights[["S_length"]],
    0.2
  )
})

test_that("preset identity distinguishes named and custom configurations", {
  for (preset_name in names(pepVet:::.pepvet_presets)) {
    preset <- pepvet_preset(preset_name)
    expect_identical(
      pepVet:::.identify_preset_used(
        gravy_range = preset$gravy_range,
        length_range = preset$length_range,
        weights = preset$weights,
        include_pI = preset$include_pI,
        has_proteome = TRUE
      ),
      preset_name,
      info = preset_name
    )
  }

  standard <- pepvet_preset("standard")
  protein_weights <- pepVet:::.validate_weights(
    standard$weights,
    has_proteome = FALSE
  )
  expect_identical(
    pepVet:::.identify_preset_used(
      standard$gravy_range,
      standard$length_range,
      protein_weights,
      standard$include_pI,
      has_proteome = FALSE
    ),
    "standard"
  )
  expect_identical(
    pepVet:::.identify_preset_used(
      c(-0.9, 0.6),
      standard$length_range,
      protein_weights,
      standard$include_pI,
      has_proteome = FALSE
    ),
    "custom"
  )

  expect_true(
    pepVet:::.same_numeric_values(c(1, 2), c(1 + 1e-9, 2))
  )
  expect_false(
    pepVet:::.same_numeric_values(c(1, 2), c(1 + 1e-5, 2))
  )
  expect_true(
    pepVet:::.same_named_weights(c(a = 1), c(a = 1 + 1e-9))
  )
  expect_false(pepVet:::.same_named_weights(c(a = 1), c(b = 1)))
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
    c(list(
      sequence = bsa_path, enzyme = "trypsin",
      missed_cleavages = 0L
    ), pepvet_preset("standard"))
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
  digest_result <- make_digest_result(
    c("PEPTIDEK", "AAAAAAAR", "AAA"),
    starts = c(1L, 9L, 17L)
  )
  baseline <- score_peptides(digest_result)
  with_pi <- score_peptides(digest_result, include_pI = TRUE)

  expect_equal(with_pi$composite_score, baseline$composite_score)
  expect_true("pI" %in% names(with_pi))
  expect_type(with_pi$pI, "list")
  expect_identical(names(with_pi$pI[[1]]), c("PEPTIDEK", "AAAAAAAR"))
  expect_identical(with_pi$preset_used, "fractionated")

  no_valid_peptides <- score_peptides(
    make_digest_result(
      c("AAA", "CCC"),
      starts = c(1L, 4L)
    ),
    include_pI = TRUE
  )
  expect_identical(no_valid_peptides$pI[[1]], numeric(0))
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

test_that("fallback peptide lengths cover every supported enzyme group", {
  expected_fallbacks <- c(
    trypsin = 12,
    `trypsin-high` = 12,
    `trypsin-low` = 12,
    `trypsin-simple` = 12,
    lysc = 24,
    `arg-c proteinase` = 24,
    `glutamyl endopeptidase` = 17,
    `asp-n endopeptidase` = 17,
    `chymotrypsin-high` = 11,
    `chymotrypsin-low` = 11,
    unknown = 12
  )

  for (enzyme_name in names(expected_fallbacks)) {
    lookup_name <- if (identical(enzyme_name, "unknown")) {
      "unsupported-but-normalized"
    } else {
      enzyme_name
    }

    expect_identical(
      pepVet:::.fallback_expected_peptide_length(lookup_name),
      as.numeric(expected_fallbacks[[enzyme_name]]),
      info = enzyme_name
    )
  }

  three_peptides <- make_digest_result(
    c(strrep("A", 8), strrep("A", 12), strrep("A", 16)),
    starts = c(1L, 9L, 21L)
  )
  two_peptides <- three_peptides[1:2, , drop = FALSE]

  expect_identical(
    pepVet:::.expected_peptide_length(three_peptides, enzyme = "lysc"),
    12
  )
  expect_identical(
    pepVet:::.expected_peptide_length(two_peptides, enzyme = "lysc"),
    24
  )
})

test_that("score_peptides warns and zeros S_count for no-cleavage digests", {
  digest_result <- digest_protein(strrep("A", 20L), enzyme = "trypsin")

  expect_warning(
    result <- score_peptides(digest_result, enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )

  expect_identical(result$S_count, 0)
  expect_identical(result$median_peptide_length, 12)
})

test_that("zero-cleavage scoring hard-fails across component partitions", {
  no_valid_peptide <- digest_protein("A", enzyme = "trypsin")
  full_length_peptide <- digest_protein(
    strrep("A", 20L),
    enzyme = "trypsin"
  )

  expect_warning(
    poor_result <- score_peptides(no_valid_peptide, enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  poor_components <- c(
    S_length = 0,
    S_coverage = 0,
    S_count = 0,
    S_hydro = 0,
    S_charge = 0
  )
  expect_equal(
    unlist(poor_result[names(poor_components)], use.names = FALSE),
    unname(poor_components)
  )
  expect_identical(poor_result$composite_score, 0)
  expect_identical(poor_result$verdict, "Poor")

  expect_warning(
    full_result <- score_peptides(full_length_peptide, enzyme = "trypsin"),
    class = "pepvet_warning_no_cleavage_sites"
  )
  expected_components <- c(
    S_length = 1,
    S_coverage = 1,
    S_count = 0,
    S_hydro = 0,
    S_charge = 0
  )
  weighted_sum <- sum(
    expected_components * expected_protein_only_weights[names(expected_components)]
  )
  actual_components <- unlist(
    full_result[names(expected_components)],
    use.names = FALSE
  )

  expect_equal(
    unname(actual_components),
    unname(expected_components),
    tolerance = 1e-12
  )
  expect_gt(weighted_sum, 0)
  expect_identical(full_result$composite_score, 0)
  expect_identical(full_result$verdict, "Poor")
})

test_that("digest validation rejects malformed tables with classed errors", {
  valid_digest <- make_digest_result(
    c("AAAAAAA", "CCCCCCC"),
    starts = c(1L, 8L)
  )

  with_extra_column <- valid_digest
  with_extra_column$extra <- TRUE
  normalized <- pepVet:::.validate_digest_result(with_extra_column)

  expect_s3_class(normalized, "tbl_df")
  expect_identical(names(normalized), pepVet:::.required_digest_columns)
  expect_type(normalized$start, "integer")
  expect_type(normalized$end, "integer")
  expect_type(normalized$length, "integer")
  expect_type(normalized$missed_cleavages, "integer")

  empty_digest <- tibble::tibble(
    protein_id = character(),
    peptide = character(),
    start = integer(),
    end = integer(),
    length = integer(),
    missed_cleavages = integer()
  )
  wrong_protein_type <- valid_digest
  wrong_protein_type$protein_id <- seq_len(nrow(valid_digest))
  wrong_peptide_type <- valid_digest
  wrong_peptide_type$peptide <- seq_len(nrow(valid_digest))
  wrong_coordinate_type <- valid_digest
  wrong_coordinate_type$start <- as.character(valid_digest$start)
  missing_value <- valid_digest
  missing_value$peptide[[1L]] <- NA_character_
  infinite_coordinate <- valid_digest
  infinite_coordinate$end[[1L]] <- Inf
  oversized_coordinate <- valid_digest
  oversized_coordinate$end[[1L]] <- .Machine$integer.max + 1
  fractional_coordinate <- valid_digest
  fractional_coordinate$start[[1L]] <- 1.5
  fractional_length <- valid_digest
  fractional_length$length[[1L]] <- 7.5
  fractional_missed <- valid_digest
  fractional_missed$missed_cleavages[[1L]] <- 0.5
  invalid_start <- valid_digest
  invalid_start$start[[1L]] <- 0L
  invalid_order <- valid_digest
  invalid_order$end[[1L]] <- invalid_order$start[[1L]] - 1L
  inconsistent_length <- valid_digest
  inconsistent_length$length[[1L]] <- inconsistent_length$length[[1L]] + 1L
  negative_missed <- valid_digest
  negative_missed$missed_cleavages[[1L]] <- -1L

  invalid_cases <- list(
    wrong_object_type = 42,
    empty = empty_digest,
    wrong_protein_type = wrong_protein_type,
    wrong_peptide_type = wrong_peptide_type,
    wrong_coordinate_type = wrong_coordinate_type,
    missing_value = missing_value,
    infinite_coordinate = infinite_coordinate,
    oversized_coordinate = oversized_coordinate,
    fractional_coordinate = fractional_coordinate,
    fractional_length = fractional_length,
    fractional_missed = fractional_missed,
    invalid_start = invalid_start,
    invalid_order = invalid_order,
    inconsistent_length = inconsistent_length,
    negative_missed = negative_missed
  )

  for (case_name in names(invalid_cases)) {
    expect_error(
      pepVet:::.validate_digest_result(invalid_cases[[case_name]]),
      class = "pepvet_error_invalid_digest",
      info = case_name
    )
  }
})

test_that("score_peptides rejects invalid digest and proteome inputs", {
  valid_digest <- make_digest_result("AAAAAAAAAAAK")
  empty_digest <- tibble::tibble(
    protein_id = character(),
    peptide = character(),
    start = integer(),
    end = integer(),
    length = integer(),
    missed_cleavages = integer()
  )

  expect_error(
    score_peptides(tibble::tibble(foo = 1)),
    class = "pepvet_error_invalid_digest"
  )
  expect_error(
    score_peptides(42),
    class = "pepvet_error_invalid_digest"
  )
  expect_error(
    score_peptides(empty_digest),
    class = "pepvet_error_invalid_digest"
  )
  expect_error(
    score_peptides(valid_digest, proteome = tibble::tibble(foo = 1)),
    class = "pepvet_error_invalid_digest"
  )
  expect_error(
    score_peptides(valid_digest, enzyme = NULL),
    class = "pepvet_error_invalid_enzyme"
  )
  expect_error(
    score_peptides(valid_digest, weights = rep(0.2, 4L)),
    class = "pepvet_error_invalid_weights"
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

test_that("score_peptides returns one stable row per input protein", {
  digest_result <- make_digest_result(
    c(
      strrep("A", 7L), strrep("A", 7L),
      strrep("C", 7L), strrep("C", 7L)
    ),
    starts = c(1L, 8L, 1L, 8L),
    protein_id = c("first", "first", "second", "second")
  )
  result <- score_peptides(digest_result)

  expect_identical(result$protein_id, c("first", "second"))
  expect_identical(nrow(result), 2L)
  expect_identical(result$median_peptide_length, c(12, 12))
  expect_true(all(result$S_length == 1))
  expect_true(all(result$S_coverage == 1))
  expect_true(all(result$S_count == 1))
  expect_true(all(result$composite_score >= 0))
  expect_true(all(result$composite_score <= 1))
})

test_that("composite equals the weighted sum and is deterministic", {
  digest_result <- digest_protein(
    reference_fasta("P02769.fasta"),
    enzyme = "trypsin"
  )
  result_one <- score_peptides(digest_result)
  result_two <- score_peptides(digest_result)
  expect_equal(result_one, result_two, tolerance = 1e-15)
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
  expect_equal(pepVet:::.score_length(exact_ratio), 0.6, tolerance = 1e-10)
  expect_equal(
    pepVet:::.valid_length_mask(boundary),
    c(TRUE, FALSE, TRUE, FALSE)
  )
})

test_that("valid peptide extraction preserves order, columns, and empty type", {
  digest_result <- make_digest_result(
    c("AAAAAAA", "A", strrep("A", 25L), strrep("A", 26L)),
    starts = c(1L, 8L, 9L, 34L)
  )

  mask <- pepVet:::.valid_length_mask(digest_result)
  extracted <- pepVet:::.extract_valid_digest(digest_result)
  empty <- pepVet:::.extract_valid_digest(
    digest_result,
    length_range = c(100L, 101L)
  )

  expect_identical(mask, c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(
    extracted$peptide,
    c("AAAAAAA", strrep("A", 25L))
  )
  expect_identical(names(extracted), names(digest_result))
  expect_identical(nrow(empty), 0L)
  expect_identical(names(empty), names(digest_result))
  expect_type(empty$peptide, "character")
  expect_type(empty$start, "integer")
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
  expect_equal(pepVet:::.score_coverage(partial_coverage), 0.6, tolerance = 1e-10)
  expect_identical(pepVet:::.score_coverage(overlapping), 1)
})

test_that("coverage scoring reduces unsorted and adjacent intervals independently", {
  digest_result <- make_digest_result(
    rep(strrep("A", 10L), 5L),
    starts = c(21L, 1L, 11L, 6L, 31L)
  )
  interval_set <- IRanges::IRanges(
    start = digest_result$start,
    end = digest_result$end
  )
  reduced_intervals <- IRanges::reduce(interval_set)
  expected_coverage <- sum(IRanges::width(reduced_intervals)) /
    max(digest_result$end)

  expect_equal(
    pepVet:::.score_coverage(digest_result),
    expected_coverage,
    tolerance = 1e-12
  )

  adjacent <- make_digest_result(
    c(strrep("A", 5L), strrep("A", 5L)),
    starts = c(1L, 6L)
  )
  expect_identical(
    pepVet:::.score_coverage(adjacent, length_range = c(1L, 25L)),
    1
  )
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
    class = "pepvet_warning_no_cleavage_sites"
  )
})

test_that("hydrophobicity scoring respects thresholds and zero-valid cases", {
  in_range <- make_digest_result(c("ALIVDEK", "ALIVDEK"))
  too_hydrophobic <- make_digest_result("IIIIIIIIK")
  too_hydrophilic <- make_digest_result("RRRRRRRRK")
  boundary_low <- make_digest_result("SWWWWYY")
  boundary_high <- make_digest_result("STVVWWW")
  zero_valid <- make_digest_result(c("AAA", "AAA"))
  missing_hydrophobicity <- make_digest_result("O")
  mixed_hydrophobicity <- make_digest_result("AO")

  expect_equal(pepVet:::.calculate_gravy("SWWWWYY"), -1)
  expect_equal(pepVet:::.calculate_gravy("STVVWWW"), 0.6)
  expect_equal(pepVet:::.calculate_gravy("AO"), 1.8, tolerance = 1e-12)
  expect_identical(pepVet:::.score_hydro(in_range), 1)
  expect_identical(pepVet:::.score_hydro(too_hydrophobic), 0)
  expect_identical(pepVet:::.score_hydro(too_hydrophilic), 0)
  expect_identical(pepVet:::.score_hydro(boundary_low), 1)
  expect_identical(pepVet:::.score_hydro(boundary_high), 1)
  expect_identical(pepVet:::.score_hydro(zero_valid), 0)
  expect_identical(
    pepVet:::.score_hydro(
      missing_hydrophobicity,
      length_range = c(1L, 25L)
    ),
    0
  )
  expect_identical(
    pepVet:::.score_hydro(
      mixed_hydrophobicity,
      gravy_range = c(-2, 2),
      length_range = c(1L, 25L)
    ),
    1
  )

  public_result <- score_peptides(
    make_digest_result(c("O", "O")),
    length_range = c(1L, 25L)
  )
  expect_identical(public_result$S_hydro, 0)
  expect_true(is.finite(public_result$composite_score))
  expect_false(anyNA(public_result))
})

test_that("charge scoring checks internal basic residues only", {
  internal_basic <- make_digest_result("AAHAAAAAAK")
  no_internal_basic <- make_digest_result("AAAAAAAAAK")
  non_tryptic_pass <- make_digest_result("AKAAAAAAAAE")
  non_tryptic_fail <- make_digest_result("AAAAAAAAAAE")
  zero_valid <- make_digest_result(c("AAA", "AAA"))
  all_basic <- make_digest_result(c("AAHAAAAAAK", "AARAAAAAAK"))
  single_basic <- make_digest_result(c("K", "R", "H"))

  expect_identical(pepVet:::.score_charge(internal_basic), 1)
  expect_identical(pepVet:::.score_charge(no_internal_basic), 0)
  expect_identical(pepVet:::.score_charge(non_tryptic_pass), 1)
  expect_identical(pepVet:::.score_charge(non_tryptic_fail), 0)
  expect_identical(pepVet:::.score_charge(zero_valid), 0)
  expect_identical(pepVet:::.score_charge(all_basic), 1)
  expect_identical(
    pepVet:::.score_charge(single_basic, length_range = c(1L, 1L)),
    0
  )
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
  expect_identical(
    pepVet:::.score_unique(
      make_digest_result(c("AAA", "CCC")),
      target_index
    ),
    0
  )
  expect_false("S_unique" %in% names(score_peptides(target)))

  unrelated <- make_digest_result(
    "CCCCCCCCK",
    protein_id = "unrelated"
  )
  expect_identical(
    score_peptides(
      target,
      proteome = unrelated,
      length_range = c(1L, 25L)
    )$S_unique,
    1
  )

  repeated <- make_digest_result(
    c("AAAAAAAAR", "AAAAAAAAR"),
    starts = c(1L, 10L),
    protein_id = "target"
  )
  expect_identical(
    pepVet:::.score_unique(
      repeated,
      pepVet:::.build_proteome_index(repeated)
    ),
    1
  )
  expect_error(
    pepVet:::.score_unique(
      combine_digest_results(target, other),
      mixed_index
    ),
    class = "pepvet_error_invalid_digest"
  )
})

test_that("verdict classification respects both decision boundaries", {
  vg <- pepVet:::.get_param("verdict_good")
  vm <- pepVet:::.get_param("verdict_moderate")
  expect_identical(
    pepVet:::.classify_verdict(c(vg, vg - 1e-4, vm, vm - 1e-4, 0, 1)),
    c("Good", "Moderate", "Moderate", "Poor", "Poor", "Good")
  )
  expect_identical(
    pepVet:::.classify_verdict(c(-Inf, Inf, NA_real_)),
    c("Poor", "Good", NA_character_)
  )
})

test_that("random valid score inputs preserve component and composite bounds", {
  withr::local_seed(20260711)
  amino_acids <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1L]]
  component_names <- names(expected_protein_only_weights)

  for (iteration in seq_len(40L)) {
    n_peptides <- sample(3L:12L, 1L)
    peptide_lengths <- sample(3L:30L, n_peptides, replace = TRUE)
    peptides <- vapply(
      peptide_lengths,
      function(peptide_length) {
        paste(sample(amino_acids, peptide_length, replace = TRUE),
          collapse = ""
        )
      },
      character(1)
    )
    starts <- c(1L, head(cumsum(peptide_lengths) + 1L, -1L))
    digest_result <- make_digest_result(peptides, starts = starts)
    weights <- stats::runif(length(component_names))
    weights <- weights / sum(weights)
    names(weights) <- component_names
    result <- score_peptides(digest_result, weights = weights)
    components <- unname(unlist(result[component_names], use.names = FALSE))

    expect_true(all(components >= 0 & components <= 1),
      info = paste("component bounds", iteration)
    )
    expect_true(
      result$composite_score >= min(components) &&
        result$composite_score <= max(components),
      info = paste("convex composite", iteration)
    )
    expect_equal(
      result$composite_score,
      sum(components * weights),
      tolerance = 1e-12,
      info = paste("weighted sum", iteration)
    )
  }
})

test_that("score_peptides reproduces boundary verdicts via weight isolation", {
  # 7 valid out of 10 peptides: S_length = 0.7 >= 0.65 => Good
  good_boundary <- make_digest_result(c(rep("AAAAAAAA", 7), rep("AAA", 3)))
  # 2 valid out of 5 peptides: S_length = 0.4 >= 0.4 => Moderate
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

# Error class tests.

test_that("score_peptides rejects invalid gravy_range", {
  d <- digest_protein(.bsa_path, enzyme = "trypsin")
  expect_error(
    score_peptides(d, gravy_range = NULL),
    class = "pepvet_error_invalid_gravy_range"
  )
  expect_error(
    score_peptides(d, gravy_range = c(1, -1)),
    class = "pepvet_error_invalid_gravy_range"
  )
})

test_that("score_peptides rejects invalid length_range", {
  d <- digest_protein(.bsa_path, enzyme = "trypsin")
  expect_error(
    score_peptides(d, length_range = NULL),
    class = "pepvet_error_invalid_length_range"
  )
  expect_error(
    score_peptides(d, length_range = c(0, 0)),
    class = "pepvet_error_invalid_length_range"
  )
})

test_that("score_peptides rejects invalid include_pI", {
  d <- digest_protein(.bsa_path, enzyme = "trypsin")
  invalid_values <- list(
    null = NULL,
    missing = NA,
    wrong_type = 1,
    multiple = c(FALSE, TRUE),
    character = "no"
  )

  for (input_name in names(invalid_values)) {
    expect_error(
      score_peptides(d, include_pI = invalid_values[[input_name]]),
      class = "pepvet_error_invalid_include_pi",
      info = input_name
    )
  }
})

test_that("pepvet_preset rejects invalid preset name", {
  invalid_values <- list(
    null = NULL,
    missing = NA_character_,
    empty = "",
    multiple = c("standard", "dia"),
    wrong_type = 42,
    unsupported = "nonexistent_preset"
  )

  for (input_name in names(invalid_values)) {
    expect_error(
      pepvet_preset(invalid_values[[input_name]]),
      class = "pepvet_error_invalid_preset",
      info = input_name
    )
  }
})
