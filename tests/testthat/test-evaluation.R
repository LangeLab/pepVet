# Cross-function consistency.

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

test_that("named input order changes presentation but not named results", {
  sequences <- Biostrings::AAStringSet(c(
    zeta = "AAAAAAAKAAAAAAARAAAAAAAK",
    alpha = "CCCCCCCKCCCCCCCRCCCCCCCK",
    mu = "PEPTIDEKAPEPTIDERAPEPTIDEK"
  ))
  reordered_sequences <- sequences[rev(seq_along(sequences))]

  forward_evaluation <- evaluate_digest(
    sequences,
    enzyme = "trypsin",
    missed_cleavages = 0L
  )
  reverse_evaluation <- evaluate_digest(
    reordered_sequences,
    enzyme = "trypsin",
    missed_cleavages = 0L
  )

  expect_identical(
    forward_evaluation$params$protein_ids,
    names(sequences)
  )
  expect_identical(
    reverse_evaluation$params$protein_ids,
    names(reordered_sequences)
  )
  expect_identical(
    forward_evaluation$scores$protein_id,
    names(sequences)
  )
  expect_identical(
    reverse_evaluation$scores$protein_id,
    names(reordered_sequences)
  )
  expect_identical(
    unique(forward_evaluation$peptides$protein_id),
    names(sequences)
  )
  expect_identical(
    unique(reverse_evaluation$peptides$protein_id),
    names(reordered_sequences)
  )

  aligned_score_index <- match(
    forward_evaluation$scores$protein_id,
    reverse_evaluation$scores$protein_id
  )
  aligned_reverse_scores <- reverse_evaluation$scores[
    aligned_score_index, , drop = FALSE
  ]
  expect_equal(
    forward_evaluation$scores,
    aligned_reverse_scores,
    tolerance = 1e-15
  )

  peptide_columns <- setdiff(
    names(forward_evaluation$peptides),
    "protein_id"
  )
  for (protein_id in names(sequences)) {
    forward_peptides <- forward_evaluation$peptides[
      forward_evaluation$peptides$protein_id == protein_id,
      peptide_columns,
      drop = FALSE
    ]
    reverse_peptides <- reverse_evaluation$peptides[
      reverse_evaluation$peptides$protein_id == protein_id,
      peptide_columns,
      drop = FALSE
    ]
    expect_identical(
      forward_peptides,
      reverse_peptides,
      info = protein_id
    )
  }

  forward_batch <- batch_evaluate(
    sequences,
    enzyme = "trypsin",
    missed_cleavages = 0L,
    cores = 1L
  )
  reverse_batch <- batch_evaluate(
    reordered_sequences,
    enzyme = "trypsin",
    missed_cleavages = 0L,
    cores = 1L
  )
  expect_identical(forward_batch$protein_id, names(sequences))
  expect_identical(reverse_batch$protein_id, names(reordered_sequences))

  aligned_batch_index <- match(
    forward_batch$protein_id,
    reverse_batch$protein_id
  )
  aligned_reverse_batch <- reverse_batch[
    aligned_batch_index, , drop = FALSE
  ]
  expect_equal(
    as.data.frame(forward_batch),
    as.data.frame(aligned_reverse_batch),
    tolerance = 1e-15
  )
  expect_identical(
    attr(forward_batch, "scoring_config"),
    attr(reverse_batch, "scoring_config")
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
  expect_identical(
    attr(batch_with, "scoring_config")$proteome_aware,
    TRUE
  )
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

# Comparison and recommendation.

test_that("compare_digests output is sorted by composite_score descending", {
  result <- compare_digests(.bsa_path, enzymes = c("trypsin", "lysc"))

  expect_s3_class(result, "tbl_df")
  expect_true(
    all(diff(result$composite_score) <= 0),
    info = "composite_score must be non-increasing across rows"
  )
  expect_identical(
    attr(result, "scoring_config")$weights,
    c(
      S_length = 0.200,
      S_coverage = 0.348,
      S_count = 0.226,
      S_hydro = 0.138,
      S_charge = 0.088
    )
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

test_that("compare_digests forwards pI output to every enzyme result", {
  result <- compare_digests(
    .bsa_path,
    enzymes = c("trypsin", "lysc"),
    include_pI = TRUE
  )

  expect_true("pI" %in% names(result))
  expect_type(result$pI, "list")
  expect_equal(length(result$pI), 2L)
})

test_that("subsetted batch_compare_enzymes objects print without error", {
  expect_warning(
    result <- suppressMessages(
      batch_compare_enzymes(.small_path, enzymes = c("trypsin", "lysc"))
    ),
    class = "pepvet_warning_no_cleavage_sites"
  )
  subsetted <- result[result$enzyme == "trypsin", c("protein_id", "composite_score")]

  expect_no_error(print(subsetted))
})

test_that("batch_compare_enzymes records class, levels, metadata, and progress", {
  sequences <- Biostrings::AAStringSet(c(
    alpha = "AKAAAAAAK",
    beta = "AKRTPK"
  ))
  expect_message(
    comparison <- batch_compare_enzymes(
      sequences,
      enzymes = c("Trypsin", " lysc "),
      missed_cleavages = 0L,
      gravy_range = c(-2.0, 2.0)
    ),
    class = "pepvet_message_batch_scoring"
  )

  expect_s3_class(comparison, "pepvet_batch_comparison")
  expect_identical(levels(comparison$enzyme), c("trypsin", "lysc"))
  expect_identical(attr(comparison, "n_proteins"), 2L)
  expect_identical(attr(comparison, "n_enzymes"), 2L)
  expect_identical(
    attr(comparison, "scoring_config")$length_range,
    c(7L, 25L)
  )
  expect_identical(
    attr(comparison, "scoring_config")$gravy_range,
    c(-2.0, 2.0)
  )
  expect_identical(nrow(comparison), 4L)
})

test_that("print.pepvet_batch_comparison exposes a stable summary contract", {
  sequences <- Biostrings::AAStringSet(c(
    alpha = "AKAAAAAAK",
    beta = "AKRTPK"
  ))
  comparison <- suppressMessages(
    batch_compare_enzymes(
      sequences,
      enzymes = c("trypsin", "lysc"),
      missed_cleavages = 0L
    )
  )
  expect_message(
    print(comparison),
    "pepVet batch enzyme comparison"
  )
  printed <- capture.output(suppressMessages(print(comparison)))
  printed_text <- paste(printed, collapse = "\n")

  expect_match(printed_text, "trypsin")
  expect_match(printed_text, "lysc")
  expect_match(printed_text, "4 rows total")
  expect_identical(
    suppressMessages(invisible(print(comparison))),
    comparison
  )

  fallback <- comparison
  attr(fallback, "enzymes") <- NULL
  attr(fallback, "n_proteins") <- NA_integer_
  attr(fallback, "n_enzymes") <- NA_integer_
  expect_message(
    suppressWarnings(print(fallback)),
    "pepVet batch enzyme comparison"
  )
})

test_that("batch_compare_enzymes serial and Unix parallel results agree", {
  skip_on_os("windows")
  sequences <- Biostrings::AAStringSet(c(
    alpha = "AKAAAAAAK",
    beta = "AKRTPK",
    gamma = "MKWVTFISLLFLFSSAYSR"
  ))
  serial <- suppressMessages(batch_compare_enzymes(
    sequences,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 0L,
    cores = 1L
  ))
  parallel <- suppressMessages(batch_compare_enzymes(
    sequences,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 0L,
    cores = 2L
  ))

  expect_equal(parallel, serial)
})

test_that("compare_digests rejects multi-protein input", {
  expect_error(
    compare_digests(.multi_path, enzymes = c("trypsin", "lysc")),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    compare_digests(character(0), enzymes = c("trypsin", "lysc")),
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
  # A poly-alanine sequence has no trypsin or lysc cut sites. Both return an
  # identical single-peptide digest and receive the same composite score.
  poly_ala <- strrep("A", 20L)
  warning_log <- new.env(parent = emptyenv())
  warning_log$classes <- character()
  result <- withCallingHandlers(
    recommend_enzyme(poly_ala, enzymes = c("trypsin", "lysc")),
    warning = function(condition) {
      warning_log$classes <- c(
        warning_log$classes,
        class(condition)[[1L]]
      )
      invokeRestart("muffleWarning")
    }
  )

  expect_type(result, "character")
  expect_identical(result, c("lysc", "trypsin"))
  expect_equal(
    sum(warning_log$classes == "pepvet_warning_no_cleavage_sites"),
    2L
  )
})

# Return structure.

test_that("evaluate_digest returns named list with scores, peptides, params", {
  result <- .fix_bsa_trypsin

  expect_type(result, "list")
  expect_identical(names(result), c("scores", "peptides", "params"))
  expect_s3_class(result$scores, "tbl_df")
  expect_s3_class(result$peptides, "tbl_df")
  expect_type(result$params, "list")
  expect_identical(
    names(result$params),
    c(
      "enzyme", "missed_cleavages", "protein_ids", "preset_used",
      "gravy_range", "length_range", "weights", "proteome_aware", "include_pI"
    )
  )
  expect_identical(result$params$preset_used, "standard")
  expect_identical(result$params$gravy_range, c(-1.0, 0.6))
  expect_identical(result$params$length_range, c(7L, 25L))
  expect_identical(
    result$params$weights,
    c(
      S_length = 0.200,
      S_coverage = 0.348,
      S_count = 0.226,
      S_hydro = 0.138,
      S_charge = 0.088
    )
  )
  expect_false(result$params$proteome_aware)
  expect_false(result$params$include_pI)
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

test_that("evaluate_digest stores custom and proteome-aware scoring metadata", {
  custom_weights <- c(
    S_length = 0.30,
    S_coverage = 0.25,
    S_count = 0.20,
    S_hydro = 0.15,
    S_charge = 0.10
  )
  custom <- evaluate_digest(
    .bsa_path,
    enzyme = "trypsin",
    weights = custom_weights,
    gravy_range = c(-2.0, 1.0),
    length_range = c(6L, 30L),
    include_pI = TRUE
  )

  expect_identical(custom$params$gravy_range, c(-2.0, 1.0))
  expect_identical(custom$params$length_range, c(6L, 30L))
  expect_identical(custom$params$weights, custom_weights)
  expect_false(custom$params$proteome_aware)
  expect_true(custom$params$include_pI)
  expect_identical(custom$params$preset_used, "custom")

  proteome_digest <- digest_protein(.multi_path, enzyme = "trypsin")
  aware_weights <- c(
    S_length = 0.16,
    S_coverage = 0.279,
    S_count = 0.181,
    S_hydro = 0.11,
    S_charge = 0.07,
    S_unique = 0.20
  )
  aware <- evaluate_digest(
    .bsa_path,
    enzyme = "trypsin",
    proteome = proteome_digest,
    weights = aware_weights
  )

  expect_true(aware$params$proteome_aware)
  expect_identical(aware$params$weights, aware_weights)
  expect_true("S_unique" %in% names(aware$scores))
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

  multi <- evaluate_digest(Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK"
  )), missed_cleavages = 0L)
  expect_identical(multi$scores$protein_id, c("first", "second"))
  expect_identical(multi$params$protein_ids, c("first", "second"))
  expect_true(all(multi$peptides$protein_id %in% multi$params$protein_ids))
})

# Error handling.

test_that("invalid sequence in batch_evaluate propagates a classed error", {
  expect_error(
    batch_evaluate("MXBZ123"),
    class = "pepvet_error_invalid_sequence"
  )
})

test_that("evaluate_digest rejects invalid input and enzyme values", {
  expect_error(
    evaluate_digest(NULL),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    evaluate_digest(42),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    evaluate_digest(character(0)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    evaluate_digest(.bsa_path, enzyme = NULL),
    class = "pepvet_error_invalid_enzyme"
  )
})

test_that("empty AAStringSet in batch_evaluate raises a classed error", {
  expect_error(
    batch_evaluate(Biostrings::AAStringSet()),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    batch_evaluate(character(0)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    batch_evaluate(42),
    class = "pepvet_error_invalid_input"
  )
})

test_that("batch_evaluate preserves names and rejects duplicate identifiers", {
  sequences <- Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK"
  ))
  one <- batch_evaluate(sequences[1L], missed_cleavages = 0L)

  expect_identical(one$protein_id, "first")
  expect_identical(nrow(one), 1L)

  duplicated <- Biostrings::AAStringSet(c(
    duplicate = "AKAAAAAAK",
    duplicate = "AKRTPK"
  ))
  expect_error(
    batch_evaluate(duplicated),
    class = "pepvet_error_invalid_input"
  )
})

test_that("batch_evaluate agrees with independent per-protein evaluations", {
  sequences <- Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK",
    third = "MKWVTFISLLFLFSSAYSR"
  ))
  batch <- batch_evaluate(sequences, missed_cleavages = 0L)
  individual <- lapply(seq_along(sequences), function(index) {
    evaluate_digest(sequences[index], missed_cleavages = 0L)
  })

  expect_identical(batch$protein_id, names(sequences))
  expect_equal(
    batch$composite_score,
    vapply(individual, function(result) result$scores$composite_score, numeric(1))
  )
  score_cols <- c(
    "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
  )
  expected_scores <- vapply(
    individual,
    function(result) as.numeric(result$scores[1L, score_cols, drop = TRUE]),
    numeric(length(score_cols))
  )
  expect_equal(
    as.numeric(as.matrix(batch[, score_cols])),
    as.numeric(t(expected_scores))
  )
  expect_equal(
    batch$median_peptide_length,
    vapply(
      individual,
      function(result) result$scores$median_peptide_length,
      numeric(1)
    )
  )
  expect_identical(
    batch$verdict,
    vapply(individual, function(result) result$scores$verdict, character(1))
  )
  expect_false("cleavage_efficiency" %in% names(batch))
})

test_that("batch_evaluate uses active scoring ranges for count and flags", {
  length_default <- batch_evaluate(
    c(length_case = "AKAAAAAAK"),
    enzyme = "trypsin",
    missed_cleavages = 0L
  )
  length_custom <- batch_evaluate(
    c(length_case = "AKAAAAAAK"),
    enzyme = "trypsin",
    missed_cleavages = 0L,
    length_range = c(8L, 25L)
  )
  expect_identical(length_default$n_valid_peptides, 1L)
  expect_identical(length_custom$n_valid_peptides, 0L)
  expect_false(length_default$flag_no_valid_peptides)
  expect_true(length_custom$flag_no_valid_peptides)

  expect_warning(
    hydro_default <- batch_evaluate(
      c(hydro_case = "VVVVVVV"),
      enzyme = "trypsin",
      missed_cleavages = 0L
    ),
    class = "pepvet_warning_no_cleavage_sites"
  )
  expect_warning(
    hydro_custom <- batch_evaluate(
      c(hydro_case = "VVVVVVV"),
      enzyme = "trypsin",
      missed_cleavages = 0L,
      gravy_range = c(-1.0, 5.0)
    ),
    class = "pepvet_warning_no_cleavage_sites"
  )
  expect_true(hydro_default$flag_hydrophobic)
  expect_false(hydro_custom$flag_hydrophobic)
  expect_identical(
    attr(hydro_custom, "scoring_config")$gravy_range,
    c(-1.0, 5.0)
  )
})

test_that("batch_evaluate serial and Unix parallel results agree", {
  skip_on_os("windows")
  sequences <- Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK",
    third = "MKWVTFISLLFLFSSAYSR",
    fourth = "AAAAAAAKAAAAAAAK"
  ))
  serial <- batch_evaluate(sequences, missed_cleavages = 0L, cores = 1L)
  parallel <- batch_evaluate(sequences, missed_cleavages = 0L, cores = 2L)

  expect_equal(parallel, serial)
})

test_that("Windows parallel dispatch selects the socket backend", {
  calls <- new.env(parent = emptyenv())
  calls$cores <- NULL
  testthat::local_mocked_bindings(
    .batch_socket_map = function(index_list, worker, cores) {
      calls$cores <- cores
      lapply(index_list, worker)
    },
    .package = "pepVet"
  )

  result <- pepVet:::.batch_parallel_map(
    as.list(1:3),
    function(index) index * 2L,
    cores = 2L,
    os_type = "windows"
  )

  expect_identical(result, as.list(c(2L, 4L, 6L)))
  expect_identical(calls$cores, 2L)
})

test_that("Windows socket workers preserve chunk order and values", {
  skip_if_not(.Platform$OS.type == "windows")
  running_from_pkgload <- "pkgload" %in% loadedNamespaces() && isTRUE(
    get("is_dev_package", asNamespace("pkgload"))("pepVet")
  )
  skip_if(running_from_pkgload, "requires an installed package")

  package_path <- normalizePath(find.package("pepVet"), mustWork = TRUE)
  result <- pepVet:::.batch_socket_map(
    as.list(1:3),
    function(index) {
      list(
        value = index * 2L,
        package_path = normalizePath(find.package("pepVet"), mustWork = TRUE)
      )
    },
    cores = 2L
  )

  expect_identical(
    vapply(result, `[[`, integer(1), "value"),
    c(2L, 4L, 6L)
  )
  expect_identical(
    vapply(result, `[[`, character(1), "package_path"),
    rep(package_path, 3L)
  )
})

test_that("Windows socket batch paths agree with serial results", {
  skip_if_not(.Platform$OS.type == "windows")
  running_from_pkgload <- "pkgload" %in% loadedNamespaces() && isTRUE(
    get("is_dev_package", asNamespace("pkgload"))("pepVet")
  )
  skip_if(running_from_pkgload, "requires an installed package")
  sequences <- Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK",
    third = "MKWVTFISLLFLFSSAYSR",
    fourth = "AAAAAAAKAAAAAAAK"
  ))

  serial <- batch_evaluate(sequences, missed_cleavages = 0L, cores = 1L)
  socket <- batch_evaluate(sequences, missed_cleavages = 0L, cores = 2L)
  expect_equal(socket, serial)

  serial_comparison <- suppressMessages(batch_compare_enzymes(
    sequences,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 0L,
    cores = 1L
  ))
  socket_comparison <- suppressMessages(batch_compare_enzymes(
    sequences,
    enzymes = c("trypsin", "lysc"),
    missed_cleavages = 0L,
    cores = 2L
  ))
  expect_equal(socket_comparison, serial_comparison)
})

test_that("batch_evaluate retries failed parallel chunks with a classed warning", {
  skip_on_os("windows")
  testthat::local_mocked_bindings(
    .batch_parallel_map = function(index_list, worker, cores) {
      lapply(index_list, function(index) {
        structure(list(), class = "try-error")
      })
    },
    .package = "pepVet"
  )
  sequences <- Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK"
  ))

  expect_warning(
    result <- batch_evaluate(sequences, missed_cleavages = 0L, cores = 2L),
    class = "pepvet_warning_parallel_retry"
  )
  baseline <- batch_evaluate(sequences, missed_cleavages = 0L, cores = 1L)
  expect_equal(result, baseline)
  expect_identical(result$protein_id, names(sequences))
})

# summarize_batch.

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

test_that("summarize_batch includes uniqueness when available", {
  proteome <- digest_protein(.multi_path)
  batch <- batch_evaluate(.multi_path, proteome = proteome)
  summary <- summarize_batch(batch)

  expect_true("S_unique" %in% names(summary$component_summary))
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

test_that(".validate_batch_result checks schema types and values", {
  valid <- .fix_batch_bsa
  expect_identical(.validate_batch_result(valid), valid)
  valid_with_unique <- valid
  valid_with_unique$S_unique <- 0.5
  expect_identical(
    .validate_batch_result(valid_with_unique),
    valid_with_unique
  )

  missing_score <- valid[, setdiff(names(valid), "S_length"), drop = FALSE]
  wrong_score_type <- valid
  wrong_score_type$S_length <- as.character(wrong_score_type$S_length)
  wrong_unique_type <- valid_with_unique
  wrong_unique_type$S_unique <- as.character(wrong_unique_type$S_unique)
  wrong_id_type <- valid
  wrong_id_type$protein_id <- 1L
  wrong_verdict_type <- valid
  wrong_verdict_type$verdict <- factor(wrong_verdict_type$verdict)
  wrong_flag_type <- valid
  wrong_flag_type$flag_short_protein <-
    as.integer(wrong_flag_type$flag_short_protein)
  invalid_score <- valid
  invalid_score$S_length[[1L]] <- Inf
  invalid_verdict <- valid
  invalid_verdict$verdict[[1L]] <- "Unknown"
  inconsistent_verdict <- valid
  inconsistent_verdict$composite_score[[1L]] <- 0.9
  inconsistent_verdict$verdict[[1L]] <- "Poor"
  blank_id <- valid
  blank_id$protein_id[[1L]] <- ""
  duplicate_id <- valid[rep(1L, 2L), , drop = FALSE]
  duplicate_id$protein_id <- c("same", "same")
  duplicate_columns <- as.data.frame(valid)
  names(duplicate_columns)[[2L]] <- names(duplicate_columns)[[1L]]
  hard_fail <- valid
  hard_fail$S_count[[1L]] <- 0
  hard_fail$S_coverage[[1L]] <- 1
  hard_fail$composite_score[[1L]] <- 0.5
  hard_fail$verdict[[1L]] <- "Moderate"
  empty <- valid[0, , drop = FALSE]

  invalid_results <- list(
    missing_score = missing_score,
    wrong_score_type = wrong_score_type,
    wrong_unique_type = wrong_unique_type,
    wrong_id_type = wrong_id_type,
    wrong_verdict_type = wrong_verdict_type,
    wrong_flag_type = wrong_flag_type,
    invalid_score = invalid_score,
    invalid_verdict = invalid_verdict,
    inconsistent_verdict = inconsistent_verdict,
    blank_id = blank_id,
    duplicate_id = duplicate_id,
    duplicate_columns = duplicate_columns,
    hard_fail = hard_fail,
    empty = empty
  )
  for (input_name in names(invalid_results)) {
    expect_error(
      .validate_batch_result(invalid_results[[input_name]]),
      class = "pepvet_error_invalid_batch_result",
      info = input_name
    )
  }

  expect_error(
    summarize_batch(duplicate_columns),
    class = "pepvet_error_invalid_batch_result"
  )
  expect_error(
    triage_proteins(duplicate_columns),
    class = "pepvet_error_invalid_batch_result"
  )
})

test_that(".batch_difficulty_flags evaluates each flag independently", {
  short <- tibble::tibble(
    protein_id = rep("short", 5L),
    peptide = rep("ACDEFGHIK", 5L),
    start = seq(1L, 37L, by = 9L),
    end = seq(9L, 45L, by = 9L),
    length = rep(9L, 5L),
    missed_cleavages = rep(0L, 5L)
  )
  no_valid <- tibble::tibble(
    protein_id = "no_valid",
    peptide = "ACDE",
    start = 1L,
    end = 200L,
    length = 4L,
    missed_cleavages = 0L
  )
  low_complexity <- tibble::tibble(
    protein_id = rep("low_complexity", 6L),
    peptide = rep("DDDDDDDDDDDDDDDDDDDD", 6L),
    start = seq(1L, 101L, by = 20L),
    end = seq(20L, 120L, by = 20L),
    length = rep(20L, 6L),
    missed_cleavages = rep(0L, 6L)
  )
  hydrophobic <- tibble::tibble(
    protein_id = rep("hydrophobic", 15L),
    peptide = rep(
      c("VVVVVVVV", "LILILILI", "FWFWFWFW", "IVIVIVIV"),
      length.out = 15L
    ),
    start = seq(1L, 113L, by = 8L),
    end = seq(8L, 120L, by = 8L),
    length = rep(8L, 15L),
    missed_cleavages = rep(0L, 15L)
  )
  unknown_gravy <- tibble::tibble(
    protein_id = "unknown_gravy",
    peptide = "OOOOOOOO",
    start = 1L,
    end = 8L,
    length = 8L,
    missed_cleavages = 0L
  )
  peptides <- .bind_rows(list(
    short, no_valid, low_complexity, hydrophobic, unknown_gravy
  ))
  ids <- c(
    "short", "no_valid", "low_complexity", "hydrophobic", "unknown_gravy"
  )
  flags <- .batch_difficulty_flags(peptides, ids)

  expect_identical(
    flags$protein_length,
    c(45L, 200L, 120L, 120L, 8L)
  )
  expect_identical(
    flags$n_valid_peptides,
    c(5L, 0L, 6L, 15L, 1L)
  )
  expect_identical(
    flags$flag_short_protein,
    c(TRUE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_identical(
    flags$flag_no_valid_peptides,
    c(FALSE, TRUE, FALSE, FALSE, FALSE)
  )
  expect_identical(
    flags$flag_low_complexity,
    c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )
  expect_identical(
    flags$flag_hydrophobic,
    c(FALSE, FALSE, FALSE, TRUE, FALSE)
  )
})

test_that("summarize_batch handles one row and constructed switch candidates", {
  one <- summarize_batch(.fix_batch_bsa)
  expect_identical(one$score_distribution[["sd"]], 0)
  expect_equal(nrow(one$problem_proteins), 1L)

  constructed <- .fix_batch_bsa[rep(1L, 4L), , drop = FALSE]
  constructed$protein_id <- c(
    "moderate_hydrophobic", "poor_short", "poor_unflagged", "good"
  )
  constructed$verdict <- c("Moderate", "Poor", "Poor", "Good")
  constructed$composite_score <- c(0.5, 0.3, 0.3, 0.8)
  constructed$flag_hydrophobic <- c(TRUE, FALSE, FALSE, FALSE)
  constructed$flag_short_protein <- c(FALSE, TRUE, FALSE, FALSE)
  constructed$flag_low_complexity <- FALSE
  constructed$flag_no_valid_peptides <- FALSE

  summary <- summarize_batch(constructed)
  expect_identical(
    summary$enzyme_switch_candidates$protein_id,
    c("moderate_hydrophobic", "poor_short")
  )
})

# triage_proteins.

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

test_that("triage_proteins covers every action partition", {
  constructed <- .fix_batch_bsa[rep(1L, 5L), , drop = FALSE]
  constructed$protein_id <- c(
    "proceed", "alternative", "other", "skip", "moderate_high"
  )
  constructed$verdict <- c("Good", "Moderate", "Poor", "Poor", "Moderate")
  constructed$S_length <- c(0.8, 0.4, 0.7, 0.7, 0.8)
  constructed$composite_score <- c(0.8, 0.5, 0.3, 0.3, 0.5)
  constructed$flag_hydrophobic <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
  constructed$flag_short_protein <- rep(FALSE, 5L)
  constructed$flag_low_complexity <- c(FALSE, FALSE, FALSE, TRUE, FALSE)
  constructed$flag_no_valid_peptides <- rep(FALSE, 5L)

  triaged <- triage_proteins(constructed)
  expect_identical(
    triaged$action,
    c(
      "proceed", "consider_alternative", "try_other_enzyme", "skip",
      "consider_alternative"
    )
  )
})

# Weight sensitivity analysis.

test_that("sensitivity_analysis returns correct structure for evaluate_digest input", {
  res <- .fix_bsa_trypsin
  withr::local_seed(42)
  sens <- sensitivity_analysis(res, n_iter = 500L)
  expect_type(sens, "list")
  expect_named(sens, c("iterations", "convergence", "summary", "settings"))
  expect_s3_class(sens$iterations, "tbl_df")
  expect_true("composite_score" %in% names(sens$iterations))
  expect_true("verdict" %in% names(sens$iterations))
  expect_equal(nrow(sens$iterations), 500L)
  expect_true(sens$summary$verdict_pct["Good"] > 0.95)
  expect_length(sens$summary$composite_ci, 2L)
  expect_true(sens$summary$composite_ci[[1L]] <= sens$summary$composite_ci[[2L]])
  expect_equal(nrow(sens$convergence), 500L)
  expect_true(all(
    sens$convergence$cumulative_stability >= 0 &
      sens$convergence$cumulative_stability <= 1
  ))
  expect_identical(sens$settings$nu, 63)
  expect_identical(sens$settings$n_iter, 500L)
  expect_identical(sens$settings$reference_weights, res$params$weights)
  expect_identical(sens$settings$distribution, "Dirichlet")
  expect_identical(
    sens$settings$verdict_thresholds,
    c(moderate = 0.40, good = 0.65)
  )
})

test_that("single sensitivity matches independent fixed-draw arithmetic", {
  weights <- stats::setNames(rep(0.2, 5L), c(
    "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
  ))
  scores <- tibble::tibble(
    protein_id = "protein_1",
    S_length = 1,
    S_coverage = 0,
    S_count = 1,
    S_hydro = 0,
    S_charge = 0,
    composite_score = 0.4,
    verdict = "Moderate"
  )
  draws <- rbind(
    c(0.70, 0.10, 0.10, 0.05, 0.05),
    c(0.10, 0.10, 0.10, 0.35, 0.35),
    rep(0.20, 5L),
    c(0.30, 0.10, 0.30, 0.15, 0.15)
  )

  sensitivity <- .sensitivity_single(
    scores,
    params = list(weights = weights),
    nu = 63,
    n_iter = 4L,
    importance = FALSE,
    corner_cases = FALSE,
    weight_draws = draws
  )

  expect_equal(
    sensitivity$iterations$composite_score,
    c(0.8, 0.2, 0.4, 0.6)
  )
  expect_identical(
    sensitivity$iterations$verdict,
    c("Good", "Poor", "Moderate", "Moderate")
  )
  expect_equal(
    sensitivity$convergence$cumulative_stability,
    c(0, 0, 1 / 3, 1 / 2)
  )
  expect_equal(
    sensitivity$summary$verdict_pct,
    c(Good = 0.25, Moderate = 0.50, Poor = 0.25)
  )
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

test_that("single sensitivity importance handles a fixed zero weight", {
  weights <- c(
    S_length = 0.200,
    S_coverage = 0.348,
    S_count = 0.226,
    S_hydro = 0.226,
    S_charge = 0
  )
  result <- evaluate_digest(.bsa_path, weights = weights)
  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(
    result,
    n_iter = 100L,
    importance = TRUE
  )

  expect_true(all(sensitivity$iterations$S_charge == 0))
  expect_identical(sensitivity$settings$fixed_zero_weights, "S_charge")
  expect_identical(sensitivity$summary$weight_importance[["S_charge"]], 0)
  expect_true(all(is.finite(sensitivity$summary$weight_importance)))
})

test_that("sensitivity_analysis corner_cases returns table", {
  res <- .fix_bsa_trypsin
  withr::local_seed(42)
  sens <- sensitivity_analysis(res, n_iter = 500L, corner_cases = TRUE)
  expect_s3_class(sens$summary$corner_cases, "tbl_df")
  expect_true("composite_at_lo" %in% names(sens$summary$corner_cases))
  expect_true("composite_at_hi" %in% names(sens$summary$corner_cases))

  components <- c(
    S_length = res$scores$S_length[[1L]],
    S_coverage = res$scores$S_coverage[[1L]],
    S_count = res$scores$S_count[[1L]],
    S_hydro = res$scores$S_hydro[[1L]],
    S_charge = res$scores$S_charge[[1L]]
  )
  reference_weights <- c(
    S_length = 0.200,
    S_coverage = 0.348,
    S_count = 0.226,
    S_hydro = 0.138,
    S_charge = 0.088
  )
  expected_bounds <- t(vapply(reference_weights, function(weight) {
    stats::qbeta(
      c(0.025, 0.975),
      shape1 = 63 * weight,
      shape2 = 63 * (1 - weight)
    )
  }, numeric(2)))
  expect_equal(sens$summary$corner_cases$lo, unname(expected_bounds[, 1L]))
  expect_equal(sens$summary$corner_cases$hi, unname(expected_bounds[, 2L]))
  expected_lo <- vapply(seq_len(nrow(sens$summary$corner_cases)), function(i) {
    weights <- reference_weights
    bound <- sens$summary$corner_cases$lo[[i]]
    weights[-i] <- weights[-i] * (1 - bound) / sum(weights[-i])
    weights[[i]] <- bound
    sum(components * weights)
  }, numeric(1))
  expected_hi <- vapply(seq_len(nrow(sens$summary$corner_cases)), function(i) {
    weights <- reference_weights
    bound <- sens$summary$corner_cases$hi[[i]]
    weights[-i] <- weights[-i] * (1 - bound) / sum(weights[-i])
    weights[[i]] <- bound
    sum(components * weights)
  }, numeric(1))
  expect_equal(
    sens$summary$corner_cases$composite_at_lo,
    expected_lo
  )
  expect_equal(
    sens$summary$corner_cases$composite_at_hi,
    expected_hi
  )
})

test_that("sensitivity_analysis verdict matches expectation for zero-cleavage protein", {
  no_cleave <- "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  expect_warning(
    res <- evaluate_digest(
      no_cleave, enzyme = "trypsin", missed_cleavages = 1L
    ),
    class = "pepvet_warning_no_cleavage_sites"
  )
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
  expect_named(
    sens$summary,
    c(
      "top1_stability", "kendall_mean", "kendall_defined_fraction",
      "reference_composites"
    )
  )
  expect_true("trypsin" %in% names(sens$summary$top1_stability))
  expect_true(sens$summary$kendall_mean > 0)
  expect_true(sens$summary$kendall_defined_fraction > 0)
})

test_that("sensitivity_analysis on batch input returns per-protein instability", {
  batch <- .fix_batch_small
  withr::local_seed(42)
  sens <- sensitivity_analysis(batch, n_iter = 500L, chunk_size = 50L)
  expect_named(sens, c("per_protein", "summary", "settings"))
  expect_s3_class(sens$per_protein, "tbl_df")
  expect_true("verdict_instability" %in% names(sens$per_protein))
  expect_true("protein_id" %in% names(sens$per_protein))
  expect_equal(nrow(sens$per_protein), nrow(batch))
  expect_named(
    sens$summary,
    c(
      "total_instability", "mean_ci_width", "ci_width_quantiles",
      "reference_weights"
    )
  )
})

test_that("sensitivity_analysis rejects invalid tuning parameters", {
  res <- .fix_bsa_trypsin

  invalid_calls <- list(
    nu_wrong_type = list(nu = "63"),
    nu_zero = list(nu = 0),
    nu_negative = list(nu = -1),
    nu_below_numeric_floor = list(nu = 0.01),
    nu_multiple = list(nu = c(1, 2)),
    nu_missing = list(nu = NA_real_),
    nu_nan = list(nu = NaN),
    nu_infinite = list(nu = Inf),
    n_iter_wrong_type = list(n_iter = "10"),
    n_iter_zero = list(n_iter = 0L),
    n_iter_negative = list(n_iter = -1L),
    n_iter_fractional = list(n_iter = 1.5),
    n_iter_missing = list(n_iter = NA_real_),
    n_iter_nan = list(n_iter = NaN),
    n_iter_infinite = list(n_iter = Inf),
    chunk_wrong_type = list(chunk_size = "10"),
    chunk_zero = list(chunk_size = 0L),
    chunk_negative = list(chunk_size = -1L),
    chunk_fractional = list(chunk_size = 1.5),
    chunk_missing = list(chunk_size = NA_real_),
    chunk_nan = list(chunk_size = NaN),
    importance_wrong_type = list(importance = "yes"),
    importance_multiple = list(importance = c(TRUE, FALSE)),
    importance_missing = list(importance = NA),
    corner_cases_wrong_type = list(corner_cases = "no"),
    corner_cases_multiple = list(corner_cases = c(TRUE, FALSE)),
    corner_cases_missing = list(corner_cases = NA)
  )

  for (input_name in names(invalid_calls)) {
    expect_error(
      do.call(
        sensitivity_analysis,
        c(list(x = res), invalid_calls[[input_name]])
      ),
      class = "pepvet_error_invalid_sensitivity_parameter",
      info = input_name
    )
  }
})

test_that("Dirichlet sampler produces correct mean and variance", {
  w0 <- c(a = 0.2, b = 0.3, c = 0.5)
  nu <- 100
  withr::local_seed(42)
  samples <- .rdirichlet(10000, nu * w0)
  expect_equal(nrow(samples), 10000)
  expect_equal(ncol(samples), 3L)
  expect_true(all(samples > 0))
  expect_equal(rowSums(samples), rep(1, 10000), tolerance = 1e-10)
  col_means <- colMeans(samples)
  expect_equal(as.numeric(col_means), as.numeric(w0), tolerance = 0.02)
  col_vars <- apply(samples, 2, var)
  expected_var <- as.numeric(w0 * (1 - w0) / (nu + 1))
  expect_equal(as.numeric(col_vars), expected_var, tolerance = 0.005)

  withr::local_seed(42)
  boundary_samples <- .rdirichlet(100L, c(20, 30, 0))
  expect_true(all(boundary_samples[, 3L] == 0))
  expect_equal(rowSums(boundary_samples), rep(1, 100L), tolerance = 1e-10)

  expect_error(
    .rdirichlet(2L, c(0, 0, 0)),
    class = "pepvet_error_invalid_sensitivity_parameter"
  )

  withr::local_seed(42)
  first <- .rdirichlet(10L, nu * w0)
  withr::local_seed(42)
  second <- .rdirichlet(10L, nu * w0)
  expect_identical(first, second)
})

test_that("plot_weight_sensitivity returns a ggplot for single-protein input", {
  skip_if_not_installed("ggplot2")
  res <- .fix_bsa_trypsin
  withr::local_seed(42)
  sens <- sensitivity_analysis(res, n_iter = 500L)
  p <- plot_weight_sensitivity(sens)
  expect_s3_class(p, "ggplot")
})

test_that("plot_weight_sensitivity marks the stored reference composite", {
  skip_if_not_installed("ggplot2")
  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(.fix_bsa_trypsin, n_iter = 80L)
  sensitivity$summary$reference_composite <- 0.123
  sensitivity$summary$composite_mean <- 0.987
  plot <- plot_weight_sensitivity(sensitivity)
  rug_layers <- vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomRug"),
    logical(1)
  )

  expect_identical(sum(rug_layers), 1L)
  expect_identical(plot$layers[[which(rug_layers)]]$data$x, 0.123)
})

test_that("sensitivity_analysis preserves custom reference weights and score", {
  weights <- c(
    S_length = 0.40,
    S_coverage = 0.20,
    S_count = 0.15,
    S_hydro = 0.15,
    S_charge = 0.10
  )
  result <- evaluate_digest(.bsa_path, weights = weights)
  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(result, n_iter = 80L)

  expect_identical(sensitivity$summary$reference_weights, weights)
  expect_equal(
    sensitivity$summary$reference_composite,
    result$scores$composite_score[[1L]]
  )
  expect_identical(
    names(sensitivity$summary$verdict_pct),
    c("Good", "Moderate", "Poor")
  )
  sampled_weights <- as.numeric(
    sensitivity$iterations[1L, names(weights), drop = TRUE]
  )
  components <- as.numeric(result$scores[1L, names(weights), drop = TRUE])
  weighted_sum <- sum(sampled_weights * components)
  expect_equal(
    sensitivity$iterations$composite_score[[1L]],
    weighted_sum
  )

  preset_result <- do.call(
    evaluate_digest,
    c(list(sequence = .bsa_path), pepvet_preset("fractionated"))
  )
  withr::local_seed(42)
  preset_sensitivity <- sensitivity_analysis(preset_result, n_iter = 20L)
  expect_identical(
    preset_sensitivity$summary$reference_weights,
    preset_result$params$weights
  )
})

test_that("sensitivity_analysis applies the zero-cleavage hard-fail rule", {
  no_cleave <- strrep("A", 20L)
  expect_warning(
    result <- evaluate_digest(no_cleave, missed_cleavages = 1L),
    class = "pepvet_warning_no_cleavage_sites"
  )
  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(result, n_iter = 80L)

  expect_identical(sensitivity$summary$reference_composite, 0)
  expect_true(all(sensitivity$iterations$composite_score == 0))
  expect_identical(sensitivity$summary$reference_verdict, "Poor")
})

test_that("multi-enzyme sensitivity uses stored reference composites", {
  comparison <- tibble::tibble(
    enzyme = factor(c("enzyme_a", "enzyme_b"), levels = c("enzyme_a", "enzyme_b")),
    protein_id = c("protein_1", "protein_1"),
    S_length = c(0.8, 0.2),
    S_coverage = c(0.8, 0.2),
    S_count = c(0.8, 0.2),
    S_hydro = c(0.8, 0.2),
    S_charge = c(0.8, 0.2),
    composite_score = c(0.9, 0.1),
    verdict = c("Good", "Poor")
  )
  attr(comparison, "scoring_config") <- list(
    weights = c(
      S_length = 0.20,
      S_coverage = 0.20,
      S_count = 0.20,
      S_hydro = 0.20,
      S_charge = 0.20
    )
  )
  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(comparison, n_iter = 80L)

  expect_equal(
    unname(sensitivity$summary$reference_composites),
    c(0.9, 0.1)
  )
  expect_identical(
    names(sensitivity$summary$top1_stability),
    c("enzyme_a", "enzyme_b")
  )

  withr::local_seed(42)
  optional <- sensitivity_analysis(
    comparison,
    n_iter = 20L,
    importance = TRUE,
    corner_cases = TRUE
  )
  expect_named(
    optional$summary$weight_importance,
    c("enzyme_a", "enzyme_b")
  )
  expect_named(optional$summary$corner_cases, c("enzyme_a", "enzyme_b"))
  expect_named(
    optional$summary$weight_importance$enzyme_a,
    c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge")
  )
  expect_s3_class(optional$summary$corner_cases$enzyme_a, "tbl_df")
})

test_that("batch sensitivity covers chunk boundaries, importance, and determinism", {
  batch <- .fix_batch_small
  withr::local_seed(42)
  first <- sensitivity_analysis(
    batch,
    n_iter = 30L,
    chunk_size = 17L,
    importance = TRUE
  )
  withr::local_seed(42)
  second <- sensitivity_analysis(
    batch,
    n_iter = 30L,
    chunk_size = 17L,
    importance = TRUE
  )
  withr::local_seed(42)
  one_chunk <- sensitivity_analysis(
    batch,
    n_iter = 30L,
    chunk_size = nrow(batch),
    importance = TRUE
  )

  expect_identical(first, second)
  expect_identical(first$per_protein, one_chunk$per_protein)
  expect_identical(first$summary, one_chunk$summary)
  expect_identical(first$settings$chunk_size, 17L)
  expect_identical(one_chunk$settings$chunk_size, nrow(batch))
  expect_equal(nrow(first$per_protein), nrow(batch))
  expect_equal(
    first$per_protein$reference_composite,
    batch$composite_score
  )
  expect_true(all(
    first$per_protein$composite_lo <= first$per_protein$composite_hi
  ))
  expect_true(all(is.finite(first$per_protein$composite_lo)))
  expect_true(all(is.finite(first$per_protein$composite_hi)))
  expect_named(
    first$summary$weight_importance,
    c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge")
  )
  expect_true(all(is.finite(first$summary$weight_importance)))
  withr::local_seed(42)
  reference_weights <- attr(batch, "scoring_config")$weights
  draws <- pepVet:::.rdirichlet(30L, 63 * reference_weights)
  components <- as.matrix(batch[, names(reference_weights), drop = FALSE])
  composites <- components %*% t(draws)
  composites[batch$S_count == 0, ] <- 0
  expected_associations <- vapply(seq_along(reference_weights), function(j) {
    mean(vapply(seq_len(nrow(batch)), function(i) {
      if (stats::sd(composites[i, ]) == 0) {
        return(0)
      }
      stats::cor(draws[, j], composites[i, ])^2
    }, numeric(1)))
  }, numeric(1))
  names(expected_associations) <- names(reference_weights)
  expect_equal(
    first$summary$weight_importance,
    expected_associations,
    tolerance = 1e-12
  )
  expect_identical(
    first$summary$reference_weights,
    attr(batch, "scoring_config")$weights
  )
})

test_that("batch sensitivity matches independent fixed-draw arithmetic", {
  weights <- stats::setNames(rep(0.2, 5L), c(
    "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
  ))
  batch <- tibble::tibble(
    protein_id = c("constant_good", "boundary_case", "hard_fail"),
    S_length = c(1, 1, 1),
    S_coverage = c(1, 0, 1),
    S_count = c(1, 1, 0),
    S_hydro = c(1, 0, 1),
    S_charge = c(1, 0, 1),
    composite_score = c(1, 0.4, 0),
    verdict = c("Good", "Moderate", "Poor")
  )
  attr(batch, "scoring_config") <- list(weights = weights)
  draws <- rbind(
    c(0.70, 0.10, 0.10, 0.05, 0.05),
    c(0.10, 0.10, 0.10, 0.35, 0.35),
    rep(0.20, 5L),
    c(0.30, 0.10, 0.30, 0.15, 0.15)
  )
  testthat::local_mocked_bindings(
    .rdirichlet = function(n, alpha) {
      expect_identical(n, 4L)
      expect_equal(alpha, 63 * weights)
      draws
    },
    .package = "pepVet"
  )

  sensitivity <- sensitivity_analysis(
    batch,
    n_iter = 4L,
    chunk_size = 2L
  )

  expect_equal(
    sensitivity$per_protein$verdict_instability,
    c(0, 0.5, 0)
  )
  expect_equal(
    sensitivity$per_protein$composite_mean,
    c(1, 0.5, 0)
  )
  expect_equal(
    sensitivity$per_protein$composite_lo,
    c(1, 0.215, 0)
  )
  expect_equal(
    sensitivity$per_protein$composite_hi,
    c(1, 0.785, 0)
  )
  expect_equal(sensitivity$summary$total_instability, 1 / 6)
})

test_that("batch sensitivity rejects unsupported corner diagnostics", {
  expect_error(
    sensitivity_analysis(
      .fix_batch_small,
      n_iter = 20L,
      corner_cases = TRUE
    ),
    class = "pepvet_error_invalid_sensitivity_parameter"
  )
})

test_that("sensitivity_analysis rejects malformed score inputs", {
  malformed_table <- tibble::tibble(
    protein_id = "protein_1",
    enzyme = "enzyme_a"
  )
  malformed_result <- list(
    scores = malformed_table,
    params = list()
  )

  expect_error(
    sensitivity_analysis(malformed_table, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    sensitivity_analysis(malformed_result, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    sensitivity_analysis(
      .fix_bsa_trypsin,
      n_iter = 1L,
      importance = TRUE
    ),
    class = "pepvet_error_invalid_sensitivity_parameter"
  )
  expect_error(
    sensitivity_analysis(42, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    sensitivity_analysis(NULL, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    sensitivity_analysis(tibble::tibble(), n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  wrong_id_type <- .fix_batch_bsa
  wrong_id_type$protein_id <- 1L
  blank_id <- .fix_batch_bsa
  blank_id$protein_id[[1L]] <- ""
  missing_id <- .fix_batch_bsa
  missing_id$protein_id[[1L]] <- NA_character_
  expect_error(
    sensitivity_analysis(wrong_id_type, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    sensitivity_analysis(blank_id, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    sensitivity_analysis(missing_id, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  wrong_score_type <- .fix_batch_bsa
  wrong_score_type$S_length <- as.character(wrong_score_type$S_length)
  expect_error(
    sensitivity_analysis(wrong_score_type, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  invalid_score <- .fix_batch_bsa
  invalid_score$S_length[[1L]] <- Inf
  expect_error(
    sensitivity_analysis(invalid_score, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  hard_fail <- .fix_batch_bsa
  hard_fail$S_count[[1L]] <- 0
  hard_fail$S_coverage[[1L]] <- 1
  hard_fail$composite_score[[1L]] <- 0.5
  hard_fail$verdict[[1L]] <- "Moderate"
  expect_error(
    sensitivity_analysis(hard_fail, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  inconsistent_verdict <- .fix_batch_bsa
  inconsistent_verdict$composite_score[[1L]] <- 0.9
  inconsistent_verdict$verdict[[1L]] <- "Poor"
  expect_error(
    sensitivity_analysis(inconsistent_verdict, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  duplicate_protein <- .fix_batch_bsa[rep(1L, 2L), , drop = FALSE]
  duplicate_protein$protein_id <- c("same", "same")
  expect_error(
    sensitivity_analysis(duplicate_protein, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  duplicate_columns <- as.data.frame(.fix_batch_bsa)
  names(duplicate_columns)[[2L]] <- names(duplicate_columns)[[1L]]
  expect_error(
    sensitivity_analysis(duplicate_columns, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  duplicate_pair <- tibble::tibble(
    enzyme = c("enzyme_a", "enzyme_a"),
    protein_id = c("protein_1", "protein_1"),
    S_length = c(0.5, 0.5),
    S_coverage = c(0.5, 0.5),
    S_count = c(0.5, 0.5),
    S_hydro = c(0.5, 0.5),
    S_charge = c(0.5, 0.5),
    composite_score = c(0.5, 0.5),
    verdict = c("Moderate", "Moderate")
  )
  expect_error(
    sensitivity_analysis(duplicate_pair, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  wrong_enzyme_type <- .fix_batch_bsa
  wrong_enzyme_type$enzyme <- 1L
  expect_error(
    sensitivity_analysis(wrong_enzyme_type, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  missing_enzyme <- .fix_batch_bsa
  missing_enzyme$enzyme <- NA_character_
  expect_error(
    sensitivity_analysis(missing_enzyme, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  blank_enzyme <- .fix_batch_bsa
  blank_enzyme$enzyme <- ""
  expect_error(
    sensitivity_analysis(blank_enzyme, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  malformed_params <- .fix_bsa_trypsin
  malformed_params$params <- "not a list"
  expect_error(
    sensitivity_analysis(malformed_params, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  malformed_metadata <- .fix_batch_bsa
  attr(malformed_metadata, "scoring_config") <- "not a list"
  expect_error(
    sensitivity_analysis(malformed_metadata, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
  multi_result <- evaluate_digest(Biostrings::AAStringSet(c(
    first = "AKAAAAAAK",
    second = "AKRTPK"
  )))
  expect_error(
    sensitivity_analysis(multi_result, n_iter = 10L),
    class = "pepvet_error_invalid_input"
  )
})

test_that("sensitivity_analysis supports proteome-aware and legacy weights", {
  proteome <- digest_protein(.multi_path)
  proteome_batch <- batch_evaluate(.multi_path, proteome = proteome)
  withr::local_seed(42)
  aware <- sensitivity_analysis(proteome_batch, n_iter = 20L)

  expect_named(
    aware$summary$reference_weights,
    c(
      "S_length", "S_coverage", "S_count",
      "S_hydro", "S_charge", "S_unique"
    )
  )

  legacy <- .fix_batch_small
  attr(legacy, "scoring_config") <- NULL
  withr::local_seed(42)
  legacy_sensitivity <- sensitivity_analysis(legacy, n_iter = 20L)
  expect_identical(
    legacy_sensitivity$summary$reference_weights,
    c(
      S_length = 0.200,
      S_coverage = 0.348,
      S_count = 0.226,
      S_hydro = 0.138,
      S_charge = 0.088
    )
  )

  legacy_aware <- proteome_batch
  attr(legacy_aware, "scoring_config") <- NULL
  withr::local_seed(42)
  legacy_aware_sensitivity <- sensitivity_analysis(
    legacy_aware,
    n_iter = 20L
  )
  expect_true(
    "S_unique" %in% names(legacy_aware_sensitivity$summary$reference_weights)
  )
})

test_that("sensitivity helpers reject invalid shared weight draws", {
  expect_error(
    .sensitivity_single(
      .fix_bsa_trypsin$scores,
      .fix_bsa_trypsin$params,
      nu = 63,
      n_iter = 2L,
      importance = FALSE,
      corner_cases = FALSE,
      weight_draws = matrix(0, nrow = 1L, ncol = 5L)
    ),
    class = "pepvet_error_invalid_sensitivity_parameter"
  )

  valid_draws <- matrix(rep(0.2, 10L), nrow = 2L, ncol = 5L)
  invalid_draws <- list(
    wrong_sum = {
      draws <- valid_draws
      draws[1L, 1L] <- 0.3
      draws
    },
    negative = {
      draws <- valid_draws
      draws[1L, 1L] <- -0.1
      draws
    },
    non_finite = {
      draws <- valid_draws
      draws[1L, 1L] <- Inf
      draws
    }
  )
  for (draws in invalid_draws) {
    expect_error(
      .sensitivity_single(
        .fix_bsa_trypsin$scores,
        .fix_bsa_trypsin$params,
        nu = 63,
        n_iter = 2L,
        importance = FALSE,
        corner_cases = FALSE,
        weight_draws = draws
      ),
      class = "pepvet_error_invalid_sensitivity_parameter"
    )
  }
})

test_that("multi-enzyme sensitivity handles tied reference ranks", {
  comparison <- tibble::tibble(
    enzyme = factor(c("enzyme_a", "enzyme_b"), levels = c("enzyme_a", "enzyme_b")),
    protein_id = c("protein_1", "protein_1"),
    S_length = c(0.5, 0.5),
    S_coverage = c(0.5, 0.5),
    S_count = c(0.5, 0.5),
    S_hydro = c(0.5, 0.5),
    S_charge = c(0.5, 0.5),
    composite_score = c(0.5, 0.5),
    verdict = c("Moderate", "Moderate")
  )
  attr(comparison, "scoring_config") <- list(
    weights = c(
      S_length = 0.20,
      S_coverage = 0.20,
      S_count = 0.20,
      S_hydro = 0.20,
      S_charge = 0.20
    )
  )
  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(comparison, n_iter = 20L)

  expect_true(is.na(sensitivity$summary$kendall_mean))
  expect_identical(sensitivity$summary$kendall_defined_fraction, 0)
  expect_equal(
    unname(sensitivity$summary$top1_stability),
    c(0.5, 0.5)
  )
})

test_that("multi-enzyme sensitivity does not call a tied reference rank stable", {
  comparison <- tibble::tibble(
    enzyme = c("enzyme_a", "enzyme_b"),
    protein_id = c("protein_1", "protein_1"),
    S_length = c(0.8, 0.2),
    S_coverage = c(0.2, 0.8),
    S_count = c(0.5, 0.5),
    S_hydro = c(0.5, 0.5),
    S_charge = c(0.5, 0.5),
    composite_score = c(0.5, 0.5),
    verdict = c("Moderate", "Moderate")
  )
  attr(comparison, "scoring_config") <- list(
    weights = stats::setNames(rep(0.2, 5L), c(
      "S_length", "S_coverage", "S_count", "S_hydro", "S_charge"
    ))
  )

  withr::local_seed(42)
  sensitivity <- sensitivity_analysis(comparison, n_iter = 500L)

  expect_true(is.na(sensitivity$summary$kendall_mean))
  expect_identical(sensitivity$summary$kendall_defined_fraction, 0)
  expect_true(all(sensitivity$summary$top1_stability > 0.4))
  expect_equal(sum(sensitivity$summary$top1_stability), 1)
})

test_that("multi-enzyme sensitivity shares one draw matrix across enzymes", {
  comparison <- tibble::tibble(
    enzyme = factor(c("enzyme_a", "enzyme_b")),
    protein_id = c("protein_1", "protein_1"),
    S_length = c(0.8, 0.2),
    S_coverage = c(0.8, 0.2),
    S_count = c(0.8, 0.2),
    S_hydro = c(0.8, 0.2),
    S_charge = c(0.8, 0.2),
    composite_score = c(0.8, 0.2),
    verdict = c("Good", "Poor")
  )
  attr(comparison, "scoring_config") <- list(
    weights = c(
      S_length = 0.20,
      S_coverage = 0.20,
      S_count = 0.20,
      S_hydro = 0.20,
      S_charge = 0.20
    )
  )
  draws_seen <- new.env(parent = emptyenv())
  draws_seen$values <- list()

  testthat::local_mocked_bindings(
    .rdirichlet = function(n, alpha) {
      matrix(rep(0.2, n * length(alpha)), nrow = n, ncol = length(alpha))
    },
    .sensitivity_single = function(scores, params, nu, n_iter,
                                   importance, corner_cases,
                                   weight_draws = NULL) {
      draws_seen$values <- c(draws_seen$values, list(weight_draws))
      list(
        iterations = tibble::tibble(
          composite_score = rep(scores$composite_score[[1L]], n_iter)
        ),
        summary = list(reference_composite = scores$composite_score[[1L]])
      )
    },
    .package = "pepVet"
  )

  .sensitivity_enzymes(
    comparison,
    nu = 63,
    n_iter = 3L,
    importance = FALSE,
    corner_cases = FALSE,
    scoring_config = attr(comparison, "scoring_config")
  )

  expect_length(draws_seen$values, 2L)
  expect_identical(draws_seen$values[[1L]], draws_seen$values[[2L]])
})

# Argument passthrough.

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
  expect_warning(
    default <- batch_compare_enzymes(small, enzymes = c("trypsin", "lysc")),
    class = "pepvet_warning_no_cleavage_sites"
  )
  expect_warning(
    wider <- batch_compare_enzymes(
      small,
      enzymes = c("trypsin", "lysc"),
      gravy_range = c(-2.0, 2.0)
    ),
    class = "pepvet_warning_no_cleavage_sites"
  )
  expect_true(wider$S_hydro[[1]] >= default$S_hydro[[1]])
})

test_that("pepvet_check passes ... to evaluate_digest", {
  bsa <- .bsa_path
  default <- pepvet_check(bsa, enzyme = "trypsin")
  wider <- pepvet_check(bsa, enzyme = "trypsin", gravy_range = c(-2.0, 2.0))
  expect_true(wider$scores$S_hydro >= default$scores$S_hydro)
})

# Error class tests.

test_that("batch_evaluate rejects invalid cores", {
  invalid_cores <- list(
    NULL, character(0), "1", 0L, -1L, 1.5, NA_real_, NaN, Inf
  )
  for (cores in invalid_cores) {
    expect_error(
      batch_evaluate(.bsa_path, cores = cores),
      class = "pepvet_error_invalid_cores"
    )
    expect_error(
      batch_compare_enzymes(.bsa_path, cores = cores),
      class = "pepvet_error_invalid_cores"
    )
  }
})

test_that("compare_digests rejects invalid enzymes", {
  invalid_enzymes <- list(
    NULL,
    character(0),
    NA_character_,
    c("trypsin", NA_character_),
    1
  )
  for (enzymes in invalid_enzymes) {
    expect_error(
      compare_digests(.bsa_path, enzymes = enzymes),
      class = "pepvet_error_invalid_enzymes"
    )
  }
  expect_error(
    compare_digests(
      .bsa_path,
      enzymes = c("Trypsin", " trypsin ")
    ),
    class = "pepvet_error_invalid_enzymes"
  )
  expect_error(
    compare_digests(.bsa_path, enzymes = "not-an-enzyme"),
    class = "pepvet_error_invalid_enzyme"
  )
})

test_that("batch_compare_enzymes rejects invalid enzyme vectors and IDs", {
  expect_error(
    batch_compare_enzymes(.small_path, enzymes = character(0)),
    class = "pepvet_error_invalid_enzymes"
  )
  expect_error(
    batch_compare_enzymes(.small_path, enzymes = NA_character_),
    class = "pepvet_error_invalid_enzymes"
  )
  duplicated <- Biostrings::AAStringSet(c(
    duplicate = "AKAAAAAAK",
    duplicate = "AKRTPK"
  ))
  expect_error(
    batch_compare_enzymes(duplicated, enzymes = c("trypsin", "lysc")),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    batch_compare_enzymes(character(0), enzymes = c("trypsin", "lysc")),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    batch_compare_enzymes(
      .small_path,
      enzymes = c("Trypsin", " trypsin ")
    ),
    class = "pepvet_error_invalid_enzymes"
  )
  expect_error(
    batch_compare_enzymes(.small_path, enzymes = "not-an-enzyme"),
    class = "pepvet_error_invalid_enzyme"
  )
})
