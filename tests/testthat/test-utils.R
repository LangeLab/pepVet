data("aa_properties", package = "pepVet", envir = environment())

test_that(
  "range validators normalize valid boundaries and reject invalid ranges",
  {
    expect_identical(
      pepVet:::.validate_gravy_range(c(-1L, 1L)),
      c(-1, 1)
    )
    expect_identical(
      pepVet:::.validate_gravy_range(c(0.4, 0.4)),
      c(0.4, 0.4)
    )
    expect_identical(
      pepVet:::.validate_length_range(c(7, 25)),
      c(7L, 25L)
    )
    expect_identical(
      pepVet:::.validate_length_range(c(10L, 10L)),
      c(10L, 10L)
    )

    invalid_cases <- list(
      list(
        name = "GRAVY wrong type",
        validator = pepVet:::.validate_gravy_range,
        value = "-1,0.6",
        class = "pepvet_error_invalid_gravy_range"
      ),
      list(
        name = "GRAVY wrong length",
        validator = pepVet:::.validate_gravy_range,
        value = c(-1, 0, 0.6),
        class = "pepvet_error_invalid_gravy_range"
      ),
      list(
        name = "GRAVY missing value",
        validator = pepVet:::.validate_gravy_range,
        value = c(-1, NA_real_),
        class = "pepvet_error_invalid_gravy_range"
      ),
      list(
        name = "GRAVY non-finite value",
        validator = pepVet:::.validate_gravy_range,
        value = c(-1, Inf),
        class = "pepvet_error_invalid_gravy_range"
      ),
      list(
        name = "GRAVY descending values",
        validator = pepVet:::.validate_gravy_range,
        value = c(0.6, -1),
        class = "pepvet_error_invalid_gravy_range"
      ),
      list(
        name = "length wrong type",
        validator = pepVet:::.validate_length_range,
        value = c("7", "25"),
        class = "pepvet_error_invalid_length_range"
      ),
      list(
        name = "length wrong length",
        validator = pepVet:::.validate_length_range,
        value = 7,
        class = "pepvet_error_invalid_length_range"
      ),
      list(
        name = "length missing value",
        validator = pepVet:::.validate_length_range,
        value = c(7, NA_real_),
        class = "pepvet_error_invalid_length_range"
      ),
      list(
        name = "length non-finite value",
        validator = pepVet:::.validate_length_range,
        value = c(7, Inf),
        class = "pepvet_error_invalid_length_range"
      ),
      list(
        name = "length descending values",
        validator = pepVet:::.validate_length_range,
        value = c(25, 7),
        class = "pepvet_error_invalid_length_range"
      ),
      list(
        name = "length non-integer value",
        validator = pepVet:::.validate_length_range,
        value = c(7.5, 25),
        class = "pepvet_error_invalid_length_range"
      ),
      list(
        name = "length below one",
        validator = pepVet:::.validate_length_range,
        value = c(0, 25),
        class = "pepvet_error_invalid_length_range"
      )
    )

    for (case in invalid_cases) {
      expect_error(
        case$validator(case$value),
        class = case$class,
        info = case$name
      )
    }

    expect_no_warning(
      expect_error(
        pepVet:::.validate_length_range(c(1e20, 1e20)),
        class = "pepvet_error_invalid_length_range"
      )
    )
  }
)

test_that("sequence helpers normalize names and supported residues", {
  expect_identical(
    pepVet:::.normalize_sequence_names(
      c(" alpha ", NA_character_, "", "beta"),
      sequence_count = 4L
    ),
    c("alpha", "sequence_2", "sequence_3", "beta")
  )
  expect_identical(
    pepVet:::.normalize_sequence_names(NULL, sequence_count = 2L),
    c("sequence_1", "sequence_2")
  )
  expect_identical(
    pepVet:::.normalize_sequence_names(
      c("duplicate", "duplicate"),
      sequence_count = 2L
    ),
    c("duplicate", "duplicate")
  )

  expect_identical(
    pepVet:::.validate_sequence("acduo", sequence_name = "example"),
    "ACDUO"
  )
  expect_identical(
    pepVet:::.normalize_peptide_sequences(
      c(first = "acduo", second = "KR")
    ),
    c("ACDUO", "KR")
  )
  expect_error(
    pepVet:::.validate_sequence("A Z", sequence_name = "example"),
    class = "pepvet_error_invalid_sequence"
  )

  invalid_inputs <- list(
    null = NULL,
    empty = character(0),
    wrong_type = 42,
    missing = NA_character_,
    blank = ""
  )
  for (input_name in names(invalid_inputs)) {
    expect_error(
      pepVet:::.normalize_peptide_sequences(invalid_inputs[[input_name]]),
      class = "pepvet_error_invalid_sequence",
      info = input_name
    )
  }
})

test_that("read_input accepts sequence classes and controlled FASTA paths", {
  character_input <- pepVet:::.read_input(
    c(alpha = "ac", beta = "duo")
  )
  expect_s4_class(character_input, "AAStringSet")
  expect_identical(
    as.character(character_input),
    c(alpha = "AC", beta = "DUO")
  )

  aa_string_input <- pepVet:::.read_input(Biostrings::AAString("acduo"))
  expect_identical(as.character(aa_string_input), c(sequence_1 = "ACDUO"))

  aa_set_input <- pepVet:::.read_input(
    Biostrings::AAStringSet(c(one = "ac", two = "duo"))
  )
  expect_identical(
    as.character(aa_set_input),
    c(one = "AC", two = "DUO")
  )

  temp_root <- withr::local_tempdir("pepVet utility path ")
  fasta_path <- file.path(temp_root, "input with spaces.fasta")
  writeLines(c(">fixture header", "acde"), fasta_path)

  absolute_input <- pepVet:::.read_input(fasta_path)
  expect_identical(
    as.character(absolute_input),
    c(`fixture header` = "ACDE")
  )

  withr::local_dir(temp_root)
  relative_input <- pepVet:::.read_input(basename(fasta_path))

  expect_identical(
    as.character(relative_input),
    c(`fixture header` = "ACDE")
  )

  expect_error(
    pepVet:::.read_input(NULL),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    pepVet:::.read_input(42),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    pepVet:::.read_input(character(0)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    pepVet:::.read_input(file.path(temp_root, "missing.fasta")),
    class = "pepvet_error_missing_file"
  )
  expect_error(
    pepVet:::.read_input(temp_root),
    class = "pepvet_error_missing_file"
  )
})

test_that("build_digest_ranges enumerates adjacent missed-cleavage intervals", {
  strict_ranges <- IRanges::IRanges(
    start = c(1L, 3L, 5L),
    end = c(2L, 4L, 6L)
  )
  expected <- list(
    start = c(1L, 1L, 1L, 3L, 3L, 5L),
    end = c(2L, 4L, 6L, 4L, 6L, 6L),
    missed_cleavages = c(0L, 1L, 2L, 0L, 1L, 0L)
  )

  expect_identical(
    pepVet:::.build_digest_ranges(strict_ranges, missed_cleavages = 0L),
    list(
      start = c(1L, 3L, 5L),
      end = c(2L, 4L, 6L),
      missed_cleavages = c(0L, 0L, 0L)
    )
  )
  expect_identical(
    pepVet:::.build_digest_ranges(strict_ranges, missed_cleavages = 2L),
    expected
  )
  expect_identical(
    pepVet:::.build_digest_ranges(strict_ranges, missed_cleavages = 10L),
    expected
  )

  empty <- pepVet:::.build_digest_ranges(
    IRanges::IRanges(),
    missed_cleavages = 2L
  )
  expect_identical(empty, list(
    start = integer(0),
    end = integer(0),
    missed_cleavages = integer(0)
  ))
})

test_that("hydrophobicity lookup caching has explicit missing-value behavior", {
  cache <- get(".pepvet_cache", envir = asNamespace("pepVet"))
  had_cached_lookup <- exists(
    "hydro_lookup", envir = cache, inherits = FALSE
  )
  previous_lookup <- if (had_cached_lookup) cache$hydro_lookup else NULL
  withr::defer({
    if (had_cached_lookup) {
      cache$hydro_lookup <- previous_lookup
    } else {
      rm(list = "hydro_lookup", envir = cache)
    }
  })

  rm(list = "hydro_lookup", envir = cache)
  first <- pepVet:::.get_hydro_lookup()
  cached_lookup <- cache$hydro_lookup
  second <- pepVet:::.get_hydro_lookup()

  expect_identical(first, cached_lookup)
  expect_identical(second, cached_lookup)
  expect_identical(names(first), aa_properties$amino_acid)
  expect_equal(
    first[c("A", "I", "U", "O")],
    c(A = 1.8, I = 4.5, U = 2.5, O = NA_real_)
  )

  observed <- pepVet:::.calculate_gravy(
    c(lower = "aliv", mixed = "AO", missing = "O")
  )
  expected <- c(
    lower = mean(c(1.8, 3.8, 4.5, 4.2)),
    mixed = 1.8,
    missing = NA_real_
  )
  expect_identical(names(observed), names(expected))
  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("ionizable composition and net charge follow independent equations", {
  sequences <- c(
    mixed = "ACDEHIKR",
    tyrosine = "YYYY",
    neutral = "AAA",
    selenocysteine = "UUUU"
  )
  residue_names <- c("C", "D", "E", "H", "K", "R", "Y", "U")
  expected_composition <- matrix(
    c(
      1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L,
      0L, 0L, 0L, 0L, 0L, 0L, 4L, 0L,
      0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L, 0L, 0L, 0L, 4L
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(NULL, residue_names)
  )

  observed_composition <- pepVet:::.ionizable_composition_matrix(sequences)
  expect_type(observed_composition, "integer")
  expect_identical(observed_composition, expected_composition)
  expect_identical(
    pepVet:::.ionizable_composition_matrix("acdehikr"),
    expected_composition[1L, , drop = FALSE]
  )
  expect_error(
    pepVet:::.ionizable_composition_matrix(character(0)),
    class = "pepvet_error_invalid_sequence"
  )
  invalid_inputs <- list(
    null = NULL,
    wrong_type = 42,
    missing = NA_character_,
    blank = ""
  )
  for (input_name in names(invalid_inputs)) {
    expect_error(
      pepVet:::.ionizable_composition_matrix(invalid_inputs[[input_name]]),
      class = "pepvet_error_invalid_sequence",
      info = input_name
    )
  }

  composition <- matrix(
    c(
      0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L,
      0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
      0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(NULL, residue_names)
  )
  pH <- 7
  terminal_charge <-
    1 / (1 + 10^(pH - 8.0)) - 1 / (1 + 10^(3.1 - pH))
  expected_charge <- c(
    terminal_charge,
    terminal_charge + 1 / (1 + 10^(pH - 10.5)),
    terminal_charge - 1 / (1 + 10^(3.9 - pH)),
    terminal_charge - 1 / (1 + 10^(5.2 - pH))
  )
  observed_charge <- pepVet:::.net_charge_at_pH(
    rep(pH, nrow(composition)),
    composition
  )

  expect_equal(observed_charge, expected_charge, tolerance = 1e-12)
  expect_true(observed_charge[[2L]] > observed_charge[[1L]])
  expect_true(observed_charge[[3L]] < observed_charge[[1L]])
  expect_true(observed_charge[[4L]] < observed_charge[[1L]])
})

test_that("row binding and proteome indexing preserve compact table contracts", {
  empty <- pepVet:::.bind_rows(list())
  expect_s3_class(empty, "tbl_df")
  expect_identical(dim(empty), c(0L, 0L))

  first <- data.frame(value = 1:2, row.names = c("a", "b"))
  second <- data.frame(value = 3:4, row.names = c("c", "d"))
  single <- pepVet:::.bind_rows(list(first))
  expect_s3_class(single, "tbl_df")
  expect_identical(single$value, 1:2)

  combined <- pepVet:::.bind_rows(list(first, second))
  expect_s3_class(combined, "tbl_df")
  expect_identical(combined$value, 1:4)
  expect_identical(rownames(combined), as.character(seq_len(4L)))

  proteome <- data.frame(
    peptide = c("AA", "BB", "AA", "AA"),
    protein_id = c("p1", "p2", "p1", "p3")
  )
  index <- pepVet:::.build_proteome_index(proteome)
  expect_true(is.environment(index))
  expect_setequal(get("AA", envir = index), c("p1", "p3"))
  expect_identical(get("BB", envir = index), "p2")
  expect_false(exists("CC", envir = index, inherits = FALSE))

  empty_index <- pepVet:::.build_proteome_index(proteome[0, , drop = FALSE])
  expect_length(ls(empty_index, all.names = TRUE), 0L)
})
