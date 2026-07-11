# export_peptide_list - Skyline.

export_fixture <- function() {
  tibble::tibble(
    protein_id = c("protein_1", "protein_2"),
    peptide = c("AAAAAAA", "ACDEFGH"),
    start = c(1L, 10L),
    end = c(7L, 16L),
    length = c(7L, 7L),
    missed_cleavages = c(0L, 1L)
  )
}

test_that("skyline format returns a tibble with the correct column names", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  result <- export_peptide_list(bsa_peps, format = "skyline")

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c("Protein", "Peptide Sequence", "Precursor Charge", "Precursor Mz")
  )
})

test_that("charges argument controls the number of rows per valid peptide", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  valid_count <- sum(bsa_peps$length >= 7L & bsa_peps$length <= 25L)

  result_two <- export_peptide_list(bsa_peps, format = "skyline", charges = 2:3)
  result_three <- export_peptide_list(bsa_peps, format = "skyline", charges = 2:4)

  expect_equal(nrow(result_two), valid_count * 2L)
  expect_equal(nrow(result_three), valid_count * 3L)
})

test_that("skyline export preserves peptide and charge order", {
  result <- export_peptide_list(
    export_fixture(),
    format = " skyline ",
    charges = c(2L, 4L),
    length_range = c(7L, 7L)
  )

  expect_identical(
    result$Protein,
    c("protein_1", "protein_1", "protein_2", "protein_2")
  )
  expect_identical(
    result$`Precursor Charge`,
    c(2L, 4L, 2L, 4L)
  )
  expect_identical(
    result$`Peptide Sequence`,
    c("AAAAAAA", "AAAAAAA", "ACDEFGH", "ACDEFGH")
  )
})

test_that("skyline m/z follows the neutral-mass and charge equation", {
  peptide <- "ACDEFGH"
  residues <- strsplit(peptide, "", fixed = TRUE)[[1L]]
  independent_residue_masses <- c(
    A = 71.03712, C = 103.00919, D = 115.02695,
    E = 129.04260, F = 147.06842, G = 57.02147,
    H = 137.05892
  )
  neutral_mass <- sum(independent_residue_masses[residues]) + 18.01056
  result <- export_peptide_list(
    export_fixture()[2, , drop = FALSE],
    format = "skyline",
    charges = c(2L, 3L)
  )
  expected <- (neutral_mass + c(2, 3) * 1.007276) / c(2, 3)

  expect_equal(result$`Precursor Mz`, expected, tolerance = 1e-10)
})

test_that("skyline mz values are positive and finite", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  result <- export_peptide_list(bsa_peps, format = "skyline")

  expect_true(all(is.finite(result$`Precursor Mz`)))
  expect_true(all(result$`Precursor Mz` > 0))
})

test_that("skyline m/z at charge 2 exceeds charge 3", {
  # For a given peptide, mz at charge 2 is generally higher than at charge 3
  peps <- digest_protein("MAAAKAAAARAAAAK")
  result <- export_peptide_list(peps, format = "skyline", charges = 2:3)

  charge2_mz <- result$`Precursor Mz`[result$`Precursor Charge` == 2L]
  charge3_mz <- result$`Precursor Mz`[result$`Precursor Charge` == 3L]
  # neutral_mass/2 + proton > neutral_mass/3 + proton for typical peptide masses
  expect_true(all(charge2_mz > charge3_mz))
})

test_that("skyline format returns empty tibble when no valid peptides exist", {
  # A very short sequence that produces no valid-length peptides
  peps <- digest_protein("MAAK")
  result <- export_peptide_list(peps, format = "skyline")

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0L)
  expect_named(
    result,
    c("Protein", "Peptide Sequence", "Precursor Charge", "Precursor Mz")
  )
})

# export_peptide_list - generic.

test_that("generic format returns all rows annotated with gravy, pI, and valid columns", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  result <- export_peptide_list(bsa_peps, format = "generic")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(bsa_peps))
  expect_true(all(c("gravy", "pI", "valid") %in% names(result)))
  expect_true(all(is.finite(result$gravy)))
  expect_true(all(is.finite(result$pI)))
})

test_that("generic valid column matches the length_range filter", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  result <- export_peptide_list(bsa_peps, format = "generic")
  expected_valid <- bsa_peps$length >= 7L & bsa_peps$length <= 25L

  expect_equal(result$valid, expected_valid)
})

test_that("generic export preserves columns and applies inclusive boundaries", {
  fixture <- export_fixture()
  result <- export_peptide_list(
    fixture,
    format = "generic",
    length_range = c(7L, 7L),
    charges = Inf
  )

  expect_named(result, c(names(fixture), "gravy", "pI", "valid"))
  expect_identical(result$valid, c(TRUE, TRUE))
  expect_true(all(is.finite(result$gravy)))
  expect_true(all(is.finite(result$pI)))
})

test_that("generic export marks values on opposite sides of the length window", {
  fixture <- tibble::tibble(
    protein_id = c("short", "boundary"),
    peptide = c("AAAAAA", "ACDEFGH"),
    start = c(1L, 1L),
    end = c(6L, 7L),
    length = c(6L, 7L),
    missed_cleavages = c(0L, 0L)
  )

  result <- export_peptide_list(
    fixture,
    format = "generic",
    length_range = c(7L, 7L)
  )

  expect_identical(result$valid, c(FALSE, TRUE))
})

test_that("generic and FASTA exports support an empty validated table", {
  empty <- tibble::tibble(
    protein_id = character(),
    peptide = character(),
    start = integer(),
    end = integer(),
    length = integer(),
    missed_cleavages = integer()
  )

  generic <- export_peptide_list(empty, format = "generic")
  fasta <- export_peptide_list(empty, format = "fasta")

  expect_identical(generic$gravy, numeric())
  expect_identical(generic$pI, numeric())
  expect_identical(generic$valid, logical())
  expect_type(fasta, "character")
  expect_length(fasta, 0L)
})

# export_peptide_list - FASTA.

test_that("fasta format returns two lines per valid peptide with > headers", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  valid_count <- sum(bsa_peps$length >= 7L & bsa_peps$length <= 25L)
  result <- export_peptide_list(bsa_peps, format = "fasta")

  expect_type(result, "character")
  expect_equal(length(result), valid_count * 2L)
  header_lines <- result[seq(1L, length(result), by = 2L)]
  expect_true(all(startsWith(header_lines, ">")))
})

test_that("fasta headers encode protein_id and start-end coordinates", {
  peps <- digest_protein("MAAAAAAAKAAAAAAAR") # two valid tryptic peptides
  result <- export_peptide_list(peps, format = "fasta")

  header_lines <- result[seq(1L, length(result), by = 2L)]
  # all headers must contain a pipe-separated coordinate suffix
  expect_true(all(grepl("\\|peptide_[0-9]+-[0-9]+$", header_lines)))
})

test_that("fasta export emits exact headers and sequences", {
  result <- export_peptide_list(export_fixture(), format = "fasta")

  expect_identical(
    result,
    c(
      ">protein_1|peptide_1-7", "AAAAAAA",
      ">protein_2|peptide_10-16", "ACDEFGH"
    )
  )
})

test_that("fasta format returns empty character vector when no valid peptides", {
  peps <- digest_protein("MAAK")
  result <- export_peptide_list(peps, format = "fasta")

  expect_type(result, "character")
  expect_length(result, 0L)
})

# export_peptide_list - file argument.

test_that("skyline format writes a file and returns file path invisibly", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  temp_root <- withr::local_tempdir()
  tmp <- file.path(temp_root, "skyline.csv")
  result <- withVisible(export_peptide_list(bsa_peps, format = "skyline", file = tmp))

  expect_false(result$visible)
  expect_identical(result$value, tmp)
  expect_true(file.exists(tmp))
})

test_that("fasta format writes a file and returns file path invisibly", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  temp_root <- withr::local_tempdir()
  tmp <- file.path(temp_root, "peptides.fasta")
  result <- withVisible(export_peptide_list(bsa_peps, format = "fasta", file = tmp))

  expect_false(result$visible)
  expect_true(file.exists(tmp))
  lines <- readLines(tmp)
  expect_true(any(startsWith(lines, ">")))
})

test_that("CSV and FASTA file output preserves the returned schema", {
  temp_root <- withr::local_tempdir()
  csv_path <- file.path(temp_root, "export.csv")
  fasta_path <- file.path(temp_root, "export.fasta")

  csv_result <- withVisible(export_peptide_list(
    export_fixture(), format = "generic", file = csv_path
  ))
  fasta_result <- withVisible(export_peptide_list(
    export_fixture(), format = "fasta", file = fasta_path
  ))

  expect_false(csv_result$visible)
  expect_false(fasta_result$visible)
  expect_identical(csv_result$value, csv_path)
  expect_identical(fasta_result$value, fasta_path)
  expect_named(utils::read.csv(csv_path, check.names = FALSE),
    c("protein_id", "peptide", "start", "end", "length",
      "missed_cleavages", "gravy", "pI", "valid")
  )
  expect_identical(readLines(fasta_path), export_peptide_list(
    export_fixture(), format = "fasta"
  ))
})

# export_peptide_list - error handling.

test_that("unsupported format raises a classed error", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")

  expect_error(
    export_peptide_list(bsa_peps, format = "xlsx"),
    class = "pepvet_error_invalid_export_format"
  )
})

test_that("invalid charges raises a classed error", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")

  expect_error(
    export_peptide_list(bsa_peps, format = "skyline", charges = 0L),
    class = "pepvet_error_invalid_charges"
  )
  expect_error(
    export_peptide_list(bsa_peps, format = "skyline", charges = integer(0)),
    class = "pepvet_error_invalid_charges"
  )

  for (invalid in list(
    NULL, "2", NA_real_, Inf, -Inf, 1.5, .Machine$integer.max + 1
  )) {
    expect_error(
      export_peptide_list(bsa_peps, format = "skyline", charges = invalid),
      class = "pepvet_error_invalid_charges"
    )
  }
})

test_that("format validation rejects missing, empty, and unsupported values", {
  fixture <- export_fixture()

  for (invalid in list(
    NULL, character(0), NA_character_, "", "   ",
    c("skyline", "fasta"), "xlsx"
  )) {
    expect_error(
      export_peptide_list(fixture, format = invalid),
      class = "pepvet_error_invalid_export_format"
    )
  }
})

test_that("non-data-frame peptides input raises a classed error", {
  for (invalid in list(NULL, NA, "not a tibble")) {
    expect_error(
      export_peptide_list(invalid),
      class = "pepvet_error_invalid_export_input"
    )
  }
})

test_that("missing required columns raises a classed error", {
  bad_df <- data.frame(peptide = "AAAAAK", stringsAsFactors = FALSE)

  expect_error(
    export_peptide_list(bad_df),
    class = "pepvet_error_invalid_export_input"
  )
})

test_that("export input validation rejects malformed schemas and values", {
  fixture <- export_fixture()

  invalid_inputs <- list(
    wrong_protein_type = within(fixture, protein_id <- factor(protein_id)),
    wrong_peptide_type = within(fixture, peptide <- factor(peptide)),
    wrong_length_type = within(fixture, length <- factor(length)),
    missing_length = fixture[, setdiff(names(fixture), "length")],
    missing_sequence = within(fixture, peptide <- NA_character_),
    fractional_length = within(fixture, length <- c(7.5, 7)),
    mismatched_length = within(fixture, length <- c(6L, 7L)),
    invalid_sequence = within(fixture, peptide <- "ZZZZZZZ"),
    empty_identifier = within(fixture, protein_id <- c("", "protein_2"))
  )

  for (invalid in invalid_inputs) {
    expect_error(
      export_peptide_list(invalid),
      class = "pepvet_error_invalid_export_input"
    )
  }

  expect_error(
    export_peptide_list(fixture[, setdiff(names(fixture), c("start", "end"))],
      format = "fasta"),
    class = "pepvet_error_invalid_export_input"
  )
  bad_coordinates <- within(fixture, end <- c(6L, 16L))
  expect_error(
    export_peptide_list(bad_coordinates, format = "fasta"),
    class = "pepvet_error_invalid_export_input"
  )
  bad_coordinate_values <- within(fixture, start <- c(0L, 10L))
  expect_error(
    export_peptide_list(bad_coordinate_values, format = "fasta"),
    class = "pepvet_error_invalid_export_input"
  )
})

test_that("duplicate export columns are rejected", {
  fixture <- export_fixture()
  names(fixture)[[2L]] <- names(fixture)[[1L]]

  expect_error(
    export_peptide_list(fixture),
    class = "pepvet_error_invalid_export_input"
  )
})

test_that("export_peptide_list rejects invalid file argument", {
  peps <- digest_protein(.bsa_path, enzyme = "trypsin")
  expect_error(
    export_peptide_list(peps, file = NA),
    class = "pepvet_error_invalid_file"
  )
  expect_error(
    export_peptide_list(peps, file = 123),
    class = "pepvet_error_invalid_file"
  )
  expect_error(
    export_peptide_list(peps, file = ""),
    class = "pepvet_error_invalid_file"
  )
  expect_error(
    export_peptide_list(peps, file = "   "),
    class = "pepvet_error_invalid_file"
  )
  invalid_parent <- file.path(
    withr::local_tempdir(), "missing-directory", "export.csv"
  )
  expect_error(
    export_peptide_list(peps, format = "generic", file = invalid_parent),
    class = "pepvet_error_invalid_file"
  )
})
