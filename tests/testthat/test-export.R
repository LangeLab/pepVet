# export_peptide_list - Skyline.

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

test_that("skyline mz values are positive and finite", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  result <- export_peptide_list(bsa_peps, format = "skyline")

  expect_true(all(is.finite(result$`Precursor Mz`)))
  expect_true(all(result$`Precursor Mz` > 0))
})

test_that("skyline mz increases with higher charge is FALSE (heavier per charge at +2 vs +3)", {
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
  peps <- digest_protein("MAAAAAAAKAAAAAAAR") # two 8-AA tryptic peptides, both valid
  result <- export_peptide_list(peps, format = "fasta")

  header_lines <- result[seq(1L, length(result), by = 2L)]
  # all headers must contain a pipe-separated coordinate suffix
  expect_true(all(grepl("\\|peptide_[0-9]+-[0-9]+$", header_lines)))
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
  tmp <- tempfile(fileext = ".csv")
  result <- withVisible(export_peptide_list(bsa_peps, format = "skyline", file = tmp))

  expect_false(result$visible)
  expect_identical(result$value, tmp)
  expect_true(file.exists(tmp))
  unlink(tmp)
})

test_that("fasta format writes a file and returns file path invisibly", {
  bsa_peps <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  tmp <- tempfile(fileext = ".fasta")
  result <- withVisible(export_peptide_list(bsa_peps, format = "fasta", file = tmp))

  expect_false(result$visible)
  expect_true(file.exists(tmp))
  lines <- readLines(tmp)
  expect_true(any(startsWith(lines, ">")))
  unlink(tmp)
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
})

test_that("non-data-frame peptides input raises a classed error", {
  expect_error(
    export_peptide_list("not a tibble"),
    class = "pepvet_error_invalid_export_input"
  )
})

test_that("missing required columns raises a classed error", {
  bad_df <- data.frame(peptide = "AAAAAK", stringsAsFactors = FALSE)

  expect_error(
    export_peptide_list(bad_df),
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
})
