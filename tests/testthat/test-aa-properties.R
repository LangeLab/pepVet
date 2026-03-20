data("aa_properties", package = "pepVet", envir = environment())

expected_amino_acids <- c(
  "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
)

expected_free_mass <- c(
  A = 89.04768, C = 121.01975, D = 133.03751, E = 147.05316,
  F = 165.07898, G = 75.03203, H = 155.06948, I = 131.09463,
  K = 146.10553, L = 131.09463, M = 149.05105, N = 132.05349,
  P = 115.06333, Q = 146.06914, R = 174.11168, S = 105.04259,
  T = 119.05824, V = 117.07898, W = 204.08988, Y = 181.07389
)

expected_residue_mass <- c(
  A = 71.03711, C = 103.00919, D = 115.02694, E = 129.04259,
  F = 147.06841, G = 57.02146, H = 137.05891, I = 113.08406,
  K = 128.09496, L = 113.08406, M = 131.04049, N = 114.04293,
  P = 97.05276, Q = 128.05858, R = 156.10111, S = 87.03203,
  T = 101.04768, V = 99.06841, W = 186.07931, Y = 163.06333
)

expected_kd <- c(
  A = 1.8, C = 2.5, D = -3.5, E = -3.5, F = 2.8,
  G = -0.4, H = -3.2, I = 4.5, K = -3.9, L = 3.8,
  M = 1.9, N = -3.5, P = -1.6, Q = -3.5, R = -4.5,
  S = -0.8, T = -0.7, V = 4.2, W = -0.9, Y = -1.3
)

expected_pka <- c(
  C = 8.3, D = 3.9, E = 4.3, H = 6.0, K = 10.5, R = 12.5, Y = 10.1
)

water_mass <- 18.01056

reference_gravy <- function(peptide_sequence) {
  residues <- strsplit(toupper(peptide_sequence), split = "", fixed = TRUE)[[1]]
  mean(unname(expected_kd[residues]))
}

test_that(
  "aa_properties is a complete reference table for the 20 standard residues",
  {
    expect_s3_class(aa_properties, "tbl_df")
    expect_equal(nrow(aa_properties), 20L)
    expect_equal(anyDuplicated(aa_properties$amino_acid), 0L)
    expect_identical(aa_properties$amino_acid, expected_amino_acids)
    expect_true(all(nchar(aa_properties$amino_acid) == 1L))
    expect_true(
      all(aa_properties$amino_acid == toupper(aa_properties$amino_acid))
    )
  }
)

test_that("aa_properties has the expected schema and missingness", {
  expect_identical(
    names(aa_properties),
    c(
      "amino_acid",
      "molecular_weight",
      "hydrophobicity",
      "pKa_side_chain",
      "is_basic"
    )
  )

  expect_type(aa_properties$amino_acid, "character")
  expect_type(aa_properties$molecular_weight, "double")
  expect_type(aa_properties$hydrophobicity, "double")
  expect_type(aa_properties$pKa_side_chain, "double")
  expect_type(aa_properties$is_basic, "logical")

  expect_false(anyNA(aa_properties$amino_acid))
  expect_false(anyNA(aa_properties$molecular_weight))
  expect_false(anyNA(aa_properties$hydrophobicity))
  expect_false(anyNA(aa_properties$is_basic))
})

test_that(
  "aa_properties molecular weights match standard free amino acid masses",
  {
    observed_mass <- setNames(
      aa_properties$molecular_weight,
      aa_properties$amino_acid
    )

    expect_equal(
      observed_mass,
      expected_free_mass,
      tolerance = 1e-4
    )
  }
)

test_that("aa_properties free masses are residue masses plus water", {
  observed_mass <- setNames(
    aa_properties$molecular_weight,
    aa_properties$amino_acid
  )

  expect_equal(
    unname(observed_mass - expected_residue_mass),
    rep(water_mass, length(expected_residue_mass)),
    tolerance = 1e-4
  )
})

test_that(
  "aa_properties hydrophobicity values match the Kyte-Doolittle scale",
  {
    observed_hydrophobicity <- setNames(
      aa_properties$hydrophobicity,
      aa_properties$amino_acid
    )

    expect_equal(observed_hydrophobicity, expected_kd)
  }
)

test_that("side-chain pKa values are present only for the ionizable residues", {
  observed_pka <- setNames(
    aa_properties$pKa_side_chain,
    aa_properties$amino_acid
  )
  ionizable_residues <- names(expected_pka)

  expect_identical(
    names(observed_pka)[!is.na(observed_pka)],
    ionizable_residues
  )
  expect_equal(observed_pka[ionizable_residues], expected_pka)
  expect_true(
    all(is.na(observed_pka[setdiff(expected_amino_acids, ionizable_residues)]))
  )
})

test_that("basic residue annotations agree with the chemistry table", {
  basic_residues <- aa_properties$amino_acid[aa_properties$is_basic]
  observed_pka <- setNames(
    aa_properties$pKa_side_chain,
    aa_properties$amino_acid
  )

  expect_identical(basic_residues, c("H", "K", "R"))
  expect_true(all(!is.na(observed_pka[basic_residues])))
  expect_true(all(observed_pka[basic_residues] > 0))
  expect_false(any(c("D", "E", "C", "Y") %in% basic_residues))
})

test_that(".calculate_gravy matches the reference hydrophobicity arithmetic", {
  test_sequences <- c(
    A = "A",
    ALIV = "ALIV",
    ACDE = "ACDE",
    WYRK = "WYRK",
    VILA = "VILA",
    GGGSS = "GGGSS"
  )

  for (sequence_name in names(test_sequences)) {
    sequence <- test_sequences[[sequence_name]]
    expect_equal(
      pepVet:::.calculate_gravy(sequence),
      reference_gravy(sequence),
      info = sequence_name
    )
  }
})

test_that(".calculate_gravy accepts lowercase sequences", {
  expect_equal(pepVet:::.calculate_gravy("aliv"), reference_gravy("ALIV"))
})

test_that(".calculate_gravy rejects non-scalar character inputs", {
  expect_error(pepVet:::.calculate_gravy(1), "single character string")
  expect_error(
    pepVet:::.calculate_gravy(c("A", "B")),
    "single character string"
  )
  expect_error(
    pepVet:::.calculate_gravy(character()),
    "single character string"
  )
})

test_that(".calculate_gravy rejects empty and missing strings", {
  expect_error(pepVet:::.calculate_gravy(""), "must not be empty")
  expect_error(pepVet:::.calculate_gravy(NA_character_), "must not be missing")
})

test_that(
  ".calculate_gravy rejects unknown amino acid codes with deduped messaging",
  {
    expect_error(
      pepVet:::.calculate_gravy("AXZAZX"),
      "Unknown amino acid code\\(s\\): X, Z\\."
    )
  }
)
