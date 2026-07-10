expected_reference_fasta_files <- c(
  "P00698.fasta",
  "P02769.fasta",
  "P0CG48.fasta",
  "P37840_isoforms.fasta",
  "P56817.fasta",
  "P68431.fasta",
  "Q8WZ42.fasta",
  "small_proteome_50_proteins.fasta"
)

expected_single_entry_widths <- c(
  P00698 = 147L,
  P02769 = 607L,
  P0CG48 = 685L,
  P56817 = 501L,
  P68431 = 136L,
  Q8WZ42 = 34350L
)

expected_isoform_accessions <- c("P37840", "P37840-2", "P37840-3")
expected_isoform_widths <- c(140L, 112L, 126L)

expected_small_proteome_widths <- c(
  O15162 = 318L,
  O15393 = 492L,
  O60502 = 916L,
  O75475 = 530L,
  O75970 = 2070L,
  O95273 = 360L,
  O95905 = 644L,
  P01889 = 362L,
  P01911 = 266L,
  P04439 = 365L,
  P05556 = 798L,
  P06400 = 928L,
  P08174 = 381L,
  P08246 = 267L,
  P08648 = 1049L,
  P10321 = 366L,
  P13667 = 645L,
  P15529 = 392L,
  P16234 = 1089L,
  P17544 = 483L,
  P25685 = 340L,
  P29320 = 983L,
  P29590 = 882L,
  P31358 = 61L,
  P31629 = 2446L,
  P31689 = 397L,
  P41226 = 1012L,
  P42694 = 1942L,
  P51784 = 963L,
  P56851 = 147L,
  P63244 = 317L,
  P63279 = 158L,
  P78310 = 365L,
  Q08648 = 103L,
  Q15326 = 602L,
  Q6FHJ7 = 346L,
  Q6P2E9 = 1401L,
  Q6ZWK4 = 172L,
  Q7L5Y9 = 396L,
  Q86V81 = 257L,
  Q86WV6 = 379L,
  Q92692 = 538L,
  Q92956 = 283L,
  Q99836 = 296L,
  Q9H3P7 = 528L,
  Q9H9K5 = 563L,
  Q9NPC3 = 277L,
  Q9NX65 = 697L,
  Q9UM44 = 414L,
  Q9Y2H6 = 1198L
)

standard_amino_acids <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1]]

reference_fasta_path <- function(file_name) {
  system.file("extdata", file_name, package = "pepVet")
}

read_reference_fasta <- function(file_name) {
  Biostrings::readAAStringSet(reference_fasta_path(file_name))
}

extract_accessions <- function(fasta_names) {
  sub("^[^|]+[|]([^|]+)[|].*$", "\\1", fasta_names)
}

extract_sequence_alphabet <- function(fasta_set) {
  sort(
    unique(
      strsplit(
        paste(as.character(fasta_set), collapse = ""),
        "",
        fixed = TRUE
      )[[1]]
    )
  )
}

fasta_widths <- function(fasta_set) {
  unname(nchar(as.character(fasta_set), type = "chars"))
}

capture_messages <- function(expr) {
  captured <- list()

  withCallingHandlers(
    expr,
    message = function(cnd) {
      captured[[length(captured) + 1L]] <<- cnd
      invokeRestart("muffleMessage")
    }
  )

  captured
}

test_that("reference FASTA inventory is exact and resolvable", {
  extdata_dir <- system.file("extdata", package = "pepVet")

  expect_true(nzchar(extdata_dir))
  expect_true(dir.exists(extdata_dir))
  expect_identical(
    sort(list.files(extdata_dir, pattern = "[.]fasta$")),
    expected_reference_fasta_files
  )

  for (file_name in expected_reference_fasta_files) {
    path <- reference_fasta_path(file_name)
    expect_true(nzchar(path), info = file_name)
    expect_true(file.exists(path), info = path)
    expect_match(basename(path), file_name, fixed = TRUE)
  }
})

test_that(
  "single-entry reference FASTAs have the expected accessions and widths",
  {
    for (accession in names(expected_single_entry_widths)) {
      fasta <- read_reference_fasta(paste0(accession, ".fasta"))

      expect_true(methods::is(fasta, "AAStringSet"), info = accession)
      expect_identical(length(fasta), 1L, info = accession)
      expect_identical(
        extract_accessions(names(fasta)),
        accession,
        info = accession
      )
      expect_identical(
        fasta_widths(fasta),
        expected_single_entry_widths[[accession]]
      )
      expect_false(anyNA(names(fasta)), info = accession)
      expect_true(nzchar(names(fasta)), info = accession)
    }
  }
)

test_that("alpha-synuclein isoform FASTA is pinned exactly", {
  fasta <- read_reference_fasta("P37840_isoforms.fasta")

  expect_s4_class(fasta, "AAStringSet")
  expect_length(fasta, 3L)
  expect_identical(
    extract_accessions(names(fasta)),
    expected_isoform_accessions
  )
  expect_identical(fasta_widths(fasta), expected_isoform_widths)
  expect_equal(anyDuplicated(names(fasta)), 0L)
})

test_that("small proteome FASTA has the expected 50 accessions and widths", {
  fasta <- read_reference_fasta("small_proteome_50_proteins.fasta")
  accessions <- extract_accessions(names(fasta))

  expect_s4_class(fasta, "AAStringSet")
  expect_length(fasta, 50L)
  expect_equal(anyDuplicated(names(fasta)), 0L)
  expect_identical(accessions, names(expected_small_proteome_widths))
  expect_identical(fasta_widths(fasta), unname(expected_small_proteome_widths))
})

test_that("all committed FASTA fixtures use the standard amino-acid alphabet", {
  for (file_name in expected_reference_fasta_files) {
    fasta <- read_reference_fasta(file_name)
    alphabet <- extract_sequence_alphabet(fasta)

    expect_true(all(fasta_widths(fasta) > 0L), info = file_name)
    expect_true(
      setequal(alphabet, intersect(alphabet, standard_amino_acids)),
      info = file_name
    )
  }
})

test_that("key reference fixture widths remain biologically plausible", {
  titin <- read_reference_fasta("Q8WZ42.fasta")
  isoforms <- read_reference_fasta("P37840_isoforms.fasta")

  expect_gt(fasta_widths(titin), 30000L)
  expect_identical(fasta_widths(isoforms), expected_isoform_widths)
})
