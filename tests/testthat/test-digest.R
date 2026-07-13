expected_enzymes <- c(
  "arg-c proteinase",
  "asp-n endopeptidase",
  "bnps-skatole",
  "caspase1",
  "caspase2",
  "caspase3",
  "caspase4",
  "caspase5",
  "caspase6",
  "caspase7",
  "caspase8",
  "caspase9",
  "caspase10",
  "chymotrypsin-high",
  "chymotrypsin-low",
  "clostripain",
  "cnbr",
  "enterokinase",
  "factor xa",
  "formic acid",
  "glutamyl endopeptidase",
  "granzyme-b",
  "hydroxylamine",
  "iodosobenzoic acid",
  "lysc",
  "lysn",
  "lysarginase",
  "neutrophil elastase",
  "ntcb",
  "pepsin1.3",
  "pepsin",
  "proline endopeptidase",
  "proteinase k",
  "staphylococcal peptidase i",
  "thermolysin",
  "thrombin",
  "trypsin",
  "trypsin-high",
  "trypsin-low",
  "trypsin-simple"
)

strict_ranges <- function(result) {
  result[result$missed_cleavages == 0L, , drop = FALSE]
}

generated_property_sequences <- function(n_records = 6L) {
  amino_acids <- strsplit(
    "ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE
  )[[1L]]
  cleavage_motif <- strsplit("AKADMAFWYRP", "", fixed = TRUE)[[1L]]
  sequence_lengths <- sample(48L:96L, n_records, replace = FALSE)

  sequences <- vapply(
    seq_len(n_records),
    function(index) {
      residues <- sample(
        amino_acids,
        sequence_lengths[[index]],
        replace = TRUE
      )
      motif_start <- sample.int(
        length(residues) - length(cleavage_motif) - 8L,
        1L
      ) + 4L
      motif_positions <- motif_start + seq_along(cleavage_motif) - 1L
      residues[motif_positions] <- cleavage_motif
      paste0(residues, collapse = "")
    },
    character(1)
  )
  names(sequences) <- sprintf("generated_%02d", seq_len(n_records))

  Biostrings::AAStringSet(sequences)
}

test_that("supported enzyme registry is pinned exactly", {
  expect_identical(pepVet:::.supported_digest_enzymes, expected_enzymes)
})

test_that("registered enzymes normalize and dispatch through digest_protein", {
  normalized <- vapply(
    expected_enzymes,
    function(enzyme) {
      pepVet:::.normalize_enzyme(paste0(" ", toupper(enzyme), " "))
    },
    character(1),
    USE.NAMES = FALSE
  )
  expect_identical(normalized, expected_enzymes)

  digests <- lapply(
    expected_enzymes,
    function(enzyme) {
      digest_protein(
        "MKWVTFISLLFLFSSAYSR",
        enzyme = enzyme,
        missed_cleavages = 0L
      )
    }
  )
  expected_columns <- c(
    "protein_id", "peptide", "start", "end", "length",
    "missed_cleavages"
  )

  expect_true(
    all(vapply(digests, nrow, integer(1)) > 0L)
  )
  expect_true(
    all(vapply(digests, function(result) {
      identical(names(result), expected_columns) && !anyNA(result)
    }, logical(1)))
  )
})

test_that("digest_protein returns the documented schema and column types", {
  result <- digest_protein("MKWVTFISLLFLFSSAYSR")

  expect_s3_class(result, "tbl_df")
  expect_identical(
    names(result),
    c("protein_id", "peptide", "start", "end", "length", "missed_cleavages")
  )
  expect_type(result$protein_id, "character")
  expect_type(result$peptide, "character")
  expect_type(result$start, "integer")
  expect_type(result$end, "integer")
  expect_type(result$length, "integer")
  expect_type(result$missed_cleavages, "integer")
  expect_false(anyNA(result))
})

test_that("annotate_cleavage_sites classifies tryptic motifs correctly", {
  bsa_path <- reference_fasta("P02769.fasta")
  annotations <- annotate_cleavage_sites(bsa_path, enzyme = "trypsin")

  expect_s3_class(annotations, "tbl_df")
  expect_identical(
    names(annotations),
    c("position", "residue", "flanking_context", "efficiency", "rule_applied")
  )

  kp_site <- annotations[annotations$position == 140L, , drop = FALSE]
  kk_site <- annotations[annotations$position == 155L, , drop = FALSE]
  default_site <- annotations[annotations$position == 2L, , drop = FALSE]

  expect_identical(kp_site$residue, "K")
  expect_identical(kp_site$efficiency, "low")
  expect_identical(kp_site$rule_applied, "proline_block")

  expect_identical(kk_site$residue, "K")
  expect_identical(kk_site$efficiency, "medium")
  expect_identical(kk_site$rule_applied, "adjacent_basic_residues")

  expect_identical(default_site$residue, "K")
  expect_identical(default_site$efficiency, "high")
  expect_identical(default_site$rule_applied, "default_trypsin_site")

  synthetic <- annotate_cleavage_sites("AKRDKPRA")
  expect_identical(synthetic$position, c(2L, 3L, 5L, 7L))
  expect_identical(synthetic$residue, c("K", "R", "K", "R"))
  expect_identical(
    synthetic$flanking_context,
    c("AKRD", "AKRDK", "RDKPR", "KPRA")
  )
  expect_identical(
    synthetic$efficiency,
    c("medium", "low", "low", "high")
  )
  expect_identical(
    synthetic$rule_applied,
    c(
      "adjacent_basic_residues", "acidic_p1_prime", "proline_block",
      "default_trypsin_site"
    )
  )

  lowercase <- annotate_cleavage_sites("akr", enzyme = " TRYPSIN ")
  expect_identical(lowercase$position, 2L)
  expect_identical(lowercase$efficiency, "medium")

  first_site <- annotate_cleavage_sites("KAAA")
  expect_identical(first_site$position, 1L)
  expect_identical(first_site$flanking_context, "KAA")
  expect_identical(first_site$efficiency, "high")
  expect_identical(first_site$rule_applied, "default_trypsin_site")

  expected_alias_output <- annotate_cleavage_sites("AKRTPK")
  for (enzyme in pepVet:::.cleavage_annotation_trypsin_enzymes) {
    expect_identical(
      annotate_cleavage_sites("AKRTPK", enzyme = enzyme),
      expected_alias_output,
      info = enzyme
    )
  }
})

test_that("annotate_cleavage_sites rejects unsupported annotation enzymes", {
  expect_error(
    annotate_cleavage_sites("MKWVTFISLLFLFSSAYSR", enzyme = "lysc"),
    class = "pepvet_error_unsupported_cleavage_annotation"
  )
})

test_that("annotate_cleavage_sites accepts input forms and empty sites", {
  sequence <- "AKRTPK"
  temp_root <- withr::local_tempdir("pepVet annotation ")
  fasta_path <- file.path(temp_root, "single protein.fasta")
  writeLines(c(">annotation_fixture", sequence), fasta_path)

  expected <- annotate_cleavage_sites(sequence)
  expect_identical(
    annotate_cleavage_sites(c(annotation_fixture = sequence)),
    expected
  )
  expect_identical(
    annotate_cleavage_sites(Biostrings::AAString(sequence)),
    expected
  )
  expect_identical(annotate_cleavage_sites(fasta_path), expected)

  empty <- annotate_cleavage_sites("AAAA")
  expect_s3_class(empty, "tbl_df")
  expect_identical(names(empty), names(expected))
  expect_identical(nrow(empty), 0L)
  expect_type(empty$position, "integer")
  expect_type(empty$efficiency, "character")
})

test_that("annotate_cleavage_sites rejects invalid sequence inputs", {
  invalid_inputs <- list(
    null = NULL,
    empty = character(0),
    wrong_type = 42,
    multiple = c(first = "AK", second = "RK")
  )
  for (input_name in names(invalid_inputs)) {
    expect_error(
      annotate_cleavage_sites(invalid_inputs[[input_name]]),
      class = "pepvet_error_invalid_input",
      info = input_name
    )
  }

  expect_error(
    annotate_cleavage_sites(NA_character_),
    class = "pepvet_error_invalid_sequence"
  )
  expect_error(
    annotate_cleavage_sites(""),
    class = "pepvet_error_invalid_sequence"
  )

  temp_root <- withr::local_tempdir("pepVet annotation errors ")
  malformed_path <- file.path(temp_root, "malformed.fasta")
  writeLines(c("not a FASTA header", "AKR"), malformed_path)
  expect_error(
    annotate_cleavage_sites(malformed_path),
    class = "pepvet_error_invalid_input",
    regexp = "FASTA"
  )

  multi_path <- file.path(temp_root, "multiple.fasta")
  writeLines(c(">first", "AKR", ">second", "TPK"), multi_path)
  expect_error(
    annotate_cleavage_sites(multi_path),
    class = "pepvet_error_invalid_input"
  )

  expect_error(
    annotate_cleavage_sites(file.path(temp_root, "missing.fasta")),
    class = "pepvet_error_missing_file"
  )
})

test_that("annotate_cleavage_sites rejects invalid enzyme values", {
  sequence <- "AKRTPK"
  invalid_enzymes <- list(
    null = NULL,
    missing = NA_character_,
    wrong_type = 42,
    empty = "",
    unrecognized = "not-an-enzyme"
  )
  for (enzyme_name in names(invalid_enzymes)) {
    expect_error(
      annotate_cleavage_sites(
        sequence,
        enzyme = invalid_enzymes[[enzyme_name]]
      ),
      class = "pepvet_error_invalid_enzyme",
      info = enzyme_name
    )
  }
})

test_that("cleavage annotation helpers preserve classification and severity", {
  annotations <- pepVet:::.annotate_trypsin_cleavage_sites("AKRDKPRA")

  expect_identical(
    vapply(
      c("trypsin", "trypsin-high", "trypsin-low", "trypsin-simple", "lysc"),
      pepVet:::.supports_cleavage_efficiency_annotations,
      logical(1),
      USE.NAMES = FALSE
    ),
    c(TRUE, TRUE, TRUE, TRUE, FALSE)
  )
  expect_identical(
    pepVet:::.least_efficient_class(c("high", "medium", "low")),
    "low"
  )
  expect_identical(
    pepVet:::.least_efficient_class(c(NA_character_, "medium")),
    "medium"
  )
  expect_identical(
    pepVet:::.least_efficient_class(character(0)),
    NA_character_
  )
  expect_identical(
    pepVet:::.least_efficient_class(NA_character_),
    NA_character_
  )

  expect_identical(
    pepVet:::.map_peptide_cleavage_efficiency(1L, 4L, annotations),
    "low"
  )
  expect_identical(
    pepVet:::.map_peptide_cleavage_efficiency(5L, 7L, annotations),
    "low"
  )
  expect_identical(
    pepVet:::.map_peptide_cleavage_efficiency(1L, 1L, annotations),
    NA_character_
  )

  expect_identical(
    pepVet:::.cleavage_efficiency_summary("AKRDKPRA", "trypsin"),
    list(n_high_efficiency_sites = 1L, n_low_efficiency_sites = 2L)
  )
  expect_identical(
    pepVet:::.cleavage_efficiency_summary("AAAA", "trypsin"),
    list(n_high_efficiency_sites = 0L, n_low_efficiency_sites = 0L)
  )
  expect_identical(
    pepVet:::.cleavage_efficiency_summary("AKRDKPRA", "lysc"),
    list(
      n_high_efficiency_sites = NA_integer_,
      n_low_efficiency_sites = NA_integer_
    )
  )
})

test_that("single protein digestion roundtrips and preserves positions", {
  fasta <- Biostrings::readAAStringSet(reference_fasta("P02769.fasta"))
  sequence <- as.character(fasta[[1]])
  result <- digest_protein(reference_fasta("P02769.fasta"), enzyme = "trypsin")
  mc0 <- strict_ranges(result)
  coverage <- IRanges::reduce(IRanges::IRanges(mc0$start, mc0$end))

  expect_identical(paste0(mc0$peptide, collapse = ""), sequence)
  expect_identical(
    mc0$peptide,
    substring(sequence, first = mc0$start, last = mc0$end)
  )
  expect_true(all(
    mc0$start >= 1L,
    mc0$start <= mc0$end,
    mc0$end <= nchar(sequence)
  ))
  expect_identical(mc0$length, as.integer(nchar(mc0$peptide)))
  expect_identical(mc0$length, mc0$end - mc0$start + 1L)
  expect_identical(IRanges::start(coverage), 1L)
  expect_identical(IRanges::end(coverage), nchar(sequence))
})

test_that("generated digests match source coordinates across enzyme families", {
  withr::local_seed(20260712)
  sequences <- generated_property_sequences()
  source_sequences <- as.character(sequences)
  enzymes <- c(
    "trypsin",
    "asp-n endopeptidase",
    "cnbr",
    "chymotrypsin-high"
  )

  for (enzyme in enzymes) {
    result <- digest_protein(
      sequences,
      enzyme = enzyme,
      missed_cleavages = 2L
    )
    source_index <- match(result$protein_id, names(source_sequences))
    row_sources <- unname(source_sequences[source_index])
    expected_peptides <- substring(
      row_sources,
      first = result$start,
      last = result$end
    )

    expect_false(anyNA(source_index), info = enzyme)
    expect_identical(unique(result$protein_id), names(sequences), info = enzyme)
    expect_true(all(result$start >= 1L), info = enzyme)
    expect_true(all(result$start <= result$end), info = enzyme)
    expect_true(
      all(result$end <= nchar(row_sources, type = "chars")),
      info = enzyme
    )
    expect_identical(result$peptide, expected_peptides, info = enzyme)
    expect_identical(
      result$length,
      as.integer(nchar(result$peptide, type = "chars")),
      info = enzyme
    )
  }
})

test_that("strict generated digests reconstruct under complete coverage", {
  withr::local_seed(20260713)
  sequences <- generated_property_sequences(n_records = 8L)
  source_sequences <- as.character(sequences)
  enzymes <- c("trypsin", "asp-n endopeptidase", "cnbr")

  for (enzyme in enzymes) {
    strict_digest <- digest_protein(
      sequences,
      enzyme = enzyme,
      missed_cleavages = 0L
    )

    for (protein_id in names(source_sequences)) {
      rows <- strict_digest[
        strict_digest$protein_id == protein_id, , drop = FALSE
      ]
      rows <- rows[order(rows$start, rows$end), , drop = FALSE]
      source_length <- nchar(
        source_sequences[[protein_id]],
        type = "chars"
      )

      ## Reconstruction requires a gap-free, non-overlapping strict partition.
      has_complete_coverage <- nrow(rows) > 0L &&
        rows$start[[1L]] == 1L &&
        rows$end[[nrow(rows)]] == source_length &&
        (
          nrow(rows) == 1L ||
            all(rows$start[-1L] == rows$end[-nrow(rows)] + 1L)
        )
      case_label <- paste(enzyme, protein_id)

      expect_true(has_complete_coverage, info = case_label)
      expect_identical(
        paste0(rows$peptide, collapse = ""),
        unname(source_sequences[[protein_id]]),
        info = case_label
      )
    }
  }
})

test_that("character and AAString input produce identical digestion results", {
  sequence <- "MKWVTFISLLFLFSSAYSR"

  expect_identical(
    digest_protein(sequence, enzyme = "trypsin"),
    digest_protein(Biostrings::AAString(sequence), enzyme = "trypsin")
  )
})

test_that("enzyme normalization accepts case and surrounding whitespace", {
  expect_identical(
    digest_protein("mkwvtfisllflfssaysr", enzyme = " Trypsin "),
    digest_protein("MKWVTFISLLFLFSSAYSR", enzyme = "trypsin")
  )
})

test_that("enzyme families follow independent cleavage oracles", {
  cases <- list(
    trypsin = list(
      sequence = "AKRTPK",
      enzyme = "trypsin",
      peptides = c("AK", "R", "TPK")
    ),
    trypsin_blocks_rp = list(
      sequence = "AKPR",
      enzyme = "trypsin",
      peptides = "AKPR"
    ),
    lysc = list(
      sequence = "AKRTPK",
      enzyme = "lysc",
      peptides = c("AK", "RTPK")
    ),
    cnbr = list(
      sequence = "AMAKMAGT",
      enzyme = "cnbr",
      peptides = c("AM", "AKM", "AGT")
    ),
    aspn = list(
      sequence = "AADDK",
      enzyme = "asp-n endopeptidase",
      peptides = c("AA", "D", "DK")
    ),
    glutamyl = list(
      sequence = "ADEEAK",
      enzyme = "glutamyl endopeptidase",
      peptides = c("ADE", "E", "AK")
    ),
    chymotrypsin = list(
      sequence = "AFYWLPA",
      enzyme = "chymotrypsin-high",
      peptides = c("AF", "Y", "W", "LPA")
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    result <- digest_protein(
      case$sequence,
      enzyme = case$enzyme,
      missed_cleavages = 0L
    )
    expect_identical(result$peptide, case$peptides, info = case_name)
    expect_identical(
      result$length,
      as.integer(nchar(result$peptide)),
      info = case_name
    )
    expect_identical(
      result$peptide,
      substring(case$sequence, result$start, result$end),
      info = case_name
    )
  }
})

test_that("named character vectors preserve their names as protein ids", {
  result <- digest_protein(
    c(foo = "AKRTPK", bar = "MKWVTFISLLFLFSSAYSR"),
    enzyme = "trypsin"
  )

  expect_identical(unique(result$protein_id), c("foo", "bar"))
})

test_that("AAStringSet input preserves sequence names across records", {
  input <- Biostrings::AAStringSet(
    c(alpha = "MKWVTFISLLFLFSSAYSR", beta = "AKRTPK")
  )
  result <- digest_protein(input, enzyme = "trypsin")

  expect_identical(unique(result$protein_id), c("alpha", "beta"))
  expect_true(all(result$protein_id %in% c("alpha", "beta")))
})

test_that("single-entry FASTA path preserves the FASTA header as protein_id", {
  fasta_path <- reference_fasta("P02769.fasta")
  fasta_header <- names(Biostrings::readAAStringSet(fasta_path))[[1]]
  result <- digest_protein(fasta_path, enzyme = "trypsin")

  expect_identical(unique(result$protein_id), fasta_header)
})

test_that("multi-entry FASTA inputs preserve all records", {
  isoform_result <- digest_protein(
    reference_fasta("P37840_isoforms.fasta"),
    enzyme = "trypsin"
  )
  proteome_result <- digest_protein(
    reference_fasta("small_proteome_50_proteins.fasta"),
    enzyme = "trypsin"
  )

  expect_identical(length(unique(isoform_result$protein_id)), 3L)
  expect_identical(length(unique(proteome_result$protein_id)), 50L)
  expect_true(all(table(isoform_result$protein_id) > 0L))
  expect_true(all(table(proteome_result$protein_id) > 0L))
})

test_that("digest_protein accepts irregular FASTA headers and extensions", {
  temp_fasta <- tempfile(fileext = ".txt")
  writeLines(c(
    ">weird header no pipes",
    "mkwvtfisllflfssaysr"
  ), temp_fasta)

  result <- digest_protein(temp_fasta, enzyme = "trypsin",
                           missed_cleavages = 0L)

  expect_identical(unique(result$protein_id), "weird header no pipes")
  expect_identical(result$peptide, c("MK", "WVTFISLLFLFSSAYSR"))
})

test_that("trypsin does not cleave a KP motif", {
  result <- digest_protein("AAAAAKPAAAAAAAR", enzyme = "trypsin",
                           missed_cleavages = 0L)

  expect_identical(nrow(result), 1L)
  expect_identical(result$peptide, "AAAAAKPAAAAAAAR")
  expect_match(result$peptide, "KP", fixed = TRUE)
})

test_that("repeated single-residue peptides are preserved without collapse", {
  result <- digest_protein("RKRKRKRK", enzyme = "trypsin",
                           missed_cleavages = 0L)

  expect_identical(nrow(result), 8L)
  expect_identical(result$peptide, c("R", "K", "R", "K", "R", "K", "R", "K"))
  expect_true(all(result$length == 1L))
})

test_that("trypsin and Lys-C produce meaningfully different digest patterns", {
  sequence <- as.character(
    Biostrings::readAAStringSet(reference_fasta("P02769.fasta"))[[1]]
  )
  trypsin_result <- strict_ranges(digest_protein(sequence, enzyme = "trypsin"))
  lysc_result <- strict_ranges(digest_protein(sequence, enzyme = "lysc"))
  lysc_nonterminal <- lysc_result$peptide[-nrow(lysc_result)]
  lysc_terminal_residue <- substring(
    lysc_nonterminal,
    nchar(lysc_nonterminal),
    nchar(lysc_nonterminal)
  )

  expect_gt(nrow(trypsin_result), nrow(lysc_result))
  expect_true(all(lysc_terminal_residue == "K"))
})

test_that("known reference digests stay pinned for BSA and lysozyme", {
  bsa <- strict_ranges(digest_protein(reference_fasta("P02769.fasta"),
                                      missed_cleavages = 0L))
  lysozyme <- strict_ranges(digest_protein(reference_fasta("P00698.fasta"),
                                           missed_cleavages = 0L))

  expect_identical(nrow(bsa), 79L)
  expect_identical(
    head(bsa$peptide, 5),
    c("MK", "WVTFISLLLLFSSAYSR", "GVFR", "R", "DTHK")
  )
  expect_identical(
    tail(bsa$peptide, 5),
    c("ATEEQLK", "TVMENFVAFVDK", "CCAADDK", "EACFAVEGPK", "LVVSTQTALA")
  )

  expect_identical(nrow(lysozyme), 19L)
  expect_identical(
    head(lysozyme$peptide, 5),
    c("MR", "SLLILVLCFLPLAALGK", "VFGR", "CELAAAMK", "R")
  )
  expect_identical(
    tail(lysozyme$peptide, 5),
    c("NR", "CK", "GTDVQAWIR", "GCR", "L")
  )
})

test_that("no-cleavage and single-residue edge cases stay stable", {
  no_cut <- digest_protein("AAAAAAAAAA", enzyme = "trypsin")
  single <- digest_protein("M", enzyme = "trypsin")

  expect_identical(no_cut$peptide, "AAAAAAAAAA")
  expect_identical(no_cut$start, 1L)
  expect_identical(no_cut$end, 10L)
  expect_identical(no_cut$length, 10L)

  expect_identical(single$peptide, "M")
  expect_identical(single$start, 1L)
  expect_identical(single$end, 1L)
  expect_identical(single$length, 1L)

  collision_root <- withr::local_tempdir("pepVet sequence collision ")
  dir.create(file.path(collision_root, "R"))
  withr::local_dir(collision_root)

  collision_digest <- digest_protein(
    "R",
    enzyme = "trypsin",
    missed_cleavages = 0L
  )
  expect_identical(collision_digest$peptide, "R")
  expect_identical(collision_digest$start, 1L)
  expect_identical(collision_digest$end, 1L)

  collision_annotations <- annotate_cleavage_sites("R")
  expect_identical(nrow(collision_annotations), 0L)
})

test_that(
  "missed cleavage output includes strict peptides and adjacent concatenations",
  {
    result <- digest_protein(
      "AKRTPK",
      enzyme = "trypsin",
      missed_cleavages = 1L
    )
    mc0 <- result[result$missed_cleavages == 0L, , drop = FALSE]
    mc1 <- result[result$missed_cleavages == 1L, , drop = FALSE]

    expect_identical(mc0$peptide, c("AK", "R", "TPK"))
    expect_identical(mc1$peptide, c("AKR", "RTPK"))
    expect_true(all(mc0$peptide %in% result$peptide))
    expect_gt(nrow(result), nrow(mc0))
    expect_identical(mc1$start, c(1L, 3L))
    expect_identical(mc1$end, c(3L, 6L))
    expect_identical(
      result$peptide,
      substring("AKRTPK", result$start, result$end)
    )
    expect_true(all(
      result$start >= 1L,
      result$start <= result$end,
      result$end <= nchar("AKRTPK")
    ))
    expect_identical(result$length, as.integer(nchar(result$peptide)))
    expect_true(all(result$missed_cleavages %in% 0:1))
    expect_identical(
      nrow(result),
      nrow(unique(result[c("start", "end", "missed_cleavages")]))
    )
  }
)

test_that("digest_protein can append peptide-level cleavage efficiency", {
  expected_columns <- c(
    "protein_id", "peptide", "start", "end", "length", "missed_cleavages",
    "cleavage_efficiency"
  )

  for (enzyme in pepVet:::.cleavage_annotation_trypsin_enzymes) {
    result <- digest_protein(
      "AKRTPK",
      enzyme = enzyme,
      missed_cleavages = 0L,
      include_cleavage_efficiency = TRUE
    )

    expect_identical(names(result), expected_columns, info = enzyme)
    expect_identical(
      result$cleavage_efficiency,
      c("medium", "medium", "high"),
      info = enzyme
    )
  }
})

test_that("unsupported enzymes get NA peptide-level cleavage efficiency", {
  result <- digest_protein(
    "AKRTPK",
    enzyme = "lysc",
    include_cleavage_efficiency = TRUE
  )

  expect_true("cleavage_efficiency" %in% names(result))
  expect_type(result$cleavage_efficiency, "character")
  expect_true(all(is.na(result$cleavage_efficiency)))
})

test_that("no-site trypsin efficiency stays typed missing", {
  result <- digest_protein(
    "AAAA",
    enzyme = "trypsin",
    include_cleavage_efficiency = TRUE
  )

  expect_identical(result$peptide, "AAAA")
  expect_type(result$cleavage_efficiency, "character")
  expect_true(all(is.na(result$cleavage_efficiency)))
})

test_that("missing and invalid FASTA paths raise a classed file error", {
  expect_error(
    digest_protein("missing_fixture.fasta"),
    class = "pepvet_error_missing_file"
  )

  temp_dir <- withr::local_tempdir("pepVet digest path ")
  expect_error(digest_protein(temp_dir), class = "pepvet_error_missing_file")

  malformed_path <- file.path(temp_dir, "malformed.fasta")
  writeLines(c("not a FASTA header", "MKWV"), malformed_path)
  expect_error(
    digest_protein(malformed_path),
    class = "pepvet_error_invalid_input",
    regexp = "FASTA"
  )
})

test_that("invalid enzyme names raise a classed error with supported names", {
  invalid_enzymes <- list(
    null = NULL,
    empty = character(0),
    missing = NA_character_,
    blank = "",
    wrong_type = 42,
    unsupported = "not-an-enzyme"
  )

  for (enzyme_name in names(invalid_enzymes)) {
    expect_error(
      digest_protein(
        "MKWVTFISLLFLFSSAYSR",
        enzyme = invalid_enzymes[[enzyme_name]]
      ),
      regexp = if (identical(enzyme_name, "unsupported")) {
        "Supported enzymes"
      } else {
        NULL
      },
      class = "pepvet_error_invalid_enzyme",
      info = enzyme_name
    )
  }
})

test_that("invalid sequence content is rejected with offending characters", {
  expect_error(
    digest_protein("MXBZ123", enzyme = "trypsin"),
    regexp = "unsupported amino acid code",
    class = "pepvet_error_invalid_sequence"
  )

  expect_error(
    digest_protein(""),
    class = "pepvet_error_invalid_sequence"
  )
})

test_that("invalid input types and missed cleavage values are rejected", {
  expect_error(digest_protein(NULL), class = "pepvet_error_invalid_input")
  expect_error(digest_protein(12345), class = "pepvet_error_invalid_input")
  expect_error(
    digest_protein(character(0)),
    class = "pepvet_error_invalid_input"
  )
  expect_error(
    digest_protein(NA_character_),
    class = "pepvet_error_invalid_sequence"
  )
  expect_error(
    digest_protein(Biostrings::AAStringSet()),
    class = "pepvet_error_invalid_input"
  )

  invalid_missed_cleavages <- list(
    negative = -1L,
    fractional = 1.5,
    missing = NA_real_,
    non_finite = Inf,
    empty = numeric(0),
    multiple = c(0L, 1L),
    wrong_type = "1"
  )
  for (value_name in names(invalid_missed_cleavages)) {
    expect_error(
      digest_protein(
        "MKWVTFISLLFLFSSAYSR",
        missed_cleavages = invalid_missed_cleavages[[value_name]]
      ),
      class = "pepvet_error_invalid_missed_cleavages",
      info = value_name
    )
  }
  expect_no_warning(
    expect_error(
      digest_protein(
        "MKWVTFISLLFLFSSAYSR",
        missed_cleavages = 1e20
      ),
      class = "pepvet_error_invalid_missed_cleavages"
    )
  )

  invalid_include_values <- list(
    null = NULL,
    missing = NA,
    wrong_type = 1,
    multiple = c(FALSE, TRUE),
    character = "yes"
  )
  for (value_name in names(invalid_include_values)) {
    expect_error(
      digest_protein(
        "MKWVTFISLLFLFSSAYSR",
        include_cleavage_efficiency = invalid_include_values[[value_name]]
      ),
      class = "pepvet_error_invalid_include_cleavage_efficiency",
      info = value_name
    )
  }
})

test_that("unnamed character inputs receive stable generated protein ids", {
  result <- digest_protein(
    c("MKWVTFISLLFLFSSAYSR", "AKRTPK"),
    enzyme = "trypsin"
  )

  expect_identical(unique(result$protein_id), c("sequence_1", "sequence_2"))
})
