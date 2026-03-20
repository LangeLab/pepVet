.get_aa_properties <- local({
  aa_properties_cache <- NULL

  function() {
    if (is.null(aa_properties_cache)) {
      data_env <- new.env(parent = emptyenv())
      utils::data("aa_properties", package = "pepVet", envir = data_env)
      aa_properties_cache <<- data_env$aa_properties
    }

    aa_properties_cache
  }
})

.supported_digest_enzymes <- c(
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

.normalize_enzyme <- function(enzyme) {
  if (!is.character(enzyme) || length(enzyme) != 1L || is.na(enzyme)) {
    cli::cli_abort(
      "{.arg enzyme} must be a single, non-missing character string.",
      class = "pepvet_error_invalid_enzyme"
    )
  }

  enzyme <- tolower(trimws(enzyme))

  if (!nzchar(enzyme)) {
    cli::cli_abort(
      "{.arg enzyme} must not be empty.",
      class = "pepvet_error_invalid_enzyme"
    )
  }

  if (!enzyme %in% .supported_digest_enzymes) {
    cli::cli_abort(
      c(
        paste(
          "{.arg enzyme} must be one of pepVet's supported",
          "cleaver-compatible enzyme names."
        ),
        "i" = paste(
          "Supported enzymes:",
          paste(.supported_digest_enzymes, collapse = ", ")
        )
      ),
      class = "pepvet_error_invalid_enzyme"
    )
  }

  enzyme
}

.validate_missed_cleavages <- function(missed_cleavages) {
  if (!is.numeric(missed_cleavages) || length(missed_cleavages) != 1L) {
    cli::cli_abort(
      "{.arg missed_cleavages} must be a single non-negative integer.",
      class = "pepvet_error_invalid_missed_cleavages"
    )
  }

  if (is.na(missed_cleavages) || missed_cleavages < 0) {
    cli::cli_abort(
      "{.arg missed_cleavages} must be a single non-negative integer.",
      class = "pepvet_error_invalid_missed_cleavages"
    )
  }

  missed_cleavages_int <- as.integer(missed_cleavages)

  if (!isTRUE(all.equal(missed_cleavages, missed_cleavages_int))) {
    cli::cli_abort(
      "{.arg missed_cleavages} must be a single non-negative integer.",
      class = "pepvet_error_invalid_missed_cleavages"
    )
  }

  missed_cleavages_int
}

.looks_like_path <- function(path) {
  grepl("[/\\\\]", path) || grepl("[.][A-Za-z0-9]+$", basename(path))
}

.normalize_sequence_names <- function(sequence_names, sequence_count) {
  normalized_names <- rep(NA_character_, sequence_count)

  if (!is.null(sequence_names)) {
    normalized_names[seq_along(sequence_names)] <- trimws(sequence_names)
  }

  missing_name <- is.na(normalized_names) | !nzchar(normalized_names)
  normalized_names[missing_name] <- paste0("sequence_", which(missing_name))
  normalized_names
}

.validate_sequence <- function(sequence, sequence_name = "sequence") {
  if (!is.character(sequence) || length(sequence) != 1L) {
    cli::cli_abort(
      "{.arg sequence} must be a single character string.",
      class = "pepvet_error_invalid_sequence"
    )
  }

  if (is.na(sequence)) {
    cli::cli_abort(
      paste0("Sequence '", sequence_name, "' must not be missing."),
      class = "pepvet_error_invalid_sequence"
    )
  }

  if (!nzchar(trimws(sequence))) {
    cli::cli_abort(
      paste0("Sequence '", sequence_name, "' must not be empty."),
      class = "pepvet_error_invalid_sequence"
    )
  }

  sequence <- toupper(sequence)
  residues <- strsplit(sequence, split = "", fixed = TRUE)[[1]]
  allowed_residues <- .get_aa_properties()$amino_acid
  invalid_residues <- unique(residues[!residues %in% allowed_residues])

  if (length(invalid_residues) > 0L) {
    cli::cli_abort(
      paste0(
        "Sequence '",
        sequence_name,
        "' contains unsupported amino acid code(s): ",
        paste(invalid_residues, collapse = ", "),
        "."
      ),
      class = "pepvet_error_invalid_sequence"
    )
  }

  sequence
}

.read_input <- function(sequence) {
  input_error <- paste(
    "{.arg sequence} must be a character sequence, named character vector,",
    paste(
      "{.cls Biostrings::AAString}, {.cls Biostrings::AAStringSet},",
      "or a FASTA file path."
    )
  )

  if (is.null(sequence)) {
    cli::cli_abort(input_error, class = "pepvet_error_invalid_input")
  }

  if (inherits(sequence, "AAStringSet")) {
    raw_sequences <- as.character(sequence)
    raw_names <- names(sequence)
  } else if (inherits(sequence, "AAString")) {
    raw_sequences <- as.character(sequence)
    raw_names <- names(Biostrings::AAStringSet(sequence))
  } else if (is.character(sequence)) {
    if (length(sequence) == 0L) {
      cli::cli_abort(input_error, class = "pepvet_error_invalid_input")
    }

    if (length(sequence) == 1L && !is.na(sequence) && file.exists(sequence)) {
      if (dir.exists(sequence)) {
        cli::cli_abort(
          paste0(
            "Expected a FASTA file path, but '",
            sequence,
            "' is a directory."
          ),
          class = "pepvet_error_missing_file"
        )
      }

      fasta_sequences <- Biostrings::readAAStringSet(sequence, format = "fasta")
      raw_sequences <- as.character(fasta_sequences)
      raw_names <- names(fasta_sequences)
    } else {
      if (
        length(sequence) == 1L &&
          !is.na(sequence) &&
          .looks_like_path(sequence)
      ) {
        cli::cli_abort(
          paste0("FASTA file not found: '", sequence, "'."),
          class = "pepvet_error_missing_file"
        )
      }

      raw_sequences <- unname(sequence)
      raw_names <- names(sequence)
    }
  } else {
    cli::cli_abort(input_error, class = "pepvet_error_invalid_input")
  }

  if (length(raw_sequences) == 0L) {
    cli::cli_abort(
      "{.arg sequence} must contain at least one sequence entry.",
      class = "pepvet_error_invalid_input"
    )
  }

  normalized_names <- .normalize_sequence_names(
    raw_names,
    length(raw_sequences)
  )
  normalized_sequences <- vapply(
    seq_along(raw_sequences),
    function(index) {
      .validate_sequence(raw_sequences[[index]], normalized_names[[index]])
    },
    character(1)
  )

  aa_string_set <- Biostrings::AAStringSet(normalized_sequences)
  names(aa_string_set) <- normalized_names
  aa_string_set
}

.build_digest_ranges <- function(strict_ranges, missed_cleavages) {
  range_starts <- IRanges::start(strict_ranges)
  range_ends <- IRanges::end(strict_ranges)
  peptide_count <- length(strict_ranges)
  digest_rows <- vector("list", peptide_count * (missed_cleavages + 1L))
  row_index <- 1L

  for (start_index in seq_len(peptide_count)) {
    max_missed <- min(missed_cleavages, peptide_count - start_index)

    for (missed in 0:max_missed) {
      end_index <- start_index + missed
      digest_rows[[row_index]] <- list(
        start = range_starts[[start_index]],
        end = range_ends[[end_index]],
        missed_cleavages = missed
      )
      row_index <- row_index + 1L
    }
  }

  digest_rows[seq_len(row_index - 1L)]
}

.calculate_gravy <- function(peptide_sequence) {
  if (!is.character(peptide_sequence) || length(peptide_sequence) != 1L) {
    cli::cli_abort("{.arg peptide_sequence} must be a single character string.")
  }

  if (is.na(peptide_sequence)) {
    cli::cli_abort("{.arg peptide_sequence} must not be missing.")
  }

  if (!nzchar(peptide_sequence)) {
    cli::cli_abort("{.arg peptide_sequence} must not be empty.")
  }

  peptide_sequence <- toupper(peptide_sequence)
  aa_properties <- .get_aa_properties()

  residues <- strsplit(peptide_sequence, split = "", fixed = TRUE)[[1]]
  residue_index <- match(residues, aa_properties$amino_acid)

  if (anyNA(residue_index)) {
    cli::cli_abort(
      paste0(
        "Unknown amino acid code(s): ",
        paste(unique(residues[is.na(residue_index)]), collapse = ", "),
        "."
      )
    )
  }

  mean(aa_properties$hydrophobicity[residue_index])
}
