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

.water_monoisotopic_mass <- 18.01056
.proton_monoisotopic_mass <- 1.007276
.terminal_pka <- c(n_term = 8.0, c_term = 3.1)
.ionizable_side_chain_pka <- c(
  C = 8.3,
  D = 3.9,
  E = 4.3,
  H = 6.0,
  K = 10.5,
  R = 12.5,
  Y = 10.1
)

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

.default_scoring_weights <- list(
  protein_only = c(
    S_length = 0.25,
    S_coverage = 0.25,
    S_count = 0.20,
    S_hydro = 0.15,
    S_charge = 0.15
  ),
  proteome_aware = c(
    S_length = 0.20,
    S_coverage = 0.20,
    S_count = 0.15,
    S_hydro = 0.15,
    S_charge = 0.10,
    S_unique = 0.20
  )
)

.preset_scoring_weights <- list(
  standard = c(
    S_length = 0.25,
    S_coverage = 0.25,
    S_count = 0.20,
    S_hydro = 0.15,
    S_charge = 0.15,
    S_unique = 0.00
  ),
  dia = c(
    S_length = 0.20,
    S_coverage = 0.30,
    S_count = 0.20,
    S_hydro = 0.10,
    S_charge = 0.10,
    S_unique = 0.10
  ),
  targeted = c(
    S_length = 0.15,
    S_coverage = 0.10,
    S_count = 0.15,
    S_hydro = 0.15,
    S_charge = 0.15,
    S_unique = 0.30
  ),
  membrane = c(
    S_length = 0.25,
    S_coverage = 0.25,
    S_count = 0.20,
    S_hydro = 0.05,
    S_charge = 0.15,
    S_unique = 0.10
  ),
  ffpe_degraded = c(
    S_length = 0.20,
    S_coverage = 0.20,
    S_count = 0.30,
    S_hydro = 0.10,
    S_charge = 0.10,
    S_unique = 0.10
  ),
  fractionated = c(
    S_length = 0.25,
    S_coverage = 0.25,
    S_count = 0.20,
    S_hydro = 0.15,
    S_charge = 0.15,
    S_unique = 0.00
  )
)

.pepvet_presets <- list(
  standard = list(
    gravy_range = c(-1.0, 0.6),
    length_range = c(7L, 25L),
    weights = .preset_scoring_weights$standard,
    include_pI = FALSE
  ),
  dia = list(
    gravy_range = c(-1.0, 0.6),
    length_range = c(7L, 30L),
    weights = .preset_scoring_weights$dia,
    include_pI = FALSE
  ),
  targeted = list(
    gravy_range = c(-0.8, 0.4),
    length_range = c(8L, 20L),
    weights = .preset_scoring_weights$targeted,
    include_pI = FALSE
  ),
  membrane = list(
    gravy_range = c(-1.0, 1.5),
    length_range = c(7L, 30L),
    weights = .preset_scoring_weights$membrane,
    include_pI = FALSE
  ),
  ffpe_degraded = list(
    gravy_range = c(-1.0, 0.8),
    length_range = c(6L, 30L),
    weights = .preset_scoring_weights$ffpe_degraded,
    include_pI = FALSE
  ),
  fractionated = list(
    gravy_range = c(-1.0, 0.6),
    length_range = c(7L, 25L),
    weights = .preset_scoring_weights$fractionated,
    include_pI = TRUE
  )
)

.validate_gravy_range <- function(gravy_range) {
  if (!is.numeric(gravy_range) || length(gravy_range) != 2L || anyNA(gravy_range)) {
    cli::cli_abort(
      "{.arg gravy_range} must be a numeric vector of length 2 with no missing values.",
      class = "pepvet_error_invalid_gravy_range"
    )
  }

  normalized_range <- as.numeric(gravy_range)

  if (!all(is.finite(normalized_range)) || normalized_range[[1]] > normalized_range[[2]]) {
    cli::cli_abort(
      "{.arg gravy_range} must contain finite values in ascending order.",
      class = "pepvet_error_invalid_gravy_range"
    )
  }

  normalized_range
}

.validate_length_range <- function(length_range) {
  if (!is.numeric(length_range) || length(length_range) != 2L || anyNA(length_range)) {
    cli::cli_abort(
      "{.arg length_range} must be a numeric vector of length 2 with no missing values.",
      class = "pepvet_error_invalid_length_range"
    )
  }

  normalized_range <- as.integer(length_range)

  if (
    any(normalized_range < 1L) ||
      normalized_range[[1]] > normalized_range[[2]] ||
      !isTRUE(all.equal(as.numeric(length_range), as.numeric(normalized_range)))
  ) {
    cli::cli_abort(
      "{.arg length_range} must contain positive integers in ascending order.",
      class = "pepvet_error_invalid_length_range"
    )
  }

  normalized_range
}

.normalize_weights <- function(weights, defaults) {
  if (!is.numeric(weights) || anyNA(weights)) {
    cli::cli_abort(
      "{.arg weights} must be a numeric vector with no missing values.",
      class = "pepvet_error_invalid_weights"
    )
  }

  expected_names <- names(defaults)

  if (is.null(names(weights))) {
    normalized_weights <- as.numeric(weights)
    names(normalized_weights) <- expected_names
    return(normalized_weights)
  }

  observed_names <- trimws(names(weights))

  if (anyNA(observed_names) || any(!nzchar(observed_names))) {
    cli::cli_abort(
      "Named {.arg weights} entries must all have non-empty names.",
      class = "pepvet_error_invalid_weights"
    )
  }

  if (!setequal(observed_names, expected_names)) {
    cli::cli_abort(
      c(
        "Named {.arg weights} must match the scoring component names.",
        "i" = paste("Expected names:", paste(expected_names, collapse = ", "))
      ),
      class = "pepvet_error_invalid_weights"
    )
  }

  normalized_weights <- as.numeric(weights[match(expected_names, observed_names)])
  names(normalized_weights) <- expected_names
  normalized_weights
}

.validate_weights <- function(weights, has_proteome) {
  defaults <- if (isTRUE(has_proteome)) {
    .default_scoring_weights$proteome_aware
  } else {
    .default_scoring_weights$protein_only
  }

  if (is.null(weights)) {
    return(defaults)
  }

  expected_lengths <- if (isTRUE(has_proteome)) c(6L) else c(5L, 6L)
  if (!length(weights) %in% expected_lengths) {
    cli::cli_abort(
      paste0(
        "{.arg weights} must contain exactly ",
        paste(expected_lengths, collapse = " or "),
        " value(s) in this scoring mode."
      ),
      class = "pepvet_error_invalid_weights"
    )
  }

  if (!isTRUE(has_proteome) && length(weights) == 6L) {
    normalized_weights <- .normalize_weights(weights, .default_scoring_weights$proteome_aware)

    if (normalized_weights[["S_unique"]] > 0) {
      cli::cli_abort(
        c(
          "{.arg weights} assigns a non-zero value to {.field S_unique} but no {.arg proteome} was supplied.",
          "i" = "Provide a proteome digest for uniqueness scoring or set S_unique to 0."
        ),
        class = "pepvet_error_invalid_weights"
      )
    }

    normalized_weights <- normalized_weights[names(defaults)]
  } else {
    normalized_weights <- .normalize_weights(weights, defaults)
  }

  if (any(normalized_weights < 0)) {
    cli::cli_abort(
      "{.arg weights} must not contain negative values.",
      class = "pepvet_error_invalid_weights"
    )
  }

  if (!isTRUE(all.equal(sum(normalized_weights), 1, tolerance = 1e-8))) {
    cli::cli_abort(
      "{.arg weights} must sum to 1.",
      class = "pepvet_error_invalid_weights"
    )
  }

  normalized_weights
}

#' Return a Named Scoring Preset
#'
#' `pepvet_preset()` returns a named list containing a GRAVY range, peptide
#' length range, and scoring weights for a supported proteomics workflow.
#' Presets are intended as editable starting points rather than hard rules.
#'
#' Presets with non-zero `S_unique` weights require a comparison proteome at
#' scoring time so uniqueness can be measured honestly.
#'
#' @param type Preset name. Supported values are `"standard"`, `"dia"`,
#'   `"targeted"`, `"membrane"`, `"ffpe_degraded"`, and `"fractionated"`.
#'
#' @return A named list with `gravy_range`, `length_range`, `weights`, and
#'   `include_pI`.
#'   The returned object can be passed directly into [score_peptides()] or
#'   [evaluate_digest()] through `do.call()` or argument splicing.
#'
#' @examples
#' pepvet_preset("standard")
#' @export
pepvet_preset <- function(type = "standard") {
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    cli::cli_abort(
      "{.arg type} must be a single, non-missing character string.",
      class = "pepvet_error_invalid_preset"
    )
  }

  normalized_type <- tolower(trimws(type))

  if (!normalized_type %in% names(.pepvet_presets)) {
    cli::cli_abort(
      c(
        "{.arg type} must be one of pepVet's supported preset names.",
        "i" = paste("Supported presets:", paste(names(.pepvet_presets), collapse = ", "))
      ),
      class = "pepvet_error_invalid_preset"
    )
  }

  preset <- .pepvet_presets[[normalized_type]]
  preset$gravy_range <- .validate_gravy_range(preset$gravy_range)
  preset$length_range <- .validate_length_range(preset$length_range)
  preset$weights <- .normalize_weights(preset$weights, .default_scoring_weights$proteome_aware)
  preset$include_pI <- .validate_include_pI(preset$include_pI)
  preset
}

.same_numeric_values <- function(x, y, tolerance = 1e-8) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
}

.same_named_weights <- function(x, y, tolerance = 1e-8) {
  identical(names(x), names(y)) && .same_numeric_values(x, y, tolerance = tolerance)
}

.identify_preset_used <- function(gravy_range,
                                  length_range,
                                  weights,
                                  include_pI,
                                  has_proteome) {
  matched_presets <- vapply(
    names(.pepvet_presets),
    function(preset_name) {
      preset <- .pepvet_presets[[preset_name]]
      effective_weights <- tryCatch(
        .validate_weights(preset$weights, has_proteome),
        pepvet_error_invalid_weights = function(error) NULL
      )

      if (is.null(effective_weights)) {
        return(FALSE)
      }

      identical(.validate_length_range(preset$length_range), length_range) &&
        .same_numeric_values(.validate_gravy_range(preset$gravy_range), gravy_range) &&
        .same_named_weights(effective_weights, weights) &&
        identical(.validate_include_pI(preset$include_pI), include_pI)
    },
    logical(1)
  )

  matches <- names(.pepvet_presets)[matched_presets]

  if (length(matches) == 1L) {
    return(matches)
  }

  "custom"
}

.build_proteome_index <- function(proteome_digests) {
  index <- new.env(hash = TRUE, parent = emptyenv())

  if (nrow(proteome_digests) == 0L) {
    return(index)
  }

  peptide_pairs <- unique(proteome_digests[c("peptide", "protein_id")])

  for (row_index in seq_len(nrow(peptide_pairs))) {
    peptide <- peptide_pairs$peptide[[row_index]]
    protein_id <- peptide_pairs$protein_id[[row_index]]
    indexed_proteins <- get0(peptide, envir = index, inherits = FALSE)

    if (is.null(indexed_proteins)) {
      assign(peptide, protein_id, envir = index)
    } else if (!protein_id %in% indexed_proteins) {
      assign(peptide, c(indexed_proteins, protein_id), envir = index)
    }
  }

  index
}

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

.normalize_peptide_sequences <- function(sequence, arg_name = "sequence") {
  if (!is.character(sequence) || length(sequence) == 0L) {
    cli::cli_abort(
      paste0(
        "{.arg ",
        arg_name,
        "} must be a non-empty character vector of peptide sequences."
      ),
      class = "pepvet_error_invalid_sequence"
    )
  }

  sequence_names <- names(sequence)

  vapply(
    seq_along(sequence),
    function(index) {
      sequence_name <- if (
        !is.null(sequence_names) &&
          !is.na(sequence_names[[index]]) &&
          nzchar(sequence_names[[index]])
      ) {
        sequence_names[[index]]
      } else {
        paste0(arg_name, "_", index)
      }

      .validate_sequence(sequence[[index]], sequence_name = sequence_name)
    },
    character(1)
  )
}

.validate_charge <- function(charge, sequence_count) {
  if (!is.numeric(charge) || anyNA(charge)) {
    cli::cli_abort(
      "{.arg charge} must be a numeric vector of non-missing integers.",
      class = "pepvet_error_invalid_charge"
    )
  }

  if (!length(charge) %in% c(1L, sequence_count)) {
    cli::cli_abort(
      "{.arg charge} must have length 1 or the same length as {.arg sequence}.",
      class = "pepvet_error_invalid_charge"
    )
  }

  normalized_charge <- as.integer(charge)

  if (
    any(normalized_charge < 0L) ||
      !isTRUE(all.equal(as.numeric(charge), as.numeric(normalized_charge)))
  ) {
    cli::cli_abort(
      "{.arg charge} must contain non-negative integers.",
      class = "pepvet_error_invalid_charge"
    )
  }

  if (length(normalized_charge) == 1L) {
    normalized_charge <- rep(normalized_charge, sequence_count)
  }

  normalized_charge
}

.validate_include_pI <- function(include_pI) {
  if (!is.logical(include_pI) || length(include_pI) != 1L || is.na(include_pI)) {
    cli::cli_abort(
      "{.arg include_pI} must be a single, non-missing logical value.",
      class = "pepvet_error_invalid_include_pi"
    )
  }

  include_pI
}

.ionizable_composition_matrix <- function(sequence) {
  normalized_sequences <- .normalize_peptide_sequences(sequence)
  ionizable_residues <- names(.ionizable_side_chain_pka)
  composition_matrix <- matrix(
    0,
    nrow = length(normalized_sequences),
    ncol = length(ionizable_residues),
    dimnames = list(NULL, ionizable_residues)
  )

  residue_indices <- lapply(
    strsplit(normalized_sequences, split = "", fixed = TRUE),
    function(residues) {
      match(residues, ionizable_residues)
    }
  )

  for (index in seq_along(residue_indices)) {
    matched_residues <- residue_indices[[index]]
    matched_residues <- matched_residues[!is.na(matched_residues)]

    if (length(matched_residues) > 0L) {
      composition_matrix[index, ] <- tabulate(
        matched_residues,
        nbins = length(ionizable_residues)
      )
    }
  }

  composition_matrix
}

.net_charge_at_pH <- function(pH, composition_matrix) {
  net_charge <-
    1 / (1 + 10^(pH - .terminal_pka[["n_term"]])) -
    1 / (1 + 10^(.terminal_pka[["c_term"]] - pH))

  for (residue in c("H", "K", "R")) {
    net_charge <- net_charge +
      composition_matrix[, residue] / (1 + 10^(pH - .ionizable_side_chain_pka[[residue]]))
  }

  for (residue in c("C", "D", "E", "Y")) {
    net_charge <- net_charge -
      composition_matrix[, residue] / (1 + 10^(.ionizable_side_chain_pka[[residue]] - pH))
  }

  net_charge
}

#' Calculate Peptide Mass or m/z
#'
#' `calculate_peptide_mass()` computes the neutral monoisotopic mass of one or
#' more unmodified peptide sequences using residue masses stored in
#' [aa_properties]. When `charge > 0`, the function returns monoisotopic m/z
#' values using the proton mass `1.007276` Da.
#'
#' @param sequence Peptide sequence supplied as a character vector. Amino-acid
#'   codes must be one-letter uppercase or lowercase symbols from the 20
#'   standard residues.
#' @param charge Optional non-negative integer scalar or integer vector with
#'   length matching `sequence`. `0L` returns neutral masses. Positive values
#'   return m/z.
#'
#' @return A numeric vector of neutral masses or m/z values. Names are
#'   preserved from `sequence` when present.
#'
#' @details This function computes masses for unmodified peptide sequences only.
#'   It does not account for chemical labels such as TMT or iTRAQ, isotopic
#'   labels such as SILAC, or post-translational modifications.
#'
#' @examples
#' calculate_peptide_mass("PEPTIDE")
#' calculate_peptide_mass(c(a = "PEPTIDE", b = "AAAAAAAR"), charge = 2L)
#' @export
calculate_peptide_mass <- function(sequence, charge = 0L) {
  normalized_sequences <- .normalize_peptide_sequences(sequence)
  normalized_charge <- .validate_charge(charge, length(normalized_sequences))
  aa_properties <- .get_aa_properties()
  residue_masses <- stats::setNames(
    aa_properties$residue_monoisotopic_mass,
    aa_properties$amino_acid
  )

  neutral_mass <- vapply(
    strsplit(normalized_sequences, split = "", fixed = TRUE),
    function(residues) {
      sum(residue_masses[residues]) + .water_monoisotopic_mass
    },
    numeric(1)
  )

  result <- neutral_mass
  charged <- normalized_charge > 0L

  if (any(charged)) {
    result[charged] <- (
      neutral_mass[charged] +
        normalized_charge[charged] * .proton_monoisotopic_mass
    ) / normalized_charge[charged]
  }

  names(result) <- names(sequence)
  result
}

#' Calculate Peptide Isoelectric Point
#'
#' `calculate_pI()` estimates peptide isoelectric points with a bisection
#' search over the Henderson-Hasselbalch net-charge equation. The calculation
#' uses hard-coded terminal pKa values of `8.0` for the N-terminus and `3.1`
#' for the C-terminus, plus side-chain pKa values from [aa_properties] for `C`,
#' `D`, `E`, `H`, `K`, `R`, and `Y`.
#'
#' @param sequence Peptide sequence supplied as a character vector. Amino-acid
#'   codes must be one-letter uppercase or lowercase symbols from the 20
#'   standard residues.
#'
#' @return A numeric vector of estimated peptide pI values. Names are preserved
#'   from `sequence` when present.
#'
#' @examples
#' calculate_pI("PEPTIDE")
#' calculate_pI(c("AAAAAAAR", "ACDEFGHIKLMNPQRSTVWY"))
#' @export
calculate_pI <- function(sequence) {
  normalized_sequences <- .normalize_peptide_sequences(sequence)

  if (length(normalized_sequences) > 5000L) {
    cli::cli_inform(
      paste0(
        "Calculating peptide pI values for ",
        length(normalized_sequences),
        " sequences."
      )
    )
  }

  composition_matrix <- .ionizable_composition_matrix(normalized_sequences)
  lower_bound <- rep(0, length(normalized_sequences))
  upper_bound <- rep(14, length(normalized_sequences))

  for (iteration in seq_len(40L)) {
    midpoint <- (lower_bound + upper_bound) / 2
    midpoint_charge <- .net_charge_at_pH(midpoint, composition_matrix)
    still_positive <- midpoint_charge > 0

    lower_bound[still_positive] <- midpoint[still_positive]
    upper_bound[!still_positive] <- midpoint[!still_positive]
  }

  result <- (lower_bound + upper_bound) / 2
  names(result) <- names(sequence)
  result
}
