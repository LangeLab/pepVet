#' @importFrom tibble tibble as_tibble add_column
#' @importFrom Biostrings AAString AAStringSet readAAStringSet
#' @importFrom IRanges start end
#' @importFrom rlang check_installed caller_env
#' @importFrom cli cli_abort cli_warn cli_inform cli_text
#' @importFrom cli cat_line cat_rule symbol
#' @importFrom cli col_blue col_green col_red col_silver col_yellow
#' @importFrom cli style_bold style_italic
#' @importFrom stats rgamma
NULL

.pepvet_cache <- new.env(parent = emptyenv())

.get_aa_properties <- function() {
  if (is.null(.pepvet_cache$aa_properties)) {
    data_env <- new.env(parent = emptyenv())
    utils::data("aa_properties", package = "pepVet", envir = data_env)
    .pepvet_cache$aa_properties <- data_env$aa_properties
  }

  .pepvet_cache$aa_properties
}

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
  Y = 10.1,
  U = 5.2
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
    S_length   = 0.200,
    S_coverage = 0.348,
    S_count    = 0.226,
    S_hydro    = 0.138,
    S_charge   = 0.088
  ),
  proteome_aware = c(
    S_length   = 0.160,
    S_coverage = 0.279,
    S_count    = 0.181,
    S_hydro    = 0.110,
    S_charge   = 0.070,
    S_unique   = 0.200
  )
)

.preset_scoring_weights <- list(
  standard = c(
    S_length   = 0.200,
    S_coverage = 0.348,
    S_count    = 0.226,
    S_hydro    = 0.138,
    S_charge   = 0.088,
    S_unique   = 0.000
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
    S_length   = 0.200,
    S_coverage = 0.348,
    S_count    = 0.226,
    S_hydro    = 0.138,
    S_charge   = 0.088,
    S_unique   = 0.000
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
    gravy_range = c(-1.0, 0.8),
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
    gravy_range = c(-1.0, 2.0),
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

.pepvet_params <- list(
  verdict_good         = 0.65,
  verdict_moderate     = 0.40,
  length_lo            = 7L,
  length_hi            = 25L,
  gravy_lo             = -2.0,
  gravy_hi             =  2.0,
  scatter_alpha        = 0.88,
  scatter_max_pts      = 5000L,
  patchwork_title_size = 15,
  patchwork_tag_size   = 14,
  sensitivity_bw       = 0.025,
  pi_binwidth          = 0.25,
  length_binwidth      = 1L,
  batch_score_binwidth = 0.05
)

.get_param <- function(name) {
  .pepvet_params[[name]]
}

.get_plot_theme_value <- function(name, default) {
  override <- .pepvet_config_env$theme_overrides[[name]]
  if (is.null(override)) default else override
}

.resolve_plot_metadata_range <- function(result, value, name, default) {
  has_params <- is.list(result) && !is.data.frame(result) &&
    is.list(result$params) && !is.null(result$params[[name]])
  uses_default <- is.numeric(value) && isTRUE(all.equal(
    unname(as.numeric(value)), unname(as.numeric(default))
  ))
  if (has_params && uses_default) result$params[[name]] else value
}

.validate_cleavage_sites_for_plot <- function(cleavage_sites,
                                              protein_length = NULL) {
  if (is.null(cleavage_sites)) {
    return(invisible(NULL))
  }

  required <- c("position", "efficiency")
  if (!is.data.frame(cleavage_sites) ||
      !all(required %in% names(cleavage_sites))) {
    .abort(
      c(
        "{.arg cleavage_sites} must be a data.frame from {.fn annotate_cleavage_sites}.",
        "i" = "Required columns: {.code position}, {.code efficiency}."
      ),
      class = "pepvet_error_invalid_cleavage_sites"
    )
  }

  positions <- cleavage_sites$position
  efficiency_input <- cleavage_sites$efficiency
  valid_efficiency_type <- is.character(efficiency_input) ||
    is.factor(efficiency_input)
  efficiency <- if (valid_efficiency_type) {
    as.character(efficiency_input)
  } else {
    character(0L)
  }
  valid_positions <- is.numeric(positions) &&
    all(is.finite(positions)) && all(positions == floor(positions)) &&
    all(positions >= 1L)
  valid_efficiency <- valid_efficiency_type &&
    !anyNA(efficiency) && all(efficiency %in% c("high", "medium", "low"))

  if (!valid_positions || !valid_efficiency ||
      (!is.null(protein_length) && any(positions > protein_length))) {
    .abort(
      "{.arg cleavage_sites} contains invalid positions or efficiency levels.",
      class = "pepvet_error_invalid_cleavage_sites"
    )
  }

  invisible(cleavage_sites)
}

.validate_domains_for_plot <- function(domains, protein_length = NULL) {
  if (is.null(domains)) {
    return(invisible(NULL))
  }

  required <- c("name", "start", "end")
  if (!is.data.frame(domains) || !all(required %in% names(domains))) {
    .abort(
      c(
        paste0(
          "{.arg domains} must be a data.frame with columns ",
          "{.code name}, {.code start}, and {.code end}."
        ),
        "i" = "Each row describes one annotated protein domain."
      ),
      class = "pepvet_error_invalid_domains"
    )
  }

  valid_names <- is.character(domains$name) && !anyNA(domains$name) &&
    all(nzchar(trimws(domains$name)))
  valid_coordinates <- all(vapply(
    domains[c("start", "end")],
    function(values) {
      is.numeric(values) && all(is.finite(values)) &&
        all(values == floor(values))
    },
    logical(1L)
  ))
  valid_order <- valid_coordinates &&
    all(domains$start >= 1L) && all(domains$end >= domains$start)
  in_protein <- is.null(protein_length) ||
    (valid_order && all(domains$end <= protein_length))

  if (!valid_names || !valid_order || !in_protein) {
    .abort(
      "{.arg domains} contains invalid names or residue coordinates.",
      class = "pepvet_error_invalid_domains"
    )
  }

  invisible(domains)
}

## Dirichlet random vector generator (base R, no dependencies)
## Returns an n x k matrix where each row sums to 1.
.rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  m <- matrix(rgamma(n * k, shape = rep(alpha, each = n)), nrow = n, ncol = k)
  totals <- rowSums(m)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    .abort(
      "Dirichlet concentration parameters produced non-finite weight draws.",
      class = "pepvet_error_invalid_sensitivity_parameter"
    )
  }
  m / totals
}

.abort <- function(message, ...,
                   class = NULL,
                   call = rlang::caller_env()) {
  cli::cli_abort(
    message, ...,
    class = c(class, "pepvet_error"),
    call = call,
    .envir = call
  )
}

.validate_unique_columns <- function(data, arg_name, class) {
  if (anyDuplicated(names(data)) > 0L) {
    .abort(
      "{.arg {arg_name}} must have unique column names.",
      class = class
    )
  }

  invisible(data)
}

.bind_rows <- function(df_list) {
  if (length(df_list) == 0L) return(tibble::tibble())
  if (length(df_list) == 1L) return(tibble::as_tibble(df_list[[1L]]))
  for (i in seq_along(df_list)) rownames(df_list[[i]]) <- NULL
  tibble::as_tibble(do.call(rbind, df_list))
}

.validate_gravy_range <- function(gravy_range) {
  if (!is.numeric(gravy_range) ||
      length(gravy_range) != 2L || anyNA(gravy_range)) {
    .abort(
      paste0(
        "{.arg gravy_range} must be a numeric vector ",
        "of length 2 with no missing values."
      ),
      class = "pepvet_error_invalid_gravy_range"
    )
  }

  normalized_range <- as.numeric(gravy_range)

  if (!all(is.finite(normalized_range)) ||
      normalized_range[[1]] > normalized_range[[2]]) {
    .abort(
      c(
        paste0(
          "{.arg gravy_range} must contain finite ",
          "values in non-decreasing order."
        ),
        "i" = "Use c(lower, upper) where lower <= upper, e.g. c(-1.0, 0.6)."
      ),
      class = "pepvet_error_invalid_gravy_range"
    )
  }

  normalized_range
}

.validate_length_range <- function(length_range) {
  if (!is.numeric(length_range) ||
      length(length_range) != 2L || anyNA(length_range)) {
    .abort(
      paste0(
        "{.arg length_range} must be a numeric vector ",
        "of length 2 with no missing values."
      ),
      class = "pepvet_error_invalid_length_range"
    )
  }

  if (!all(is.finite(length_range))) {
    .abort(
      paste0(
        "{.arg length_range} must contain finite positive integers in ",
        "non-decreasing order."
      ),
      class = "pepvet_error_invalid_length_range"
    )
  }

  if (
    any(length_range < 1) ||
      any(length_range > .Machine$integer.max) ||
      length_range[[1]] > length_range[[2]] ||
      !isTRUE(all.equal(as.numeric(length_range), floor(length_range)))
  ) {
    .abort(
      paste0(
        "{.arg length_range} must contain positive integers in ",
        "non-decreasing order."
      ),
      class = "pepvet_error_invalid_length_range"
    )
  }

  as.integer(length_range)
}

.normalize_weights <- function(weights, defaults) {
  if (!is.numeric(weights) || anyNA(weights) || !all(is.finite(weights))) {
    .abort(
      "{.arg weights} must be a finite numeric vector with no missing values.",
      class = "pepvet_error_invalid_weights"
    )
  }

  if (length(weights) != length(defaults)) {
    .abort(
      "{.arg weights} must have the expected number of scoring components.",
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
    .abort(
      "Named {.arg weights} entries must all have non-empty names.",
      class = "pepvet_error_invalid_weights"
    )
  }

  if (anyDuplicated(observed_names) > 0L) {
    .abort(
      "Named {.arg weights} must have unique names.",
      class = "pepvet_error_invalid_weights"
    )
  }

  if (!setequal(observed_names, expected_names)) {
    .abort(
      c(
        "Named {.arg weights} must match the scoring component names.",
        "i" = "Expected names: {.val {expected_names}}"
      ),
      class = "pepvet_error_invalid_weights"
    )
  }

  normalized_weights <- as.numeric(
    weights[match(expected_names, observed_names)]
  )
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
    .abort(
      c(
        paste0(
          "{.arg weights} must contain exactly ",
          paste(expected_lengths, collapse = " or "),
          " value(s) in this scoring mode."
        ),
        "i" = "Pass 5 values for protein-only mode or 6 for proteome-aware mode."
      ),
      class = "pepvet_error_invalid_weights"
    )
  }

  if (!isTRUE(has_proteome) && length(weights) == 6L) {
    normalized_weights <- .normalize_weights(
      weights, .default_scoring_weights$proteome_aware
    )

    if (normalized_weights[["S_unique"]] > 0) {
      .abort(
        c(
          paste0(
            "{.arg weights} assigns a non-zero value to {.field S_unique} ",
            "but no {.arg proteome} was supplied."
          ),
          "i" = paste0(
            "Provide a proteome digest for uniqueness ",
            "scoring or set S_unique to 0."
          )
        ),
        class = "pepvet_error_invalid_weights"
      )
    }

    normalized_weights <- normalized_weights[names(defaults)]
  } else {
    normalized_weights <- .normalize_weights(weights, defaults)
  }

  if (any(normalized_weights < 0)) {
    .abort(
      "{.arg weights} must not contain negative values.",
      class = "pepvet_error_invalid_weights"
    )
  }

  if (!isTRUE(all.equal(sum(normalized_weights), 1, tolerance = 1e-8))) {
    .abort(
      "{.arg weights} must sum to 1.",
      class = "pepvet_error_invalid_weights"
    )
  }

  normalized_weights
}

#' Return a named scoring preset
#'
#' `pepvet_preset()` returns a named list containing a GRAVY range, peptide
#' length range, and scoring weights for a supported proteomics workflow.
#' Presets are intended as editable starting points rather than hard rules.
#' Their exact ranges and weights are conservative package choices, not
#' empirically calibrated boundaries.
#'
#' @param type Preset name. Defaults to `"standard"`. Supported values are
#'   `"standard"`, `"dia"`, `"targeted"`, `"membrane"`, `"ffpe_degraded"`, and
#'   `"fractionated"`. If `NULL`, raises an error.
#'
#' @section Presets:
#'
#' **`"standard"`** : Default scoring configuration for routine DDA examples.
#' Uses the `[7, 25]` length range, `[-1.0, 0.6]` GRAVY range, and default
#' protein-only weights.
#'
#' **`"dia"`** : A wider starting configuration for DIA or SWATH examples.
#' Uses the `[7, 30]` length range, `[-1.0, 0.8]` GRAVY range, and a larger
#' coverage weight than the standard preset.
#'
#' **`"targeted"`** : A narrower starting configuration for SRM, PRM, or MRM
#' examples. Uses the `[8, 20]` length range and `[-0.8, 0.4]` GRAVY range.
#' When a background proteome digest is supplied, uniqueness receives 30 percent
#' of the composite weight.
#'
#' **`"membrane"`** : A wider hydrophobicity configuration for membrane-protein
#' review. Uses the `[7, 30]` length range, extends the upper GRAVY boundary to
#' `2.0`, and assigns 5 percent of the composite weight to \code{S_hydro}.
#'
#' **`"ffpe_degraded"`** : A broader length configuration for exploratory work
#' with degraded or FFPE-derived material. Uses the `[6, 30]` length range and
#' assigns more weight to \code{S_count} than the standard preset.
#'
#' **`"fractionated"`** : SCX / high-pH RP fractionation planning. Same
#' scoring parameters as \code{"standard"} but with \code{include_pI = TRUE}
#' to append peptide-level pI values for fractionation-aware analysis.
#'
#' @family utils
#' @section Limitations:
#'   The six presets encode package priors. Inspect the returned values and use
#'   explicit arguments when the experimental context calls for other choices.
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
    .abort(
      "{.arg type} must be a single, non-missing character string.",
      class = "pepvet_error_invalid_preset"
    )
  }

  normalized_type <- tolower(trimws(type))

  if (!normalized_type %in% names(.pepvet_presets)) {
    .abort(
      c(
        "{.arg type} must be one of pepVet's supported preset names.",
        "i" = "Supported presets: {.val {names(.pepvet_presets)}}"
      ),
      class = "pepvet_error_invalid_preset"
    )
  }

  preset <- .pepvet_presets[[normalized_type]]
  preset$gravy_range <- .validate_gravy_range(preset$gravy_range)
  preset$length_range <- .validate_length_range(preset$length_range)
  preset$weights <- .normalize_weights(
    preset$weights, .default_scoring_weights$proteome_aware
  )
  preset$include_pI <- .validate_include_pI(preset$include_pI)
  preset
}

.same_numeric_values <- function(x, y, tolerance = 1e-8) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
}

.same_named_weights <- function(x, y, tolerance = 1e-8) {
  identical(names(x), names(y)) &&
    .same_numeric_values(x, y, tolerance = tolerance)
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
        .same_numeric_values(
          .validate_gravy_range(preset$gravy_range), gravy_range
        ) &&
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

  ## unique() removes duplicate (peptide, protein_id) pairs before grouping,
  ## so each list element already contains the deduplicated protein_id set.
  peptide_pairs <- unique(proteome_digests[c("peptide", "protein_id")])
  grouped <- split(peptide_pairs$protein_id, peptide_pairs$peptide)
  list2env(grouped, envir = index)
  index
}

.normalize_enzyme <- function(enzyme) {
  if (!is.character(enzyme) || length(enzyme) != 1L || is.na(enzyme)) {
    .abort(
      "{.arg enzyme} must be a single, non-missing character string.",
      class = "pepvet_error_invalid_enzyme"
    )
  }

  enzyme <- tolower(trimws(enzyme))

  if (!nzchar(enzyme)) {
    .abort(
      "{.arg enzyme} must not be empty.",
      class = "pepvet_error_invalid_enzyme"
    )
  }

  if (!enzyme %in% .supported_digest_enzymes) {
    .abort(
      c(
        paste0(
          "{.arg enzyme} must be one of pepVet's supported ",
          "cleaver-compatible enzyme names."
        ),
        "i" = "Supported enzymes: {.val { .supported_digest_enzymes}}"
      ),
      class = "pepvet_error_invalid_enzyme"
    )
  }

  enzyme
}

.validate_missed_cleavages <- function(missed_cleavages) {
  if (!is.numeric(missed_cleavages) || length(missed_cleavages) != 1L) {
    .abort(
      "{.arg missed_cleavages} must be a single non-negative integer.",
      class = "pepvet_error_invalid_missed_cleavages"
    )
  }

  if (
    is.na(missed_cleavages) ||
      !is.finite(missed_cleavages) ||
      missed_cleavages < 0
  ) {
    .abort(
      "{.arg missed_cleavages} must be a single non-negative integer.",
      class = "pepvet_error_invalid_missed_cleavages"
    )
  }

  if (missed_cleavages > .Machine$integer.max ||
      !isTRUE(all.equal(missed_cleavages, floor(missed_cleavages)))) {
    .abort(
      "{.arg missed_cleavages} must be a single non-negative integer.",
      class = "pepvet_error_invalid_missed_cleavages"
    )
  }

  as.integer(missed_cleavages)
}

.validate_include_cleavage_efficiency <- function(include_cleavage_efficiency) {
  if (
    !is.logical(include_cleavage_efficiency) ||
      length(include_cleavage_efficiency) != 1L ||
      is.na(include_cleavage_efficiency)
  ) {
    .abort(
      paste0(
        "{.arg include_cleavage_efficiency} must be a single, ",
        "non-missing logical value."
      ),
      class = "pepvet_error_invalid_include_cleavage_efficiency"
    )
  }

  include_cleavage_efficiency
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
  if (is.na(sequence)) {
    .abort(
      "Sequence '{sequence_name}' must not be missing.",
      class = "pepvet_error_invalid_sequence"
    )
  }

  if (!nzchar(trimws(sequence))) {
    .abort(
      "Sequence '{sequence_name}' must not be empty.",
      class = "pepvet_error_invalid_sequence"
    )
  }

  sequence <- toupper(sequence)
  residues <- strsplit(sequence, split = "", fixed = TRUE)[[1]]
  allowed_residues <- .get_aa_properties()$amino_acid
  invalid_residues <- unique(residues[!residues %in% allowed_residues])

  if (length(invalid_residues) > 0L) {
    .abort(
      paste0(
        "Sequence {.val {sequence_name}} contains unsupported ",
        "amino acid code(s): {.val {invalid_residues}}."
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
    .abort(input_error, class = "pepvet_error_invalid_input")
  }

  if (inherits(sequence, "AAStringSet")) {
    raw_sequences <- as.character(sequence)
    raw_names <- names(sequence)
  } else if (inherits(sequence, "AAString")) {
    raw_sequences <- as.character(sequence)
    raw_names <- names(Biostrings::AAStringSet(sequence))
  } else if (is.character(sequence)) {
    if (length(sequence) == 0L) {
      .abort(input_error, class = "pepvet_error_invalid_input")
    }

    directory_as_sequence <- FALSE
    if (length(sequence) == 1L && !is.na(sequence) && dir.exists(sequence)) {
      residues <- strsplit(toupper(sequence), split = "", fixed = TRUE)[[1]]
      directory_as_sequence <- nzchar(trimws(sequence)) &&
        all(residues %in% .get_aa_properties()$amino_acid)
    }

    if (
      length(sequence) == 1L &&
        !is.na(sequence) &&
        file.exists(sequence) &&
        !directory_as_sequence
    ) {
      if (dir.exists(sequence)) {
        .abort(
          c(
            "Expected a FASTA file path, but '{sequence}' is a directory.",
            "i" = "Provide a path to a FASTA file, not a directory."
          ),
          class = "pepvet_error_missing_file"
        )
      }

      fasta_sequences <- tryCatch(
        Biostrings::readAAStringSet(sequence, format = "fasta"),
        error = function(error) {
          .abort(
            c(
              "Could not parse {.arg sequence} as a FASTA file.",
              "i" = paste0(
                "Check the FASTA header and sequence records. Parser detail: ",
                conditionMessage(error)
              )
            ),
            class = "pepvet_error_invalid_input"
          )
        }
      )
      raw_sequences <- as.character(fasta_sequences)
      raw_names <- names(fasta_sequences)
    } else {
      if (
        length(sequence) == 1L &&
          !is.na(sequence) &&
          .looks_like_path(sequence)
      ) {
        .abort(
          c(
            "FASTA file not found: '{sequence}'.",
            "i" = "Check that the file path exists and is readable."
          ),
          class = "pepvet_error_missing_file"
        )
      }

      raw_sequences <- unname(sequence)
      raw_names <- names(sequence)
    }
  } else {
    .abort(input_error, class = "pepvet_error_invalid_input")
  }

  if (length(raw_sequences) == 0L) {
    .abort(
      "{.arg sequence} must contain at least one sequence entry.",
      class = "pepvet_error_invalid_input"
    )
  }

  normalized_names <- .normalize_sequence_names(
    raw_names,
    length(raw_sequences)
  )
  if (anyDuplicated(normalized_names) > 0L) {
    .abort(
      "{.arg sequence} must have unique protein identifiers.",
      class = "pepvet_error_invalid_input"
    )
  }
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

  ## Fast path: no missed cleavages, return vectors directly with no loop.
  if (missed_cleavages == 0L) {
    return(list(
      start            = range_starts,
      end              = range_ends,
      missed_cleavages = integer(peptide_count)
    ))
  }

  ## Pre-allocate output vectors (avoids repeated list reallocation).
  total_rows <- sum(vapply(
    seq_len(peptide_count),
    function(i) min(missed_cleavages, peptide_count - i) + 1L,
    integer(1)
  ))
  out_starts <- integer(total_rows)
  out_ends <- integer(total_rows)
  out_mc <- integer(total_rows)
  k <- 1L

  for (si in seq_len(peptide_count)) {
    max_mc <- min(missed_cleavages, peptide_count - si)
    for (mc in 0L:max_mc) {
      out_starts[k] <- range_starts[[si]]
      out_ends[k] <- range_ends[[si + mc]]
      out_mc[k] <- mc
      k <- k + 1L
    }
  }

  list(start = out_starts, end = out_ends, missed_cleavages = out_mc)
}

## Cached hydrophobicity lookup table (AA -> Kyte-Doolittle score).
## Built once from .get_aa_properties() and reused for all subsequent calls.
.get_hydro_lookup <- function() {
  if (is.null(.pepvet_cache$hydro_lookup)) {
    aa_props <- .get_aa_properties()
    .pepvet_cache$hydro_lookup <- stats::setNames(
      aa_props$hydrophobicity, aa_props$amino_acid
    )
  }
  .pepvet_cache$hydro_lookup
}

## Compute GRAVY scores for a character vector of peptide sequences.
## Caller is responsible for passing a validated, non-empty character vector of
## supported amino-acid sequences (as produced by digest_protein output). Case
## is normalized here, and a sequence with no known hydrophobicity returns NA.
## Mixed sequences average known residue values and ignore unknown values.
## Uses the cached hydrophobicity lookup; skips per-call validation.
.calculate_gravy <- function(peptide_vector) {
  hydro_lookup <- .get_hydro_lookup()
  res_lists <- strsplit(peptide_vector, "", fixed = TRUE)
  vapply(
    res_lists,
    function(res) {
      values <- hydro_lookup[toupper(res)]
      if (all(is.na(values))) {
        return(NA_real_)
      }
      mean(values, na.rm = TRUE)
    },
    numeric(1)
  )
}

.normalize_peptide_sequences <- function(sequence, arg_name = "sequence") {
  if (!is.character(sequence) || length(sequence) == 0L) {
    .abort(
      paste0(
        "{.arg {arg_name}} must be a non-empty character vector ",
        "of peptide sequences."
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
  if (!is.numeric(charge) || anyNA(charge) || !all(is.finite(charge))) {
    .abort(
      "{.arg charge} must be a numeric vector of non-missing integers.",
      class = "pepvet_error_invalid_charge"
    )
  }

  if (!length(charge) %in% c(1L, sequence_count)) {
    .abort(
      "{.arg charge} must have length 1 or the same length as {.arg sequence}.",
      class = "pepvet_error_invalid_charge"
    )
  }

  if (
    any(charge < 0) ||
      any(charge > .Machine$integer.max) ||
      !isTRUE(all.equal(as.numeric(charge), floor(charge)))
  ) {
    .abort(
      "{.arg charge} must contain non-negative integers.",
      class = "pepvet_error_invalid_charge"
    )
  }

  normalized_charge <- as.integer(charge)

  if (length(normalized_charge) == 1L) {
    normalized_charge <- rep(normalized_charge, sequence_count)
  }

  normalized_charge
}

.validate_include_pI <- function(include_pI) {
  if (!is.logical(include_pI) ||
      length(include_pI) != 1L || is.na(include_pI)) {
    .abort(
      "{.arg include_pI} must be a single, non-missing logical value.",
      class = "pepvet_error_invalid_include_pi"
    )
  }

  include_pI
}

## Internal: sequences already uppercase supported-amino-acid strings
## (no re-validation).
## Counts each of the 8 ionizable residues per peptide using vectorized gsub
## (8 gsub passes over the whole vector) instead of strsplit + a for-loop over
## every sequence.  This is >20x faster for large peptide sets.
.ionizable_composition_matrix_internal <- function(normalized_sequences) {
  ionizable_residues <- names(.ionizable_side_chain_pka)
  seq_lengths <- nchar(normalized_sequences)
  composition_matrix <- matrix(
    0L,
    nrow = length(normalized_sequences),
    ncol = length(ionizable_residues),
    dimnames = list(NULL, ionizable_residues)
  )

  for (j in seq_along(ionizable_residues)) {
    res <- ionizable_residues[[j]]
    composition_matrix[, j] <- seq_lengths -
      nchar(gsub(res, "", normalized_sequences, fixed = TRUE))
  }

  composition_matrix
}

.ionizable_composition_matrix <- function(sequence) {
  normalized_sequences <- .normalize_peptide_sequences(sequence)
  .ionizable_composition_matrix_internal(normalized_sequences)
}

.net_charge_at_pH <- function(pH, composition_matrix) {
  net_charge <-
    1 / (1 + 10^(pH - .terminal_pka[["n_term"]])) -
    1 / (1 + 10^(.terminal_pka[["c_term"]] - pH))

  for (residue in c("H", "K", "R")) {
    net_charge <- net_charge +
      composition_matrix[, residue] /
        (1 + 10^(pH - .ionizable_side_chain_pka[[residue]]))
  }

  for (residue in c("C", "D", "E", "Y", "U")) {
    net_charge <- net_charge -
      composition_matrix[, residue] /
        (1 + 10^(.ionizable_side_chain_pka[[residue]] - pH))
  }

  net_charge
}

#' Calculate peptide mass or m/z
#'
#' `calculate_peptide_mass()` computes the neutral monoisotopic mass of one or
#' more unmodified peptide sequences using residue masses stored in
#' [aa_properties]. When `charge > 0`, the function returns monoisotopic m/z
#' values using the proton mass `1.007276` Da.
#'
#' @param sequence Peptide sequence supplied as a character vector. Amino-acid
#'   codes must be one-letter uppercase or lowercase symbols from the 22
#'   supported residues in [aa_properties]. If `NULL`, raises an error.
#' @param charge Optional non-negative integer scalar or integer vector with
#'   length matching `sequence`. Defaults to `0L`. `0L` returns neutral masses.
#'   Positive values return m/z. If `NULL`, raises an error.
#'
#' @details This function computes masses for unmodified peptide sequences only.
#'   It does not account for chemical labels such as TMT or iTRAQ, isotopic
#'   labels such as SILAC, or post-translational modifications.
#'
#' @family utils
#' @section Limitations:
#'   This function computes monoisotopic mass only. Average mass is not
#'   supported.
#'
#' @return A numeric vector of neutral masses or m/z values. Names are
#'   preserved from `sequence` when present.
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

#' Calculate peptide isoelectric point
#'
#' `calculate_pI()` estimates peptide isoelectric points with a bisection
#' search over the Henderson-Hasselbalch net-charge equation. The calculation
#' uses hard-coded terminal pKa values of `8.0` for the N-terminus and `3.1`
#' for the C-terminus, plus side-chain pKa values from [aa_properties] for `C`,
#' `D`, `E`, `H`, `K`, `R`, `Y`, and `U`.
#'
#' @param sequence Peptide sequence supplied as a character vector. Amino-acid
#'   codes must be one-letter uppercase or lowercase symbols from the 22
#'   supported residues in [aa_properties]. If `NULL`, raises an error.
#'
#' @family utils
#' @section Limitations:
#'   pI estimation uses the Lehninger pKa set. Calculation may be slow for
#'   more than 5000 sequences.
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

  if (length(normalized_sequences) > .get_param("scatter_max_pts")) {
    cli::cli_inform(
      paste0(
        "Calculating peptide pI values for ",
        "{length(normalized_sequences)} sequences."
      ),
      class = "pepvet_message_calculating_pi"
    )
  }

  ## Use the internal helper to avoid re-normalizing sequences that were just
  ## validated above.
  composition_matrix <-
    .ionizable_composition_matrix_internal(normalized_sequences)
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
