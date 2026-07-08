.required_digest_columns <- c(
  "protein_id",
  "peptide",
  "start",
  "end",
  "length",
  "missed_cleavages"
)

.validate_digest_result <- function(digest_result, arg_name = "digest_result") {
  if (!inherits(digest_result, "data.frame")) {
    .abort(
      paste0(
        "{.arg {arg_name}} must be a digest tibble or data frame ",
        "with pepVet digest columns."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  missing_columns <- setdiff(.required_digest_columns, names(digest_result))

  if (length(missing_columns) > 0L) {
    .abort(
      c(
        "{.arg {arg_name}} is missing required digest columns.",
        "i" = "Missing columns: {.val {missing_columns}}"
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  if (nrow(digest_result) == 0L) {
    .abort(
      "{.arg {arg_name}} must contain at least one peptide row.",
      class = "pepvet_error_invalid_digest"
    )
  }

  digest_result <- digest_result[, .required_digest_columns, drop = FALSE]

  if (
    !is.character(digest_result$protein_id) ||
      !is.character(digest_result$peptide)
  ) {
    .abort(
      paste0(
        "{.arg {arg_name}} must store {.field protein_id} and ",
        "{.field peptide} as character columns."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  numeric_columns <- c("start", "end", "length", "missed_cleavages")

  if (!all(vapply(digest_result[numeric_columns], is.numeric, logical(1)))) {
    .abort(
      paste0(
        "{.arg {arg_name}} must store coordinate and count ",
        "columns as numeric values."
      ),
      class = "pepvet_error_invalid_digest"
    )
  }

  if (anyNA(digest_result)) {
    .abort(
      "{.arg {arg_name}} must not contain missing values.",
      class = "pepvet_error_invalid_digest"
    )
  }

  peptide_lengths <- nchar(digest_result$peptide, type = "chars")
  row_lengths <- as.integer(digest_result$length)
  row_starts <- as.integer(digest_result$start)
  row_ends <- as.integer(digest_result$end)
  row_missed <- as.integer(digest_result$missed_cleavages)

  if (any(row_starts < 1L) || any(row_ends < row_starts)) {
    .abort(
      "{.arg {arg_name}} contains invalid start/end coordinates.",
      class = "pepvet_error_invalid_digest"
    )
  }

  if (
    any(row_lengths != peptide_lengths) ||
      any(row_lengths != (row_ends - row_starts + 1L))
  ) {
    .abort(
      "{.arg {arg_name}} must have internally consistent peptide lengths.",
      class = "pepvet_error_invalid_digest"
    )
  }

  if (any(row_missed < 0L)) {
    .abort(
      "{.arg {arg_name}} must not contain negative missed-cleavage counts.",
      class = "pepvet_error_invalid_digest"
    )
  }

  digest_result$start <- row_starts
  digest_result$end <- row_ends
  digest_result$length <- row_lengths
  digest_result$missed_cleavages <- row_missed
  tibble::as_tibble(digest_result)
}

.valid_length_mask <- function(protein_digest, length_range = c(7L, 25L)) {
  protein_digest$length >= length_range[[1]] &
    protein_digest$length <= length_range[[2]]
}

.extract_valid_digest <- function(protein_digest, length_range = c(7L, 25L)) {
  protein_digest[
    .valid_length_mask(protein_digest, length_range), , drop = FALSE
  ]
}

.fallback_expected_peptide_length <- function(enzyme = "trypsin") {
  fallback_groups <- list(
    `12` = c("trypsin", "trypsin-high", "trypsin-low", "trypsin-simple"),
    `24` = c("lysc", "arg-c proteinase"),
    `17` = c("glutamyl endopeptidase", "asp-n endopeptidase"),
    `11` = c("chymotrypsin-high", "chymotrypsin-low")
  )

  for (fallback_value in names(fallback_groups)) {
    if (enzyme %in% fallback_groups[[fallback_value]]) {
      return(as.numeric(fallback_value))
    }
  }

  12
}

.expected_peptide_length <- function(protein_digest, enzyme = "trypsin") {
  peptide_count <- nrow(protein_digest)

  if (peptide_count < 3L) {
    return(as.numeric(.fallback_expected_peptide_length(enzyme)))
  }

  as.numeric(stats::median(protein_digest$length))
}

.has_no_cleavage_sites <- function(protein_digest) {
  nrow(protein_digest) == 1L &&
    protein_digest$start[[1]] == 1L &&
    protein_digest$end[[1]] == max(protein_digest$end)
}

.score_length <- function(protein_digest, length_range = c(7L, 25L)) {
  sum(.valid_length_mask(protein_digest, length_range)) / nrow(protein_digest)
}

.score_coverage <- function(protein_digest,
                             length_range = c(7L, 25L),
                             valid_digest = NULL) {
  if (is.null(valid_digest)) {
    valid_digest <- .extract_valid_digest(protein_digest, length_range)
  }

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  protein_length <- max(protein_digest$end)

  # Compute covered bases without S4 IRanges dispatch: sort intervals by start,
  # then use cummax of prior ends to merge overlaps in one vectorized pass.
  s <- valid_digest$start
  e <- valid_digest$end
  ord <- order(s)
  s <- s[ord]
  e <- e[ord]
  prev_max_end <- c(0L, cummax(e[-length(e)]))
  sum(pmax(0L, e - pmax(s - 1L, prev_max_end))) / protein_length
}

.score_count <- function(protein_digest,
                         enzyme = "trypsin",
                         length_range = c(7L, 25L),
                         valid_digest = NULL) {
  if (.has_no_cleavage_sites(protein_digest)) {
    cli::cli_warn(
      paste0(
        "Protein {.val {protein_digest$protein_id[[1]]}} has no cleavage ",
        "sites for {.val {enzyme}}. S_count set to 0."
      )
    )

    return(0)
  }

  if (is.null(valid_digest)) {
    valid_digest <- .extract_valid_digest(protein_digest, length_range)
  }
  valid_count <- nrow(valid_digest)
  protein_length <- max(protein_digest$end)
  expected_count <- protein_length /
    .expected_peptide_length(protein_digest, enzyme)

  min(valid_count / expected_count, 1)
}

.score_hydro <- function(protein_digest,
                         gravy_range = c(-1.0, 0.6),
                         length_range = c(7L, 25L),
                         valid_digest = NULL) {
  if (is.null(valid_digest)) {
    valid_digest <- .extract_valid_digest(protein_digest, length_range)
  }

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  # Batch GRAVY: builds the hydrophobicity lookup once and calls strsplit on
  # the whole vector instead of once per peptide.
  gravy_values <- .calculate_gravy_vec(valid_digest$peptide)
  tolerance <- sqrt(.Machine$double.eps)

  mean(
    gravy_values >= (gravy_range[[1]] - tolerance) &
      gravy_values <= (gravy_range[[2]] + tolerance)
  )
}

.score_charge <- function(protein_digest,
                           length_range = c(7L, 25L),
                           valid_digest = NULL) {
  if (is.null(valid_digest)) {
    valid_digest <- .extract_valid_digest(protein_digest, length_range)
  }

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  # Vectorized: check for any internal basic residue (K, R, H) using grepl
  # with a character class instead of strsplit + %in% per peptide.
  # Trim the last character to check only internal positions (non-C-terminal).
  # Single-residue peptides have no internal positions; treat as lacking charge.
  peptide_lengths <- nchar(valid_digest$peptide)
  internal_seqs <- ifelse(peptide_lengths > 1L,
    substr(valid_digest$peptide, 1L, peptide_lengths - 1L), "")
  mean(grepl("[KRH]", internal_seqs, perl = TRUE))
}

.score_unique <- function(protein_digest,
                          proteome_index,
                          length_range = c(7L, 25L),
                          valid_digest = NULL) {
  if (is.null(valid_digest)) {
    valid_digest <- .extract_valid_digest(protein_digest, length_range)
  }

  if (nrow(valid_digest) == 0L) {
    return(0)
  }

  protein_id <- unique(protein_digest$protein_id)

  if (length(protein_id) != 1L) {
    .abort(
      "{.arg protein_digest} must contain rows from exactly one protein.",
      class = "pepvet_error_invalid_digest"
    )
  }

  is_unique <- vapply(
    valid_digest$peptide,
    function(peptide) {
      indexed_proteins <- get0(
        peptide,
        envir = proteome_index,
        inherits = FALSE
      )

      if (is.null(indexed_proteins)) {
        return(TRUE)
      }

      length(indexed_proteins) == 1L && identical(indexed_proteins, protein_id)
    },
    logical(1)
  )

  mean(is_unique)
}

.classify_verdict <- function(composite_score) {
  ifelse(
    composite_score >= .get_param("verdict_good"),
    "Good",
    ifelse(composite_score >= .get_param("verdict_moderate"),
      "Moderate", "Poor")
  )
}

.score_components <- function(protein_digest,
                              proteome_index = NULL,
                              enzyme = "trypsin",
                              length_range = c(7L, 25L),
                              gravy_range = c(-1.0, 0.6)) {
  # Extract valid digest once to avoid 4 independent [.data.frame subsets.
  valid_digest <- .extract_valid_digest(protein_digest, length_range)

  component_scores <- c(
    S_length   = .score_length(protein_digest, length_range),
    S_coverage = .score_coverage(protein_digest, length_range, valid_digest),
    S_count    = .score_count(
      protein_digest,
      enzyme = enzyme,
      length_range = length_range,
      valid_digest = valid_digest
    ),
    S_hydro    = .score_hydro(
      protein_digest, gravy_range, length_range, valid_digest
    ),
    S_charge   = .score_charge(protein_digest, length_range, valid_digest)
  )

  if (!is.null(proteome_index)) {
    component_scores <- c(
      component_scores,
      S_unique = .score_unique(
        protein_digest, proteome_index, length_range, valid_digest
      )
    )
  }

  component_scores
}

#' Score digested peptides for proteomics suitability
#'
#' `score_peptides()` summarizes a pepVet digest tibble into per-protein scoring
#' components and a weighted composite suitability score. The scoring model is
#' designed for digest planning and enzyme comparison, not for post-search
#' peptide detectability prediction.
#'
#' @param digest_result A digest tibble produced by [digest_protein()] or an
#'   equivalent table containing the columns `protein_id`, `peptide`, `start`,
#'   `end`, `length`, and `missed_cleavages`. If `NULL` or missing required
#'   columns, raises an error.
#' @param proteome Optional digest tibble representing the comparison proteome
#'   used for peptide uniqueness scoring. When `NULL` (default), `S_unique` is
#'   excluded and protein-only default weights are used.
#' @param weights Optional numeric weight vector. In protein-only mode the
#'   default weights are `c(S_length = 0.200, S_coverage = 0.348,
#'   S_count = 0.226, S_hydro = 0.138, S_charge = 0.088)`.
#'   In proteome-aware mode the weights are
#'   `c(S_length = 0.160, S_coverage = 0.279, S_count = 0.181,
#'   S_hydro = 0.110, S_charge = 0.070, S_unique = 0.200)`.
#'   Weights were derived via analytical hierarchy process (AHP) with
#'   pairwise comparisons grounded in proteomics literature: peptide length
#'   is the primary MS detectability filter, sequence coverage carries the
#'   main biological signal, peptide count supports statistical confidence,
#'   and GRAVY/charge capture independent LC-MS dimensions.
#' @param gravy_range Numeric vector of length 2 defining the inclusive GRAVY
#'   range used by `S_hydro`. Defaults to `c(-1.0, 0.6)`. If `NULL`, raises an
#'   error.
#' @param length_range Integer vector of length 2 defining the inclusive valid
#'   peptide length range used by `S_length`, `S_coverage`, `S_count`,
#'   `S_hydro`, `S_charge`, and `S_unique`. Defaults to `c(7L, 25L)`. If
#'   `NULL`, raises an error.
#' @param enzyme Cleavage rule name used to choose the fallback expected peptide
#'   length when the digest contains fewer than three peptides. Defaults to
#'   `"trypsin"`. If `NULL`, raises an error.
#' @param include_pI Logical flag indicating whether to append a `pI` list
#'   column containing peptide-level pI values for valid peptides. Defaults to
#'   `FALSE`. If `NULL`, raises an error.
#'
#' @details Valid peptides are defined as peptides with lengths between 7 and
#'   25 residues inclusive by default, but this window can be changed with
#'   `length_range`. Coverage is computed from valid peptide coordinates with
#'   overlapping intervals reduced before coverage is summed. `S_hydro` uses
#'   the inclusive `gravy_range`, and `S_unique` is only computed when a
#'   proteome digest is supplied. `S_charge` measures the fraction of
#'   valid peptides that contain at least one non-terminal basic residue,
#'   capturing extra charge-state richness rather than baseline ionizability.
#'   Composite verdicts are classified as `Good` for scores >= 0.65,
#'   `Moderate` for scores >= 0.4, and `Poor` otherwise.
#'
#' @section Limitations:
#' pepVet is a rule-based digest-ranking model, not a peptide detectability
#' predictor. The verdict thresholds are heuristic, the default and preset
#' weights are expert priors rather than empirically fit coefficients, and the
#' scoring windows assume conventional C18 reversed-phase LC with ESI. pepVet
#' does not model PTMs, chemical labels, chromatographic gradients, or
#' instrument-specific fragmentation behavior. Composite scores are
#' interpretable rankings with no physical unit; they are not calibrated
#' probabilities. Cross-workflow comparisons are only meaningful when the
#' resolved scoring configuration matches, which is why `preset_used` is
#' recorded in the output.
#'
#' @return A tibble with one row per `protein_id` and the component score
#'   columns `S_length`, `S_coverage`, `S_count`, `S_hydro`,
#'   `S_charge`,
#'   optional `S_unique`, plus `composite_score`, `verdict`, and
#'   `median_peptide_length`, and `preset_used`. The
#'   `median_peptide_length` column records the digest-level denominator used in
#'   the enzyme-aware `S_count` calculation. The `preset_used` column records
#'   the named preset whose resolved scoring configuration exactly matches the
#'   current call, or `"custom"` otherwise. When `include_pI = TRUE`, the
#'   output also includes a `pI` list column with one named numeric vector per
#'   protein, storing valid-peptide pI values keyed by peptide sequence.
#'
#' @examples
#' digest_result <- digest_protein("MKWVTFISLLFLFSSAYSR")
#' score_peptides(digest_result)
#'
#' target_and_background <- digest_protein(
#'   c(target = "AAAAAAARAAAAAAAK", background = "AAAAAAARGGGGGGGK")
#' )
#' score_peptides(
#'   target_and_background[target_and_background$protein_id == "target", ],
#'   proteome = target_and_background
#' )
#' @family scoring
#' @export
score_peptides <- function(digest_result,
                           proteome = NULL,
                           weights = NULL,
                           gravy_range = c(-1.0, 0.6),
                           length_range = c(7L, 25L),
                           enzyme = "trypsin",
                           include_pI = FALSE) {
  validated_digest <- .validate_digest_result(digest_result)
  has_proteome <- !is.null(proteome)
  proteome_index <- NULL
  normalized_enzyme <- .normalize_enzyme(enzyme)
  normalized_gravy_range <- .validate_gravy_range(gravy_range)
  normalized_length_range <- .validate_length_range(length_range)
  normalized_include_pI <- .validate_include_pI(include_pI)

  if (has_proteome) {
    validated_proteome <- .validate_digest_result(
      proteome,
      arg_name = "proteome"
    )
    proteome_index <- .build_proteome_index(validated_proteome)
  }

  normalized_weights <- .validate_weights(weights, has_proteome)
  preset_used <- .identify_preset_used(
    gravy_range = normalized_gravy_range,
    length_range = normalized_length_range,
    weights = normalized_weights,
    include_pI = normalized_include_pI,
    has_proteome = has_proteome
  )
  protein_levels <- unique(validated_digest$protein_id)
  protein_groups <- split(
    seq_len(nrow(validated_digest)),
    factor(validated_digest$protein_id, levels = protein_levels)
  )

  # Convert to base data.frame once so that per-protein subsetting uses
  # [.data.frame instead of the slower [.tbl_df throughout the inner loop.
  validated_digest <- as.data.frame(validated_digest, stringsAsFactors = FALSE)

  scored_rows <- lapply(
    protein_groups,
    function(row_indices) {
      protein_digest <- validated_digest[row_indices, , drop = FALSE]
      component_scores <- .score_components(
        protein_digest,
        proteome_index,
        enzyme = normalized_enzyme,
        length_range = normalized_length_range,
        gravy_range = normalized_gravy_range
      )
      weighted_components <- component_scores[names(normalized_weights)]
      composite_score <- sum(weighted_components * normalized_weights)

      # Zero-cleavage hard-fail: S_count == 0 and S_coverage == 0 means
      # the protein was not digested by this enzyme.  Override composite
      # and verdict so an undigestible protein is always "Poor".
      if (component_scores[["S_count"]] == 0 &&
            component_scores[["S_coverage"]] == 0) {
        composite_score <- 0
      }

      median_peptide_length <- .expected_peptide_length(
        protein_digest,
        normalized_enzyme
      )
      row_values <- c(
        list(protein_id = protein_digest$protein_id[[1]]),
        as.list(component_scores),
        list(
          composite_score = composite_score,
          verdict = .classify_verdict(composite_score),
          median_peptide_length = median_peptide_length,
          preset_used = preset_used
        )
      )

      if (isTRUE(normalized_include_pI)) {
        valid_digest <- .extract_valid_digest(
          protein_digest, normalized_length_range
        )
        peptide_pI <- if (nrow(valid_digest) == 0L) {
          numeric(0)
        } else {
          stats::setNames(
            calculate_pI(valid_digest$peptide), valid_digest$peptide
          )
        }

        row_values$pI <- list(peptide_pI)
      }

      tibble::as_tibble(row_values)
    }
  )

  .bind_rows(scored_rows)
}
