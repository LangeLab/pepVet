#' Evaluate a Proteolytic Digest
#'
#' `evaluate_digest()` combines [digest_protein()] and [score_peptides()] into
#' a single call and returns a named list containing the peptide table, the
#' score table, and the resolved input parameters. Use it when you want a full
#' digest object for one protein and one enzyme without manually wiring the two
#' lower-level functions together.
#'
#' @param sequence Protein input. Accepts the same forms as [digest_protein()]:
#'   a character sequence, named character vector, `Biostrings::AAString`,
#'   `Biostrings::AAStringSet`, or a FASTA file path.
#' @param enzyme Enzyme name passed to [digest_protein()]. Defaults to
#'   `"trypsin"`.
#' @param missed_cleavages Maximum missed cleavages passed to
#'   [digest_protein()]. Defaults to `0L`.
#' @param include_cleavage_efficiency Logical flag passed to [digest_protein()].
#'   When `TRUE`, the returned peptide table gains a `cleavage_efficiency`
#'   column. This does not affect the score components.
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for peptide uniqueness scoring.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
#'   When scoring a non-tryptic digest directly, [evaluate_digest()] forwards
#'   the selected `enzyme` so enzyme-aware S_count denominators stay aligned
#'   with the digest.
#' @param ... Additional scoring arguments passed to [score_peptides()], such
#'   as `gravy_range` and `length_range`. This makes workflow presets from
#'   [pepvet_preset()] directly compatible with [evaluate_digest()] through
#'   `do.call()` or argument splicing.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{`scores`}{A tibble from [score_peptides()] with one row per
#'       protein, plus the informational columns `n_high_efficiency_sites` and
#'       `n_low_efficiency_sites`.}
#'     \item{`peptides`}{A tibble from [digest_protein()] with one row per
#'       peptide.}
#'     \item{`params`}{A list recording the resolved `enzyme` name,
#'       `missed_cleavages` count, `protein_ids` found in the input, and the
#'       resolved `preset_used` label from [score_peptides()].}
#'   }
#'
#' @details `evaluate_digest()` preserves pepVet's scoring metadata so the
#' returned object can be interpreted honestly outside the immediate scoring
#' call. In particular, `params$preset_used` records whether the resolved
#' scoring configuration matches one of pepVet's named presets or should be
#' treated as `"custom"`. The cleavage-efficiency counts summarize annotated
#' trypsin-family cleavage sites only; unsupported enzymes currently receive
#' `NA` in these informational fields.
#'
#' @seealso [digest_protein()], [score_peptides()], [compare_digests()],
#'   [batch_evaluate()]
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' result <- evaluate_digest(bsa_path, enzyme = "trypsin")
#' result$scores
#' result$params$enzyme
#' result$params$preset_used
#' @export
# nolint start: object_usage_linter.
evaluate_digest <- function(sequence,
                            enzyme = "trypsin",
                            missed_cleavages = 0L,
                            include_cleavage_efficiency = FALSE,
                            proteome = NULL,
                            weights = NULL,
                            ...) {
  normalized_input <- .read_input(sequence)

  peptides <- digest_protein(normalized_input,
    enzyme = enzyme,
    missed_cleavages = missed_cleavages,
    include_cleavage_efficiency = include_cleavage_efficiency
  )
  scores <- score_peptides(
    peptides,
    proteome = proteome,
    weights = weights,
    ...,
    enzyme = enzyme
  )
  normalized_enzyme <- .normalize_enzyme(enzyme)
  cleavage_counts <- lapply(
    seq_along(normalized_input),
    function(index) {
      counts <- .cleavage_efficiency_summary(
        as.character(normalized_input[[index]]),
        normalized_enzyme
      )

      tibble::tibble(
        protein_id = names(normalized_input)[[index]],
        n_high_efficiency_sites = counts$n_high_efficiency_sites,
        n_low_efficiency_sites = counts$n_low_efficiency_sites
      )
    }
  )
  cleavage_counts <- tibble::as_tibble(do.call(rbind, cleavage_counts))
  score_index <- match(scores$protein_id, cleavage_counts$protein_id)
  scores <- tibble::add_column(
    scores,
    n_high_efficiency_sites = cleavage_counts$n_high_efficiency_sites[score_index],
    n_low_efficiency_sites = cleavage_counts$n_low_efficiency_sites[score_index],
    .after = "preset_used"
  )

  list(
    scores = scores,
    peptides = peptides,
    params = list(
      enzyme = normalized_enzyme,
      missed_cleavages = as.integer(missed_cleavages),
      protein_ids = unique(peptides$protein_id),
      preset_used = scores$preset_used[[1L]]
    )
  )
}
# nolint end

#' Compare Multiple Enzymes on a Single Protein
#'
#' `compare_digests()` runs [evaluate_digest()] for each enzyme in `enzymes`
#' and returns a tibble of scores sorted by `composite_score` descending. It is
#' the main ranking function for pre-experimental enzyme selection.
#'
#' @param sequence A single-protein input. Accepts the same forms as
#'   [digest_protein()] but must resolve to exactly one protein.
#' @param enzymes Character vector of enzyme names to compare. Each name must
#'   be one of pepVet's supported cleaver-compatible enzyme names.
#' @param missed_cleavages Maximum missed cleavages passed to
#'   [digest_protein()] for every enzyme. Defaults to `0L`.
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for all enzyme evaluations.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
#' @param ... Additional arguments passed to [evaluate_digest()]. This includes
#'   scoring arguments such as `gravy_range`, `length_range`, and
#'   `include_pI`, plus `include_cleavage_efficiency` when peptide-level
#'   cleavage annotations are requested during comparison.
#'
#' @return A tibble with one row per enzyme and columns `enzyme` followed by
#'   the score columns returned by [evaluate_digest()], sorted by
#'   `composite_score` descending.
#'
#' @seealso [evaluate_digest()], [recommend_enzyme()]
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' compare_digests(bsa_path, enzymes = c("trypsin", "lysc"))
#' @export
# nolint start: object_usage_linter.
compare_digests <- function(sequence,
                            enzymes = c("trypsin", "lysc"),
                            missed_cleavages = 0L,
                            proteome = NULL,
                            weights = NULL,
                            ...) {
  if (!is.character(enzymes) || length(enzymes) == 0L || anyNA(enzymes)) {
    cli::cli_abort(
      paste(
        "{.arg enzymes} must be a non-empty character vector",
        "with no missing values."
      ),
      class = "pepvet_error_invalid_enzymes"
    )
  }

  normalized_input <- .read_input(sequence)

  if (length(normalized_input) != 1L) {
    cli::cli_abort(
      paste(
        "{.arg sequence} must resolve to exactly one protein",
        "for enzyme comparison."
      ),
      class = "pepvet_error_invalid_input"
    )
  }

  scored_rows <- lapply(enzymes, function(enzyme) {
    ev <- evaluate_digest(
      normalized_input,
      enzyme = enzyme,
      missed_cleavages = missed_cleavages,
      proteome = proteome,
      weights = weights,
      ...
    )
    tibble::add_column(ev$scores, enzyme = ev$params$enzyme, .before = 1L)
  })

  result <- tibble::as_tibble(do.call(rbind, scored_rows))
  result[order(result$composite_score, decreasing = TRUE), , drop = FALSE]
}
# nolint end

#' Recommend the Best Enzyme for a Single Protein
#'
#' `recommend_enzyme()` calls [compare_digests()] and returns the name of the
#' enzyme with the highest composite score. When two or more enzymes are tied,
#' all tied enzyme names are returned in alphabetical order. This function is
#' useful in scripted triage pipelines where you need a compact recommendation
#' but still want ranking logic that stays aligned with [compare_digests()].
#'
#' @param sequence A single-protein input passed to [compare_digests()].
#' @param enzymes Character vector of enzyme names to compare.
#' @param missed_cleavages Maximum missed cleavages. Defaults to `0L`.
#' @param proteome Optional proteome digest tibble for uniqueness scoring.
#' @param weights Optional scoring weight vector.
#' @param ... Additional scoring arguments passed to [compare_digests()] and
#'   ultimately to [evaluate_digest()] and [score_peptides()].
#'
#' @return A character vector of one or more enzyme names. Length greater than
#'   one only when top scores are tied within floating-point tolerance.
#'
#' @seealso [compare_digests()], [evaluate_digest()]
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' recommend_enzyme(bsa_path, enzymes = c("trypsin", "lysc"))
#' @export
# nolint start: object_usage_linter.
recommend_enzyme <- function(sequence,
                             enzymes = c("trypsin", "lysc"),
                             missed_cleavages = 0L,
                             proteome = NULL,
                             weights = NULL,
                             ...) {
  comparison <- compare_digests(
    sequence,
    enzymes = enzymes,
    missed_cleavages = missed_cleavages,
    proteome = proteome,
    weights = weights,
    ...
  )

  top_score <- max(comparison$composite_score)
  tied <- comparison$enzyme[
    abs(comparison$composite_score - top_score) < sqrt(.Machine$double.eps)
  ]
  sort(tied)
}
# nolint end

#' Batch-Evaluate Multiple Proteins
#'
#' `batch_evaluate()` calls [evaluate_digest()] independently for each protein
#' in `sequences` and returns a named list of results. The output for each
#' protein is bit-identical to calling [evaluate_digest()] on that protein
#' individually. Use it for small proteomes, panels, or fixture-level quality
#' assessment before moving to a larger parallel workflow.
#'
#' @param sequences Multi-protein input. Accepts the same forms as
#'   [digest_protein()]. Must resolve to at least one protein.
#' @param enzyme Enzyme name passed to [digest_protein()]. Defaults to
#'   `"trypsin"`.
#' @param missed_cleavages Maximum missed cleavages. Defaults to `0L`.
#' @param include_cleavage_efficiency Logical flag passed to
#'   [evaluate_digest()] and ultimately [digest_protein()]. When `TRUE`, each
#'   per-protein peptide table includes a `cleavage_efficiency` column.
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for every protein evaluation.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
#' @param ... Additional scoring arguments passed to [score_peptides()], such
#'   as `gravy_range` and `length_range`.
#'
#' @return A named list where each element is the result of [evaluate_digest()]
#'   for the corresponding protein. Names match the `protein_id` values from
#'   the input.
#'
#' @seealso [evaluate_digest()], [compare_digests()]
#'
#' @examples
#' small_proteome <- system.file(
#'   "extdata", "small_proteome_50_proteins.fasta",
#'   package = "pepVet"
#' )
#' results <- batch_evaluate(small_proteome, enzyme = "trypsin")
#' length(results)
#' results[[1]]$scores
#' @export
# nolint start: object_usage_linter.
batch_evaluate <- function(sequences,
                           enzyme = "trypsin",
                           missed_cleavages = 0L,
                           include_cleavage_efficiency = FALSE,
                           proteome = NULL,
                           weights = NULL,
                           ...) {
  normalized_input <- .read_input(sequences)

  results <- lapply(seq_along(normalized_input), function(index) {
    evaluate_digest(
      normalized_input[index],
      enzyme = enzyme,
      missed_cleavages = missed_cleavages,
      include_cleavage_efficiency = include_cleavage_efficiency,
      proteome = proteome,
      weights = weights,
      ...
    )
  })

  names(results) <- names(normalized_input)
  results
}
# nolint end
