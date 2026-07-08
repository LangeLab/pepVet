#' Evaluate a proteolytic digest
#'
#' `evaluate_digest()` combines [digest_protein()] and [score_peptides()] into
#' a single call and returns a named list containing the peptide table, the
#' score table, and the resolved input parameters. Use it when you want a full
#' digest object for one protein and one enzyme without manually wiring the two
#' lower-level functions together.
#'
#' @param sequence Protein input. Accepts the same forms as [digest_protein()]:
#'   a character sequence, named character vector, `Biostrings::AAString`,
#'   `Biostrings::AAStringSet`, or a FASTA file path. If `NULL` or empty,
#'   raises an error.
#' @param enzyme Enzyme name passed to [digest_protein()]. Defaults to
#'   `"trypsin"`.
#' @param missed_cleavages Maximum missed cleavages passed to
#'   [digest_protein()]. Defaults to `1L`.
#' @param include_cleavage_efficiency Logical flag passed to [digest_protein()].
#'   When `TRUE`, the returned peptide table gains a `cleavage_efficiency`
#'   column. This does not affect the score components.
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for peptide uniqueness scoring. When `NULL` (default), no uniqueness
#'   scoring is performed.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
#'   When `NULL` (default), uses pepVet's default scoring weights.
#'   When scoring a non-tryptic digest directly, [evaluate_digest()] forwards
#'   the selected `enzyme` so enzyme-aware S_count denominators stay aligned
#'   with the digest.
#' @param ... Additional scoring arguments passed to [score_peptides()], such
#'   as `gravy_range` and `length_range`. This makes workflow presets from
#'   [pepvet_preset()] directly compatible with [evaluate_digest()] through
#'   `do.call()` or argument splicing.
#'
#' @details `evaluate_digest()` preserves pepVet's scoring metadata so the
#' returned object can be interpreted honestly outside the immediate scoring
#' call. In particular, `params$preset_used` records whether the resolved
#' scoring configuration matches one of pepVet's named presets or should be
#' treated as `"custom"`. The cleavage-efficiency counts summarize annotated
#' trypsin-family cleavage sites only; unsupported enzymes receive `NA` in
#' these informational fields.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{\code{scores}}{A tibble from [score_peptides()] with one row per
#'       protein, plus the informational columns \code{n_high_efficiency_sites}
#'       and \code{n_low_efficiency_sites}.}
#'     \item{\code{peptides}}{A tibble from [digest_protein()] with one row per
#'       peptide.}
#'     \item{\code{params}}{A list recording the resolved \code{enzyme} name,
#'       \code{missed_cleavages} count, \code{protein_ids} found in the input,
#'       and the resolved \code{preset_used} label from [score_peptides()].}
#'   }
#'
#' @section Limitations:
#' Cleavage efficiency annotations only work for trypsin-family enzymes.
#' Unsupported enzymes receive `NA` counts. The scoring model is rule-based;
#' see [score_peptides()] for scope details.
#'
#' @family evaluation
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
                            missed_cleavages = 1L,
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
  cleavage_counts <- .bind_rows(cleavage_counts)
  score_index <- match(scores$protein_id, cleavage_counts$protein_id)
  scores <- tibble::add_column(
    scores,
    n_high_efficiency_sites = cleavage_counts$n_high_efficiency_sites[
      score_index
    ],
    n_low_efficiency_sites = cleavage_counts$n_low_efficiency_sites[
      score_index
    ],
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

#' Compare multiple enzymes on a single protein
#'
#' `compare_digests()` runs [evaluate_digest()] for each enzyme in `enzymes`
#' and returns a tibble of scores sorted by `composite_score` descending. Main
#' ranking function for pre-experimental enzyme selection.
#'
#' @param sequence A single-protein input. Accepts the same forms as
#'   [digest_protein()] but must resolve to exactly one protein. If `NULL` or
#'   empty, raises an error.
#' @param enzymes Character vector of enzyme names to compare. Defaults to
#'   `c("trypsin", "lysc")`. Each name must be one of pepVet's supported
#'   cleaver-compatible enzyme names.
#' @param missed_cleavages Maximum missed cleavages passed to
#'   [digest_protein()] for every enzyme. Defaults to `1L`.
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for all enzyme evaluations. When `NULL` (default), no uniqueness scoring.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
#'   When `NULL` (default), uses pepVet's default scoring weights.
#' @param ... Additional arguments passed to [evaluate_digest()]. This includes
#'   scoring arguments such as `gravy_range`, `length_range`, and
#'   `include_pI`, plus `include_cleavage_efficiency` when peptide-level
#'   cleavage annotations are requested during comparison.
#'
#' @inheritParams evaluate_digest
#'
#' @return A tibble with one row per enzyme and columns `enzyme` followed by
#'   the score columns returned by [evaluate_digest()], sorted by
#'   `composite_score` descending.
#'
#' @section Limitations:
#' Single-protein only. Use [batch_compare_enzymes()] for proteome-wide
#' enzyme comparison.
#'
#' @family evaluation
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
                            missed_cleavages = 1L,
                            proteome = NULL,
                            weights = NULL,
                            ...) {
  if (!is.character(enzymes) || length(enzymes) == 0L || anyNA(enzymes)) {
    .abort(
      paste0("{.arg enzymes} must be a non-empty character vector ",
        "with no missing values."),
      class = "pepvet_error_invalid_enzymes"
    )
  }

  normalized_input <- .read_input(sequence)

  if (length(normalized_input) != 1L) {
    .abort(
      paste0("{.arg sequence} must resolve to exactly one protein ",
        "for enzyme comparison."),
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

  result <- .bind_rows(scored_rows)
  result[order(result$composite_score, decreasing = TRUE), , drop = FALSE]
}
# nolint end

#' Recommend the best enzyme for a single protein
#'
#' `recommend_enzyme()` calls [compare_digests()] and returns the name of the
#' enzyme with the highest composite score. When two or more enzymes are tied,
#' all tied enzyme names are returned in alphabetical order. Compact
#' recommendation for scripted triage pipelines that stays aligned with
#' [compare_digests()].
#'
#' @param sequence A single-protein input passed to [compare_digests()]. If
#'   `NULL` or empty, raises an error.
#' @param enzymes Character vector of enzyme names to compare. Defaults to
#'   `c("trypsin", "lysc")`.
#' @param missed_cleavages Maximum missed cleavages. Defaults to `1L`.
#' @param proteome Optional proteome digest tibble for uniqueness scoring.
#'   When `NULL` (default), no uniqueness scoring.
#' @param weights Optional scoring weight vector. When `NULL` (default), uses
#'   pepVet's default scoring weights.
#' @param ... Additional scoring arguments passed to [compare_digests()] and
#'   ultimately to [evaluate_digest()] and [score_peptides()].
#'
#' @inheritParams compare_digests
#'
#' @return A character vector of one or more enzyme names. Length greater than
#'   one only when top scores are tied within floating-point tolerance.
#'
#' @section Limitations:
#' Single-protein only. When multiple enzymes tie within tolerance, all are
#' returned in alphabetical order with no further tie-breaking.
#'
#' @family evaluation
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
                             missed_cleavages = 1L,
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
    abs(comparison$composite_score - top_score) < 1e-5
  ]
  sort(tied)
}
# nolint end

#' Batch-evaluate multiple proteins
#'
#' `batch_evaluate()` calls [evaluate_digest()] independently for each protein
#' in `sequences` and returns a flat tibble with one row per protein. Columns
#' include `protein_id`, `protein_length`, all component scores,
#'   `composite_score`,
#' `verdict`, `n_peptides`, `n_valid_peptides`,
#'   `median_peptide_length`, and four
#' sequence-level difficulty flags. Pass the result to [summarize_batch()] for
#' aggregate statistics or to [triage_proteins()] for action recommendations.
#'
#' @param sequences Multi-protein input. Accepts the same forms as
#'   [digest_protein()]. Must resolve to at least one protein. If `NULL` or
#'   empty, raises an error.
#' @param enzyme Enzyme name passed to [digest_protein()]. Defaults to
#'   `"trypsin"`.
#' @param missed_cleavages Maximum missed cleavages. Defaults to `1L`.
#' @param include_cleavage_efficiency Logical flag passed to
#'   [evaluate_digest()] and ultimately [digest_protein()]. Defaults to `FALSE`.
#'   When `TRUE`, each per-protein peptide table includes a
#'   `cleavage_efficiency`
#'   column (does not affect the flat batch tibble columns).
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for every protein evaluation. When `NULL` (default), no uniqueness
#'   scoring is performed and the `S_unique` column is omitted.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
#'   When `NULL` (default), uses pepVet's default scoring weights.
#' @param cores Number of parallel workers for protein-level chunking.
#'   `1L` (default) runs sequentially with no extra dependencies. Values
#'   greater than `1L` split the input into equal chunks and process each
#'   chunk in a forked worker via `parallel::mclapply` (Unix only; on Windows
#'   `cores` is silently capped to `1L`).
#' @param ... Additional scoring arguments passed to [score_peptides()], such
#'   as `gravy_range` and `length_range`.
#'
#' @return A tibble with one row per protein. Fixed columns:
#'   `protein_id`, `protein_length`, `n_peptides`, `n_valid_peptides`,
#'   all available component scores (`S_length`, `S_coverage`, `S_count`,
#'   `S_hydro`, `S_charge`; `S_unique` when proteome is provided),
#'   `composite_score`, `verdict`, `median_peptide_length`,
#'   `flag_short_protein`, `flag_hydrophobic`, `flag_low_complexity`,
#'   `flag_no_valid_peptides`.
#'
#' @section Limitations:
#' Parallel execution only on Unix (`parallel::mclapply`). Windows falls
#' back to serial silently. If a worker crashes, the failed chunk is retried
#' sequentially, which is slower but produces correct results.
#'
#' @family evaluation
#'
#' @seealso [evaluate_digest()], [compare_digests()], [summarize_batch()],
#'   [triage_proteins()]
#'
#' @examples
#' small_proteome <- system.file(
#'   "extdata", "small_proteome_50_proteins.fasta",
#'   package = "pepVet"
#' )
#' batch <- batch_evaluate(small_proteome, enzyme = "trypsin")
#' nrow(batch)
#' batch[, c("protein_id", "composite_score", "verdict")]
#' @export
# nolint start: object_usage_linter.
batch_evaluate <- function(sequences,
                           enzyme = "trypsin",
                           missed_cleavages = 1L,
                           include_cleavage_efficiency = FALSE,
                           proteome = NULL,
                           weights = NULL,
                           cores = 1L,
                           ...) {
  cores <- suppressWarnings(as.integer(cores))
  if (is.na(cores) || cores < 1L) {
    .abort(
      "{.arg cores} must be a positive integer.",
      class = "pepvet_error_invalid_cores"
    )
  }

  normalized_input <- .read_input(sequences)
  extra_args <- list(...)
  n_proteins <- length(normalized_input)

  ## On Windows mclapply is unavailable; silently run serial.
  effective_cores <- if (.Platform$OS.type == "windows") {
    1L
  } else {
    min(cores, n_proteins)
  }

  if (effective_cores > 1L) {
    idx_list <- parallel::splitIndices(n_proteins, effective_cores)
    results <- parallel::mclapply(idx_list, function(idx) {
      .batch_evaluate_inner(
        normalized_input[idx], enzyme, missed_cleavages,
        include_cleavage_efficiency, proteome, weights, extra_args
      )
    }, mc.cores = effective_cores)

    ## Check for worker failures (mclapply returns try-error on crash).
    ## Retry failed chunks sequentially. Slower but correct.
    failed <- vapply(results, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      n_failed <- sum(failed)
      cli::cli_warn(
        paste0("{n_failed} parallel worker{?s} failed. Retrying failed ",
          "chunk{?s} sequentially."),
        class = "pepvet_warning_parallel_retry"
      )
      for (i in which(failed)) {
        results[[i]] <- .batch_evaluate_inner(
          normalized_input[idx_list[[i]]], enzyme, missed_cleavages,
          include_cleavage_efficiency, proteome, weights, extra_args
        )
      }
    }

    return(.bind_rows(results))
  }

  .batch_evaluate_inner(
    normalized_input, enzyme, missed_cleavages,
    include_cleavage_efficiency, proteome, weights, extra_args
  )
}
# nolint end

## Private batch helpers

## Core pipeline for a pre-parsed AAStringSet, called by batch_evaluate() both
## in serial mode and from within each mclapply worker. Accepts extra scoring
## arguments as a captured list (extra_args) so they survive fork serialization.
.batch_evaluate_inner <- function(normalized_input, enzyme, missed_cleavages,
                                  include_cleavage_efficiency, proteome,
                                  weights, extra_args) {
  protein_ids <- names(normalized_input)

  all_peptides <- digest_protein(
    normalized_input,
    enzyme = enzyme,
    missed_cleavages = missed_cleavages,
    include_cleavage_efficiency = include_cleavage_efficiency
  )

  all_scores <- do.call(
    score_peptides,
    c(
      list(
        all_peptides,
        proteome = proteome,
        weights  = weights,
        enzyme   = enzyme
      ),
      extra_args
    )
  )

  score_cols <- intersect(
    names(all_scores),
    c(
      "S_length", "S_coverage", "S_count", "S_hydro", "S_charge",
      "S_unique", "composite_score", "verdict", "median_peptide_length"
    )
  )

  flags <- .batch_difficulty_flags(all_peptides, protein_ids)

  scores_reordered <- all_scores[
    match(protein_ids, all_scores$protein_id),
    score_cols,
    drop = FALSE
  ]

  tibble::as_tibble(c(
    list(
      protein_id       = protein_ids,
      protein_length   = flags$protein_length,
      n_peptides       = flags$n_peptides,
      n_valid_peptides = flags$n_valid_peptides
    ),
    as.list(scores_reordered),
    list(
      flag_short_protein     = flags$flag_short_protein,
      flag_hydrophobic       = flags$flag_hydrophobic,
      flag_low_complexity    = flags$flag_low_complexity,
      flag_no_valid_peptides = flags$flag_no_valid_peptides
    )
  ))
}

.validate_batch_result <- function(batch_result) {
  if (!inherits(batch_result, "data.frame")) {
    .abort(
      "{.arg batch_result} must be a tibble returned by {.fn batch_evaluate}.",
      class = "pepvet_error_invalid_batch_result"
    )
  }

  if (nrow(batch_result) == 0L) {
    .abort(
      "{.arg batch_result} must contain at least one protein row.",
      class = "pepvet_error_invalid_batch_result"
    )
  }

  required_cols <- c(
    "protein_id", "composite_score", "verdict",
    "flag_short_protein", "flag_hydrophobic",
    "flag_low_complexity", "flag_no_valid_peptides"
  )
  missing_cols <- setdiff(required_cols, names(batch_result))

  if (length(missing_cols) > 0L) {
    .abort(
      c(
        paste0("{.arg batch_result} is missing required columns from ",
          "{.fn batch_evaluate}."),
        "i" = "Missing: {.val {missing_cols}}"
      ),
      class = "pepvet_error_invalid_batch_result"
    )
  }

  batch_result
}

## Vectorized difficulty flags for a full multi-protein peptide tibble.
## Returns a named list of vectors, each of length == length(protein_ids),
## in the same order as protein_ids.
.batch_difficulty_flags <- function(all_peptides, protein_ids) {
  pid_factor <- factor(all_peptides$protein_id, levels = protein_ids)

  ## protein_length and n_peptides
  protein_length <- as.integer(tapply(all_peptides$end, pid_factor, max))
  n_peptides <- as.integer(tabulate(pid_factor))

  ## valid peptide mask (length 7-25)
  valid_mask <-
    all_peptides$length >= .get_param("length_lo") &
      all_peptides$length <= .get_param("length_hi")
  n_valid_peptides <- as.integer(tabulate(pid_factor[valid_mask],
    nbins = length(protein_ids)))

  ## flags derivable from counts
  flag_short_protein <- protein_length < 100L
  flag_no_valid_peptides <- n_valid_peptides == 0L

  ## flag_hydrophobic: median GRAVY of valid peptides > 0.6 per protein
  flag_hydrophobic <- logical(length(protein_ids))
  if (any(valid_mask)) {
    gravy_vals <- .calculate_gravy(all_peptides$peptide[valid_mask])
    median_gravy <- tapply(gravy_vals, pid_factor[valid_mask], stats::median)
    flag_hydrophobic[match(names(median_gravy), protein_ids)] <-
      median_gravy > 0.6
  }

  ## flag_low_complexity: dominant AA > 50% in reconstructed MC=0 sequence.
  ## Build one concatenated sequence per protein from MC=0 peptides (sorted by
  ## start), then check character frequencies.
  mc0_mask <- all_peptides$missed_cleavages == 0L
  mc0_peps <- all_peptides[mc0_mask, c("protein_id", "start", "peptide"),
    drop = FALSE
  ]
  mc0_peps <- mc0_peps[order(mc0_peps$protein_id, mc0_peps$start), ,
    drop = FALSE
  ]
  prot_seqs <- tapply(mc0_peps$peptide, mc0_peps$protein_id,
    paste,
    collapse = "", simplify = FALSE
  )

  flag_low_complexity <- logical(length(protein_ids))
  names(flag_low_complexity) <- protein_ids
  for (pid in names(prot_seqs)) {
    s <- prot_seqs[[pid]]
    if (nchar(s) > 0L) {
      chars <- strsplit(s, "", fixed = TRUE)[[1L]]
      flag_low_complexity[[pid]] <-
        max(tabulate(match(chars, unique(chars)))) / length(chars) > 0.5
    }
  }

  list(
    protein_length = protein_length,
    n_peptides = n_peptides,
    n_valid_peptides = n_valid_peptides,
    flag_short_protein = flag_short_protein,
    flag_no_valid_peptides = flag_no_valid_peptides,
    flag_hydrophobic = flag_hydrophobic,
    flag_low_complexity = unname(flag_low_complexity)
  )
}

#' Summarize a batch digest evaluation
#'
#' `summarize_batch()` extracts aggregate statistics from a [batch_evaluate()]
#' result tibble. Returns a named list with verdict distribution, score
#' distribution, per-component averages, the lowest-scoring proteins, and a
#' heuristic set of enzyme-switch candidates.
#'
#' @param batch_result A tibble returned by [batch_evaluate()]. If `NULL` or
#'   empty, raises an error.
#'
#' @details `enzyme_switch_candidates` is a heuristic flag list derived from
#'   sequence-level difficulty flags, not from running alternative enzymes.
#'   Use [compare_digests()] to confirm whether a specific alternative enzyme
#'   improves the verdict for a flagged protein.
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{\code{verdict_counts}}{A tibble with columns \code{verdict},
#'       \code{n}, and \code{pct} covering the three verdict categories.}
#'     \item{\code{score_distribution}}{A named numeric vector with \code{mean},
#'       \code{median}, \code{sd}, \code{q25}, \code{q75}, \code{min}, and
#'       \code{max} of composite scores.}
#'     \item{\code{component_summary}}{A named numeric vector of per-component
#'       mean scores. The lowest values identify the weakest scoring dimension
#'       across the proteome.}
#'     \item{\code{problem_proteins}}{A tibble of proteins in the bottom 10% by
#'       composite score, ordered ascending, with all score and flag columns.}
#'     \item{\code{enzyme_switch_candidates}}{A tibble of Moderate or Poor
#'       proteins where \code{flag_hydrophobic} or \code{flag_short_protein} is
#'       \code{TRUE}, indicating that enzyme or preset choice is the likely
#'       limiting factor.}
#'   }
#'
#' @section Limitations:
#' The `enzyme_switch_candidates` are heuristic flags based on sequence
#' difficulty, not actual re-evaluation with alternative enzymes. Use
#' [compare_digests()] to confirm whether a switch improves the verdict.
#'
#' @family evaluation
#'
#' @seealso [batch_evaluate()], [triage_proteins()]
#'
#' @examples
#' small_path <- system.file(
#'   "extdata", "small_proteome_50_proteins.fasta",
#'   package = "pepVet"
#' )
#' batch <- batch_evaluate(small_path, enzyme = "trypsin")
#' summary <- summarize_batch(batch)
#' summary$verdict_counts
#' summary$component_summary
#' @export
# nolint start: object_usage_linter.
summarize_batch <- function(batch_result) {
  .validate_batch_result(batch_result)
  flat <- batch_result

  verdict_levels <- c("Good", "Moderate", "Poor")
  counts <- vapply(
    verdict_levels,
    function(v) sum(flat$verdict == v),
    integer(1)
  )
  verdict_counts <- tibble::tibble(
    verdict = verdict_levels,
    n       = counts,
    pct     = round(counts / nrow(flat) * 100, 1)
  )

  scores <- flat$composite_score
  score_distribution <- c(
    mean   = mean(scores),
    median = stats::median(scores),
    sd     = stats::sd(scores),
    q25    = stats::quantile(scores, 0.25, names = FALSE),
    q75    = stats::quantile(scores, 0.75, names = FALSE),
    min    = min(scores),
    max    = max(scores)
  )

  score_cols <- intersect(
    names(flat),
    c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge", "S_unique")
  )
  component_summary <- vapply(flat[score_cols], mean, numeric(1))

  threshold <- stats::quantile(flat$composite_score, 0.1, names = FALSE)
  problem_mask <- flat$composite_score <= threshold
  problem_proteins <- flat[problem_mask, , drop = FALSE]
  problem_proteins <- problem_proteins[
    order(problem_proteins$composite_score), ,
    drop = FALSE
  ]

  switch_mask <- flat$verdict %in% c("Moderate", "Poor") &
    (flat$flag_hydrophobic | flat$flag_short_protein)
  keep_cols <- intersect(
    names(flat),
    c(
      "protein_id", "verdict", "composite_score",
      "flag_hydrophobic", "flag_short_protein"
    )
  )
  enzyme_switch_candidates <- flat[switch_mask, keep_cols, drop = FALSE]

  list(
    verdict_counts           = verdict_counts,
    score_distribution       = score_distribution,
    component_summary        = component_summary,
    problem_proteins         = problem_proteins,
    enzyme_switch_candidates = enzyme_switch_candidates
  )
}
# nolint end

#' Compare multiple enzymes across a full proteome
#'
#' `batch_compare_enzymes()` scores every protein in `sequences` against each
#' enzyme in `enzymes` and returns a tidy tibble with one row per
#' protein-enzyme pair. The input proteome is parsed once; each enzyme then
#' calls [batch_evaluate()] with `cores` workers that split proteins into
#' equal chunks processed via fork copy-on-write (`parallel::mclapply`).
#' This means `cores` workers are fully utilised regardless of how many enzymes
#' are compared, and the speedup applies equally to single-enzyme calls.
#'
#' @param sequences Multi-protein input. Accepts the same forms as
#'   [digest_protein()]. Must resolve to at least one protein. If `NULL` or
#'   empty, raises an error.
#' @param enzymes Character vector of enzyme names to compare. Each name must
#'   be one of pepVet's supported enzyme names. Defaults to a panel of five
#'   commonly compared enzymes: trypsin, lysc, chymotrypsin-high,
#'   asp-n endopeptidase, and arg-c proteinase.
#' @param cores Number of parallel workers passed to [batch_evaluate()] for
#'   each enzyme. Proteins are split into `cores` equal chunks and processed
#'   in parallel via `parallel::mclapply` (Unix only; silently serial on
#'   Windows). Enzymes are always processed sequentially.
#' @param missed_cleavages Maximum missed cleavages passed to
#'   [batch_evaluate()] for every enzyme. Defaults to `1L`.
#' @param proteome Optional proteome digest tibble passed to [batch_evaluate()]
#'   for every enzyme. When `NULL` (default), no uniqueness scoring and the
#'   `S_unique` column is omitted.
#' @param weights Optional scoring weight vector passed to [batch_evaluate()].
#'   When `NULL` (default), uses pepVet's default scoring weights.
#' @param ... Additional scoring arguments passed to [batch_evaluate()], such
#'   as `gravy_range` and `length_range`.
#'
#' @inheritParams batch_evaluate
#'
#' @return A tibble of class `pepvet_batch_comparison` with one row per
#'   protein-enzyme pair. Columns: `protein_id`, `enzyme` (factor ordered
#'   by the input `enzymes` vector, so ggplot2 axis and facet order matches
#'   your specification), then all columns returned by [batch_evaluate()]:
#'   `protein_length`, `n_peptides`, `n_valid_peptides`, component scores,
#'   `composite_score`, `verdict`, `median_peptide_length`, and the four
#'   difficulty flags. Printing shows a per-enzyme summary table before the
#'   tibble rows.
#'
#' @section Limitations:
#' Enzymes run sequentially even when `cores > 1`; only proteins within each
#' enzyme are parallelized. Parallel execution requires Unix.
#'
#' @family evaluation
#'
#' @seealso [batch_evaluate()], [compare_digests()], [summarize_batch()]
#'
#' @examples
#' small <- system.file(
#'   "extdata", "small_proteome_50_proteins.fasta",
#'   package = "pepVet"
#' )
#' result <- batch_compare_enzymes(small, enzymes = c("trypsin", "lysc"))
#' result
#' result[result$enzyme == "trypsin", c("protein_id", "composite_score")]
#' @export
# nolint start: object_usage_linter.
batch_compare_enzymes <- function(
  sequences,
  enzymes = c(
    "trypsin", "lysc", "chymotrypsin-high",
    "asp-n endopeptidase", "arg-c proteinase"
  ),
  cores = 1L,
  missed_cleavages = 1L,
  proteome = NULL,
  weights = NULL,
  ...
) {
  if (!is.character(enzymes) || length(enzymes) == 0L || anyNA(enzymes)) {
    .abort(
      "{.arg enzymes} must be a non-empty character vector with no NAs.",
      class = "pepvet_error_invalid_enzymes"
    )
  }

  cores <- suppressWarnings(as.integer(cores))
  if (is.na(cores) || cores < 1L) {
    .abort(
      "{.arg cores} must be a positive integer.",
      class = "pepvet_error_invalid_cores"
    )
  }

  ## Parse input once. Each per-enzyme batch_evaluate() call receives the
  ## same in-memory AAStringSet, which fork workers share via copy-on-write.
  normalized_input <- .read_input(sequences)
  normalized_enzymes <- vapply(enzymes, .normalize_enzyme, character(1L),
    USE.NAMES = FALSE
  )
  n_proteins <- length(normalized_input)
  n_enzymes <- length(normalized_enzymes)

  cli::cli_inform(
    "Scoring {n_proteins} protein{?s} against {n_enzymes} enzyme{?s}.",
    class = "pepvet_message_batch_scoring"
  )

  extra_args <- list(...)

  results <- lapply(normalized_enzymes, function(enz) {
    row <- do.call(
      batch_evaluate,
      c(
        list(
          sequences        = normalized_input,
          enzyme           = enz,
          missed_cleavages = missed_cleavages,
          proteome         = proteome,
          weights          = weights,
          cores            = cores
        ),
        extra_args
      )
    )
    tibble::add_column(row, enzyme = enz, .after = "protein_id")
  })

  combined <- .bind_rows(results)
  combined$enzyme <- factor(combined$enzyme, levels = normalized_enzymes)

  result <- combined
  class(result) <- c("pepvet_batch_comparison", class(result))
  attr(result, "n_proteins") <- n_proteins
  attr(result, "n_enzymes") <- n_enzymes
  attr(result, "enzymes") <- normalized_enzymes
  result
}
# nolint end

#' @export
print.pepvet_batch_comparison <- function(x, ...) {
  has_summary_shape <- all(
    c("enzyme", "composite_score", "verdict") %in% names(x)
  )

  if (!has_summary_shape) {
    plain_x <- x
    class(plain_x) <- setdiff(class(x), "pepvet_batch_comparison")
    print(plain_x, ...)
    return(invisible(x))
  }

  enzymes <- attr(x, "enzymes")
  if (is.null(enzymes) || length(enzymes) == 0L) {
    enzymes <- unique(as.character(x$enzyme))
  }

  n_prot <- attr(x, "n_proteins")
  if (is.null(n_prot) || length(n_prot) != 1L || is.na(n_prot)) {
    n_prot <- length(unique(x$protein_id))
  }

  n_enz <- attr(x, "n_enzymes")
  if (is.null(n_enz) || length(n_enz) != 1L || is.na(n_enz)) {
    n_enz <- length(enzymes)
  }

  cli::cli_text(
    "pepVet batch enzyme comparison: {n_prot} protein{?s} x {n_enz} enzyme{?s}"
  )

  summary_rows <- lapply(enzymes, function(enz) {
    sub <- x[as.character(x$enzyme) == enz, , drop = FALSE]
    data.frame(
      enzyme = enz,
      n = nrow(sub),
      med_score = round(stats::median(sub$composite_score), 3),
      pct_good = round(100 * mean(sub$verdict == "Good"), 1),
      pct_moderate = round(100 * mean(sub$verdict == "Moderate"), 1),
      pct_poor = round(100 * mean(sub$verdict == "Poor"), 1),
      stringsAsFactors = FALSE
    )
  })
  summary_tbl <- .bind_rows(summary_rows)

  cat("\n")
  print(summary_tbl, n = n_enz)
  cat(sprintf(
    "\n%s rows total (%s proteins x %s enzymes).\n",
    format(n_prot * n_enz, big.mark = ","),
    format(n_prot, big.mark = ","),
    n_enz
  ))
  invisible(x)
}

#' Triage proteins from a batch evaluation
#'
#' `triage_proteins()` appends an `action` column to the flat tibble returned
#' by [batch_evaluate()] with deterministic recommendations based on each
#' protein's verdict and difficulty flags.
#'
#' @param batch_result A tibble returned by [batch_evaluate()]. If `NULL` or
#'   empty, raises an error.
#'
#' @return A tibble with one row per protein containing all score and flag
#'   columns from the flat batch summary, plus an `action` column. Possible
#'   values:
#'   \describe{
#'     \item{\code{"proceed"}}{Good verdict. No intervention indicated.}
#'     \item{\code{"consider_alternative"}}{Moderate verdict with at least one
#'       component score below 0.5. A preset change or missed-cleavage
#'       adjustment may improve results.}
#'     \item{\code{"try_other_enzyme"}}{Moderate or Poor verdict with
#'       \code{flag_hydrophobic} or \code{flag_short_protein}, or any Poor
#'       verdict without an intrinsic complexity flag. An alternative enzyme
#'       is the most likely improvement path.}
#'     \item{\code{"skip"}}{No valid peptides or low-complexity sequence. No
#'       standard enzyme choice is expected to substantially improve the
#'       score.}
#'   }
#'
#' @section Limitations:
#' Triage actions are advisory, based on heuristic difficulty flags from the
#' batch score columns. They do not re-evaluate the protein with alternative
#' enzymes. Use [compare_digests()] for that.
#'
#' @family evaluation
#'
#' @seealso [batch_evaluate()], [summarize_batch()]
#'
#' @examples
#' small_path <- system.file(
#'   "extdata", "small_proteome_50_proteins.fasta",
#'   package = "pepVet"
#' )
#' batch <- batch_evaluate(small_path, enzyme = "trypsin")
#' triaged <- triage_proteins(batch)
#' table(triaged$action)
#' @export
# nolint start: object_usage_linter.
triage_proteins <- function(batch_result) {
  .validate_batch_result(batch_result)
  flat <- batch_result

  score_cols <- intersect(
    names(flat),
    c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge", "S_unique")
  )

  ## Fully vectorized: avoids row-by-row subsetting via vapply.
  ## any_component_low: TRUE when at least one component score < 0.5.
  any_component_low <- if (length(score_cols) > 0L) {
    rowSums(as.matrix(flat[score_cols]) < 0.5, na.rm = TRUE) > 0L
  } else {
    rep(FALSE, nrow(flat))
  }

  action <- ifelse(
    flat$verdict == "Good",
    "proceed",
    ifelse(
      flat$flag_no_valid_peptides | flat$flag_low_complexity,
      "skip",
      ifelse(
        flat$flag_hydrophobic | flat$flag_short_protein |
          flat$verdict == "Poor",
        "try_other_enzyme",
        ifelse(
          any_component_low,
          "consider_alternative",
          "consider_alternative"
        )
      )
    )
  )

  tibble::add_column(flat, action = action)
}
# nolint end


## Weight sensitivity analysis

#' Weight sensitivity analysis
#'
#' `sensitivity_analysis()` perturbs the default scoring weights using a
#' Dirichlet distribution and reports how often the verdict (or enzyme ranking)
#' changes. This answers the question: "If the weights were slightly different,
#' would the recommendation change?"
#'
#' @param x An [evaluate_digest()], [compare_digests()], or [batch_evaluate()]
#'   result.
#' @param nu Dirichlet concentration scaling factor.  Controls how much the
#'   perturbed weights are allowed to vary from the defaults.  Default `63`
#'   gives a standard deviation of approximately 0.05 for a weight of 0.2.
#'   Larger values produce tighter distributions; smaller values allow more
#'   variation.
#' @param n_iter Number of Monte Carlo iterations.  Default `10000L`.
#' @param chunk_size Number of proteins to process per chunk in batch mode.
#'   Default `10000L`.  Set lower on memory-constrained machines.
#' @param importance Logical.  If `TRUE`, compute the squared Pearson
#'   correlation (R squared) between each z-scored weight and the composite
#'   score across iterations, indicating which weight drives the most variance.
#' @param corner_cases Logical.  If `TRUE`, report the composite score when
#'   each weight is at its 95 percent Dirichlet interval bound (others held at
#'   default and renormalised).
#'
#' @return A list.  For single-protein input: `iterations` (tibble of per-draw
#'   weights and composites), `convergence` (cumulative stability trace),
#'   `summary` (probabilistic verdict, composite CI, optional R squared values
#'   and corner-case table).  For multi-enzyme input: additionally reports
#'   `top1_stability` and Kendall rank correlation.  For batch input: returns
#'   `per_protein` (instability per protein) and `summary` aggregates.
#'
#' @section Limitations:
#' Monte Carlo estimates are approximate. Stability depends on `n_iter`.
#' The analysis only perturbs weights within the Dirichlet framework; it does
#' not test alternative scoring functions or parameter ranges.
#'
#' @family evaluation
#'
#' @examples
#' bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
#' res <- evaluate_digest(bsa_path, enzyme = "trypsin")
#' sens <- sensitivity_analysis(res, n_iter = 1000L)
#' sens$summary$verdict_pct
#'
#' @export
sensitivity_analysis <- function(x, nu = 63, n_iter = 10000L,
                                 chunk_size = 10000L,
                                 importance = FALSE,
                                 corner_cases = FALSE) {
  if (is.data.frame(x)) {
    if ("enzyme" %in% names(x) && length(unique(x$protein_id)) == 1L) {
      .sensitivity_enzymes(x, nu, n_iter, importance, corner_cases)
    } else {
      .sensitivity_batch(x, nu, n_iter, chunk_size, importance)
    }
  } else if (is.list(x) && "scores" %in% names(x)) {
    .sensitivity_single(
      x$scores, x$params, nu, n_iter, importance, corner_cases
    )
  } else {
    .abort(
      paste0("{.arg x} must be an {.fn evaluate_digest}, ",
        "{.fn compare_digests}, or {.fn batch_evaluate} result."),
      class = "pepvet_error_invalid_input"
    )
  }
}


## Single-protein sensitivity

.sensitivity_single <- function(scores, params, nu, n_iter,
                                importance, corner_cases) {
  score_row <- as.list(scores[1L, , drop = FALSE])
  has_unique <- "S_unique" %in% names(score_row)
  w0 <- if (has_unique) {
    .default_scoring_weights$proteome_aware
  } else {
    .default_scoring_weights$protein_only
  }

  comp_names <- names(w0)
  s_vec <- unlist(score_row[comp_names], use.names = FALSE)

  alpha <- nu * w0
  W <- .rdirichlet(as.integer(n_iter), alpha)
  colnames(W) <- comp_names

  composites <- drop(W %*% s_vec)
  verdicts <- ifelse(
    composites >= .get_param("verdict_good"), "Good",
    ifelse(composites >= .get_param("verdict_moderate"), "Moderate", "Poor")
  )

  default_comp <- sum(w0 * s_vec)
  default_verdict <- ifelse(
    default_comp >= .get_param("verdict_good"), "Good",
    ifelse(default_comp >= .get_param("verdict_moderate"), "Moderate", "Poor")
  )

  iter_df <- as.data.frame(W)
  iter_df$iteration <- seq_len(n_iter)
  iter_df$composite_score <- composites
  iter_df$verdict <- verdicts
  iter_df <- tibble::as_tibble(iter_df)

  stab <- cumsum(verdicts == default_verdict) / seq_len(n_iter)
  conv <- tibble::tibble(
    iteration = seq_len(n_iter),
    cumulative_stability = stab,
    mc_se = sqrt(stab * (1 - stab) / seq_len(n_iter))
  )

  vtab <- table(verdicts)
  verdict_pct <- as.numeric(prop.table(vtab))
  names(verdict_pct) <- names(vtab)
  qi <- stats::quantile(composites, c(0.025, 0.975), na.rm = TRUE)

  out <- list(
    iterations = iter_df,
    convergence = conv,
    summary = list(
      verdict_pct      = verdict_pct,
      composite_ci     = unname(qi),
      composite_mean   = mean(composites, na.rm = TRUE)
    )
  )

  if (importance) {
    zW <- scale(W)
    r2 <- vapply(seq_len(ncol(zW)), function(j) {
      stats::cor(zW[, j], composites)^2
    }, numeric(1))
    names(r2) <- comp_names
    out$summary$weight_importance <- r2
  }

  if (corner_cases) {
    half_span <- sqrt(w0 * (1 - w0) / (nu + 1)) * stats::qnorm(0.975)
    lo <- pmax(0, w0 - half_span)
    hi <- pmin(1, w0 + half_span)
    cc <- data.frame(
      weight = comp_names,
      default = w0,
      lo = lo,
      hi = hi,
      composite_at_lo = NA_real_,
      composite_at_hi = NA_real_,
      stringsAsFactors = FALSE
    )
    for (i in seq_along(comp_names)) {
      w_lo <- w0
      w_lo[i] <- lo[i]
      w_lo <- w_lo / sum(w_lo)
      cc$composite_at_lo[i] <- sum(w_lo * s_vec)

      w_hi <- w0
      w_hi[i] <- hi[i]
      w_hi <- w_hi / sum(w_hi)
      cc$composite_at_hi[i] <- sum(w_hi * s_vec)
    }
    out$summary$corner_cases <- tibble::as_tibble(cc)
  }

  out
}


## Multi-enzyme sensitivity

.sensitivity_enzymes <- function(enzyme_tbl, nu, n_iter,
                                 importance, corner_cases) {
  enzymes <- unique(enzyme_tbl$enzyme)
  res_list <- lapply(enzymes, function(enz) {
    sub <- enzyme_tbl[enzyme_tbl$enzyme == enz, , drop = FALSE]
    params <- list(enzyme = enz, protein_ids = sub$protein_id[[1L]])
    .sensitivity_single(sub, params, nu, n_iter, importance, corner_cases)
  })

  n_enz <- length(enzymes)
  rank_matrix <- matrix(NA_real_, nrow = n_iter, ncol = n_enz)
  for (k in seq_len(n_enz)) {
    rank_matrix[, k] <- res_list[[k]]$iterations$composite_score
  }
  top1 <- apply(rank_matrix, 1, function(row) enzymes[which.max(row)])
  top1_stab <- prop.table(table(top1))

  default_composites <- vapply(
    res_list, function(r) r$summary$composite_mean,
    numeric(1)
  )
  default_rank <- order(default_composites, decreasing = TRUE)
  kendalls <- apply(rank_matrix, 1, function(row) {
    iter_rank <- order(row, decreasing = TRUE)
    stats::cor(default_rank, iter_rank, method = "kendall")
  })
  kendall_mean <- mean(kendalls, na.rm = TRUE)

  list(
    summary = list(
      top1_stability = top1_stab,
      kendall_mean   = kendall_mean
    )
  )
}


## Batch sensitivity

.sensitivity_batch <- function(batch_tbl, nu, n_iter, chunk_size,
                               importance) {
  has_unique <- "S_unique" %in% names(batch_tbl)
  w0 <- if (has_unique) {
    .default_scoring_weights$proteome_aware
  } else {
    .default_scoring_weights$protein_only
  }
  comp_names <- names(w0)

  S <- as.matrix(batch_tbl[, comp_names, drop = FALSE])
  n_prot <- nrow(S)

  alpha <- nu * w0
  W <- .rdirichlet(as.integer(n_iter), alpha)

  n_chunks <- ceiling(n_prot / chunk_size)
  chunk_starts <- seq(1L, n_prot, by = chunk_size)

  per_protein_list <- vector("list", n_chunks)

  for (ci in seq_len(n_chunks)) {
    idx <- chunk_starts[ci]:min(chunk_starts[ci] + chunk_size - 1L, n_prot)
    S_chunk <- S[idx, , drop = FALSE]
    C_chunk <- S_chunk %*% t(W)

    def <- drop(S_chunk %*% w0)
    def_verdict <- ifelse(
      def >= .get_param("verdict_good"), "Good",
      ifelse(def >= .get_param("verdict_moderate"), "Moderate", "Poor")
    )

    ## Compare each iteration verdict with the default verdict using full 3-level
    ## classification (not just Good vs not-Good like the previous
    ## single-threshold check).
    iter_good <- C_chunk >= .get_param("verdict_good")
    iter_mod  <- C_chunk >= .get_param("verdict_moderate")
    iter_verdict <- ifelse(iter_good, "Good",
      ifelse(iter_mod, "Moderate", "Poor"))
    same <- iter_verdict == def_verdict
    instability <- 1 - rowMeans(same, na.rm = TRUE)

    ## Empirical quantiles (one pass per row for both tails)
    ci_mat <- apply(C_chunk, 1L, stats::quantile,
      probs = c(0.025, 0.975), na.rm = TRUE)
    ci_lo <- ci_mat[1L, ]
    ci_hi <- ci_mat[2L, ]
    comp_mean <- rowMeans(C_chunk, na.rm = TRUE)

    per_protein_list[[ci]] <- tibble::tibble(
      protein_id          = batch_tbl$protein_id[idx],
      verdict_instability = instability,
      composite_mean      = comp_mean,
      composite_lo        = ci_lo,
      composite_hi        = ci_hi
    )

    ## Keep composites for importance calculation (shared across chunks)
    if (importance) {
      if (ci == 1L) {
        C_all <- C_chunk
      } else {
        C_all <- rbind(C_all, C_chunk)
      }
    }

    rm(C_chunk, same)
  }

  per_protein <- .bind_rows(per_protein_list)

  ## Preserve extra grouping columns from the input (e.g. enzyme)
  extra_cols <- setdiff(names(batch_tbl), names(per_protein))
  for (ec in extra_cols) {
    per_protein[[ec]] <- batch_tbl[[ec]]
  }

  total_inst <- mean(per_protein$verdict_instability, na.rm = TRUE)
  ci_widths <- per_protein$composite_hi - per_protein$composite_lo

  out <- list(
    per_protein = per_protein,
    summary = list(
      total_instability = total_inst,
      mean_ci_width = mean(ci_widths, na.rm = TRUE),
      ci_width_quantiles = stats::quantile(ci_widths,
        c(0, 0.25, 0.5, 0.75, 1),
        na.rm = TRUE
      )
    )
  )

  if (importance) {
    ## C_all is n_prot x n_iter, W is n_iter x n_comp
    r2_per_prot <- vapply(seq_len(n_prot), function(i) {
      vapply(seq_along(comp_names), function(j) {
        stats::cor(W[, j], C_all[i, ])^2
      }, numeric(1))
    }, numeric(length(comp_names)))

    r2_mean <- rowMeans(r2_per_prot, na.rm = TRUE)
    names(r2_mean) <- comp_names
    out$summary$weight_importance <- r2_mean
  }

  out
}
