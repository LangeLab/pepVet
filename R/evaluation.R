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
#' in `sequences` and returns a flat tibble with one row per protein. Columns
#' include `protein_id`, `protein_length`, all component scores, `composite_score`,
#' `verdict`, `n_peptides`, `n_valid_peptides`, `median_peptide_length`, and four
#' sequence-level difficulty flags. Pass the result to [summarize_batch()] for
#' aggregate statistics or to [triage_proteins()] for action recommendations.
#'
#' @param sequences Multi-protein input. Accepts the same forms as
#'   [digest_protein()]. Must resolve to at least one protein.
#' @param enzyme Enzyme name passed to [digest_protein()]. Defaults to
#'   `"trypsin"`.
#' @param missed_cleavages Maximum missed cleavages. Defaults to `0L`.
#' @param include_cleavage_efficiency Logical flag passed to
#'   [evaluate_digest()] and ultimately [digest_protein()]. When `TRUE`, each
#'   per-protein peptide table includes a `cleavage_efficiency` column (does
#'   not affect the flat batch tibble columns).
#' @param proteome Optional proteome digest tibble passed to [score_peptides()]
#'   for every protein evaluation. When supplied, an `S_unique` column appears
#'   in the returned tibble.
#' @param weights Optional scoring weight vector passed to [score_peptides()].
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
                           missed_cleavages = 0L,
                           include_cleavage_efficiency = FALSE,
                           proteome = NULL,
                           weights = NULL,
                           ...) {
  normalized_input <- .read_input(sequences)
  protein_ids <- names(normalized_input)

  # Digest all proteins in one call — eliminates 20K individual evaluate_digest
  # invocations and the extra .read_input() overhead inside each one.
  all_peptides <- digest_protein(
    normalized_input,
    enzyme = enzyme,
    missed_cleavages = missed_cleavages,
    include_cleavage_efficiency = include_cleavage_efficiency
  )

  # Score all proteins in one call — eliminates 20K separate score_peptides
  # setup overheads and the .cleavage_efficiency_summary() call that
  # evaluate_digest runs per protein even though batch_evaluate discards it.
  all_scores <- score_peptides(
    all_peptides,
    proteome = proteome,
    weights  = weights,
    ...,
    enzyme = enzyme
  )

  score_cols <- intersect(
    names(all_scores),
    c(
      "S_length", "S_coverage", "S_count", "S_hydro", "S_charge",
      "S_unique", "composite_score", "verdict", "median_peptide_length"
    )
  )

  # Compute all difficulty flags across all proteins at once — eliminates
  # 20K per-protein tibble subsets, 20K .calculate_gravy_vec() calls, and
  # 20K tibble::as_tibble() constructions.
  flags <- .batch_difficulty_flags(all_peptides, protein_ids)

  # Reorder scores to match protein_ids order and extract score_cols.
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
# nolint end

# ---- Private batch helpers ----

.validate_batch_result <- function(batch_result) {
  if (!inherits(batch_result, "data.frame")) {
    cli::cli_abort(
      "{.arg batch_result} must be a tibble returned by {.fn batch_evaluate}.",
      class = "pepvet_error_invalid_batch_result"
    )
  }

  if (nrow(batch_result) == 0L) {
    cli::cli_abort(
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
    cli::cli_abort(
      c(
        "{.arg batch_result} is missing required columns from {.fn batch_evaluate}.",
        "i" = paste("Missing:", paste(missing_cols, collapse = ", "))
      ),
      class = "pepvet_error_invalid_batch_result"
    )
  }

  batch_result
}

# Vectorized difficulty flags for a full multi-protein peptide tibble.
# Returns a named list of vectors, each of length == length(protein_ids),
# in the same order as protein_ids.
.batch_difficulty_flags <- function(all_peptides, protein_ids) {
  pid_factor <- factor(all_peptides$protein_id, levels = protein_ids)

  # protein_length and n_peptides
  protein_length   <- as.integer(tapply(all_peptides$end, pid_factor, max))
  n_peptides       <- as.integer(tabulate(pid_factor))

  # valid peptide mask (length 7–25)
  valid_mask       <- all_peptides$length >= 7L & all_peptides$length <= 25L
  n_valid_peptides <- as.integer(tabulate(pid_factor[valid_mask]))

  # flags derivable from counts
  flag_short_protein    <- protein_length < 100L
  flag_no_valid_peptides <- n_valid_peptides == 0L

  # flag_hydrophobic: median GRAVY of valid peptides > 0.6 per protein
  flag_hydrophobic <- logical(length(protein_ids))
  if (any(valid_mask)) {
    gravy_vals   <- .calculate_gravy_vec(all_peptides$peptide[valid_mask])
    median_gravy <- tapply(gravy_vals, pid_factor[valid_mask], stats::median)
    flag_hydrophobic[match(names(median_gravy), protein_ids)] <- median_gravy > 0.6
  }

  # flag_low_complexity: dominant AA > 50% in reconstructed MC=0 sequence.
  # Build one concatenated sequence per protein from MC=0 peptides (sorted by
  # start), then check character frequencies.
  mc0_mask  <- all_peptides$missed_cleavages == 0L
  mc0_peps  <- all_peptides[mc0_mask, c("protein_id", "start", "peptide"),
                             drop = FALSE]
  mc0_peps  <- mc0_peps[order(mc0_peps$protein_id, mc0_peps$start), ,
                         drop = FALSE]
  prot_seqs <- tapply(mc0_peps$peptide, mc0_peps$protein_id,
                      paste, collapse = "", simplify = FALSE)

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
    protein_length        = protein_length,
    n_peptides            = n_peptides,
    n_valid_peptides      = n_valid_peptides,
    flag_short_protein    = flag_short_protein,
    flag_no_valid_peptides = flag_no_valid_peptides,
    flag_hydrophobic      = flag_hydrophobic,
    flag_low_complexity   = unname(flag_low_complexity)
  )
}

.compute_difficulty_flags <- function(peptides) {
  protein_length <- max(peptides$end)
  valid_mask <- peptides$length >= 7L & peptides$length <= 25L
  n_valid <- sum(valid_mask)
  valid_peps <- peptides[valid_mask, , drop = FALSE]

  flag_short_protein <- protein_length < 100L
  flag_no_valid_peptides <- n_valid == 0L

  flag_hydrophobic <- if (n_valid > 0L) {
    gravy_vals <- .calculate_gravy_vec(valid_peps$peptide)
    stats::median(gravy_vals) > 0.6
  } else {
    FALSE
  }

  mc0_peps <- peptides[peptides$missed_cleavages == 0L, , drop = FALSE]
  mc0_peps <- mc0_peps[order(mc0_peps$start), , drop = FALSE]
  protein_seq <- paste(mc0_peps$peptide, collapse = "")
  flag_low_complexity <- if (nchar(protein_seq) > 0L) {
    aa_counts <- table(strsplit(protein_seq, "", fixed = TRUE)[[1]])
    max(aa_counts) / sum(aa_counts) > 0.5
  } else {
    FALSE
  }

  list(
    protein_length        = protein_length,
    n_valid_peptides      = n_valid,
    flag_short_protein    = flag_short_protein,
    flag_no_valid_peptides = flag_no_valid_peptides,
    flag_hydrophobic      = flag_hydrophobic,
    flag_low_complexity   = flag_low_complexity
  )
}

#' Summarize a Batch Digest Evaluation
#'
#' `summarize_batch()` extracts aggregate statistics from a [batch_evaluate()]
#' result tibble. It returns a named list covering verdict distribution, score
#' distribution, per-component averages, the lowest-scoring proteins, and a
#' heuristic set of enzyme-switch candidates.
#'
#' @param batch_result A tibble returned by [batch_evaluate()].
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{`verdict_counts`}{A tibble with columns `verdict`, `n`, and `pct`
#'       covering the three verdict categories.}
#'     \item{`score_distribution`}{A named numeric vector with `mean`,
#'       `median`, `sd`, `q25`, `q75`, `min`, and `max` of composite scores.}
#'     \item{`component_summary`}{A named numeric vector of per-component mean
#'       scores. The lowest values identify the weakest scoring dimension
#'       across the proteome.}
#'     \item{`problem_proteins`}{A tibble of proteins in the bottom 10% by
#'       composite score, ordered ascending, with all score and flag columns.}
#'     \item{`enzyme_switch_candidates`}{A tibble of Moderate or Poor proteins
#'       where `flag_hydrophobic` or `flag_short_protein` is `TRUE`, indicating
#'       that enzyme or preset choice is the likely limiting factor.}
#'   }
#'
#' @details `enzyme_switch_candidates` is a heuristic flag list derived from
#'   sequence-level difficulty flags, not from running alternative enzymes.
#'   Use [compare_digests()] to confirm whether a specific alternative enzyme
#'   improves the verdict for a flagged protein.
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
    order(problem_proteins$composite_score), , drop = FALSE
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

#' Triage Proteins from a Batch Evaluation
#'
#' `triage_proteins()` appends an `action` column to the flat tibble returned
#' by [batch_evaluate()] with deterministic recommendations based on each
#' protein's verdict and difficulty flags.
#'
#' @param batch_result A tibble returned by [batch_evaluate()].
#'
#' @return A tibble with one row per protein containing all score and flag
#'   columns from the flat batch summary, plus an `action` column. Possible
#'   values:
#'   \describe{
#'     \item{`"proceed"`}{Good verdict. No intervention indicated.}
#'     \item{`"consider_alternative"`}{Moderate verdict with at least one
#'       component score below 0.5. A preset change or missed-cleavage
#'       adjustment may improve results.}
#'     \item{`"try_other_enzyme"`}{Moderate or Poor verdict with
#'       `flag_hydrophobic` or `flag_short_protein`, or any Poor verdict
#'       without an intrinsic complexity flag. An alternative enzyme is the
#'       most likely improvement path.}
#'     \item{`"skip"`}{No valid peptides or low-complexity sequence. No
#'       standard enzyme choice is expected to substantially improve the
#'       score.}
#'   }
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

  # Fully vectorized: avoids row-by-row subsetting via vapply.
  # any_component_low: TRUE when at least one component score < 0.5.
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
        flat$flag_hydrophobic | flat$flag_short_protein | flat$verdict == "Poor",
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
