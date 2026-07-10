## data-raw/peptideatlas-concordance.R
## PeptideAtlas concordance analysis
##
## Tests whether pepVet-valid peptides are observed at higher rates in
## PeptideAtlas than pepVet-invalid peptides. Grid-searches length and GRAVY
## boundaries to find data-optimal ranges. Calibrates the composite threshold
## via logistic regression + Youden's J.
##
## Produces pre-computed tables in inst/extdata/peptideatlas-concordance/.
##
## Dependencies: Biostrings, tibble, pepVet
## Data: scratch/peptideatlas-human-2026-01.fasta (PeptideAtlas build)
##       scratch/human-ref.fasta (UniProt human reference proteome)

library(Biostrings)
library(tibble)
devtools::load_all(".", quiet = TRUE)

set.seed(42)
cat("=== PeptideAtlas Concordance Analysis ===\n\n")

## ---- Configuration ----------------------------------------------------------

pa_fasta     <- "scratch/peptideatlas-human-2026-01.fasta"
ref_fasta    <- "scratch/human-ref.fasta"
out_dir      <- "inst/extdata/peptideatlas-concordance"
n_proteins   <- 500L
n_bootstrap  <- 2000L

default_length_range <- c(7L, 25L)
default_gravy_range  <- c(-1.0, 0.6)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Step 1: Parse PeptideAtlas observed peptides ---------------------------
cat("Step 1: Parsing observed peptides from PeptideAtlas...\n")
pa_seqs <- readAAStringSet(pa_fasta)
all_observed <- unique(toupper(as.character(pa_seqs)))
cat(sprintf("  Total entries: %d\n", length(pa_seqs)))
cat(sprintf("  Unique peptides: %d\n", length(all_observed)))
cat(sprintf("  Length range: %d - %d\n", min(nchar(all_observed)), max(nchar(all_observed))))

rm(pa_seqs)
cat(sprintf("  %d unique peptides loaded\n\n", length(all_observed)))

## ---- Step 2: Select proteins from reference proteome ------------------------
cat("Step 2: Selecting proteins from human reference proteome...\n")
ref_seqs <- readAAStringSet(ref_fasta)

aa_regex <- "^[ACDEFGHIKLMNPQRSTVWY]+$"
seqs_char <- as.character(ref_seqs)
nchar_vec <- nchar(seqs_char)
has_standard_aa <- grepl(aa_regex, seqs_char)

valid_idx <- which(nchar_vec >= 50 & has_standard_aa)
cat(sprintf("  Proteins with length >= 50 AA and standard AAs: %d / %d\n",
            length(valid_idx), length(ref_seqs)))

selected_idx <- sample(valid_idx, min(n_proteins, length(valid_idx)))
selected_proteins <- ref_seqs[selected_idx]
cat(sprintf("  Selected %d proteins\n\n", length(selected_proteins)))

## ---- Step 3: Digest + classify + check observation --------------------------
cat("Step 3: Digesting with pepVet and checking observation status...\n")

## Bulk digest with pepVet's digest_protein (uses cleaver under the hood)
all_peptides <- digest_protein(
  selected_proteins,
  enzyme = "trypsin",
  missed_cleavages = 1L
)
cat(sprintf("  Total peptides from digest_protein: %d\n", nrow(all_peptides)))
cat(sprintf("  Unique proteins in digest: %d\n",
            length(unique(all_peptides$protein_id))))

## Bulk scoring with batch_evaluate
batch_scores <- batch_evaluate(
  selected_proteins,
  enzyme = "trypsin",
  missed_cleavages = 1L
)
cat(sprintf("  batch_evaluate returned %d proteins\n", nrow(batch_scores)))

## Add GRAVY for every peptide using pepVet's internal function
all_peptides$gravy <- .calculate_gravy(all_peptides$peptide)

## Check observation status for each peptide (fast C-level matching)
all_peptides$is_observed <- all_peptides$peptide %in% all_observed
rm(all_observed)

## ---- Aggregate per-protein results ------------------------------------------
cat("\nStep 4: Aggregating per-protein results...\n")

protein_ids <- unique(all_peptides$protein_id)
n_prot <- length(protein_ids)

## Vectorized aggregation via tapply (single pass over all peptides)
prot_idx <- match(all_peptides$protein_id, protein_ids)
valid_mask <- all_peptides$length >= default_length_range[1] &
              all_peptides$length <= default_length_range[2] &
              all_peptides$gravy >= default_gravy_range[1] &
              all_peptides$gravy <= default_gravy_range[2]

per_protein <- data.frame(
  protein_id               = protein_ids,
  n_theoretical            = as.integer(tapply(seq_len(nrow(all_peptides)), prot_idx, length)),
  n_observed               = as.integer(tapply(all_peptides$is_observed, prot_idx, sum)),
  n_pepVet_valid           = as.integer(tapply(valid_mask, prot_idx, sum)),
  n_observed_and_valid     = as.integer(
    tapply(valid_mask & all_peptides$is_observed, prot_idx, sum)
  ),
  n_observed_and_invalid   = as.integer(
    tapply(!valid_mask & all_peptides$is_observed, prot_idx, sum)
  ),
  stringsAsFactors = FALSE
)

per_protein$detection_rate_all     <- per_protein$n_observed / per_protein$n_theoretical
per_protein$detection_rate_valid   <- ifelse(per_protein$n_pepVet_valid > 0L,
  per_protein$n_observed_and_valid / per_protein$n_pepVet_valid, NA_real_)
per_protein$detection_rate_invalid <- ifelse(
  per_protein$n_theoretical > per_protein$n_pepVet_valid,
  per_protein$n_observed_and_invalid / (
    per_protein$n_theoretical - per_protein$n_pepVet_valid
  ),
  NA_real_
)

## Merge the original sequence length, score, and verdict from batch_evaluate.
## The first peptide length is not the protein length, especially for missed
## cleavage digests.
scores_sub <- batch_scores[, c(
  "protein_id", "protein_length", "composite_score", "verdict"
)]
per_protein <- merge(per_protein, scores_sub, by = "protein_id", all.x = TRUE)

## Validate the generated research artifact before writing it. These checks
## belong here because the PeptideAtlas table is a data product, not part of
## pepVet's runtime API.
if (
  nrow(per_protein) != length(selected_proteins) ||
    anyNA(per_protein$protein_length) ||
    any(per_protein$protein_length < 50L)
) {
  stop(
    paste(
      "PeptideAtlas artifact invariant failed:",
      "expected one complete protein-length record per selected protein,",
      "with all lengths >= 50 AA."
    ),
    call. = FALSE
  )
}

cat(sprintf("  Proteins with observed peptides: %d / %d (%.1f%%)\n",
            sum(per_protein$n_observed > 0), n_prot,
            100 * mean(per_protein$n_observed > 0)))

## Save intermediate results
saveRDS(per_protein, file.path(out_dir, "per-protein-results.rds"))
write.csv(per_protein, file.path(out_dir, "per-protein-results.csv"),
          row.names = FALSE)
cat("  Saved per-protein results\n")

## ---- Step 4b: Enrichment measurement with bootstrap -------------------------
cat("\nStep 4b: Measuring enrichment with bootstrap CI...\n")

valid_obs_rate   <- mean(per_protein$detection_rate_valid, na.rm = TRUE)
invalid_obs_rate <- mean(per_protein$detection_rate_invalid, na.rm = TRUE)
enrichment       <- valid_obs_rate - invalid_obs_rate

cat(sprintf("  Mean detection rate (pepVet-valid peptides):   %.4f\n", valid_obs_rate))
cat(sprintf("  Mean detection rate (pepVet-invalid peptides): %.4f\n", invalid_obs_rate))
cat(sprintf("  Enrichment (valid - invalid):                 %.4f\n", enrichment))

boot_enrichment <- numeric(n_bootstrap)
for (b in seq_len(n_bootstrap)) {
  idx <- sample(n_prot, replace = TRUE)
  boot_valid   <- mean(per_protein$detection_rate_valid[idx], na.rm = TRUE)
  boot_invalid <- mean(per_protein$detection_rate_invalid[idx], na.rm = TRUE)
  boot_enrichment[b] <- boot_valid - boot_invalid
}
ci_95 <- quantile(boot_enrichment, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("  95%% CI for enrichment: [%.4f, %.4f]\n", ci_95[1], ci_95[2]))
cat(sprintf("  Enrichment is positive: %s\n\n", if (ci_95[1] > 0) "YES" else "NO"))

enrichment_table <- tibble(
  metric = c("detection_rate_valid", "detection_rate_invalid",
             "enrichment", "ci_lower", "ci_upper"),
  value  = c(valid_obs_rate, invalid_obs_rate, enrichment, ci_95[1], ci_95[2])
)

## ---- Step 5: Grid search for optimal length/GRAVY boundaries ----------------
cat("Step 5: Grid-searching length and GRAVY boundaries...\n")

## Vectorized approach: precompute per-peptide attributes once, then evaluate
## all parameter combos in bulk using vectorized masks over all 60K+ peptides.
peptide_lengths <- all_peptides$length
peptide_gravy   <- all_peptides$gravy
peptide_observed <- all_peptides$is_observed
peptide_protein  <- all_peptides$protein_id

## Map protein_id to index in per_protein table
prot_idx <- match(peptide_protein, per_protein$protein_id)

length_lo_candidates <- c(5L, 6L, 7L, 8L, 9L, 10L)
length_hi_candidates <- c(20L, 25L, 30L, 35L, 40L)
gravy_lo_candidates  <- c(-1.5, -1.2, -1.0, -0.8, -0.5)
gravy_hi_candidates  <- c(0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

build_combo <- function(len_lo, len_hi, g_lo, g_hi) {
  valid_mask <- peptide_lengths >= len_lo & peptide_lengths <= len_hi &
    peptide_gravy >= g_lo & peptide_gravy <= g_hi

  ## Aggregate per-protein using tapply (fast, base R)
  n_valid     <- as.integer(tapply(valid_mask, prot_idx, sum))
  obs_valid   <- as.integer(tapply(valid_mask & peptide_observed, prot_idx, sum))
  obs_invalid <- as.integer(tapply(!valid_mask & peptide_observed, prot_idx, sum))
  n_invalid   <- per_protein$n_theoretical - n_valid

  dr_valid   <- ifelse(n_valid > 0, obs_valid / n_valid, NA_real_)
  dr_invalid <- ifelse(n_invalid > 0, obs_invalid / n_invalid, NA_real_)

  list(
    n_valid_total   = sum(n_valid),
    n_invalid_total = sum(n_invalid),
    dr_valid_mean   = mean(dr_valid, na.rm = TRUE),
    dr_invalid_mean = mean(dr_invalid, na.rm = TRUE),
    enrichment      = mean(dr_valid, na.rm = TRUE) - mean(dr_invalid, na.rm = TRUE)
  )
}

## Build all valid parameter combos first
param_grid <- expand.grid(
  length_min = length_lo_candidates,
  length_max = length_hi_candidates,
  gravy_min  = gravy_lo_candidates,
  gravy_max  = gravy_hi_candidates,
  stringsAsFactors = FALSE
)
param_grid <- param_grid[param_grid$length_max > param_grid$length_min &
                           param_grid$gravy_max > param_grid$gravy_min, ]

n_combos <- nrow(param_grid)
cat(sprintf("  Evaluating %d parameter combinations...\n", n_combos))

grid_list <- vector("list", n_combos)
for (i in seq_len(n_combos)) {
  res <- build_combo(
    param_grid$length_min[i], param_grid$length_max[i],
    param_grid$gravy_min[i],  param_grid$gravy_max[i]
  )
  grid_list[[i]] <- tibble(
    length_min = param_grid$length_min[i],
    length_max = param_grid$length_max[i],
    gravy_min  = param_grid$gravy_min[i],
    gravy_max  = param_grid$gravy_max[i],
    n_valid_total   = res$n_valid_total,
    n_invalid_total = res$n_invalid_total,
    detection_rate_valid   = res$dr_valid_mean,
    detection_rate_invalid = res$dr_invalid_mean,
    enrichment = res$enrichment
  )
}

grid_tbl <- .bind_rows(grid_list)
grid_tbl <- grid_tbl[order(grid_tbl$enrichment, decreasing = TRUE), ]
rownames(grid_tbl) <- NULL

cat(sprintf("  Grid search complete: %d combinations\n", nrow(grid_tbl)))
cat("  Top 10 parameter combinations by enrichment:\n")
print(head(grid_tbl, 10), row.names = FALSE)

best_params <- grid_tbl[1, ]
cat(sprintf("\n  Best parameters:\n"))
cat(sprintf("    length_range: [%d, %d]\n", best_params$length_min, best_params$length_max))
cat(sprintf("    gravy_range:  [%.1f, %.1f]\n", best_params$gravy_min, best_params$gravy_max))
cat(sprintf("    Enrichment:    %.4f\n", best_params$enrichment))
cat(sprintf("  Default parameters enrichment: %.4f\n\n", enrichment))

## Check if optimal differs from default by meaningful amount
len_diff <- abs(best_params$length_min - default_length_range[1]) +
  abs(best_params$length_max - default_length_range[2])
gravy_diff <- abs(best_params$gravy_min - default_gravy_range[1]) +
  abs(best_params$gravy_max - default_gravy_range[2])

if (len_diff > 2 || gravy_diff > 0.2) {
  cat(sprintf(
    "  *** Optimal boundaries differ from defaults (len delta=%d, gravy delta=%.1f).\n",
    len_diff, gravy_diff
  ))
  cat("      Consider updating default parameters.\n")
} else {
  cat("  Optimal boundaries are close to defaults. Defaults validated.\n")
}
cat("\n")

saveRDS(grid_tbl, file.path(out_dir, "grid-search-results.rds"))
write.csv(grid_tbl, file.path(out_dir, "grid-search-results.csv"),
          row.names = FALSE)
cat("  Saved grid search results\n\n")

## ---- Step 6: Logistic regression + Youden's J --------------------------------
cat("Step 6: Logistic regression for threshold calibration...\n")

per_protein$has_observed_peptide <- per_protein$n_observed > 0

logit_fit <- glm(
  has_observed_peptide ~ composite_score,
  data = per_protein,
  family = binomial
)
cat("  Logistic regression summary:\n")
print(summary(logit_fit)$coefficients)

## Compute Youden's J across a grid of thresholds
threshold_seq <- seq(0.1, 0.95, by = 0.005)
pred_probs <- predict(logit_fit, type = "response")
pred_class_obs <- per_protein$has_observed_peptide

youden_results <- data.frame(
  threshold = threshold_seq,
  sensitivity = NA_real_,
  specificity = NA_real_,
  youden_j = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(threshold_seq)) {
  thresh <- threshold_seq[i]
  pred <- pred_probs >= thresh

  tp <- sum(pred & pred_class_obs)
  fn <- sum(!pred & pred_class_obs)
  fp <- sum(pred & !pred_class_obs)
  tn <- sum(!pred & !pred_class_obs)

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_

  youden_results$sensitivity[i] <- sensitivity
  youden_results$specificity[i] <- specificity
  youden_results$youden_j[i] <- sensitivity + specificity - 1
}

best_idx <- which.max(youden_results$youden_j)
best_threshold <- youden_results$threshold[best_idx]

cat(sprintf("\n  Youden-optimal threshold: %.3f\n", best_threshold))
cat(sprintf("    Sensitivity: %.3f\n", youden_results$sensitivity[best_idx]))
cat(sprintf("    Specificity: %.3f\n", youden_results$specificity[best_idx]))
cat(sprintf("    Youden's J:  %.3f\n", youden_results$youden_j[best_idx]))
cat(sprintf("  Current default threshold (Good): 0.65\n"))
cat(sprintf("  Difference: %.3f\n", best_threshold - 0.65))

if (abs(best_threshold - 0.65) < 0.05) {
  cat("  Default threshold validated (within 0.05 of Youden-optimal).\n")
} else {
  cat("  *** Default threshold differs from Youden-optimal by > 0.05.\n")
  cat("      Consider updating .pepvet_params$verdict_good.\n")
}

saveRDS(youden_results, file.path(out_dir, "threshold-calibration.rds"))
write.csv(youden_results, file.path(out_dir, "threshold-calibration.csv"),
          row.names = FALSE)
cat("  Saved threshold calibration results\n\n")

## ---- Summary -----------------------------------------------------------------
cat("========================================\n")
cat("   PEPTIDEATLAS CONCORDANCE SUMMARY\n")
cat("========================================\n")
cat(sprintf("  Proteins analyzed:           %d\n", n_prot))
cat(sprintf("  Total theoretical peptides:  %d\n", sum(per_protein$n_theoretical)))
cat(sprintf("  Proteins w/ observed peptides: %d (%.1f%%)\n",
            sum(per_protein$has_observed_peptide),
            100 * mean(per_protein$has_observed_peptide)))
cat(sprintf("  Default enrichment:          %.4f [%.4f, %.4f]\n",
            enrichment, ci_95[1], ci_95[2]))
cat(sprintf("  Best length_range:           [%d, %d]\n",
            best_params$length_min, best_params$length_max))
cat(sprintf("  Best gravy_range:            [%.1f, %.1f]\n",
            best_params$gravy_min, best_params$gravy_max))
cat(sprintf("  Best enrichment:             %.4f\n", best_params$enrichment))
cat(sprintf("  Youden-optimal threshold:    %.3f\n", best_threshold))
cat(sprintf("  Default Good threshold:      0.65\n"))
cat(sprintf("\n  Output files in: %s/\n", out_dir))
cat("========================================\n")
