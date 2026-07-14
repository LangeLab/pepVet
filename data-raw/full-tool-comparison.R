## data-raw/full-tool-comparison.R
## Descriptive tool comparison for package documentation.
## Retained outputs and documented interfaces are compared where available.
## Output: inst/extdata/comparison-data/ with CSVs, vignettes/ with figures.
##
## pepVet stands for "pre-acquisition proteolytic digest evaluation".
## The comparison is not an experimental validation or performance benchmark.

library(ggplot2)
library(patchwork)
pkgload::load_all(quiet = TRUE)

out_dir <- "inst/extdata/comparison-data"
fig_dir <- "vignettes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Protein panel
proteins <- list(
  BSA = list(
    path = system.file("extdata", "P02769.fasta", package = "pepVet"),
    name = "Serum albumin", species = "Bos taurus", len = 607
  ),
  H3 = list(
    path = system.file("extdata", "P68431.fasta", package = "pepVet"),
    name = "Histone H3.1", species = "Homo sapiens", len = 136
  ),
  BACE1 = list(
    path = system.file("extdata", "P56817.fasta", package = "pepVet"),
    name = "Beta-secretase 1", species = "Homo sapiens", len = 501
  ),
  LYSO = list(
    path = system.file("extdata", "P00698.fasta", package = "pepVet"),
    name = "Lysozyme C", species = "Gallus gallus", len = 147
  ),
  UBIQ = list(
    path = system.file("extdata", "P0CG48.fasta", package = "pepVet"),
    name = "Ubiquitin-40S fusion", species = "Homo sapiens", len = 685
  )
)

enzymes <- c(
  "trypsin", "lysc", "chymotrypsin-high",
  "glutamyl endopeptidase", "asp-n endopeptidase"
)
presets <- c(
  "standard", "dia", "targeted", "membrane",
  "ffpe_degraded", "fractionated"
)

cat("============================================================\n")
cat("FULL TOOL COMPARISON: pepVet vs MS-Digest vs ExPASy\n")
cat("                       vs PeptideRanger vs Protein Cleaver vs ProteaseGuru\n")
cat("============================================================\n\n")

## ======================================================================
## SECTION A: BASELINE OVERLAP
## Retained peptide outputs for one matched digest configuration.
## ======================================================================
cat("Section A: Baseline overlap of retained peptide lists\n\n")

pepvet_results <- list()
for (pn in names(proteins)) {
  path <- proteins[[pn]]$path
  for (enz in enzymes) {
    ev <- evaluate_digest(path, enzyme = enz, missed_cleavages = 1L)
    pepvet_results[[paste(pn, enz, sep = "_")]] <- list(
      peptides = ev$peptides,
      scores = ev$scores
    )
  }
}

msd_raw <- read.csv("inst/extdata/tool-data/tool-compare-msdigest.csv",
  stringsAsFactors = FALSE
)
exp_raw <- read.csv("inst/extdata/tool-data/tool-compare-expasy.csv",
  stringsAsFactors = FALSE
)
pvp_raw <- read.csv("inst/extdata/tool-data/pepvet-bsa-trypsin.csv",
  stringsAsFactors = FALSE
)

msd_peps <- unique(msd_raw$peptide)
exp_peps <- unique(exp_raw$peptide)
pvp_peps <- unique(pvp_raw$peptide)

overlap <- data.frame(
  Tool = c("MS-Digest", "ExPASy PeptideMass", "pepVet"),
  N_unique = c(length(msd_peps), length(exp_peps), length(pvp_peps)),
  Common_with_pepVet = c(
    sum(msd_peps %in% pvp_peps),
    sum(exp_peps %in% pvp_peps),
    length(pvp_peps)
  ),
  Shared_all_three = sum(msd_peps %in% pvp_peps & msd_peps %in% exp_peps)
)
overlap$Pct_overlap <- round(overlap$Common_with_pepVet / overlap$N_unique * 100, 1)
cat("BSA + trypsin MC=1 peptide overlap:\n")
print(overlap, row.names = FALSE)
cat("\nThese counts describe only the retained BSA configuration.\n\n")

write.csv(overlap, file.path(out_dir, "sectionA-overlap.csv"), row.names = FALSE)

## ======================================================================
## SECTION B: TOOL CAPABILITY CATALOG
## What each tool adds beyond the peptide list.
## ======================================================================
cat("Section B: Tool capability catalog\n\n")

capabilities <- data.frame(
  Capability = c(
    "Peptide list + masses",
    "pI output",
    "Hydrophobicity metric",
    "Composite digest score",
    "Score-band label",
    "Multi-enzyme comparison",
    "Enzyme ranking",
    "Workflow presets",
    "Batch input",
    "Sensitivity analysis",
    "Peptide-level ML score",
    "3D structure mapping",
    "Retention time prediction",
    "Skyline/FASTA export",
    "Non-GUI interface",
    "R package"
  ),
  pepVet = c(
    "Yes", "Peptide", "GRAVY", "Weighted model score", "Heuristic bands",
    "Yes", "Highest model score", "6 presets", "Yes", "Weight perturbation",
    "No", "No", "No", "Yes", "Yes", "Yes"
  ),
  MS_Digest = c(
    "Yes", "No", "Bull-Breese / HPLC index", "No", "No",
    "One rule per run", "No", "No", "Multiple proteins", "Not assessed",
    "No", "No", "No", "No", "Command line / XML", "No"
  ),
  ExPASy = c(
    "Yes", "Protein", "No", "No", "No",
    "One rule per run", "No", "No", "Single protein", "Not assessed",
    "No", "No", "No", "No", "No", "No"
  ),
  Protein_Cleaver = c(
    "Yes", "No", "No", "No", "No",
    "Yes", "Application ranking", "No", "Yes", "Not assessed",
    "No (rules-based flag)", "Yes", "No documented output", "No",
    "No documented API", "No"
  ),
  ProteaseGuru = c(
    "Yes", "No", "SSRCalc", "No", "No",
    "Yes", "Comparison output", "No", "Protein database", "Not assessed",
    "No", "No", "No", "No", "No", "No"
  ),
  PeptideRanger = c(
    "No digest output", "No", "Not an output", "No", "No",
    "No", "No", "No", "Peptide vectors", "Not assessed",
    "RF score (0-1)", "No", "No", "No", "Yes", "Yes"
  )
)
cat("Tool capability matrix:\n\n")
for (i in seq_len(nrow(capabilities))) {
  cat(sprintf("  %-35s\n", capabilities$Capability[i]))
  cat(sprintf("    pepVet:          %s\n", capabilities$pepVet[i]))
  cat(sprintf("    MS-Digest:        %s\n", capabilities$MS_Digest[i]))
  cat(sprintf("    ExPASy:           %s\n", capabilities$ExPASy[i]))
  cat(sprintf("    Protein Cleaver:  %s\n", capabilities$Protein_Cleaver[i]))
  cat(sprintf("    ProteaseGuru:     %s\n", capabilities$ProteaseGuru[i]))
  cat(sprintf("    PeptideRanger:    %s\n", capabilities$PeptideRanger[i]))
}
write.csv(capabilities, file.path(out_dir, "sectionB-capabilities.csv"),
  row.names = FALSE
)

## ======================================================================
## SECTION C: pepVet SCORING MODEL
## Multi-component quality score across proteins and enzymes.
## ======================================================================
cat("\nSection C: pepVet scoring model across all proteins\n\n")

score_rows <- list()
for (pn in names(proteins)) {
  for (enz in enzymes) {
    r <- pepvet_results[[paste(pn, enz, sep = "_")]]
    s <- r$scores
    score_rows[[length(score_rows) + 1]] <- data.frame(
      protein = pn,
      enzyme = enz,
      n_total = nrow(r$peptides),
      n_valid = sum(r$peptides$length >= 7 & r$peptides$length <= 25),
      S_length = s$S_length,
      S_coverage = s$S_coverage,
      S_count = s$S_count,
      S_hydro = s$S_hydro,
      S_charge = s$S_charge,
      composite = s$composite_score,
      verdict = s$verdict,
      stringsAsFactors = FALSE
    )
  }
}
scores_tbl <- .bind_rows(score_rows)
cat("pepVet scores across 5 proteins x 5 enzymes:\n\n")
print(scores_tbl, row.names = FALSE)

write.csv(scores_tbl, file.path(out_dir, "sectionC-scores-all.csv"),
  row.names = FALSE
)

cat("\nHighest model score per protein:\n")
for (pn in names(proteins)) {
  sub <- scores_tbl[scores_tbl$protein == pn, ]
  best <- sub[which.max(sub$composite), ]
  cat(sprintf(
    "  %-6s -> %-25s (composite=%.3f, %s)\n",
    pn, best$enzyme, best$composite, best$verdict
  ))
}

## ======================================================================
## SECTION C2: WORKFLOW PRESET EFFECT
## How presets change scores on different protein types.
## ======================================================================
cat("\nSection C2: Workflow preset effects\n\n")
preset_sequences <- lapply(
  proteins[c("BSA", "H3", "BACE1")],
  function(protein) Biostrings::readAAStringSet(protein$path)
)
preset_background <- digest_protein(
  do.call(c, unname(preset_sequences)),
  enzyme = "trypsin",
  missed_cleavages = 1L
)
preset_rows <- list()
for (pn in c("BSA", "H3", "BACE1")) {
  path <- proteins[[pn]]$path
  for (pr_name in presets) {
    pr <- pepvet_preset(pr_name)
    ev <- do.call(
      evaluate_digest,
      c(
        list(
          sequence = path,
          enzyme = "trypsin",
          missed_cleavages = 1L,
          proteome = preset_background
        ),
        pr
      )
    )
    lr <- pr$length_range
    preset_rows[[length(preset_rows) + 1]] <- data.frame(
      protein = pn,
      preset = pr_name,
      composite = ev$scores$composite_score,
      verdict = ev$scores$verdict,
      n_total = nrow(ev$peptides),
      n_length_valid = sum(
        ev$peptides$length >= lr[1] & ev$peptides$length <= lr[2]
      ),
      S_hydro = ev$scores$S_hydro,
      S_unique = ev$scores$S_unique,
      preset_used = ev$scores$preset_used,
      include_pI = ev$params$include_pI,
      len_lo = lr[1],
      len_hi = lr[2],
      stringsAsFactors = FALSE
    )
  }
}
preset_tbl <- .bind_rows(preset_rows)
print(preset_tbl, row.names = FALSE)
write.csv(preset_tbl, file.path(out_dir, "sectionC2-presets.csv"),
  row.names = FALSE
)

## ======================================================================
## SECTION D: pepVet + PeptideRanger PIPELINE
## Descriptive comparison of a selected window and an ML score.
## ======================================================================
cat("\nSection D: pepVet + PeptideRanger pipeline\n\n")

peptideranger_output <- file.path(out_dir, "sectionD-peptideranger.csv")
if (requireNamespace("PeptideRanger", quietly = TRUE)) {
  library(PeptideRanger)
  rf_model <- RFmodel_ProteomicsDB

  pr_rows <- list()
  for (pn in names(proteins)) {
    path <- proteins[[pn]]$path
    for (enz in enzymes) {
      ev <- evaluate_digest(path, enzyme = enz, missed_cleavages = 1L)
      peps <- ev$peptides
      names(peps)[2] <- "peptide"
      peps$gravy <- .calculate_gravy(peps$peptide)
      peps$window_pass <- peps$length >= 7 & peps$length <= 25 &
        peps$gravy >= -1.0 & peps$gravy <= 0.6

      pr_res <- peptide_predictions(peps$peptide,
        prediction_model = rf_model
      )
      peps$pr_score <- pr_res$RF_score[match(
        peps$peptide,
        pr_res$sequence
      )]

      pr_rows[[length(pr_rows) + 1]] <- data.frame(
        protein = pn,
        enzyme = enz,
        verdict = ev$scores$verdict,
        n_total = nrow(peps),
        n_window_pass = sum(peps$window_pass),
        PR_mean_all = mean(peps$pr_score, na.rm = TRUE),
        PR_mean_window_pass = mean(
          peps$pr_score[peps$window_pass], na.rm = TRUE
        ),
        PR_mean_window_fail = mean(
          peps$pr_score[!peps$window_pass], na.rm = TRUE
        ),
        mean_score_difference = mean(
          peps$pr_score[peps$window_pass], na.rm = TRUE
        ) - mean(peps$pr_score[!peps$window_pass], na.rm = TRUE),
        cor_PR_length = cor(peps$pr_score, peps$length,
          use = "complete.obs"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  pr_tbl <- .bind_rows(pr_rows)
  print(pr_tbl, row.names = FALSE)
  write.csv(pr_tbl, peptideranger_output,
    row.names = FALSE
  )

  cat("\nDescriptive comparison:\n")
  cat("  pepVet supplies digest-level model scores.\n")
  cat("  PeptideRanger supplies peptide-level random-forest scores.\n")
  cat("  The reported difference compares selected-window group means.\n")
} else {
  stop(
    "PeptideRanger is not installed; Section D cannot be regenerated.",
    call. = FALSE
  )
}

## ======================================================================
## SECTION E: PROTEIN CLEAVER ANALYSIS
## Based on source code analysis (gkoulouras/ProteinCleaver).
## ======================================================================
cat("\nSection E: Protein Cleaver source analysis\n\n")
cat("Source: https://github.com/gkoulouras/ProteinCleaver (revision not recorded)\n")
cat("Engine: cleaver::cleave() + cleaver::cleavageRanges()\n")
cat("Identifiable = length in window AND mass in window (binary flag)\n")
cat("No GRAVY, no graded score, no verdict, no presets\n")
cat("Documented outputs include sequence and 3D structure mapping.\n")
cat("The assessed revision did not include a GRAVY-based digest score.\n\n")

pc_mass_table <- c(
  "A" = 71.0788, "C" = 103.1388, "D" = 115.0886, "E" = 129.1155,
  "F" = 147.1766, "G" = 57.0519, "H" = 137.1411, "I" = 113.1594,
  "K" = 128.1741, "L" = 113.1594, "M" = 131.1926, "N" = 114.1038,
  "O" = 237.3018, "P" = 97.1167, "Q" = 128.1307, "R" = 156.1875,
  "S" = 87.0782, "T" = 101.1051, "U" = 150.0388, "V" = 99.1326,
  "W" = 186.2132, "Y" = 163.1760, "H2O" = 18.01524
)

simulate_pc_identifiable <- function(pep_seq, min_len = 7, max_len = 35,
                                     min_mass = 400, max_mass = 4000) {
  if (!is.character(pep_seq) || length(pep_seq) != 1L ||
      is.na(pep_seq) || !grepl("^[ACDEFGHIKLMNPQRSTVWYOU]+$", pep_seq)) {
    return(FALSE)
  }
  len <- nchar(pep_seq)
  if (len < min_len || len > max_len) {
    return(FALSE)
  }
  aa <- strsplit(pep_seq, "")[[1]]
  mass <- sum(pc_mass_table[aa]) + pc_mass_table["H2O"]
  mass >= min_mass && mass <= max_mass
}

pc_comparison <- list()
for (pn in c("BSA", "H3", "BACE1")) {
  ev <- pepvet_results[[paste(pn, "trypsin", sep = "_")]]
  peps <- ev$peptides
  names(peps)[2] <- "peptide"
  peps$gravy <- .calculate_gravy(peps$peptide)
  peps$pc_identifiable <- vapply(peps$peptide, simulate_pc_identifiable,
                                  logical(1))
  peps$pepvet_window_pass <- peps$length >= 7 & peps$length <= 25 &
    peps$gravy >= -1.0 & peps$gravy <= 0.6

  pc_comparison[[pn]] <- data.frame(
    protein = pn,
    n_peptides = nrow(peps),
    pc_identifiable = sum(peps$pc_identifiable),
    pepvet_window_pass = sum(peps$pepvet_window_pass),
    both_pass = sum(peps$pc_identifiable & peps$pepvet_window_pass),
    pc_only = sum(peps$pc_identifiable & !peps$pepvet_window_pass),
    pepvet_window_only = sum(
      peps$pepvet_window_pass & !peps$pc_identifiable
    ),
    stringsAsFactors = FALSE
  )
}
pc_tbl <- .bind_rows(pc_comparison)
cat("Protein Cleaver simulation vs pepVet (trypsin MC=1):\n")
print(pc_tbl, row.names = FALSE)
write.csv(pc_tbl, file.path(out_dir, "sectionE-protein-cleaver.csv"),
  row.names = FALSE
)

## ======================================================================
## GENERATION SUMMARY
## ======================================================================
cat("\nGenerated descriptive tables for the documentation vignette.\n")

##
## FIGURES
##
cat("\nGenerating figures\n\n")

enzyme_labels <- c(
  "trypsin" = "Trypsin", "lysc" = "Lys-C",
  "chymotrypsin-high" = "Chymotrypsin",
  "glutamyl endopeptidase" = "Glu-C",
  "asp-n endopeptidase" = "Asp-N"
)
verdict_fill <- c(
  "Good" = "#3A8C5F", "Moderate" = "#D4A76A",
  "Poor" = "#C46A6A"
)
enzyme_fill <- c("#3B7A9E", "#6BA292", "#D4A76A", "#C46A6A", "#8B7EB5")
names(enzyme_fill) <- names(enzyme_labels)

.label_color <- function(v) ifelse(v > 0.5, "white", "grey20")

## Figure 1: Faceted scoring heatmap
score_cols <- c("S_length", "S_coverage", "S_count", "S_hydro", "S_charge")
scores_long <- tidyr::pivot_longer(
  scores_tbl[, c("protein", "enzyme", score_cols)],
  cols = all_of(score_cols),
  names_to = "component",
  values_to = "value"
)
scores_long$enzyme <- factor(scores_long$enzyme, levels = enzymes)
scores_long$component <- factor(scores_long$component, levels = score_cols)

p1 <- ggplot(
  scores_long,
  aes(x = enzyme, y = component, fill = value)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(
    label = sprintf("%.2f", value),
    color = value <= 0.5
  ), size = 2.8) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "grey20"),
    guide = "none"
  ) +
  scale_fill_viridis_c(option = "D", limits = c(0, 1), name = "Score") +
  scale_x_discrete() +
  facet_wrap(~protein, nrow = 1) +
  labs(
    title = "pepVet scoring dimensions",
    subtitle = "Component scores across 5 proteins and 5 enzymes, MC=1",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig1_scoring_heatmap.png"), p1,
  width = 12, height = 4, dpi = 150
)
cat("  fig1_scoring_heatmap.png\n")

## Figure 2: Enzyme comparison with verdict threshold lines
p2 <- ggplot(scores_tbl, aes(x = protein, y = composite, fill = enzyme)) +
  geom_col(position = "dodge", color = "black", alpha = 0.85) +
  geom_hline(yintercept = c(0.40, 0.65),
    linetype = "dashed", linewidth = 0.5,
    color = c("#D4A76A", "#3A8C5F")
  ) +
  scale_x_discrete(expand = expansion(add = c(0.6, 1.0))) +
  annotate("label", x = 5.6, y = 0.65,
    label = "Good", hjust = 0.5, vjust = -0.4, size = 3.2,
    fill = "#3A8C5F", color = "white", fontface = "bold"
  ) +
  annotate("label", x = 5.6, y = 0.40,
    label = "Moderate", hjust = 0.5, vjust = 1.4, size = 3.2,
    fill = "#D4A76A", color = "white", fontface = "bold"
  ) +
  annotate("label", x = 5.6, y = 0.15,
    label = "Poor", hjust = 0.5, vjust = 0.5, size = 3.2,
    fill = "#C46A6A", color = "white", fontface = "bold"
  ) +
  scale_fill_manual(
    values = enzyme_fill,
    labels = enzyme_labels[names(enzyme_fill)]
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.02)) +
  labs(
    title = "Enzyme comparison across proteins",
    y = "Composite score", x = NULL, fill = "Enzyme"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fig2_enzyme_comparison.png"), p2,
  width = 8, height = 5, dpi = 150
)
cat("  fig2_enzyme_comparison.png\n")

## Figure 3: PR mean-score difference faceted by protein
if (exists("pr_tbl")) {
  pr_enrich <- pr_tbl[!is.na(pr_tbl$mean_score_difference), ]
  pr_enrich$sign <- ifelse(
    pr_enrich$mean_score_difference >= 0, "positive", "negative"
  )
  pr_enrich$enzyme_label <- enzyme_labels[pr_enrich$enzyme]

  mean_enr <- tapply(
    pr_enrich$mean_score_difference,
    pr_enrich$enzyme_label,
    mean,
    na.rm = TRUE
  )
  enz_order <- names(sort(mean_enr))
  pr_enrich$enzyme_label <- factor(pr_enrich$enzyme_label, levels = enz_order)

  p3 <- ggplot(
    pr_enrich,
    aes(x = mean_score_difference, y = enzyme_label, fill = sign)
  ) +
    geom_col(color = "black", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = sprintf("%+.2f", mean_score_difference)),
      hjust = ifelse(pr_enrich$mean_score_difference >= 0, -0.2, 1.2),
      size = 2.5
    ) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey50") +
    scale_fill_manual(values = c(
      "positive" = "#3A8C5F",
      "negative" = "#C46A6A"
    )) +
    facet_grid(protein ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "PeptideRanger score difference by selected window",
      subtitle = paste0(
        "Difference in mean PR score between peptides inside and outside ",
        "the selected length-and-GRAVY window"
      ),
      x = "Mean score difference", y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 1, face = "bold", size = 9),
      strip.background = element_rect(fill = "#F0F0F0", color = NA)
    )

  ggsave(file.path(fig_dir, "fig3_pr_enrichment.png"), p3,
    width = 9, height = 6, dpi = 150
  )
  cat("  fig3_pr_enrichment.png\n")
}

## Figure 4: Overlap ordered by count
overlap_plot <- data.frame(
  category = c(
    "All three tools",
    "pepVet + MS-Digest", "pepVet + ExPASy",
    "MS-Digest + ExPASy",
    "pepVet only", "MS-Digest only", "ExPASy only"
  ),
  count = c(
    sum(pvp_peps %in% msd_peps & pvp_peps %in% exp_peps),
    sum(pvp_peps %in% msd_peps & !(pvp_peps %in% exp_peps)),
    sum(pvp_peps %in% exp_peps & !(pvp_peps %in% msd_peps)),
    sum(msd_peps %in% exp_peps & !(msd_peps %in% pvp_peps)),
    sum(!(pvp_peps %in% msd_peps) & !(pvp_peps %in% exp_peps)),
    sum(!(msd_peps %in% pvp_peps) & !(msd_peps %in% exp_peps)),
    sum(!(exp_peps %in% pvp_peps) & !(exp_peps %in% msd_peps))
  )
)
overlap_plot$category <- reorder(overlap_plot$category, overlap_plot$count)

p4 <- ggplot(overlap_plot, aes(x = category, y = count)) +
  geom_col(fill = "#3B7A9E", color = "black", alpha = 0.85) +
  geom_text(aes(label = count), hjust = -0.3, size = 3.5) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, max(overlap_plot$count) * 1.15),
    expand = c(0, 0)
  ) +
  labs(
    title = "Peptide overlap: BSA trypsin MC=1",
    subtitle = paste0(
      "pepVet returns the most peptides (", length(pvp_peps), ") ",
      "as it applies no default mass or length filters"
    ),
    x = NULL, y = "Number of peptides"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(fig_dir, "fig4_overlap.png"), p4,
  width = 8, height = 5, dpi = 150
)
cat("  fig4_overlap.png\n")

## Figure 5: Preset effects faceted, with criteria annotations
presets_detail <- data.frame(
  preset = names(.pepvet_presets),
  detail = vapply(names(.pepvet_presets), function(n) {
    p <- .pepvet_presets[[n]]
    sprintf(
      "len %d-%d\nGRAVY %.1f-%.1f",
      p$length_range[1], p$length_range[2],
      p$gravy_range[1], p$gravy_range[2]
    )
  }, character(1)),
  stringsAsFactors = FALSE
)

p5 <- ggplot(preset_tbl, aes(x = preset)) +
  geom_col(aes(y = n_total), fill = "#D8DDE6", alpha = 0.7, color = NA, width = 0.85) +
  geom_col(aes(y = n_length_valid, fill = verdict), color = "black", alpha = 0.85, width = 0.85) +
  geom_text(aes(y = n_length_valid, label = n_length_valid), vjust = -0.3, size = 3.2) +
  geom_text(aes(y = n_length_valid, label = sprintf("%.0f%%", 100 * n_length_valid / n_total)),
    vjust = 1.8, size = 2.5, color = "grey40"
  ) +
  geom_text(
    data = presets_detail, aes(x = preset, y = 0, label = detail),
    vjust = 1.3, size = 2.5, color = "grey40", lineheight = 0.9,
    inherit.aes = FALSE
  ) +
  facet_wrap(~protein, nrow = 1) +
  scale_fill_manual(values = verdict_fill) +
  labs(
    title = "Preset effects on length-valid peptide count",
    subtitle = paste0(
      "Grey bar = total digested peptides. ",
      "Coloured bar = peptides within the preset length range. ",
      "Percent = length-valid / total."
    ),
    x = NULL, y = "Number of peptides", fill = "Verdict"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    legend.position = "bottom",
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "fig5_preset_effect.png"), p5,
  width = 10, height = 5, dpi = 150
)
cat("  fig5_preset_effect.png\n")

## Figure 6: Protein Cleaver vs pepVet with annotations
pc_long <- tidyr::pivot_longer(
  pc_tbl[, c("protein", "both_pass", "pc_only")],
  cols = c("both_pass", "pc_only"),
  names_to = "classifier",
  values_to = "count"
)
pc_long$classifier <- ifelse(pc_long$classifier == "both_pass",
  "Both filters", "PC identifiable only"
)

p6 <- ggplot(pc_long, aes(x = protein, y = count, fill = classifier)) +
  geom_col(position = "dodge", color = "black", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = count),
    position = position_dodge(width = 0.6),
    vjust = -0.3, size = 3.5
  ) +
  scale_fill_manual(values = c(
    "Both filters" = "#3A8C5F",
    "PC identifiable only" = "#C46A6A"
  )) +
  labs(
    title = "pepVet vs Protein Cleaver (PC)",
    subtitle = paste0(
      "pepVet applies a GRAVY filter (KD -1.0 to 0.6) ",
      "that PC lacks. PC uses a different mass table."
    ),
    x = NULL, y = "Number of peptides", fill = "Classification"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "fig6_pc_comparison.png"), p6,
  width = 8, height = 5, dpi = 150
)
cat("  fig6_pc_comparison.png\n")

cat("\nCOMPARISON COMPLETE\n")
cat(sprintf("Data: %s/*.csv\n", out_dir))
cat(sprintf("Figures: %s/*.png\n", fig_dir))
cat(sprintf("Total figures: %d\n", length(list.files(fig_dir, pattern = "\\.png$"))))
