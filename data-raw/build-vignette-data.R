## data-raw/build-vignette-data.R
## Compile all comparison data into vignette-ready RDS datasets.
## Run after fetch-external-data.R, full-tool-comparison.R

pkgload::load_all(quiet = TRUE)

tool_dir <- file.path("inst", "extdata", "tool-data")
comp_dir <- file.path("inst", "extdata", "comparison-data")
out_dir  <- "vignettes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##
## 1. Load raw CSV data
##
pepvet_raw   <- read.csv(file.path(tool_dir, "pepvet-bsa-trypsin.csv"),
                          stringsAsFactors = FALSE)
msdigest_raw <- read.csv(file.path(tool_dir, "tool-compare-msdigest.csv"),
                          stringsAsFactors = FALSE)
expasy_raw   <- read.csv(file.path(tool_dir, "tool-compare-expasy.csv"),
                          stringsAsFactors = FALSE)
presets      <- read.csv(file.path(comp_dir, "sectionC2-presets.csv"),
                          stringsAsFactors = FALSE)
pr_summary   <- read.csv(file.path(tool_dir, "tool-compare-peptideranger.csv"),
                          stringsAsFactors = FALSE)
pr_peptides  <- read.csv(file.path(tool_dir, "peptideranger-all-peptides.csv"),
                          stringsAsFactors = FALSE)

cat(sprintf("pepVet:   %d rows\n", nrow(pepvet_raw)))
cat(sprintf("MS-Digest: %d rows\n", nrow(msdigest_raw)))
cat(sprintf("ExPASy:   %d rows\n", nrow(expasy_raw)))

##
## 2. Compute pepVet-derived properties (mass, mz, gravy)
##
peptides_seqs <- pepvet_raw$peptide

pepvet_raw$gravy <- .calculate_gravy_vec(peptides_seqs)

mass_mi <- calculate_peptide_mass(peptides_seqs, charge = 0L)
mass_mh <- calculate_peptide_mass(peptides_seqs, charge = 1L)

pepvet_raw$mass_neutral_mi <- mass_mi
pepvet_raw$mz_mh_plus      <- mass_mh

##
## 3. Build tool-level summary table
##
all_pepvet <- unique(pepvet_raw$peptide)
all_msd    <- unique(msdigest_raw$peptide)
all_exp    <- unique(expasy_raw$peptide)

tool_summary <- data.frame(
  Tool = c("pepVet", "MS-Digest", "ExPASy PeptideMass"),
  Total_peptides = c(length(all_pepvet), length(all_msd), length(all_exp)),
  Common_with_pepVet = c(
    length(all_pepvet),
    sum(all_msd %in% all_pepvet),
    sum(all_exp %in% all_pepvet)
  ),
  stringsAsFactors = FALSE
)
tool_summary$Pct_common <- round(
  tool_summary$Common_with_pepVet / tool_summary$Total_peptides * 100, 1
)

cat("\nTool summary:\n")
print(tool_summary, row.names = FALSE)

##
## 4. Build master peptide comparison table
##
master <- pepvet_raw[, c("peptide", "start", "end", "length",
                          "missed_cleavages", "gravy",
                          "mass_neutral_mi", "mz_mh_plus")]

master$in_pepvet   <- TRUE
master$in_msdigest <- master$peptide %in% all_msd
master$in_expasy   <- master$peptide %in% all_exp
master$n_tools     <- 1L + master$in_msdigest + master$in_expasy

idx_msd <- match(master$peptide, msdigest_raw$peptide)
master$msdigest_mz_mi <- msdigest_raw$mz_mi[idx_msd]
master$msdigest_mz_av <- msdigest_raw$mz_av[idx_msd]

idx_exp <- match(master$peptide, expasy_raw$peptide)
master$expasy_mass <- expasy_raw$mass[idx_exp]

idx_pr <- match(master$peptide,
                pr_peptides$peptide[pr_peptides$protein == "BSA" &
                                    pr_peptides$enzyme == "trypsin"])
master$pr_score <- pr_peptides$pr_score[idx_pr]

master <- master[order(master$start), ]

saveRDS(master, file.path(out_dir, "vignette-peptides.rds"))
cat(sprintf("\nvignette-peptides.rds: %d rows, %d cols\n",
            nrow(master), ncol(master)))

##
## 5. Preset comparison
##
saveRDS(presets, file.path(out_dir, "vignette-presets.rds"))
cat(sprintf("vignette-presets.rds: %d rows\n", nrow(presets)))

##
## 6. PeptideRanger summary
##
saveRDS(pr_summary, file.path(out_dir, "vignette-peptideranger.rds"))
cat(sprintf("vignette-peptideranger.rds: %d rows\n", nrow(pr_summary)))

##
## 7. Mass correlation summary
##
common_pep <- master[master$in_msdigest & master$in_expasy, ]
cor_msd_exp <- cor(common_pep$msdigest_mz_mi, common_pep$expasy_mass,
                   use = "complete.obs")
cor_pep_msd <- cor(master$mz_mh_plus[master$in_msdigest],
                   master$msdigest_mz_mi[master$in_msdigest],
                   use = "complete.obs")
cor_pep_exp <- cor(master$mz_mh_plus[master$in_expasy],
                   master$expasy_mass[master$in_expasy],
                   use = "complete.obs")

mass_corr <- data.frame(
  tool_pair = c("pepVet vs MS-Digest", "pepVet vs ExPASy",
                "MS-Digest vs ExPASy"),
  correlation = round(c(cor_pep_msd, cor_pep_exp, cor_msd_exp), 6),
  n_common = c(sum(master$in_msdigest), sum(master$in_expasy),
               sum(master$in_msdigest & master$in_expasy)),
  stringsAsFactors = FALSE
)
cat("\nMass correlations:\n")
print(mass_corr, row.names = FALSE)

saveRDS(mass_corr, file.path(out_dir, "vignette-mass-corr.rds"))

##
## 8. Overlap summary
##
overlap_summary <- data.frame(
  category = c("All three tools", "pepVet only", "MS-Digest only",
               "ExPASy only", "pepVet + MS-Digest", "pepVet + ExPASy",
               "MS-Digest + ExPASy"),
  n_peptides = c(
    sum(master$n_tools == 3),
    sum(master$n_tools == 1 & master$in_pepvet),
    sum(!(all_msd %in% master$peptide[master$in_msdigest])),
    sum(!(all_exp %in% master$peptide[master$in_expasy])),
    sum(master$in_msdigest & !master$in_expasy),
    sum(master$in_expasy & !master$in_msdigest),
    sum(master$in_msdigest & master$in_expasy)
  ),
  stringsAsFactors = FALSE
)
cat("\nOverlap summary:\n")
print(overlap_summary, row.names = FALSE)

saveRDS(overlap_summary, file.path(out_dir, "vignette-overlap.rds"))

cat("\nAll vignette data written to inst/extdata/\n")
