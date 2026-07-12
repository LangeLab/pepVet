# Generate example plots for README
# Run from repo root: Rscript man/figures/generate-readme-plots.R

library(ggplot2)
library(patchwork)

if (!file.exists("DESCRIPTION") ||
    !identical(read.dcf("DESCRIPTION", fields = "Package")[[1L]], "pepVet")) {
  stop("Run this script from the pepVet repository root.", call. = FALSE)
}

devtools::load_all(".", quiet = TRUE)

out_dir <- Sys.getenv("PEPVET_FIGURE_DIR", unset = "man/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
small_path <- system.file("extdata", "small_proteome_50_proteins.fasta", package = "pepVet")

cat("Generating digest_profile_bsa_trypsin.png ...\n")
res <- evaluate_digest(bsa_path, enzyme = "trypsin", missed_cleavages = 1L)
p1 <- plot_digest_profile(res)
ggplot2::ggsave(file.path(out_dir, "digest_profile_bsa_trypsin.png"),
  p1, device = ragg::agg_png, width = 12, height = 10, dpi = 150,
  bg = "white")

cat("Generating batch_comparison_10_enzymes_50_proteins.png ...\n")
enzymes <- c("trypsin", "lysc", "chymotrypsin-high", "asp-n endopeptidase",
  "glutamyl endopeptidase", "arg-c proteinase", "thermolysin",
  "pepsin", "staphylococcal peptidase i", "proteinase k")
batch_comp <- batch_compare_enzymes(small_path, enzymes = enzymes)
p2 <- plot_batch_comparison(batch_comp)
ggplot2::ggsave(file.path(out_dir, "batch_comparison_10_enzymes_50_proteins.png"),
  p2, device = ragg::agg_png, width = 18, height = 12, dpi = 150,
  bg = "white")

cat("Done.\n")
