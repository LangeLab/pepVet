# Shared test fixtures.
#
# testthat sources this file before any test file runs. Expensive
# evaluate_digest() and batch_evaluate() calls belong here so they run once per
# test-suite run rather than once per test block.
#
# Fixture names use .fix_<protein>_<enzyme>[_mc<n>]. A leading dot keeps them
# out of ordinary test output and reduces name collisions.

# Path helper shared across all test files.
reference_fasta <- function(file_name) {
  system.file("extdata", file_name, package = "pepVet")
}

.bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
.h3_path <- system.file("extdata", "P68431.fasta", package = "pepVet")
.small_path <- system.file("extdata", "small_proteome_50_proteins.fasta",
  package = "pepVet"
)
.multi_path <- system.file("extdata", "P37840_isoforms.fasta",
  package = "pepVet"
)

# Single-protein evaluate_digest results.
.fix_bsa_trypsin <- evaluate_digest(.bsa_path, enzyme = "trypsin")
.fix_h3_trypsin <- evaluate_digest(.h3_path, enzyme = "trypsin")
.fix_bsa_lysc <- evaluate_digest(.bsa_path, enzyme = "lysc")
.fix_h3_lysc <- evaluate_digest(.h3_path, enzyme = "lysc")
# Chymotrypsin-high cuts after F/Y/W/L and creates a different size profile.
.fix_bsa_chymotryp <- evaluate_digest(.bsa_path, enzyme = "chymotrypsin-high")
.fix_h3_chymotryp <- evaluate_digest(.h3_path, enzyme = "chymotrypsin-high")
# Asp-N cuts N-terminal to D and exercises the opposite direction.
.fix_bsa_aspn <- evaluate_digest(.bsa_path, enzyme = "asp-n endopeptidase")
# Glutamyl endopeptidase cuts after D/E and differs from trypsin.
.fix_bsa_glute <- evaluate_digest(.bsa_path, enzyme = "glutamyl endopeptidase")
.fix_bsa_mc0 <- evaluate_digest(.bsa_path,
  enzyme = "trypsin",
  missed_cleavages = 0L
)
.fix_bsa_mc1 <- evaluate_digest(.bsa_path,
  enzyme = "trypsin",
  missed_cleavages = 1L
)
.fix_bsa_mc2 <- evaluate_digest(.bsa_path,
  enzyme = "trypsin",
  missed_cleavages = 2L
)
.fix_bsa_lysc_mc1 <- evaluate_digest(.bsa_path,
  enzyme = "lysc",
  missed_cleavages = 1L
)
.fix_bsa_chymotryp_mc1 <- evaluate_digest(.bsa_path,
  enzyme = "chymotrypsin-high",
  missed_cleavages = 1L
)

# Batch results.
.fix_batch_trypsin <- batch_evaluate(
  c(
    Biostrings::readAAStringSet(.bsa_path),
    Biostrings::readAAStringSet(.h3_path)
  ),
  enzyme = "trypsin"
)
.fix_batch_small <- batch_evaluate(.small_path, enzyme = "trypsin")
.fix_batch_bsa <- batch_evaluate(.bsa_path, enzyme = "trypsin")
.fix_batch_bsa_mc1 <- batch_evaluate(.bsa_path,
  enzyme = "trypsin",
  missed_cleavages = 1L
)
.fix_batch_h3 <- batch_evaluate(.h3_path, enzyme = "trypsin")
.fix_batch_multi <- batch_evaluate(.multi_path)
.fix_batch_chymotryp <- batch_evaluate(
  c(
    Biostrings::readAAStringSet(.bsa_path),
    Biostrings::readAAStringSet(.h3_path)
  ),
  enzyme = "chymotrypsin-high"
)

# Plotting wrappers return cached fixtures for common missed-cleavage levels.
.bsa_result <- function(mc = 0L) {
  if (mc == 0L) {
    return(.fix_bsa_mc0)
  }
  if (mc == 1L) {
    return(.fix_bsa_mc1)
  }
  if (mc == 2L) {
    return(.fix_bsa_mc2)
  }
  evaluate_digest(.bsa_path, enzyme = "trypsin", missed_cleavages = mc)
}

.h3_result <- function(mc = 0L) {
  if (mc == 0L) {
    return(.fix_h3_trypsin)
  }
  evaluate_digest(.h3_path, enzyme = "trypsin", missed_cleavages = mc)
}

.bsa_cs <- function() {
  annotate_cleavage_sites(.bsa_path, enzyme = "trypsin")
}

.h3_cs <- function() {
  annotate_cleavage_sites(.h3_path, enzyme = "trypsin")
}

.bsa_comparison <- function(
  enzymes = c(
    "trypsin", "lysc",
    "glutamyl endopeptidase"
  )
) {
  compare_digests(.bsa_path, enzymes = enzymes)
}
