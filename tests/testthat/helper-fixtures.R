# ── Shared test fixtures ──────────────────────────────────────────────────────
#
# This file is sourced automatically by testthat before any test file runs.
# All expensive evaluate_digest() and batch_evaluate() calls go here so they
# execute ONCE per test-suite run, not once per test_that() block.
#
# Naming convention:
#   .fix_<protein>_<enzyme>[_mc<n>]
#
# Objects are prefixed with "." so they are invisible in the workspace and
# unlikely to clash with test-local variables.
# ─────────────────────────────────────────────────────────────────────────────

# Path helpers — shared across all test files (previously local to each file)
reference_fasta <- function(file_name) {
  system.file("extdata", file_name, package = "pepVet")
}

.bsa_path    <- system.file("extdata", "P02769.fasta", package = "pepVet")
.h3_path     <- system.file("extdata", "P68431.fasta", package = "pepVet")
.small_path  <- system.file("extdata", "small_proteome_50_proteins.fasta",
                             package = "pepVet")
.multi_path  <- system.file("extdata", "P37840_isoforms.fasta",
                             package = "pepVet")

# Single-protein evaluate_digest results ──────────────────────────────────────
.fix_bsa_trypsin    <- evaluate_digest(.bsa_path, enzyme = "trypsin")
.fix_h3_trypsin     <- evaluate_digest(.h3_path,  enzyme = "trypsin")
.fix_bsa_lysc       <- evaluate_digest(.bsa_path, enzyme = "lysc")
.fix_h3_lysc        <- evaluate_digest(.h3_path,  enzyme = "lysc")
# chymotrypsin-high cuts after F/Y/W/L — very different peptide size profile
.fix_bsa_chymotryp  <- evaluate_digest(.bsa_path, enzyme = "chymotrypsin-high")
.fix_h3_chymotryp   <- evaluate_digest(.h3_path,  enzyme = "chymotrypsin-high")
# asp-n endopeptidase cuts N-terminal to D — unique directionality
.fix_bsa_aspn       <- evaluate_digest(.bsa_path, enzyme = "asp-n endopeptidase")
# glutamyl endopeptidase cuts after D/E — moderate count, different than trypsin
.fix_bsa_glute      <- evaluate_digest(.bsa_path, enzyme = "glutamyl endopeptidase")
.fix_bsa_mc0        <- evaluate_digest(.bsa_path, enzyme = "trypsin",
                                        missed_cleavages = 0L)
.fix_bsa_mc1        <- evaluate_digest(.bsa_path, enzyme = "trypsin",
                                        missed_cleavages = 1L)
.fix_bsa_mc2        <- evaluate_digest(.bsa_path, enzyme = "trypsin",
                                        missed_cleavages = 2L)
.fix_bsa_lysc_mc1   <- evaluate_digest(.bsa_path, enzyme = "lysc",
                                        missed_cleavages = 1L)
.fix_bsa_chymotryp_mc1 <- evaluate_digest(.bsa_path, enzyme = "chymotrypsin-high",
                                            missed_cleavages = 1L)

# Batch result ────────────────────────────────────────────────────────────────
if (requireNamespace("Biostrings", quietly = TRUE)) {
  .fix_batch_trypsin <- batch_evaluate(
    c(Biostrings::readAAStringSet(.bsa_path),
      Biostrings::readAAStringSet(.h3_path)),
    enzyme = "trypsin"
  )
  # Batch for evaluation tests
  .fix_batch_small        <- batch_evaluate(.small_path, enzyme = "trypsin")
  .fix_batch_bsa          <- batch_evaluate(.bsa_path,   enzyme = "trypsin")
  .fix_batch_bsa_mc1      <- batch_evaluate(.bsa_path,   enzyme = "trypsin",
                                             missed_cleavages = 1L)
  .fix_batch_h3           <- batch_evaluate(.h3_path,    enzyme = "trypsin")
  .fix_batch_multi        <- batch_evaluate(.multi_path)
  # Non-trypsin batch fixtures for enzyme-diversity testing
  .fix_batch_chymotryp    <- batch_evaluate(
    c(Biostrings::readAAStringSet(.bsa_path),
      Biostrings::readAAStringSet(.h3_path)),
    enzyme = "chymotrypsin-high"
  )
} else {
  .fix_batch_trypsin   <- NULL
  .fix_batch_small     <- NULL
  .fix_batch_bsa       <- NULL
  .fix_batch_bsa_mc1   <- NULL
  .fix_batch_h3        <- NULL
  .fix_batch_multi     <- NULL
  .fix_batch_chymotryp <- NULL
}
