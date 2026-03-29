# pepVet

pepVet is an R package for proteolytic digest simulation, peptide-set scoring, enzyme comparison, and workflow preset selection for bottom-up proteomics. Give it a protein sequence and an enzyme. It returns peptide coordinates, score components, a ranked comparison, and enough detail to explain why one enzyme is a better starting point than another.

## What pepVet does

The choice of proteolytic enzyme shapes every downstream result in a proteomics experiment. Cut too aggressively and you get thousands of tiny fragments below the detection threshold. Cut too conservatively and overlong peptides fail to fly or resolve on the column.

pepVet quantifies this trade-off with five orthogonal scoring components and a weighted composite score.

## Version 0.1.0

Version 0.1.0 completes the current Phase 2 model pass and turns pepVet into a more usable planning tool.

- `S_count` uses an enzyme-aware expected peptide length instead of a fixed trypsin-derived denominator.
- Workflow presets configure peptide-length windows, GRAVY windows, weights, and preset tracking for common experiment types.
- `calculate_peptide_mass()` and `calculate_pI()` provide peptide-level mass and pI utilities for fractionation-aware planning.
- `annotate_cleavage_sites()` and optional digest-level `cleavage_efficiency` flags expose sequence-local missed-cleavage risk for trypsin-family digests.

## Installation

pepVet depends on Bioconductor packages. `cleaver` is required and should be installed through `BiocManager` before the GitHub install.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("Biostrings", "IRanges", "S4Vectors", "cleaver"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("LangeLab/pepVet", dependencies = TRUE)
```

## Core workflow

```r
library(pepVet)

bsa <- system.file("extdata", "P02769.fasta", package = "pepVet")

ev <- evaluate_digest(bsa, enzyme = "trypsin", missed_cleavages = 1L)
ev$scores

comp <- compare_digests(
  bsa,
  enzymes = c("trypsin", "lysc", "glutamyl endopeptidase", "asp-n endopeptidase")
)
comp

digest_report(comp)
recommend_enzyme(bsa, enzymes = c("trypsin", "lysc", "glutamyl endopeptidase"))
```

`digest_protein()` returns a tibble with one row per peptide. `score_peptides()` returns a tibble with one row per protein. `compare_digests()` returns a tibble with one row per enzyme. `evaluate_digest()` returns a named list with `scores`, `peptides`, and `params`.

## Cleavage efficiency annotations

pepVet now distinguishes between cleavage-site location and cleavage-site confidence for trypsin-family digests. `annotate_cleavage_sites()` classifies each candidate K/R site as `high`, `medium`, or `low` efficiency using local P1-P1' sequence context.

```r
annotate_cleavage_sites(bsa, enzyme = "trypsin")

ev_eff <- evaluate_digest(
  bsa,
  enzyme = "trypsin",
  missed_cleavages = 1L,
  include_cleavage_efficiency = TRUE
)

ev_eff$peptides
ev_eff$scores[, c("protein_id", "n_high_efficiency_sites", "n_low_efficiency_sites")]
```

These annotations are informational. They are not score components, and they are intentionally limited to local sequence context. pepVet does not yet model structural accessibility, PTMs, or extended subsite preferences.

## Workflow presets

Each preset returns a named list with `gravy_range`, `length_range`, and `weights`, so you can pass it straight into `evaluate_digest()` or `score_peptides()`.

```r
pepvet_preset("standard")

syn_path <- system.file("extdata", "P37840_isoforms.fasta", package = "pepVet")
syn_proteome <- digest_protein(syn_path, enzyme = "trypsin")
targeted <- pepvet_preset("targeted")
do.call(
  evaluate_digest,
  c(list(sequence = syn_path, enzyme = "trypsin", proteome = syn_proteome), targeted)
)
```

| Preset          | Best fit                   | Main scoring shift                                        |
| --------------- | -------------------------- | --------------------------------------------------------- |
| `standard`      | Routine DDA                | Package defaults                                          |
| `dia`           | DIA and SWATH              | Higher coverage weight                                    |
| `targeted`      | SRM, PRM, MRM              | Higher uniqueness weight                                  |
| `membrane`      | Hydrophobic proteins       | Wider GRAVY window                                        |
| `ffpe_degraded` | Degraded samples           | Wider peptide-length window                               |
| `fractionated`  | High-pH RP or SCX planning | Standard score with fractionation-oriented interpretation |

Presets with non-zero `S_unique` weights require a proteome digest. pepVet rejects a non-zero uniqueness weight when no proteome is supplied.

## Scoring model

pepVet scores each digest with five core components and one optional proteome-aware component.

| Score        | What it measures                                                        | Why it matters                                                                                   |
| ------------ | ----------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `S_length`   | Fraction of peptides inside the active peptide-length window            | Very short and very long peptides lower identification rates                                     |
| `S_coverage` | Fraction of the parent protein covered by valid peptides                | Dark regions weaken protein-level interpretation                                                 |
| `S_count`    | Valid peptide count relative to the enzyme-aware expected density       | Too few peptides weaken evidence. Too many reflect over-digestion                                |
| `S_hydro`    | Fraction of valid peptides inside the active GRAVY window               | Extreme hydrophobicity or hydrophilicity hurts LC behavior                                       |
| `S_charge`   | Fraction of valid peptides with at least one non-terminal basic residue | Higher values indicate more opportunities for multi-charge states and richer fragment ion series |
| `S_unique`   | Fraction of valid peptides unique in a supplied proteome                | Shared peptides weaken protein-level attribution                                                 |

`S_charge` does not mean a peptide can or cannot ionize. Tryptic peptides still carry the free N-terminus and often a terminal Lys or Arg. The score is meant to distinguish baseline ionizability from extra internal basic-residue richness, which tends to support higher charge states and richer b/y ion series.

The composite score is a weighted sum. Verdict thresholds remain `Good` at `>= 0.70`, `Moderate` at `>= 0.40`, and `Poor` below `0.40`. These thresholds are heuristic. They are not calibrated probabilities. [REF]

## Positioning and scope

pepVet is not a peptide detectability predictor. It is a rule-based, multi-criteria digest-ranking model for pre-acquisition planning.

- ML detectability tools estimate per-peptide detection from experimental training data.
- pepVet ranks per-protein digest quality from an in silico digest and explicit workflow assumptions.
- pepVet scores are interpretable rankings, not calibrated probabilities.

Cross-workflow comparisons are only meaningful when the resolved scoring configuration matches. `score_peptides()` records this in the `preset_used` column, and `evaluate_digest()` carries the same label in `params$preset_used` so downstream reporting can distinguish named presets from custom scoring setups.

## Fixtures

The package ships pinned FASTA fixtures in `inst/extdata/` for regression tests and worked examples.

| File                               | Protein                  | Why it is useful                                  |
| ---------------------------------- | ------------------------ | ------------------------------------------------- |
| `P02769.fasta`                     | BSA                      | Stable positive-control digest                    |
| `P68431.fasta`                     | Histone H3.1             | Basic protein that exposes trypsin over-digestion |
| `P56817.fasta`                     | BACE1                    | Mixed composition with hydrophobic segments       |
| `Q8WZ42.fasta`                     | Titin                    | Very large protein for scale checks               |
| `P0CG48.fasta`                     | Ubiquitin                | Short protein edge case                           |
| `P37840_isoforms.fasta`            | Alpha-synuclein isoforms | Proteome-aware uniqueness example                 |
| `small_proteome_50_proteins.fasta` | 50-protein set           | Batch workflow fixture                            |

## Website content

The pkgdown site now covers four entry points.

- Getting started: the end-to-end pipeline and return shapes
- Enzyme selection: worked comparisons on real proteins
- Workflow presets: when to use each preset and how to modify it
- Scoring model: equations, thresholds, and current limits

## License

MIT. See [LICENSE.md](LICENSE.md).
