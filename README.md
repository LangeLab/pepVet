# pepVet

pepVet is an R package for proteolytic digest simulation, peptide-set scoring, enzyme comparison, and workflow preset selection for bottom-up proteomics. Give it a protein sequence and an enzyme. It returns peptide coordinates, score components, a ranked comparison, and enough detail to explain why one enzyme is a better starting point than another.

## What pepVet does

The choice of proteolytic enzyme shapes every downstream result in a proteomics experiment. Cut too aggressively and you get thousands of tiny fragments below the detection threshold. Cut too conservatively and overlong peptides fail to fly or resolve on the column.

pepVet quantifies this trade-off with five orthogonal scoring components and a weighted composite score.

## Version 0.0.4

Version 0.0.4 changes two parts of the scoring model.

- `S_count` now uses an enzyme-aware expected peptide length instead of a fixed trypsin-derived denominator.
- Workflow presets now configure peptide-length windows, GRAVY windows, and weights for common experiment types.

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

## Workflow presets

Preset support is part of `v0.0.4`. Each preset returns a named list with `gravy_range`, `length_range`, and `weights`, so you can pass it straight into `evaluate_digest()` or `score_peptides()`.

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

| Score        | What it measures                                                  | Why it matters                                                               |
| ------------ | ----------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `S_length`   | Fraction of peptides inside the active peptide-length window      | Very short and very long peptides lower identification rates                 |
| `S_coverage` | Fraction of the parent protein covered by valid peptides          | Dark regions weaken protein-level interpretation                             |
| `S_count`    | Valid peptide count relative to the enzyme-aware expected density | Too few peptides weaken evidence. Too many reflect over-digestion            |
| `S_hydro`    | Fraction of valid peptides inside the active GRAVY window         | Extreme hydrophobicity or hydrophilicity hurts LC behavior                   |
| `S_charge`   | Fraction of valid peptides with at least one non-terminal basic residue | Higher values indicate more opportunities for multi-charge states and richer fragment ion series |
| `S_unique`   | Fraction of valid peptides unique in a supplied proteome          | Shared peptides weaken protein-level attribution                             |

`S_charge` does not mean a peptide can or cannot ionize. Tryptic peptides still carry the free N-terminus and often a terminal Lys or Arg. The score is meant to distinguish baseline ionizability from extra internal basic-residue richness, which tends to support higher charge states and richer b/y ion series.

The composite score is a weighted sum. Verdict thresholds remain `Good` at `>= 0.70`, `Moderate` at `>= 0.40`, and `Poor` below `0.40`. These thresholds are heuristic. They are not calibrated probabilities. [REF]

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
