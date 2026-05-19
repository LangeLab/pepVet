<!-- markdownlint-disable MD033 MD036 MD041 -->
<br/>
<p align="center">
  <img src="man/figures/logo.svg" alt="pepVet logo" width="96">
</p>

<h1 align="center">pepVet</h1>

<p align="center">
  Proteolytic digest evaluation for bottom-up proteomics. Score, compare, and triage enzyme choices before any sample reaches the instrument.
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-v0.1.3-2C5F8A?style=flat-square" alt="v0.1.3">
  <img src="https://img.shields.io/badge/R-%3E%3D4.6-276DC3?style=flat-square&logo=r&logoColor=white" alt="R >= 4.6">
  <img src="https://img.shields.io/badge/Bioconductor-3.23-87B13F?style=flat-square" alt="Bioconductor 3.23">
  <a href="https://github.com/LangeLab/pepVet/actions/workflows/R-CMD-check.yaml">
    <img src="https://img.shields.io/github/actions/workflow/status/LangeLab/pepVet/R-CMD-check.yaml?label=R%20CMD%20check&style=flat-square" alt="R CMD check">
  </a>
  <a href="LICENSE.md">
    <img src="https://img.shields.io/badge/license-MIT-4B9D6E?style=flat-square" alt="MIT">
  </a>
  <a href="https://langelab.github.io/pepVet/">
    <img src="https://img.shields.io/badge/docs-pkgdown-2C5F8A?style=flat-square" alt="pkgdown docs">
  </a>
</p>

---

## What pepVet does

The choice of proteolytic enzyme shapes every downstream result in a proteomics experiment. Cut too aggressively and you get thousands of tiny fragments below the detection threshold. Cut too conservatively and overlong peptides fail to fly or resolve on the column.

pepVet quantifies this trade-off before the bench work begins. Give it a protein sequence and an enzyme. It returns peptide coordinates, five orthogonal score components, a composite verdict, and a ranked enzyme comparison that makes the reasoning explicit.

## Features

**Digest simulation**

- `digest_protein()` cleaves any protein sequence with any of 40 cleaver-compatible enzyme rules and returns a peptide tibble with coordinates and missed-cleavage counts.
- `annotate_cleavage_sites()` labels each trypsin-family cleavage site as `high`, `medium`, or `low` efficiency using local P1-P1' sequence context.

**Scoring**

- `score_peptides()` summarises a peptide set into five component scores (`S_length`, `S_coverage`, `S_count`, `S_hydro`, `S_charge`) plus an optional sixth (`S_unique`) when a background proteome digest is supplied.
- `pepvet_preset()` returns workflow-specific parameter sets for DDA, DIA, targeted, membrane, FFPE/degraded, and fractionated workflows. Pass the preset list directly into `evaluate_digest()` or `score_peptides()`.

**Evaluation and comparison**

- `evaluate_digest()` wraps digest and scoring into one call and returns a named list with `$scores`, `$peptides`, and `$params`.
- `compare_digests()` runs `evaluate_digest()` across a vector of enzymes for a single protein and returns a ranked tibble.
- `recommend_enzyme()` returns the name of the best-scoring enzyme from a `compare_digests()` run.

**Batch workflows**

- `batch_evaluate()` evaluates every protein in a multi-FASTA independently and returns a flat tibble with one row per protein. All score columns, verdicts, and difficulty flags are directly accessible.
- `summarize_batch()` computes proteome-level verdict distribution, composite score statistics, per-component means, the bottom-10% proteins, and heuristic enzyme-switch candidates.
- `triage_proteins()` appends an `action` column (`proceed`, `consider_alternative`, `try_other_enzyme`, `skip`) to the batch tibble for downstream filtering.

**Reporting and export**

- `pepvet_check()` evaluates a protein and immediately prints a styled console report. One call for interactive exploration; result returned invisibly.
- `digest_report()` renders a colour-coded console summary for a single `evaluate_digest()` result or a `compare_digests()` tibble.
- `export_peptide_list()` filters valid peptides and exports them as a Skyline-compatible transition list (with computed precursor m/z), a generic annotated CSV (with GRAVY, pI, and validity), or a FASTA file.

**Peptide properties**

- `calculate_peptide_mass()` computes monoisotopic neutral mass.
- `calculate_pI()` computes isoelectric point using a Lehninger-style pKa set.

## Scoring model

Five core components, one optional proteome-aware component, one weighted composite.

| Score        | What it measures                                            | Why it matters                                              |
| ------------ | ----------------------------------------------------------- | ----------------------------------------------------------- |
| `S_length`   | Fraction of peptides inside the active length window        | Short and long peptides lower identification rates          |
| `S_coverage` | Fraction of the protein covered by valid peptides           | Dark regions weaken protein-level interpretation            |
| `S_count`    | Valid count relative to enzyme-aware expected density       | Too few weakens evidence; too many reflects over-digestion  |
| `S_hydro`    | Fraction of valid peptides inside the active GRAVY window   | Extreme hydrophobicity or hydrophilicity hurts LC behaviour |
| `S_charge`   | Valid peptides with at least one non-terminal basic residue | Proxy for multi-charge potential and richer fragment series |
| `S_unique`   | Fraction of valid peptides unique in a supplied proteome    | Shared peptides weaken protein-level attribution            |

Verdict thresholds: `Good` >= 0.65, `Moderate` >= 0.40, `Poor` < 0.40. These are heuristic ranking labels, not calibrated probabilities.

## Installation

pepVet depends on Bioconductor packages. Install the Bioconductor dependencies first.

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

## Quick start

```r
library(pepVet)

bsa <- system.file("extdata", "P02769.fasta", package = "pepVet")

# Evaluate and print in one call
pepvet_check(bsa, enzyme = "trypsin")

# Full evaluation object
ev <- evaluate_digest(bsa, enzyme = "trypsin", missed_cleavages = 1L)
ev$scores

# Compare enzymes
comp <- compare_digests(
  bsa,
  enzymes = c("trypsin", "lysc", "glutamyl endopeptidase", "asp-n endopeptidase")
)
digest_report(comp)
recommend_enzyme(bsa, enzymes = c("trypsin", "lysc", "glutamyl endopeptidase"))
```

## Batch workflow

```r
small_proteome <- system.file(
  "extdata", "small_proteome_50_proteins.fasta", package = "pepVet"
)

# One row per protein
batch <- batch_evaluate(small_proteome, enzyme = "trypsin", missed_cleavages = 1L)

# Proteome-level statistics
s <- summarize_batch(batch)
s$verdict_counts
s$component_summary

# Per-protein action recommendations
triaged <- triage_proteins(batch)
table(triaged$action)
```

## Export

```r
peps <- digest_protein(bsa, enzyme = "trypsin", missed_cleavages = 1L)

# Skyline transition list
export_peptide_list(peps, format = "skyline", charges = 2:3)
export_peptide_list(peps, format = "skyline", file = "bsa_transitions.csv")

# Generic annotated table (GRAVY, pI, valid flag)
export_peptide_list(peps, format = "generic")

# FASTA for valid peptides only
export_peptide_list(peps, format = "fasta", file = "bsa_peptides.fasta")
```

## Workflow presets

Each preset returns a list with `gravy_range`, `length_range`, and `weights` that can be passed directly into `evaluate_digest()` or `score_peptides()`.

```r
pepvet_preset("standard")

targeted <- pepvet_preset("targeted")
do.call(evaluate_digest, c(list(sequence = bsa, enzyme = "trypsin"), targeted))
```

| Preset          | Best fit                   | Key parameters                                      | Literature basis         |
| --------------- | -------------------------- | --------------------------------------------------- | ------------------------ |
| `standard`      | Routine DDA                | `[7,25]` aa, GRAVY `[-1,0.6]`, AHP weights          | Tabb 2008                |
| `dia`           | DIA and SWATH              | `[7,30]` aa, GRAVY `[-1,0.8]`, high coverage weight | Ludwig 2018              |
| `targeted`      | SRM, PRM, MRM              | `[8,20]` aa, GRAVY `[-0.8,0.4]`, S_unique 30%       | Lange 2008, Picotti 2012 |
| `membrane`      | Hydrophobic proteins       | GRAVY `[-1.0,2.0]`, S_hydro 5%                      | Vit & Petrak 2017        |
| `ffpe_degraded` | Degraded material          | `[6,30]` aa, high S_count weight                    | Coscia 2020, Buczak 2023 |
| `fractionated`  | SCX or high-pH RP planning | Same as standard, `include_pI = TRUE`               | —                        |

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

## Scope

pepVet is not a peptide detectability predictor. It is a rule-based, multi-criteria digest-ranking model for pre-acquisition planning. Scores are interpretable rankings within a given enzyme-workflow combination, not calibrated probabilities.

## Documentation and support

- Website: [langelab.github.io/pepVet](https://langelab.github.io/pepVet/)
- Changelog: [NEWS.md](NEWS.md)
- Bug reports and questions: [GitHub Issues](https://github.com/LangeLab/pepVet/issues)

## Citation

```r
citation("pepVet")
```

## License

MIT. See [LICENSE.md](LICENSE.md).

## Contributing

Pull requests, bug reports, and documentation fixes are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for the review workflow and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) for community standards.
