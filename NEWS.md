<!-- markdownlint-disable MD025 -->

# pepVet 0.0.4

## Scoring

* Reworked `S_count` to use an enzyme-aware expected peptide length instead of a fixed trypsin-derived constant.
* Added `median_peptide_length` to scoring output so enzyme comparisons show the denominator behind `S_count`.
* Added configurable `gravy_range` and `length_range` arguments across the scoring and evaluation helpers.
* Added `pepvet_preset()` with presets for standard DDA, DIA, targeted assays, membrane proteomics, FFPE-style degraded samples, and fractionated workflows.

## Documentation

* Rewrote the package website copy, docstrings, README, and articles for `v0.0.4` with detailed workflow guidance and preset-specific examples.
* Added a dedicated article on workflow presets, including practical examples and use-case notes for each preset.

## Site

* Updated pkgdown navigation to surface the preset guide alongside the core workflow articles.

# pepVet 0.0.3

## New functions

* `evaluate_digest()` wraps `digest_protein()` and `score_peptides()` into a single call.
* `compare_digests()` runs multi-enzyme comparison for a single protein, sorted by composite score.
* `recommend_enzyme()` returns the top-scoring enzyme name from a comparison.
* `batch_evaluate()` runs `evaluate_digest()` across every protein in a multi-FASTA file.
* `digest_report()` prints styled console output for evaluation and comparison results, with colour-coded bar charts and ranked tables.

## Documentation

* Rewrote `README.md` with scoring component table, reference fixture table, and workflow diagram.
* Expanded the getting-started vignette to cover all seven exported functions end-to-end.
* Added *Choosing a Proteolytic Enzyme* article covering enzyme biology, worked comparisons on BSA, Histone H3.1, and alpha-synuclein isoforms, and guidance for membrane proteins, phosphoproteomics, and IDPs.
* Added *Understanding the Scoring Model* article with mathematical definitions for all components, weight customisation guidance, verdict calibration notes, and known limitations.

## Site

* Configured pkgdown site with structured navbar, grouped reference index, and `flatly` Bootstrap 5 theme.

# pepVet 0.0.2

## Core engine

* `digest_protein()` performs cleaver-compatible digestion with validated input handling for character sequences, `AAString`, `AAStringSet`, and FASTA paths.
* `score_peptides()` computes five component scores (`S_length`, `S_coverage`, `S_count`, `S_hydro`, `S_charge`) and optional proteome-aware `S_unique`.

## Data

* Added `aa_properties` reference tibble with Kyte-Doolittle hydrophobicity, molecular weight, and side-chain pKa values.
* Added eight reference FASTA fixtures in `inst/extdata/`.

## CI

* Established GitHub Actions workflows for R CMD check and lint.
