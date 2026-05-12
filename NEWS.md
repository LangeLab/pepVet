<!-- markdownlint-disable MD025 MD024 -->

# pepVet (development version)

## Batch evaluation, triage, and export

* `batch_evaluate()` docs updated to reference `summarize_batch()` and `triage_proteins()` as the aggregate and triage layer.
* Added `summarize_batch()` to extract aggregate statistics from a `batch_evaluate()` result: verdict distribution, score distribution, per-component averages, bottom-10% problem proteins, and heuristic enzyme-switch candidates.
* Added `triage_proteins()` to convert a `batch_evaluate()` result into a flat per-protein tibble with an `action` column (`"proceed"`, `"consider_alternative"`, `"try_other_enzyme"`, `"skip"`) derived from verdict and sequence-level difficulty flags.
* Added `pepvet_check()` as a single-call convenience wrapper that evaluates a digest and immediately prints a styled console report.
* Added `export_peptide_list()` in a new `R/export.R` module. Exports valid peptides in `"skyline"` (precursor-charge CSV for Skyline import), `"generic"` (annotated CSV with GRAVY and validity), or `"fasta"` format. Writes to a file when `file` is specified; returns a tibble or character vector otherwise.

# pepVet 0.1.0

## Digestion and annotations

* Added `annotate_cleavage_sites()` for sequence-local cleavage-efficiency annotation of trypsin-family digestion sites.
* Added optional `cleavage_efficiency` output in `digest_protein()` so peptide tables can carry high/medium/low cleavage-risk context without changing the default schema.
* Added `n_high_efficiency_sites` and `n_low_efficiency_sites` to `evaluate_digest()` output as informational protein-level summaries.

## Scoring and utilities

* Added peptide mass and pI utilities plus residue-level monoisotopic mass data for fractionation-aware planning workflows.
* Retained preset tracking through `preset_used` so score tables can distinguish shipped workflows from custom configurations.

## Documentation

* Reworked the README and core articles for the `0.1.0` release so digestion, scoring, presets, and cleavage annotations are explained as one coherent planning workflow.

# pepVet 0.0.4

## Scoring

* Reworked `S_count` to use an enzyme-aware expected peptide length instead of a fixed trypsin-derived constant.
* Added `median_peptide_length` to scoring output so enzyme comparisons show the denominator behind `S_count`.
* Added configurable `gravy_range` and `length_range` arguments across the scoring and evaluation helpers.
* Added `pepvet_preset()` with presets for standard DDA, DIA, targeted assays, membrane proteomics, FFPE-style degraded samples, and fractionated workflows.
* Clarified that `S_charge` tracks extra internal basic-residue richness, not baseline peptide ionizability.
* Added `preset_used` metadata to scoring output and `evaluate_digest()` params so named presets and custom scoring configurations can be distinguished explicitly.

## Documentation

* Rewrote the package website copy, docstrings, README, and articles for `v0.0.4` with detailed workflow guidance and preset-specific examples.
* Added a dedicated article on workflow presets, including practical examples and use-case notes for each preset.
* Added scoring-model positioning, evidence-basis, scope, and known-limitations guidance so heuristic and literature-backed assumptions are documented in-package.

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
