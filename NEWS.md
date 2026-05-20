<!-- markdownlint-disable MD025 MD024 -->

# pepVet 0.1.4

> Current working version. Not yet released.

## Housekeeping

* Standardised cli error messages across the package. Replaced all `paste0("{.arg ", var, "}")` patterns with `"{.arg {var}}"` cli glue syntax in `.abort()`, `cli_warn()`, and `cli_inform()` calls (G-06). Also converted `"i"` bullet `paste()` calls to `{.val {x}}` formatters.

* Added `.bind_rows()` helper and replaced all bare `do.call(rbind, ...)` calls with it (G-01). The helper returns an empty tibble for empty input instead of failing silently.

# pepVet 0.1.3

## Defaults

* Changed default `missed_cleavages` from `0L` to `1L` across all evaluation functions. MC=1 reflects standard bottom-up proteomics practice; MC=0 produces unrealistically poor scores for every enzyme.
* Changed default scoring weights to AHP-derived values: `S_length = 0.200, S_coverage = 0.348, S_count = 0.226, S_hydro = 0.138, S_charge = 0.088`. Weights were derived via Analytic Hierarchy Process with pairwise comparisons grounded in the proteomics literature (CR = 0.028).
* Lowered the Good verdict threshold from `0.70` to `0.65` to better align with the observed score distribution under the new weights. Centralised verdict thresholds in `.pepvet_params` accessed via `.get_param()`.

## Performance

* `digest_protein()` rewritten for batch workloads. `.cleavage_ranges()` is now called once on the full `AAStringSet` instead of once per protein, eliminating repeated S4 dispatch overhead. The non-efficiency path pre-allocates six output vectors and fills them in a single loop before constructing one tibble, replacing ~20 K individual tibble builds followed by `do.call(rbind, ...)`.
* `batch_evaluate()` restructured around two bulk calls — one `digest_protein()` and one `score_peptides()` for the entire input — instead of a per-protein loop over `evaluate_digest()`. A new internal helper `.batch_difficulty_flags()` computes all four difficulty flags across all proteins simultaneously via `tabulate()`, `tapply()`, and `.calculate_gravy_vec()`, replacing ~20 K per-protein subset/GRAVY/`tibble::as_tibble()` cycles.
* `triage_proteins()` fully vectorized. Row-by-row `vapply` logic replaced with nested `ifelse()` and a `rowSums(matrix < 0.5)` component check, making the function O(n) in a single pass.
* `S_coverage` scoring no longer uses `IRanges::reduce()`. Overlapping intervals are now merged inline via `order()` + `cummax()` on sorted starts and ends — same result with no S4 dispatch.
* Scoring helpers (`S_coverage`, `S_count`, `S_hydro`, `S_charge`) accept a pre-computed `valid_digest` argument so callers can extract valid peptides once and share it across all four components instead of calling `.extract_valid_digest()` four times per protein.

## New functions

* Added `batch_compare_enzymes()` to score an entire proteome against multiple enzymes in one call, returning a tidy tibble of class `pepvet_batch_comparison` with one row per protein–enzyme pair. Parallel execution is supported via `parallel::mclapply` on Unix (fork copy-on-write, zero serialization overhead) and `parallel::parLapply` on Windows; both are part of base R. Proteins are split into equal chunks, so all `cores` workers are utilised regardless of the number of enzymes.

## Parallel robustness

* Fixed a bug in `batch_evaluate()` where `parallel::mclapply` worker crashes produced silently corrupted results. Failed chunks are now detected and retried sequentially.

## Amino acid data

* Added `U` (selenocysteine, Sec) to the `aa_properties` reference table. 25 human proteins contain selenocysteine and were previously rejected during input validation.
* Added `O` (pyrrolysine, Pyl) to the `aa_properties` reference table (row 22). Monoisotopic mass verified from PubChem CID 119813 and ChEBI CHEBI:91273 (C₁₂H₂₁N₃O₃, 255.15829 Da). Hydrophobicity and pKa are `NA`: the Kyte–Doolittle scale predates the discovery of pyrrolysine (1982 vs 2002), and the ε-amino group is tied up in an amide bond and is not titratable. `calculate_gravy()` now passes `na.rm = TRUE` so sequences containing O return a valid GRAVY score computed from the remaining residues.

## Removals

* Removed `plot_enzyme_protein_heatmap()`. The 2D tile matrix did not add enough over the other comparison functions (`plot_enzyme_comparison()`, `plot_batch_comparison()`) to justify maintenance.
* Removed `plot_component_scatter()`. The 2D scatter of component scores was redundant with the information already visible in `plot_proteome_overview()` and `plot_batch_comparison()`.
* Removed `plot_batch_summary()`. The two-panel overview was superseded by the richer `plot_proteome_overview()` and `plot_batch_comparison()`.
* Removed `plot_protein_comparison()`. The grouped bar chart did not add enough insight beyond `plot_enzyme_comparison()` and the batch-level comparison functions.

# pepVet 0.1.2

## Visualization

* Added `plot_digest_profile()` as the flagship single-protein diagnostic. Four-panel layout showing length distribution, GRAVY scatter, coverage map, and component scores.
* Added `plot_coverage_map()` for horizontal sequence coverage visualisation with valid/invalid peptide segments and missed-cleavage expansion lanes.
* Added `plot_cleavage_map()` for vertical cleavage-site ticks with fragment blocks and optional efficiency coloring from `annotate_cleavage_sites()`.
* Added `plot_enzyme_comparison()` for comparing component scores across multiple enzymes with sorted bars and recommendation badge.
* Added `plot_protein_comparison()` for comparing component scores across multiple proteins under a single enzyme with verdict badges.
* Added `plot_enzyme_protein_heatmap()` for a 2D tile matrix of proteins versus enzymes with composite-score gradient and verdict overlay.
* Added `plot_length_distribution()` for peptide-length histograms with valid-range shading and multi-input faceted mode.
* Added `plot_gravy_landscape()` for 2D scatter of length versus GRAVY with valid-region rectangle and multi-input comparison.
* Added `plot_pI_distribution()` for isoelectric-point histograms with fraction-bin coloring and multi-input density overlay.
* Added `plot_missed_cleavage_impact()` for showing how component scores change across MC=0/1/2.
* Added `plot_batch_summary()` for two-panel proteome overview: verdict histogram and score-versus-length scatter.
* Added `plot_proteome_heatmap()` for hierarchical clustering of component scores across a batch. Requires `pheatmap` (guarded).
* Added `plot_component_scatter()` for 2D scatter of any two score components across all batch proteins.
* Added shared infrastructure in `R/plot_utils.R`: `.pepvet_pal` colour palette, `.pepvet_theme()` base theme, and reusable panel builders.

## Housekeeping

* Established coding conventions and standardisation strategy for consistency across the package.

# pepVet 0.1.1

## Batch evaluation, triage, and export

* `batch_evaluate()` now returns a flat tibble with one row per protein instead of a named list of `evaluate_digest()` results. Columns: `protein_id`, `protein_length`, `n_peptides`, `n_valid_peptides`, all component scores, `composite_score`, `verdict`, `median_peptide_length`, and four sequence-level difficulty flags. This is a breaking change for any code that accessed results via `batch[[protein_id]]$scores`.
* Added `summarize_batch()` to extract aggregate statistics from the batch tibble: verdict distribution, composite score distribution, per-component means, bottom-10% problem proteins, and heuristic enzyme-switch candidates.
* Added `triage_proteins()` to append an `action` column (`"proceed"`, `"consider_alternative"`, `"try_other_enzyme"`, `"skip"`) to the batch tibble based on verdict and sequence-level difficulty flags.
* Added `pepvet_check()` as a single-call convenience wrapper: evaluates a protein digest and immediately prints a styled console report. Returns the `evaluate_digest()` result invisibly.
* Added `export_peptide_list()` in a new `R/export.R` module. Exports valid peptides in `"skyline"` (precursor-charge CSV with computed m/z for Skyline import), `"generic"` (annotated CSV with GRAVY, pI, and validity), or `"fasta"` format. Writes to a file when `file` is specified; returns a tibble or character vector otherwise.
* Fixed `export_peptide_list()` generic format to include a `pI` column alongside `gravy` and `valid`.

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
