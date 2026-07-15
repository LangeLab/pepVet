<!-- markdownlint-disable MD025 MD024 -->

# pepVet 0.99.0

## User-visible changes

* Redesigned interactive digest reports with ASCII-safe score profiles, explicit scoring settings and peptide counts, wrapped identifiers, and width-aware enzyme comparisons.
* Preserved resolved scoring metadata through single-protein, comparison, batch, sensitivity, and plotting workflows.
* Restored Windows protein-chunk parallelism through base R socket workers. Unix uses fork workers, and failed chunks on either backend are retried sequentially with a classed warning.

## Validation and fixes

* Standardized malformed-input handling across digestion, scoring, evaluation, export, reporting, diagnostics, and plotting. Invalid inputs now fail through package conditions before reaching base R errors.
* Rejected duplicate protein identifiers and duplicate data-frame columns before they can collapse or corrupt downstream results.
* Applied the zero-count digest rule consistently: digests without peptides inside the active length range receive a zero composite score and a Poor verdict.
* Corrected the workflow-preset comparison to apply complete preset settings against an explicit three-protein digest background.
* Distinguished theoretical products, length-valid peptides, combined comparison windows, model score bands, and model-based rankings throughout package-facing documentation.
* Fixed missed-cleavage filtering in peptide-overlap plots when all levels are requested, and validated plot titles and figure dimensions before rendering.
* Corrected weight sensitivity tie handling, corner diagnostics, batch aggregation, and single-result plot semantics. Results now record simulation settings, and batch memory is bounded by a smaller default chunk.

## Testing

* Reworked the test suite around public contracts, scientific invariants, classed conditions, metadata, plot semantics, and platform-specific batch execution.
* Added generated-sequence properties and a malformed-input audit covering all 37 exported functions.
* Added fast contract tests for shipped FASTA, tool-comparison, and PeptideAtlas artifacts. Data-generation workflows now run as separate offline audits and fail clearly when required inputs are missing.
* Removed an architecture-dependent threshold collision from the diagnostics ablation oracle so installed-package tests retain the same verdict-flip expectation across BLAS implementations.
* Corrected Windows-specific test assumptions about normalized path separators and avoided reinstalling an already installed package during `R CMD check`.

## Documentation and metadata

* Prepared the package metadata and installation guidance for initial Bioconductor submission.
* Reconciled the README, vignettes, help pages, package metadata, and generated site with the current scoring and validation contracts.
* Removed the redundant custom package citation because pepVet does not yet have an associated publication DOI. Installed citations now come from DESCRIPTION metadata.
* Described the scoring weights as documented expert priors because their original derivation record is not retained. The numeric weights are unchanged.

## Release infrastructure

* Added Linux coverage reporting with a 90 percent floor, retained coverage artifacts, and Codecov upload support.
* Added pkgdown builds for pull requests and GitHub Pages deployment from `main`, with bounded build stages and cancellation of obsolete runs.
* Added a source-package artifact build alongside parallel release-R checks on Linux, macOS, and Windows, with unfinished checks cancelled after the first check failure.
* Split CI installation into role-specific `pak` dependency plans and preserved resolved-library caches when later workflow stages fail.
* Bounded dependency installation in every workflow so a stalled package resolver fails promptly instead of consuming an entire job timeout.
* Consolidated routine Linux build, check, coverage, and pkgdown work behind one uncached dependency installation while retaining parallel macOS and Windows checks. This avoids intermittent cache-service stalls across duplicated Linux setup jobs.
* Added version-tag release-candidate builds and tag/published-release BiocCheck workflows. Release artifacts are built in temporary storage and are not published automatically.

# pepVet 0.1.7

## Empirical concordance

* Added a reproducible PeptideAtlas concordance analysis. In the sampled human tryptic digests, peptides inside the default length and GRAVY windows had a 38.0 percentage-point higher mean observation rate than peptides outside them. The result is limited to the sampled data and is not an experimental calibration of pepVet scores.
* Added bounded grid-search and threshold-calibration analyses, documented their limits, and retained the current defaults as strict, conservative settings by design.
* Corrected the concordance artifacts to store full protein lengths and regenerated the committed CSV and RDS files.

## Package quality

* Standardized docstrings, limitations, error handling, plot conventions, package-qualified calls, comments, and line formatting across the source.
* Unified GRAVY calculation under one vectorized helper and removed stale helper references from package and data-generation code.
* Removed inactive GitHub workflow files and generated vignette data from version control so the repository matches the current release process.

## Fixes

* Fixed edge cases in row binding, cleavage-efficiency summaries, GRAVY missing-value handling, and sensitivity plot validation.
* Fixed `recommend_enzyme()` tie handling so all enzymes within tolerance are returned in alphabetical order.
* Corrected charge scoring for single-residue peptides so terminal basic residues are not counted as internal charge sites.
* Fixed plot configuration updates so invalid changes leave the active configuration unchanged.

# pepVet 0.1.6

## New functions

* `score_diagnostics()` runs VIF, PCA, and ablation analyses on a `batch_evaluate()` result to quantify multicollinearity, dimensionality, and component importance in the scoring model.
* `plot_score_diagnostics()` visualises the diagnostics result as a three-panel figure: VIF bar chart with severity-colored bars, PCA scree plot with cumulative variance, and ablation waterfall with error bars and verdict-flip counts.

## New documentation

* Added *Score Diagnostics for pepVet Scoring Models* article covering VIF, PCA, and ablation analysis with worked examples.
* Added *pepVet in the Tool Landscape* article comparing pepVet against MS-Digest, ExPASy PeptideMass, Protein Cleaver, ProteaseGuru, and PeptideRanger.
* Added diagnostics cross-reference to the introduction vignette.

## Bug fixes

* Fixed `sensitivity_analysis()` batch-mode verdict instability: now uses full 3-level classification (Good/Moderate/Poor) instead of checking only the Good boundary.
* Fixed vignette pipe example where `recommend_enzyme()` received a comparison tibble instead of a protein sequence.
* Fixed tool-comparison vignette: replaced fragile `../inst/extdata/` paths with `system.file()`.
* Fixed verdict threshold table in introduction vignette: said `0.7 / 0.4`, code uses `0.65 / 0.40`.

## Housekeeping

* Refactored `zzz.R` to avoid the `lockBinding` R CMD check NOTE while preserving active bindings for mutable plot configuration.
* Fixed README installation instructions (removed `S4Vectors` from BiocManager install; pulled in automatically by Biostrings).
* Enabled automatic CI triggers for R CMD check on push and PR across ubuntu, macOS, and Windows.
* Added `.gitattributes`, `.Rbuildignore` updates, `.onAttach` removal, and cli cleanup.
* Docstring standardisation and `.bind_rows()` type-safety improvements.

# pepVet 0.1.5

## New functions

* `sensitivity_analysis()` performs Monte Carlo weight perturbation via Dirichlet sampling. Reports simulated verdict frequencies, composite intervals, rank stability (top-1, Kendall tau), per-protein instability in batch mode, weight importance (R squared), and corner-case composites. Optional `importance` and `corner_cases` diagnostics.
* `plot_weight_sensitivity()` generates a verdict-coloured density ridge plot from a sensitivity analysis result, with threshold lines and a rug mark at the default composite score.

## Scoring

* Zero-cleavage hard-fail: proteins with no cleavage sites for a given enzyme now receive composite = 0 and verdict = "Poor" regardless of other scores. Previously an uncuttable protein with ideal length/hydro/charge could score up to 0.426, crossing the Moderate threshold. The fix belongs in the scoring engine so sensitivity analysis treats these as deterministic negative controls.

## Plot fixes

* `plot_digest_profile()` coverage panel now filters to the exact requested MC level instead of showing all levels. Overlapping peptides are stacked into sub-tracks via greedy interval packing, and the panel title includes the MC level. Gap overlays remain but the subtitle no longer counts uncovered regions; a small italic caption notes when some gaps are too narrow to render at plot scale.

## Site

* Updated pkgdown color palette to match the new pepVet logo colors. Primary/secondary swapped for better link contrast. Favicons regenerated.
* Replaced all logo references with the new `pepVet-logo.png`. Navbar shadow, rounded code blocks, sticky TOC sidebar added.
* Flattened the Articles dropdown. Removed redundant "Getting Started" entry and section headers.
* Added prev/next navigation buttons to reference pages (client-side JS).

# pepVet 0.1.4

## Housekeeping

* Standardised cli error messages across the package. Replaced all `paste0("{.arg ", var, "}")` patterns with `"{.arg {var}}"` cli glue syntax in `.abort()`, `cli_warn()`, and `cli_inform()` calls. Also converted `"i"` bullet `paste()` calls to `{.val {x}}` formatters.
* Added `.bind_rows()` helper and replaced all bare `do.call(rbind, ...)` calls with it. The helper returns an empty tibble for empty input instead of failing silently.
* Replaced `import(cleaver)` with `@importFrom cleaver cleavageRanges`. Removed the `.cleavage_ranges` workaround via `get("cleavageRanges", envir = asNamespace("cleaver"))` and now calls `cleaver::cleavageRanges()` directly.
* Removed dead code: `.compute_difficulty_flags()` in `evaluation.R` was defined but never called. Its logic was superseded by `.batch_difficulty_flags()`.
* Made `plot_coverage_map()` gradient stops (`color_by = "hydrophobicity"`) data-driven instead of hardcoding the GRAVY valid-range boundary at 0.6. The 4 color stops are now evenly spaced across the actual GRAVY range of the displayed peptides.
* Suppressed three R CMD check notes by adding `.lintr`, `tmp/`, and `paper/` to `.Rbuildignore`.
* Added missing `@return` roxygen tags to 10 internal helper functions across the source.
* Removed duplicate `.classify_verdict()` from `plot_utils.R`. The function was defined in both `scoring.R` and `plot_utils.R` with identical logic; the duplicate silently shadowed the original due to file-sourcing order.
* Standardised patchwork title and tag sizes across all plot functions. Added `patchwork_tag_size` to `.pepvet_params`. Fixed `plot_gravy_landscape()` title size 13 -> 15. Changed `plot_enzyme_comparison()`, `plot_digest_profile()`, `plot_proteome_overview()`, `plot_batch_comparison()` to use `.get_param()` for tag sizes instead of hardcoded 14.
* Fixed `plot_batch_comparison()` heatmap colorbar to use explicit `unit(..., "pt")` instead of numeric (lines) values.

## Bug fixes

* Fixed `.onLoad()` failure when called on an already-loaded namespace. The `rm()` call on locked namespace bindings was guarded with `bindingIsActive()` checks so the active-binding installation runs exactly once.

# pepVet 0.1.3

## Defaults

* Changed default `missed_cleavages` from `0L` to `1L` across all evaluation functions. MC=1 reflects standard bottom-up proteomics practice; MC=0 produces unrealistically poor scores for every enzyme.
* Changed default scoring weights to documented expert-prior values: `S_length = 0.200, S_coverage = 0.348, S_count = 0.226, S_hydro = 0.138, S_charge = 0.088`. The derivation record needed to support the earlier method and consistency-ratio description was not retained, so those claims have been removed.
* Lowered the Good verdict threshold from `0.70` to `0.65` to better align with the observed score distribution under the new weights. Centralised verdict thresholds in `.pepvet_params` accessed via `.get_param()`.

## Performance

* `digest_protein()` rewritten for batch workloads. `.cleavage_ranges()` is now called once on the full `AAStringSet` instead of once per protein, eliminating repeated S4 dispatch overhead. The non-efficiency path pre-allocates six output vectors and fills them in a single loop before constructing one tibble, replacing ~20 K individual tibble builds followed by `do.call(rbind, ...)`.
* `batch_evaluate()` restructured around two bulk calls, one `digest_protein()` and one `score_peptides()` for the entire input, instead of a per-protein loop over `evaluate_digest()`. A new internal helper `.batch_difficulty_flags()` computes all four difficulty flags across all proteins simultaneously via `tabulate()`, `tapply()`, and a vectorized GRAVY calculation, replacing ~20 K per-protein subset/GRAVY/`tibble::as_tibble()` cycles.
* `triage_proteins()` fully vectorized. Row-by-row `vapply` logic replaced with nested `ifelse()` and a `rowSums(matrix < 0.5)` component check, making the function O(n) in a single pass.
* `S_coverage` scoring no longer uses `IRanges::reduce()`. Overlapping intervals are now merged inline via `order()` + `cummax()` on sorted starts and ends. Same result with no S4 dispatch.
* Scoring helpers (`S_coverage`, `S_count`, `S_hydro`, `S_charge`) accept a pre-computed `valid_digest` argument so callers can extract valid peptides once and share it across all four components instead of calling `.extract_valid_digest()` four times per protein.

## New functions

* Added `batch_compare_enzymes()` to score an entire proteome against multiple enzymes in one call, returning a tidy tibble of class `pepvet_batch_comparison` with one row per protein-enzyme pair. Parallel execution is supported via `parallel::mclapply` on Unix (fork copy-on-write, zero serialization overhead) and `parallel::parLapply` on Windows; both are part of base R. Proteins are split into equal chunks, so all `cores` workers are utilised regardless of the number of enzymes.

## Parallel robustness

* Fixed a bug in `batch_evaluate()` where `parallel::mclapply` worker crashes produced silently corrupted results. Failed chunks are now detected and retried sequentially.

## Amino acid data

* Added `U` (selenocysteine, Sec) to the `aa_properties` reference table. 25 human proteins contain selenocysteine and were previously rejected during input validation.
* Added `O` (pyrrolysine, Pyl) to the `aa_properties` reference table (row 22). Monoisotopic mass verified from PubChem CID 119813 and ChEBI CHEBI:91273 (C₁₂H₂₁N₃O₃, 255.15829 Da). Hydrophobicity and pKa are `NA`: the Kyte-Doolittle scale predates the discovery of pyrrolysine (1982 vs 2002), and the ε-amino group is tied up in an amide bond and is not titratable. `calculate_gravy()` now passes `na.rm = TRUE` so sequences containing O return a valid GRAVY score computed from the remaining residues.

## Removals

* Removed `plot_enzyme_protein_heatmap()`. The 2D tile matrix did not add enough over the other comparison functions (`plot_enzyme_comparison()`, `plot_batch_comparison()`) to justify maintenance.
* Removed `plot_component_scatter()`. The 2D scatter of component scores was redundant with the information already visible in `plot_proteome_overview()` and `plot_batch_comparison()`.
* Removed `plot_batch_summary()`. The two-panel overview was superseded by the richer `plot_proteome_overview()` and `plot_batch_comparison()`.
* Removed `plot_protein_comparison()`. The grouped bar chart did not add enough insight beyond `plot_enzyme_comparison()` and the batch-level comparison functions.

# pepVet 0.1.2

## Visualization

* Added `plot_digest_profile()` as a four-panel single-protein diagnostic showing length distribution, GRAVY scatter, coverage map, and component scores.
* Added `plot_coverage_map()` for horizontal sequence coverage visualisation with valid/invalid peptide segments and missed-cleavage expansion lanes.
* Added `plot_cleavage_map()` for vertical cleavage-site ticks with fragment blocks and optional efficiency coloring from `annotate_cleavage_sites()`.
* Added `plot_enzyme_comparison()` for comparing component scores across multiple enzymes with sorted bars and a top-score badge.
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
