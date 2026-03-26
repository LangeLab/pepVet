# pepVet 0.0.3

* Added `digest_report()` — styled console output for `evaluate_digest()` and
  `compare_digests()` results, with colour-coded bar charts and ranked
  comparison tables.
* Added `evaluate_digest()` — single-call wrapper combining `digest_protein()`
  and `score_peptides()`.
* Added `compare_digests()` — multi-enzyme comparison for a single protein,
  sorted by composite score.
* Added `recommend_enzyme()` — returns the top-scoring enzyme name from a
  comparison.
* Added `batch_evaluate()` — runs `evaluate_digest()` across every protein in
  a multi-FASTA file.
* Rewrote `README.md` with a clear introduction, scoring component table,
  reference fixture table, and workflow diagram.
* Substantially expanded the getting-started vignette to cover all seven
  exported functions end-to-end.
* Added *Choosing a Proteolytic Enzyme* article covering biology, worked
  comparisons on BSA, Histone H3.1, and alpha-synuclein isoforms, and guidance
  for membrane proteins, phosphoproteomics, and IDPs.
* Added *Understanding the Scoring Model* article with full mathematical
  definitions for all five (six) components, weight customisation guidance,
  verdict calibration notes, and known limitations.
* Configured pkgdown site with a structured navbar, grouped reference index,
  and `flatly` Bootstrap 5 theme.

# pepVet 0.0.2

* Added protein-level scoring engine: `score_peptides()` with five component
  scores (`S_length`, `S_coverage`, `S_count`, `S_hydro`, `S_charge`) and
  optional proteome-aware `S_unique`.
* Added amino acid reference data (`aa_properties`) with Kyte-Doolittle
  hydrophobicity, molecular weight, and side-chain pKa values.
* Added `digest_protein()` — cleaver-compatible digestion with validated input
  handling for character sequences, `AAString`, `AAStringSet`, and FASTA paths.
* Added eight reference FASTA fixtures in `inst/extdata/`.
* Established GitHub Actions CI workflows for R CMD check and lint.
