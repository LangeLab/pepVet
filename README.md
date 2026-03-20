# pepVet

`pepVet` is a Bioconductor-oriented R package for simulating proteolytic digests and evaluating peptide sets for proteomics workflows.

## Overview

Version `0.0.2` provides a usable core workflow for peptide-centric method development. The package can normalize protein input, simulate cleaver-compatible digests, and score peptide sets at the protein level with validated component metrics.

Current capabilities:

- digest simulation via `digest_protein()`
- peptide-set scoring via `score_peptides()`
- input handling for character sequences, `AAString`, `AAStringSet`, and FASTA paths
- pinned reference FASTA fixtures in `inst/extdata/`
- exact start/end coordinate tracking, including missed cleavages
- amino-acid property reference data in `aa_properties`
- component scores for length, coverage, count, hydrophobicity, charge, and optional proteome uniqueness
- full validation through `devtools::test()`, `lintr`, `pkgdown`, and `R CMD check`

## Quick Start

```r
library(pepVet)

digest_result <- digest_protein("MKWVTFISLLFLFSSAYSR")
score_peptides(digest_result)
```

```r
library(Biostrings)

digest_protein(AAString("MKWVTFISLLFLFSSAYSR"), enzyme = "trypsin")
```

```r
digest_protein(
 system.file("extdata", "P02769.fasta", package = "pepVet"),
 enzyme = "lysc",
 missed_cleavages = 1L
)
```

```r
bsa_digest <- digest_protein(
 system.file("extdata", "P02769.fasta", package = "pepVet"),
 enzyme = "trypsin"
)

score_peptides(bsa_digest)
```

```r
proteome_digest <- digest_protein(
 c(target = "AAAAAAARAAAAAAAK", background = "AAAAAAARGGGGGGGK")
)

score_peptides(
 proteome_digest[proteome_digest$protein_id == "target", ],
 proteome = proteome_digest
)
```

## Supported Inputs

- single character sequence
- named character vector of sequences
- `Biostrings::AAString`
- `Biostrings::AAStringSet`
- FASTA file path, including multi-entry fixtures and irregular but valid FASTA
 headers

## Digest Output

`digest_protein()` returns a tibble with:

- `protein_id`
- `peptide`
- `start`
- `end`
- `length`
- `missed_cleavages`

The implementation uses cleaver-compatible strict cut rules and expands missed cleavages inside pepVet so repeated peptides and overlapping coordinates remain exact.

## Scoring Output

`score_peptides()` returns one row per `protein_id` with:

- `S_length`
- `S_coverage`
- `S_count`
- `S_hydro`
- `S_charge`
- optional `S_unique` when a proteome digest is supplied
- `composite_score`
- `verdict`

Protein-only scoring uses default weights of `0.25/0.25/0.20/0.15/0.15` for length, coverage, count, hydrophobicity, and charge. Proteome-aware scoring adds `S_unique` and switches to `0.20/0.20/0.15/0.15/0.10/0.20`.

## Reference Fixtures

The package ships pinned FASTA fixtures in `inst/extdata/` for:

- BSA (`P02769.fasta`)
- lysozyme C (`P00698.fasta`)
- beta-secretase 1 (`P56817.fasta`)
- ubiquitin (`P0CG48.fasta`)
- histone H3.1 (`P68431.fasta`)
- titin (`Q8WZ42.fasta`)
- alpha-synuclein isoforms (`P37840_isoforms.fasta`)
- a 50-protein human fixture (`small_proteome_50_proteins.fasta`)

## Roadmap

1. Add evaluation, comparison, and recommendation helpers on top of the digest and scoring layers.
2. Add batch workflows and proteome-aware summaries for multi-protein inputs.
3. Expand console reporting and documentation around score interpretation.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details.

---
