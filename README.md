# pepVet

`pepVet` is a Bioconductor-oriented R package for simulating proteolytic digests and evaluating peptide sets for proteomics workflows.

## Status

This repository now includes a working digestion layer plus pinned reference fixtures and validation tests. The scoring and evaluation layers are still to come, but the package can already normalize protein input, digest with cleaver-compatible rules, and return exact peptide coordinates.

Current implementation highlights:

- digest simulation via `digest_protein()`
- input handling for character sequences, `AAString`, `AAStringSet`, and FASTA
 paths
- pinned reference FASTA fixtures in `inst/extdata/`
- exact start/end coordinate tracking, including missed cleavages
- amino-acid property reference data in `aa_properties`
- installed-package validation with `devtools::test()`, `lintr`, `pkgdown`,
 and `R CMD check`

## Quick Start

```r
library(pepVet)

digest_protein("MKWVTFISLLFLFSSAYSR")
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

## Next steps

1. Implement the scoring engine on top of the digest output.
2. Add digest evaluation and comparison helpers.
3. Expand the pkgdown site with scoring and evaluation walkthroughs.
