# pepVet

**pepVet** is a Bioconductor-oriented R package for simulating proteolytic
digests and evaluating peptide sets for proteomics workflows. Given a protein
sequence and an enzyme, pepVet tells you how suitable the resulting peptides
are for downstream LC-MS/MS detection — and helps you pick the best enzyme
before you ever run a gel.

## What pepVet does

Choosing a proteolytic enzyme is one of the earliest and most consequential
decisions in a proteomics experiment. Cut too aggressively and you get
thousands of tiny, undetectable fragments. Cut too conservatively and large
peptides fail to fly or resolve on the column. pepVet quantifies this
trade-off with five orthogonal scoring components and a weighted composite
score that you can act on immediately.

```r
library(pepVet)

bsa <- system.file("extdata", "P02769.fasta", package = "pepVet")

# One-step evaluate
evaluate_digest(bsa, enzyme = "trypsin", missed_cleavages = 1L)

# Compare all candidate enzymes side by side
comp <- compare_digests(bsa,
  enzymes = c("trypsin", "lysc", "glutamyl endopeptidase",
              "asp-n endopeptidase", "chymotrypsin-high"))

# Print a styled report to the console
digest_report(comp)

# Get the winner programmatically
recommend_enzyme(bsa, enzymes = c("trypsin", "lysc", "glutamyl endopeptidase"))
```

## Installation

pepVet depends on Bioconductor infrastructure. Install those first, then
install pepVet from GitHub.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Biostrings", "IRanges", "S4Vectors"))

install.packages("remotes")
remotes::install_github("LangeLab/pepVet", dependencies = TRUE)
```

## The scored workflow

```
digest_protein()  →  score_peptides()  →  verdict
                                               │
evaluate_digest() ─────────────────────────────┤
compare_digests() ─── rank enzymes ─────────────┤
recommend_enzyme()─── top pick ─────────────────┘
batch_evaluate()  ─── whole proteome ──────────►  digest_report()
```

Every function accepts the same flexible input: a raw character sequence,
a named character vector, a `Biostrings::AAString` or `AAStringSet`, or a
FASTA file path.

## Scoring components

| Score        | What it measures                                       | Why it matters                                   |
| ------------ | ------------------------------------------------------ | ------------------------------------------------ |
| `S_length`   | Fraction of peptides in the 7–25 aa detection window   | Short peptides vanish; long ones don't fly       |
| `S_coverage` | Fraction of protein sequence covered by valid peptides | Blind spots hurt quantification confidence       |
| `S_count`    | Normalised peptide count                               | Too few peptides = fragile identification        |
| `S_hydro`    | Mean GRAVY score of valid peptides                     | Extreme hydrophobics stick to columns            |
| `S_charge`   | Fraction of peptides with ≥ 1 basic residue            | Charge enables ESI ionisation                    |
| `S_unique`   | Fraction unique within a proteome _(optional)_         | Shared peptides confound protein-level inference |

The **composite score** is a weighted sum (default weights 0.25/0.25/0.20/0.15/0.15).
Verdicts: **Good** ≥ 0.70 · **Moderate** ≥ 0.40 · **Poor** < 0.40.

## Reference fixtures

The package ships pinned FASTA fixtures in `inst/extdata/` covering a range of
challenging proteins:

| File                               | Protein              | Why it's interesting                         |
| ---------------------------------- | -------------------- | -------------------------------------------- |
| `P02769.fasta`                     | BSA                  | Workhorse standard; trypsin scores Good      |
| `P68431.fasta`                     | Histone H3.1         | Very basic tail; trypsin scores Poor         |
| `Q8WZ42.fasta`                     | Titin                | Largest human protein; stress-tests coverage |
| `P0CG48.fasta`                     | Ubiquitin            | Very small; tests edge-case peptide count    |
| `P37840_isoforms.fasta`            | α-Synuclein          | Multi-isoform; tests proteome-aware scoring  |
| `small_proteome_50_proteins.fasta` | 50-protein human set | Batch workflow fixture                       |

## Learn more

- **Getting Started** — full walkthrough from raw sequence to console report
- **Choosing an Enzyme** — the biology behind each scoring component and a
  worked multi-enzyme comparison
- **Scoring Deep Dive** — weight arithmetic, boundary conditions, and how to
  customise the score for your experiment
- **Reference** — complete function documentation

## License

MIT — see [LICENSE.md](LICENSE.md).

---
