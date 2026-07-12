# README figure maintenance

`generate-readme-plots.R` creates the two tracked plots used in the main README:

- `digest_profile_bsa_trypsin.png`
- `batch_comparison_10_enzymes_50_proteins.png`

## Regenerating the plots

Run the script from the repository root when changes to plotting, themes, scoring, or input fixtures affect the README examples. It uses the current package source together with the shipped P02769 and 50-protein FASTA fixtures. The maintenance environment requires `devtools`, `ggplot2`, `patchwork`, and `ragg`.

Set `PEPVET_FIGURE_DIR` to inspect newly generated plots without replacing the tracked files:

```sh
PEPVET_FIGURE_DIR=/tmp/pepvet-readme-figures \
  Rscript man/figures/generate-readme-plots.R
```

Check the image dimensions and inspect the visible content before replacing either tracked file. Byte-for-byte PNG identity is not expected across R versions, graphics devices, fonts, or dependency versions.

## Logo files

`logo.png` and `pepVet-logo.png` are maintained separately from the README plots. They currently contain identical image data because both filenames are still referenced. If those references are consolidated, remove the unused alias in the same change. The README plot script must not replace either logo.
