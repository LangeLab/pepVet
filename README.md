# pepVet

`pepVet` is a Bioconductor-oriented R package for simulating proteolytic
digests and evaluating peptide sets for proteomics workflows.

## Status

This repository currently contains a lean implementation scaffold:

- package metadata in `DESCRIPTION`
- build exclusions in `.Rbuildignore`
- container bootstrap in `Dockerfile` pinned to Bioconductor `RELEASE_3_22`
- dependency isolation via `renv`
- project environment pinned to R `4.5.x` and Bioconductor `3.22`
- a minimal exported API stub in `R/digest.R`
- one initial test file in `tests/testthat/`
- one starter vignette in `vignettes/`

## Next steps

1. Implement the digestion and scoring pipeline in vertical slices.
2. Add runnable CI workflows.
3. Add package documentation, checks, and site automation.
