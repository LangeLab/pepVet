#' Amino acid properties
#'
#' Curates a reference table for the 22 amino acids (20 standard + U + O) used
#' by pepVet scoring and validation utilities. The `molecular_weight` column
#' stores free amino acid monoisotopic masses rather than peptide residue
#' masses. Subtract the mass of water (18.01056 Da) to get residue-level
#' masses. The `pKa_side_chain` column stores conventional reference values for
#' the ionizable side chains C, D, E, H, K, R, and Y. Non-ionizable residues
#' are recorded as `NA`.
#'
#' @format A tibble with 22 rows and 6 variables:
#' \describe{
#'   \item{amino_acid}{Single-letter amino acid code.}
#'   \item{molecular_weight}{Free amino acid monoisotopic mass in daltons.}
#'   \item{residue_monoisotopic_mass}{Residue monoisotopic mass in daltons,
#'   equal to \code{molecular_weight - 18.01056}.}
#'   \item{hydrophobicity}{Kyte-Doolittle hydrophobicity value.}
#'   \item{pKa_side_chain}{Conventional side-chain reference pKa, or
#'   \code{NA} for non-ionizable residues.}
#'   \item{is_basic}{Logical flag for basic residues H, K, and R.}
#' }
#'
#' @source
#' Hydrophobicity values follow Kyte J, Doolittle RF (1982). "A simple method
#' for displaying the hydropathic character of a protein." *Journal of
#' Molecular Biology*, 157(1), 105-132.
#' [doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0).
#'
#' Monoisotopic masses follow the free amino acid convention used by ExPASy
#' FindMod: <https://web.expasy.org/findmod/findmod_masses.html>.
#'
#' Side-chain pKa values follow conventional reference values summarized by
#' Thurlkill RL, Grimsley GR, Scholtz JM, Pace CN (2006). "pK values of the
#' ionizable groups of proteins." *Protein Science*, 15(5), 1214-1218.
#' [doi:10.1110/ps.051840806](https://doi.org/10.1110/ps.051840806); and
#' Pace CN, Grimsley GR, Scholtz JM (2009). "Protein ionizable groups: pK
#' values and their contribution to protein stability and solubility."
#' *Journal of Biological Chemistry*, 284(20), 13285-13289.
#' [doi:10.1074/jbc.R800080200](https://doi.org/10.1074/jbc.R800080200).
#'
#' @keywords datasets
aa_properties <- NULL
