aa_properties <- tibble::tibble(
  amino_acid = c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    # 21st genetically encoded amino acid (UGA codon + SECIS element).
    # Present in ~25 human selenoproteins (GPX1, GPX2, SELENOP, etc.).
    # Hydrophobicity treated as equivalent to Cys (same scale position).
    # pKa of selenol group ~5.2 (vs Cys thiol ~8.3); is_basic = FALSE.
    "U",
    # 22nd genetically encoded amino acid (UAG amber codon).
    # Found in methanogenic archaea (Methanosarcina spp.) and
    # Desulfitobacterium hafniense; not present in humans.
    # Hydrophobicity: no Kyte-Doolittle value exists (scale published 1982;
    # pyrrolysine discovered 2002). NA used, do not substitute a guess.
    # pKa: epsilon-amino is in an amide bond; not titratable. NA.
    # is_basic: epsilon-amino tied up in amide bond. FALSE.
    "O"
  ),
  molecular_weight = c(
    89.04768, 121.01975, 133.03751, 147.05316, 165.07898,
    75.03203, 155.06948, 131.09463, 146.10553, 131.09463,
    149.05105, 132.05349, 115.06333, 146.06914, 174.11168,
    105.04259, 119.05824, 117.07898, 204.08988, 181.07389,
    168.96420,  # C3H7NO2Se monoisotopic (80Se = 79.9165196; PubChem/Unimod)
    255.15829   # C12H21N3O3 monoisotopic; PubChem exact mass 255.15829154
  ),
  hydrophobicity = c(
    1.8, 2.5, -3.5, -3.5, 2.8,
    -0.4, -3.2, 4.5, -3.9, 3.8,
    1.9, -3.5, -1.6, -3.5, -4.5,
    -0.8, -0.7, 4.2, -0.9, -1.3,
    2.5,  # treated as Cys equivalent (Kyte-Doolittle scale)
    NA_real_  # no validated Kyte-Doolittle value; scale predates discovery
  ),
  pKa_side_chain = c(
    NA_real_, 8.3, 3.9, 4.3, NA_real_,
    NA_real_, 6.0, NA_real_, 10.5, NA_real_,
    NA_real_, NA_real_, NA_real_, NA_real_, 12.5,
    NA_real_, NA_real_, NA_real_, NA_real_, 10.1,
    5.2,  # selenol pKa
    NA_real_  # epsilon-amino in amide bond; not titratable
  ),
  is_basic = c(
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, TRUE, FALSE, TRUE, FALSE,
    FALSE, FALSE, FALSE, FALSE, TRUE,
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE,
    FALSE
  )
)

aa_properties$residue_monoisotopic_mass <-
  aa_properties$molecular_weight - 18.01056

aa_properties <- aa_properties[
  c(
    "amino_acid",
    "molecular_weight",
    "residue_monoisotopic_mass",
    "hydrophobicity",
    "pKa_side_chain",
    "is_basic"
  )
]

if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

save(aa_properties, file = "data/aa_properties.rda", compress = "xz")
