aa_properties <- tibble::tibble(
  amino_acid = c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  ),
  molecular_weight = c(
    89.04768, 121.01975, 133.03751, 147.05316, 165.07898,
    75.03203, 155.06948, 131.09463, 146.10553, 131.09463,
    149.05105, 132.05349, 115.06333, 146.06914, 174.11168,
    105.04259, 119.05824, 117.07898, 204.08988, 181.07389
  ),
  hydrophobicity = c(
    1.8, 2.5, -3.5, -3.5, 2.8,
    -0.4, -3.2, 4.5, -3.9, 3.8,
    1.9, -3.5, -1.6, -3.5, -4.5,
    -0.8, -0.7, 4.2, -0.9, -1.3
  ),
  pKa_side_chain = c(
    NA_real_, 8.3, 3.9, 4.3, NA_real_,
    NA_real_, 6.0, NA_real_, 10.5, NA_real_,
    NA_real_, NA_real_, NA_real_, NA_real_, 12.5,
    NA_real_, NA_real_, NA_real_, NA_real_, 10.1
  ),
  is_basic = c(
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, TRUE, FALSE, TRUE, FALSE,
    FALSE, FALSE, FALSE, FALSE, TRUE,
    FALSE, FALSE, FALSE, FALSE, FALSE
  )
)

if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

save(aa_properties, file = "data/aa_properties.rda", compress = "xz")
