## data-raw/fetch-external-data.R
## Parse raw tool outputs from MS-Digest and ExPASy PeptideMass into clean CSVs.

library(tibble)

out_dir <- file.path("inst", "extdata", "tool-data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##
## 1. MS-Digest
##
raw_msd <- file.path("inst", "extdata", "tool-data", "msdigest-bsa-trypsin.csv")
if (file.exists(raw_msd)) {
  msd <- read.delim(raw_msd,
    skip = 1, header = FALSE, sep = "\t",
    quote = "", strip.white = TRUE
  )
  names(msd) <- c(
    "num", "mz_mi", "mz_av", "modifications", "start", "end",
    "missed_cleavages", "sequence_raw"
  )

  seq_pattern <- "^\\([A-Z-]+\\)([A-Z]+)\\([A-Z-]+\\)$"
  msd$peptide <- gsub(seq_pattern, "\\1", msd$sequence_raw)
  msd$length <- nchar(msd$peptide)
  msd$tool <- "MS-Digest"
  msd$modifications <- ifelse(nchar(msd$modifications) == 0, NA, msd$modifications)

  msd_unique <- msd[!duplicated(msd$peptide), ]
  msd_unique <- msd_unique[order(msd_unique$start), ]

  write.csv(msd_unique, file.path(out_dir, "tool-compare-msdigest.csv"),
    row.names = FALSE
  )
  cat(sprintf(
    "MS-Digest: %d rows (%d unique peptides)\n", nrow(msd),
    nrow(msd_unique)
  ))
} else {
  cat("MS-Digest CSV not found at", raw_msd, "\n")
}

##
## 2. ExPASy PeptideMass
##
raw_expasy <- file.path("inst", "extdata", "tool-data", "expasy-bsa-trypsin.html")
if (file.exists(raw_expasy)) {
  html <- readLines(raw_expasy, warn = FALSE)

  pat <- "<!-- ([0-9.]+)\\|([0-9]+)-([0-9]+)\\|([0-9]+)\\|([^|]*)\\|*\\|*\\|*\\|*\\|*([A-Z ]+)?.*-->"
  lines <- grep("<!-- [0-9.]+\\|[0-9]+-[0-9]+\\|[0-9]+\\|", html, value = TRUE)

  parse_one <- function(line) {
    m <- regmatches(line, regexec(
      "<!-- ([0-9.]+)\\|([0-9]+)-([0-9]+)\\|([0-9]+)\\|([^|]*)\\|{3,}([A-Z ]+)?.*-->",
      line
    ))[[1]]
    if (length(m) < 7) {
      return(NULL)
    }
    data.frame(
      mass = as.numeric(m[2]),
      start = as.integer(m[3]),
      end = as.integer(m[4]),
      missed_cleavages = as.integer(m[5]),
      modifications = ifelse(nchar(m[6]) == 0, NA, m[6]),
      peptide = gsub(" ", "", m[7]),
      stringsAsFactors = FALSE
    )
  }

  expasy_list <- lapply(lines, parse_one)
  expasy_list <- expasy_list[!vapply(expasy_list, is.null, logical(1))]
  expasy <- do.call(rbind, expasy_list)
  expasy$length <- nchar(expasy$peptide)
  expasy$tool <- "PeptideMass"

  write.csv(expasy, file.path(out_dir, "tool-compare-expasy.csv"),
    row.names = FALSE, na = ""
  )
  cat(sprintf("ExPASy: %d peptides\n", nrow(expasy)))
} else {
  cat("ExPASy HTML not found at", raw_expasy, "\n")
}

##
## 3. pepVet reference data
##
library(pepVet)
bsa_path <- system.file("extdata", "P02769.fasta", package = "pepVet")
bsa <- evaluate_digest(bsa_path, enzyme = "trypsin", missed_cleavages = 1L)
pepvet <- bsa$peptides
names(pepvet)[2] <- "peptide"
pepvet$tool <- "pepVet"
write.csv(pepvet, file.path(out_dir, "pepvet-bsa-trypsin.csv"),
  row.names = FALSE
)
cat(sprintf("pepVet: %d peptides\n", nrow(pepvet)))

cat("\nDone. Files written to", out_dir, "/\n")
