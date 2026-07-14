artifact_read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

artifact_expect_schema <- function(data, columns) {
  testthat::expect_identical(names(data), columns)
  testthat::expect_gt(nrow(data), 0L)
}

artifact_expect_finite_or_na <- function(values) {
  testthat::expect_true(all(is.na(values) | is.finite(values)))
}

artifact_expect_range <- function(values, lower, upper, tolerance = 0) {
  observed <- values[!is.na(values)]
  testthat::expect_true(all(is.finite(observed)))
  testthat::expect_true(all(observed >= lower - tolerance))
  testthat::expect_true(all(observed <= upper + tolerance))
}

artifact_mean_or_na <- function(values) {
  if (length(values) == 0L) NA_real_ else mean(values)
}

artifact_median_or_na <- function(values) {
  if (length(values) == 0L) NA_real_ else median(values)
}

artifact_expect_absolute <- function(actual, expected, tolerance, info) {
  if (is.na(expected)) {
    testthat::expect_true(is.na(actual), info = info)
  } else {
    testthat::expect_true(abs(actual - expected) <= tolerance, info = info)
  }
}

artifact_key <- function(data, columns) {
  do.call(paste, c(data[columns], sep = "\r"))
}

artifact_expect_csv_rds_parity <- function(csv_path, rds_path) {
  csv_data <- artifact_read_csv(csv_path)
  rds_data <- as.data.frame(readRDS(rds_path))

  testthat::expect_identical(names(csv_data), names(rds_data))
  testthat::expect_identical(nrow(csv_data), nrow(rds_data))

  for (column in names(csv_data)) {
    if (is.numeric(csv_data[[column]])) {
      testthat::expect_equal(
        csv_data[[column]],
        rds_data[[column]],
        tolerance = 1e-12,
        info = column
      )
    } else {
      testthat::expect_identical(
        as.character(csv_data[[column]]),
        as.character(rds_data[[column]]),
        info = column
      )
    }
  }
}

test_that("tool-data artifacts preserve raw and derived contracts", {
  tool_dir <- system.file("extdata", "tool-data", package = "pepVet")
  expected_files <- c(
    "expasy-bsa-trypsin.html",
    "msdigest-bsa-trypsin.csv",
    "peptideranger-all-peptides.csv",
    "pepvet-bsa-trypsin.csv",
    "tool-compare-expasy.csv",
    "tool-compare-msdigest.csv",
    "tool-compare-peptideranger.csv"
  )

  expect_identical(sort(list.files(tool_dir)), sort(expected_files))

  raw_msd_path <- file.path(tool_dir, "msdigest-bsa-trypsin.csv")
  raw_msd_lines <- readLines(raw_msd_path, warn = FALSE)
  expect_identical(
    raw_msd_lines[[1L]],
    paste(
      c(
        "Number", "m/z (mi)", "m/z (av)", "Modifications", "Start",
        "End", "Missed Cleavages", "Sequence"
      ),
      collapse = "\t"
    )
  )

  raw_msd <- read.delim(
    raw_msd_path,
    skip = 1L,
    header = FALSE,
    sep = "\t",
    quote = "",
    strip.white = TRUE
  )
  names(raw_msd) <- c(
    "num", "mz_mi", "mz_av", "modifications", "start", "end",
    "missed_cleavages", "sequence_raw"
  )
  expect_identical(nrow(raw_msd), 153L)
  raw_pattern <- "^\\([A-Z-]+\\)([A-Z]+)\\([A-Z-]+\\)$"
  expect_true(all(grepl(raw_pattern, raw_msd$sequence_raw)))
  expect_true(all(is.finite(raw_msd$mz_mi)))
  expect_true(all(is.finite(raw_msd$mz_av)))
  expect_true(all(raw_msd$start > 0L & raw_msd$end >= raw_msd$start))
  expect_true(all(raw_msd$missed_cleavages >= 0L))

  pepvet <- artifact_read_csv(file.path(tool_dir, "pepvet-bsa-trypsin.csv"))
  artifact_expect_schema(
    pepvet,
    c(
      "protein_id", "peptide", "start", "end", "length",
      "missed_cleavages", "tool"
    )
  )
  expect_identical(nrow(pepvet), 157L)
  expect_true(all(pepvet$tool == "pepVet"))
  expect_identical(length(unique(pepvet$protein_id)), 1L)
  expect_true(all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", pepvet$peptide)))
  expect_true(all(pepvet$start > 0L & pepvet$end >= pepvet$start))
  expect_true(all(pepvet$end - pepvet$start + 1L == pepvet$length))
  expect_true(all(nchar(pepvet$peptide) == pepvet$length))
  expect_true(all(pepvet$length > 0L & pepvet$missed_cleavages >= 0L))

  msdigest <- artifact_read_csv(
    file.path(tool_dir, "tool-compare-msdigest.csv")
  )
  artifact_expect_schema(
    msdigest,
    c(
      "num", "mz_mi", "mz_av", "modifications", "start", "end",
      "missed_cleavages", "sequence_raw", "peptide", "length", "tool"
    )
  )
  expect_identical(nrow(msdigest), 130L)
  expect_identical(anyDuplicated(msdigest$peptide), 0L)
  expect_true(all(msdigest$tool == "MS-Digest"))
  expect_true(all(is.finite(msdigest$mz_mi) & msdigest$mz_mi > 0))
  expect_true(all(is.finite(msdigest$mz_av) & msdigest$mz_av > 0))
  expect_true(all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", msdigest$peptide)))
  raw_peptides <- sub(raw_pattern, "\\1", raw_msd$sequence_raw)
  expect_identical(sort(unique(raw_peptides)), sort(msdigest$peptide))
  expect_true(all(msdigest$end - msdigest$start + 1L == msdigest$length))
  expect_true(all(nchar(msdigest$peptide) == msdigest$length))
  expect_true(all(msdigest$missed_cleavages >= 0L))

  expasy <- artifact_read_csv(file.path(tool_dir, "tool-compare-expasy.csv"))
  artifact_expect_schema(
    expasy,
    c(
      "mass", "start", "end", "missed_cleavages", "modifications",
      "peptide", "length", "tool"
    )
  )
  expect_identical(nrow(expasy), 155L)
  expect_true(all(expasy$tool == "PeptideMass"))
  expect_true(all(is.finite(expasy$mass) & expasy$mass > 0))
  expect_true(all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", expasy$peptide)))
  expect_true(all(expasy$end - expasy$start + 1L == expasy$length))
  expect_true(all(nchar(expasy$peptide) == expasy$length))
  expect_true(all(expasy$missed_cleavages >= 0L))

  html <- readLines(file.path(tool_dir, "expasy-bsa-trypsin.html"), warn = FALSE)
  html_records <- grep(
    "^<!-- [0-9.]+\\|[0-9]+-[0-9]+\\|[0-9]+\\|",
    html,
    value = TRUE
  )
  expect_identical(length(html_records), nrow(expasy))
  html_matches <- regmatches(
    html_records,
    regexec(
      "^<!-- ([0-9.]+)\\|([0-9]+)-([0-9]+)\\|([0-9]+)\\|",
      html_records
    )
  )
  expect_true(all(vapply(html_matches, length, integer(1)) == 5L))

  pr_all <- artifact_read_csv(
    file.path(tool_dir, "peptideranger-all-peptides.csv")
  )
  artifact_expect_schema(
    pr_all,
    c(
      "protein", "enzyme", "peptide", "length", "gravy", "window_pass",
      "pr_score", "composite_score", "verdict"
    )
  )
  expect_true(all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", pr_all$peptide)))
  expect_true(all(nchar(pr_all$peptide) == pr_all$length))
  expect_true(all(is.finite(pr_all$gravy)))
  artifact_expect_range(pr_all$pr_score, 0, 1, tolerance = 1e-12)
  artifact_expect_range(pr_all$composite_score, 0, 1)
  expect_true(all(pr_all$window_pass %in% c(TRUE, FALSE)))
  expect_true(all(pr_all$verdict %in% c("Good", "Moderate", "Poor")))

  pr_summary <- artifact_read_csv(
    file.path(tool_dir, "tool-compare-peptideranger.csv")
  )
  artifact_expect_schema(
    pr_summary,
    c(
      "protein", "enzyme", "n_total", "n_window_pass", "composite",
      "verdict", "mean_PR_all", "mean_PR_window_pass",
      "median_PR_window_pass", "cor_PR_length", "cor_PR_gravy",
      "mean_score_difference"
    )
  )
  expect_identical(nrow(pr_summary), 9L)
  pr_all_key <- artifact_key(pr_all, c("protein", "enzyme"))
  pr_summary_key <- artifact_key(pr_summary, c("protein", "enzyme"))
  expect_equal(sort(unique(pr_all_key)), sort(pr_summary_key))
  expect_identical(anyDuplicated(pr_summary_key), 0L)
  expect_equal(
    as.integer(table(pr_all_key)[pr_summary_key]),
    pr_summary$n_total
  )
  pass_counts <- tapply(pr_all$window_pass, pr_all_key, sum)
  expect_equal(
    as.integer(pass_counts[pr_summary_key]),
    pr_summary$n_window_pass
  )
  expect_true(all(pr_summary$n_total > 0L))
  expect_true(all(pr_summary$n_window_pass >= 0L &
    pr_summary$n_window_pass <= pr_summary$n_total))
  artifact_expect_range(pr_summary$composite, 0, 1)
  expect_true(all(is.finite(pr_summary$mean_PR_all)))
  for (row_index in seq_len(nrow(pr_summary))) {
    group <- pr_all[pr_all_key == pr_summary_key[[row_index]], , drop = FALSE]
    pass_scores <- group$pr_score[group$window_pass]
    fail_scores <- group$pr_score[!group$window_pass]
    expect_equal(
      pr_summary$mean_PR_all[[row_index]],
      mean(group$pr_score),
      tolerance = 1e-12,
      info = pr_summary_key[[row_index]]
    )
    expect_equal(
      pr_summary$mean_PR_window_pass[[row_index]],
      artifact_mean_or_na(pass_scores),
      tolerance = 1e-12,
      info = pr_summary_key[[row_index]]
    )
    expect_equal(
      pr_summary$median_PR_window_pass[[row_index]],
      artifact_median_or_na(pass_scores),
      tolerance = 1e-12,
      info = pr_summary_key[[row_index]]
    )
    expect_equal(
      pr_summary$mean_score_difference[[row_index]],
      artifact_mean_or_na(pass_scores) - artifact_mean_or_na(fail_scores),
      tolerance = 1e-12,
      info = pr_summary_key[[row_index]]
    )
    artifact_expect_absolute(
      pr_summary$cor_PR_length[[row_index]],
      cor(group$pr_score, group$length),
      tolerance = 5.1e-5,
      info = pr_summary_key[[row_index]]
    )
    artifact_expect_absolute(
      pr_summary$cor_PR_gravy[[row_index]],
      cor(group$pr_score, group$gravy),
      tolerance = 5.1e-5,
      info = pr_summary_key[[row_index]]
    )
  }
  for (column in c(
    "mean_PR_all", "mean_PR_window_pass", "median_PR_window_pass",
    "cor_PR_length", "cor_PR_gravy", "mean_score_difference"
  )) {
    artifact_expect_finite_or_na(pr_summary[[column]])
  }
})

test_that("comparison-data artifacts preserve grids, categories, and overlap math", {
  comparison_dir <- system.file(
    "extdata", "comparison-data", package = "pepVet"
  )
  schemas <- list(
    "sectionA-overlap.csv" = c(
      "Tool", "N_unique", "Common_with_pepVet", "Shared_all_three",
      "Pct_overlap"
    ),
    "sectionB-capabilities.csv" = c(
      "Capability", "pepVet", "MS_Digest", "ExPASy", "Protein_Cleaver",
      "ProteaseGuru", "PeptideRanger"
    ),
    "sectionC-scores-all.csv" = c(
      "protein", "enzyme", "n_total", "n_valid", "S_length",
      "S_coverage", "S_count", "S_hydro", "S_charge", "composite",
      "verdict"
    ),
    "sectionC2-presets.csv" = c(
      "protein", "preset", "composite", "verdict", "n_total",
      "n_length_valid", "S_hydro", "S_unique", "preset_used",
      "include_pI", "len_lo", "len_hi"
    ),
    "sectionD-peptideranger.csv" = c(
      "protein", "enzyme", "verdict", "n_total", "n_window_pass",
      "PR_mean_all", "PR_mean_window_pass", "PR_mean_window_fail",
      "mean_score_difference", "cor_PR_length"
    ),
    "sectionE-protein-cleaver.csv" = c(
      "protein", "n_peptides", "pc_identifiable", "pepvet_window_pass",
      "both_pass", "pc_only", "pepvet_window_only"
    )
  )
  expected_rows <- c(
    "sectionA-overlap.csv" = 3L,
    "sectionB-capabilities.csv" = 16L,
    "sectionC-scores-all.csv" = 25L,
    "sectionC2-presets.csv" = 18L,
    "sectionD-peptideranger.csv" = 25L,
    "sectionE-protein-cleaver.csv" = 3L
  )
  expect_identical(sort(list.files(comparison_dir)), sort(names(schemas)))

  tables <- lapply(names(schemas), function(name) {
    data <- artifact_read_csv(file.path(comparison_dir, name))
    artifact_expect_schema(data, schemas[[name]])
    expect_identical(nrow(data), expected_rows[[name]])
    data
  })
  names(tables) <- names(schemas)

  pepvet <- artifact_read_csv(
    system.file("extdata", "tool-data", "pepvet-bsa-trypsin.csv", package = "pepVet")
  )
  msdigest <- artifact_read_csv(
    system.file("extdata", "tool-data", "tool-compare-msdigest.csv", package = "pepVet")
  )
  expasy <- artifact_read_csv(
    system.file("extdata", "tool-data", "tool-compare-expasy.csv", package = "pepVet")
  )
  peptide_sets <- list(
    "MS-Digest" = unique(msdigest$peptide),
    "ExPASy PeptideMass" = unique(expasy$peptide),
    pepVet = unique(pepvet$peptide)
  )
  pepvet_set <- peptide_sets$pepVet
  expected_overlap <- data.frame(
    Tool = unname(names(peptide_sets)),
    N_unique = unname(vapply(peptide_sets, length, integer(1))),
    Common_with_pepVet = unname(vapply(
      peptide_sets,
      function(peptides) sum(peptides %in% pepvet_set),
      integer(1)
    )),
    Shared_all_three = unname(vapply(
      peptide_sets,
      function(peptides) {
        sum(peptides %in% peptide_sets[["MS-Digest"]] &
          peptides %in% peptide_sets[["ExPASy PeptideMass"]] &
          peptides %in% pepvet_set)
      },
      integer(1)
    )),
    stringsAsFactors = FALSE
  )
  expected_overlap$Pct_overlap <- round(
    expected_overlap$Common_with_pepVet / expected_overlap$N_unique * 100,
    1
  )
  expect_equal(tables[["sectionA-overlap.csv"]], expected_overlap)

  expected_capabilities <- c(
    "Peptide list + masses", "pI output", "Hydrophobicity metric",
    "Composite digest score", "Score-band label",
    "Multi-enzyme comparison", "Enzyme ranking", "Workflow presets",
    "Batch input", "Sensitivity analysis", "Peptide-level ML score",
    "3D structure mapping",
    "Retention time prediction", "Skyline/FASTA export",
    "Non-GUI interface", "R package"
  )
  capabilities <- tables[["sectionB-capabilities.csv"]]
  expect_identical(capabilities$Capability, expected_capabilities)
  capability_values <- unlist(capabilities[-1L], use.names = FALSE)
  expect_false(anyNA(capability_values))
  expect_true(all(nzchar(capability_values)))

  scores <- tables[["sectionC-scores-all.csv"]]
  expected_score_grid <- expand.grid(
    protein = c("BSA", "H3", "BACE1", "LYSO", "UBIQ"),
    enzyme = c(
      "trypsin", "lysc", "chymotrypsin-high", "glutamyl endopeptidase",
      "asp-n endopeptidase"
    ),
    stringsAsFactors = FALSE
  )
  expect_equal(
    sort(artifact_key(scores, c("protein", "enzyme"))),
    sort(artifact_key(expected_score_grid, c("protein", "enzyme")))
  )
  expect_identical(
    anyDuplicated(artifact_key(scores, c("protein", "enzyme"))),
    0L
  )
  score_columns <- c(
    "S_length", "S_coverage", "S_count", "S_hydro", "S_charge", "composite"
  )
  for (column in score_columns) {
    artifact_expect_range(scores[[column]], 0, 1)
  }
  expect_true(all(scores$n_total > 0L))
  expect_true(all(scores$n_valid >= 0L & scores$n_valid <= scores$n_total))
  verdict_good <- get(".pepvet_params", asNamespace("pepVet"))$verdict_good
  verdict_moderate <- get(".pepvet_params", asNamespace("pepVet"))$verdict_moderate
  expected_verdict <- ifelse(
    scores$composite >= verdict_good,
    "Good",
    ifelse(scores$composite >= verdict_moderate, "Moderate", "Poor")
  )
  expect_identical(scores$verdict, expected_verdict)

  comparison_fastas <- c(
    BSA = "P02769.fasta",
    H3 = "P68431.fasta",
    BACE1 = "P56817.fasta",
    LYSO = "P00698.fasta",
    UBIQ = "P0CG48.fasta"
  )
  comparison_enzymes <- c(
    "trypsin", "lysc", "chymotrypsin-high",
    "glutamyl endopeptidase", "asp-n endopeptidase"
  )
  expected_scores <- do.call(rbind, lapply(names(comparison_fastas), function(protein) {
    do.call(rbind, lapply(comparison_enzymes, function(enzyme) {
      result <- evaluate_digest(
        system.file(
          "extdata", comparison_fastas[[protein]], package = "pepVet"
        ),
        enzyme = enzyme,
        missed_cleavages = 1L
      )
      score <- result$scores
      data.frame(
        protein = protein,
        enzyme = enzyme,
        n_total = nrow(result$peptides),
        n_valid = sum(result$peptides$length >= 7L &
          result$peptides$length <= 25L),
        S_length = score$S_length,
        S_coverage = score$S_coverage,
        S_count = score$S_count,
        S_hydro = score$S_hydro,
        S_charge = score$S_charge,
        composite = score$composite_score,
        verdict = score$verdict,
        stringsAsFactors = FALSE
      )
    }))
  }))
  expected_score_key <- artifact_key(expected_scores, c("protein", "enzyme"))
  observed_score_key <- artifact_key(scores, c("protein", "enzyme"))
  score_index <- match(expected_score_key, observed_score_key)
  expect_false(anyNA(score_index))
  expect_identical(anyDuplicated(observed_score_key), 0L)
  for (column in c(
    "n_total", "n_valid", "S_length", "S_coverage", "S_count",
    "S_hydro", "S_charge", "composite"
  )) {
    expect_equal(
      scores[[column]][score_index],
      expected_scores[[column]],
      tolerance = 1e-12,
      info = column
    )
  }
  expect_identical(scores$verdict[score_index], expected_scores$verdict)

  presets <- tables[["sectionC2-presets.csv"]]
  preset_registry <- get(".pepvet_presets", asNamespace("pepVet"))
  expect_identical(sort(unique(presets$preset)), sort(names(preset_registry)))
  for (preset_name in names(preset_registry)) {
    rows <- presets[presets$preset == preset_name, , drop = FALSE]
    preset <- preset_registry[[preset_name]]
    expect_identical(nrow(rows), 3L, info = preset_name)
    expect_true(all(rows$len_lo == preset$length_range[[1L]]))
    expect_true(all(rows$len_hi == preset$length_range[[2L]]))
    expect_true(all(rows$n_total > 0L & rows$n_length_valid >= 0L))
    expect_true(all(rows$n_length_valid <= rows$n_total))
    artifact_expect_range(rows$composite, 0, 1)
    artifact_expect_range(rows$S_hydro, 0, 1)
    artifact_expect_range(rows$S_unique, 0, 1)
    expect_true(all(rows$preset_used == preset_name))
    expect_true(all(rows$include_pI == preset$include_pI))
  }

  optional_rows <- tables[["sectionD-peptideranger.csv"]]
  expect_identical(
    sort(unique(optional_rows$protein)),
    sort(c("BSA", "H3", "BACE1", "LYSO", "UBIQ"))
  )
  expect_identical(
    sort(unique(optional_rows$enzyme)),
    sort(c(
      "trypsin", "lysc", "chymotrypsin-high", "glutamyl endopeptidase",
      "asp-n endopeptidase"
    ))
  )
  optional_key <- artifact_key(optional_rows, c("protein", "enzyme"))
  expect_identical(anyDuplicated(optional_key), 0L)
  expect_equal(
    sort(optional_key),
    sort(artifact_key(expected_score_grid, c("protein", "enzyme")))
  )
  expect_true(all(optional_rows$n_total > 0L))
  expect_true(all(optional_rows$n_window_pass >= 0L &
    optional_rows$n_window_pass <= optional_rows$n_total))
  artifact_expect_range(optional_rows$PR_mean_all, 0, 1, tolerance = 1e-12)
  artifact_expect_range(
    optional_rows$PR_mean_window_pass, 0, 1, tolerance = 1e-12
  )
  artifact_expect_range(
    optional_rows$PR_mean_window_fail, 0, 1, tolerance = 1e-12
  )
  artifact_expect_range(
    optional_rows$mean_score_difference, -1, 1, tolerance = 1e-12
  )
  expect_true(all(is.finite(optional_rows$PR_mean_all)))
  for (column in c(
    "PR_mean_all", "PR_mean_window_pass", "PR_mean_window_fail",
    "mean_score_difference", "cor_PR_length"
  )) {
    artifact_expect_finite_or_na(optional_rows[[column]])
  }
  expect_true(all(optional_rows$verdict %in% c("Good", "Moderate", "Poor")))

  protein_cleaver <- tables[["sectionE-protein-cleaver.csv"]]
  expect_identical(protein_cleaver$protein, c("BSA", "H3", "BACE1"))
  count_columns <- c(
    "n_peptides", "pc_identifiable", "pepvet_window_pass", "both_pass",
    "pc_only", "pepvet_window_only"
  )
  for (column in count_columns) {
    expect_true(all(protein_cleaver[[column]] >= 0L))
  }
  expect_true(all(protein_cleaver$pc_identifiable <= protein_cleaver$n_peptides))
  expect_true(all(
    protein_cleaver$pepvet_window_pass <= protein_cleaver$n_peptides
  ))
  expect_true(all(protein_cleaver$both_pass <= protein_cleaver$pc_identifiable))
  expect_true(all(
    protein_cleaver$both_pass <= protein_cleaver$pepvet_window_pass
  ))
  expect_true(all(
    protein_cleaver$pc_identifiable ==
      protein_cleaver$both_pass + protein_cleaver$pc_only
  ))
  expect_true(all(
    protein_cleaver$pepvet_window_pass ==
      protein_cleaver$both_pass + protein_cleaver$pepvet_window_only
  ))
  expect_true(all(
    protein_cleaver$n_peptides >=
      protein_cleaver$both_pass + protein_cleaver$pc_only +
        protein_cleaver$pepvet_window_only
  ))
})

test_that("preset comparison rows use complete preset configurations", {
  comparison_dir <- system.file(
    "extdata", "comparison-data", package = "pepVet"
  )
  observed <- artifact_read_csv(
    file.path(comparison_dir, "sectionC2-presets.csv")
  )
  fasta_names <- c(
    BSA = "P02769.fasta",
    H3 = "P68431.fasta",
    BACE1 = "P56817.fasta"
  )
  fasta_paths <- vapply(
    fasta_names,
    function(name) system.file("extdata", name, package = "pepVet"),
    character(1)
  )
  background_sequences <- lapply(
    fasta_paths,
    Biostrings::readAAStringSet
  )
  background <- digest_protein(
    do.call(c, unname(background_sequences)),
    enzyme = "trypsin",
    missed_cleavages = 1L
  )

  for (row_index in seq_len(nrow(observed))) {
    row <- observed[row_index, , drop = FALSE]
    preset <- pepvet_preset(row$preset)
    direct <- do.call(
      evaluate_digest,
      c(
        list(
          sequence = fasta_paths[[row$protein]],
          enzyme = "trypsin",
          missed_cleavages = 1L,
          proteome = background
        ),
        preset
      )
    )
    score <- direct$scores
    info <- paste(row$protein, row$preset, sep = " / ")
    expect_equal(row$composite, score$composite_score, tolerance = 1e-12,
      info = info
    )
    expect_identical(row$verdict, score$verdict, info = info)
    expect_equal(row$S_hydro, score$S_hydro, tolerance = 1e-12, info = info)
    expect_equal(row$S_unique, score$S_unique, tolerance = 1e-12, info = info)
    expect_identical(row$preset_used, score$preset_used, info = info)
    expect_identical(row$include_pI, direct$params$include_pI, info = info)
    expect_equal(row$n_total, nrow(direct$peptides), info = info)
    expect_equal(
      row$n_length_valid,
      sum(direct$peptides$length >= preset$length_range[[1L]] &
        direct$peptides$length <= preset$length_range[[2L]]),
      info = info
    )
  }
})

test_that("PeptideAtlas artifacts preserve parity and scientific invariants", {
  artifact_dir <- system.file(
    "extdata", "peptideatlas-concordance", package = "pepVet"
  )
  expect_identical(
    sort(list.files(artifact_dir)),
    sort(c(
      "grid-search-results.csv", "grid-search-results.rds",
      "per-protein-results.csv", "per-protein-results.rds",
      "threshold-calibration.csv", "threshold-calibration.rds"
    ))
  )
  artifact_expect_csv_rds_parity(
    file.path(artifact_dir, "per-protein-results.csv"),
    file.path(artifact_dir, "per-protein-results.rds")
  )
  artifact_expect_csv_rds_parity(
    file.path(artifact_dir, "grid-search-results.csv"),
    file.path(artifact_dir, "grid-search-results.rds")
  )
  artifact_expect_csv_rds_parity(
    file.path(artifact_dir, "threshold-calibration.csv"),
    file.path(artifact_dir, "threshold-calibration.rds")
  )

  per_protein <- artifact_read_csv(
    file.path(artifact_dir, "per-protein-results.csv")
  )
  artifact_expect_schema(
    per_protein,
    c(
      "protein_id", "n_theoretical", "n_observed", "n_pepVet_valid",
      "n_observed_and_valid", "n_observed_and_invalid", "detection_rate_all",
      "detection_rate_valid", "detection_rate_invalid", "protein_length",
      "composite_score", "verdict"
    )
  )
  expect_identical(nrow(per_protein), 500L)
  expect_identical(anyDuplicated(per_protein$protein_id), 0L)
  expect_type(per_protein$protein_length, "integer")
  expect_true(all(per_protein$protein_length >= 50L))
  count_columns <- c(
    "n_theoretical", "n_observed", "n_pepVet_valid", "n_observed_and_valid",
    "n_observed_and_invalid"
  )
  for (column in count_columns) {
    expect_true(all(!is.na(per_protein[[column]]) & per_protein[[column]] >= 0L))
  }
  expect_true(all(per_protein$n_observed <= per_protein$n_theoretical))
  expect_true(all(per_protein$n_pepVet_valid <= per_protein$n_theoretical))
  expect_true(all(
    per_protein$n_observed ==
      per_protein$n_observed_and_valid + per_protein$n_observed_and_invalid
  ))
  expect_true(all(
    per_protein$n_observed_and_valid <=
      pmin(per_protein$n_observed, per_protein$n_pepVet_valid)
  ))
  expect_true(all(
    per_protein$n_observed_and_invalid <= per_protein$n_observed
  ))
  for (column in c(
    "detection_rate_all", "detection_rate_valid", "detection_rate_invalid"
  )) {
    artifact_expect_finite_or_na(per_protein[[column]])
    artifact_expect_range(per_protein[[column]], 0, 1)
  }
  expect_true(all(is.finite(per_protein$detection_rate_all)))
  artifact_expect_range(per_protein$composite_score, 0, 1)
  expect_true(all(per_protein$verdict %in% c("Good", "Moderate", "Poor")))

  grid <- artifact_read_csv(file.path(artifact_dir, "grid-search-results.csv"))
  artifact_expect_schema(
    grid,
    c(
      "length_min", "length_max", "gravy_min", "gravy_max",
      "n_valid_total", "n_invalid_total", "detection_rate_valid",
      "detection_rate_invalid", "enrichment"
    )
  )
  expect_identical(nrow(grid), 900L)
  expect_identical(
    anyDuplicated(artifact_key(
      grid,
      c("length_min", "length_max", "gravy_min", "gravy_max")
    )),
    0L
  )
  expect_true(all(grid$length_min %in% c(5L, 6L, 7L, 8L, 9L, 10L)))
  expect_true(all(grid$length_max %in% c(20L, 25L, 30L, 35L, 40L)))
  expect_true(all(grid$gravy_min %in% c(-1.5, -1.2, -1.0, -0.8, -0.5)))
  expect_true(all(grid$gravy_max %in% c(0.3, 0.4, 0.5, 0.6, 0.8, 1.0)))
  expect_true(all(grid$length_max > grid$length_min))
  expect_true(all(grid$gravy_max > grid$gravy_min))
  expected_theoretical <- sum(per_protein$n_theoretical)
  expect_true(all(
    grid$n_valid_total + grid$n_invalid_total == expected_theoretical
  ))
  artifact_expect_range(grid$detection_rate_valid, 0, 1)
  artifact_expect_range(grid$detection_rate_invalid, 0, 1)
  artifact_expect_range(grid$enrichment, -1, 1)
  expect_true(all(is.finite(grid$detection_rate_valid)))
  expect_true(all(is.finite(grid$detection_rate_invalid)))
  expect_true(all(is.finite(grid$enrichment)))
  expect_true(all(diff(grid$enrichment) <= 1e-12))
  expect_true(all(
    abs(grid$enrichment -
      (grid$detection_rate_valid - grid$detection_rate_invalid)) < 1e-12
  ))

  threshold <- artifact_read_csv(
    file.path(artifact_dir, "threshold-calibration.csv")
  )
  artifact_expect_schema(
    threshold,
    c("threshold", "sensitivity", "specificity", "youden_j")
  )
  expect_identical(nrow(threshold), 171L)
  expect_equal(threshold$threshold, seq(0.1, 0.95, by = 0.005), tolerance = 1e-12)
  artifact_expect_range(threshold$threshold, 0.1, 0.95)
  artifact_expect_range(threshold$sensitivity, 0, 1)
  artifact_expect_range(threshold$specificity, 0, 1)
  artifact_expect_range(threshold$youden_j, -1, 1)
  expect_true(all(is.finite(threshold$sensitivity)))
  expect_true(all(is.finite(threshold$specificity)))
  expect_true(all(is.finite(threshold$youden_j)))
  expect_true(all(diff(threshold$sensitivity) <= 1e-12))
  expect_true(all(diff(threshold$specificity) >= -1e-12))
  expect_equal(
    threshold$youden_j,
    threshold$sensitivity + threshold$specificity - 1,
    tolerance = 1e-12
  )
})
