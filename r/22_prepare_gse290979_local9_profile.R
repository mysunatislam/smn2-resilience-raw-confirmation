source(file.path("r", "common.R"))

arguments <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% arguments

full_sheet_path <- file.path(
  ROOT,
  "config",
  "GSE290979_raw_sample_sheet.tsv"
)
profile_sheet_path <- file.path(
  ROOT,
  "config",
  "GSE290979_local9_sample_sheet.tsv"
)
selection_path <- file.path(
  ROOT,
  "config",
  "GSE290979_local9_selection.tsv"
)
profile_directory <- file.path(
  ROOT,
  "config",
  "local_raw",
  "GSE290979",
  "local9"
)
summary_path <- file.path(profile_directory, "profile_summary.tsv")
benchmark_path <- file.path(profile_directory, "benchmark_sample.txt")

full <- read.delim(
  full_sheet_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_columns <- c(
  "array_index", "sample_id", "run_accession", "gsm", "donor_line",
  "genotype", "treatment", "replicate", "analysis_group", "read_pairs",
  "base_count", "read_length", "fastq_bytes_total", "fastq_r1_url",
  "fastq_r2_url", "fastq_r1_md5", "fastq_r2_md5", "ena_instrument",
  "paper_instrument"
)
stopifnot(
  identical(names(full), required_columns),
  nrow(full) == 31L,
  identical(full$array_index, seq_len(31L)),
  !anyDuplicated(full$sample_id),
  !anyDuplicated(full$run_accession),
  all(full$read_length == 151L)
)

select_central_library <- function(frame) {
  distance <- abs(frame$read_pairs - stats::median(frame$read_pairs))
  frame[order(distance, frame$run_accession)[1L], , drop = FALSE]
}

untreated <- full[full$treatment == "NT", , drop = FALSE]
selected_untreated <- do.call(
  rbind,
  lapply(split(untreated, untreated$donor_line), select_central_library)
)

treatment <- full[
  full$treatment %in% c("Scramble", "R6-Mo"),
  ,
  drop = FALSE
]
treatment_key <- interaction(
  treatment$donor_line,
  treatment$treatment,
  drop = TRUE,
  lex.order = TRUE
)
selected_treatment <- do.call(
  rbind,
  lapply(split(treatment, treatment_key), select_central_library)
)

selected <- rbind(selected_untreated, selected_treatment)
analysis_order <- c("CTRL_NT", "SMA_NT", "SMA_Scramble", "SMA_R6-Mo")
line_order <- c("C1", "C2", "C3", "S2", "S3")
selected <- selected[order(
  match(selected$analysis_group, analysis_order),
  match(selected$donor_line, line_order),
  selected$replicate,
  selected$run_accession
), , drop = FALSE]
selected$array_index <- seq_len(nrow(selected))
rownames(selected) <- NULL

stopifnot(
  nrow(selected) == 9L,
  identical(
    selected$sample_id,
    c(
      "BULK-SAM-109", "BULK-SAM-104", "BULK-SAM-107",
      "BULK-SAM-160", "BULK-SAM-148", "BULK-SAM-163",
      "BULK-SAM-151", "BULK-SAM-166", "BULK-SAM-157"
    )
  ),
  sum(selected$treatment == "NT") == 5L,
  sum(selected$treatment == "Scramble") == 2L,
  sum(selected$treatment == "R6-Mo") == 2L,
  setequal(
    unique(selected$donor_line[selected$treatment == "NT"]),
    line_order
  ),
  setequal(
    unique(selected$donor_line[selected$treatment != "NT"]),
    c("S2", "S3")
  )
)

selection <- selected[, c(
  "array_index", "sample_id", "run_accession", "gsm", "donor_line",
  "genotype", "treatment", "replicate", "analysis_group", "read_pairs",
  "fastq_bytes_total"
)]
selection$selection_basis <- ifelse(
  selection$treatment == "NT",
  "closest_to_within_line_median_read_pairs",
  "closest_to_within_line_condition_median_read_pairs"
)

sample_lists <- list(
  disease_group1_sma = selected$sample_id[
    selected$genotype == "SMA" & selected$treatment == "NT"
  ],
  disease_group2_control = selected$sample_id[
    selected$genotype == "CTRL" & selected$treatment == "NT"
  ],
  treatment_group1_r6 = selected$sample_id[
    selected$treatment == "R6-Mo"
  ],
  treatment_group2_scramble = selected$sample_id[
    selected$treatment == "Scramble"
  ]
)
stopifnot(
  identical(unname(lengths(sample_lists)), c(2L, 3L, 2L, 2L)),
  all(vapply(sample_lists, function(values) !anyDuplicated(values), logical(1)))
)

benchmark_row <- selected[
  order(selected$fastq_bytes_total, selected$run_accession)[1L],
  ,
  drop = FALSE
]
profile_summary <- data.frame(
  metric = c(
    "profile", "selection_schema_version", "deposited_libraries",
    "selected_libraries", "selected_untreated_libraries",
    "selected_treatment_libraries", "disease_independent_lines",
    "treatment_independent_lines", "selected_fastq_bytes",
    "selected_fastq_gib", "deposited_fastq_gib", "benchmark_sample_id",
    "benchmark_run_accession", "benchmark_fastq_gib",
    "untreated_selection_rule", "treatment_selection_rule",
    "benchmark_selection_rule", "raw_scope"
  ),
  value = c(
    "local9", "1", nrow(full), nrow(selected),
    sum(selected$treatment == "NT"),
    sum(selected$treatment != "NT"),
    length(unique(selected$donor_line[selected$treatment == "NT"])),
    length(unique(selected$donor_line[selected$treatment != "NT"])),
    format(sum(selected$fastq_bytes_total), scientific = FALSE, trim = TRUE),
    sprintf("%.3f", sum(selected$fastq_bytes_total) / 1024^3),
    sprintf("%.3f", sum(full$fastq_bytes_total) / 1024^3),
    benchmark_row$sample_id,
    benchmark_row$run_accession,
    sprintf("%.3f", benchmark_row$fastq_bytes_total / 1024^3),
    "closest_to_within_line_median_read_pairs_then_run_accession",
    paste0(
      "closest_to_within_line_condition_median_read_pairs_",
      "then_run_accession"
    ),
    "smallest_complete_fastq_pair_inside_frozen_local9",
    paste0(
      "pre_specified_local_raw_read_sensitivity_not_validation21_",
      "and_not_full_reanalysis"
    )
  ),
  stringsAsFactors = FALSE
)

assert_generated_profile <- function() {
  required_files <- c(
    profile_sheet_path,
    selection_path,
    summary_path,
    benchmark_path,
    file.path(profile_directory, paste0(names(sample_lists), ".txt"))
  )
  if (!all(file.exists(required_files))) {
    stop("Generated local9 profile files are missing")
  }
  observed_sheet <- read.delim(
    profile_sheet_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  observed_selection <- read.delim(
    selection_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  observed_summary <- read.delim(
    summary_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(
    identical(observed_sheet, selected),
    identical(observed_selection, selection),
    identical(observed_summary, profile_summary),
    identical(readLines(benchmark_path, warn = FALSE), benchmark_row$sample_id)
  )
  for (name in names(sample_lists)) {
    observed <- readLines(
      file.path(profile_directory, paste0(name, ".txt")),
      warn = FALSE
    )
    stopifnot(identical(observed, sample_lists[[name]]))
  }
  invisible(TRUE)
}

if (check_only) {
  assert_generated_profile()
  cat(
    "local9 profile check passed:",
    nrow(selected), "libraries,",
    sprintf("%.1f GiB", sum(selected$fastq_bytes_total) / 1024^3),
    "compressed FASTQ; benchmark", benchmark_row$sample_id, "\n"
  )
  quit(save = "no", status = 0L)
}

dir.create(profile_directory, recursive = TRUE, showWarnings = FALSE)
write.table(
  selected,
  profile_sheet_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
write.table(
  selection,
  selection_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
for (name in names(sample_lists)) {
  writeLines(
    sample_lists[[name]],
    file.path(profile_directory, paste0(name, ".txt")),
    useBytes = TRUE
  )
}
writeLines(benchmark_row$sample_id, benchmark_path, useBytes = TRUE)
write.table(
  profile_summary,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
assert_generated_profile()

cat(
  "Prepared local9 profile:",
  nrow(selected), "libraries,",
  sprintf("%.1f GiB", sum(selected$fastq_bytes_total) / 1024^3),
  "compressed FASTQ; benchmark", benchmark_row$sample_id, "\n"
)
