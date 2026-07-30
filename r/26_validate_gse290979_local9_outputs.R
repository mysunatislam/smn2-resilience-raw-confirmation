script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) != 1L) {
  stop("This script must be run with Rscript")
}
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

parse_options <- function(arguments) {
  parsed <- list()
  for (argument in arguments) {
    pieces <- regmatches(
      argument,
      regexec("^--([A-Za-z0-9-]+)=(.*)$", argument)
    )[[1L]]
    if (length(pieces) != 3L) {
      stop("Expected --name=value argument, received: ", argument)
    }
    parsed[[pieces[2L]]] <- pieces[3L]
  }
  parsed
}

option_value <- function(options, name, default = NULL) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) default else value
}

read_tsv <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_tsv <- function(frame, path) {
  write.table(
    frame,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

metric_value <- function(frame, metric) {
  value <- frame$value[match(metric, frame$metric)]
  if (length(value) != 1L || is.na(value)) {
    stop("Required summary metric is missing: ", metric)
  }
  as.character(value)
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
output_root <- option_value(
  options,
  "output-root",
  file.path(ROOT, "results", "r", "raw_confirmation", "local9")
)
sample_sheet_path <- option_value(
  options,
  "sample-sheet",
  file.path(ROOT, "config", "GSE290979_local9_sample_sheet.tsv")
)
expected_candidates <- suppressWarnings(as.integer(option_value(
  options,
  "expected-candidates",
  "37"
)))
minimum_matched_genes <- suppressWarnings(as.integer(option_value(
  options,
  "minimum-matched-genes",
  "1000"
)))
completion_marker <- option_value(
  options,
  "completion-marker",
  "LOCAL9_RAW_ANALYSIS_COMPLETE.tsv"
)
expected_quantification_method <- option_value(
  options,
  "expected-quantification-method",
  "genome_alignment_featurecounts"
)
expected_whole_genome_alignment <- toupper(option_value(
  options,
  "expected-whole-genome-alignment",
  "TRUE"
))
if (is.na(expected_candidates) || expected_candidates < 1L) {
  stop("--expected-candidates must be a positive integer")
}
if (is.na(minimum_matched_genes) || minimum_matched_genes < 1L) {
  stop("--minimum-matched-genes must be a positive integer")
}
if (basename(completion_marker) != completion_marker) {
  stop("--completion-marker must be a file name")
}
if (!expected_whole_genome_alignment %in% c("TRUE", "FALSE")) {
  stop("--expected-whole-genome-alignment must be TRUE or FALSE")
}

file_names <- c(
  raw_counts = "GSE290979_local9_raw_gene_counts.tsv",
  analysis_summary = "GSE290979_local9_analysis_summary.tsv",
  disease = "GSE290979_local9_SMA_vs_CTRL_edgeR.tsv",
  disease_lolo = "GSE290979_local9_disease_LOLO_summary.tsv",
  treatment = "GSE290979_local9_R6_vs_scramble_paired_edgeR.tsv",
  treatment_linewise = "GSE290979_local9_treatment_linewise_summary.tsv",
  concordance =
    "GSE290979_local9_raw_vs_processed_concordance_summary.tsv",
  candidates = "GSE290979_local9_raw_candidate_confirmation.tsv",
  candidate_summary =
    "GSE290979_local9_raw_candidate_confirmation_summary.tsv"
)
paths <- setNames(file.path(output_root, file_names), names(file_names))
required_paths <- c(sample_sheet_path, paths)
if (!all(file.exists(required_paths))) {
  stop(
    "Local9 validation inputs are missing: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}

sample_sheet <- read_tsv(sample_sheet_path)
raw_counts <- read_tsv(paths[["raw_counts"]])
analysis_summary <- read_tsv(paths[["analysis_summary"]])
disease <- read_tsv(paths[["disease"]])
disease_lolo <- read_tsv(paths[["disease_lolo"]])
treatment <- read_tsv(paths[["treatment"]])
treatment_linewise <- read_tsv(paths[["treatment_linewise"]])
concordance <- read_tsv(paths[["concordance"]])
candidates <- read_tsv(paths[["candidates"]])
candidate_summary <- read_tsv(paths[["candidate_summary"]])

checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = paste(observed, collapse = ","),
    expected = expected,
    stringsAsFactors = FALSE
  )
}

sample_columns <- setdiff(
  names(raw_counts),
  c("gene_id", "ensembl_gene_id", "gene_name")
)
add_check(
  "frozen_sample_sheet",
  nrow(sample_sheet) == 9L &&
    identical(sample_sheet$array_index, seq_len(9L)),
  nrow(sample_sheet),
  "9 libraries in frozen array order"
)
add_check(
  "raw_count_sample_order",
  identical(sample_columns, sample_sheet$sample_id),
  sample_columns,
  paste(sample_sheet$sample_id, collapse = ",")
)
add_check(
  "raw_count_gene_rows",
  nrow(raw_counts) >= minimum_matched_genes &&
    !anyDuplicated(raw_counts$gene_id),
  nrow(raw_counts),
  paste0(">=", minimum_matched_genes, " unique gene IDs")
)
add_check(
  "analysis_unit_counts",
  metric_value(analysis_summary, "libraries") == "9" &&
    metric_value(analysis_summary, "untreated_donor_lines") == "5" &&
    metric_value(
      analysis_summary,
      "paired_treatment_donor_lines"
    ) == "2",
  c(
    metric_value(analysis_summary, "libraries"),
    metric_value(analysis_summary, "untreated_donor_lines"),
    metric_value(analysis_summary, "paired_treatment_donor_lines")
  ),
  "9 libraries, 5 untreated lines, 2 treatment pairs"
)
add_check(
  "all_quantification_decisions_pass",
  metric_value(analysis_summary, "all_sample_decisions_pass") == "TRUE",
  metric_value(analysis_summary, "all_sample_decisions_pass"),
  "TRUE"
)
add_check(
  "quantification_scope",
  metric_value(analysis_summary, "quantification_method") ==
    expected_quantification_method &&
    toupper(metric_value(
      analysis_summary,
      "whole_genome_alignment"
    )) == expected_whole_genome_alignment,
  c(
    metric_value(analysis_summary, "quantification_method"),
    metric_value(analysis_summary, "whole_genome_alignment")
  ),
  paste(
    expected_quantification_method,
    expected_whole_genome_alignment,
    sep = ","
  )
)
add_check(
  "random_split_absent",
  metric_value(analysis_summary, "random_split_used") == "FALSE" &&
    metric_value(candidate_summary, "random_split_used") == "FALSE",
  c(
    metric_value(analysis_summary, "random_split_used"),
    metric_value(candidate_summary, "random_split_used")
  ),
  "FALSE in count and concordance summaries"
)
add_check(
  "disease_lolo_complete",
  nrow(disease) >= minimum_matched_genes &&
    !anyDuplicated(disease$gene_id) &&
    !anyDuplicated(disease_lolo$gene_id) &&
    all(disease$gene_id %in% disease_lolo$gene_id) &&
    all(disease_lolo$lolo_folds_tested >= 1L) &&
    all(disease_lolo$lolo_folds_tested <= 5L) &&
    sum(disease_lolo$lolo_folds_tested == 5L) >= minimum_matched_genes,
  c(
    nrow(disease),
    nrow(disease_lolo),
    sum(disease_lolo$lolo_folds_tested == 5L)
  ),
  paste0(
    ">=", minimum_matched_genes,
    " disease genes; all represented in LOLO; >=",
    minimum_matched_genes, " genes tested in all 5 folds"
  )
)
add_check(
  "treatment_linewise_complete",
  nrow(treatment) >= minimum_matched_genes &&
    !anyDuplicated(treatment$gene_id) &&
    setequal(treatment$gene_id, treatment_linewise$gene_id) &&
    all(treatment_linewise$donor_lines == 2L),
  c(nrow(treatment), nrow(treatment_linewise)),
  paste0(">=", minimum_matched_genes, " genes with 2 donor-line effects")
)
add_check(
  "concordance_contrasts",
  setequal(
    concordance$contrast,
    c("SMA_vs_CTRL", "R6-Mo_vs_Scramble")
  ) &&
    nrow(concordance) == 2L,
  concordance$contrast,
  "SMA_vs_CTRL,R6-Mo_vs_Scramble"
)
add_check(
  "concordance_metrics_finite",
  all(concordance$matched_genes >= minimum_matched_genes) &&
    all(is.finite(concordance$spearman_effect_correlation)) &&
    all(is.finite(concordance$pearson_effect_correlation)) &&
    all(
      concordance$direction_agreement_fraction >= 0 &
        concordance$direction_agreement_fraction <= 1
    ),
  concordance$matched_genes,
  paste0(
    ">=", minimum_matched_genes,
    " matched genes; finite correlations; direction fraction in [0,1]"
  )
)
allowed_tiers <- c(
  "raw_direction_and_unit_robust",
  "raw_direction_only",
  "raw_tested_not_full_pattern",
  "raw_not_mapped"
)
add_check(
  "candidate_table_integrity",
  nrow(candidates) == expected_candidates &&
    !anyDuplicated(candidates$gene_symbol) &&
    all(candidates$local9_raw_confirmation_tier %in% allowed_tiers),
  nrow(candidates),
  paste0(expected_candidates, " unique candidates with allowed tiers")
)
add_check(
  "candidate_summary_consistency",
  metric_value(candidate_summary, "candidates") ==
    as.character(nrow(candidates)) &&
    metric_value(candidate_summary, "inference_role") ==
      "raw_read_sensitivity_confirmation_not_independent_validation",
  c(
    metric_value(candidate_summary, "candidates"),
    metric_value(candidate_summary, "inference_role")
  ),
  paste0(
    nrow(candidates),
    ",raw_read_sensitivity_confirmation_not_independent_validation"
  )
)

validation <- do.call(rbind, checks)
validation_path <- file.path(
  output_root,
  "GSE290979_local9_raw_analysis_validation.tsv"
)
write_tsv(validation, validation_path)
if (!all(validation$passed)) {
  stop(
    "Local9 raw-output validation failed: ",
    paste(validation$check[!validation$passed], collapse = ", ")
  )
}

completion <- data.frame(
  metric = c(
    "status", "completed_utc", "libraries", "candidate_count",
    "minimum_matched_genes", "random_split_used", "inference_role",
    "quantification_method", "whole_genome_alignment"
  ),
  value = c(
    "COMPLETE",
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    nrow(sample_sheet),
    nrow(candidates),
    min(concordance$matched_genes),
    "FALSE",
    "raw_read_sensitivity_confirmation_not_independent_validation",
    metric_value(analysis_summary, "quantification_method"),
    metric_value(analysis_summary, "whole_genome_alignment")
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  completion,
  file.path(output_root, completion_marker)
)
message(
  "Validated complete local9 ",
  expected_quantification_method, " analysis: ",
  nrow(validation), " integrity checks passed."
)
