script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

analysis_script <- file.path(
  ROOT,
  "r",
  "25_compare_gse290979_local9_to_processed.R"
)
validation_script <- file.path(
  ROOT,
  "r",
  "26_validate_gse290979_local9_outputs.R"
)
stopifnot(length(parse(file = analysis_script)) > 0L)
stopifnot(length(parse(file = validation_script)) > 0L)
analysis_lines <- readLines(analysis_script, warn = FALSE)
required_patterns <- c(
  "stats::cor",
  "raw_read_sensitivity_confirmation_not_independent_validation",
  "human_cross_model_holdout_candidates",
  "raw_direction_and_unit_robust",
  "minimum-matched-genes"
)
stopifnot(all(vapply(
  required_patterns,
  function(pattern) any(grepl(pattern, analysis_lines, fixed = TRUE)),
  logical(1)
)))
forbidden_patterns <- c(
  "(^|[^A-Za-z0-9_])sample[[:space:]]*[(]",
  "(^|[^A-Za-z0-9_])runif[[:space:]]*[(]",
  "(^|[^A-Za-z0-9_])rnorm[[:space:]]*[(]",
  "(^|[^A-Za-z0-9_.])set[.]seed[[:space:]]*[(]"
)
stopifnot(!any(vapply(
  forbidden_patterns,
  function(pattern) any(grepl(pattern, analysis_lines)),
  logical(1)
)))

fixture_root <- tempfile("local9_concordance_fixture_")
raw_root <- file.path(fixture_root, "raw")
processed_root <- file.path(fixture_root, "processed")
output_root <- file.path(fixture_root, "output")
dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

gene_count <- 240L
gene_index <- seq_len(gene_count)
gene_ids <- sprintf("ENSG%011d.1", gene_index)
gene_names <- sprintf("GENE%03d", gene_index)
natural_effect <- ifelse(gene_index %% 2L == 0L, 1, -1) *
  (0.25 + gene_index / 300)
raw_disease_effect <- -natural_effect +
  ((gene_index %% 7L) - 3L) / 100
raw_treatment_effect <- natural_effect +
  ((gene_index %% 5L) - 2L) / 100
processed_disease_effect <- raw_disease_effect * 0.92 +
  ((gene_index %% 11L) - 5L) / 80
processed_treatment_effect <- raw_treatment_effect * 1.04 +
  ((gene_index %% 13L) - 6L) / 90
nominal_p <- pmin(0.99, 0.001 + gene_index / 10000)

write_tsv <- function(frame, path) {
  write.table(
    frame,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

raw_template <- function(effect) {
  data.frame(
    gene_id = gene_ids,
    ensembl_gene_id = sub("[.].*$", "", gene_ids),
    gene_name = gene_names,
    logFC = effect,
    logCPM = 4 + gene_index / 100,
    PValue = nominal_p,
    FDR = pmin(1, nominal_p * 1.5),
    stringsAsFactors = FALSE
  )
}
write_tsv(
  raw_template(raw_disease_effect),
  file.path(raw_root, "GSE290979_local9_SMA_vs_CTRL_edgeR.tsv")
)
write_tsv(
  raw_template(raw_treatment_effect),
  file.path(
    raw_root,
    "GSE290979_local9_R6_vs_scramble_paired_edgeR.tsv"
  )
)
write_tsv(
  data.frame(
    gene_id = c(gene_ids, "ENSG99999999999.1"),
    lolo_folds_tested = c(
      ifelse(gene_index == 1L, 4L, 5L),
      1L
    ),
    lolo_all_same_direction = TRUE,
    stringsAsFactors = FALSE
  ),
  file.path(raw_root, "GSE290979_local9_disease_LOLO_summary.tsv")
)
write_tsv(
  data.frame(
    gene_id = gene_ids,
    donor_lines = 2L,
    same_direction = TRUE,
    stringsAsFactors = FALSE
  ),
  file.path(
    raw_root,
    "GSE290979_local9_treatment_linewise_summary.tsv"
  )
)

processed_disease_path <- file.path(
  processed_root,
  "GSE290979_SMA_vs_CTRL_pseudobulk.tsv"
)
processed_treatment_path <- file.path(
  processed_root,
  "GSE290979_R6_vs_scramble_pseudobulk.tsv"
)
candidate_path <- file.path(
  processed_root,
  "human_cross_model_holdout_candidates.tsv"
)
write_tsv(
  data.frame(
    gene_symbol = gene_names,
    sma_vs_control_log2_effect = processed_disease_effect,
    logCPM = 5 + gene_index / 120,
    p_value = nominal_p,
    q_value = pmin(1, nominal_p * 1.4),
    stringsAsFactors = FALSE
  ),
  processed_disease_path
)
write_tsv(
  data.frame(
    gene_symbol = gene_names,
    r6_vs_scramble_log2_effect = processed_treatment_effect,
    logCPM = 5 + gene_index / 130,
    p_value = nominal_p,
    q_value = pmin(1, nominal_p * 1.3),
    stringsAsFactors = FALSE
  ),
  processed_treatment_path
)
write_tsv(
  data.frame(
    gene_symbol = gene_names[seq_len(12L)],
    omn_log2_effect_omn_vs_sc = natural_effect[seq_len(12L)],
    prior_candidate_rank = seq_len(12L),
    stringsAsFactors = FALSE
  ),
  candidate_path
)

output <- suppressWarnings(system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(
    shQuote(analysis_script),
    paste0(
      "--raw-output-root=",
      normalizePath(raw_root, winslash = "/", mustWork = TRUE)
    ),
    paste0(
      "--processed-disease=",
      normalizePath(processed_disease_path, winslash = "/", mustWork = TRUE)
    ),
    paste0(
      "--processed-treatment=",
      normalizePath(
        processed_treatment_path,
        winslash = "/",
        mustWork = TRUE
      )
    ),
    paste0(
      "--candidate-table=",
      normalizePath(candidate_path, winslash = "/", mustWork = TRUE)
    ),
    paste0(
      "--output-root=",
      normalizePath(output_root, winslash = "/", mustWork = TRUE)
    ),
    "--minimum-matched-genes=100"
  ),
  stdout = TRUE,
  stderr = TRUE
))
status <- attr(output, "status")
if (!is.null(status)) {
  stop(
    "Synthetic local9 raw concordance failed with status ",
    status,
    ":\n",
    paste(output, collapse = "\n")
  )
}

expected_outputs <- c(
  "GSE290979_local9_raw_vs_processed_disease.tsv",
  "GSE290979_local9_raw_vs_processed_treatment.tsv",
  "GSE290979_local9_raw_vs_processed_concordance_summary.tsv",
  "GSE290979_local9_raw_candidate_confirmation.tsv",
  "GSE290979_local9_raw_candidate_confirmation_summary.tsv",
  "GSE290979_local9_raw_concordance_R_sessionInfo.txt"
)
stopifnot(all(file.exists(file.path(output_root, expected_outputs))))

concordance <- read.delim(
  file.path(
    output_root,
    "GSE290979_local9_raw_vs_processed_concordance_summary.tsv"
  ),
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(concordance) == 2L,
  all(concordance$matched_genes == gene_count),
  all(concordance$spearman_effect_correlation > 0.95),
  all(concordance$direction_agreement_fraction == 1)
)
candidate_confirmation <- read.delim(
  file.path(
    output_root,
    "GSE290979_local9_raw_candidate_confirmation.tsv"
  ),
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(candidate_confirmation) == 12L,
  all(candidate_confirmation$local9_raw_full_directional_pattern),
  !candidate_confirmation$local9_raw_direction_and_unit_robust[1L],
  all(candidate_confirmation$local9_raw_direction_and_unit_robust[-1L]),
  candidate_confirmation$local9_raw_confirmation_tier[1L] ==
    "raw_direction_only",
  all(
    candidate_confirmation$local9_raw_confirmation_tier[-1L] ==
      "raw_direction_and_unit_robust"
  )
)
candidate_summary <- read.delim(
  file.path(
    output_root,
    "GSE290979_local9_raw_candidate_confirmation_summary.tsv"
  ),
  stringsAsFactors = FALSE
)
summary_value <- function(metric) {
  candidate_summary$value[match(metric, candidate_summary$metric)]
}
stopifnot(
  summary_value("candidates") == "12",
  summary_value("raw_direction_and_unit_robust") == "11",
  summary_value("random_split_used") == "FALSE",
  summary_value("inference_role") ==
    "raw_read_sensitivity_confirmation_not_independent_validation"
)

raw_input_names <- c(
  "GSE290979_local9_SMA_vs_CTRL_edgeR.tsv",
  "GSE290979_local9_disease_LOLO_summary.tsv",
  "GSE290979_local9_R6_vs_scramble_paired_edgeR.tsv",
  "GSE290979_local9_treatment_linewise_summary.tsv"
)
stopifnot(all(file.copy(
  file.path(raw_root, raw_input_names),
  file.path(output_root, raw_input_names)
)))
sample_ids <- sprintf("SAMPLE%d", seq_len(9L))
sample_sheet_path <- file.path(fixture_root, "sample_sheet.tsv")
write_tsv(
  data.frame(
    array_index = seq_len(9L),
    sample_id = sample_ids,
    stringsAsFactors = FALSE
  ),
  sample_sheet_path
)
raw_count_fixture <- data.frame(
  gene_id = gene_ids,
  ensembl_gene_id = sub("[.].*$", "", gene_ids),
  gene_name = gene_names,
  stringsAsFactors = FALSE
)
for (sample_index in seq_along(sample_ids)) {
  raw_count_fixture[[sample_ids[sample_index]]] <-
    100L + gene_index + sample_index
}
write_tsv(
  raw_count_fixture,
  file.path(output_root, "GSE290979_local9_raw_gene_counts.tsv")
)
write_tsv(
  data.frame(
    metric = c(
      "libraries", "untreated_donor_lines",
      "paired_treatment_donor_lines", "random_split_used",
      "all_sample_decisions_pass", "quantification_method",
      "whole_genome_alignment"
    ),
    value = c(
      "9", "5", "2", "FALSE", "TRUE",
      "genome_alignment_featurecounts", "TRUE"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_root, "GSE290979_local9_analysis_summary.tsv")
)
validation_output <- suppressWarnings(system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(
    shQuote(validation_script),
    paste0(
      "--output-root=",
      normalizePath(output_root, winslash = "/", mustWork = TRUE)
    ),
    paste0(
      "--sample-sheet=",
      normalizePath(sample_sheet_path, winslash = "/", mustWork = TRUE)
    ),
    "--expected-candidates=12",
    "--minimum-matched-genes=100"
  ),
  stdout = TRUE,
  stderr = TRUE
))
validation_status <- attr(validation_output, "status")
if (!is.null(validation_status)) {
  stop(
    "Synthetic local9 output validation failed with status ",
    validation_status,
    ":\n",
    paste(validation_output, collapse = "\n")
  )
}
validation <- read.delim(
  file.path(
    output_root,
    "GSE290979_local9_raw_analysis_validation.tsv"
  ),
  stringsAsFactors = FALSE
)
completion <- read.delim(
  file.path(output_root, "LOCAL9_RAW_ANALYSIS_COMPLETE.tsv"),
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(validation) >= 10L,
  all(validation$passed),
  completion$value[completion$metric == "status"] == "COMPLETE",
  completion$value[completion$metric == "random_split_used"] == "FALSE",
  completion$value[completion$metric == "quantification_method"] ==
    "genome_alignment_featurecounts",
  completion$value[completion$metric == "whole_genome_alignment"] == "TRUE"
)

unlink(fixture_root, recursive = TRUE, force = TRUE)
cat("local9 raw concordance and output validation tests passed\n")
