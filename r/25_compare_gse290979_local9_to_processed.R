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
source(file.path(ROOT, "r", "common.R"))

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

write_output <- function(frame, name) {
  write.table(
    frame,
    file.path(output_root, name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

collapse_raw_symbols <- function(frame) {
  frame <- frame[
    !is.na(frame$gene_name) & nzchar(frame$gene_name),
    ,
    drop = FALSE
  ]
  ordering <- order(
    frame$gene_name,
    -frame$logCPM,
    frame$PValue,
    frame$gene_id
  )
  frame <- frame[ordering, , drop = FALSE]
  frame[!duplicated(frame$gene_name), , drop = FALSE]
}

contrast_concordance <- function(
  raw,
  processed,
  raw_effect,
  processed_effect,
  label
) {
  joined <- merge(
    raw,
    processed,
    by = "gene_symbol",
    all = FALSE,
    sort = FALSE
  )
  complete <- is.finite(joined[[raw_effect]]) &
    is.finite(joined[[processed_effect]])
  joined <- joined[complete, , drop = FALSE]
  if (nrow(joined) < minimum_matched_genes) {
    stop(
      label, " has only ", nrow(joined),
      " matched genes; minimum is ", minimum_matched_genes
    )
  }
  nonzero <- joined[[raw_effect]] != 0 & joined[[processed_effect]] != 0
  summary <- data.frame(
    contrast = label,
    matched_genes = nrow(joined),
    spearman_effect_correlation = stats::cor(
      joined[[raw_effect]],
      joined[[processed_effect]],
      method = "spearman"
    ),
    pearson_effect_correlation = stats::cor(
      joined[[raw_effect]],
      joined[[processed_effect]],
      method = "pearson"
    ),
    direction_agreement_fraction = mean(
      sign(joined[[raw_effect]][nonzero]) ==
        sign(joined[[processed_effect]][nonzero])
    ),
    raw_nominal_p_lt_0_05 = sum(joined$raw_p_value < 0.05),
    processed_nominal_p_lt_0_05 = sum(
      joined$processed_p_value < 0.05
    ),
    both_nominal_p_lt_0_05 = sum(
      joined$raw_p_value < 0.05 &
        joined$processed_p_value < 0.05
    ),
    stringsAsFactors = FALSE
  )
  list(joined = joined, summary = summary)
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
raw_output_root <- option_value(
  options,
  "raw-output-root",
  file.path(ROOT, "results", "r", "raw_confirmation", "local9")
)
processed_disease_path <- option_value(
  options,
  "processed-disease",
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE290979_SMA_vs_CTRL_pseudobulk.tsv"
  )
)
processed_treatment_path <- option_value(
  options,
  "processed-treatment",
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE290979_R6_vs_scramble_pseudobulk.tsv"
  )
)
candidate_path <- option_value(
  options,
  "candidate-table",
  file.path(
    ROOT,
    "results",
    "r",
    "robustness",
    "human_cross_model_holdout_candidates.tsv"
  )
)
output_root <- option_value(options, "output-root", raw_output_root)
minimum_matched_genes <- suppressWarnings(as.integer(option_value(
  options,
  "minimum-matched-genes",
  "1000"
)))
if (is.na(minimum_matched_genes) || minimum_matched_genes < 1L) {
  stop("--minimum-matched-genes must be a positive integer")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  file.path(raw_output_root, "GSE290979_local9_SMA_vs_CTRL_edgeR.tsv"),
  file.path(
    raw_output_root,
    "GSE290979_local9_R6_vs_scramble_paired_edgeR.tsv"
  ),
  file.path(
    raw_output_root,
    "GSE290979_local9_disease_LOLO_summary.tsv"
  ),
  file.path(
    raw_output_root,
    "GSE290979_local9_treatment_linewise_summary.tsv"
  ),
  processed_disease_path,
  processed_treatment_path,
  candidate_path
)
if (!all(file.exists(required_paths))) {
  stop(
    "Raw concordance inputs are missing: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}

raw_disease <- read.delim(
  required_paths[1L],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw_treatment <- read.delim(
  required_paths[2L],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw_disease_lolo <- read.delim(
  required_paths[3L],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
raw_treatment_linewise <- read.delim(
  required_paths[4L],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
processed_disease <- read.delim(
  processed_disease_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
processed_treatment <- read.delim(
  processed_treatment_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
candidates <- read.delim(
  candidate_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

raw_disease <- collapse_raw_symbols(raw_disease)
raw_treatment <- collapse_raw_symbols(raw_treatment)
raw_disease <- merge(
  raw_disease,
  raw_disease_lolo,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
raw_treatment <- merge(
  raw_treatment,
  raw_treatment_linewise,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

raw_disease_concordance <- data.frame(
  gene_symbol = raw_disease$gene_name,
  raw_gene_id = raw_disease$gene_id,
  raw_log2_effect = raw_disease$logFC,
  raw_logCPM = raw_disease$logCPM,
  raw_p_value = raw_disease$PValue,
  raw_q_value = raw_disease$FDR,
  raw_lolo_folds_tested = raw_disease$lolo_folds_tested,
  raw_lolo_all_same_direction = raw_disease$lolo_all_same_direction,
  stringsAsFactors = FALSE
)
processed_disease_concordance <- data.frame(
  gene_symbol = processed_disease$gene_symbol,
  processed_log2_effect =
    processed_disease$sma_vs_control_log2_effect,
  processed_logCPM = processed_disease$logCPM,
  processed_p_value = processed_disease$p_value,
  processed_q_value = processed_disease$q_value,
  stringsAsFactors = FALSE
)
disease <- contrast_concordance(
  raw_disease_concordance,
  processed_disease_concordance,
  "raw_log2_effect",
  "processed_log2_effect",
  "SMA_vs_CTRL"
)
write_output(
  disease$joined,
  "GSE290979_local9_raw_vs_processed_disease.tsv"
)

raw_treatment_concordance <- data.frame(
  gene_symbol = raw_treatment$gene_name,
  raw_gene_id = raw_treatment$gene_id,
  raw_log2_effect = raw_treatment$logFC,
  raw_logCPM = raw_treatment$logCPM,
  raw_p_value = raw_treatment$PValue,
  raw_q_value = raw_treatment$FDR,
  raw_treatment_lines = raw_treatment$donor_lines,
  raw_treatment_lines_same_direction = raw_treatment$same_direction,
  stringsAsFactors = FALSE
)
processed_treatment_concordance <- data.frame(
  gene_symbol = processed_treatment$gene_symbol,
  processed_log2_effect =
    processed_treatment$r6_vs_scramble_log2_effect,
  processed_logCPM = processed_treatment$logCPM,
  processed_p_value = processed_treatment$p_value,
  processed_q_value = processed_treatment$q_value,
  stringsAsFactors = FALSE
)
treatment <- contrast_concordance(
  raw_treatment_concordance,
  processed_treatment_concordance,
  "raw_log2_effect",
  "processed_log2_effect",
  "R6-Mo_vs_Scramble"
)
write_output(
  treatment$joined,
  "GSE290979_local9_raw_vs_processed_treatment.tsv"
)
concordance_summary <- rbind(disease$summary, treatment$summary)
write_output(
  concordance_summary,
  "GSE290979_local9_raw_vs_processed_concordance_summary.tsv"
)

raw_candidate_disease <- raw_disease_concordance
names(raw_candidate_disease)[names(raw_candidate_disease) != "gene_symbol"] <-
  paste0(
    "local9_disease_",
    names(raw_candidate_disease)[names(raw_candidate_disease) != "gene_symbol"]
  )
raw_candidate_treatment <- raw_treatment_concordance
names(
  raw_candidate_treatment
)[names(raw_candidate_treatment) != "gene_symbol"] <- paste0(
  "local9_treatment_",
  names(raw_candidate_treatment)[
    names(raw_candidate_treatment) != "gene_symbol"
  ]
)
candidate_raw <- merge(
  candidates,
  raw_candidate_disease,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
candidate_raw <- merge(
  candidate_raw,
  raw_candidate_treatment,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
candidate_raw$local9_raw_disease_opposes_natural <- (
  is.finite(candidate_raw$local9_disease_raw_log2_effect) &
    is.finite(candidate_raw$omn_log2_effect_omn_vs_sc) &
    sign(candidate_raw$local9_disease_raw_log2_effect) ==
      -sign(candidate_raw$omn_log2_effect_omn_vs_sc)
)
candidate_raw$local9_raw_r6_aligns_natural <- (
  is.finite(candidate_raw$local9_treatment_raw_log2_effect) &
    is.finite(candidate_raw$omn_log2_effect_omn_vs_sc) &
    sign(candidate_raw$local9_treatment_raw_log2_effect) ==
      sign(candidate_raw$omn_log2_effect_omn_vs_sc)
)
candidate_raw$local9_raw_disease_lolo_robust <- (
  !is.na(candidate_raw$local9_disease_raw_lolo_folds_tested) &
    candidate_raw$local9_disease_raw_lolo_folds_tested == 5L &
  !is.na(candidate_raw$local9_disease_raw_lolo_all_same_direction) &
    candidate_raw$local9_disease_raw_lolo_all_same_direction
)
candidate_raw$local9_raw_treatment_line_robust <- (
  !is.na(
    candidate_raw$
      local9_treatment_raw_treatment_lines_same_direction
  ) &
    candidate_raw$local9_treatment_raw_treatment_lines_same_direction
)
candidate_raw$local9_raw_full_directional_pattern <- (
  candidate_raw$local9_raw_disease_opposes_natural &
    candidate_raw$local9_raw_r6_aligns_natural
)
candidate_raw$local9_raw_direction_and_unit_robust <- (
  candidate_raw$local9_raw_full_directional_pattern &
    candidate_raw$local9_raw_disease_lolo_robust &
    candidate_raw$local9_raw_treatment_line_robust
)
candidate_raw$local9_raw_confirmation_tier <- ifelse(
  candidate_raw$local9_raw_direction_and_unit_robust,
  "raw_direction_and_unit_robust",
  ifelse(
    candidate_raw$local9_raw_full_directional_pattern,
    "raw_direction_only",
    ifelse(
      is.finite(candidate_raw$local9_disease_raw_log2_effect) |
        is.finite(candidate_raw$local9_treatment_raw_log2_effect),
      "raw_tested_not_full_pattern",
      "raw_not_mapped"
    )
  )
)
write_output(
  candidate_raw,
  "GSE290979_local9_raw_candidate_confirmation.tsv"
)
candidate_summary <- data.frame(
  metric = c(
    "candidates", "raw_disease_mapped", "raw_treatment_mapped",
    "raw_full_directional_pattern",
    "raw_direction_and_unit_robust",
    "raw_disease_nominal_p_lt_0_05",
    "raw_treatment_nominal_p_lt_0_05",
    "random_split_used", "inference_role"
  ),
  value = c(
    nrow(candidate_raw),
    sum(is.finite(candidate_raw$local9_disease_raw_log2_effect)),
    sum(is.finite(candidate_raw$local9_treatment_raw_log2_effect)),
    sum(candidate_raw$local9_raw_full_directional_pattern),
    sum(candidate_raw$local9_raw_direction_and_unit_robust),
    sum(candidate_raw$local9_disease_raw_p_value < 0.05, na.rm = TRUE),
    sum(candidate_raw$local9_treatment_raw_p_value < 0.05, na.rm = TRUE),
    "FALSE",
    "raw_read_sensitivity_confirmation_not_independent_validation"
  ),
  stringsAsFactors = FALSE
)
write_output(
  candidate_summary,
  "GSE290979_local9_raw_candidate_confirmation_summary.tsv"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_root,
    "GSE290979_local9_raw_concordance_R_sessionInfo.txt"
  )
)
message(
  "Completed raw-versus-processed concordance: ",
  disease$summary$matched_genes, " disease genes and ",
  treatment$summary$matched_genes, " treatment genes."
)
