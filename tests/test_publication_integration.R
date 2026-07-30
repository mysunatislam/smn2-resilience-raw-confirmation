script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) != 1L) {
  stop("This test must be run with Rscript")
}
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

write_fixture <- function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

fixture_root <- tempfile("publication_integration_")
raw_root <- file.path(fixture_root, "raw")
output_root <- file.path(fixture_root, "output")
manuscript_root <- file.path(fixture_root, "manuscript")
dir.create(raw_root, recursive = TRUE)

completion <- data.frame(
  metric = c(
    "status", "libraries", "candidate_count", "random_split_used",
    "inference_role", "quantification_method", "whole_genome_alignment"
  ),
  value = c(
    "COMPLETE", "9", "6", "FALSE",
    "raw_read_sensitivity_confirmation_not_independent_validation",
    "salmon_transcriptome_quasimapping", "FALSE"
  ),
  stringsAsFactors = FALSE
)
write_fixture(
  completion,
  file.path(raw_root, "LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv")
)

validation <- data.frame(
  check = paste0("check_", seq_len(13L)),
  passed = TRUE,
  observed = "pass",
  expected = "pass",
  stringsAsFactors = FALSE
)
write_fixture(
  validation,
  file.path(raw_root, "GSE290979_local9_raw_analysis_validation.tsv")
)

concordance <- data.frame(
  contrast = c("SMA_vs_CTRL", "R6-Mo_vs_Scramble"),
  matched_genes = c(20L, 20L),
  spearman_effect_correlation = c(0.8, 0.6),
  pearson_effect_correlation = c(0.8, 0.6),
  direction_agreement_fraction = c(0.8, 0.7),
  raw_nominal_p_lt_0_05 = c(2L, 1L),
  processed_nominal_p_lt_0_05 = c(2L, 1L),
  both_nominal_p_lt_0_05 = c(1L, 0L),
  stringsAsFactors = FALSE
)
write_fixture(
  concordance,
  file.path(
    raw_root,
    "GSE290979_local9_raw_vs_processed_concordance_summary.tsv"
  )
)

gene_symbols <- LETTERS[1:6]
processed_robust <- c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)
raw_robust <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
raw_direction <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
candidates <- data.frame(
  gene_symbol = gene_symbols,
  cross_model_resilience_rank = seq_along(gene_symbols),
  exploratory_cross_model_rank = seq_along(gene_symbols),
  cross_model_resilience_score = seq(0.9, 0.4, length.out = 6),
  omn_log2_effect_omn_vs_sc = rep(1, 6),
  omn_p_value = rep(0.01, 6),
  omn_q_value = rep(0.5, 6),
  disease_sma_vs_control_log2_effect = rep(-1, 6),
  treatment_r6_vs_scramble_log2_effect = rep(1, 6),
  gse108094_sma_vs_control_log2_effect = rep(-1, 6),
  strict_splicing = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  all_estimable_biological_unit_checks_robust = processed_robust,
  gse69175_disease_opposes_natural = rep(FALSE, 6),
  gse290980_mn_opposes_natural_resistance = rep(FALSE, 6),
  gse243076_adult_mn_detected = rep(TRUE, 6),
  gse243076_adult_mn_localized = rep(FALSE, 6),
  local9_disease_raw_gene_id = c("g1", "g2", "g3", "g4", "g5", ""),
  local9_disease_raw_log2_effect = c(-1, -1, -1, -1, 1, NA),
  local9_disease_raw_p_value = rep(0.2, 6),
  local9_disease_raw_q_value = rep(0.8, 6),
  local9_disease_raw_lolo_folds_tested = c(5, 5, 5, 5, 5, NA),
  local9_disease_raw_lolo_all_same_direction =
    c(TRUE, TRUE, TRUE, FALSE, FALSE, NA),
  local9_treatment_raw_gene_id = c("g1", "g2", "g3", "g4", "g5", ""),
  local9_treatment_raw_log2_effect = c(1, 1, 1, 1, -1, NA),
  local9_treatment_raw_p_value = rep(0.2, 6),
  local9_treatment_raw_q_value = rep(0.8, 6),
  local9_treatment_raw_treatment_lines = c(2, 2, 2, 2, 2, NA),
  local9_treatment_raw_treatment_lines_same_direction =
    c(TRUE, TRUE, FALSE, FALSE, FALSE, NA),
  local9_raw_disease_opposes_natural = raw_direction,
  local9_raw_r6_aligns_natural = raw_direction,
  local9_raw_disease_lolo_robust = raw_robust,
  local9_raw_treatment_line_robust = raw_robust,
  local9_raw_full_directional_pattern = raw_direction,
  local9_raw_direction_and_unit_robust = raw_robust,
  local9_raw_confirmation_tier = c(
    "raw_direction_and_unit_robust",
    "raw_direction_and_unit_robust",
    "raw_direction_only",
    "raw_direction_only",
    "raw_tested_not_full_pattern",
    "raw_not_mapped"
  ),
  stringsAsFactors = FALSE
)
write_fixture(
  candidates,
  file.path(raw_root, "GSE290979_local9_raw_candidate_confirmation.tsv")
)

set.seed(1)
for (contrast_name in c("disease", "treatment")) {
  x <- seq(-2, 2, length.out = 20)
  joined <- data.frame(
    gene_symbol = paste0("gene", seq_along(x)),
    raw_log2_effect = x,
    processed_log2_effect = x + stats::rnorm(length(x), sd = 0.2),
    stringsAsFactors = FALSE
  )
  file_name <- if (contrast_name == "disease") {
    "GSE290979_local9_raw_vs_processed_disease.tsv"
  } else {
    "GSE290979_local9_raw_vs_processed_treatment.tsv"
  }
  write_fixture(joined, file.path(raw_root, file_name))
}

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
}
arguments <- c(
  file.path(ROOT, "r", "28_integrate_local9_salmon_publication.R"),
  paste0("--raw-root=", raw_root),
  paste0("--output-root=", output_root),
  paste0("--manuscript-root=", manuscript_root),
  "--expected-candidates=6"
)
status <- system2(rscript, arguments)
if (!identical(status, 0L)) {
  stop("Publication integration script failed with status ", status)
}

result <- read.delim(
  file.path(output_root, "human_raw_confirmed_publication_candidates.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_tiers <- c(
  "tier_1_raw_and_processed_unit_robust",
  "tier_2_raw_unit_robust",
  "tier_3_processed_unit_robust",
  "tier_4_raw_direction_only",
  "tier_5_other_cross_model",
  "tier_5_other_cross_model"
)
if (!identical(result$publication_evidence_tier, expected_tiers)) {
  stop(
    "Unexpected tier ordering: ",
    paste(result$publication_evidence_tier, collapse = ",")
  )
}
if (!identical(result$gene_symbol, gene_symbols)) {
  stop("Frozen-rank tie-breaking changed unexpectedly")
}
if (!identical(result$frozen_discovery_rank_unchanged, seq_len(6L))) {
  stop("Frozen discovery ranks were modified")
}
if (any(result$raw_confirmation_is_independent_validation)) {
  stop("Raw confirmation was incorrectly labeled independent")
}

summary <- read.delim(
  file.path(output_root, "publication_analysis_summary.tsv"),
  stringsAsFactors = FALSE
)
summary_value <- setNames(summary$value, summary$metric)
expected_counts <- c(
  tier_1_candidates = "1",
  tier_2_candidates = "1",
  tier_3_candidates = "1",
  tier_4_candidates = "1",
  tier_5_candidates = "2"
)
if (!identical(
  unname(summary_value[names(expected_counts)]),
  unname(expected_counts)
)) {
  stop("Synthetic publication tier counts are incorrect")
}

expected_outputs <- c(
  file.path(manuscript_root, "supplementary_table_S1_candidates.tsv"),
  file.path(manuscript_root, "supplementary_table_S2_summary.tsv"),
  file.path(manuscript_root, "supplementary_table_S3_concordance.tsv"),
  file.path(output_root, "R_sessionInfo.txt"),
  file.path(
    manuscript_root, "figures",
    "figure_1_raw_vs_processed_concordance.png"
  ),
  file.path(
    manuscript_root, "figures",
    "figure_2_candidate_evidence_matrix.png"
  )
)
if (!all(file.exists(expected_outputs))) {
  stop(
    "Publication outputs are missing: ",
    paste(expected_outputs[!file.exists(expected_outputs)], collapse = ", ")
  )
}
if (any(file.info(expected_outputs)$size <= 0)) {
  stop("One or more publication outputs are empty")
}

cat("Publication integration test passed\n")
