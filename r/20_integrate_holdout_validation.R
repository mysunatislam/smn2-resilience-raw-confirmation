source(file.path("r", "common.R"))

integrated <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "human_cross_model_cell_resolved_gene_table.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gse69175 <- read.delim(
  file.path(
    ROOT, "results", "r", "external_validation",
    "GSE69175_expression_SMA_vs_control.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
omn_lodo <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "GSE93939_primary_LODO_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
omn_same_lodo <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "GSE93939_HiSeq2000_LODO_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
disease_lolo <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "GSE290979_disease_LOLO_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
treatment_lines <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "GSE290979_treatment_linewise_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(integrated$gene_symbol),
  !anyDuplicated(gse69175$gene_symbol),
  !anyDuplicated(omn_lodo$gene_symbol),
  !anyDuplicated(omn_same_lodo$gene_symbol),
  !anyDuplicated(disease_lolo$gene_symbol),
  !anyDuplicated(treatment_lines$gene_symbol)
)

names(gse69175)[names(gse69175) != "gene_symbol"] <- paste0(
  "gse69175_",
  names(gse69175)[names(gse69175) != "gene_symbol"]
)
names(treatment_lines)[names(treatment_lines) != "gene_symbol"] <- paste0(
  "treatment_linewise_",
  names(treatment_lines)[names(treatment_lines) != "gene_symbol"]
)
holdout <- Reduce(function(left, right) {
  merge(left, right, by = "gene_symbol", all.x = TRUE, sort = FALSE)
}, list(
  integrated,
  gse69175,
  omn_lodo,
  omn_same_lodo,
  disease_lolo,
  treatment_lines
))
stopifnot(nrow(holdout) == nrow(integrated))

finite <- function(values) !is.na(values) & is.finite(values)
natural <- holdout$omn_log2_effect_omn_vs_sc
gse69175_effect <- holdout$gse69175_sma_vs_control_log2_effect
gse69175_usable <-
  !is.na(holdout$gse69175_external_direction_usable) &
  holdout$gse69175_external_direction_usable & finite(gse69175_effect)
holdout$gse69175_direction_usable <- gse69175_usable
holdout$gse69175_disease_opposes_natural <-
  gse69175_usable & finite(natural) & natural * gse69175_effect < 0

logical_or_false <- function(values) !is.na(values) & values
holdout$gse93939_primary_lodo_robust <- logical_or_false(
  holdout$primary_lodo_all_folds_same_direction
)
holdout$gse93939_same_platform_lodo_robust <- logical_or_false(
  holdout$same_platform_lodo_all_folds_same_direction
)
holdout$gse290979_disease_lolo_robust <- logical_or_false(
  holdout$disease_lolo_all_folds_same_direction
)
holdout$gse290979_treatment_line_robust <- logical_or_false(
  holdout$treatment_linewise_both_lines_match_full_direction
)
holdout$biological_unit_robustness_count <- rowSums(cbind(
  holdout$gse93939_primary_lodo_robust,
  holdout$gse290979_disease_lolo_robust,
  holdout$gse290979_treatment_line_robust
))
holdout$all_estimable_biological_unit_checks_robust <-
  holdout$biological_unit_robustness_count == 3L

# Dataset-level sensitivity is deterministic. GSE290979 is omitted as one
# source, so its disease and treatment components leave together.
eligible <- logical_or_false(holdout$cross_model_ranking_eligible)
components <- cbind(
  GSE93939 = holdout$natural_resistance_percentile,
  GSE290979_disease = holdout$gse290979_sma_depletion_percentile,
  GSE290979_treatment = holdout$gse290979_r6_alignment_percentile,
  GSE108094 = holdout$gse108094_sma_depletion_percentile
)
holdout$lodo_dataset_omit_GSE93939_score <- rowMeans(
  components[, c("GSE290979_disease", "GSE290979_treatment", "GSE108094")],
  na.rm = FALSE
)
holdout$lodo_dataset_omit_GSE290979_score <- rowMeans(
  components[, c("GSE93939", "GSE108094")],
  na.rm = FALSE
)
holdout$lodo_dataset_omit_GSE108094_score <- rowMeans(
  components[, c("GSE93939", "GSE290979_disease", "GSE290979_treatment")],
  na.rm = FALSE
)

rank_eligible <- function(values, eligible) {
  result <- rep(NA_integer_, length(values))
  use <- which(eligible & finite(values))
  result[use] <- rank(-values[use], ties.method = "min")
  result
}
holdout$lodo_dataset_omit_GSE93939_rank <- rank_eligible(
  holdout$lodo_dataset_omit_GSE93939_score,
  eligible
)
holdout$lodo_dataset_omit_GSE290979_rank <- rank_eligible(
  holdout$lodo_dataset_omit_GSE290979_score,
  eligible
)
holdout$lodo_dataset_omit_GSE108094_rank <- rank_eligible(
  holdout$lodo_dataset_omit_GSE108094_score,
  eligible
)
dataset_ranks <- cbind(
  holdout$lodo_dataset_omit_GSE93939_rank,
  holdout$lodo_dataset_omit_GSE290979_rank,
  holdout$lodo_dataset_omit_GSE108094_rank
)
holdout$lodo_dataset_median_rank <- rep(NA_real_, nrow(holdout))
holdout$lodo_dataset_worst_rank <- rep(NA_real_, nrow(holdout))
holdout$lodo_dataset_median_rank[eligible] <- apply(
  dataset_ranks[eligible, , drop = FALSE],
  1,
  median
)
holdout$lodo_dataset_worst_rank[eligible] <- apply(
  dataset_ranks[eligible, , drop = FALSE],
  1,
  max
)
ranked_genes <- sum(eligible)
holdout$lodo_dataset_worst_rank_stability <- ifelse(
  eligible,
  1 - (holdout$lodo_dataset_worst_rank - 1) / (ranked_genes - 1),
  NA_real_
)

holdout$holdout_evidence_tier <- "not_cross_model_ranked"
holdout$holdout_evidence_tier[eligible] <- "ranked_context"
holdout$holdout_evidence_tier[
  eligible & holdout$all_estimable_biological_unit_checks_robust
] <- "biological_unit_direction_robust"
holdout$holdout_evidence_tier[
  eligible & holdout$all_estimable_biological_unit_checks_robust &
    holdout$gse69175_disease_opposes_natural
] <- "unit_robust_plus_GSE69175_direction"

holdout <- holdout[order(
  -as.integer(eligible),
  -holdout$biological_unit_robustness_count,
  -as.integer(holdout$gse69175_disease_opposes_natural),
  holdout$lodo_dataset_worst_rank,
  holdout$gene_symbol,
  na.last = TRUE
), , drop = FALSE]
rownames(holdout) <- NULL

candidates <- holdout[
  logical_or_false(holdout$exploratory_natural_plus_full_cross_model),
  ,
  drop = FALSE
]
candidates <- candidates[order(
  -candidates$biological_unit_robustness_count,
  -as.integer(candidates$gse69175_disease_opposes_natural),
  candidates$lodo_dataset_worst_rank,
  candidates$cross_model_resilience_rank
), , drop = FALSE]
stopifnot(nrow(candidates) == 37L)

write_tsv(
  holdout,
  "results/r/robustness/human_cross_model_holdout_gene_table.tsv"
)
write_tsv(
  candidates,
  "results/r/robustness/human_cross_model_holdout_candidates.tsv"
)

summary <- data.frame(
  metric = c(
    "ranked_genes", "leave_one_dataset_out_folds",
    "exploratory_candidates",
    "candidates_primary_GSE93939_LODO_robust",
    "candidates_GSE290979_disease_LOLO_robust",
    "candidates_GSE290979_treatment_both_lines_robust",
    "candidates_all_biological_unit_checks_robust",
    "candidates_GSE69175_direction_usable",
    "candidates_GSE69175_opposes_natural",
    "candidates_unit_robust_plus_GSE69175_direction",
    "GSE69175_LOLO_interpretation"
  ),
  value = c(
    ranked_genes, 3L, nrow(candidates),
    sum(candidates$gse93939_primary_lodo_robust),
    sum(candidates$gse290979_disease_lolo_robust),
    sum(candidates$gse290979_treatment_line_robust),
    sum(candidates$all_estimable_biological_unit_checks_robust),
    sum(candidates$gse69175_direction_usable),
    sum(candidates$gse69175_disease_opposes_natural),
    sum(
      candidates$all_estimable_biological_unit_checks_robust &
        candidates$gse69175_disease_opposes_natural
    ),
    paste(
      "Not estimable: one control line and one SMA line; GSE69175 is",
      "directional sensitivity evidence and never a discovery fold."
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "results/r/robustness/human_holdout_validation_summary.tsv")

top <- head(candidates, 25L)
plot_matrix <- cbind(
  `GSE93939 donor LODO` = top$primary_lodo_holdout_sign_agreement_fraction,
  `GSE290979 disease LOLO` = top$disease_lolo_holdout_sign_agreement_fraction,
  `GSE290979 treatment lines` =
    top$treatment_linewise_line_sign_agreement_fraction,
  `GSE69175 direction` = ifelse(
    top$gse69175_direction_usable,
    as.numeric(top$gse69175_disease_opposes_natural),
    NA_real_
  ),
  `Leave-dataset-out stability` = top$lodo_dataset_worst_rank_stability
)
rownames(plot_matrix) <- top$gene_symbol
png(
  file.path(
    ROOT, "results", "r", "figures",
    "human_holdout_robustness_top25.png"
  ),
  width = 6333,
  height = 5833,
  res = 600
)
par(mar = c(13, 12, 5, 3))
image(
  x = seq_len(ncol(plot_matrix)),
  y = seq_len(nrow(plot_matrix)),
  z = t(plot_matrix[nrow(plot_matrix):1, , drop = FALSE]),
  col = hcl.colors(60, "YlGnBu"),
  zlim = c(0, 1),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = ""
)
axis(1, at = seq_len(ncol(plot_matrix)), labels = colnames(plot_matrix), las = 2)
axis(
  2,
  at = seq_len(nrow(plot_matrix)),
  labels = rev(rownames(plot_matrix)),
  las = 1
)
for (column in seq_len(ncol(plot_matrix))) {
  for (row in seq_len(nrow(plot_matrix))) {
    value <- plot_matrix[nrow(plot_matrix) - row + 1L, column]
    if (is.na(value)) {
      rect(
        column - 0.5,
        row - 0.5,
        column + 0.5,
        row + 0.5,
        col = "#D9D9D9",
        border = NA
      )
      text(column, row, "NA", col = "#444444", cex = 0.70)
      next
    }
    text(
      column,
      row,
      sprintf("%.2f", value),
      col = ifelse(value < 0.55, "white", "black"),
      cex = 0.70
    )
  }
}
box()
dev.off()

write_session_info(
  "results/r/robustness/holdout_integration_R_sessionInfo.txt"
)
cat(
  "Holdout validation integration complete:",
  sum(candidates$all_estimable_biological_unit_checks_robust),
  "of 37 candidates pass all estimable unit-direction checks and",
  sum(candidates$gse69175_disease_opposes_natural),
  "have GSE69175 directional support\n"
)
