source(file.path("r", "common.R"))

core <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "integration",
    "human_R_two_track_gene_table.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
external <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "external_validation",
    "GSE108094_expression_SMA_vs_control.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  !anyDuplicated(core$gene_symbol),
  !anyDuplicated(external$gene_symbol),
  all(c(
    "omn_log2_effect_omn_vs_sc",
    "disease_sma_vs_control_log2_effect",
    "treatment_r6_vs_scramble_log2_effect"
  ) %in% names(core)),
  all(c(
    "sma_vs_control_log2_effect", "external_direction_usable"
  ) %in% names(external))
)

names(external)[names(external) != "gene_symbol"] <- paste0(
  "gse108094_",
  names(external)[names(external) != "gene_symbol"]
)
cross_model <- merge(
  core,
  external,
  by = "gene_symbol",
  all = TRUE,
  sort = FALSE
)

finite_effect <- function(values) !is.na(values) & is.finite(values)
natural_effect <- cross_model$omn_log2_effect_omn_vs_sc
organoid_disease_effect <-
  cross_model$disease_sma_vs_control_log2_effect
organoid_treatment_effect <-
  cross_model$treatment_r6_vs_scramble_log2_effect
external_disease_effect <-
  cross_model$gse108094_sma_vs_control_log2_effect
external_usable <-
  !is.na(cross_model$gse108094_external_direction_usable) &
  cross_model$gse108094_external_direction_usable

cross_model$complete_core_effects <-
  finite_effect(natural_effect) &
  finite_effect(organoid_disease_effect) &
  finite_effect(organoid_treatment_effect)
cross_model$complete_cross_model_effects <-
  cross_model$complete_core_effects &
  finite_effect(external_disease_effect) & external_usable
cross_model$approved_hgnc_gene <-
  !is.na(cross_model$hgnc_status) & cross_model$hgnc_status == "Approved"
cross_model$cross_model_ranking_eligible <-
  cross_model$complete_cross_model_effects &
  cross_model$approved_hgnc_gene

percentile_rank <- function(values, eligible) {
  output <- rep(NA_real_, length(values))
  use <- eligible & finite_effect(values)
  count <- sum(use)
  if (count > 0L) {
    output[use] <- (rank(values[use], ties.method = "average") - 0.5) /
      count
  }
  output
}

eligible <- cross_model$cross_model_ranking_eligible
cross_model$natural_resistance_percentile <- percentile_rank(
  natural_effect,
  eligible
)
cross_model$gse290979_sma_depletion_percentile <- percentile_rank(
  -organoid_disease_effect,
  eligible
)
cross_model$gse290979_r6_alignment_percentile <- percentile_rank(
  organoid_treatment_effect,
  eligible
)
cross_model$gse108094_sma_depletion_percentile <- percentile_rank(
  -external_disease_effect,
  eligible
)
cross_model$cross_model_resilience_score <- rowMeans(cbind(
  cross_model$natural_resistance_percentile,
  cross_model$gse290979_sma_depletion_percentile,
  cross_model$gse290979_r6_alignment_percentile,
  cross_model$gse108094_sma_depletion_percentile
), na.rm = FALSE)

opposite_sign <- function(left, right) {
  finite_effect(left) & finite_effect(right) & left * right < 0
}
same_sign <- function(left, right) {
  finite_effect(left) & finite_effect(right) & left * right > 0
}
cross_model$gse290979_disease_opposes_natural <- opposite_sign(
  natural_effect,
  organoid_disease_effect
)
cross_model$gse290979_r6_aligns_natural <- same_sign(
  natural_effect,
  organoid_treatment_effect
)
cross_model$gse108094_disease_opposes_natural <- opposite_sign(
  natural_effect,
  external_disease_effect
)
cross_model$directional_support_count <- rowSums(cbind(
  cross_model$gse290979_disease_opposes_natural,
  cross_model$gse290979_r6_aligns_natural,
  cross_model$gse108094_disease_opposes_natural
))
cross_model$full_directional_resilience_pattern <-
  cross_model$complete_cross_model_effects &
  natural_effect > 0 &
  organoid_disease_effect < 0 &
  organoid_treatment_effect > 0 &
  external_disease_effect < 0
cross_model$full_directional_vulnerability_pattern <-
  cross_model$complete_cross_model_effects &
  natural_effect < 0 &
  organoid_disease_effect > 0 &
  organoid_treatment_effect < 0 &
  external_disease_effect > 0
cross_model$exploratory_natural_resistance_support <-
  !is.na(cross_model$exploratory_omn_positive_p05) &
  cross_model$exploratory_omn_positive_p05
cross_model$exploratory_natural_plus_full_cross_model <-
  cross_model$exploratory_natural_resistance_support &
  cross_model$full_directional_resilience_pattern

rank_order <- order(
  -cross_model$cross_model_resilience_score,
  cross_model$gene_symbol,
  na.last = NA
)
cross_model$cross_model_resilience_rank <- NA_integer_
cross_model$cross_model_resilience_rank[rank_order] <- seq_along(rank_order)
exploratory_rank_order <- which(
  cross_model$exploratory_natural_plus_full_cross_model
)
exploratory_rank_order <- exploratory_rank_order[order(
  -cross_model$cross_model_resilience_score[exploratory_rank_order],
  cross_model$gene_symbol[exploratory_rank_order]
)]
cross_model$exploratory_cross_model_rank <- NA_integer_
cross_model$exploratory_cross_model_rank[exploratory_rank_order] <-
  seq_along(exploratory_rank_order)

cross_model$cross_model_evidence_class <- "context_only"
cross_model$cross_model_evidence_class[
  cross_model$complete_cross_model_effects
] <- "four_model_effect_context"
cross_model$cross_model_evidence_class[
  cross_model$full_directional_resilience_pattern
] <- "full_directional_resilience_pattern"
cross_model$cross_model_evidence_class[
  cross_model$exploratory_natural_plus_full_cross_model
] <- "exploratory_natural_plus_full_directional_pattern"
cross_model$cross_model_evidence_class[
  cross_model$full_directional_vulnerability_pattern
] <- "full_directional_vulnerability_pattern"

cross_model <- cross_model[order(
  -as.integer(cross_model$cross_model_ranking_eligible),
  cross_model$cross_model_resilience_rank,
  cross_model$gene_symbol,
  na.last = TRUE
), , drop = FALSE]
rownames(cross_model) <- NULL

pairwise_correlation <- function(
  comparison,
  x,
  y,
  x_label,
  y_label,
  additional_filter = rep(TRUE, length(x))
) {
  use <- finite_effect(x) & finite_effect(y) & additional_filter
  test <- suppressWarnings(cor.test(
    x[use], y[use],
    method = "spearman",
    exact = FALSE
  ))
  data.frame(
    comparison = comparison,
    x_effect = x_label,
    y_effect = y_label,
    genes = sum(use),
    spearman_rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

correlations <- do.call(rbind, list(
  pairwise_correlation(
    "natural_vs_organoid_disease_opposition",
    natural_effect, -organoid_disease_effect,
    "GSE93939 OMN minus spinal",
    "negative GSE290979 SMA minus control"
  ),
  pairwise_correlation(
    "natural_vs_organoid_smn_response",
    natural_effect, organoid_treatment_effect,
    "GSE93939 OMN minus spinal",
    "GSE290979 R6 minus scramble"
  ),
  pairwise_correlation(
    "natural_vs_external_disease_opposition",
    natural_effect, -external_disease_effect,
    "GSE93939 OMN minus spinal",
    "negative GSE108094 SMA minus control",
    external_usable
  ),
  pairwise_correlation(
    "organoid_vs_external_sma_disease",
    organoid_disease_effect, external_disease_effect,
    "GSE290979 SMA minus control",
    "GSE108094 SMA minus control",
    external_usable
  ),
  pairwise_correlation(
    "organoid_disease_vs_inverse_smn_response",
    organoid_disease_effect, -organoid_treatment_effect,
    "GSE290979 SMA minus control",
    "negative GSE290979 R6 minus scramble"
  )
))

strict_splicing <-
  !is.na(cross_model$splicing_strict_event_count) &
  cross_model$splicing_strict_event_count > 0
strict_cross_model <- cross_model[strict_splicing, , drop = FALSE]
exploratory_shortlist <- cross_model[
  cross_model$exploratory_natural_plus_full_cross_model,
  ,
  drop = FALSE
]
exploratory_shortlist <- exploratory_shortlist[order(
  exploratory_shortlist$exploratory_cross_model_rank
), , drop = FALSE]

primary_panel <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_primary_validation_panel.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
primary_panel_cross_model <- merge(
  primary_panel[, c(
    "panel_order", "current_gene_symbol", "event_type",
    "validation_priority_score"
  )],
  cross_model,
  by.x = "current_gene_symbol",
  by.y = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
primary_panel_cross_model <- primary_panel_cross_model[order(
  primary_panel_cross_model$panel_order
), , drop = FALSE]

bridge_panel <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_exploratory_OMN_bridge_panel.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
bridge_cross_model <- merge(
  bridge_panel[, c(
    "current_gene_symbol", "event_type", "validation_priority_score"
  )],
  cross_model,
  by.x = "current_gene_symbol",
  by.y = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
bridge_cross_model <- bridge_cross_model[order(
  bridge_cross_model$cross_model_resilience_rank
), , drop = FALSE]

summary <- data.frame(
  metric = c(
    "genes_in_union", "genes_with_three_core_effects",
    "genes_with_four_model_effects", "approved_genes_ranked",
    "full_directional_resilience_pattern",
    "full_directional_vulnerability_pattern",
    "exploratory_natural_plus_full_cross_model",
    "strict_splicing_genes_with_external_direction",
    "primary_panel_genes_with_external_direction",
    "bridge_genes_with_external_direction"
  ),
  value = c(
    nrow(cross_model), sum(cross_model$complete_core_effects),
    sum(cross_model$complete_cross_model_effects),
    sum(cross_model$cross_model_ranking_eligible),
    sum(cross_model$full_directional_resilience_pattern),
    sum(cross_model$full_directional_vulnerability_pattern),
    sum(cross_model$exploratory_natural_plus_full_cross_model),
    sum(strict_cross_model$complete_cross_model_effects),
    sum(primary_panel_cross_model$complete_cross_model_effects),
    sum(bridge_cross_model$complete_cross_model_effects)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  cross_model,
  "results/r/integration/human_cross_model_resilience_gene_table.tsv"
)
write_tsv(
  strict_cross_model,
  "results/r/integration/human_cross_model_strict_splicing_annotation.tsv"
)
write_tsv(
  exploratory_shortlist,
  "results/r/integration/human_cross_model_exploratory_resilience_candidates.tsv"
)
write_tsv(
  primary_panel_cross_model,
  "results/r/integration/human_cross_model_primary_validation_panel.tsv"
)
write_tsv(
  bridge_cross_model,
  "results/r/integration/human_cross_model_exploratory_bridge.tsv"
)
write_tsv(
  correlations,
  "results/r/integration/human_cross_model_effect_correlations.tsv"
)
write_tsv(
  summary,
  "results/r/integration/human_cross_model_resilience_summary.tsv"
)

top <- head(exploratory_shortlist, 25L)
component_matrix <- as.matrix(top[, c(
  "natural_resistance_percentile",
  "gse290979_sma_depletion_percentile",
  "gse290979_r6_alignment_percentile",
  "gse108094_sma_depletion_percentile"
)])
strict_top <-
  !is.na(top$splicing_strict_event_count) &
  top$splicing_strict_event_count > 0
rownames(component_matrix) <- paste0(
  top$gene_symbol,
  ifelse(strict_top, " *", "")
)

png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "human_cross_model_resilience_top25.png"
  ),
  width = 6333,
  height = 5667,
  res = 600
)
par(mar = c(12, 12, 5, 3))
image(
  x = seq_len(ncol(component_matrix)),
  y = seq_len(nrow(component_matrix)),
  z = t(component_matrix[nrow(component_matrix):1, , drop = FALSE]),
  col = hcl.colors(60, "YlGnBu", rev = FALSE),
  zlim = c(0, 1),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = ""
)
axis(
  1,
  at = seq_len(ncol(component_matrix)),
  labels = c(
    "Natural resistance\nGSE93939",
    "SMA depletion\nGSE290979",
    "SMN-response alignment\nGSE290979",
    "External SMA depletion\nGSE108094"
  ),
  las = 2
)
axis(
  2,
  at = seq_len(nrow(component_matrix)),
  labels = rev(rownames(component_matrix)),
  las = 1
)
for (column in seq_len(ncol(component_matrix))) {
  for (row in seq_len(nrow(component_matrix))) {
    value <- component_matrix[nrow(component_matrix) - row + 1L, column]
    text(
      column,
      row,
      sprintf("%.2f", value),
      col = ifelse(value < 0.60, "white", "black"),
      cex = 0.72
    )
  }
}
mtext(
  "* GSE290979 strict splice-restoration gene",
  side = 1,
  line = 10,
  adj = 0,
  cex = 0.8
)
box()
dev.off()

plot_rows <- cross_model$cross_model_ranking_eligible
x <- cross_model$natural_resistance_percentile[plot_rows]
y <- cross_model$gse108094_sma_depletion_percentile[plot_rows]
strict_plot <- strict_splicing[plot_rows]
full_plot <- cross_model$full_directional_resilience_pattern[plot_rows]
colors <- ifelse(
  strict_plot,
  "#C44E52",
  ifelse(full_plot, "#009E73", "#77777755")
)
png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "human_cross_model_natural_vs_external.png"
  ),
  width = 5667,
  height = 5000,
  res = 600
)
par(mar = c(6, 6, 4, 2))
plot(
  x,
  y,
  pch = 16,
  cex = ifelse(strict_plot, 1.0, 0.45),
  col = colors,
  xlab = "GSE93939 natural-resistance effect percentile",
  ylab = "GSE108094 SMA-depletion effect percentile",
  main = ""
)
abline(h = 0.5, v = 0.5, col = "grey75", lty = 3)
abline(0, 1, col = "grey70", lty = 2)
label_genes <- c(
  "KCNAB3", "LINC00665", "TSPOAP1",
  primary_panel$current_gene_symbol
)
label_rows <- plot_rows & cross_model$gene_symbol %in% label_genes
text(
  cross_model$natural_resistance_percentile[label_rows],
  cross_model$gse108094_sma_depletion_percentile[label_rows],
  labels = cross_model$gene_symbol[label_rows],
  pos = 3,
  cex = 0.65
)
legend(
  "bottomright",
  legend = c(
    "Strict splice-restoration gene",
    "Full four-model resilience direction",
    "Other ranked gene"
  ),
  col = c("#C44E52", "#009E73", "#777777"),
  pch = 16,
  bty = "n"
)
dev.off()

write_session_info(
  "results/r/integration/human_cross_model_R_sessionInfo.txt"
)

cat(
  "Cross-model human motor-neuron resilience analysis complete:",
  sum(cross_model$cross_model_ranking_eligible),
  "approved genes ranked;",
  sum(cross_model$full_directional_resilience_pattern),
  "full directional resilience patterns;",
  sum(cross_model$exploratory_natural_plus_full_cross_model),
  "also meet the exploratory natural-resistance criterion\n"
)
