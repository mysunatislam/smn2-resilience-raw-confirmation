source(file.path("r", "common.R"))

cross_model <- read.delim(
  file.path(
    ROOT, "results", "r", "integration",
    "human_cross_model_resilience_gene_table.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
atlas <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE243076_C20_motor_neuron_context.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
deg <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_cell_type_pseudobulk_DEGs.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
mn_deg <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_motor_neuron_pseudobulk_DEGs.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
mn_markers <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_motor_neuron_marker_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  !anyDuplicated(cross_model$gene_symbol),
  !anyDuplicated(atlas$current_gene_symbol),
  !anyDuplicated(mn_deg$current_gene_symbol),
  !anyDuplicated(mn_markers$current_gene_symbol)
)

atlas_columns <- c(
  "current_gene_symbol", "source_gene_symbol", "C20",
  "c20_other_neuron_mean", "c20_rank_of_21",
  "c20_log2_enrichment_pc_0_001", "c20_log2_enrichment_pc_0_01",
  "c20_log2_enrichment_pc_0_1", "c20_expression_percentile",
  "c20_is_max_neuronal_cluster", "c20_robust_localization"
)
atlas_context <- atlas[atlas_columns]
names(atlas_context)[names(atlas_context) == "current_gene_symbol"] <- "gene_symbol"
names(atlas_context)[names(atlas_context) != "gene_symbol"] <- paste0(
  "gse243076_",
  names(atlas_context)[names(atlas_context) != "gene_symbol"]
)

mn_columns <- c(
  "current_gene_symbol", "source_gene_symbol", "p_value", "q_value",
  "published_log2fc_ctrl_vs_sma", "sma_vs_control_log2_effect"
)
mn_context <- mn_deg[mn_columns]
names(mn_context)[names(mn_context) == "current_gene_symbol"] <- "gene_symbol"
names(mn_context)[names(mn_context) != "gene_symbol"] <- paste0(
  "gse290980_mn_",
  names(mn_context)[names(mn_context) != "gene_symbol"]
)

marker_context <- mn_markers
names(marker_context)[names(marker_context) == "current_gene_symbol"] <- "gene_symbol"

cell_type_split <- split(deg, deg$current_gene_symbol)
cell_type_context <- do.call(rbind, lapply(cell_type_split, function(frame) {
  data.frame(
    gene_symbol = frame$current_gene_symbol[1],
    gse290980_significant_cell_type_count = length(unique(frame$cell_type)),
    gse290980_significant_cell_types = paste(
      sort(unique(frame$cell_type)),
      collapse = ","
    ),
    gse290980_cell_types_sma_higher = sum(frame$sma_vs_control_log2_effect > 0),
    gse290980_cell_types_sma_lower = sum(frame$sma_vs_control_log2_effect < 0),
    stringsAsFactors = FALSE
  )
}))
rownames(cell_type_context) <- NULL

integrated <- Reduce(function(left, right) {
  merge(left, right, by = "gene_symbol", all.x = TRUE, sort = FALSE)
}, list(
  cross_model,
  atlas_context,
  mn_context,
  marker_context,
  cell_type_context
))
stopifnot(nrow(integrated) == nrow(cross_model))

finite <- function(values) !is.na(values) & is.finite(values)
natural <- integrated$omn_log2_effect_omn_vs_sc
mn_disease <- integrated$gse290980_mn_sma_vs_control_log2_effect
integrated$gse290980_mn_reported_deg <- finite(mn_disease)
integrated$gse290980_mn_sma_depleted <- finite(mn_disease) & mn_disease < 0
integrated$gse290980_mn_opposes_natural_resistance <-
  finite(natural) & finite(mn_disease) & sign(natural) == -sign(mn_disease)
integrated$gse243076_adult_mn_detected <-
  !is.na(integrated$gse243076_C20) & integrated$gse243076_C20 > 0
integrated$gse243076_adult_mn_localized <-
  !is.na(integrated$gse243076_c20_robust_localization) &
  integrated$gse243076_c20_robust_localization
integrated$cell_resolved_support_tier <- "no_cell_resolved_support"
integrated$cell_resolved_support_tier[
  integrated$gse243076_adult_mn_detected
] <- "independent_adult_MN_detected"
integrated$cell_resolved_support_tier[
  integrated$gse290980_mn_opposes_natural_resistance
] <- "same_study_MN_disease_direction"
integrated$cell_resolved_support_tier[
  integrated$gse290980_mn_opposes_natural_resistance &
    integrated$gse243076_adult_mn_detected
] <- "same_study_MN_direction_plus_independent_adult_MN_detection"
integrated$cell_resolved_support_tier[
  integrated$gse290980_mn_opposes_natural_resistance &
    integrated$gse243076_adult_mn_localized
] <- "same_study_MN_direction_plus_independent_adult_MN_localization"

primary_panel <- read.delim(
  file.path(
    ROOT, "results", "r", "integration",
    "human_cross_model_primary_validation_panel.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
bridge <- read.delim(
  file.path(
    ROOT, "results", "r", "integration",
    "human_cross_model_exploratory_bridge.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
integrated$in_primary_splicing_panel <-
  integrated$gene_symbol %in% primary_panel$current_gene_symbol
integrated$in_exploratory_splicing_bridge <-
  integrated$gene_symbol %in% bridge$current_gene_symbol

write_tsv(
  integrated,
  "results/r/cell_resolved/human_cross_model_cell_resolved_gene_table.tsv"
)

shortlist <- integrated[
  !is.na(integrated$exploratory_natural_plus_full_cross_model) &
    integrated$exploratory_natural_plus_full_cross_model,
  ,
  drop = FALSE
]
shortlist <- shortlist[order(
  -as.integer(shortlist$gse290980_mn_opposes_natural_resistance),
  -as.integer(shortlist$gse243076_adult_mn_localized),
  shortlist$cross_model_resilience_rank
), ]
panel_context <- integrated[
  match(primary_panel$current_gene_symbol, integrated$gene_symbol),
  ,
  drop = FALSE
]
bridge_context <- integrated[
  match(bridge$current_gene_symbol, integrated$gene_symbol),
  ,
  drop = FALSE
]
stopifnot(
  nrow(shortlist) == 37L,
  nrow(panel_context) == 12L,
  nrow(bridge_context) == 3L,
  !anyNA(panel_context$gene_symbol),
  !anyNA(bridge_context$gene_symbol)
)
write_tsv(
  shortlist,
  "results/r/cell_resolved/human_cross_model_cell_resolved_candidates.tsv"
)
write_tsv(
  panel_context,
  "results/r/cell_resolved/human_primary_panel_cell_resolved_context.tsv"
)
write_tsv(
  bridge_context,
  "results/r/cell_resolved/human_bridge_cell_resolved_context.tsv"
)

summary <- data.frame(
  metric = c(
    "cross_model_ranked_genes",
    "ranked_genes_in_GSE290980_MN_DEG_list",
    "exploratory_candidates",
    "candidates_in_GSE290980_MN_DEG_list",
    "candidates_with_MN_disease_opposition",
    "candidates_detected_in_independent_adult_MN",
    "candidates_with_independent_adult_MN_localization",
    "candidates_with_both_cell_resolved_supports",
    "primary_panel_genes_with_MN_disease_opposition",
    "primary_panel_genes_with_independent_adult_MN_localization",
    "statistical_interpretation"
  ),
  value = c(
    sum(integrated$cross_model_ranking_eligible),
    sum(integrated$cross_model_ranking_eligible & integrated$gse290980_mn_reported_deg),
    nrow(shortlist),
    sum(shortlist$gse290980_mn_reported_deg),
    sum(shortlist$gse290980_mn_opposes_natural_resistance),
    sum(shortlist$gse243076_adult_mn_detected),
    sum(shortlist$gse243076_adult_mn_localized),
    sum(
      shortlist$gse290980_mn_opposes_natural_resistance &
        shortlist$gse243076_adult_mn_localized
    ),
    sum(panel_context$gse290980_mn_opposes_natural_resistance),
    sum(panel_context$gse243076_adult_mn_localized),
    paste(
      "GSE290980 is same-study pseudobulk context; GSE243076 is independent",
      "adult localization without an SMA contrast; neither increases donor n"
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  summary,
  "results/r/cell_resolved/human_cell_resolved_summary.tsv"
)

plot_rows <- shortlist$gse290980_mn_reported_deg &
  !is.na(shortlist$gse243076_c20_log2_enrichment_pc_0_01)
png(
  file.path(
    ROOT, "results", "r", "figures",
    "human_cell_resolved_candidate_context.png"
  ),
  width = 1750,
  height = 1450,
  res = 180
)
par(mar = c(6, 6, 4, 2))
plot(
  shortlist$gse243076_c20_log2_enrichment_pc_0_01[plot_rows],
  shortlist$gse290980_mn_sma_vs_control_log2_effect[plot_rows],
  pch = 16,
  col = ifelse(
    shortlist$gse290980_mn_opposes_natural_resistance[plot_rows],
    "#167D8D",
    "#B44B4B"
  ),
  xlab = "Independent adult C20 log2 enrichment",
  ylab = "Same-study MN SMA minus control log2 effect",
  main = "Cell-resolved context for cross-model candidates"
)
abline(h = 0, v = 0, col = "grey70", lty = 2)
label_rows <- plot_rows
if (any(label_rows)) {
  text(
    shortlist$gse243076_c20_log2_enrichment_pc_0_01[label_rows],
    shortlist$gse290980_mn_sma_vs_control_log2_effect[label_rows],
    labels = shortlist$gene_symbol[label_rows],
    pos = 3,
    cex = 0.7
  )
}
legend(
  "topright",
  legend = c("Opposes natural-resistance direction", "Other MN DEG direction"),
  col = c("#167D8D", "#B44B4B"),
  pch = 16,
  bty = "n"
)
dev.off()

write_session_info(
  "results/r/cell_resolved/cell_resolved_integration_R_sessionInfo.txt"
)
cat(
  "Cell-resolved integration complete:",
  sum(shortlist$gse290980_mn_opposes_natural_resistance),
  "of 37 candidates have same-study MN disease opposition and",
  sum(shortlist$gse243076_adult_mn_localized),
  "show independent adult C20 localization\n"
)
