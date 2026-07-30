source(file.path("r", "common.R"))

hgnc <- read.delim(
  file.path(ROOT, "data", "metadata", "hgnc_complete_set.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)
resolver <- build_hgnc_resolver(hgnc)
stopifnot(
  identical(unname(resolver$aliases[["GAS8"]]), "DRC4"),
  identical(unname(resolver$aliases[["ENTHD2"]]), "TEPSIN"),
  identical(unname(resolver$aliases[["ASUN"]]), "INTS13")
)

omn <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE93939_OMN_vs_SC_limma_voom.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
disease <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE290979_SMA_vs_CTRL_pseudobulk.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
treatment <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE290979_R6_vs_scramble_pseudobulk.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
splicing <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_splicing_gene_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

omn_selected <- omn[, c(
  "gene_symbol",
  "log2_effect_omn_vs_sc",
  "p_value",
  "q_value",
  "same_platform_log2_effect",
  "same_platform_p_value",
  "same_platform_q_value",
  "direction_consistent"
)]
omn_mapped <- remap_unique(
  omn_selected,
  "gene_symbol",
  resolver,
  "omn_"
)

disease_selected <- disease[, c(
  "gene_symbol",
  "sma_vs_control_log2_effect",
  "p_value",
  "q_value"
)]
disease_mapped <- remap_unique(
  disease_selected,
  "gene_symbol",
  resolver,
  "disease_"
)

treatment_selected <- treatment[, c(
  "gene_symbol",
  "r6_vs_scramble_log2_effect",
  "p_value",
  "q_value",
  "library_sensitivity_log2_effect",
  "library_sensitivity_q_value"
)]
treatment_mapped <- remap_unique(
  treatment_selected,
  "gene_symbol",
  resolver,
  "treatment_"
)

splicing$source_symbol <- splicing$gene_symbol
splicing$gene_symbol <- resolve_symbols(splicing$gene_symbol, resolver)
stopifnot(!anyDuplicated(splicing$gene_symbol))
names(splicing)[names(splicing) != "gene_symbol"] <- paste0(
  "splicing_",
  names(splicing)[names(splicing) != "gene_symbol"]
)

integrated <- Reduce(
  function(left, right) merge(
    left,
    right,
    by = "gene_symbol",
    all = TRUE,
    sort = FALSE
  ),
  list(omn_mapped, disease_mapped, treatment_mapped, splicing)
)

annotation <- hgnc[!duplicated(hgnc$symbol), c(
  "symbol",
  "name",
  "locus_group",
  "locus_type",
  "location",
  "status"
)]
names(annotation) <- c(
  "gene_symbol",
  "hgnc_name",
  "hgnc_locus_group",
  "hgnc_locus_type",
  "hgnc_location",
  "hgnc_status"
)
integrated <- merge(
  integrated,
  annotation,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

integrated$strict_splicing <-
  !is.na(integrated$splicing_strict_event_count) &
  integrated$splicing_strict_event_count > 0
integrated$robust_omn_positive_fdr05 <-
  !is.na(integrated$omn_log2_effect_omn_vs_sc) &
  integrated$omn_log2_effect_omn_vs_sc > 0 &
  integrated$omn_q_value < 0.05 &
  integrated$omn_same_platform_log2_effect > 0 &
  integrated$omn_same_platform_q_value < 0.05
integrated$exploratory_omn_positive_p05 <-
  !is.na(integrated$omn_log2_effect_omn_vs_sc) &
  integrated$omn_log2_effect_omn_vs_sc > 0 &
  integrated$omn_p_value < 0.05 &
  integrated$omn_same_platform_log2_effect > 0 &
  integrated$omn_same_platform_p_value < 0.10
integrated$sma_expression_depleted <-
  !is.na(integrated$disease_sma_vs_control_log2_effect) &
  integrated$disease_sma_vs_control_log2_effect < 0
integrated$r6_expression_increased <-
  !is.na(integrated$treatment_r6_vs_scramble_log2_effect) &
  integrated$treatment_r6_vs_scramble_log2_effect > 0

integrated$R_evidence_class <- "context_only"
integrated$R_evidence_class[integrated$strict_splicing] <-
  "strict_splicing"
integrated$R_evidence_class[
  integrated$strict_splicing &
    integrated$exploratory_omn_positive_p05
] <- "strict_splicing_plus_exploratory_OMN"
integrated$R_evidence_class[
  integrated$strict_splicing &
    integrated$robust_omn_positive_fdr05
] <- "strict_splicing_plus_robust_OMN"
integrated$independent_validation_required <- TRUE

integrated <- integrated[order(
  -as.integer(integrated$strict_splicing),
  -ifelse(
    is.na(integrated$splicing_strict_event_count),
    0,
    integrated$splicing_strict_event_count
  ),
  -as.integer(integrated$exploratory_omn_positive_p05),
  ifelse(is.na(integrated$omn_p_value), 1, integrated$omn_p_value)
), ]

strict <- integrated[integrated$strict_splicing, , drop = FALSE]
exploratory_bridge <- strict[
  strict$exploratory_omn_positive_p05,
  ,
  drop = FALSE
]
robust_bridge <- strict[
  strict$robust_omn_positive_fdr05,
  ,
  drop = FALSE
]

write_tsv(
  integrated,
  "results/r/integration/human_R_two_track_gene_table.tsv"
)
write_tsv(
  strict,
  "results/r/integration/human_R_strict_splicing_candidates.tsv"
)
write_tsv(
  exploratory_bridge,
  "results/r/integration/human_R_exploratory_bridge_candidates.tsv"
)

python_omn <- read.delim(
  file.path(
    ROOT,
    "results",
    "differential_expression",
    "GSE93939_OMN_vs_SC_donor_aware.tsv"
  ),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
python_omn$gene_symbol <- rownames(python_omn)
omn_comparison <- merge(
  omn[, c("gene_symbol", "log2_effect_omn_vs_sc")],
  python_omn[, c("gene_symbol", "log2_effect_omn_vs_sc")],
  by = "gene_symbol",
  suffixes = c("_R", "_Python")
)

python_disease <- read.delim(
  file.path(
    ROOT,
    "results",
    "differential_expression",
    "GSE290979_SMA_vs_CTRL.tsv"
  ),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
python_disease$gene_symbol <- rownames(python_disease)
disease_comparison <- merge(
  disease[, c("gene_symbol", "sma_vs_control_log2_effect")],
  python_disease[, c("gene_symbol", "sma_vs_ctrl_log2_effect")],
  by = "gene_symbol"
)

python_treatment <- read.delim(
  file.path(
    ROOT,
    "results",
    "differential_expression",
    "GSE290979_R6Mo_vs_Scramble.tsv"
  ),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
python_treatment$gene_symbol <- rownames(python_treatment)
treatment_comparison <- merge(
  treatment[, c("gene_symbol", "r6_vs_scramble_log2_effect")],
  python_treatment[, c("gene_symbol", "r6mo_vs_scramble_log2_effect")],
  by = "gene_symbol"
)

r_strict_events <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_strict_corrected_events.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
python_strict_events <- read.delim(
  file.path(
    ROOT,
    "results",
    "splicing",
    "GSE290979_strict_corrected_rescue_events.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
r_event_pairs <- paste(
  r_strict_events$event_type,
  r_strict_events$disease_event_id,
  r_strict_events$treatment_event_id,
  sep = "|"
)
python_event_pairs <- paste(
  python_strict_events$event_type,
  python_strict_events$disease_event_id,
  python_strict_events$treatment_event_id,
  sep = "|"
)
stopifnot(setequal(r_event_pairs, python_event_pairs))

method_comparison <- data.frame(
  comparison = c(
    "GSE93939_R_vs_Python_effect",
    "GSE290979_disease_R_pseudobulk_vs_Python_effect",
    "GSE290979_treatment_R_pseudobulk_vs_Python_effect",
    "strict_splicing_event_pair_identity"
  ),
  shared_features = c(
    nrow(omn_comparison),
    nrow(disease_comparison),
    nrow(treatment_comparison),
    length(r_event_pairs)
  ),
  pearson_correlation = c(
    cor(
      omn_comparison$log2_effect_omn_vs_sc_R,
      omn_comparison$log2_effect_omn_vs_sc_Python
    ),
    cor(
      disease_comparison$sma_vs_control_log2_effect,
      disease_comparison$sma_vs_ctrl_log2_effect
    ),
    cor(
      treatment_comparison$r6_vs_scramble_log2_effect,
      treatment_comparison$r6mo_vs_scramble_log2_effect
    ),
    NA_real_
  ),
  exact_set_match = c(NA, NA, NA, TRUE)
)
write_tsv(
  method_comparison,
  "results/r/integration/R_vs_Python_method_comparison.tsv"
)

summary <- data.frame(
  metric = c(
    "R_OMN_positive_genes_FDR05",
    "R_strict_splicing_genes",
    "R_strict_plus_robust_OMN_FDR05",
    "R_strict_plus_exploratory_OMN_p05",
    "R_strict_plus_exploratory_OMN_p05_genes"
  ),
  value = c(
    sum(integrated$robust_omn_positive_fdr05, na.rm = TRUE),
    nrow(strict),
    nrow(robust_bridge),
    nrow(exploratory_bridge),
    paste(sort(exploratory_bridge$gene_symbol), collapse = ";")
  )
)
write_tsv(
  summary,
  "results/r/integration/human_R_two_track_summary.tsv"
)

png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "R_vs_Python_GSE93939_effects.png"
  ),
  width = 1600,
  height = 1500,
  res = 180
)
plot(
  omn_comparison$log2_effect_omn_vs_sc_Python,
  omn_comparison$log2_effect_omn_vs_sc_R,
  pch = 16,
  cex = 0.45,
  col = grDevices::adjustcolor("#3C5488", alpha.f = 0.35),
  xlab = "Earlier Python adjusted effect",
  ylab = "R limma-voom blocked effect",
  main = sprintf(
    "GSE93939 effect agreement (r = %.3f)",
    method_comparison$pearson_correlation[1]
  )
)
abline(a = 0, b = 1, col = "#C44E52", lwd = 2)
abline(h = 0, v = 0, col = "grey55")
dev.off()

cat(
  "R integration complete:",
  nrow(strict),
  "strict splicing genes;",
  nrow(robust_bridge),
  "robust and",
  nrow(exploratory_bridge),
  "exploratory OMN bridge genes; strict R/Python events identical\n"
)
