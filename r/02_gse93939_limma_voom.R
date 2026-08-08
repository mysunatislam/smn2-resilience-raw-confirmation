source(file.path("r", "common.R"))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

make_design <- function(metadata, include_platform) {
  design <- cbind(
    `(Intercept)` = 1,
    OMN_vs_SC = as.integer(metadata$group == "OMN"),
    sex_M = as.integer(metadata$sex == "M"),
    age_z = zscore(metadata$age_at_death),
    postmortem_delay_z = zscore(metadata$postmortem_delay_hours),
    source_NDRI = as.integer(metadata$tissue_source == "NDRI"),
    source_NIH = as.integer(metadata$tissue_source == "NIH")
  )
  if (include_platform) {
    design <- cbind(
      design,
      platform_HiSeq2500 = as.integer(metadata$platform == "HiSeq2500")
    )
  }
  variable_columns <- vapply(
    seq_len(ncol(design)),
    function(index) {
      colnames(design)[index] == "(Intercept)" ||
        stats::sd(design[, index]) > 0
    },
    logical(1)
  )
  design <- design[, variable_columns, drop = FALSE]
  rownames(design) <- rownames(metadata)
  design
}

fit_blocked_voom <- function(
  counts,
  metadata,
  keep,
  include_platform,
  analysis_name
) {
  counts <- counts[keep, rownames(metadata), drop = FALSE]
  design <- make_design(metadata, include_platform)
  check_full_rank(design, analysis_name)

  y <- DGEList(counts = counts)
  y <- normLibSizes(y, method = "TMM")

  initial_voom <- voom(y, design, plot = FALSE)
  initial_correlation <- duplicateCorrelation(
    initial_voom,
    design,
    block = metadata$donor_id
  )
  weighted_voom <- voomWithQualityWeights(
    y,
    design,
    plot = FALSE,
    block = metadata$donor_id,
    correlation = initial_correlation$consensus.correlation
  )
  final_correlation <- duplicateCorrelation(
    weighted_voom,
    design,
    block = metadata$donor_id
  )
  fit <- lmFit(
    weighted_voom,
    design,
    block = metadata$donor_id,
    correlation = final_correlation$consensus.correlation
  )
  fit <- eBayes(fit, robust = TRUE)
  table <- topTable(
    fit,
    coef = "OMN_vs_SC",
    number = Inf,
    sort.by = "none"
  )
  table$gene_symbol <- rownames(table)
  table$n_samples <- nrow(metadata)
  table$n_donors <- length(unique(metadata$donor_id))
  table$consensus_donor_correlation <-
    final_correlation$consensus.correlation

  list(
    table = table,
    y = y,
    voom = weighted_voom,
    fit = fit,
    design = design,
    donor_correlation = final_correlation$consensus.correlation
  )
}

counts <- read_count_matrix("GSE93939")
metadata <- read_sample_metadata("GSE93939")
assert_identical_samples(counts, metadata, "GSE93939")

metadata <- metadata[metadata$group %in% c("OMN", "SC"), , drop = FALSE]
counts <- counts[, rownames(metadata), drop = FALSE]
stopifnot(
  nrow(metadata) == 32L,
  length(unique(metadata$donor_id)) == 19L,
  sum(metadata$group == "OMN") == 20L,
  sum(metadata$group == "SC") == 12L
)

raw_dge <- DGEList(counts = counts)
keep <- rowSums(cpm(raw_dge) >= 1) >= 7
stopifnot(sum(keep) > 10000L)

primary <- fit_blocked_voom(
  counts,
  metadata,
  keep,
  include_platform = TRUE,
  analysis_name = "GSE93939 adjusted OMN versus spinal"
)

same_platform_metadata <- metadata[
  metadata$platform == "HiSeq2000",
  ,
  drop = FALSE
]
stopifnot(
  nrow(same_platform_metadata) == 25L,
  length(unique(same_platform_metadata$donor_id)) == 13L
)
same_platform <- fit_blocked_voom(
  counts,
  same_platform_metadata,
  keep,
  include_platform = FALSE,
  analysis_name = "GSE93939 HiSeq2000 sensitivity"
)

primary_table <- primary$table
names(primary_table)[names(primary_table) == "logFC"] <-
  "log2_effect_omn_vs_sc"
names(primary_table)[names(primary_table) == "P.Value"] <- "p_value"
names(primary_table)[names(primary_table) == "adj.P.Val"] <- "q_value"

sensitivity_columns <- same_platform$table[, c(
  "gene_symbol",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "n_samples",
  "n_donors",
  "consensus_donor_correlation"
)]
names(sensitivity_columns) <- c(
  "gene_symbol",
  "same_platform_log2_effect",
  "same_platform_AveExpr",
  "same_platform_t",
  "same_platform_p_value",
  "same_platform_q_value",
  "same_platform_n_samples",
  "same_platform_n_donors",
  "same_platform_donor_correlation"
)

result <- merge(
  primary_table,
  sensitivity_columns,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
annotation <- read.delim(
  file.path(ROOT, "data", "processed", "HGNC_gene_annotation.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotation <- annotation[!duplicated(annotation$symbol), , drop = FALSE]
result <- merge(
  result,
  annotation,
  by.x = "gene_symbol",
  by.y = "symbol",
  all.x = TRUE,
  sort = FALSE
)
result$sex_chromosome <- result$chromosome %in% c("X", "Y")
result$direction_consistent <-
  sign(result$log2_effect_omn_vs_sc) ==
  sign(result$same_platform_log2_effect)
result <- result[order(result$q_value, -abs(result$log2_effect_omn_vs_sc)), ]

write_tsv(
  result,
  "results/r/differential_expression/GSE93939_OMN_vs_SC_limma_voom.tsv"
)

comparison <- result[
  is.finite(result$log2_effect_omn_vs_sc) &
    is.finite(result$same_platform_log2_effect),
]
effect_correlation <- cor(
  comparison$log2_effect_omn_vs_sc,
  comparison$same_platform_log2_effect,
  method = "pearson"
)
summary <- data.frame(
  analysis = c("primary_adjusted", "HiSeq2000_sensitivity"),
  samples = c(nrow(metadata), nrow(same_platform_metadata)),
  donors = c(
    length(unique(metadata$donor_id)),
    length(unique(same_platform_metadata$donor_id))
  ),
  genes_tested = c(nrow(primary$table), nrow(same_platform$table)),
  genes_fdr_005 = c(
    sum(primary$table$adj.P.Val < 0.05),
    sum(same_platform$table$adj.P.Val < 0.05)
  ),
  donor_correlation = c(
    primary$donor_correlation,
    same_platform$donor_correlation
  ),
  primary_sensitivity_effect_correlation = c(effect_correlation, NA_real_)
)
summary$model_terms <- c(
  paste(colnames(primary$design), collapse = ";"),
  paste(colnames(same_platform$design), collapse = ";")
)
write_tsv(
  summary,
  "results/r/differential_expression/GSE93939_limma_voom_summary.tsv"
)

png(
  file.path(ROOT, "results", "r", "figures", "GSE93939_R_MDS.png"),
  width = 6000,
  height = 4667,
  res = 600
)
colors <- ifelse(metadata$group == "OMN", "#C44E52", "#3C5488")
plotMDS(
  primary$y,
  col = colors,
  pch = 19,
  main = ""
)
legend(
  "topright",
  legend = c("Oculomotor", "Spinal"),
  col = c("#C44E52", "#3C5488"),
  pch = 19,
  bty = "n"
)
dev.off()

png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "GSE93939_R_primary_vs_sensitivity.png"
  ),
  width = 5333,
  height = 5000,
  res = 600
)
plot(
  comparison$log2_effect_omn_vs_sc,
  comparison$same_platform_log2_effect,
  pch = 16,
  cex = 0.45,
  col = grDevices::adjustcolor("#3C5488", alpha.f = 0.35),
  xlab = "Adjusted OMN - spinal log2 effect",
  ylab = "HiSeq2000 sensitivity log2 effect",
  main = ""
)
abline(h = 0, v = 0, col = "grey55")
abline(a = 0, b = 1, col = "#C44E52", lwd = 2)
dev.off()

cat(
  "GSE93939 R analysis complete:",
  nrow(primary$table),
  "genes;",
  sum(primary$table$adj.P.Val < 0.05),
  "at FDR < 0.05; sensitivity r =",
  sprintf("%.3f", effect_correlation),
  "\n"
)
