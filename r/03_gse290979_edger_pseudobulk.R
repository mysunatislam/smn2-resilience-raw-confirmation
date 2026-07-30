source(file.path("r", "common.R"))

suppressPackageStartupMessages({
  library(edgeR)
})

aggregate_counts <- function(counts, groups) {
  aggregated <- rowsum(
    t(counts),
    group = groups,
    reorder = FALSE
  )
  t(aggregated)
}

fit_ql <- function(
  counts,
  design,
  coefficient,
  min_libraries,
  analysis_name
) {
  check_full_rank(design, analysis_name)
  initial <- DGEList(counts = counts)
  keep <- rowSums(cpm(initial) >= 1) >= min_libraries
  y <- DGEList(counts = counts[keep, , drop = FALSE])
  y <- normLibSizes(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = coefficient)
  table <- topTags(test, n = Inf, sort.by = "none")$table
  table$gene_symbol <- rownames(table)
  table$analysis <- analysis_name
  list(
    table = table,
    y = y,
    fit = fit,
    design = design,
    coefficient = coefficient
  )
}

counts <- read_count_matrix("GSE290979")
metadata <- read_sample_metadata("GSE290979")
assert_identical_samples(counts, metadata, "GSE290979")
metadata <- metadata[colnames(counts), , drop = FALSE]

# Disease inference uses one summed pseudobulk per donor line. Genotype and
# donor line are perfectly confounded, so individual organoid libraries cannot
# create additional independent disease donors.
untreated_metadata <- metadata[metadata$treatment == "NT", , drop = FALSE]
untreated_counts <- counts[, rownames(untreated_metadata), drop = FALSE]
disease_groups <- untreated_metadata$`cell line`
disease_counts <- aggregate_counts(untreated_counts, disease_groups)
disease_line_order <- colnames(disease_counts)
disease_metadata <- do.call(
  rbind,
  lapply(disease_line_order, function(cell_line) {
    rows <- untreated_metadata[
      untreated_metadata$`cell line` == cell_line,
      ,
      drop = FALSE
    ]
    stopifnot(length(unique(rows$genotype)) == 1L)
    data.frame(
      cell_line = cell_line,
      genotype = unique(rows$genotype),
      contributing_libraries = nrow(rows),
      row.names = cell_line,
      stringsAsFactors = FALSE
    )
  })
)
disease_metadata$genotype <- factor(
  disease_metadata$genotype,
  levels = c("CTRL", "SMA")
)
disease_design <- model.matrix(~genotype, disease_metadata)
stopifnot(
  ncol(disease_counts) == 5L,
  sum(disease_metadata$genotype == "CTRL") == 3L,
  sum(disease_metadata$genotype == "SMA") == 2L
)
disease_fit <- fit_ql(
  disease_counts,
  disease_design,
  coefficient = "genotypeSMA",
  min_libraries = 2L,
  analysis_name = "SMA_minus_control_donor_line_pseudobulk"
)

# Treatment inference is paired by SMA donor line. Four preparation libraries
# are summed within each line-condition cell, yielding two biological pairs.
treatment_metadata <- metadata[
  metadata$treatment %in% c("Scramble", "R6-Mo"),
  ,
  drop = FALSE
]
treatment_counts <- counts[, rownames(treatment_metadata), drop = FALSE]
treatment_groups <- paste(
  treatment_metadata$`cell line`,
  treatment_metadata$treatment,
  sep = "__"
)
treatment_pseudobulk <- aggregate_counts(
  treatment_counts,
  treatment_groups
)
treatment_group_order <- colnames(treatment_pseudobulk)
treatment_pseudobulk_metadata <- do.call(
  rbind,
  lapply(treatment_group_order, function(group_id) {
    pieces <- strsplit(group_id, "__", fixed = TRUE)[[1]]
    rows <- treatment_metadata[
      treatment_metadata$`cell line` == pieces[1] &
        treatment_metadata$treatment == pieces[2],
      ,
      drop = FALSE
    ]
    data.frame(
      cell_line = pieces[1],
      treatment = pieces[2],
      contributing_libraries = nrow(rows),
      row.names = group_id,
      stringsAsFactors = FALSE
    )
  })
)
treatment_pseudobulk_metadata$cell_line <- factor(
  treatment_pseudobulk_metadata$cell_line
)
treatment_pseudobulk_metadata$treatment <- factor(
  treatment_pseudobulk_metadata$treatment,
  levels = c("Scramble", "R6-Mo")
)
treatment_design <- model.matrix(
  ~cell_line + treatment,
  treatment_pseudobulk_metadata
)
treatment_coefficient <- grep(
  "^treatment",
  colnames(treatment_design),
  value = TRUE
)
stopifnot(
  ncol(treatment_pseudobulk) == 4L,
  length(treatment_coefficient) == 1L,
  all(treatment_pseudobulk_metadata$contributing_libraries == 4L)
)
treatment_fit <- fit_ql(
  treatment_pseudobulk,
  treatment_design,
  coefficient = treatment_coefficient,
  min_libraries = 2L,
  analysis_name = "R6_minus_scramble_paired_line_pseudobulk"
)

# A library-level treatment model is retained only as a sensitivity analysis.
# It estimates preparation-level consistency, not additional donor replication.
library_treatment_metadata <- treatment_metadata
library_treatment_metadata$`cell line` <- factor(
  library_treatment_metadata$`cell line`
)
library_treatment_metadata$treatment <- factor(
  library_treatment_metadata$treatment,
  levels = c("Scramble", "R6-Mo")
)
library_treatment_design <- model.matrix(
  ~`cell line` + treatment,
  library_treatment_metadata
)
library_treatment_coefficient <- grep(
  "^treatment",
  colnames(library_treatment_design),
  value = TRUE
)
library_treatment_fit <- fit_ql(
  treatment_counts,
  library_treatment_design,
  coefficient = library_treatment_coefficient,
  min_libraries = 4L,
  analysis_name = "R6_minus_scramble_library_level_sensitivity"
)

disease_result <- disease_fit$table
names(disease_result)[names(disease_result) == "logFC"] <-
  "sma_vs_control_log2_effect"
names(disease_result)[names(disease_result) == "PValue"] <- "p_value"
names(disease_result)[names(disease_result) == "FDR"] <- "q_value"
disease_result <- disease_result[order(
  disease_result$q_value,
  disease_result$p_value
), ]

treatment_result <- treatment_fit$table
names(treatment_result)[names(treatment_result) == "logFC"] <-
  "r6_vs_scramble_log2_effect"
names(treatment_result)[names(treatment_result) == "PValue"] <- "p_value"
names(treatment_result)[names(treatment_result) == "FDR"] <- "q_value"

library_sensitivity <- library_treatment_fit$table[, c(
  "gene_symbol",
  "logFC",
  "PValue",
  "FDR"
)]
names(library_sensitivity) <- c(
  "gene_symbol",
  "library_sensitivity_log2_effect",
  "library_sensitivity_p_value",
  "library_sensitivity_q_value"
)
treatment_result <- merge(
  treatment_result,
  library_sensitivity,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
treatment_result <- treatment_result[order(
  treatment_result$q_value,
  treatment_result$p_value
), ]

write_tsv(
  disease_result,
  "results/r/differential_expression/GSE290979_SMA_vs_CTRL_pseudobulk.tsv"
)
write_tsv(
  treatment_result,
  "results/r/differential_expression/GSE290979_R6_vs_scramble_pseudobulk.tsv"
)
write_tsv(
  cbind(
    sample_id = rownames(disease_metadata),
    disease_metadata
  ),
  "results/r/dataset_audit/GSE290979_disease_pseudobulk_samples.tsv"
)
write_tsv(
  cbind(
    sample_id = rownames(treatment_pseudobulk_metadata),
    treatment_pseudobulk_metadata
  ),
  "results/r/dataset_audit/GSE290979_treatment_pseudobulk_samples.tsv"
)

shared_treatment <- merge(
  treatment_fit$table[, c("gene_symbol", "logFC")],
  library_treatment_fit$table[, c("gene_symbol", "logFC")],
  by = "gene_symbol",
  suffixes = c("_pseudobulk", "_library")
)
treatment_effect_correlation <- cor(
  shared_treatment$logFC_pseudobulk,
  shared_treatment$logFC_library
)
summary <- data.frame(
  analysis = c(
    "SMA_minus_control_donor_line_pseudobulk",
    "R6_minus_scramble_paired_line_pseudobulk",
    "R6_minus_scramble_library_level_sensitivity"
  ),
  libraries_in_source = c(15L, 16L, 16L),
  modeled_samples = c(5L, 4L, 16L),
  independent_donor_lines = c(5L, 2L, 2L),
  genes_tested = c(
    nrow(disease_fit$table),
    nrow(treatment_fit$table),
    nrow(library_treatment_fit$table)
  ),
  genes_fdr_005 = c(
    sum(disease_fit$table$FDR < 0.05),
    sum(treatment_fit$table$FDR < 0.05),
    sum(library_treatment_fit$table$FDR < 0.05)
  ),
  residual_df = c(
    disease_fit$fit$df.residual[1],
    treatment_fit$fit$df.residual[1],
    library_treatment_fit$fit$df.residual[1]
  ),
  pseudobulk_library_effect_correlation = c(
    NA_real_,
    treatment_effect_correlation,
    NA_real_
  )
)
write_tsv(
  summary,
  "results/r/differential_expression/GSE290979_edger_summary.tsv"
)

png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "GSE290979_R_disease_pseudobulk_MDS.png"
  ),
  width = 1600,
  height = 1400,
  res = 180
)
colors <- ifelse(
  disease_metadata$genotype == "SMA",
  "#C44E52",
  "#3C5488"
)
plotMDS(
  disease_fit$y,
  col = colors,
  pch = 19,
  labels = rownames(disease_metadata),
  main = "GSE290979 donor-line pseudobulks"
)
legend(
  x = -2.8,
  y = -1.45,
  legend = c("Control", "SMA"),
  col = c("#3C5488", "#C44E52"),
  pch = 19,
  bty = "n"
)
dev.off()

cat(
  "GSE290979 R pseudobulk analysis complete:",
  sum(disease_fit$table$FDR < 0.05),
  "disease genes and",
  sum(treatment_fit$table$FDR < 0.05),
  "treatment genes at FDR < 0.05; treatment sensitivity r =",
  sprintf("%.3f", treatment_effect_correlation),
  "\n"
)
