source(file.path("r", "common.R"))

suppressPackageStartupMessages({
  library(edgeR)
})

counts <- read_count_matrix("GSE93939")
metadata <- read_sample_metadata("GSE93939")
assert_identical_samples(counts, metadata, "GSE93939")

# Donor-level summation is used only within HiSeq2000. Summing a donor's
# libraries across two sequencing platforms would confound biological
# aggregation with platform-specific measurement.
metadata <- metadata[
  metadata$group %in% c("OMN", "SC") &
    metadata$platform == "HiSeq2000",
  ,
  drop = FALSE
]
counts <- counts[, rownames(metadata), drop = FALSE]
donor_counts <- t(rowsum(
  t(counts),
  group = metadata$donor_id,
  reorder = FALSE
))
donor_order <- colnames(donor_counts)

donor_metadata <- do.call(
  rbind,
  lapply(donor_order, function(donor_id) {
    rows <- metadata[metadata$donor_id == donor_id, , drop = FALSE]
    constant_columns <- c(
      "group",
      "sex",
      "age_at_death",
      "postmortem_delay_hours",
      "tissue_source"
    )
    stopifnot(all(vapply(
      rows[, constant_columns, drop = FALSE],
      function(values) length(unique(values)) == 1L,
      logical(1)
    )))
    data.frame(
      donor_id = donor_id,
      group = unique(rows$group),
      sex = unique(rows$sex),
      age_at_death = unique(rows$age_at_death),
      postmortem_delay_hours = unique(rows$postmortem_delay_hours),
      tissue_source = unique(rows$tissue_source),
      contributing_libraries = nrow(rows),
      row.names = donor_id,
      stringsAsFactors = FALSE
    )
  })
)
stopifnot(
  ncol(donor_counts) == 13L,
  sum(donor_metadata$group == "OMN") == 6L,
  sum(donor_metadata$group == "SC") == 7L
)

design <- cbind(
  `(Intercept)` = 1,
  OMN_vs_SC = as.integer(donor_metadata$group == "OMN"),
  sex_M = as.integer(donor_metadata$sex == "M"),
  age_z = zscore(donor_metadata$age_at_death),
  postmortem_delay_z = zscore(donor_metadata$postmortem_delay_hours),
  source_NDRI = as.integer(donor_metadata$tissue_source == "NDRI")
)
rownames(design) <- rownames(donor_metadata)
check_full_rank(design, "GSE93939 HiSeq2000 donor pseudobulk")

initial <- DGEList(counts = donor_counts)
keep <- rowSums(cpm(initial) >= 1) >= 4L
y <- DGEList(counts = donor_counts[keep, , drop = FALSE])
y <- normLibSizes(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "OMN_vs_SC")
result <- topTags(test, n = Inf, sort.by = "none")$table
result$gene_symbol <- rownames(result)
names(result)[names(result) == "logFC"] <- "log2_effect_omn_vs_sc"
names(result)[names(result) == "PValue"] <- "p_value"
names(result)[names(result) == "FDR"] <- "q_value"
result <- result[order(result$q_value, result$p_value), ]

blocked <- read.delim(
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
comparison <- merge(
  result[, c("gene_symbol", "log2_effect_omn_vs_sc")],
  blocked[, c("gene_symbol", "same_platform_log2_effect")],
  by = "gene_symbol"
)
effect_correlation <- cor(
  comparison$log2_effect_omn_vs_sc,
  comparison$same_platform_log2_effect
)

write_tsv(
  result,
  "results/r/differential_expression/GSE93939_HiSeq2000_donor_pseudobulk.tsv"
)
write_tsv(
  cbind(sample_id = rownames(donor_metadata), donor_metadata),
  "results/r/dataset_audit/GSE93939_HiSeq2000_donor_pseudobulk_samples.tsv"
)
summary <- data.frame(
  analysis = "HiSeq2000_donor_pseudobulk",
  source_libraries = nrow(metadata),
  donor_pseudobulks = nrow(donor_metadata),
  design_parameters = ncol(design),
  residual_df = fit$df.residual[1],
  genes_tested = nrow(result),
  genes_fdr_005 = sum(result$q_value < 0.05),
  blocked_pseudobulk_effect_correlation = effect_correlation,
  model_terms = paste(colnames(design), collapse = ";")
)
write_tsv(
  summary,
  "results/r/differential_expression/GSE93939_donor_pseudobulk_summary.tsv"
)

cat(
  "GSE93939 donor pseudobulk complete:",
  nrow(result),
  "genes;",
  sum(result$q_value < 0.05),
  "at FDR < 0.05; blocked-model effect r =",
  sprintf("%.3f", effect_correlation),
  "\n"
)
