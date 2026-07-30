source(file.path("r", "common.R"))

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

rank_safe_design <- function(metadata, include_platform) {
  design <- cbind(
    `(Intercept)` = 1,
    OMN_vs_SC = as.integer(metadata$group == "OMN")
  )
  nuisance <- cbind(
    sex_M = as.integer(metadata$sex == "M"),
    age_z = zscore(metadata$age_at_death),
    postmortem_delay_z = zscore(metadata$postmortem_delay_hours),
    source_NDRI = as.integer(metadata$tissue_source == "NDRI"),
    source_NIH = as.integer(metadata$tissue_source == "NIH")
  )
  if (include_platform) {
    nuisance <- cbind(
      nuisance,
      platform_HiSeq2500 = as.integer(metadata$platform == "HiSeq2500")
    )
  }

  dropped <- character()
  for (column in colnames(nuisance)) {
    candidate <- cbind(design, nuisance[, column])
    colnames(candidate)[ncol(candidate)] <- column
    if (
      is.finite(stats::sd(nuisance[, column])) &&
      stats::sd(nuisance[, column]) > 0 &&
      qr(candidate)$rank == ncol(candidate)
    ) {
      design <- candidate
    } else {
      dropped <- c(dropped, column)
    }
  }
  rownames(design) <- rownames(metadata)
  check_full_rank(design, "LODO fold")
  list(design = design, dropped = dropped)
}

fit_gse93939_fold <- function(
  counts,
  metadata,
  global_keep,
  held_out_donor,
  include_platform,
  analysis,
  fixed_donor_correlation
) {
  fold_metadata <- metadata[
    metadata$donor_id != held_out_donor,
    ,
    drop = FALSE
  ]
  fold_counts <- counts[global_keep, rownames(fold_metadata), drop = FALSE]
  nonzero <- rowSums(fold_counts) > 0
  fold_counts <- fold_counts[nonzero, , drop = FALSE]
  design_info <- rank_safe_design(fold_metadata, include_platform)
  design <- design_info$design

  y <- DGEList(counts = fold_counts)
  y <- normLibSizes(y, method = "TMM")
  # Fold-specific duplicateCorrelation and array-quality-weight optimization
  # are prohibitively slow for all-gene jackknifing. The consensus donor
  # correlation is locked from the corresponding complete-cohort model. Each
  # fold still re-estimates TMM factors, voom weights, covariate coefficients,
  # and empirical-Bayes variance moderation.
  blocked_voom <- voom(
    y,
    design,
    plot = FALSE,
    block = fold_metadata$donor_id,
    correlation = fixed_donor_correlation
  )
  fit <- lmFit(
    blocked_voom,
    design,
    block = fold_metadata$donor_id,
    correlation = fixed_donor_correlation
  )
  fit <- eBayes(fit, robust = TRUE)
  result <- topTable(
    fit,
    coef = "OMN_vs_SC",
    number = Inf,
    sort.by = "none"
  )
  data.frame(
    analysis = analysis,
    held_out_unit = held_out_donor,
    held_out_group = unique(metadata$group[metadata$donor_id == held_out_donor]),
    gene_symbol = rownames(result),
    log2_effect = result$logFC,
    p_value = result$P.Value,
    donors_in_fold = length(unique(fold_metadata$donor_id)),
    retained_terms = paste(colnames(design), collapse = ";"),
    dropped_terms = paste(design_info$dropped, collapse = ";"),
    donor_correlation = fixed_donor_correlation,
    fold_weighting = paste(
      "blocked_voom_fixed_complete_model_correlation",
      "without_array_quality_weights",
      sep = "_"
    ),
    stringsAsFactors = FALSE
  )
}

summarize_jackknife <- function(
  long,
  full,
  full_effect_column,
  expected_folds,
  prefix
) {
  split_effects <- split(long$log2_effect, long$gene_symbol)
  genes <- names(split_effects)
  full_effect <- full[[full_effect_column]][match(genes, full$gene_symbol)]
  fold_count <- lengths(split_effects)
  sign_fraction <- vapply(seq_along(genes), function(index) {
    values <- split_effects[[index]]
    mean(sign(values) == sign(full_effect[index]))
  }, numeric(1))
  output <- data.frame(
    gene_symbol = genes,
    full_log2_effect = full_effect,
    estimable_folds = fold_count,
    median_holdout_log2_effect = vapply(
      split_effects,
      median,
      numeric(1),
      na.rm = TRUE
    ),
    minimum_holdout_log2_effect = vapply(
      split_effects,
      min,
      numeric(1),
      na.rm = TRUE
    ),
    maximum_holdout_log2_effect = vapply(
      split_effects,
      max,
      numeric(1),
      na.rm = TRUE
    ),
    minimum_absolute_holdout_effect = vapply(
      split_effects,
      function(values) min(abs(values), na.rm = TRUE),
      numeric(1)
    ),
    holdout_sign_agreement_fraction = sign_fraction,
    all_folds_same_direction =
      fold_count == expected_folds & sign_fraction == 1,
    effect_range_crosses_zero = vapply(
      split_effects,
      function(values) min(values, na.rm = TRUE) <= 0 &&
        max(values, na.rm = TRUE) >= 0,
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
  names(output)[names(output) != "gene_symbol"] <- paste0(
    prefix,
    names(output)[names(output) != "gene_symbol"]
  )
  output
}

# GSE93939: leave all libraries from one donor out together.
omn_counts <- read_count_matrix("GSE93939")
omn_metadata <- read_sample_metadata("GSE93939")
assert_identical_samples(omn_counts, omn_metadata, "GSE93939")
omn_metadata <- omn_metadata[
  omn_metadata$group %in% c("OMN", "SC"),
  ,
  drop = FALSE
]
omn_counts <- omn_counts[, rownames(omn_metadata), drop = FALSE]
omn_keep <- rowSums(cpm(DGEList(counts = omn_counts)) >= 1) >= 7L
omn_donors <- unique(omn_metadata$donor_id)
stopifnot(length(omn_donors) == 19L)
omn_model_summary <- read.delim(
  file.path(
    ROOT, "results", "r", "differential_expression",
    "GSE93939_limma_voom_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
primary_fixed_correlation <- omn_model_summary$donor_correlation[
  omn_model_summary$analysis == "primary_adjusted"
]
same_fixed_correlation <- omn_model_summary$donor_correlation[
  omn_model_summary$analysis == "HiSeq2000_sensitivity"
]
stopifnot(
  length(primary_fixed_correlation) == 1L,
  is.finite(primary_fixed_correlation),
  length(same_fixed_correlation) == 1L,
  is.finite(same_fixed_correlation)
)

omn_lodo <- do.call(rbind, lapply(omn_donors, function(donor) {
  cat("GSE93939 primary LODO:", donor, "\n")
  fit_gse93939_fold(
    omn_counts,
    omn_metadata,
    omn_keep,
    donor,
    include_platform = TRUE,
    analysis = "GSE93939_primary_LODO",
    fixed_donor_correlation = primary_fixed_correlation
  )
}))

same_metadata <- omn_metadata[
  omn_metadata$platform == "HiSeq2000",
  ,
  drop = FALSE
]
same_donors <- unique(same_metadata$donor_id)
stopifnot(length(same_donors) == 13L)
same_lodo <- do.call(rbind, lapply(same_donors, function(donor) {
  cat("GSE93939 HiSeq2000 LODO:", donor, "\n")
  fit_gse93939_fold(
    omn_counts,
    same_metadata,
    omn_keep,
    donor,
    include_platform = FALSE,
    analysis = "GSE93939_HiSeq2000_LODO",
    fixed_donor_correlation = same_fixed_correlation
  )
}))

omn_full <- read.delim(
  file.path(
    ROOT, "results", "r", "differential_expression",
    "GSE93939_OMN_vs_SC_limma_voom.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
omn_lodo_summary <- summarize_jackknife(
  omn_lodo,
  omn_full,
  "log2_effect_omn_vs_sc",
  19L,
  "primary_lodo_"
)
same_lodo_summary <- summarize_jackknife(
  same_lodo,
  omn_full,
  "same_platform_log2_effect",
  13L,
  "same_platform_lodo_"
)

# GSE290979 disease: leave one complete donor line out per refit.
aggregate_counts <- function(counts, groups) {
  t(rowsum(t(counts), group = groups, reorder = FALSE))
}
organoid_counts <- read_count_matrix("GSE290979")
organoid_metadata <- read_sample_metadata("GSE290979")
assert_identical_samples(organoid_counts, organoid_metadata, "GSE290979")
organoid_metadata <- organoid_metadata[colnames(organoid_counts), , drop = FALSE]
untreated_metadata <- organoid_metadata[
  organoid_metadata$treatment == "NT",
  ,
  drop = FALSE
]
disease_counts <- aggregate_counts(
  organoid_counts[, rownames(untreated_metadata), drop = FALSE],
  untreated_metadata$`cell line`
)
disease_lines <- colnames(disease_counts)
disease_groups <- vapply(disease_lines, function(line) {
  unique(untreated_metadata$genotype[untreated_metadata$`cell line` == line])
}, character(1))
disease_keep <- rowSums(cpm(DGEList(counts = disease_counts)) >= 1) >= 2L
stopifnot(
  length(disease_lines) == 5L,
  sum(disease_groups == "CTRL") == 3L,
  sum(disease_groups == "SMA") == 2L
)

disease_lolo <- do.call(rbind, lapply(disease_lines, function(held_line) {
  keep_lines <- disease_lines != held_line
  fold_counts <- disease_counts[disease_keep, keep_lines, drop = FALSE]
  fold_metadata <- data.frame(
    genotype = factor(
      disease_groups[keep_lines],
      levels = c("CTRL", "SMA")
    ),
    row.names = disease_lines[keep_lines]
  )
  design <- model.matrix(~genotype, fold_metadata)
  check_full_rank(design, paste("GSE290979 LOLO", held_line))
  y <- DGEList(counts = fold_counts)
  y <- normLibSizes(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "genotypeSMA")
  table <- topTags(test, n = Inf, sort.by = "none")$table
  cat("GSE290979 disease LOLO:", held_line, "\n")
  data.frame(
    analysis = "GSE290979_disease_LOLO",
    held_out_unit = held_line,
    held_out_group = unname(disease_groups[held_line]),
    gene_symbol = rownames(table),
    log2_effect = table$logFC,
    p_value = table$PValue,
    donor_lines_in_fold = sum(keep_lines),
    stringsAsFactors = FALSE
  )
}))

disease_full <- read.delim(
  file.path(
    ROOT, "results", "r", "differential_expression",
    "GSE290979_SMA_vs_CTRL_pseudobulk.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
disease_lolo_summary <- summarize_jackknife(
  disease_lolo,
  disease_full,
  "sma_vs_control_log2_effect",
  5L,
  "disease_lolo_"
)

# Two treatment lines cannot support a formal LOLO refit. Report each line's
# normalized paired effect so replication is assessed without library leakage.
treatment_metadata <- organoid_metadata[
  organoid_metadata$treatment %in% c("Scramble", "R6-Mo"),
  ,
  drop = FALSE
]
treatment_groups <- paste(
  treatment_metadata$`cell line`,
  treatment_metadata$treatment,
  sep = "__"
)
treatment_counts <- aggregate_counts(
  organoid_counts[, rownames(treatment_metadata), drop = FALSE],
  treatment_groups
)
treatment_y <- DGEList(counts = treatment_counts)
treatment_y <- normLibSizes(treatment_y, method = "TMM")
treatment_log_cpm <- cpm(treatment_y, log = TRUE, prior.count = 2)
treatment_lines <- sort(unique(treatment_metadata$`cell line`))
stopifnot(length(treatment_lines) == 2L)
treatment_line_effects <- do.call(rbind, lapply(treatment_lines, function(line) {
  r6 <- paste(line, "R6-Mo", sep = "__")
  scramble <- paste(line, "Scramble", sep = "__")
  data.frame(
    analysis = "GSE290979_treatment_linewise",
    cell_line = line,
    gene_symbol = rownames(treatment_log_cpm),
    log2_effect = treatment_log_cpm[, r6] - treatment_log_cpm[, scramble],
    stringsAsFactors = FALSE
  )
}))
treatment_full <- read.delim(
  file.path(
    ROOT, "results", "r", "differential_expression",
    "GSE290979_R6_vs_scramble_pseudobulk.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
treatment_genes <- treatment_full$gene_symbol[
  is.finite(treatment_full$r6_vs_scramble_log2_effect)
]
treatment_full_effect <- treatment_full$r6_vs_scramble_log2_effect[
  match(treatment_genes, treatment_full$gene_symbol)
]
line_effect <- function(line) {
  rows <- treatment_line_effects$cell_line == line
  values <- setNames(
    treatment_line_effects$log2_effect[rows],
    treatment_line_effects$gene_symbol[rows]
  )
  unname(values[treatment_genes])
}
s2_effect <- line_effect("S2")
s3_effect <- line_effect("S3")
stopifnot(all(is.finite(s2_effect)), all(is.finite(s3_effect)))
treatment_line_summary <- data.frame(
  gene_symbol = treatment_genes,
  full_log2_effect = treatment_full_effect,
  S2_log2_effect = s2_effect,
  S3_log2_effect = s3_effect,
  line_sign_agreement_fraction = (
    as.integer(sign(s2_effect) == sign(treatment_full_effect)) +
      as.integer(sign(s3_effect) == sign(treatment_full_effect))
  ) / 2,
  both_lines_match_full_direction =
    sign(s2_effect) == sign(treatment_full_effect) &
      sign(s3_effect) == sign(treatment_full_effect),
  lines_agree_with_each_other = sign(s2_effect) == sign(s3_effect),
  formal_LOLO_estimable = FALSE,
  stringsAsFactors = FALSE
)

write_tsv(omn_lodo, "results/r/robustness/GSE93939_primary_LODO_long.tsv")
write_tsv(
  same_lodo,
  "results/r/robustness/GSE93939_HiSeq2000_LODO_long.tsv"
)
write_tsv(
  omn_lodo_summary,
  "results/r/robustness/GSE93939_primary_LODO_summary.tsv"
)
write_tsv(
  same_lodo_summary,
  "results/r/robustness/GSE93939_HiSeq2000_LODO_summary.tsv"
)
write_tsv(
  disease_lolo,
  "results/r/robustness/GSE290979_disease_LOLO_long.tsv"
)
write_tsv(
  disease_lolo_summary,
  "results/r/robustness/GSE290979_disease_LOLO_summary.tsv"
)
write_tsv(
  treatment_line_effects,
  "results/r/robustness/GSE290979_treatment_linewise_long.tsv"
)
write_tsv(
  treatment_line_summary,
  "results/r/robustness/GSE290979_treatment_linewise_summary.tsv"
)

policy <- read.delim(
  file.path(ROOT, "config", "resampling_policy.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(policy, "results/r/robustness/resampling_policy.tsv")
summary <- data.frame(
  analysis = c(
    "GSE93939_primary_LODO", "GSE93939_HiSeq2000_LODO",
    "GSE290979_disease_LOLO", "GSE290979_treatment_linewise"
  ),
  biological_units = c(19L, 13L, 5L, 2L),
  folds_or_line_effects = c(19L, 13L, 5L, 2L),
  formal_holdout_estimable = c(TRUE, TRUE, TRUE, FALSE),
  genes_summarized = c(
    nrow(omn_lodo_summary), nrow(same_lodo_summary),
    nrow(disease_lolo_summary), nrow(treatment_line_summary)
  ),
  all_unit_directions_match_full = c(
    sum(omn_lodo_summary$primary_lodo_all_folds_same_direction),
    sum(same_lodo_summary$same_platform_lodo_all_folds_same_direction),
    sum(disease_lolo_summary$disease_lolo_all_folds_same_direction),
    sum(treatment_line_summary$both_lines_match_full_direction)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "results/r/robustness/biological_unit_holdout_summary.tsv")
write_session_info("results/r/robustness/holdout_R_sessionInfo.txt")

cat(
  "Biological-unit holdout analysis complete:",
  "19 GSE93939 donor folds, 13 same-platform donor folds,",
  "5 GSE290979 disease line folds, and 2 treatment line effects\n"
)
