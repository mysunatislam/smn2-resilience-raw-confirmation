source(file.path("r", "common.R"))

checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(passed)) {
    stop(name, " failed: observed ", observed, ", expected ", expected)
  }
}

dataset_audit <- read.delim(
  file.path(ROOT, "results", "r", "dataset_audit", "core_dataset_audit.tsv"),
  stringsAsFactors = FALSE
)
add_check(
  "core_library_counts",
  identical(dataset_audit$libraries, c(39L, 31L)),
  paste(dataset_audit$libraries, collapse = ","),
  "39,31"
)
add_check(
  "core_independent_units",
  identical(dataset_audit$primary_independent_units, c(19L, 5L)),
  paste(dataset_audit$primary_independent_units, collapse = ","),
  "19,5"
)

omn_summary <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE93939_limma_voom_summary.tsv"
  ),
  stringsAsFactors = FALSE
)
add_check(
  "GSE93939_genes_tested",
  omn_summary$genes_tested[1] == 16443L,
  omn_summary$genes_tested[1],
  16443L
)
add_check(
  "GSE93939_primary_FDR_hits",
  omn_summary$genes_fdr_005[1] == 2L,
  omn_summary$genes_fdr_005[1],
  2L
)

organoid_summary <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "differential_expression",
    "GSE290979_edger_summary.tsv"
  ),
  stringsAsFactors = FALSE
)
add_check(
  "GSE290979_primary_FDR_hits",
  all(organoid_summary$genes_fdr_005[1:2] == 0L),
  paste(organoid_summary$genes_fdr_005[1:2], collapse = ","),
  "0,0"
)

event_summary <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_event_summary.tsv"
  ),
  stringsAsFactors = FALSE
)
add_check(
  "splicing_exact_common",
  sum(event_summary$exact_common) == 182L,
  sum(event_summary$exact_common),
  182L
)
add_check(
  "splicing_pooled_corrected",
  sum(event_summary$pooled_corrected) == 108L,
  sum(event_summary$pooled_corrected),
  108L
)
add_check(
  "splicing_strict_events",
  sum(event_summary$strict_both_lines) == 83L,
  sum(event_summary$strict_both_lines),
  83L
)

integration_summary <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "integration",
    "human_R_two_track_summary.tsv"
  ),
  stringsAsFactors = FALSE
)
summary_value <- setNames(
  integration_summary$value,
  integration_summary$metric
)
add_check(
  "strict_splicing_genes",
  as.integer(summary_value[["R_strict_splicing_genes"]]) == 74L,
  summary_value[["R_strict_splicing_genes"]],
  74L
)
add_check(
  "robust_OMN_bridge_genes",
  as.integer(summary_value[["R_strict_plus_robust_OMN_FDR05"]]) == 0L,
  summary_value[["R_strict_plus_robust_OMN_FDR05"]],
  0L
)
add_check(
  "exploratory_OMN_bridge_genes",
  as.integer(summary_value[["R_strict_plus_exploratory_OMN_p05"]]) == 3L,
  summary_value[["R_strict_plus_exploratory_OMN_p05"]],
  3L
)

method_comparison <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "integration",
    "R_vs_Python_method_comparison.tsv"
  ),
  stringsAsFactors = FALSE
)
event_match <- method_comparison[
  method_comparison$comparison == "strict_splicing_event_pair_identity",
  "exact_set_match"
]
add_check(
  "R_Python_strict_event_identity",
  identical(event_match, TRUE),
  event_match,
  TRUE
)

raw_sample_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_raw_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE290979_raw_libraries",
  nrow(raw_sample_sheet) == 31L,
  nrow(raw_sample_sheet),
  31L
)
add_check(
  "GSE290979_raw_read_length",
  identical(sort(unique(raw_sample_sheet$read_length)), 151L),
  paste(sort(unique(raw_sample_sheet$read_length)), collapse = ","),
  151L
)
add_check(
  "GSE290979_raw_read_pairs",
  sum(raw_sample_sheet$read_pairs) == 2377570407,
  format(sum(raw_sample_sheet$read_pairs), scientific = FALSE),
  2377570407
)

raw_list_directory <- file.path(ROOT, "config", "rmats", "GSE290979")
raw_list_lengths <- vapply(
  c(
    "disease_group1_sma", "disease_group2_control",
    "treatment_group1_r6", "treatment_group2_scramble"
  ),
  function(name) length(readLines(
    file.path(raw_list_directory, paste0(name, ".txt")),
    warn = FALSE
  )),
  integer(1)
)
add_check(
  "GSE290979_rmats_contrast_order",
  identical(unname(raw_list_lengths), c(6L, 9L, 8L, 8L)),
  paste(raw_list_lengths, collapse = ","),
  "6,9,8,8"
)

fixed_manifest <- read.delim(
  file.path(raw_list_directory, "fixed_events_manifest.tsv"),
  stringsAsFactors = FALSE
)
add_check(
  "GSE290979_fixed_event_classes",
  identical(unname(fixed_manifest$events), c(35L, 6L, 5L, 22L, 15L)),
  paste(fixed_manifest$events, collapse = ","),
  "35,6,5,22,15"
)

validation_ranking <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_validation_event_ranking.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "validation_ranking_dimensions",
  nrow(validation_ranking) == 83L &&
    length(unique(validation_ranking$current_gene_symbol)) == 74L,
  paste(
    nrow(validation_ranking),
    length(unique(validation_ranking$current_gene_symbol)),
    sep = ","
  ),
  "83,74"
)
ranking_is_ordered <-
  identical(validation_ranking$validation_rank, seq_len(83L)) &&
  all(diff(validation_ranking$validation_priority_score) <= 1e-12)
scores_are_bounded <-
  all(is.finite(validation_ranking$validation_priority_score)) &&
  all(validation_ranking$validation_priority_score >= 0) &&
  all(validation_ranking$validation_priority_score <= 1)
add_check(
  "validation_scores_bounded_and_ranked",
  ranking_is_ordered && scores_are_bounded,
  paste(ranking_is_ordered, scores_are_bounded, sep = ","),
  "TRUE,TRUE"
)

validation_panel <- read.delim(
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
panel_design_pass <-
  nrow(validation_panel) == 12L &&
  length(unique(validation_panel$current_gene_symbol)) == 12L &&
  setequal(validation_panel$event_type, c("ES", "A5SS", "A3SS", "MXE", "RI"))
add_check(
  "primary_validation_panel_design",
  panel_design_pass,
  paste(
    nrow(validation_panel),
    length(unique(validation_panel$current_gene_symbol)),
    length(unique(validation_panel$event_type)),
    sep = ","
  ),
  "12,12,5"
)
add_check(
  "primary_validation_panel_HGNC_status",
  all(validation_panel$hgnc_status == "Approved"),
  paste(sort(unique(validation_panel$hgnc_status)), collapse = ","),
  "Approved"
)

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
add_check(
  "exploratory_OMN_bridge_panel",
  nrow(bridge_panel) == 3L && setequal(
    bridge_panel$current_gene_symbol,
    c("KCNAB3", "LINC00665", "TSPOAP1")
  ),
  paste(sort(bridge_panel$current_gene_symbol), collapse = ","),
  "KCNAB3,LINC00665,TSPOAP1"
)

support_exception <- validation_ranking$current_gene_symbol[
  !validation_ranking$junction_support_exact
]
add_check(
  "junction_support_source_exception",
  identical(support_exception, "ECD"),
  paste(support_exception, collapse = ","),
  "ECD"
)

score_specification <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "splicing",
    "GSE290979_R_validation_score_specification.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "validation_score_weights",
  nrow(score_specification) == 5L &&
    isTRUE(all.equal(sum(score_specification$weight), 1)),
  paste(nrow(score_specification), sum(score_specification$weight), sep = ","),
  "5,1"
)

selection <- read.delim(
  file.path(ROOT, "config", "dataset_selection.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
external_selection <- selection[
  selection$decision == "EXTERNAL_VALIDATION",
  ,
  drop = FALSE
]
add_check(
  "external_validation_dataset_registry",
  nrow(external_selection) == 1L &&
    identical(external_selection$accession, "GSE108094") &&
    identical(
      external_selection$role,
      "independent_ipsc_sma_motor_neuron_validation"
    ),
  paste(external_selection$accession, external_selection$role, sep = ","),
  "GSE108094,independent_ipsc_sma_motor_neuron_validation"
)

external_summary <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "external_validation",
    "GSE108094_audit_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
external_value <- setNames(external_summary$value, external_summary$metric)
add_check(
  "GSE108094_external_design",
  identical(
    as.integer(external_value[c(
      "libraries", "independent_lines", "control_lines", "sma_lines"
    )]),
    c(8L, 4L, 2L, 2L)
  ),
  paste(external_value[c(
    "libraries", "independent_lines", "control_lines", "sma_lines"
  )], collapse = ","),
  "8,4,2,2"
)
add_check(
  "GSE108094_expression_audit",
  identical(
    as.integer(external_value[c(
      "official_expression_rows", "official_expression_OK_rows",
      "official_expression_FDR05_rows"
    )]),
    c(63657L, 24691L, 1892L)
  ),
  paste(external_value[c(
    "official_expression_rows", "official_expression_OK_rows",
    "official_expression_FDR05_rows"
  )], collapse = ","),
  "63657,24691,1892"
)
add_check(
  "GSE108094_splicing_audit",
  identical(
    as.integer(external_value[c(
      "external_ES_rows", "external_ES_FDR05_abs_dPSI_gt_005"
    )]),
    c(29237L, 1947L)
  ),
  paste(external_value[c(
    "external_ES_rows", "external_ES_FDR05_abs_dPSI_gt_005"
  )], collapse = ","),
  "29237,1947"
)

external_manifest <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "external_validation",
    "GSE108094_resource_manifest.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE108094_resource_hashes",
  nrow(external_manifest) == 3L && all(external_manifest$sha256_pass),
  paste(nrow(external_manifest), sum(external_manifest$sha256_pass), sep = ","),
  "3,3"
)

cross_summary <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "integration",
    "human_cross_model_resilience_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cross_value <- setNames(cross_summary$value, cross_summary$metric)
add_check(
  "cross_model_ranked_universe",
  as.integer(cross_value[["approved_genes_ranked"]]) == 11326L,
  cross_value[["approved_genes_ranked"]],
  11326L
)
add_check(
  "cross_model_direction_patterns",
  identical(
    as.integer(cross_value[c(
      "full_directional_resilience_pattern",
      "full_directional_vulnerability_pattern",
      "exploratory_natural_plus_full_cross_model"
    )]),
    c(626L, 468L, 37L)
  ),
  paste(cross_value[c(
    "full_directional_resilience_pattern",
    "full_directional_vulnerability_pattern",
    "exploratory_natural_plus_full_cross_model"
  )], collapse = ","),
  "626,468,37"
)
add_check(
  "cross_model_splicing_coverage",
  identical(
    as.integer(cross_value[c(
      "strict_splicing_genes_with_external_direction",
      "primary_panel_genes_with_external_direction",
      "bridge_genes_with_external_direction"
    )]),
    c(65L, 12L, 3L)
  ),
  paste(cross_value[c(
    "strict_splicing_genes_with_external_direction",
    "primary_panel_genes_with_external_direction",
    "bridge_genes_with_external_direction"
  )], collapse = ","),
  "65,12,3"
)

cross_shortlist <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "integration",
    "human_cross_model_exploratory_resilience_candidates.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
shortlist_strict <- cross_shortlist$gene_symbol[
  !is.na(cross_shortlist$splicing_strict_event_count) &
    cross_shortlist$splicing_strict_event_count > 0
]
add_check(
  "cross_model_unique_strict_shortlist_gene",
  nrow(cross_shortlist) == 37L &&
    identical(shortlist_strict, "LINC00665"),
  paste(nrow(cross_shortlist), paste(shortlist_strict, collapse = ","), sep = ","),
  "37,LINC00665"
)

cross_correlations <- read.delim(
  file.path(
    ROOT,
    "results",
    "r",
    "integration",
    "human_cross_model_effect_correlations.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "cross_model_global_effect_correlations",
  nrow(cross_correlations) == 5L &&
    max(abs(cross_correlations$spearman_rho)) < 0.12,
  paste(
    nrow(cross_correlations),
    sprintf("%.6f", max(abs(cross_correlations$spearman_rho))),
    sep = ","
  ),
  "5,<0.12"
)

selection <- read.delim(
  file.path(ROOT, "config", "dataset_selection.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cell_decisions <- setNames(selection$decision, selection$accession)
add_check(
  "cell_resolved_dataset_roles",
  identical(
    unname(cell_decisions[c("GSE290980", "GSE243076", "EGAD00001011259")]),
    c(
      "CELL_RESOLVED_INTERNAL", "INDEPENDENT_CELL_LOCALIZATION",
      "CONTROLLED_ACCESS_FUTURE_VALIDATION"
    )
  ),
  paste(
    cell_decisions[c("GSE290980", "GSE243076", "EGAD00001011259")],
    collapse = ","
  ),
  paste(
    c(
      "CELL_RESOLVED_INTERNAL", "INDEPENDENT_CELL_LOCALIZATION",
      "CONTROLLED_ACCESS_FUTURE_VALIDATION"
    ),
    collapse = ","
  )
)

cell_resource_registry <- read.delim(
  file.path(ROOT, "config", "cell_resolved_resources.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "cell_resolved_resource_registry",
  nrow(cell_resource_registry) == 7L &&
    !anyDuplicated(cell_resource_registry[c("dataset", "resource")]) &&
    all(grepl("^[0-9a-f]{64}$", cell_resource_registry$sha256)) &&
    identical(
      as.integer(table(factor(
        cell_resource_registry$dataset,
        levels = c("GSE243076", "GSE290980")
      ))),
      c(3L, 4L)
    ),
  paste(
    nrow(cell_resource_registry),
    paste(table(cell_resource_registry$dataset), collapse = ","),
    sep = ","
  ),
  "7,3,4"
)

atlas_resources <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE243076_resource_manifest.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE243076_resource_hashes",
  nrow(atlas_resources) == 3L && all(atlas_resources$sha256_pass),
  paste(nrow(atlas_resources), sum(atlas_resources$sha256_pass), sep = ","),
  "3,3"
)

atlas_raw_qc <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE243076_raw_10x_donor_qc.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE243076_raw_10x_dimensions",
  nrow(atlas_raw_qc) == 9L &&
    sum(atlas_raw_qc$filtered_nuclei) == 74711L &&
    sum(atlas_raw_qc$nuclei_passing_threshold_only_qc_proxy) == 72834L &&
    all(atlas_raw_qc$features == 36601L),
  paste(
    nrow(atlas_raw_qc), sum(atlas_raw_qc$filtered_nuclei),
    sum(atlas_raw_qc$nuclei_passing_threshold_only_qc_proxy),
    paste(unique(atlas_raw_qc$features), collapse = ","),
    sep = ","
  ),
  "9,74711,72834,36601"
)

atlas_summary <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved", "GSE243076_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
atlas_value <- setNames(atlas_summary$value, atlas_summary$metric)
add_check(
  "GSE243076_author_atlas_dimensions",
  identical(
    as.integer(atlas_value[c(
      "independent_donors", "author_reported_post_qc_nuclei",
      "features", "neuronal_clusters"
    )]),
    c(9L, 64021L, 36601L, 21L)
  ) && atlas_value[["author_motor_neuron_cluster"]] == "C20",
  paste(atlas_value[c(
    "independent_donors", "author_reported_post_qc_nuclei",
    "features", "neuronal_clusters", "author_motor_neuron_cluster"
  )], collapse = ","),
  "9,64021,36601,21,C20"
)

atlas_markers <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE243076_C20_marker_validation.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE243076_C20_marker_identity",
  nrow(atlas_markers) == 8L &&
    all(atlas_markers$C20 > 0) &&
    sum(atlas_markers$c20_robust_localization) == 4L &&
    all(c("CHAT", "SLC5A7", "MNX1", "ISL1") %in%
      atlas_markers$current_gene_symbol[atlas_markers$c20_robust_localization]),
  paste(
    nrow(atlas_markers), sum(atlas_markers$c20_robust_localization),
    collapse = ","
  ),
  "8,4"
)

single_cell_resources <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_resource_manifest.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE290980_resource_hashes",
  nrow(single_cell_resources) == 4L && all(single_cell_resources$sha256_pass),
  paste(
    nrow(single_cell_resources), sum(single_cell_resources$sha256_pass),
    sep = ","
  ),
  "4,4"
)

single_cell_samples <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_sample_metadata.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE290980_replication_units",
  nrow(single_cell_samples) == 8L &&
    length(unique(single_cell_samples$donor_line)) == 4L &&
    identical(
      as.integer(table(factor(
        single_cell_samples$genotype,
        levels = c("CTRL", "SMA")
      ))),
      c(4L, 4L)
    ),
  paste(
    nrow(single_cell_samples), length(unique(single_cell_samples$donor_line)),
    paste(table(single_cell_samples$genotype), collapse = ","),
    sep = ","
  ),
  "8,4,4,4"
)

single_cell_clusters <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_cluster_annotation.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
single_cell_markers <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_cluster_markers.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE290980_cluster_marker_structure",
  nrow(single_cell_clusters) == 19L &&
    sum(single_cell_clusters$summary_cell_type == "MN") == 4L &&
    nrow(single_cell_markers) == 5700L,
  paste(
    nrow(single_cell_clusters),
    sum(single_cell_clusters$summary_cell_type == "MN"),
    nrow(single_cell_markers),
    sep = ","
  ),
  "19,4,5700"
)

single_cell_deg <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "GSE290980_cell_type_DEG_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_cell_degs <- c(2031L, 834L, 385L, 303L, 141L, 128L, 121L)
add_check(
  "GSE290980_cell_type_DEG_counts",
  identical(single_cell_deg$significant_degs, expected_cell_degs),
  paste(single_cell_deg$significant_degs, collapse = ","),
  paste(expected_cell_degs, collapse = ",")
)

cell_summary <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "human_cell_resolved_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cell_value <- setNames(cell_summary$value, cell_summary$metric)
cell_metrics <- c(
  "cross_model_ranked_genes",
  "ranked_genes_in_GSE290980_MN_DEG_list",
  "exploratory_candidates",
  "candidates_in_GSE290980_MN_DEG_list",
  "candidates_with_MN_disease_opposition",
  "candidates_detected_in_independent_adult_MN",
  "candidates_with_independent_adult_MN_localization",
  "candidates_with_both_cell_resolved_supports",
  "primary_panel_genes_with_MN_disease_opposition",
  "primary_panel_genes_with_independent_adult_MN_localization"
)
expected_cell_values <- c(11326L, 1385L, 37L, 5L, 2L, 35L, 0L, 0L, 2L, 0L)
add_check(
  "cell_resolved_integration_counts",
  identical(as.integer(cell_value[cell_metrics]), expected_cell_values),
  paste(cell_value[cell_metrics], collapse = ","),
  paste(expected_cell_values, collapse = ",")
)

cell_candidates <- read.delim(
  file.path(
    ROOT, "results", "r", "cell_resolved",
    "human_cross_model_cell_resolved_candidates.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
opposition_genes <- sort(cell_candidates$gene_symbol[
  cell_candidates$gse290980_mn_opposes_natural_resistance
])
add_check(
  "cell_resolved_candidate_identity",
  nrow(cell_candidates) == 37L &&
    identical(opposition_genes, c("PFKFB2", "SHQ1")),
  paste(nrow(cell_candidates), paste(opposition_genes, collapse = ","), sep = ","),
  "37,PFKFB2,SHQ1"
)

gse69175_resources <- read.delim(
  file.path(
    ROOT, "results", "r", "external_validation",
    "GSE69175_resource_manifest.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "GSE69175_resource_hashes",
  nrow(gse69175_resources) == 2L && all(gse69175_resources$sha256_pass),
  paste(nrow(gse69175_resources), sum(gse69175_resources$sha256_pass), sep = ","),
  "2,2"
)

gse69175_summary <- read.delim(
  file.path(
    ROOT, "results", "r", "external_validation",
    "GSE69175_audit_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gse69175_value <- setNames(gse69175_summary$value, gse69175_summary$metric)
gse69175_metrics <- c(
  "libraries", "independent_lines", "control_lines", "sma_lines",
  "official_expression_rows", "official_expression_OK_rows",
  "collapsed_approved_genes", "direction_usable_genes",
  "official_expression_FDR05_rows"
)
expected_gse69175 <- c(4L, 2L, 1L, 1L, 19286L, 5394L, 17430L, 5093L, 968L)
add_check(
  "GSE69175_audit_counts",
  identical(as.integer(gse69175_value[gse69175_metrics]), expected_gse69175),
  paste(gse69175_value[gse69175_metrics], collapse = ","),
  paste(expected_gse69175, collapse = ",")
)
add_check(
  "GSE69175_LOLO_not_estimable",
  identical(gse69175_value[["LOLO_estimable"]], "FALSE"),
  gse69175_value[["LOLO_estimable"]],
  "FALSE"
)

symbol_examples <- normalize_excel_gene_symbol(c(
  "1-Mar", "11-Sep", "Mar-6", "Sep-14", "KCNAB3"
))
add_check(
  "legacy_month_gene_symbol_repair",
  identical(
    symbol_examples,
    c("MARCHF1", "SEPTIN11", "MARCHF6", "SEPTIN14", "KCNAB3")
  ),
  paste(symbol_examples, collapse = ","),
  "MARCHF1,SEPTIN11,MARCHF6,SEPTIN14,KCNAB3"
)

primary_lodo <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "GSE93939_primary_LODO_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
same_lodo <- read.delim(
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
add_check(
  "biological_unit_holdout_fold_completeness",
  nrow(primary_lodo) == 16443L &&
    all(primary_lodo$primary_lodo_estimable_folds == 19L) &&
    nrow(same_lodo) == 16443L &&
    all(same_lodo$same_platform_lodo_estimable_folds == 13L) &&
    nrow(disease_lolo) == 17014L &&
    all(disease_lolo$disease_lolo_estimable_folds == 5L),
  paste(
    nrow(primary_lodo), unique(primary_lodo$primary_lodo_estimable_folds),
    nrow(same_lodo), unique(same_lodo$same_platform_lodo_estimable_folds),
    nrow(disease_lolo), unique(disease_lolo$disease_lolo_estimable_folds),
    sep = ","
  ),
  "16443,19,16443,13,17014,5"
)

holdout_summary <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "human_holdout_validation_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
holdout_value <- setNames(holdout_summary$value, holdout_summary$metric)
holdout_metrics <- c(
  "ranked_genes", "leave_one_dataset_out_folds",
  "exploratory_candidates",
  "candidates_primary_GSE93939_LODO_robust",
  "candidates_GSE290979_disease_LOLO_robust",
  "candidates_GSE290979_treatment_both_lines_robust",
  "candidates_all_biological_unit_checks_robust",
  "candidates_GSE69175_direction_usable",
  "candidates_GSE69175_opposes_natural",
  "candidates_unit_robust_plus_GSE69175_direction"
)
expected_holdout <- c(11326L, 3L, 37L, 35L, 15L, 14L, 7L, 9L, 6L, 0L)
add_check(
  "holdout_integration_counts",
  identical(as.integer(holdout_value[holdout_metrics]), expected_holdout),
  paste(holdout_value[holdout_metrics], collapse = ","),
  paste(expected_holdout, collapse = ",")
)

holdout_candidates <- read.delim(
  file.path(
    ROOT, "results", "r", "robustness",
    "human_cross_model_holdout_candidates.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
unit_robust_genes <- sort(holdout_candidates$gene_symbol[
  holdout_candidates$all_estimable_biological_unit_checks_robust
])
gse69175_support_genes <- sort(holdout_candidates$gene_symbol[
  holdout_candidates$gse69175_disease_opposes_natural
])
add_check(
  "holdout_candidate_identity",
  identical(
    unit_robust_genes,
    sort(c("LY6H", "HS3ST5", "ZNF853", "PNCK", "IL17D", "CLPTM1", "PDPR"))
  ) && identical(
    gse69175_support_genes,
    sort(c("ROBO2", "ATP9B", "HEY1", "LRRC4", "SPIN1", "EIF4ENIF1"))
  ),
  paste(
    paste(unit_robust_genes, collapse = ";"),
    paste(gse69175_support_genes, collapse = ";"),
    sep = ","
  ),
  paste(
    paste(sort(c("LY6H", "HS3ST5", "ZNF853", "PNCK", "IL17D", "CLPTM1", "PDPR")), collapse = ";"),
    paste(sort(c("ROBO2", "ATP9B", "HEY1", "LRRC4", "SPIN1", "EIF4ENIF1")), collapse = ";"),
    sep = ","
  )
)

publication_summary <- read.delim(
  file.path(
    ROOT, "results", "r", "publication",
    "publication_analysis_summary.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
publication_value <- setNames(
  publication_summary$value,
  publication_summary$metric
)
publication_count_metrics <- c(
  "frozen_candidate_count",
  "raw_libraries",
  "raw_integrity_checks_passed",
  "raw_candidate_genes_mapped",
  "raw_full_directional_pattern",
  "processed_unit_robust_candidates",
  "raw_unit_robust_candidates",
  "raw_and_processed_unit_robust_candidates",
  "raw_or_processed_unit_robust_candidates",
  "tier_1_candidates",
  "tier_2_candidates",
  "tier_3_candidates",
  "tier_4_candidates",
  "tier_5_candidates"
)
expected_publication_counts <- c(
  37L, 9L, 13L, 35L, 21L, 7L, 9L, 4L, 12L,
  4L, 5L, 3L, 10L, 15L
)
add_check(
  "raw_confirmed_publication_counts",
  identical(
    as.integer(publication_value[publication_count_metrics]),
    expected_publication_counts
  ),
  paste(publication_value[publication_count_metrics], collapse = ","),
  paste(expected_publication_counts, collapse = ",")
)
add_check(
  "raw_confirmed_publication_scope",
  identical(
    unname(publication_value[c(
      "raw_quantification_method",
      "whole_genome_alignment",
      "random_split_used",
      "raw_confirmation_independent_validation"
    )]),
    c(
      "salmon_transcriptome_quasimapping",
      "FALSE", "FALSE", "FALSE"
    )
  ),
  paste(publication_value[c(
    "raw_quantification_method",
    "whole_genome_alignment",
    "random_split_used",
    "raw_confirmation_independent_validation"
  )], collapse = ","),
  "salmon_transcriptome_quasimapping,FALSE,FALSE,FALSE"
)

publication_candidates <- read.delim(
  file.path(
    ROOT, "results", "r", "publication",
    "human_raw_confirmed_publication_candidates.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
publication_tier_1 <- sort(publication_candidates$gene_symbol[
  publication_candidates$publication_evidence_tier ==
    "tier_1_raw_and_processed_unit_robust"
])
publication_tier_2 <- sort(publication_candidates$gene_symbol[
  publication_candidates$publication_evidence_tier ==
    "tier_2_raw_unit_robust"
])
publication_tier_3 <- sort(publication_candidates$gene_symbol[
  publication_candidates$publication_evidence_tier ==
    "tier_3_processed_unit_robust"
])
add_check(
  "raw_confirmed_publication_candidate_identity",
  identical(
    publication_tier_1,
    sort(c("LY6H", "HS3ST5", "ZNF853", "IL17D"))
  ) && identical(
    publication_tier_2,
    sort(c("LSAMP", "FAAH", "C12orf60", "KCNJ4", "CABLES2"))
  ) && identical(
    publication_tier_3,
    sort(c("PNCK", "CLPTM1", "PDPR"))
  ),
  paste(
    paste(publication_tier_1, collapse = ";"),
    paste(publication_tier_2, collapse = ";"),
    paste(publication_tier_3, collapse = ";"),
    sep = ","
  ),
  paste(
    paste(sort(c("LY6H", "HS3ST5", "ZNF853", "IL17D")), collapse = ";"),
    paste(sort(c(
      "LSAMP", "FAAH", "C12orf60", "KCNJ4", "CABLES2"
    )), collapse = ";"),
    paste(sort(c("PNCK", "CLPTM1", "PDPR")), collapse = ";"),
    sep = ","
  )
)
holdout_rank_match <- match(
  publication_candidates$gene_symbol,
  holdout_candidates$gene_symbol
)
add_check(
  "publication_preserves_frozen_discovery_rank",
  nrow(publication_candidates) == 37L &&
    identical(
      publication_candidates$publication_priority_order,
      seq_len(37L)
    ) &&
    all(!is.na(holdout_rank_match)) &&
    identical(
      publication_candidates$frozen_discovery_rank_unchanged,
      holdout_candidates$cross_model_resilience_rank[holdout_rank_match]
    ) &&
    all(!publication_candidates$raw_confirmation_is_independent_validation),
  paste(
    nrow(publication_candidates),
    all(!is.na(holdout_rank_match)),
    all(!publication_candidates$raw_confirmation_is_independent_validation),
    sep = ","
  ),
  "37,TRUE,TRUE"
)

publication_concordance <- read.delim(
  file.path(
    ROOT, "results", "r", "publication",
    "raw_vs_processed_concordance.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_concordance <- data.frame(
  contrast = c("SMA_vs_CTRL", "R6-Mo_vs_Scramble"),
  matched_genes = c(14636L, 14385L),
  spearman_effect_correlation = c(
    0.778292840776912,
    0.587905715385299
  ),
  direction_agreement_fraction = c(
    0.80657283410768,
    0.732082029892249
  ),
  stringsAsFactors = FALSE
)
concordance_match <- match(
  expected_concordance$contrast,
  publication_concordance$contrast
)
add_check(
  "raw_vs_processed_publication_concordance",
  all(!is.na(concordance_match)) &&
    identical(
      publication_concordance$matched_genes[concordance_match],
      expected_concordance$matched_genes
    ) &&
    isTRUE(all.equal(
      publication_concordance$spearman_effect_correlation[
        concordance_match
      ],
      expected_concordance$spearman_effect_correlation,
      tolerance = 1e-12
    )) &&
    isTRUE(all.equal(
      publication_concordance$direction_agreement_fraction[
        concordance_match
      ],
      expected_concordance$direction_agreement_fraction,
      tolerance = 1e-12
    )),
  paste(
    publication_concordance$matched_genes[concordance_match],
    sprintf(
      "%.6f",
      publication_concordance$spearman_effect_correlation[
        concordance_match
      ]
    ),
    sep = ":",
    collapse = ","
  ),
  "14636:0.778293,14385:0.587906"
)

publication_files <- c(
  file.path(
    ROOT, "manuscript", "supplementary_table_S1_candidates.tsv"
  ),
  file.path(
    ROOT, "manuscript", "supplementary_table_S2_summary.tsv"
  ),
  file.path(
    ROOT, "manuscript", "supplementary_table_S3_concordance.tsv"
  ),
  file.path(
    ROOT, "manuscript", "figures",
    "figure_1_raw_vs_processed_concordance.png"
  ),
  file.path(
    ROOT, "manuscript", "figures",
    "figure_2_candidate_evidence_matrix.png"
  )
)
publication_files_exist <- all(file.exists(publication_files))
publication_files_nonempty <- publication_files_exist &&
  all(file.info(publication_files)$size > 0)
add_check(
  "publication_artifacts_exist",
  publication_files_nonempty,
  paste(
    basename(publication_files),
    file.exists(publication_files),
    sep = ":",
    collapse = ","
  ),
  "all present and nonempty"
)

policy <- read.delim(
  file.path(ROOT, "config", "resampling_policy.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
add_check(
  "resampling_policy_constraints",
  policy$status[policy$accession == "GSE69175"] ==
    "STRUCTURALLY_NOT_ESTIMABLE" &&
    policy$status[
      policy$accession == "GSE290979" &
        policy$analysis == "R6_vs_scramble_treatment"
    ] == "NOT_FORMALLY_ESTIMABLE",
  paste(policy$status[policy$accession %in% c("GSE69175", "GSE290979")], collapse = ","),
  "GSE69175 structurally unavailable; GSE290979 treatment not formally estimable"
)

analytic_files <- c(
  setdiff(
    list.files(file.path(ROOT, "r"), pattern = "[.]R$", full.names = TRUE),
    file.path(ROOT, "r", "06_validate_outputs.R")
  ),
  list.files(file.path(ROOT, "scripts"), pattern = "[.]py$", full.names = TRUE)
)
random_split_pattern <- paste(
  "(^|[^A-Za-z0-9_.])set[.]seed[[:space:]]*[(]|",
  "(^|[^A-Za-z0-9_])sample[[:space:]]*[(]|",
  "initial_split|train_test_split|",
  "createDataPartition|vfold_cv|random_state|np[.]random",
  sep = ""
)
random_split_hits <- unlist(lapply(analytic_files, function(path) {
  lines <- readLines(path, warn = FALSE)
  matches <- grep(random_split_pattern, lines, ignore.case = TRUE, value = TRUE)
  if (length(matches)) paste(basename(path), matches, sep = ":") else character()
}))
add_check(
  "no_random_analytic_splits",
  length(random_split_hits) == 0L,
  paste(random_split_hits, collapse = ";"),
  "none"
)

r_scripts <- list.files(
  file.path(ROOT, "r"),
  pattern = "[.]R$",
  full.names = TRUE
)
parse_pass <- all(vapply(r_scripts, function(path) {
  !inherits(try(parse(file = path), silent = TRUE), "try-error")
}, logical(1)))
add_check(
  "R_scripts_parse",
  parse_pass,
  sum(vapply(r_scripts, function(path) {
    !inherits(try(parse(file = path), silent = TRUE), "try-error")
  }, logical(1))),
  length(r_scripts)
)

package_versions <- data.frame(
  package = c("R", "edgeR", "limma", "readxl", "digest", "Matrix"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("edgeR")),
    as.character(packageVersion("limma")),
    as.character(packageVersion("readxl")),
    as.character(packageVersion("digest")),
    as.character(packageVersion("Matrix"))
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  do.call(rbind, checks),
  "results/r/R_output_validation.tsv"
)
write_tsv(
  package_versions,
  "results/r/R_package_versions.tsv"
)
write_session_info("results/r/R_sessionInfo.txt")

cat(length(checks), "R output checks passed\n")
