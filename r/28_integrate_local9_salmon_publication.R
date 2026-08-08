script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) != 1L) {
  stop("This script must be run with Rscript")
}
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))
source(file.path(ROOT, "r", "common.R"))

parse_options <- function(arguments) {
  parsed <- list()
  for (argument in arguments) {
    pieces <- regmatches(
      argument,
      regexec("^--([A-Za-z0-9-]+)=(.*)$", argument)
    )[[1L]]
    if (length(pieces) != 3L) {
      stop("Expected --name=value argument, received: ", argument)
    }
    parsed[[pieces[2L]]] <- pieces[3L]
  }
  parsed
}

option_value <- function(options, name, default = NULL) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) default else value
}

read_required <- function(path) {
  if (!file.exists(path)) {
    stop("Required input is missing: ", path)
  }
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

require_columns <- function(frame, columns, label) {
  missing <- setdiff(columns, names(frame))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

as_flag <- function(values) {
  if (is.logical(values)) {
    return(values)
  }
  text <- toupper(trimws(as.character(values)))
  result <- rep(NA, length(text))
  result[text %in% c("TRUE", "T", "1")] <- TRUE
  result[text %in% c("FALSE", "F", "0")] <- FALSE
  result
}

flag_or_false <- function(values) {
  result <- as_flag(values)
  result[is.na(result)] <- FALSE
  result
}

write_output <- function(frame, path) {
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

options <- parse_options(commandArgs(trailingOnly = TRUE))
raw_root <- option_value(
  options,
  "raw-root",
  file.path(
    ROOT, "results", "r", "raw_confirmation", "local9_salmon"
  )
)
output_root <- option_value(
  options,
  "output-root",
  file.path(ROOT, "results", "r", "publication")
)
manuscript_root <- option_value(
  options,
  "manuscript-root",
  file.path(ROOT, "manuscript")
)
expected_candidates <- suppressWarnings(as.integer(option_value(
  options,
  "expected-candidates",
  "37"
)))
if (is.na(expected_candidates) || expected_candidates < 1L) {
  stop("--expected-candidates must be a positive integer")
}

completion <- read_required(file.path(
  raw_root, "LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv"
))
validation <- read_required(file.path(
  raw_root, "GSE290979_local9_raw_analysis_validation.tsv"
))
concordance <- read_required(file.path(
  raw_root, "GSE290979_local9_raw_vs_processed_concordance_summary.tsv"
))
candidates <- read_required(file.path(
  raw_root, "GSE290979_local9_raw_candidate_confirmation.tsv"
))
disease <- read_required(file.path(
  raw_root, "GSE290979_local9_raw_vs_processed_disease.tsv"
))
treatment <- read_required(file.path(
  raw_root, "GSE290979_local9_raw_vs_processed_treatment.tsv"
))

require_columns(completion, c("metric", "value"), "Completion marker")
completion_value <- setNames(completion$value, completion$metric)
required_completion <- c(
  status = "COMPLETE",
  libraries = "9",
  candidate_count = as.character(expected_candidates),
  random_split_used = "FALSE",
  inference_role =
    "raw_read_sensitivity_confirmation_not_independent_validation",
  quantification_method = "salmon_transcriptome_quasimapping",
  whole_genome_alignment = "FALSE"
)
observed_completion <- completion_value[names(required_completion)]
if (
  any(is.na(observed_completion)) ||
    !identical(unname(observed_completion), unname(required_completion))
) {
  stop(
    "Completion marker does not satisfy the publication gate. Observed: ",
    paste(
      names(required_completion),
      observed_completion,
      sep = "=",
      collapse = "; "
    )
  )
}

require_columns(validation, c("check", "passed"), "Validation table")
validation_flags <- as_flag(validation$passed)
if (
  nrow(validation) != 13L ||
    any(is.na(validation_flags)) ||
    !all(validation_flags)
) {
  stop(
    "Expected exactly 13 passing local9 integrity checks; observed ",
    sum(validation_flags %in% TRUE), " of ", nrow(validation)
  )
}

require_columns(
  concordance,
  c(
    "contrast", "matched_genes", "spearman_effect_correlation",
    "pearson_effect_correlation", "direction_agreement_fraction"
  ),
  "Concordance summary"
)
expected_contrasts <- c("SMA_vs_CTRL", "R6-Mo_vs_Scramble")
if (
  !setequal(concordance$contrast, expected_contrasts) ||
    any(!is.finite(concordance$matched_genes)) ||
    any(!is.finite(concordance$spearman_effect_correlation)) ||
    any(!is.finite(concordance$pearson_effect_correlation)) ||
    any(!is.finite(concordance$direction_agreement_fraction))
) {
  stop("Concordance summary is incomplete or non-finite")
}

candidate_columns <- c(
  "gene_symbol",
  "cross_model_resilience_rank",
  "exploratory_cross_model_rank",
  "cross_model_resilience_score",
  "omn_log2_effect_omn_vs_sc",
  "omn_p_value",
  "omn_q_value",
  "disease_sma_vs_control_log2_effect",
  "treatment_r6_vs_scramble_log2_effect",
  "gse108094_sma_vs_control_log2_effect",
  "strict_splicing",
  "all_estimable_biological_unit_checks_robust",
  "gse69175_disease_opposes_natural",
  "gse290980_mn_opposes_natural_resistance",
  "gse243076_adult_mn_detected",
  "gse243076_adult_mn_localized",
  "local9_disease_raw_gene_id",
  "local9_disease_raw_log2_effect",
  "local9_disease_raw_p_value",
  "local9_disease_raw_q_value",
  "local9_disease_raw_lolo_folds_tested",
  "local9_disease_raw_lolo_all_same_direction",
  "local9_treatment_raw_gene_id",
  "local9_treatment_raw_log2_effect",
  "local9_treatment_raw_p_value",
  "local9_treatment_raw_q_value",
  "local9_treatment_raw_treatment_lines",
  "local9_treatment_raw_treatment_lines_same_direction",
  "local9_raw_disease_opposes_natural",
  "local9_raw_r6_aligns_natural",
  "local9_raw_disease_lolo_robust",
  "local9_raw_treatment_line_robust",
  "local9_raw_full_directional_pattern",
  "local9_raw_direction_and_unit_robust",
  "local9_raw_confirmation_tier"
)
require_columns(candidates, candidate_columns, "Candidate table")
if (
  nrow(candidates) != expected_candidates ||
    anyNA(candidates$gene_symbol) ||
    any(!nzchar(candidates$gene_symbol)) ||
    anyDuplicated(candidates$gene_symbol)
) {
  stop(
    "Candidate table must contain ", expected_candidates,
    " unique, nonempty gene symbols"
  )
}

processed_robust <- flag_or_false(
  candidates$all_estimable_biological_unit_checks_robust
)
raw_robust <- flag_or_false(
  candidates$local9_raw_direction_and_unit_robust
)
raw_direction <- flag_or_false(
  candidates$local9_raw_full_directional_pattern
)
if (any(raw_robust & !raw_direction)) {
  stop("Raw unit-robust candidates must also satisfy the raw direction pattern")
}

tier <- rep("tier_5_other_cross_model", nrow(candidates))
tier[raw_direction & !raw_robust & !processed_robust] <-
  "tier_4_raw_direction_only"
tier[processed_robust & !raw_robust] <-
  "tier_3_processed_unit_robust"
tier[raw_robust & !processed_robust] <-
  "tier_2_raw_unit_robust"
tier[raw_robust & processed_robust] <-
  "tier_1_raw_and_processed_unit_robust"
tier_rank <- match(
  tier,
  c(
    "tier_1_raw_and_processed_unit_robust",
    "tier_2_raw_unit_robust",
    "tier_3_processed_unit_robust",
    "tier_4_raw_direction_only",
    "tier_5_other_cross_model"
  )
)

rank_value <- suppressWarnings(as.numeric(
  candidates$cross_model_resilience_rank
))
rank_value[!is.finite(rank_value)] <- Inf
publication_order <- order(
  tier_rank,
  rank_value,
  candidates$gene_symbol
)
candidates$publication_evidence_tier <- tier
candidates$publication_evidence_tier_rank <- tier_rank
candidates$publication_priority_order <- match(
  seq_len(nrow(candidates)),
  publication_order
)
candidates$frozen_discovery_rank_unchanged <-
  candidates$cross_model_resilience_rank
candidates$raw_confirmation_is_independent_validation <- FALSE

publication_columns <- c(
  "publication_priority_order",
  "publication_evidence_tier",
  "publication_evidence_tier_rank",
  "gene_symbol",
  "frozen_discovery_rank_unchanged",
  "exploratory_cross_model_rank",
  "cross_model_resilience_score",
  "omn_log2_effect_omn_vs_sc",
  "omn_p_value",
  "omn_q_value",
  "disease_sma_vs_control_log2_effect",
  "treatment_r6_vs_scramble_log2_effect",
  "gse108094_sma_vs_control_log2_effect",
  "all_estimable_biological_unit_checks_robust",
  "local9_disease_raw_log2_effect",
  "local9_disease_raw_p_value",
  "local9_disease_raw_q_value",
  "local9_disease_raw_lolo_folds_tested",
  "local9_disease_raw_lolo_all_same_direction",
  "local9_treatment_raw_log2_effect",
  "local9_treatment_raw_p_value",
  "local9_treatment_raw_q_value",
  "local9_treatment_raw_treatment_lines",
  "local9_treatment_raw_treatment_lines_same_direction",
  "local9_raw_disease_opposes_natural",
  "local9_raw_r6_aligns_natural",
  "local9_raw_disease_lolo_robust",
  "local9_raw_treatment_line_robust",
  "local9_raw_full_directional_pattern",
  "local9_raw_direction_and_unit_robust",
  "local9_raw_confirmation_tier",
  "gse69175_disease_opposes_natural",
  "gse290980_mn_opposes_natural_resistance",
  "gse243076_adult_mn_detected",
  "gse243076_adult_mn_localized",
  "strict_splicing",
  "raw_confirmation_is_independent_validation"
)
publication_candidates <- candidates[
  publication_order,
  publication_columns,
  drop = FALSE
]

genes_in <- function(mask) {
  paste(sort(candidates$gene_symbol[mask]), collapse = ";")
}
summary_frame <- data.frame(
  metric = c(
    "analysis_status",
    "frozen_candidate_count",
    "raw_libraries",
    "raw_integrity_checks_passed",
    "raw_quantification_method",
    "whole_genome_alignment",
    "random_split_used",
    "raw_confirmation_independent_validation",
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
    "tier_5_candidates",
    "tier_1_genes",
    "tier_2_genes",
    "tier_3_genes",
    "raw_unit_robust_genes",
    "processed_unit_robust_genes"
  ),
  value = c(
    "COMPLETE",
    nrow(candidates),
    completion_value[["libraries"]],
    sum(validation_flags),
    completion_value[["quantification_method"]],
    completion_value[["whole_genome_alignment"]],
    completion_value[["random_split_used"]],
    "FALSE",
    sum(
      !is.na(candidates$local9_disease_raw_gene_id) &
        nzchar(candidates$local9_disease_raw_gene_id)
    ),
    sum(raw_direction),
    sum(processed_robust),
    sum(raw_robust),
    sum(raw_robust & processed_robust),
    sum(raw_robust | processed_robust),
    sum(tier_rank == 1L),
    sum(tier_rank == 2L),
    sum(tier_rank == 3L),
    sum(tier_rank == 4L),
    sum(tier_rank == 5L),
    genes_in(tier_rank == 1L),
    genes_in(tier_rank == 2L),
    genes_in(tier_rank == 3L),
    genes_in(raw_robust),
    genes_in(processed_robust)
  ),
  stringsAsFactors = FALSE
)

concordance_publication <- concordance[
  match(expected_contrasts, concordance$contrast),
  ,
  drop = FALSE
]

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_root, recursive = TRUE, showWarnings = FALSE)
figure_root <- file.path(manuscript_root, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

write_output(
  publication_candidates,
  file.path(output_root, "human_raw_confirmed_publication_candidates.tsv")
)
write_output(
  summary_frame,
  file.path(output_root, "publication_analysis_summary.tsv")
)
write_output(
  concordance_publication,
  file.path(output_root, "raw_vs_processed_concordance.tsv")
)
write_output(
  publication_candidates,
  file.path(manuscript_root, "supplementary_table_S1_candidates.tsv")
)
write_output(
  summary_frame,
  file.path(manuscript_root, "supplementary_table_S2_summary.tsv")
)
write_output(
  concordance_publication,
  file.path(manuscript_root, "supplementary_table_S3_concordance.tsv")
)

require_columns(
  disease,
  c("raw_log2_effect", "processed_log2_effect"),
  "Disease concordance table"
)
require_columns(
  treatment,
  c("raw_log2_effect", "processed_log2_effect"),
  "Treatment concordance table"
)

draw_concordance_panel <- function(frame, summary_row, title, point_color) {
  complete <- is.finite(frame$raw_log2_effect) &
    is.finite(frame$processed_log2_effect)
  x <- frame$raw_log2_effect[complete]
  y <- frame$processed_log2_effect[complete]
  bounds <- stats::quantile(
    c(x, y),
    probs = c(0.005, 0.995),
    na.rm = TRUE,
    names = FALSE
  )
  if (!all(is.finite(bounds)) || diff(bounds) <= 0) {
    bounds <- range(c(x, y), finite = TRUE)
  }
  graphics::plot(
    x,
    y,
    pch = 16,
    cex = 0.35,
    col = point_color,
    xlim = bounds,
    ylim = bounds,
    xlab = "Raw FASTQ Salmon effect (log2)",
    ylab = "Deposited processed-matrix effect (log2)",
    main = title,
    bty = "l"
  )
  graphics::abline(a = 0, b = 1, col = "#555555", lty = 2, lwd = 1.2)
  fit <- stats::lm(y ~ x)
  graphics::abline(fit, col = "#111111", lwd = 1.4)
  label <- sprintf(
    "n = %s\nSpearman rho = %.3f\nDirection agreement = %.3f",
    format(summary_row$matched_genes, big.mark = ","),
    summary_row$spearman_effect_correlation,
    summary_row$direction_agreement_fraction
  )
  graphics::legend(
    "topleft",
    legend = label,
    bty = "n",
    text.col = "#222222",
    cex = 0.85
  )
}

draw_concordance_figure <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 2),
    mar = c(4.5, 4.8, 3.5, 1.2),
    oma = c(0, 0, 0, 0)
  )
  draw_concordance_panel(
    disease,
    concordance_publication[1L, , drop = FALSE],
    "SMA versus control",
    grDevices::adjustcolor("#177E89", alpha.f = 0.28)
  )
  draw_concordance_panel(
    treatment,
    concordance_publication[2L, , drop = FALSE],
    "R6-MO versus scramble",
    grDevices::adjustcolor("#C45A2D", alpha.f = 0.28)
  )
}

concordance_png <- file.path(
  figure_root, "figure_1_raw_vs_processed_concordance.png"
)
grDevices::png(
  concordance_png,
  width = 6545,
  height = 3409,
  res = 600
)
draw_concordance_figure()
grDevices::dev.off()
concordance_pdf <- file.path(
  figure_root, "figure_1_raw_vs_processed_concordance.pdf"
)
grDevices::pdf(
  concordance_pdf,
  width = 10.8,
  height = 5.7,
  useDingbats = FALSE
)
draw_concordance_figure()
grDevices::dev.off()

evidence_columns <- c(
  "Processed unit\nrobust",
  "Raw disease\nopposes",
  "Raw disease\nLOLO",
  "Raw R6-MO\naligns",
  "Raw treatment\nboth lines",
  "GSE69175\ndirection",
  "GSE290980\nMN direction",
  "Adult MN\ndetected",
  "Strict splice\nrestoration"
)
evidence <- cbind(
  processed_robust,
  flag_or_false(candidates$local9_raw_disease_opposes_natural),
  flag_or_false(candidates$local9_raw_disease_lolo_robust),
  flag_or_false(candidates$local9_raw_r6_aligns_natural),
  flag_or_false(candidates$local9_raw_treatment_line_robust),
  flag_or_false(candidates$gse69175_disease_opposes_natural),
  flag_or_false(candidates$gse290980_mn_opposes_natural_resistance),
  flag_or_false(candidates$gse243076_adult_mn_detected),
  flag_or_false(candidates$strict_splicing)
)
colnames(evidence) <- evidence_columns
evidence <- evidence[publication_order, , drop = FALSE]
evidence_genes <- candidates$gene_symbol[publication_order]
evidence_tiers <- tier_rank[publication_order]

draw_evidence_figure <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(9.5, 5.8, 2.8, 1.2), xpd = FALSE)
  nr <- nrow(evidence)
  nc <- ncol(evidence)
  graphics::plot(
    NA,
    xlim = c(0.5, nc + 6.2),
    ylim = c(0.5, nr + 0.5),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "",
    bty = "n"
  )
  support_colors <- c(
    "#177E89", "#2A6FBB", "#2A6FBB", "#3A8D5D", "#3A8D5D",
    "#D08A25", "#B44C5E", "#7A5AA6", "#222222"
  )
  y_positions <- nr:1
  for (row_index in seq_len(nr)) {
    for (column_index in seq_len(nc)) {
      fill <- if (isTRUE(evidence[row_index, column_index])) {
        support_colors[column_index]
      } else {
        "#E7E9ED"
      }
      graphics::rect(
        column_index - 0.43,
        y_positions[row_index] - 0.43,
        column_index + 0.43,
        y_positions[row_index] + 0.43,
        col = fill,
        border = "white",
        lwd = 0.6
      )
    }
  }
  tier_colors <- c("#0C5962", "#2878A8", "#8E6C1F", "#8A536A", "#8B9099")
  for (row_index in seq_len(nr)) {
    graphics::rect(
      nc + 0.53,
      y_positions[row_index] - 0.43,
      nc + 0.72,
      y_positions[row_index] + 0.43,
      col = tier_colors[evidence_tiers[row_index]],
      border = NA
    )
  }
  graphics::axis(
    1,
    at = seq_len(nc),
    labels = colnames(evidence),
    las = 2,
    tick = FALSE,
    cex.axis = 0.68,
    line = -0.2
  )
  graphics::axis(
    2,
    at = y_positions,
    labels = evidence_genes,
    las = 1,
    tick = FALSE,
    cex.axis = 0.72
  )
  graphics::text(
    nc + 0.62,
    nr + 1.1,
    "Tier",
    srt = 90,
    cex = 0.7,
    font = 2
  )
  graphics::legend(
    x = nc + 1.05,
    y = nr + 0.35,
    legend = c(
      "Tier 1: raw + processed robust",
      "Tier 2: raw robust",
      "Tier 3: processed robust",
      "Tier 4: raw direction only",
      "Tier 5: other cross-model"
    ),
    fill = tier_colors,
    border = NA,
    bty = "n",
    cex = 0.68,
    xjust = 0,
    yjust = 1
  )
}

evidence_png <- file.path(
  figure_root, "figure_2_candidate_evidence_matrix.png"
)
grDevices::png(
  evidence_png,
  width = 6500,
  height = 6875,
  res = 600
)
draw_evidence_figure()
grDevices::dev.off()
evidence_pdf <- file.path(
  figure_root, "figure_2_candidate_evidence_matrix.pdf"
)
grDevices::pdf(
  evidence_pdf,
  width = 11.2,
  height = 11.5,
  useDingbats = FALSE
)
draw_evidence_figure()
grDevices::dev.off()

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "R_sessionInfo.txt")
)

cat(
  "Publication integration complete:",
  nrow(publication_candidates), "frozen candidates;",
  sum(tier_rank == 1L), "tier 1;",
  sum(tier_rank == 2L), "tier 2;",
  sum(tier_rank == 3L), "tier 3\n"
)
