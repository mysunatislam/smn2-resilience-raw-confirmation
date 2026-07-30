source(file.path("r", "common.R"))

suppressPackageStartupMessages({
  suppressWarnings(library(readxl))
})

event_types <- c("ES", "A5SS", "A3SS", "MXE", "RI")
supplement_directory <- file.path(
  ROOT,
  "data",
  "raw",
  "GSE290979",
  "publication_supplement"
)
disease_workbook <- file.path(
  supplement_directory,
  "41467_2025_67725_MOESM9_ESM.xlsx"
)
treatment_workbook <- file.path(
  supplement_directory,
  "41467_2025_67725_MOESM10_ESM.xlsx"
)
strict_path <- file.path(
  ROOT,
  "results",
  "r",
  "splicing",
  "GSE290979_R_strict_corrected_events.tsv"
)
annotation_path <- file.path(
  ROOT,
  "results",
  "r",
  "integration",
  "human_R_strict_splicing_candidates.tsv"
)

strict <- read.delim(
  strict_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
annotations <- read.delim(
  annotation_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(strict) == 83L,
  length(unique(strict$gene_symbol)) == 74L,
  all(strict$strict_corrected_rescue),
  all(strict$event_type %in% event_types)
)

parse_count_vector <- function(value, expected_length) {
  pieces <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
  if (length(pieces) != expected_length) {
    stop("Expected ", expected_length, " junction counts, found ", length(pieces))
  }
  counts <- suppressWarnings(as.numeric(pieces))
  if (anyNA(counts) || any(counts < 0)) {
    stop("Invalid junction-count vector")
  }
  counts
}

parse_junction_support <- function(row, group, expected_length) {
  ijc_value <- as.character(row[[paste0("IJC_SAMPLE_", group)]])
  sjc <- parse_count_vector(
    row[[paste0("SJC_SAMPLE_", group)]],
    expected_length
  )
  ijc_pieces <- strsplit(ijc_value, ",", fixed = TRUE)[[1]]
  if (length(ijc_pieces) == expected_length) {
    return(list(
      counts = parse_count_vector(ijc_value, expected_length) + sjc,
      exact = TRUE
    ))
  }
  if (
    length(ijc_pieces) == 1L &&
      grepl("^[0-9.]+[Ee][+]?[0-9]+$", ijc_value)
  ) {
    return(list(counts = sjc, exact = FALSE))
  }
  stop(
    "Unrecognized IJC vector for event ", row$ID,
    ", group ", group, ": ", ijc_value
  )
}

event_coordinates <- function(row, event_type) {
  empty <- data.frame(
    primary_feature = NA_character_,
    primary_start_0base = NA_real_,
    primary_end = NA_real_,
    secondary_feature = NA_character_,
    secondary_start_0base = NA_real_,
    secondary_end = NA_real_,
    upstream_start_0base = NA_real_,
    upstream_end = NA_real_,
    downstream_start_0base = NA_real_,
    downstream_end = NA_real_,
    stringsAsFactors = FALSE
  )
  if (event_type == "ES") {
    empty$primary_feature <- "skipped_exon"
    empty$primary_start_0base <- row$exonStart_0base
    empty$primary_end <- row$exonEnd
    empty$upstream_start_0base <- row$upstreamES
    empty$upstream_end <- row$upstreamEE
    empty$downstream_start_0base <- row$downstreamES
    empty$downstream_end <- row$downstreamEE
  } else if (event_type %in% c("A5SS", "A3SS")) {
    empty$primary_feature <- "long_exon"
    empty$primary_start_0base <- row$longExonStart_0base
    empty$primary_end <- row$longExonEnd
    empty$secondary_feature <- "short_exon"
    empty$secondary_start_0base <- row$shortES
    empty$secondary_end <- row$shortEE
    empty$downstream_start_0base <- row$flankingES
    empty$downstream_end <- row$flankingEE
  } else if (event_type == "MXE") {
    empty$primary_feature <- "first_exon"
    empty$primary_start_0base <- row$X1stExonStart_0base
    empty$primary_end <- row$X1stExonEnd
    empty$secondary_feature <- "second_exon"
    empty$secondary_start_0base <- row$X2ndExonStart_0base
    empty$secondary_end <- row$X2ndExonEnd
    empty$upstream_start_0base <- row$upstreamES
    empty$upstream_end <- row$upstreamEE
    empty$downstream_start_0base <- row$downstreamES
    empty$downstream_end <- row$downstreamEE
  } else if (event_type == "RI") {
    empty$primary_feature <- "retained_intron"
    empty$primary_start_0base <- row$riExonStart_0base
    empty$primary_end <- row$riExonEnd
    empty$upstream_start_0base <- row$upstreamES
    empty$upstream_end <- row$upstreamEE
    empty$downstream_start_0base <- row$downstreamES
    empty$downstream_end <- row$downstreamEE
  }
  empty
}

read_sheets <- function(path) {
  frames <- setNames(vector("list", length(event_types)), event_types)
  for (event_type in event_types) {
    frame <- as.data.frame(
      suppressWarnings(read_excel(path, sheet = event_type)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    frames[[event_type]] <- frame
  }
  frames
}

disease_frames <- read_sheets(disease_workbook)
treatment_frames <- read_sheets(treatment_workbook)
records <- vector("list", nrow(strict))

for (index in seq_len(nrow(strict))) {
  event <- strict[index, , drop = FALSE]
  event_type <- event$event_type
  disease <- disease_frames[[event_type]]
  treatment <- treatment_frames[[event_type]]
  disease_row <- disease[
    as.character(disease$ID) == as.character(event$disease_event_id),
    ,
    drop = FALSE
  ]
  treatment_row <- treatment[
    as.character(treatment$ID) == as.character(event$treatment_event_id),
    ,
    drop = FALSE
  ]
  if (nrow(disease_row) != 1L || nrow(treatment_row) != 1L) {
    stop("Could not uniquely recover source rows for ", event$event_key)
  }

  sma_support <- parse_junction_support(disease_row, 1L, 6L)
  control_support <- parse_junction_support(disease_row, 2L, 9L)
  scramble_support <- parse_junction_support(treatment_row, 1L, 8L)
  r6_support <- parse_junction_support(treatment_row, 2L, 8L)
  sma_counts <- sma_support$counts
  control_counts <- control_support$counts
  scramble_counts <- scramble_support$counts
  r6_counts <- r6_support$counts
  support_exact <- all(c(
    sma_support$exact, control_support$exact,
    scramble_support$exact, r6_support$exact
  ))

  group_medians <- c(
    sma = median(sma_counts),
    control = median(control_counts),
    scramble = median(scramble_counts),
    r6 = median(r6_counts)
  )
  all_counts <- c(sma_counts, control_counts, scramble_counts, r6_counts)
  coordinates <- event_coordinates(disease_row, event_type)
  records[[index]] <- cbind(
    data.frame(
      event_type = event_type,
      event_key = event$event_key,
      source_gene_symbol = event$gene_symbol,
      chromosome = event$chromosome,
      strand = event$strand,
      disease_event_id = event$disease_event_id,
      treatment_event_id = event$treatment_event_id,
      disease_fdr = event$disease_fdr,
      treatment_fdr = event$treatment_fdr,
      disease_delta_sma_minus_control =
        event$disease_delta_sma_minus_control,
      treatment_delta_r6_minus_scramble =
        event$treatment_delta_r6_minus_scramble,
      S2_distance_improvement = event$S2_distance_improvement,
      S3_distance_improvement = event$S3_distance_improvement,
      sma_median_informative_jc = group_medians[["sma"]],
      control_median_informative_jc = group_medians[["control"]],
      scramble_median_informative_jc = group_medians[["scramble"]],
      r6_median_informative_jc = group_medians[["r6"]],
      minimum_group_median_informative_jc = min(group_medians),
      minimum_sample_informative_jc = min(all_counts),
      fraction_samples_ge_10_informative_jc = mean(all_counts >= 10),
      fraction_samples_ge_20_informative_jc = mean(all_counts >= 20),
      junction_support_exact = support_exact,
      junction_support_basis = if (support_exact) {
        "exact_ijc_plus_sjc"
      } else {
        "conservative_sjc_only_for_excel_coerced_ijc"
      },
      stringsAsFactors = FALSE
    ),
    coordinates
  )
}

ranking <- do.call(rbind, records)
rownames(ranking) <- NULL

annotation <- annotations[, c(
  "gene_symbol", "splicing_source_symbol", "hgnc_name", "hgnc_locus_group",
  "hgnc_locus_type", "hgnc_location", "hgnc_status",
  "omn_log2_effect_omn_vs_sc", "omn_p_value", "omn_q_value",
  "exploratory_omn_positive_p05"
)]
names(annotation)[1:2] <- c("current_gene_symbol", "source_gene_symbol")
annotation <- annotation[!duplicated(annotation$source_gene_symbol), , drop = FALSE]
ranking <- merge(
  ranking,
  annotation,
  by = "source_gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
ranking$current_gene_symbol[is.na(ranking$current_gene_symbol)] <-
  ranking$source_gene_symbol[is.na(ranking$current_gene_symbol)]

event_counts <- table(ranking$current_gene_symbol)
ranking$strict_events_for_gene <- as.integer(
  event_counts[ranking$current_gene_symbol]
)
assay_weight <- c(ES = 1.00, A5SS = 0.90, A3SS = 0.90, MXE = 0.75, RI = 0.85)
ranking$effect_component <- pmin(
  pmin(
    abs(ranking$disease_delta_sma_minus_control),
    abs(ranking$treatment_delta_r6_minus_scramble)
  ) / 0.20,
  1
)
ranking$significance_component <- (
  pmin(-log10(pmax(ranking$disease_fdr, .Machine$double.xmin)) / 5, 1) +
    pmin(-log10(pmax(ranking$treatment_fdr, .Machine$double.xmin)) / 5, 1)
) / 2
ranking$minimum_line_distance_improvement <- pmin(
  ranking$S2_distance_improvement,
  ranking$S3_distance_improvement
)
ranking$consistency_component <- pmin(
  pmax(ranking$minimum_line_distance_improvement, 0) / 0.10,
  1
)
ranking$support_component <- 0.60 * pmin(
  log1p(ranking$minimum_group_median_informative_jc) / log1p(50),
  1
) + 0.40 * ranking$fraction_samples_ge_10_informative_jc
ranking$assay_component <- unname(assay_weight[ranking$event_type])
ranking$validation_priority_score <-
  0.30 * ranking$effect_component +
  0.25 * ranking$significance_component +
  0.25 * ranking$consistency_component +
  0.15 * ranking$support_component +
  0.05 * ranking$assay_component
ranking$support_class <- ifelse(
  ranking$minimum_group_median_informative_jc >= 50 &
    ranking$fraction_samples_ge_10_informative_jc >= 0.80,
  "HIGH",
  ifelse(
    ranking$minimum_group_median_informative_jc >= 20 &
      ranking$fraction_samples_ge_10_informative_jc >= 0.70,
    "MODERATE",
    "LOW"
  )
)

ranking <- ranking[order(
  -ranking$validation_priority_score,
  ranking$disease_fdr,
  ranking$treatment_fdr,
  ranking$event_key
), , drop = FALSE]
rownames(ranking) <- NULL
ranking$validation_rank <- seq_len(nrow(ranking))

approved <- !is.na(ranking$hgnc_status) & ranking$hgnc_status == "Approved"
selected_indices <- integer()
selection_reasons <- character()
for (event_type in event_types) {
  candidates <- which(
    ranking$event_type == event_type & approved &
      !ranking$current_gene_symbol %in% ranking$current_gene_symbol[selected_indices]
  )
  if (length(candidates) == 0L) stop("No approved candidate for ", event_type)
  selected_indices <- c(selected_indices, candidates[1])
  selection_reasons <- c(
    selection_reasons,
    paste0("highest_scoring_", event_type, "_anchor")
  )
}
for (index in which(approved)) {
  if (length(selected_indices) >= 12L) break
  if (
    index %in% selected_indices ||
    ranking$current_gene_symbol[index] %in%
      ranking$current_gene_symbol[selected_indices]
  ) {
    next
  }
  selected_indices <- c(selected_indices, index)
  selection_reasons <- c(selection_reasons, "highest_remaining_unique_gene")
}
stopifnot(length(selected_indices) == 12L)
panel <- ranking[selected_indices, , drop = FALSE]
panel$selection_reason <- selection_reasons
panel <- panel[order(-panel$validation_priority_score), , drop = FALSE]
panel$panel_order <- seq_len(nrow(panel))

stopifnot(
  nrow(ranking) == 83L,
  length(unique(ranking$current_gene_symbol)) == 74L,
  identical(
    as.integer(table(factor(ranking$event_type, levels = event_types))),
    c(35L, 6L, 5L, 22L, 15L)
  ),
  all(is.finite(ranking$validation_priority_score)),
  all(ranking$validation_priority_score >= 0 &
    ranking$validation_priority_score <= 1),
  sum(!ranking$junction_support_exact) == 1L,
  identical(
    ranking$current_gene_symbol[!ranking$junction_support_exact],
    "ECD"
  ),
  nrow(panel) == 12L,
  length(unique(panel$current_gene_symbol)) == 12L,
  all(event_types %in% panel$event_type)
)

bridge_panel <- ranking[
  !is.na(ranking$exploratory_omn_positive_p05) &
    ranking$exploratory_omn_positive_p05,
  ,
  drop = FALSE
]
stopifnot(
  nrow(bridge_panel) == 3L,
  setequal(
    bridge_panel$current_gene_symbol,
    c("KCNAB3", "LINC00665", "TSPOAP1")
  )
)

score_specification <- data.frame(
  component = c(
    "effect", "significance", "both_line_consistency",
    "junction_support", "assay_feasibility"
  ),
  weight = c(0.30, 0.25, 0.25, 0.15, 0.05),
  definition = c(
    "min(abs(disease_delta),abs(treatment_delta))/0.20, capped at 1",
    "mean of disease and treatment -log10(FDR)/5, each capped at 1",
    "min(S2_distance_improvement,S3_distance_improvement)/0.10, bounded 0 to 1",
    "0.60*log1p(minimum group-median informative junction count)/log1p(50), capped at 1, plus 0.40*fraction of samples with at least 10 counts",
    "event-class feasibility: ES=1.00, A5SS=0.90, A3SS=0.90, RI=0.85, MXE=0.75"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  ranking,
  "results/r/splicing/GSE290979_R_validation_event_ranking.tsv"
)
write_tsv(
  panel,
  "results/r/splicing/GSE290979_R_primary_validation_panel.tsv"
)
write_tsv(
  bridge_panel,
  "results/r/splicing/GSE290979_R_exploratory_OMN_bridge_panel.tsv"
)
write_tsv(
  score_specification,
  "results/r/splicing/GSE290979_R_validation_score_specification.tsv"
)

top <- head(ranking, 20L)
colors <- c(
  ES = "#3C5488",
  A5SS = "#009E73",
  A3SS = "#E69F00",
  MXE = "#C44E52",
  RI = "#7A5195"
)
labels <- paste0(top$current_gene_symbol, " (", top$event_type, ")")
png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "GSE290979_R_validation_priority_top20.png"
  ),
  width = 1800,
  height = 1500,
  res = 180
)
par(mar = c(5, 12, 4, 2))
y <- rev(seq_len(nrow(top)))
plot(
  top$validation_priority_score,
  y,
  pch = 19,
  cex = 1.4,
  col = unname(colors[top$event_type]),
  yaxt = "n",
  ylab = "",
  xlab = "Validation-priority score",
  xlim = c(0, 1),
  main = "GSE290979 strict splicing events: top 20 for validation"
)
segments(0, y, top$validation_priority_score, y, col = "grey80", lwd = 2)
points(
  top$validation_priority_score,
  y,
  pch = 19,
  cex = 1.4,
  col = unname(colors[top$event_type])
)
axis(2, at = y, labels = labels, las = 1)
abline(v = seq(0, 1, by = 0.2), col = "grey92", lty = 3)
legend(
  "bottomright",
  legend = names(colors),
  col = unname(colors),
  pch = 19,
  bty = "n",
  ncol = 2
)
dev.off()

component_matrix <- as.matrix(panel[, c(
  "effect_component", "significance_component", "consistency_component",
  "support_component", "assay_component"
)])
rownames(component_matrix) <- paste0(
  panel$current_gene_symbol,
  " (",
  panel$event_type,
  ")"
)
png(
  file.path(
    ROOT,
    "results",
    "r",
    "figures",
    "GSE290979_R_validation_panel_components.png"
  ),
  width = 1800,
  height = 1300,
  res = 180
)
par(mar = c(10, 12, 4, 2))
image(
  x = seq_len(ncol(component_matrix)),
  y = seq_len(nrow(component_matrix)),
  z = t(component_matrix[nrow(component_matrix):1, , drop = FALSE]),
  col = hcl.colors(50, "YlGnBu", rev = FALSE),
  zlim = c(0, 1),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = "Transparent score components for the 12-event panel"
)
axis(
  1,
  at = seq_len(ncol(component_matrix)),
  labels = c("Effect", "FDR", "Both-line", "Junction support", "Assay"),
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
    text(column, row, sprintf("%.2f", value), cex = 0.8)
  }
}
box()
dev.off()

summary <- data.frame(
  metric = c(
    "strict_events_ranked", "strict_genes", "primary_panel_events",
    "primary_panel_unique_genes", "event_classes_represented",
    "exploratory_OMN_bridge_events", "high_support_events",
    "moderate_support_events", "low_support_events",
    "source_excel_coerced_support_events"
  ),
  value = c(
    nrow(ranking), length(unique(ranking$current_gene_symbol)), nrow(panel),
    length(unique(panel$current_gene_symbol)),
    length(unique(panel$event_type)),
    nrow(bridge_panel),
    sum(ranking$support_class == "HIGH"),
    sum(ranking$support_class == "MODERATE"),
    sum(ranking$support_class == "LOW"),
    sum(!ranking$junction_support_exact)
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  summary,
  "results/r/splicing/GSE290979_R_validation_ranking_summary.tsv"
)

cat(
  "Validation ranking complete:",
  nrow(ranking),
  "events ranked; 12-event panel =",
  paste(panel$current_gene_symbol, collapse = ", "),
  "\n"
)
