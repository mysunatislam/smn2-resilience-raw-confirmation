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

arguments <- commandArgs(trailingOnly = TRUE)
argument_value <- function(prefix, default = NULL) {
  matches <- arguments[startsWith(arguments, prefix)]
  if (length(matches) == 0L) {
    return(default)
  }
  if (length(matches) != 1L) {
    stop("Expected one argument beginning with ", prefix)
  }
  sub(prefix, "", matches, fixed = TRUE)
}

rmats_root <- argument_value("--rmats-root=")
output_root <- argument_value("--output-root=")
if (is.null(rmats_root) || is.null(output_root)) {
  stop("Required arguments: --rmats-root=PATH --output-root=PATH")
}
rmats_root <- normalizePath(rmats_root, winslash = "/", mustWork = TRUE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)

frozen_path <- file.path(
  ROOT,
  "config",
  "frozen_83_splice_events_2026-08-01.tsv"
)
primary_path <- file.path(
  ROOT,
  "config",
  "frozen_12_primary_splice_events_2026-08-01.tsv"
)
frozen <- read.delim(
  frozen_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
primary <- read.delim(
  primary_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(frozen) == 83L,
  nrow(primary) == 12L,
  !anyDuplicated(frozen$event_key),
  !anyDuplicated(primary$event_key),
  all(primary$event_key %in% frozen$event_key)
)

event_types <- c("ES", "A5SS", "A3SS", "MXE", "RI")
rmats_types <- c(ES = "SE", A5SS = "A5SS", A3SS = "A3SS", MXE = "MXE", RI = "RI")
event_specs <- list(
  ES = c(
    "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  ),
  A5SS = c(
    "geneSymbol", "chr", "strand", "longExonStart_0base",
    "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"
  ),
  A3SS = c(
    "geneSymbol", "chr", "strand", "longExonStart_0base",
    "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"
  ),
  MXE = c(
    "geneSymbol", "chr", "strand", "1stExonStart_0base",
    "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  ),
  RI = c(
    "geneSymbol", "chr", "strand", "riExonStart_0base", "riExonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  )
)

format_key_value <- function(values) {
  if (is.numeric(values) || is.integer(values)) {
    return(ifelse(
      is.na(values),
      "",
      format(values, scientific = FALSE, trim = TRUE, digits = 15)
    ))
  }
  ifelse(is.na(values), "", trimws(as.character(values)))
}

make_event_key <- function(frame, event_type) {
  columns <- lapply(
    frame[, event_specs[[event_type]], drop = FALSE],
    format_key_value
  )
  do.call(paste, c(list(event_type), columns, sep = "|"))
}

parse_number_list <- function(value, expected_length) {
  pieces <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
  pieces[pieces %in% c("", "NA", "NaN", "nan", "None")] <- NA_character_
  if (length(pieces) != expected_length) {
    stop(
      "Expected ", expected_length, " comma-separated values, found ",
      length(pieces)
    )
  }
  suppressWarnings(as.numeric(pieces))
}

number_matrix <- function(values, expected_length) {
  result <- t(vapply(
    values,
    parse_number_list,
    numeric(expected_length),
    expected_length = expected_length
  ))
  storage.mode(result) <- "numeric"
  result
}

finite_row_mean <- function(matrix) {
  result <- rowMeans(matrix, na.rm = TRUE)
  result[!is.finite(result)] <- NA_real_
  result
}

finite_row_median <- function(matrix) {
  result <- apply(matrix, 1L, median, na.rm = TRUE)
  result[!is.finite(result)] <- NA_real_
  result
}

read_contrast <- function(contrast, group_lengths) {
  frames <- vector("list", length(event_types))
  names(frames) <- event_types
  for (event_type in event_types) {
    path <- file.path(
      rmats_root,
      contrast,
      paste0(rmats_types[[event_type]], ".MATS.JC.txt")
    )
    if (!file.exists(path)) {
      stop("Missing rMATS output: ", path)
    }
    frame <- read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = ""
    )
    required <- c(
      event_specs[[event_type]], "PValue", "FDR", "IncLevelDifference",
      "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "SJC_SAMPLE_1",
      "IJC_SAMPLE_2", "SJC_SAMPLE_2"
    )
    missing <- setdiff(required, names(frame))
    if (length(missing) > 0L) {
      stop(
        "Missing columns in ", basename(path), ": ",
        paste(missing, collapse = ", ")
      )
    }
    frame$event_type <- event_type
    frame$event_key <- make_event_key(frame, event_type)
    frame$exact_match_multiplicity <- ave(
      seq_len(nrow(frame)),
      frame$event_key,
      FUN = length
    )
    frames[[event_type]] <- frame
  }
  frames
}

disease_frames <- read_contrast("disease", c(2L, 3L))
treatment_frames <- read_contrast("treatment", c(2L, 2L))

empty_record <- function(event) {
  data.frame(
    event_type = event$event_type,
    event_key = event$event_key,
    gene_symbol = event$gene_symbol,
    chromosome = event$chromosome,
    strand = event$strand,
    disease_exact_matches = 0L,
    treatment_exact_matches = 0L,
    disease_raw_delta = NA_real_,
    treatment_raw_delta = NA_real_,
    disease_p_value = NA_real_,
    disease_rmats_fdr_fixed_set = NA_real_,
    treatment_p_value = NA_real_,
    treatment_rmats_fdr_fixed_set = NA_real_,
    raw_control_mean_psi = NA_real_,
    raw_sma_mean_psi = NA_real_,
    raw_scramble_mean_psi = NA_real_,
    raw_r6_mean_psi = NA_real_,
    S2_raw_disease_psi = NA_real_,
    S3_raw_disease_psi = NA_real_,
    S2_raw_scramble_psi = NA_real_,
    S3_raw_scramble_psi = NA_real_,
    S2_raw_r6_psi = NA_real_,
    S3_raw_r6_psi = NA_real_,
    S2_raw_treatment_delta = NA_real_,
    S3_raw_treatment_delta = NA_real_,
    S2_raw_distance_improvement = NA_real_,
    S3_raw_distance_improvement = NA_real_,
    sma_median_informative_jc = NA_real_,
    control_median_informative_jc = NA_real_,
    r6_median_informative_jc = NA_real_,
    scramble_median_informative_jc = NA_real_,
    minimum_group_median_informative_jc = NA_real_,
    structurally_recovered = FALSE,
    adequate_junction_support = FALSE,
    disease_direction_reproduced = FALSE,
    treatment_reversal_reproduced = FALSE,
    S2_corrected = FALSE,
    S3_corrected = FALSE,
    both_lines_corrected = FALSE,
    strong_raw_confirmation = FALSE,
    structural_recovery_reason = "",
    raw_confirmation_limiting_reason = "",
    stringsAsFactors = FALSE
  )
}

extract_unique <- function(frames, event) {
  frame <- frames[[event$event_type]]
  matches <- which(frame$event_key == event$event_key)
  list(
    frame = frame,
    matches = matches,
    row = if (length(matches) == 1L) frame[matches, , drop = FALSE] else NULL
  )
}

records <- vector("list", nrow(frozen))
support_records <- vector("list", nrow(frozen))
disease_sample_names <- c(
  "BULK-SAM-160", "BULK-SAM-148",
  "BULK-SAM-109", "BULK-SAM-104", "BULK-SAM-107"
)
treatment_sample_names <- c(
  "BULK-SAM-166", "BULK-SAM-157",
  "BULK-SAM-163", "BULK-SAM-151"
)

for (index in seq_len(nrow(frozen))) {
  event <- frozen[index, , drop = FALSE]
  disease <- extract_unique(disease_frames, event)
  treatment <- extract_unique(treatment_frames, event)
  record <- empty_record(event)
  record$disease_exact_matches <- length(disease$matches)
  record$treatment_exact_matches <- length(treatment$matches)

  if (length(disease$matches) == 1L && length(treatment$matches) == 1L) {
    disease_row <- disease$row
    treatment_row <- treatment$row
    disease_psi1 <- number_matrix(disease_row$IncLevel1, 2L)
    disease_psi2 <- number_matrix(disease_row$IncLevel2, 3L)
    treatment_psi1 <- number_matrix(treatment_row$IncLevel1, 2L)
    treatment_psi2 <- number_matrix(treatment_row$IncLevel2, 2L)
    disease_counts1 <-
      number_matrix(disease_row$IJC_SAMPLE_1, 2L) +
      number_matrix(disease_row$SJC_SAMPLE_1, 2L)
    disease_counts2 <-
      number_matrix(disease_row$IJC_SAMPLE_2, 3L) +
      number_matrix(disease_row$SJC_SAMPLE_2, 3L)
    treatment_counts1 <-
      number_matrix(treatment_row$IJC_SAMPLE_1, 2L) +
      number_matrix(treatment_row$SJC_SAMPLE_1, 2L)
    treatment_counts2 <-
      number_matrix(treatment_row$IJC_SAMPLE_2, 2L) +
      number_matrix(treatment_row$SJC_SAMPLE_2, 2L)

    record$disease_raw_delta <- as.numeric(disease_row$IncLevelDifference)
    record$treatment_raw_delta <- as.numeric(treatment_row$IncLevelDifference)
    record$disease_p_value <- as.numeric(disease_row$PValue)
    record$disease_rmats_fdr_fixed_set <- as.numeric(disease_row$FDR)
    record$treatment_p_value <- as.numeric(treatment_row$PValue)
    record$treatment_rmats_fdr_fixed_set <- as.numeric(treatment_row$FDR)
    record$raw_sma_mean_psi <- finite_row_mean(disease_psi1)
    record$raw_control_mean_psi <- finite_row_mean(disease_psi2)
    record$raw_r6_mean_psi <- finite_row_mean(treatment_psi1)
    record$raw_scramble_mean_psi <- finite_row_mean(treatment_psi2)
    record$S2_raw_disease_psi <- disease_psi1[1, 1]
    record$S3_raw_disease_psi <- disease_psi1[1, 2]
    record$S2_raw_r6_psi <- treatment_psi1[1, 1]
    record$S3_raw_r6_psi <- treatment_psi1[1, 2]
    record$S2_raw_scramble_psi <- treatment_psi2[1, 1]
    record$S3_raw_scramble_psi <- treatment_psi2[1, 2]
    record$S2_raw_treatment_delta <-
      record$S2_raw_r6_psi - record$S2_raw_scramble_psi
    record$S3_raw_treatment_delta <-
      record$S3_raw_r6_psi - record$S3_raw_scramble_psi
    record$S2_raw_distance_improvement <-
      abs(record$S2_raw_scramble_psi - record$raw_control_mean_psi) -
      abs(record$S2_raw_r6_psi - record$raw_control_mean_psi)
    record$S3_raw_distance_improvement <-
      abs(record$S3_raw_scramble_psi - record$raw_control_mean_psi) -
      abs(record$S3_raw_r6_psi - record$raw_control_mean_psi)
    record$sma_median_informative_jc <- finite_row_median(disease_counts1)
    record$control_median_informative_jc <- finite_row_median(disease_counts2)
    record$r6_median_informative_jc <- finite_row_median(treatment_counts1)
    record$scramble_median_informative_jc <- finite_row_median(treatment_counts2)
    record$minimum_group_median_informative_jc <- min(
      record$sma_median_informative_jc,
      record$control_median_informative_jc,
      record$r6_median_informative_jc,
      record$scramble_median_informative_jc,
      na.rm = TRUE
    )
    finite_required <- c(
      record$disease_raw_delta,
      record$treatment_raw_delta,
      record$raw_control_mean_psi,
      record$raw_sma_mean_psi,
      record$raw_scramble_mean_psi,
      record$raw_r6_mean_psi,
      record$S2_raw_scramble_psi,
      record$S3_raw_scramble_psi,
      record$S2_raw_r6_psi,
      record$S3_raw_r6_psi
    )
    record$structurally_recovered <- all(is.finite(finite_required))
    record$adequate_junction_support <-
      record$structurally_recovered &&
      is.finite(record$minimum_group_median_informative_jc) &&
      record$minimum_group_median_informative_jc >= 10
    record$disease_direction_reproduced <-
      record$structurally_recovered &&
      sign(record$disease_raw_delta) ==
        sign(event$disease_delta_sma_minus_control) &&
      abs(record$disease_raw_delta) >= 0.05
    record$treatment_reversal_reproduced <-
      record$structurally_recovered &&
      sign(record$treatment_raw_delta) ==
        sign(event$treatment_delta_r6_minus_scramble) &&
      sign(record$treatment_raw_delta) == -sign(record$disease_raw_delta) &&
      abs(record$treatment_raw_delta) >= 0.05
    record$S2_corrected <-
      record$structurally_recovered &&
      sign(record$S2_raw_treatment_delta) == -sign(record$disease_raw_delta) &&
      record$S2_raw_distance_improvement >= 0.05
    record$S3_corrected <-
      record$structurally_recovered &&
      sign(record$S3_raw_treatment_delta) == -sign(record$disease_raw_delta) &&
      record$S3_raw_distance_improvement >= 0.05
    record$both_lines_corrected <- record$S2_corrected && record$S3_corrected

    support_values <- c(
      disease_counts1[1, ],
      disease_counts2[1, ],
      treatment_counts1[1, ],
      treatment_counts2[1, ]
    )
    support_record <- data.frame(
      event_key = rep(event$event_key, 9L),
      event_type = rep(event$event_type, 9L),
      gene_symbol = rep(event$gene_symbol, 9L),
      sample_id = c(disease_sample_names, treatment_sample_names),
      contrast = c(rep("disease", 5L), rep("treatment", 4L)),
      group = c(
        rep("SMA", 2L), rep("control", 3L),
        rep("R6-MO", 2L), rep("scramble", 2L)
      ),
      informative_junction_count = support_values,
      stringsAsFactors = FALSE
    )
  } else {
    support_record <- data.frame(
      event_key = character(),
      event_type = character(),
      gene_symbol = character(),
      sample_id = character(),
      contrast = character(),
      group = character(),
      informative_junction_count = numeric(),
      stringsAsFactors = FALSE
    )
  }

  if (
    length(disease$matches) > 1L ||
      length(treatment$matches) > 1L
  ) {
    record$structural_recovery_reason <- "ambiguous_exact_match"
  } else if (
    length(disease$matches) == 0L ||
      length(treatment$matches) == 0L
  ) {
    same_gene_disease <- any(
      disease$frame$geneSymbol == event$gene_symbol
    )
    same_gene_treatment <- any(
      treatment$frame$geneSymbol == event$gene_symbol
    )
    record$structural_recovery_reason <- if (
      same_gene_disease || same_gene_treatment
    ) {
      "annotation_or_chromosome_mismatch"
    } else {
      "no_exact_structural_match"
    }
  } else if (!record$structurally_recovered) {
    record$structural_recovery_reason <- "exact_structure_without_finite_psi"
  } else {
    record$structural_recovery_reason <- "recovered"
  }

  records[[index]] <- record
  support_records[[index]] <- support_record
}

results <- do.call(rbind, records)
rownames(results) <- NULL
support <- do.call(rbind, support_records)
rownames(support) <- NULL

results$disease_within83_q <- p.adjust(results$disease_p_value, method = "BH")
results$treatment_within83_q <- p.adjust(
  results$treatment_p_value,
  method = "BH"
)
results$strong_raw_confirmation <-
  results$structurally_recovered &
  results$adequate_junction_support &
  results$disease_direction_reproduced &
  results$treatment_reversal_reproduced &
  results$both_lines_corrected &
  !is.na(results$disease_within83_q) &
  results$disease_within83_q < 0.05 &
  !is.na(results$treatment_within83_q) &
  results$treatment_within83_q < 0.05

limiting_reason <- function(index) {
  row <- results[index, , drop = FALSE]
  if (row$strong_raw_confirmation) {
    return("strong_raw_confirmation")
  }
  if (!row$structurally_recovered) {
    return(row$structural_recovery_reason)
  }
  if (!row$adequate_junction_support) {
    return("insufficient_junction_support")
  }
  if (!row$disease_direction_reproduced) {
    return("disease_direction_not_reproduced")
  }
  if (!row$treatment_reversal_reproduced) {
    return("treatment_reversal_not_reproduced")
  }
  if (!row$both_lines_corrected) {
    return("correction_not_reproduced_in_both_lines")
  }
  if (
    is.na(row$disease_within83_q) ||
      row$disease_within83_q >= 0.05
  ) {
    return("disease_within83_q_not_below_0.05")
  }
  "treatment_within83_q_not_below_0.05"
}
results$raw_confirmation_limiting_reason <- vapply(
  seq_len(nrow(results)),
  limiting_reason,
  character(1)
)

frozen_columns <- frozen[, c(
  "event_key", "disease_delta_sma_minus_control",
  "treatment_delta_r6_minus_scramble", "disease_fdr", "treatment_fdr"
)]
names(frozen_columns)[2:5] <- c(
  "processed_disease_delta",
  "processed_treatment_delta",
  "processed_disease_fdr",
  "processed_treatment_fdr"
)
results <- merge(
  results,
  frozen_columns,
  by = "event_key",
  all.x = TRUE,
  sort = FALSE
)
results <- results[match(frozen$event_key, results$event_key), , drop = FALSE]
rownames(results) <- NULL
stopifnot(identical(results$event_key, frozen$event_key))

finite_disease <- results$structurally_recovered &
  is.finite(results$disease_raw_delta) &
  is.finite(results$processed_disease_delta)
finite_treatment <- results$structurally_recovered &
  is.finite(results$treatment_raw_delta) &
  is.finite(results$processed_treatment_delta)

safe_correlation <- function(x, y, method) {
  if (length(x) < 3L) {
    return(NA_real_)
  }
  suppressWarnings(cor(x, y, method = method))
}

disease_spearman <- safe_correlation(
  results$processed_disease_delta[finite_disease],
  results$disease_raw_delta[finite_disease],
  "spearman"
)
disease_pearson <- safe_correlation(
  results$processed_disease_delta[finite_disease],
  results$disease_raw_delta[finite_disease],
  "pearson"
)
treatment_spearman <- safe_correlation(
  results$processed_treatment_delta[finite_treatment],
  results$treatment_raw_delta[finite_treatment],
  "spearman"
)
treatment_pearson <- safe_correlation(
  results$processed_treatment_delta[finite_treatment],
  results$treatment_raw_delta[finite_treatment],
  "pearson"
)

primary_results <- results[
  match(primary$event_key, results$event_key),
  ,
  drop = FALSE
]
primary_results$panel_order <- primary$panel_order
primary_results$source_gene_symbol <- primary$source_gene_symbol
primary_results <- primary_results[
  order(primary_results$panel_order),
  ,
  drop = FALSE
]

panel_recovered <- sum(primary_results$structurally_recovered)
panel_disease <- sum(primary_results$disease_direction_reproduced)
panel_full_correction <- sum(
  primary_results$treatment_reversal_reproduced &
    primary_results$both_lines_corrected &
    primary_results$adequate_junction_support
)
panel_attempted <- nrow(primary_results) == 12L
panel_correlation_eligible <- sum(finite_disease) >= 30L
panel_success <- all(
  panel_attempted,
  panel_recovered >= 8L,
  panel_disease >= 6L,
  panel_full_correction >= 4L,
  panel_correlation_eligible,
  is.finite(disease_spearman),
  disease_spearman >= 0.50
)

summary_table <- data.frame(
  metric = c(
    "frozen_events",
    "primary_events_attempted",
    "structurally_recovered_83",
    "adequate_junction_support_83",
    "disease_direction_reproduced_83",
    "treatment_reversal_reproduced_83",
    "both_lines_corrected_83",
    "strong_raw_confirmation_83",
    "primary_structurally_recovered",
    "primary_disease_direction_reproduced",
    "primary_treatment_reversal_both_lines_adequate_support",
    "disease_raw_processed_pairs",
    "disease_raw_processed_spearman",
    "disease_raw_processed_pearson",
    "treatment_raw_processed_pairs",
    "treatment_raw_processed_spearman",
    "treatment_raw_processed_pearson",
    "primary_panel_success"
  ),
  value = c(
    83L,
    panel_attempted,
    sum(results$structurally_recovered),
    sum(results$adequate_junction_support),
    sum(results$disease_direction_reproduced),
    sum(results$treatment_reversal_reproduced),
    sum(results$both_lines_corrected),
    sum(results$strong_raw_confirmation),
    panel_recovered,
    panel_disease,
    panel_full_correction,
    sum(finite_disease),
    disease_spearman,
    disease_pearson,
    sum(finite_treatment),
    treatment_spearman,
    treatment_pearson,
    panel_success
  ),
  stringsAsFactors = FALSE
)

concordance <- rbind(
  data.frame(
    contrast = "disease_SMA_minus_control",
    event_key = results$event_key[finite_disease],
    event_type = results$event_type[finite_disease],
    gene_symbol = results$gene_symbol[finite_disease],
    processed_delta_psi = results$processed_disease_delta[finite_disease],
    raw_delta_psi = results$disease_raw_delta[finite_disease],
    direction_agreement =
      sign(results$processed_disease_delta[finite_disease]) ==
        sign(results$disease_raw_delta[finite_disease]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    contrast = "treatment_R6_minus_scramble",
    event_key = results$event_key[finite_treatment],
    event_type = results$event_type[finite_treatment],
    gene_symbol = results$gene_symbol[finite_treatment],
    processed_delta_psi = results$processed_treatment_delta[finite_treatment],
    raw_delta_psi = results$treatment_raw_delta[finite_treatment],
    direction_agreement =
      sign(results$processed_treatment_delta[finite_treatment]) ==
        sign(results$treatment_raw_delta[finite_treatment]),
    stringsAsFactors = FALSE
  )
)

write_tsv <- function(frame, filename) {
  write.table(
    frame,
    file.path(output_root, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}
write_tsv(results, "frozen_83_raw_splice_confirmation.tsv")
write_tsv(primary_results, "frozen_12_primary_raw_splice_confirmation.tsv")
write_tsv(support, "frozen_83_junction_read_support.tsv")
write_tsv(
  results[!results$strong_raw_confirmation, c(
    "event_key", "event_type", "gene_symbol", "structurally_recovered",
    "adequate_junction_support", "disease_direction_reproduced",
    "treatment_reversal_reproduced", "both_lines_corrected",
    "structural_recovery_reason", "raw_confirmation_limiting_reason"
  )],
  "frozen_83_unconfirmed_reasons.tsv"
)
write_tsv(concordance, "raw_vs_processed_delta_psi.tsv")
write_tsv(summary_table, "raw_splice_confirmation_summary.tsv")

plot_path <- file.path(output_root, "raw_vs_processed_delta_psi.pdf")
grDevices::pdf(plot_path, width = 8.5, height = 4.25, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
plot(
  results$processed_disease_delta[finite_disease],
  results$disease_raw_delta[finite_disease],
  pch = 19,
  col = "#176B87",
  xlab = "Processed disease delta PSI",
  ylab = "Raw disease delta PSI",
  main = sprintf("Disease: Spearman %.2f", disease_spearman)
)
abline(h = 0, v = 0, col = "grey75")
abline(0, 1, lty = 2, col = "grey35")
plot(
  results$processed_treatment_delta[finite_treatment],
  results$treatment_raw_delta[finite_treatment],
  pch = 19,
  col = "#B04A5A",
  xlab = "Processed treatment delta PSI",
  ylab = "Raw treatment delta PSI",
  main = sprintf("Treatment: Spearman %.2f", treatment_spearman)
)
abline(h = 0, v = 0, col = "grey75")
abline(0, 1, lty = 2, col = "grey35")
grDevices::dev.off()

provenance <- data.frame(
  metric = c(
    "analysis_completed",
    "analysis_role",
    "frozen_event_manifest",
    "frozen_primary_manifest",
    "event_matching",
    "disease_group1",
    "disease_group2",
    "treatment_group1",
    "treatment_group2",
    "rmats_fdr_scope",
    "confirmatory_multiplicity",
    "primary_panel_success"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    "prospective_targeted_raw_splice_junction_confirmation",
    basename(frozen_path),
    basename(primary_path),
    "exact_event_class_chromosome_strand_and_all_coordinates",
    "SMA_untreated_S2_S3",
    "control_untreated_C1_C2_C3",
    "R6-MO_S2_S3",
    "scramble_S2_S3",
    "rMATS_fixed_83_event_set_not_genome_wide",
    "BH_within_frozen_83",
    as.character(panel_success)
  ),
  stringsAsFactors = FALSE
)
write_tsv(provenance, "raw_splice_confirmation_provenance.tsv")

mandatory <- c(
  "frozen_83_raw_splice_confirmation.tsv",
  "frozen_12_primary_raw_splice_confirmation.tsv",
  "frozen_83_junction_read_support.tsv",
  "frozen_83_unconfirmed_reasons.tsv",
  "raw_vs_processed_delta_psi.tsv",
  "raw_splice_confirmation_summary.tsv",
  "raw_vs_processed_delta_psi.pdf",
  "raw_splice_confirmation_provenance.tsv"
)
mandatory_paths <- file.path(output_root, mandatory)
stopifnot(
  all(file.exists(mandatory_paths)),
  all(file.info(mandatory_paths)$size > 0L),
  nrow(results) == 83L,
  nrow(primary_results) == 12L,
  nrow(summary_table) == 18L
)

write_tsv(
  data.frame(
    status = "COMPLETE",
    completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    frozen_events = 83L,
    primary_events = 12L,
    primary_panel_success = panel_success,
    note = paste(
      "Statistical confirmation tables complete;",
      "12 sashimi plots are a separate mandatory visualization phase."
    ),
    stringsAsFactors = FALSE
  ),
  "RAW_RMATS_TABLE_ANALYSIS_COMPLETE.tsv"
)

cat(
  "Raw rMATS table analysis complete:",
  sum(results$structurally_recovered), "of 83 structurally recovered;",
  panel_recovered, "of 12 primary recovered;",
  "panel success =", panel_success, "\n"
)
