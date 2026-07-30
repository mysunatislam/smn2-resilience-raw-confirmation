source(file.path("r", "common.R"))

suppressPackageStartupMessages({
  suppressWarnings(library(readxl))
})

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
published_workbook <- file.path(
  supplement_directory,
  "41467_2025_67725_MOESM11_ESM.xlsx"
)
source_workbook <- file.path(
  supplement_directory,
  "41467_2025_67725_MOESM15_ESM.xlsx"
)

event_types <- c("ES", "A5SS", "A3SS", "MXE", "RI")
event_specs <- list(
  ES = c(
    "gene_symbol", "chr", "strand", "exonStart_0base", "exonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  ),
  A5SS = c(
    "gene_symbol", "chr", "strand", "longExonStart_0base",
    "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"
  ),
  A3SS = c(
    "gene_symbol", "chr", "strand", "longExonStart_0base",
    "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"
  ),
  MXE = c(
    "gene_symbol", "chr", "strand", "X1stExonStart_0base",
    "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  ),
  RI = c(
    "gene_symbol", "chr", "strand", "riExonStart_0base", "riExonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  )
)
published_specs <- list(
  ES = c("gene_symbol", "chr", "strand", "exonStart_0base", "exonEnd"),
  A5SS = c(
    "gene_symbol", "chr", "strand", "longExonStart_0base", "longExonEnd"
  ),
  A3SS = c(
    "gene_symbol", "chr", "strand", "longExonStart_0base", "longExonEnd"
  ),
  MXE = c(
    "gene_symbol", "chr", "strand", "X1stExonStart_0base", "X1stExonEnd"
  ),
  RI = c(
    "gene_symbol", "chr", "strand", "riExonStart_0base", "riExonEnd"
  )
)

scramble_samples <- list(
  S2 = paste0("BULK-SAM-", 161:164),
  S3 = paste0("BULK-SAM-", 150:153)
)
r6_samples <- list(
  S2 = paste0("BULK-SAM-", 165:168),
  S3 = paste0("BULK-SAM-", 154:157)
)

normalize_gene_symbol <- normalize_excel_gene_symbol

format_key_values <- function(values) {
  if (is.numeric(values)) {
    return(ifelse(
      is.na(values),
      "",
      format(values, scientific = FALSE, trim = TRUE, digits = 15)
    ))
  }
  ifelse(is.na(values), "", trimws(as.character(values)))
}

make_event_key <- function(frame, event_type, specs = event_specs) {
  columns <- lapply(
    frame[, specs[[event_type]], drop = FALSE],
    format_key_values
  )
  do.call(
    paste,
    c(list(event_type), columns, sep = "|")
  )
}

parse_psi <- function(value, expected_length) {
  pieces <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1]])
  if (length(pieces) != expected_length) {
    stop(
      "Expected ",
      expected_length,
      " PSI values, found ",
      length(pieces)
    )
  }
  pieces[pieces %in% c("NA", "nan", "", "None")] <- NA_character_
  as.numeric(pieces)
}

psi_matrix <- function(values, expected_length) {
  t(vapply(
    values,
    parse_psi,
    numeric(expected_length),
    expected_length = expected_length
  ))
}

read_event_sheets <- function(path) {
  frames <- setNames(vector("list", length(event_types)), event_types)
  for (event_type in event_types) {
    frame <- as.data.frame(
      suppressWarnings(read_excel(path, sheet = event_type)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(frame)[names(frame) == "ID"] <- "source_event_id"
    names(frame)[names(frame) == "geneSymbol"] <- "gene_symbol"
    frame$gene_symbol <- normalize_gene_symbol(frame$gene_symbol)
    frame$event_type <- event_type
    frame$event_key <- make_event_key(frame, event_type)
    if (anyDuplicated(frame$event_key)) {
      stop("Duplicated full event key in ", basename(path), " ", event_type)
    }
    frames[[event_type]] <- frame
  }
  frames
}

add_disease_metrics <- function(frame) {
  sma_psi <- psi_matrix(frame$IncLevel1, 6L)
  control_psi <- psi_matrix(frame$IncLevel2, 9L)
  frame$sma_mean_inclusion <- rowMeans(sma_psi, na.rm = TRUE)
  frame$control_mean_inclusion <- rowMeans(control_psi, na.rm = TRUE)
  frame$disease_delta_sma_minus_control <- as.numeric(
    frame$IncLevelDifference
  )
  calculated <- frame$sma_mean_inclusion - frame$control_mean_inclusion
  stopifnot(max(
    abs(calculated - frame$disease_delta_sma_minus_control),
    na.rm = TRUE
  ) < 0.0021)
  frame
}

add_treatment_metrics <- function(frame) {
  scramble_psi <- psi_matrix(frame$IncLevel1, 8L)
  r6_psi <- psi_matrix(frame$IncLevel2, 8L)
  frame$scramble_mean_inclusion <- rowMeans(scramble_psi, na.rm = TRUE)
  frame$r6_mean_inclusion <- rowMeans(r6_psi, na.rm = TRUE)
  frame$rmats_delta_scramble_minus_r6 <- as.numeric(
    frame$IncLevelDifference
  )
  frame$treatment_delta_r6_minus_scramble <-
    -frame$rmats_delta_scramble_minus_r6
  calculated <- frame$scramble_mean_inclusion - frame$r6_mean_inclusion
  stopifnot(max(
    abs(calculated - frame$rmats_delta_scramble_minus_r6),
    na.rm = TRUE
  ) < 0.0021)

  for (cell_line in c("S2", "S3")) {
    scramble_columns <- paste0(
      scramble_samples[[cell_line]],
      "_IncLevel"
    )
    r6_columns <- paste0(r6_samples[[cell_line]], "_IncLevel")
    stopifnot(
      all(scramble_columns %in% names(frame)),
      all(r6_columns %in% names(frame))
    )
    frame[[paste0(cell_line, "_scramble_mean")]] <- rowMeans(
      data.matrix(frame[, scramble_columns, drop = FALSE]),
      na.rm = TRUE
    )
    frame[[paste0(cell_line, "_r6_mean")]] <- rowMeans(
      data.matrix(frame[, r6_columns, drop = FALSE]),
      na.rm = TRUE
    )
    frame[[paste0(cell_line, "_treatment_delta")]] <-
      frame[[paste0(cell_line, "_r6_mean")]] -
      frame[[paste0(cell_line, "_scramble_mean")]]
  }
  frame
}

validate_named_treatment_order <- function(frames) {
  expected_order <- NULL
  for (event_type in event_types) {
    frame <- frames[[event_type]]
    named_columns <- names(frame)[
      grepl("^BULK-SAM-[0-9]+_IncLevel$", names(frame))
    ]
    stopifnot(length(named_columns) == 16L)
    sample_order <- sub("_IncLevel$", "", named_columns)
    if (is.null(expected_order)) {
      expected_order <- sample_order
    } else {
      stopifnot(identical(sample_order, expected_order))
    }
    stopifnot(
      isTRUE(all.equal(
        psi_matrix(frame$IncLevel1, 8L),
        data.matrix(frame[, named_columns[1:8], drop = FALSE]),
        tolerance = 1e-12,
        check.attributes = FALSE
      )),
      isTRUE(all.equal(
        psi_matrix(frame$IncLevel2, 8L),
        data.matrix(frame[, named_columns[9:16], drop = FALSE]),
        tolerance = 1e-12,
        check.attributes = FALSE
      ))
    )
  }
  expected_order
}

merge_event_type <- function(disease, treatment) {
  disease <- add_disease_metrics(disease)
  treatment <- add_treatment_metrics(treatment)
  disease_selected <- data.frame(
    event_type = disease$event_type,
    event_key = disease$event_key,
    gene_symbol = disease$gene_symbol,
    chromosome = disease$chr,
    strand = disease$strand,
    disease_event_id = disease$source_event_id,
    disease_fdr = disease$FDR,
    sma_mean_inclusion = disease$sma_mean_inclusion,
    control_mean_inclusion = disease$control_mean_inclusion,
    disease_delta_sma_minus_control =
      disease$disease_delta_sma_minus_control,
    stringsAsFactors = FALSE
  )
  treatment_selected <- data.frame(
    event_type = treatment$event_type,
    event_key = treatment$event_key,
    treatment_event_id = treatment$source_event_id,
    treatment_fdr = treatment$FDR,
    scramble_mean_inclusion = treatment$scramble_mean_inclusion,
    r6_mean_inclusion = treatment$r6_mean_inclusion,
    rmats_delta_scramble_minus_r6 =
      treatment$rmats_delta_scramble_minus_r6,
    treatment_delta_r6_minus_scramble =
      treatment$treatment_delta_r6_minus_scramble,
    S2_scramble_mean = treatment$S2_scramble_mean,
    S2_r6_mean = treatment$S2_r6_mean,
    S2_treatment_delta = treatment$S2_treatment_delta,
    S3_scramble_mean = treatment$S3_scramble_mean,
    S3_r6_mean = treatment$S3_r6_mean,
    S3_treatment_delta = treatment$S3_treatment_delta,
    stringsAsFactors = FALSE
  )
  common <- merge(
    disease_selected,
    treatment_selected,
    by = c("event_type", "event_key"),
    all = FALSE,
    sort = FALSE
  )
  common$pooled_distance_improvement <-
    abs(common$scramble_mean_inclusion - common$control_mean_inclusion) -
    abs(common$r6_mean_inclusion - common$control_mean_inclusion)
  for (cell_line in c("S2", "S3")) {
    common[[paste0(cell_line, "_distance_improvement")]] <-
      abs(
        common[[paste0(cell_line, "_scramble_mean")]] -
          common$control_mean_inclusion
      ) -
      abs(
        common[[paste0(cell_line, "_r6_mean")]] -
          common$control_mean_inclusion
      )
    common[[paste0(cell_line, "_direction_reversal")]] <-
      common$disease_delta_sma_minus_control *
      common[[paste0(cell_line, "_treatment_delta")]] < 0
  }
  common$biological_direction_reversal <-
    common$disease_delta_sma_minus_control *
    common$treatment_delta_r6_minus_scramble < 0
  common$corrected_rescue <-
    common$biological_direction_reversal &
    common$pooled_distance_improvement > 0
  common$strict_corrected_rescue <-
    common$corrected_rescue &
    common$S2_direction_reversal &
    common$S3_direction_reversal &
    common$S2_distance_improvement > 0 &
    common$S3_distance_improvement > 0
  common
}

audit_published <- function(disease_frames, treatment_frames) {
  records <- list()
  for (event_type in event_types) {
    published <- as.data.frame(
      suppressWarnings(read_excel(published_workbook, sheet = event_type)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(published)[names(published) == "ID"] <- "source_event_id"
    names(published)[names(published) == "geneSymbol"] <- "gene_symbol"
    published$gene_symbol <- normalize_gene_symbol(published$gene_symbol)
    published$published_key <- make_event_key(
      published,
      event_type,
      published_specs
    )

    disease_source <- disease_frames[[event_type]][, c(
      "source_event_id",
      "event_key",
      "IncLevelDifference",
      "gene_symbol"
    )]
    names(disease_source) <- c(
      "source_event_id",
      "disease_event_key",
      "disease_source_delta",
      "disease_source_gene"
    )
    treatment_source <- treatment_frames[[event_type]][, c(
      "source_event_id",
      "event_key",
      "IncLevelDifference",
      "gene_symbol"
    )]
    names(treatment_source) <- c(
      "source_event_id",
      "treatment_event_key",
      "treatment_source_delta_scramble_minus_r6",
      "treatment_source_gene"
    )

    disease_published <- published[
      published$comp == "SMA-CTRL",
      c(
        "source_event_id",
        "published_key",
        "gene_symbol",
        "IncLevelDifference"
      )
    ]
    names(disease_published)[4] <- "disease_published_delta"
    disease_published <- merge(
      disease_published,
      disease_source,
      by = "source_event_id",
      all.x = TRUE,
      sort = FALSE
    )
    treatment_published <- published[
      published$comp == "SMAmo-SMAscr",
      c(
        "source_event_id",
        "published_key",
        "gene_symbol",
        "IncLevelDifference"
      )
    ]
    names(treatment_published)[4] <- "treatment_published_delta"
    treatment_published <- merge(
      treatment_published,
      treatment_source,
      by = "source_event_id",
      all.x = TRUE,
      sort = FALSE
    )
    stopifnot(
      max(abs(
        disease_published$disease_published_delta -
          disease_published$disease_source_delta
      )) < 1e-12,
      max(abs(
        treatment_published$treatment_published_delta -
          treatment_published$treatment_source_delta_scramble_minus_r6
      )) < 1e-12,
      !anyDuplicated(disease_published$published_key),
      !anyDuplicated(treatment_published$published_key)
    )

    paired <- merge(
      disease_published,
      treatment_published,
      by = "published_key",
      suffixes = c("_disease", "_treatment"),
      all = FALSE,
      sort = FALSE
    )
    paired$event_type <- event_type
    paired$biological_treatment_delta_r6_minus_scramble <-
      -paired$treatment_source_delta_scramble_minus_r6
    paired$same_biological_direction <-
      paired$disease_source_delta *
      paired$biological_treatment_delta_r6_minus_scramble > 0
    paired$biological_direction_reversal <-
      paired$disease_source_delta *
      paired$biological_treatment_delta_r6_minus_scramble < 0
    paired$structurally_identical_full_event <-
      paired$disease_event_key == paired$treatment_event_key
    records[[event_type]] <- paired
  }
  do.call(rbind, records)
}

parse_smn_ratios <- function() {
  source <- as.data.frame(
    suppressMessages(suppressWarnings(read_excel(
        source_workbook,
        sheet = "Supp. Fig. 13",
        col_names = FALSE
      ))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  records <- list()
  record_index <- 1L
  group_offsets <- list(
    Untreated = 1:3,
    Scramble = 4:7,
    `R6-MO` = 8:11
  )
  source_character <- as.matrix(data.frame(lapply(source, as.character)))
  label_matrix <- source_character %in% c("S1", "S2", "S3")
  dim(label_matrix) <- dim(source_character)
  label_hits <- which(label_matrix, arr.ind = TRUE)
  candidate_indices <- which(vapply(
    seq_len(nrow(label_hits)),
    function(index) {
      row_index <- label_hits[index, "row"]
      column_index <- label_hits[index, "col"]
      if (column_index + 11L > ncol(source)) {
        return(FALSE)
      }
      values <- suppressWarnings(as.numeric(unlist(
        source[row_index, column_index + 1:11, drop = FALSE]
      )))
      sum(is.finite(values)) == 11L
    },
    logical(1)
  ))
  candidates <- label_hits[candidate_indices, , drop = FALSE]
  candidate_labels <- source_character[candidates]
  stopifnot(
    nrow(candidates) == 3L,
    setequal(candidate_labels, c("S1", "S2", "S3"))
  )
  for (candidate_index in seq_len(nrow(candidates))) {
    row_index <- candidates[candidate_index, "row"]
    column_index <- candidates[candidate_index, "col"]
    cell_line <- source_character[row_index, column_index]
    stopifnot(cell_line %in% c("S1", "S2", "S3"))
    for (group_name in names(group_offsets)) {
      group_columns <- column_index + group_offsets[[group_name]]
      values <- as.numeric(unlist(
        source[row_index, group_columns, drop = FALSE]
      ))
      values <- values[is.finite(values)]
      for (replicate in seq_along(values)) {
        records[[record_index]] <- data.frame(
          cell_line = cell_line,
          group = group_name,
          replicate = replicate,
          full_length_to_delta7_ratio = values[replicate],
          stringsAsFactors = FALSE
        )
        record_index <- record_index + 1L
      }
    }
  }
  do.call(rbind, records)
}

disease_frames <- read_event_sheets(disease_workbook)
treatment_frames <- read_event_sheets(treatment_workbook)
named_treatment_order <- validate_named_treatment_order(treatment_frames)

# Group sizes in the disease workbook uniquely match the six SMA and nine
# control untreated libraries. Treatment group identities are independently
# fixed by the named PSI columns validated above.
metadata <- read_sample_metadata("GSE290979")
stopifnot(
  sum(metadata$genotype == "SMA" & metadata$treatment == "NT") == 6L,
  sum(metadata$genotype == "CTRL" & metadata$treatment == "NT") == 9L,
  sum(metadata$treatment == "Scramble") == 8L,
  sum(metadata$treatment == "R6-Mo") == 8L,
  setequal(
    named_treatment_order[1:8],
    rownames(metadata)[metadata$treatment == "Scramble"]
  ),
  setequal(
    named_treatment_order[9:16],
    rownames(metadata)[metadata$treatment == "R6-Mo"]
  )
)

common_by_type <- lapply(event_types, function(event_type) {
  merge_event_type(
    disease_frames[[event_type]],
    treatment_frames[[event_type]]
  )
})
names(common_by_type) <- event_types
common <- do.call(rbind, common_by_type)
rownames(common) <- NULL
strict <- common[common$strict_corrected_rescue, , drop = FALSE]
published_audit <- audit_published(disease_frames, treatment_frames)
rownames(published_audit) <- NULL
smn_ratios <- parse_smn_ratios()

stopifnot(
  sum(vapply(disease_frames, nrow, integer(1))) == 10553L,
  sum(vapply(treatment_frames, nrow, integer(1))) == 2614L,
  nrow(common) == 182L,
  sum(common$biological_direction_reversal) == 112L,
  sum(common$corrected_rescue) == 108L,
  nrow(strict) == 83L,
  length(unique(strict$gene_symbol)) == 74L,
  all(c("SEPTIN10", "SEPTIN11") %in% strict$gene_symbol),
  !any(grepl("^[0-9]{5}(\\.0+)?$", strict$gene_symbol)),
  nrow(published_audit) == 111L,
  sum(published_audit$same_biological_direction) == 111L,
  sum(published_audit$structurally_identical_full_event) == 67L
)

strict_pairs <- paste(
  strict$event_type,
  strict$disease_event_id,
  strict$treatment_event_id,
  sep = "|"
)
published_pairs <- paste(
  published_audit$event_type,
  published_audit$source_event_id_disease,
  published_audit$source_event_id_treatment,
  sep = "|"
)
stopifnot(length(intersect(strict_pairs, published_pairs)) == 0L)

event_summary <- do.call(
  rbind,
  lapply(event_types, function(event_type) {
    subset <- common[common$event_type == event_type, , drop = FALSE]
    data.frame(
      event_type = event_type,
      exact_common = nrow(subset),
      biological_direction_reversal =
        sum(subset$biological_direction_reversal),
      pooled_corrected = sum(subset$corrected_rescue),
      strict_both_lines = sum(subset$strict_corrected_rescue)
    )
  })
)

gene_groups <- split(common, common$gene_symbol)
gene_summary <- do.call(
  rbind,
  lapply(names(gene_groups), function(gene_symbol) {
    group <- gene_groups[[gene_symbol]]
    strict_group <- group[group$strict_corrected_rescue, , drop = FALSE]
    data.frame(
      gene_symbol = gene_symbol,
      common_event_count = nrow(group),
      corrected_event_count = sum(group$corrected_rescue),
      strict_event_count = sum(group$strict_corrected_rescue),
      strict_event_types = paste(
        sort(unique(strict_group$event_type)),
        collapse = ";"
      ),
      max_abs_disease_delta = max(
        abs(group$disease_delta_sma_minus_control)
      ),
      max_abs_treatment_delta = max(
        abs(group$treatment_delta_r6_minus_scramble)
      ),
      stringsAsFactors = FALSE
    )
  })
)
gene_summary <- gene_summary[order(
  -gene_summary$strict_event_count,
  -gene_summary$max_abs_treatment_delta
), ]

ratio_means <- aggregate(
  full_length_to_delta7_ratio ~ cell_line + group,
  smn_ratios,
  mean
)
ratio_wide <- reshape(
  ratio_means,
  idvar = "cell_line",
  timevar = "group",
  direction = "wide"
)
ratio_wide$r6_over_scramble <-
  ratio_wide$`full_length_to_delta7_ratio.R6-MO` /
  ratio_wide$full_length_to_delta7_ratio.Scramble

write_tsv(
  common,
  "results/r/splicing/GSE290979_R_common_splicing_events.tsv"
)
write_tsv(
  strict,
  "results/r/splicing/GSE290979_R_strict_corrected_events.tsv"
)
write_tsv(
  gene_summary,
  "results/r/splicing/GSE290979_R_splicing_gene_summary.tsv"
)
write_tsv(
  published_audit,
  "results/r/splicing/GSE290979_R_published_sign_audit.tsv"
)
write_tsv(
  event_summary,
  "results/r/splicing/GSE290979_R_event_summary.tsv"
)
write_tsv(
  data.frame(
    rmats_position = seq_along(named_treatment_order),
    sample_id = named_treatment_order,
    group = rep(c("Scramble", "R6-MO"), each = 8L)
  ),
  "results/r/splicing/GSE290979_R_named_treatment_order.tsv"
)
write_tsv(
  smn_ratios,
  "results/r/splicing/GSE290979_R_SMN_full_length_delta7_ratio.tsv"
)
write_tsv(
  ratio_wide,
  "results/r/splicing/GSE290979_R_SMN_ratio_by_line.tsv"
)

cat(
  "GSE290979 R splicing audit complete:",
  nrow(common),
  "exact common,",
  sum(common$corrected_rescue),
  "corrected,",
  nrow(strict),
  "strict events across",
  length(unique(strict$gene_symbol)),
  "genes; published same-direction pairs =",
  sum(published_audit$same_biological_direction),
  "\n"
)
