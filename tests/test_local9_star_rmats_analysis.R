script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) != 1L) {
  stop("This test must be run with Rscript")
}
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

frozen <- read.delim(
  file.path(
    ROOT,
    "config",
    "frozen_83_splice_events_2026-08-01.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(nrow(frozen) == 83L)

scratch <- tempfile("local9_star_rmats_test_")
rmats_root <- file.path(scratch, "rmats")
output_root <- file.path(scratch, "analysis")
dir.create(file.path(rmats_root, "disease"), recursive = TRUE)
dir.create(file.path(rmats_root, "treatment"), recursive = TRUE)

event_types <- c(ES = "SE", A5SS = "A5SS", A3SS = "A3SS", MXE = "MXE", RI = "RI")
fixed_root <- file.path(
  ROOT,
  "config",
  "rmats",
  "GSE290979",
  "fixed_events"
)
gene_mapping <- read.delim(
  file.path(
    ROOT,
    "config",
    "rmats",
    "GSE290979",
    "fixed_event_gene_id_mapping.tsv"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

format_values <- function(matrix) {
  apply(matrix, 1L, paste, collapse = ",")
}

for (event_type in names(event_types)) {
  rmats_type <- event_types[[event_type]]
  frame <- read.delim(
    file.path(fixed_root, paste0("fromGTF.", rmats_type, ".txt")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  event_keys <- frozen$event_key[frozen$event_type == event_type]
  frozen_type <- frozen[match(event_keys, frozen$event_key), , drop = FALSE]
  stopifnot(nrow(frame) == nrow(frozen_type))

  disease_sign <- sign(frozen_type$disease_delta_sma_minus_control)
  control <- matrix(0.50, nrow(frame), 3L)
  sma <- matrix(0.50 + 0.20 * disease_sign, nrow(frame), 2L)
  scramble <- matrix(0.50 + 0.20 * disease_sign, nrow(frame), 2L)
  r6 <- matrix(0.50, nrow(frame), 2L)

  add_columns <- function(source, group1, group2) {
    source$IJC_SAMPLE_1 <- format_values(matrix(
      12L,
      nrow(source),
      ncol(group1)
    ))
    source$SJC_SAMPLE_1 <- format_values(matrix(
      8L,
      nrow(source),
      ncol(group1)
    ))
    source$IJC_SAMPLE_2 <- format_values(matrix(
      12L,
      nrow(source),
      ncol(group2)
    ))
    source$SJC_SAMPLE_2 <- format_values(matrix(
      8L,
      nrow(source),
      ncol(group2)
    ))
    source$IncFormLen <- 300L
    source$SkipFormLen <- 149L
    source$PValue <- 1e-8
    source$FDR <- 1e-7
    source$IncLevel1 <- format_values(group1)
    source$IncLevel2 <- format_values(group2)
    source$IncLevelDifference <- rowMeans(group1) - rowMeans(group2)
    source
  }

  disease <- add_columns(frame, sma, control)
  treatment <- add_columns(frame, r6, scramble)
  mapping_type <- gene_mapping[
    gene_mapping$event_type == event_type,
    ,
    drop = FALSE
  ]
  mapping_index <- match(
    as.character(disease[[1L]]),
    as.character(mapping_type$event_id)
  )
  current_symbols <- mapping_type$gencode_gene_name[mapping_index]
  missing_current_symbols <- is.na(current_symbols) | !nzchar(current_symbols)
  current_symbols[missing_current_symbols] <-
    disease$geneSymbol[missing_current_symbols]
  disease$geneSymbol <- paste0('"', current_symbols, '"')
  treatment$geneSymbol <- paste0('"', current_symbols, '"')
  disease$GeneID <- paste0('"', disease$GeneID, '"')
  treatment$GeneID <- paste0('"', treatment$GeneID, '"')
  write.table(
    disease,
    file.path(
      rmats_root,
      "disease",
      paste0(rmats_type, ".MATS.JC.txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  write.table(
    treatment,
    file.path(
      rmats_root,
      "treatment",
      paste0(rmats_type, ".MATS.JC.txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
status <- system2(
  rscript,
  c(
    file.path(ROOT, "r", "30_analyze_local9_star_rmats.R"),
    paste0("--rmats-root=", rmats_root),
    paste0("--output-root=", output_root)
  )
)
stopifnot(status == 0L)

results <- read.delim(
  file.path(output_root, "frozen_83_raw_splice_confirmation.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
panel <- read.delim(
  file.path(output_root, "frozen_12_primary_raw_splice_confirmation.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
summary <- read.delim(
  file.path(output_root, "raw_splice_confirmation_summary.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
summary_value <- setNames(summary$value, summary$metric)

stopifnot(
  nrow(results) == 83L,
  nrow(panel) == 12L,
  identical(results$event_key, frozen$event_key),
  all(results$structurally_recovered),
  all(results$adequate_junction_support),
  all(results$disease_direction_reproduced),
  all(results$treatment_reversal_reproduced),
  all(results$both_lines_corrected),
  all(results$strong_raw_confirmation),
  all(nzchar(results$disease_raw_event_id)),
  all(nzchar(results$treatment_raw_event_id)),
  all(is.finite(results$disease_rmats_fdr_target_locus_discovery)),
  all(is.finite(results$treatment_rmats_fdr_target_locus_discovery)),
  all(results$raw_confirmation_limiting_reason ==
    "strong_raw_confirmation"),
  as.integer(summary_value[["structurally_recovered_83"]]) == 83L,
  as.integer(summary_value[["strong_raw_confirmation_83"]]) == 83L,
  as.integer(summary_value[["primary_events_attempted"]]) == 12L,
  as.integer(summary_value[["primary_structurally_recovered"]]) == 12L,
  as.logical(summary_value[["primary_panel_success"]]),
  file.exists(
    file.path(output_root, "RAW_RMATS_TABLE_ANALYSIS_COMPLETE.tsv")
  ),
  file.info(
    file.path(output_root, "raw_vs_processed_delta_psi.pdf")
  )$size > 0L
)

cat("local9 STAR/rMATS synthetic analysis test passed\n")

empty_class_root <- file.path(scratch, "analysis_empty_class")
for (contrast in c("disease", "treatment")) {
  path <- file.path(rmats_root, contrast, "A3SS.MATS.JC.txt")
  frame <- read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = ""
  )
  write.table(
    frame[0, , drop = FALSE],
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}
status <- system2(
  rscript,
  c(
    file.path(ROOT, "r", "30_analyze_local9_star_rmats.R"),
    paste0("--rmats-root=", rmats_root),
    paste0("--output-root=", empty_class_root)
  )
)
stopifnot(status == 0L)
empty_class_results <- read.delim(
  file.path(empty_class_root, "frozen_83_raw_splice_confirmation.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(empty_class_results) == 83L,
  !any(empty_class_results$structurally_recovered[
    empty_class_results$event_type == "A3SS"
  ]),
  all(empty_class_results$structurally_recovered[
    empty_class_results$event_type != "A3SS"
  ])
)

cat("local9 STAR/rMATS empty-class regression test passed\n")
