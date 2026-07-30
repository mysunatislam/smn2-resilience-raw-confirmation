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

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is required; run r/00_install_packages.R")
}

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

marker_value <- function(path, metric) {
  marker <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  index <- match(metric, marker$metric)
  if (is.na(index)) {
    stop("Metric '", metric, "' is missing from ", path)
  }
  marker$value[index]
}

write_output <- function(frame, name, row_names = FALSE) {
  path <- file.path(output_root, name)
  write.table(
    frame,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = row_names,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

fit_edger <- function(counts, metadata, design, coefficient) {
  y <- edgeR::DGEList(counts = counts, samples = metadata)
  keep <- edgeR::filterByExpr(y, design = design)
  if (sum(keep) < 100L) {
    stop("Fewer than 100 genes passed edgeR filtering")
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  table$gene_id <- rownames(table)
  rownames(table) <- NULL
  list(y = y, table = table)
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
default_work_root <- if (
  .Platform$OS.type == "windows" && dir.exists("E:/")
) {
  "E:/smn2_gse290979_local9"
} else {
  file.path(ROOT, "local_work", "gse290979_local9")
}
work_root <- option_value(options, "work-root", default_work_root)
work_root <- normalizePath(work_root, winslash = "/", mustWork = TRUE)
output_root <- option_value(
  options,
  "output-root",
  file.path(ROOT, "results", "r", "raw_confirmation", "local9")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
count_root <- option_value(
  options,
  "count-root",
  file.path(work_root, "output")
)
count_root <- normalizePath(count_root, winslash = "/", mustWork = TRUE)
decision_marker <- option_value(
  options,
  "decision-marker",
  "BENCHMARK_DECISION.tsv"
)
count_file <- option_value(options, "count-file", "gene_counts.tsv")
quantification_method <- option_value(
  options,
  "quantification-method",
  "genome_alignment_featurecounts"
)
whole_genome_alignment <- identical(
  quantification_method,
  "genome_alignment_featurecounts"
)
if (
  basename(decision_marker) != decision_marker ||
    basename(count_file) != count_file
) {
  stop("--decision-marker and --count-file must be file names")
}

sample_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_local9_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_sheet <- sample_sheet[order(sample_sheet$array_index), , drop = FALSE]
stopifnot(
  nrow(sample_sheet) == 9L,
  length(unique(sample_sheet$sample_id)) == 9L,
  length(unique(sample_sheet$donor_line[sample_sheet$treatment == "NT"])) == 5L
)

decision_paths <- file.path(
  count_root,
  sample_sheet$sample_id,
  decision_marker
)
if (!all(file.exists(decision_paths))) {
  stop(
    "All nine sample decision markers are required; missing: ",
    paste(sample_sheet$sample_id[!file.exists(decision_paths)], collapse = ", ")
  )
}
decisions <- vapply(
  decision_paths,
  marker_value,
  character(1),
  metric = "status"
)
if (any(decisions != "PASS")) {
  stop(
    "Every local9 sample must be PASS; observed: ",
    paste(sample_sheet$sample_id, decisions, sep = "=", collapse = ", ")
  )
}

count_paths <- file.path(
  count_root,
  sample_sheet$sample_id,
  count_file
)
if (!all(file.exists(count_paths))) {
  stop(
    "All nine raw count files are required; missing: ",
    paste(sample_sheet$sample_id[!file.exists(count_paths)], collapse = ", ")
  )
}
count_tables <- lapply(count_paths, function(path) {
  table <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("gene_id", "count")
  if (!all(required %in% names(table))) {
    stop("Raw count file lacks required columns: ", path)
  }
  if (!"gene_name" %in% names(table)) {
    table$gene_name <- NA_character_
  }
  table <- table[c("gene_id", "gene_name", "count")]
  if (
    anyDuplicated(table$gene_id) ||
      any(!is.finite(table$count)) ||
      any(table$count < 0)
  ) {
    stop("Invalid gene identifiers or counts in ", path)
  }
  table
})
reference_gene_ids <- count_tables[[1L]]$gene_id
if (!all(vapply(
  count_tables,
  function(table) identical(table$gene_id, reference_gene_ids),
  logical(1)
))) {
  stop("Gene identifier order differs across raw count files")
}

counts <- do.call(cbind, lapply(count_tables, `[[`, "count"))
storage.mode(counts) <- "numeric"
rownames(counts) <- reference_gene_ids
colnames(counts) <- sample_sheet$sample_id
gene_names <- vapply(seq_along(reference_gene_ids), function(index) {
  values <- unique(na.omit(vapply(
    count_tables,
    function(table) table$gene_name[index],
    character(1)
  )))
  values <- values[nzchar(values)]
  if (length(values)) values[1L] else NA_character_
}, character(1))
gene_annotation <- data.frame(
  gene_id = reference_gene_ids,
  ensembl_gene_id = sub("[.][0-9]+$", "", reference_gene_ids),
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

raw_matrix <- cbind(
  gene_annotation,
  as.data.frame(counts, check.names = FALSE)
)
write_output(raw_matrix, "GSE290979_local9_raw_gene_counts.tsv")

disease_metadata <- sample_sheet[
  sample_sheet$treatment == "NT",
  ,
  drop = FALSE
]
disease_metadata$genotype <- factor(
  disease_metadata$genotype,
  levels = c("CTRL", "SMA")
)
disease_counts <- counts[, disease_metadata$sample_id, drop = FALSE]
disease_design <- model.matrix(~ genotype, data = disease_metadata)
rownames(disease_design) <- disease_metadata$sample_id
disease_fit <- fit_edger(
  disease_counts,
  disease_metadata,
  disease_design,
  "genotypeSMA"
)
disease_result <- merge(
  gene_annotation,
  disease_fit$table,
  by = "gene_id",
  all.y = TRUE,
  sort = FALSE
)
disease_result$contrast <- "SMA_minus_CTRL_untreated_donor_lines"
write_output(
  disease_result,
  "GSE290979_local9_SMA_vs_CTRL_edgeR.tsv"
)

disease_folds <- lapply(unique(disease_metadata$donor_line), function(omit) {
  keep_samples <- disease_metadata$donor_line != omit
  metadata <- droplevels(disease_metadata[keep_samples, , drop = FALSE])
  fold_counts <- disease_counts[, metadata$sample_id, drop = FALSE]
  design <- model.matrix(~ genotype, data = metadata)
  rownames(design) <- metadata$sample_id
  fit <- fit_edger(fold_counts, metadata, design, "genotypeSMA")
  data.frame(
    omitted_donor_line = omit,
    fit$table,
    stringsAsFactors = FALSE
  )
})
disease_lolo <- do.call(rbind, disease_folds)
rownames(disease_lolo) <- NULL
write_output(
  disease_lolo,
  "GSE290979_local9_disease_LOLO_long.tsv"
)
full_direction <- setNames(
  sign(disease_result$logFC),
  disease_result$gene_id
)
lolo_split <- split(disease_lolo, disease_lolo$gene_id)
disease_lolo_summary <- do.call(rbind, lapply(
  lolo_split,
  function(frame) {
    direction <- full_direction[frame$gene_id[1L]]
    data.frame(
      gene_id = frame$gene_id[1L],
      lolo_folds_tested = nrow(frame),
      lolo_median_logFC = stats::median(frame$logFC),
      lolo_min_logFC = min(frame$logFC),
      lolo_max_logFC = max(frame$logFC),
      lolo_same_direction_fraction = mean(sign(frame$logFC) == direction),
      lolo_all_same_direction = all(sign(frame$logFC) == direction),
      stringsAsFactors = FALSE
    )
  }
))
rownames(disease_lolo_summary) <- NULL
write_output(
  disease_lolo_summary,
  "GSE290979_local9_disease_LOLO_summary.tsv"
)

treatment_metadata <- sample_sheet[
  sample_sheet$treatment %in% c("Scramble", "R6-Mo"),
  ,
  drop = FALSE
]
treatment_metadata$donor_line <- factor(treatment_metadata$donor_line)
treatment_metadata$treatment <- factor(
  treatment_metadata$treatment,
  levels = c("Scramble", "R6-Mo")
)
treatment_counts <- counts[, treatment_metadata$sample_id, drop = FALSE]
treatment_design <- model.matrix(
  ~ donor_line + treatment,
  data = treatment_metadata
)
rownames(treatment_design) <- treatment_metadata$sample_id
treatment_fit <- fit_edger(
  treatment_counts,
  treatment_metadata,
  treatment_design,
  "treatmentR6-Mo"
)
treatment_result <- merge(
  gene_annotation,
  treatment_fit$table,
  by = "gene_id",
  all.y = TRUE,
  sort = FALSE
)
treatment_result$contrast <- "R6-Mo_minus_Scramble_paired_donor_lines"
write_output(
  treatment_result,
  "GSE290979_local9_R6_vs_scramble_paired_edgeR.tsv"
)

treatment_logcpm <- edgeR::cpm(
  treatment_fit$y,
  log = TRUE,
  prior.count = 2
)
linewise <- do.call(rbind, lapply(
  levels(treatment_metadata$donor_line),
  function(line) {
    line_metadata <- treatment_metadata[
      treatment_metadata$donor_line == line,
      ,
      drop = FALSE
    ]
    scramble <- line_metadata$sample_id[
      line_metadata$treatment == "Scramble"
    ]
    r6 <- line_metadata$sample_id[line_metadata$treatment == "R6-Mo"]
    stopifnot(length(scramble) == 1L, length(r6) == 1L)
    data.frame(
      donor_line = line,
      gene_id = rownames(treatment_logcpm),
      r6_minus_scramble_logCPM = (
        treatment_logcpm[, r6] - treatment_logcpm[, scramble]
      ),
      stringsAsFactors = FALSE
    )
  }
))
write_output(
  linewise,
  "GSE290979_local9_treatment_linewise_logCPM.tsv"
)
linewise_split <- split(linewise, linewise$gene_id)
treatment_linewise_summary <- do.call(rbind, lapply(
  linewise_split,
  function(frame) {
    effects <- frame$r6_minus_scramble_logCPM
    data.frame(
      gene_id = frame$gene_id[1L],
      donor_lines = nrow(frame),
      median_logCPM_effect = stats::median(effects),
      min_logCPM_effect = min(effects),
      max_logCPM_effect = max(effects),
      same_direction = length(unique(sign(effects))) == 1L,
      stringsAsFactors = FALSE
    )
  }
))
rownames(treatment_linewise_summary) <- NULL
write_output(
  treatment_linewise_summary,
  "GSE290979_local9_treatment_linewise_summary.tsv"
)

sample_qc <- data.frame(
  sample_sheet,
  raw_library_size = colSums(counts[, sample_sheet$sample_id, drop = FALSE]),
  stringsAsFactors = FALSE
)
write_output(sample_qc, "GSE290979_local9_sample_qc.tsv")
analysis_summary <- data.frame(
  metric = c(
    "libraries", "untreated_donor_lines", "control_donor_lines",
    "sma_donor_lines", "paired_treatment_donor_lines",
    "disease_lolo_folds", "treatment_lolo_estimable",
    "random_split_used", "all_sample_decisions_pass",
    "quantification_method", "whole_genome_alignment",
    "count_decision_marker"
  ),
  value = c(
    nrow(sample_sheet),
    nrow(disease_metadata),
    sum(disease_metadata$genotype == "CTRL"),
    sum(disease_metadata$genotype == "SMA"),
    length(unique(treatment_metadata$donor_line)),
    length(disease_folds),
    "FALSE_two_pairs_leave_one_pair_has_zero_residual_df",
    "FALSE",
    all(decisions == "PASS"),
    quantification_method,
    whole_genome_alignment,
    decision_marker
  ),
  stringsAsFactors = FALSE
)
write_output(analysis_summary, "GSE290979_local9_analysis_summary.tsv")
writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "GSE290979_local9_R_sessionInfo.txt")
)

message(
  "Completed local9 ", quantification_method,
  " analysis: 5 untreated donor lines, 5 disease LOLO folds, ",
  "and 2 paired treatment lines."
)
