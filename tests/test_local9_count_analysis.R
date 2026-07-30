script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

analysis_script <- file.path(
  ROOT,
  "r",
  "24_analyze_gse290979_local9_counts.R"
)
stopifnot(length(parse(file = analysis_script)) > 0L)
analysis_lines <- readLines(analysis_script, warn = FALSE)
required_patterns <- c(
  "edgeR::glmQLFit",
  "edgeR::glmQLFTest",
  "SMA_minus_CTRL_untreated_donor_lines",
  "R6-Mo_minus_Scramble_paired_donor_lines",
  "disease_LOLO",
  "treatment_linewise",
  "FALSE_two_pairs_leave_one_pair_has_zero_residual_df",
  "random_split_used",
  "quantification_method",
  "whole_genome_alignment"
)
stopifnot(all(vapply(
  required_patterns,
  function(pattern) any(grepl(pattern, analysis_lines, fixed = TRUE)),
  logical(1)
)))
forbidden_patterns <- c(
  "(^|[^A-Za-z0-9_])sample[[:space:]]*[(]",
  "(^|[^A-Za-z0-9_])runif[[:space:]]*[(]",
  "(^|[^A-Za-z0-9_])rnorm[[:space:]]*[(]",
  "(^|[^A-Za-z0-9_.])set[.]seed[[:space:]]*[(]"
)
stopifnot(!any(vapply(
  forbidden_patterns,
  function(pattern) any(grepl(pattern, analysis_lines)),
  logical(1)
)))

sample_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_local9_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fixture_root <- tempfile("local9_count_fixture_")
output_root <- tempfile("local9_count_output_")
dir.create(fixture_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

gene_count <- 240L
gene_index <- seq_len(gene_count)
gene_ids <- sprintf("ENSG%011d.1", gene_index)
gene_names <- sprintf("GENE%03d", gene_index)
for (sample_index in seq_len(nrow(sample_sheet))) {
  sample <- sample_sheet[sample_index, , drop = FALSE]
  sample_output <- file.path(
    fixture_root,
    "output",
    sample$sample_id
  )
  dir.create(sample_output, recursive = TRUE, showWarnings = FALSE)
  write.table(
    data.frame(
      metric = "status",
      value = "PASS",
      stringsAsFactors = FALSE
    ),
    file.path(sample_output, "BENCHMARK_DECISION.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  counts <- (
    80L +
      3L * gene_index +
      ((gene_index * sample_index) %% 23L)
  )
  if (sample$treatment == "NT" && sample$genotype == "SMA") {
    counts[gene_index <= 30L] <- counts[gene_index <= 30L] + 90L
  }
  if (sample$treatment == "R6-Mo") {
    counts[gene_index > 30L & gene_index <= 60L] <-
      counts[gene_index > 30L & gene_index <= 60L] + 75L
  }
  write.table(
    data.frame(
      gene_id = gene_ids,
      gene_name = gene_names,
      count = counts,
      stringsAsFactors = FALSE
    ),
    file.path(sample_output, "gene_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

output <- suppressWarnings(system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(
    shQuote(analysis_script),
    paste0(
      "--work-root=",
      normalizePath(fixture_root, winslash = "/", mustWork = TRUE)
    ),
    paste0(
      "--output-root=",
      normalizePath(output_root, winslash = "/", mustWork = TRUE)
    )
  ),
  stdout = TRUE,
  stderr = TRUE
))
status <- attr(output, "status")
if (!is.null(status)) {
  stop(
    "Synthetic local9 count analysis failed with status ",
    status,
    ":\n",
    paste(output, collapse = "\n")
  )
}

expected_outputs <- c(
  "GSE290979_local9_raw_gene_counts.tsv",
  "GSE290979_local9_SMA_vs_CTRL_edgeR.tsv",
  "GSE290979_local9_disease_LOLO_long.tsv",
  "GSE290979_local9_disease_LOLO_summary.tsv",
  "GSE290979_local9_R6_vs_scramble_paired_edgeR.tsv",
  "GSE290979_local9_treatment_linewise_logCPM.tsv",
  "GSE290979_local9_treatment_linewise_summary.tsv",
  "GSE290979_local9_sample_qc.tsv",
  "GSE290979_local9_analysis_summary.tsv",
  "GSE290979_local9_R_sessionInfo.txt"
)
stopifnot(all(file.exists(file.path(output_root, expected_outputs))))
summary <- read.delim(
  file.path(output_root, "GSE290979_local9_analysis_summary.tsv"),
  stringsAsFactors = FALSE
)
summary_value <- function(metric) {
  summary$value[match(metric, summary$metric)]
}
stopifnot(
  summary_value("libraries") == "9",
  summary_value("untreated_donor_lines") == "5",
  summary_value("paired_treatment_donor_lines") == "2",
  summary_value("disease_lolo_folds") == "5",
  summary_value("random_split_used") == "FALSE",
  summary_value("all_sample_decisions_pass") == "TRUE",
  summary_value("quantification_method") ==
    "genome_alignment_featurecounts",
  summary_value("whole_genome_alignment") == "TRUE"
)
lolo <- read.delim(
  file.path(output_root, "GSE290979_local9_disease_LOLO_long.tsv"),
  stringsAsFactors = FALSE
)
linewise <- read.delim(
  file.path(
    output_root,
    "GSE290979_local9_treatment_linewise_logCPM.tsv"
  ),
  stringsAsFactors = FALSE
)
stopifnot(
  length(unique(lolo$omitted_donor_line)) == 5L,
  length(unique(linewise$donor_line)) == 2L
)

unlink(fixture_root, recursive = TRUE, force = TRUE)
unlink(output_root, recursive = TRUE, force = TRUE)
cat("local9 deterministic raw-count analysis tests passed\n")
