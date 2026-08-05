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
manuscript_root <- file.path(ROOT, "manuscript")

read_table <- function(filename) {
  path <- file.path(manuscript_root, filename)
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    stop("Missing or empty publication output: ", filename)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

results <- read_table("supplementary_table_S4_raw_splice_confirmation.tsv")
primary <- read_table("supplementary_table_S5_primary_raw_splice_panel.tsv")
reasons <- read_table(
  "supplementary_table_S6_raw_splice_unconfirmed_reasons.tsv"
)
support <- read_table("supplementary_table_S7_raw_junction_support.tsv")
summary <- read_table("supplementary_table_S8_raw_splice_summary.tsv")
provenance <- read_table(
  "supplementary_table_S9_raw_splice_provenance.tsv"
)
manifest <- read_table("sashimi_manifest.tsv")

summary_value <- setNames(summary$value, summary$metric)
provenance_value <- setNames(provenance$value, provenance$metric)
expected_recovered <- c(
  "OTUD3", "TXNDC11", "KIZ", "EIF4G3", "COL5A2", "ASAP1"
)

stopifnot(
  nrow(results) == 83L,
  nrow(primary) == 12L,
  nrow(reasons) == 83L,
  nrow(support) == 63L,
  nrow(summary) == 18L,
  nrow(provenance) == 12L,
  nrow(manifest) == 12L,
  identical(
    results$gene_symbol[results$structurally_recovered],
    expected_recovered
  ),
  sum(results$adequate_junction_support) == 5L,
  sum(results$disease_direction_reproduced) == 4L,
  sum(results$treatment_reversal_reproduced) == 3L,
  sum(results$both_lines_corrected) == 0L,
  sum(results$strong_raw_confirmation) == 0L,
  sum(primary$structurally_recovered) == 1L,
  identical(
    primary$gene_symbol[primary$structurally_recovered],
    "COL5A2"
  ),
  as.integer(summary_value[["primary_events_attempted"]]) == 12L,
  as.integer(summary_value[["structurally_recovered_83"]]) == 6L,
  as.integer(summary_value[["primary_structurally_recovered"]]) == 1L,
  as.integer(summary_value[["strong_raw_confirmation_83"]]) == 0L,
  identical(
    provenance_value[["rmats_fdr_scope"]],
    "rMATS_target_locus_de_novo_not_genome_wide"
  ),
  all(manifest$plot_event_source == "frozen_plot_only_definition"),
  identical(manifest$panel_order, seq_len(12L))
)

pdf_paths <- file.path(manuscript_root, manifest$pdf_path)
png_paths <- sub("\\.pdf$", ".png", pdf_paths)
figure_3 <- file.path(
  manuscript_root,
  "figures",
  c(
    "figure_3_raw_splice_confirmation.pdf",
    "figure_3_raw_splice_confirmation.png"
  )
)
all_artifacts <- c(pdf_paths, png_paths, figure_3)
stopifnot(
  all(file.exists(all_artifacts)),
  all(file.info(all_artifacts)$size > 0L)
)

pdf_signature <- charToRaw("%PDF")
valid_pdf <- vapply(
  c(pdf_paths, figure_3[1L]),
  function(path) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    identical(readBin(connection, what = "raw", n = 4L), pdf_signature)
  },
  logical(1)
)
png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
valid_png <- vapply(
  c(png_paths, figure_3[2L]),
  function(path) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    identical(readBin(connection, what = "raw", n = 8L), png_signature)
  },
  logical(1)
)
stopifnot(all(valid_pdf), all(valid_png))

manuscript <- paste(
  readLines(file.path(manuscript_root, "manuscript.md"), warn = FALSE),
  collapse = "\n"
)
stopifnot(
  grepl("six of[[:space:]]+83 frozen structures", manuscript),
  grepl("none met the full", manuscript, ignore.case = TRUE),
  grepl("Supplementary Figures S12-S23", manuscript, fixed = TRUE)
)

cat("Raw splice publication package test passed\n")
