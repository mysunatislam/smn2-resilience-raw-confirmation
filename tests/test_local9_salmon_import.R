script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))
import_script <- file.path(
  ROOT,
  "r",
  "27_import_gse290979_local9_salmon.R"
)
salmon_runner_path <- file.path(
  ROOT,
  "local",
  "gse290979_raw",
  "run_local9_salmon.ps1"
)
stopifnot(length(parse(file = import_script)) > 0L)
script_lines <- readLines(import_script, warn = FALSE)
salmon_runner <- readLines(salmon_runner_path, warn = FALSE)
required_patterns <- c(
  "tximport::tximport",
  "countsFromAbundance = \"lengthScaledTPM\"",
  "salmon_transcriptome_quasimapping",
  "whole_genome_alignment",
  "SALMON_QUANT_COMPLETE.tsv"
)
stopifnot(all(vapply(
  required_patterns,
  function(pattern) any(grepl(pattern, script_lines, fixed = TRUE)),
  logical(1)
)))
verified_quant_patterns <- c(
  '"quantify-verified"',
  "param([switch]$VerifiedOnly)",
  "Invoke-QuantifyPhase -VerifiedOnly",
  "Skipping $SampleId because its FASTQs are not yet",
  "Get-SalmonQuantMutexName",
  "$QuantMutex.WaitOne(0)",
  "another Salmon worker",
  "Test-CompletedSalmonQuantification"
)
stopifnot(all(vapply(
  verified_quant_patterns,
  function(pattern) any(grepl(pattern, salmon_runner, fixed = TRUE)),
  logical(1)
)))

sample_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_local9_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fixture_root <- tempfile("local9_salmon_fixture_")
quant_root <- file.path(fixture_root, "salmon", "quant")
count_root <- file.path(fixture_root, "salmon", "gene_counts")
reference_root <- file.path(
  fixture_root,
  "reference",
  "gencode_v47_salmon"
)
dir.create(quant_root, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_root, recursive = TRUE, showWarnings = FALSE)

gene_count <- 1200L
transcripts_per_gene <- 2L
transcript_count <- gene_count * transcripts_per_gene
gene_index <- rep(seq_len(gene_count), each = transcripts_per_gene)
transcript_index <- seq_len(transcript_count)
gene_ids <- sprintf("ENSG%011d.1", gene_index)
transcript_ids <- sprintf("ENST%011d.1", transcript_index)
transcript_names <- sprintf("TRANSCRIPT%04d", transcript_index)
gene_names <- sprintf("GENE%03d", gene_index)
headers <- paste0(
  ">", transcript_ids, "|", gene_ids,
  "|OTTHUMT", transcript_index,
  "|HAVANAT", transcript_index,
  "|", transcript_names,
  "|", gene_names,
  "|", nchar(rep("A", transcript_count)),
  "|1"
)
fasta_lines <- as.vector(rbind(headers, rep("A", transcript_count)))
transcriptome_path <- file.path(
  reference_root,
  "gencode.v47.transcripts.fa.gz"
)
connection <- gzfile(transcriptome_path, open = "wt")
writeLines(fasta_lines, connection)
close(connection)

write_tsv <- function(frame, path) {
  write.table(
    frame,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
for (sample_index in seq_len(nrow(sample_sheet))) {
  sample_id <- sample_sheet$sample_id[sample_index]
  sample_root <- file.path(quant_root, sample_id)
  quant_path <- file.path(sample_root, "salmon_quant.attempt-001")
  dir.create(quant_path, recursive = TRUE, showWarnings = FALSE)
  num_reads <- (
    40 + transcript_index + ((sample_index * transcript_index) %% 17)
  )
  write_tsv(
    data.frame(
      Name = transcript_ids,
      Length = 1000,
      EffectiveLength = 850,
      TPM = num_reads / sum(num_reads) * 1e6,
      NumReads = num_reads,
      stringsAsFactors = FALSE
    ),
    file.path(quant_path, "quant.sf")
  )
  write_tsv(
    data.frame(
      metric = c("status", "quant_path"),
      value = c(
        "PASS",
        normalizePath(quant_path, winslash = "/", mustWork = TRUE)
      ),
      stringsAsFactors = FALSE
    ),
    file.path(sample_root, "SALMON_QUANT_COMPLETE.tsv")
  )
}

output <- suppressWarnings(system2(
  file.path(R.home("bin"), "Rscript.exe"),
  c(
    shQuote(import_script),
    paste0(
      "--work-root=",
      normalizePath(fixture_root, winslash = "/", mustWork = TRUE)
    ),
    paste0(
      "--quant-root=",
      normalizePath(quant_root, winslash = "/", mustWork = TRUE)
    ),
    paste0("--count-root=", normalizePath(
      count_root,
      winslash = "/",
      mustWork = FALSE
    )),
    paste0(
      "--transcriptome=",
      normalizePath(transcriptome_path, winslash = "/", mustWork = TRUE)
    )
  ),
  stdout = TRUE,
  stderr = TRUE
))
status <- attr(output, "status")
if (!is.null(status)) {
  stop(
    "Synthetic Salmon import failed with status ",
    status,
    ":\n",
    paste(output, collapse = "\n")
  )
}

summary_path <- file.path(count_root, "SALMON_GENE_IMPORT_COMPLETE.tsv")
stopifnot(file.exists(summary_path))
summary <- read.delim(summary_path, stringsAsFactors = FALSE)
summary_value <- function(metric) {
  summary$value[match(metric, summary$metric)]
}
stopifnot(
  summary_value("status") == "PASS",
  summary_value("libraries") == "9",
  summary_value("genes") == as.character(gene_count),
  summary_value("quantification_method") ==
    "salmon_transcriptome_quasimapping",
  summary_value("counts_from_abundance") == "lengthScaledTPM",
  summary_value("whole_genome_alignment") == "FALSE",
  summary_value("random_split_used") == "FALSE"
)
for (sample_id in sample_sheet$sample_id) {
  count_path <- file.path(count_root, sample_id, "gene_counts.tsv")
  marker_path <- file.path(
    count_root,
    sample_id,
    "SALMON_QUANT_COMPLETE.tsv"
  )
  stopifnot(file.exists(count_path), file.exists(marker_path))
  counts <- read.delim(count_path, stringsAsFactors = FALSE)
  marker <- read.delim(marker_path, stringsAsFactors = FALSE)
  stopifnot(
    nrow(counts) == gene_count,
    !anyDuplicated(counts$gene_id),
    identical(as.character(counts$gene_name), unique(gene_names)),
    all(is.finite(counts$count)),
    all(counts$count >= 0),
    marker$value[marker$metric == "status"] == "PASS"
  )
}

unlink(fixture_root, recursive = TRUE, force = TRUE)
cat("local9 Salmon tximport tests passed\n")
