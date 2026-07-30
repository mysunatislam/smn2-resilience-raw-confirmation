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

if (!requireNamespace("tximport", quietly = TRUE)) {
  stop("tximport is required; run r/00_install_packages.R")
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

write_tsv <- function(frame, path) {
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

read_gencode_headers <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  chunks <- list()
  chunk_index <- 0L
  repeat {
    lines <- readLines(connection, n = 100000L, warn = FALSE)
    if (!length(lines)) {
      break
    }
    headers <- lines[startsWith(lines, ">")]
    if (!length(headers)) {
      next
    }
    fields <- strsplit(sub("^>", "", headers), "|", fixed = TRUE)
    field_count <- lengths(fields)
    if (any(field_count < 6L)) {
      stop("Unexpected GENCODE transcript FASTA header format")
    }
    chunk_index <- chunk_index + 1L
    chunks[[chunk_index]] <- data.frame(
      transcript_id = vapply(fields, `[[`, character(1), 1L),
      gene_id = vapply(fields, `[[`, character(1), 2L),
      gene_name = vapply(fields, `[[`, character(1), 6L),
      stringsAsFactors = FALSE
    )
  }
  if (!length(chunks)) {
    stop("No transcript headers were read from ", path)
  }
  mapping <- unique(do.call(rbind, chunks))
  if (
    any(!nzchar(mapping$transcript_id)) ||
      any(!nzchar(mapping$gene_id)) ||
      anyDuplicated(mapping$transcript_id)
  ) {
    stop("GENCODE transcript-to-gene mapping is invalid")
  }
  mapping
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
default_work_root <- if (
  .Platform$OS.type == "windows" && dir.exists("E:/")
) {
  "E:/smn2_gse290979_local9"
} else {
  file.path(ROOT, "local_work", "gse290979_local9")
}
work_root <- normalizePath(
  option_value(options, "work-root", default_work_root),
  winslash = "/",
  mustWork = TRUE
)
quant_root <- normalizePath(
  option_value(
    options,
    "quant-root",
    file.path(work_root, "salmon", "quant")
  ),
  winslash = "/",
  mustWork = TRUE
)
count_root <- option_value(
  options,
  "count-root",
  file.path(work_root, "salmon", "gene_counts")
)
dir.create(count_root, recursive = TRUE, showWarnings = FALSE)
count_root <- normalizePath(count_root, winslash = "/", mustWork = TRUE)
transcriptome_path <- normalizePath(
  option_value(
    options,
    "transcriptome",
    file.path(
      work_root,
      "reference",
      "gencode_v47_salmon",
      "gencode.v47.transcripts.fa.gz"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)

sample_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_local9_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
sample_sheet <- sample_sheet[order(sample_sheet$array_index), , drop = FALSE]
stopifnot(
  nrow(sample_sheet) == 9L,
  !anyDuplicated(sample_sheet$sample_id)
)

quant_markers <- file.path(
  quant_root,
  sample_sheet$sample_id,
  "SALMON_QUANT_COMPLETE.tsv"
)
if (!all(file.exists(quant_markers))) {
  stop(
    "Salmon completion markers are missing: ",
    paste(
      sample_sheet$sample_id[!file.exists(quant_markers)],
      collapse = ", "
    )
  )
}
statuses <- vapply(
  quant_markers,
  marker_value,
  character(1),
  metric = "status"
)
if (any(statuses != "PASS")) {
  stop("Every Salmon quantification marker must have status PASS")
}
quant_paths <- vapply(
  quant_markers,
  marker_value,
  character(1),
  metric = "quant_path"
)
quant_files <- file.path(quant_paths, "quant.sf")
if (!all(file.exists(quant_files))) {
  stop(
    "Salmon quant.sf files are missing: ",
    paste(quant_files[!file.exists(quant_files)], collapse = ", ")
  )
}
names(quant_files) <- sample_sheet$sample_id

tx_annotation <- read_gencode_headers(transcriptome_path)
tx2gene <- unique(tx_annotation[c("transcript_id", "gene_id")])
txi <- tximport::tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)
counts <- txi$counts
if (
  !identical(colnames(counts), sample_sheet$sample_id) ||
    nrow(counts) < 1000L ||
    any(!is.finite(counts)) ||
    any(counts < 0)
) {
  stop("Imported Salmon gene-count matrix failed integrity checks")
}

gene_annotation <- tx_annotation[
  !duplicated(tx_annotation$gene_id),
  c("gene_id", "gene_name"),
  drop = FALSE
]
gene_names <- gene_annotation$gene_name[
  match(rownames(counts), gene_annotation$gene_id)
]
if (sum(!is.na(gene_names) & nzchar(gene_names)) < 1000L) {
  stop("Fewer than 1000 imported genes have GENCODE gene names")
}

for (sample_id in sample_sheet$sample_id) {
  sample_root <- file.path(count_root, sample_id)
  dir.create(sample_root, recursive = TRUE, showWarnings = FALSE)
  count_path <- file.path(sample_root, "gene_counts.tsv")
  write_tsv(
    data.frame(
      gene_id = rownames(counts),
      gene_name = gene_names,
      count = counts[, sample_id],
      stringsAsFactors = FALSE
    ),
    count_path
  )
  source_marker <- quant_markers[
    match(sample_id, sample_sheet$sample_id)
  ]
  marker <- data.frame(
    metric = c(
      "status", "sample_id", "quantification_method",
      "counts_from_abundance", "genes", "gene_count_sum",
      "source_quant_marker", "source_quant_file"
    ),
    value = c(
      "PASS", sample_id, "salmon_transcriptome_quasimapping",
      "lengthScaledTPM", nrow(counts), sum(counts[, sample_id]),
      source_marker, quant_files[[sample_id]]
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    marker,
    file.path(sample_root, "SALMON_QUANT_COMPLETE.tsv")
  )
}

write_tsv(
  data.frame(
    metric = c(
      "status", "libraries", "genes", "quantification_method",
      "counts_from_abundance", "whole_genome_alignment",
      "random_split_used", "transcriptome"
    ),
    value = c(
      "PASS", ncol(counts), nrow(counts),
      "salmon_transcriptome_quasimapping", "lengthScaledTPM",
      "FALSE", "FALSE", transcriptome_path
    ),
    stringsAsFactors = FALSE
  ),
  file.path(count_root, "SALMON_GENE_IMPORT_COMPLETE.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(count_root, "SALMON_GENE_IMPORT_R_sessionInfo.txt")
)
message(
  "Imported ", ncol(counts), " Salmon libraries into ",
  nrow(counts), " GENCODE gene rows with lengthScaledTPM counts."
)
