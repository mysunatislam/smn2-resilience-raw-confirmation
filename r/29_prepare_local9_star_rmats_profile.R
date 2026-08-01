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
check_only <- "--check-only" %in% arguments
argument_value <- function(prefix, default) {
  matches <- arguments[startsWith(arguments, prefix)]
  if (length(matches) == 0L) {
    return(default)
  }
  if (length(matches) != 1L) {
    stop("Expected one argument beginning with ", prefix)
  }
  sub(prefix, "", matches, fixed = TRUE)
}

output_directory <- argument_value(
  "--output-dir=",
  file.path(ROOT, "config", "rmats", "GSE290979", "local9")
)
sample_sheet_path <- file.path(
  ROOT,
  "config",
  "GSE290979_local9_sample_sheet.tsv"
)
fixed_event_directory <- file.path(
  ROOT,
  "config",
  "rmats",
  "GSE290979",
  "fixed_events"
)
frozen_event_path <- file.path(
  ROOT,
  "config",
  "frozen_83_splice_events_2026-08-01.tsv"
)

samples <- read.delim(
  sample_sheet_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
frozen_events <- read.delim(
  frozen_event_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(samples) == 9L,
  identical(samples$array_index, seq_len(9L)),
  !anyDuplicated(samples$sample_id),
  !anyDuplicated(samples$run_accession),
  all(samples$read_length == 151L),
  nrow(frozen_events) == 83L,
  !anyDuplicated(frozen_events$event_key)
)

line_order <- c("C1", "C2", "C3", "S2", "S3")
samples <- samples[order(
  match(samples$analysis_group, c(
    "CTRL_NT", "SMA_NT", "SMA_Scramble", "SMA_R6-Mo"
  )),
  match(samples$donor_line, line_order)
), , drop = FALSE]
rownames(samples) <- NULL

sample_lists <- list(
  disease_group1_sma = samples$sample_id[
    samples$genotype == "SMA" & samples$treatment == "NT"
  ],
  disease_group2_control = samples$sample_id[
    samples$genotype == "CTRL" & samples$treatment == "NT"
  ],
  treatment_group1_r6 = samples$sample_id[
    samples$treatment == "R6-Mo"
  ],
  treatment_group2_scramble = samples$sample_id[
    samples$treatment == "Scramble"
  ]
)
stopifnot(
  identical(
    sample_lists$disease_group1_sma,
    c("BULK-SAM-160", "BULK-SAM-148")
  ),
  identical(
    sample_lists$disease_group2_control,
    c("BULK-SAM-109", "BULK-SAM-104", "BULK-SAM-107")
  ),
  identical(
    sample_lists$treatment_group1_r6,
    c("BULK-SAM-166", "BULK-SAM-157")
  ),
  identical(
    sample_lists$treatment_group2_scramble,
    c("BULK-SAM-163", "BULK-SAM-151")
  )
)

coordinate_columns <- list(
  ES = c(
    "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ),
  A5SS = c(
    "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
    "flankingES", "flankingEE"
  ),
  A3SS = c(
    "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
    "flankingES", "flankingEE"
  ),
  MXE = c(
    "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
    "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES",
    "downstreamEE"
  ),
  RI = c(
    "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  )
)
file_event_type <- c(ES = "SE", A5SS = "A5SS", A3SS = "A3SS", MXE = "MXE", RI = "RI")
padding <- 1000L

event_intervals <- do.call(rbind, lapply(
  names(coordinate_columns),
  function(event_type) {
    path <- file.path(
      fixed_event_directory,
      paste0("fromGTF.", file_event_type[[event_type]], ".txt")
    )
    frame <- read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    columns <- coordinate_columns[[event_type]]
    stopifnot(all(c("chr", columns) %in% names(frame)))
    coordinates <- as.matrix(frame[, columns, drop = FALSE])
    storage.mode(coordinates) <- "numeric"
    data.frame(
      chromosome = frame$chr,
      start = pmax(0L, apply(coordinates, 1L, min) - padding),
      end = apply(coordinates, 1L, max) + padding,
      stringsAsFactors = FALSE
    )
  }
))
stopifnot(nrow(event_intervals) == 83L)
event_intervals <- event_intervals[order(
  event_intervals$chromosome,
  event_intervals$start,
  event_intervals$end
), , drop = FALSE]

merge_intervals <- function(frame) {
  merged <- vector("list", 0L)
  for (index in seq_len(nrow(frame))) {
    current <- frame[index, , drop = FALSE]
    previous <- length(merged)
    if (
      previous > 0L &&
        merged[[previous]]$chromosome == current$chromosome &&
        current$start <= merged[[previous]]$end
    ) {
      merged[[previous]]$end <- max(merged[[previous]]$end, current$end)
    } else {
      merged[[previous + 1L]] <- current
    }
  }
  result <- do.call(rbind, merged)
  rownames(result) <- NULL
  result
}
target_intervals <- merge_intervals(event_intervals)
stopifnot(
  nrow(target_intervals) <= 83L,
  all(target_intervals$start >= 0L),
  all(target_intervals$end > target_intervals$start)
)

sample_order <- samples[, c(
  "sample_id", "run_accession", "donor_line", "genotype", "treatment",
  "analysis_group", "fastq_r1_md5", "fastq_r2_md5"
)]
profile_summary <- data.frame(
  metric = c(
    "profile", "freeze_date", "libraries", "compressed_fastq_gib",
    "disease_group1_sma", "disease_group2_control",
    "treatment_group1_r6", "treatment_group2_scramble",
    "fixed_events", "primary_events", "merged_target_intervals",
    "target_padding_bp", "reference", "star_version", "star_twopass_mode",
    "star_read_map_number", "rmats_version",
    "rmats2sashimiplot_version",
    "bam_retention", "fastq_deletion", "bam_deletion", "inference_role"
  ),
  value = c(
    "local9_star_rmats", "2026-08-01", nrow(samples),
    sprintf("%.3f", sum(samples$fastq_bytes_total) / 1024^3),
    length(sample_lists$disease_group1_sma),
    length(sample_lists$disease_group2_control),
    length(sample_lists$treatment_group1_r6),
    length(sample_lists$treatment_group2_scramble),
    nrow(frozen_events), 12L, nrow(target_intervals), padding,
    "GENCODE_v47_GRCh38_primary_assembly",
    "2.7.10a", "None", "-1", "4.3.0", "4.0.0",
    "streamed_event_locus_coordinate_sorted_bam",
    "FALSE", "FALSE",
    "prospective_targeted_raw_splice_junction_confirmation"
  ),
  stringsAsFactors = FALSE
)

write_tsv <- function(frame, path, column_names = TRUE) {
  write.table(
    frame,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = column_names,
    na = ""
  )
}

assert_profile <- function() {
  required <- c(
    "local9_sample_order.tsv",
    "fixed_events_padded_1kb_merged.bed",
    "profile_summary.tsv",
    paste0(names(sample_lists), ".txt")
  )
  paths <- file.path(output_directory, required)
  if (!all(file.exists(paths))) {
    stop(
      "Generated local9 STAR/rMATS profile files are missing: ",
      paste(required[!file.exists(paths)], collapse = ", ")
    )
  }
  observed_order <- read.delim(
    file.path(output_directory, "local9_sample_order.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  observed_bed <- read.delim(
    file.path(output_directory, "fixed_events_padded_1kb_merged.bed"),
    header = FALSE,
    col.names = names(target_intervals),
    stringsAsFactors = FALSE
  )
  observed_summary <- read.delim(
    file.path(output_directory, "profile_summary.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stopifnot(
    identical(observed_order, sample_order),
    isTRUE(all.equal(observed_bed, target_intervals, check.attributes = FALSE)),
    identical(observed_summary, profile_summary)
  )
  for (name in names(sample_lists)) {
    observed <- readLines(
      file.path(output_directory, paste0(name, ".txt")),
      warn = FALSE
    )
    stopifnot(identical(observed, sample_lists[[name]]))
  }
  invisible(TRUE)
}

if (check_only) {
  assert_profile()
  cat(
    "local9 STAR/rMATS profile check passed:",
    nrow(samples), "libraries,",
    nrow(frozen_events), "events,",
    nrow(target_intervals), "target intervals\n"
  )
  quit(save = "no", status = 0L)
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write_tsv(
  sample_order,
  file.path(output_directory, "local9_sample_order.tsv")
)
write_tsv(
  target_intervals,
  file.path(output_directory, "fixed_events_padded_1kb_merged.bed"),
  column_names = FALSE
)
write_tsv(
  profile_summary,
  file.path(output_directory, "profile_summary.tsv")
)
for (name in names(sample_lists)) {
  writeLines(
    sample_lists[[name]],
    file.path(output_directory, paste0(name, ".txt")),
    useBytes = TRUE
  )
}
assert_profile()

cat(
  "Prepared local9 STAR/rMATS profile:",
  nrow(samples), "libraries,",
  nrow(frozen_events), "events,",
  nrow(target_intervals), "target intervals\n"
)
