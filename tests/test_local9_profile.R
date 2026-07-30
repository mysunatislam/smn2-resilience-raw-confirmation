source(file.path("r", "common.R"))

profile_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_local9_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selection <- read.delim(
  file.path(ROOT, "config", "GSE290979_local9_selection.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
profile_directory <- file.path(
  ROOT,
  "config",
  "local_raw",
  "GSE290979",
  "local9"
)
summary <- read.delim(
  file.path(profile_directory, "profile_summary.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
metric <- setNames(summary$value, summary$metric)
benchmark_sample <- readLines(
  file.path(profile_directory, "benchmark_sample.txt"),
  warn = FALSE
)

stopifnot(
  nrow(profile_sheet) == 9L,
  identical(profile_sheet$array_index, seq_len(9L)),
  !anyDuplicated(profile_sheet$sample_id),
  identical(profile_sheet$sample_id, selection$sample_id),
  sum(profile_sheet$treatment == "NT") == 5L,
  sum(profile_sheet$treatment == "Scramble") == 2L,
  sum(profile_sheet$treatment == "R6-Mo") == 2L,
  length(unique(
    profile_sheet$donor_line[profile_sheet$treatment == "NT"]
  )) == 5L,
  length(unique(
    profile_sheet$donor_line[profile_sheet$treatment != "NT"]
  )) == 2L,
  length(benchmark_sample) == 1L,
  benchmark_sample %in% profile_sheet$sample_id,
  abs(sum(profile_sheet$fastq_bytes_total) / 1024^3 - 82.193) < 0.001,
  benchmark_sample == profile_sheet$sample_id[
    order(
      profile_sheet$fastq_bytes_total,
      profile_sheet$run_accession
    )[1L]
  ],
  as.numeric(metric[["selected_fastq_bytes"]]) ==
    sum(profile_sheet$fastq_bytes_total),
  metric[["raw_scope"]] ==
    paste0(
      "pre_specified_local_raw_read_sensitivity_not_validation21_",
      "and_not_full_reanalysis"
    )
)

full_sheet <- read.delim(
  file.path(ROOT, "config", "GSE290979_raw_sample_sheet.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
selection_groups <- split(
  profile_sheet,
  paste(profile_sheet$donor_line, profile_sheet$treatment, sep = "__")
)
for (group_name in names(selection_groups)) {
  selected <- selection_groups[[group_name]]
  candidates <- full_sheet[
    full_sheet$donor_line == selected$donor_line &
      full_sheet$treatment == selected$treatment,
    ,
    drop = FALSE
  ]
  distance <- abs(
    candidates$read_pairs - stats::median(candidates$read_pairs)
  )
  expected <- candidates[
    order(distance, candidates$run_accession)[1L],
    "sample_id"
  ]
  stopifnot(identical(selected$sample_id, expected))
}

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
status <- system2(
  rscript,
  c(
    file.path(ROOT, "r", "22_prepare_gse290979_local9_profile.R"),
    "--check-only"
  )
)
stopifnot(status == 0L)

cat("local9 deterministic selection tests passed\n")
