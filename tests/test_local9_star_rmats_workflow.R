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

profile_directory <- file.path(
  ROOT,
  "config",
  "rmats",
  "GSE290979",
  "local9"
)
profile <- read.delim(
  file.path(profile_directory, "profile_summary.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
profile_value <- setNames(profile$value, profile$metric)
order_table <- read.delim(
  file.path(profile_directory, "local9_sample_order.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
target_bed <- read.delim(
  file.path(profile_directory, "fixed_events_padded_1kb_merged.bed"),
  header = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(order_table) == 9L,
  nrow(target_bed) == 73L,
  profile_value[["fixed_events"]] == "83",
  profile_value[["primary_events"]] == "12",
  profile_value[["reference"]] == "GENCODE_v47_GRCh38_primary_assembly",
  profile_value[["star_version"]] == "2.7.10a",
  profile_value[["rmats_version"]] == "4.3.0",
  profile_value[["rmats2sashimiplot_version"]] == "4.0.0",
  profile_value[["fastq_deletion"]] == "FALSE",
  profile_value[["bam_deletion"]] == "FALSE"
)

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
status <- system2(
  rscript,
  c(
    file.path(ROOT, "r", "29_prepare_local9_star_rmats_profile.R"),
    "--check-only"
  )
)
stopifnot(status == 0L)

workflow_path <- file.path(
  ROOT,
  "local",
  "gse290979_star_rmats",
  "run_local9_star_rmats.sh"
)
bootstrap_path <- file.path(
  ROOT,
  "local",
  "gse290979_star_rmats",
  "bootstrap_wsl_tools.sh"
)
configuration_path <- file.path(
  ROOT,
  "local",
  "gse290979_star_rmats",
  "config.env.example"
)
workflow <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")
bootstrap <- paste(readLines(bootstrap_path, warn = FALSE), collapse = "\n")
configuration <- readLines(configuration_path, warn = FALSE)
deviation <- file.path(
  ROOT,
  "docs",
  "raw_splice_confirmation_execution_deviation_2026-08-01.md"
)

stopifnot(
  !grepl("(^|\\n)[[:space:]]*rm[[:space:]]", workflow),
  !grepl("(^|\\n)[[:space:]]*rm[[:space:]]", bootstrap),
  grepl('DELETE_FASTQ:-0.*!= "0"', workflow),
  grepl('DELETE_BAM:-0.*!= "0"', workflow),
  grepl("--fixed-event-set", workflow, fixed = TRUE),
  grepl("--outStd SAM", workflow, fixed = TRUE),
  grepl('samtools view -@ 1 -u -L "${TARGET_BED}"', workflow, fixed = TRUE),
  grepl("phase_sashimi", workflow, fixed = TRUE),
  any(configuration == "EXPECTED_STAR_VERSION=2.7.10a"),
  any(configuration == "EXPECTED_RMATS_VERSION=4.3.0"),
  any(configuration == "EXPECTED_RMATS2SASHIMIPLOT_VERSION=4.0.0"),
  any(configuration == "DELETE_FASTQ=0"),
  any(configuration == "DELETE_BAM=0"),
  file.exists(deviation)
)

cat("local9 STAR/rMATS workflow safety test passed\n")
