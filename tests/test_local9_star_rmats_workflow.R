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
fixed_event_directory <- file.path(
  ROOT,
  "config",
  "rmats",
  "GSE290979",
  "fixed_events"
)
fixed_event_files <- file.path(
  fixed_event_directory,
  paste0("fromGTF.", c("SE", "A5SS", "A3SS", "MXE", "RI"), ".txt")
)
fixed_events <- lapply(
  fixed_event_files,
  read.delim,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
fixed_gene_ids <- unlist(lapply(fixed_events, `[[`, "GeneID"), use.names = FALSE)
gene_id_manifest <- read.delim(
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

stopifnot(
  nrow(order_table) == 9L,
  nrow(target_bed) == 73L,
  profile_value[["fixed_events"]] == "83",
  profile_value[["primary_events"]] == "12",
  profile_value[["reference"]] == "GENCODE_v47_GRCh38_primary_assembly",
  profile_value[["star_version"]] == "2.7.10a",
  profile_value[["star_twopass_mode"]] == "None",
  profile_value[["star_read_map_number"]] == "-1",
  profile_value[["rmats_version"]] == "4.3.0",
  profile_value[["rmats2sashimiplot_version"]] == "4.0.0",
  profile_value[["fastq_deletion"]] == "FALSE",
  profile_value[["bam_deletion"]] == "FALSE",
  sum(vapply(fixed_events, nrow, integer(1))) == 83L,
  nrow(gene_id_manifest) == 83L,
  sum(grepl("^ENSG[0-9]+\\.[0-9]+$", fixed_gene_ids)) +
    sum(gene_id_manifest$match_method == "unmapped_in_gencode_v47") == 83L
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
status <- system2(
  rscript,
  c(
    file.path(ROOT, "r", "31_reconcile_fixed_event_gene_ids.R"),
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
windows_launcher_path <- file.path(
  ROOT,
  "local",
  "gse290979_star_rmats",
  "run_remaining_cohort.cmd"
)
configuration_path <- file.path(
  ROOT,
  "local",
  "gse290979_star_rmats",
  "config.env.example"
)
workflow <- paste(readLines(workflow_path, warn = FALSE), collapse = "\n")
bootstrap <- paste(readLines(bootstrap_path, warn = FALSE), collapse = "\n")
windows_launcher <- paste(
  readLines(windows_launcher_path, warn = FALSE),
  collapse = "\n"
)
configuration <- readLines(configuration_path, warn = FALSE)
deviation <- file.path(
  ROOT,
  "docs",
  "raw_splice_confirmation_execution_deviation_2026-08-01.md"
)

stopifnot(
  !grepl("(^|\\n)[[:space:]]*rm[[:space:]]", workflow),
  !grepl("(^|\\n)[[:space:]]*rm[[:space:]]", bootstrap),
  !grepl("(^|\\n)[[:space:]]*(del|erase)[[:space:]]", windows_launcher),
  grepl('DELETE_FASTQ:-0.*!= "0"', workflow),
  grepl('DELETE_BAM:-0.*!= "0"', workflow),
  grepl("--novelSS", workflow, fixed = TRUE),
  !grepl("--fixed-event-set", workflow, fixed = TRUE),
  grepl("target_locus_discovery", workflow, fixed = TRUE),
  grepl("rmats_result_row_count", workflow, fixed = TRUE),
  grepl("--outStd SAM", workflow, fixed = TRUE),
  grepl('--twopassMode "${STAR_TWOPASS_MODE}"', workflow, fixed = TRUE),
  grepl('--readMapNumber "${STAR_READ_MAP_NUMBER}"', workflow, fixed = TRUE),
  grepl('--outTmpDir "${star_linux_tmp}"', workflow, fixed = TRUE),
  grepl('samtools view -@ 1 -u -L "${TARGET_BED}"', workflow, fixed = TRUE),
  grepl("phase_sashimi", workflow, fixed = TRUE),
  grepl("frozen_plot_only_definition", workflow, fixed = TRUE),
  grepl(
    "run_local9_star_rmats.sh align",
    windows_launcher,
    fixed = TRUE
  ),
  any(configuration == "EXPECTED_STAR_VERSION=2.7.10a"),
  any(configuration == "STAR_TWOPASS_MODE=None"),
  any(configuration == "STAR_READ_MAP_NUMBER=-1"),
  any(configuration == "EXPECTED_RMATS_VERSION=4.3.0"),
  any(configuration == "EXPECTED_RMATS2SASHIMIPLOT_VERSION=4.0.0"),
  any(configuration ==
    "STAR_LINUX_TMP_ROOT=/tmp/smn2_gse290979_star_rmats"),
  any(configuration == "DELETE_FASTQ=0"),
  any(configuration == "DELETE_BAM=0"),
  file.exists(deviation)
)

cat("local9 STAR/rMATS workflow safety test passed\n")
