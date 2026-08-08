script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
stopifnot(length(script_file) == 1L)
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

audit_path <- file.path(
  ROOT, "manuscript", "submission", "figure_asset_audit.tsv"
)
marker_path <- file.path(
  ROOT, "manuscript", "submission", "FIGURE_ASSET_AUDIT_COMPLETE.tsv"
)
stopifnot(file.exists(audit_path), file.exists(marker_path))
audit <- read.delim(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
marker <- read.delim(marker_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(
  nrow(audit) == 28L,
  all(audit$status == "PASS"),
  identical(marker$status, "COMPLETE"),
  marker$assets_checked == 28L,
  marker$assets_failed == 0L,
  sum(audit$publication_basis == "png_600_dpi") == 11L,
  sum(audit$publication_basis == "vector_pdf") == 17L
)
cat("Figure asset audit test passed\n")
