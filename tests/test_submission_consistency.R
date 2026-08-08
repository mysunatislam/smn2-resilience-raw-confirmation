script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
if (length(script_file) != 1L) {
  stop("Run this test with Rscript")
}

ROOT <- dirname(dirname(normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)))
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
}

status <- system2(
  rscript,
  file.path(ROOT, "r", "35_audit_submission_consistency.R"),
  stdout = TRUE,
  stderr = TRUE
)
exit_status <- attr(status, "status")
if (is.null(exit_status)) {
  exit_status <- 0L
}
if (exit_status != 0L) {
  stop(paste(status, collapse = "\n"))
}

audit_path <- file.path(
  ROOT, "manuscript", "submission", "numerical_consistency_audit.tsv"
)
marker_path <- file.path(
  ROOT, "manuscript", "submission", "NUMERICAL_CONSISTENCY_AUDIT_COMPLETE.tsv"
)
if (!file.exists(audit_path) || !file.exists(marker_path)) {
  stop("Submission consistency outputs were not created")
}

audit <- read.delim(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
marker <- read.delim(marker_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(audit) < 15L || !all(audit$pass)) {
  stop("Submission consistency audit contains a failed or missing check")
}
if (nrow(marker) != 1L || marker$status != "COMPLETE" || marker$failed != 0L) {
  stop("Submission consistency completion marker is invalid")
}

cat("Submission consistency regression test passed\n")
