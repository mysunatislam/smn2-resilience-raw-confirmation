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

output <- system2(
  rscript,
  file.path(ROOT, "r", "36_audit_manuscript_structure.R"),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(output, "status")
if (is.null(status)) {
  status <- 0L
}
if (status != 0L) {
  stop(paste(output, collapse = "\n"))
}

audit <- read.delim(
  file.path(ROOT, "manuscript", "submission", "manuscript_structure_audit.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (nrow(audit) != 10L || !all(audit$pass)) {
  stop("Manuscript structure audit contains a failed or missing check")
}

cat("Manuscript structure regression test passed\n")
