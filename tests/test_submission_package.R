script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
stopifnot(length(script_file) == 1L)
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

pdf_path <- file.path(
  ROOT, "output", "pdf", "supplementary_material_review.pdf"
)
zip_path <- file.path(
  ROOT, "output", "submission", "supplementary_tables_S1-S15.zip"
)
manifest_path <- file.path(
  ROOT, "output", "submission", "submission_artifact_manifest.tsv"
)
stopifnot(
  file.exists(pdf_path),
  file.info(pdf_path)$size > 1000000,
  file.exists(zip_path),
  file.info(zip_path)$size > 1000,
  file.exists(manifest_path)
)

connection <- file(pdf_path, open = "rb")
signature <- rawToChar(readBin(connection, what = "raw", n = 5L))
close(connection)
stopifnot(identical(signature, "%PDF-"))

archive <- unzip(zip_path, list = TRUE)
expected_tables <- paste0(
  "supplementary_table_S",
  1:15,
  c(
    "_candidates.tsv",
    "_summary.tsv",
    "_concordance.tsv",
    "_raw_splice_confirmation.tsv",
    "_primary_raw_splice_panel.tsv",
    "_raw_splice_unconfirmed_reasons.tsv",
    "_raw_junction_support.tsv",
    "_raw_splice_summary.tsv",
    "_raw_splice_provenance.tsv",
    "_permutation_summary.tsv",
    "_negative_controls.tsv",
    "_score_sensitivity.tsv",
    "_candidate_rank_stability.tsv",
    "_bootstrap_candidate_stability.tsv",
    "_bootstrap_summary.tsv"
  )
)
stopifnot(
  nrow(archive) == 17L,
  setequal(
    archive$Name,
    c(expected_tables, "supplementary_table_index.tsv", "README.txt")
  )
)

manifest <- read.delim(
  manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(manifest) == 2L,
  all(manifest$submission_status == "REVIEW_ONLY"),
  all(nchar(manifest$sha256) == 64L),
  identical(
    as.numeric(manifest$bytes),
    c(file.info(pdf_path)$size, file.info(zip_path)$size)
  )
)

readiness <- read.delim(
  file.path(ROOT, "manuscript", "submission", "submission_readiness.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  readiness$status[readiness$item == "Supplementary PDF"] ==
    "PASS_REVIEW_BUILD",
  readiness$status[readiness$item == "Overall submission status"] ==
    "NOT_READY"
)
cat("Submission package test passed\n")
