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

index_path <- file.path(
  ROOT,
  "manuscript",
  "supplementary_figure_index.tsv"
)
if (!file.exists(index_path)) {
  stop("Supplementary figure index is missing")
}

figure_index <- read.delim(
  index_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_columns <- c(
  "figure_id",
  "filename",
  "title",
  "evidence_role",
  "generating_script"
)
if (!all(required_columns %in% names(figure_index))) {
  stop("Supplementary figure index has missing columns")
}
if (!identical(figure_index$figure_id, paste0("S", seq_len(23L)))) {
  stop("Supplementary figure identifiers must run from S1 through S23")
}
if (anyDuplicated(figure_index$filename)) {
  stop("Supplementary figure filenames must be unique")
}

figure_paths <- file.path(ROOT, "manuscript", figure_index$filename)
missing <- !file.exists(figure_paths)
if (any(missing)) {
  stop(
    "Supplementary figures are missing: ",
    paste(figure_index$filename[missing], collapse = ", ")
  )
}
if (any(file.info(figure_paths)$size <= 0)) {
  stop("One or more supplementary figures are empty")
}

png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
valid_png <- vapply(
  figure_paths,
  function(path) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    identical(readBin(connection, what = "raw", n = 8L), png_signature)
  },
  logical(1)
)
if (!all(valid_png)) {
  stop(
    "Invalid PNG files: ",
    paste(figure_index$filename[!valid_png], collapse = ", ")
  )
}

manuscript <- paste(
  readLines(
    file.path(ROOT, "manuscript", "manuscript.md"),
    warn = FALSE
  ),
  collapse = "\n"
)
unlinked <- !vapply(
  figure_index$filename[seq_len(11L)],
  grepl,
  logical(1),
  x = manuscript,
  fixed = TRUE
)
if (any(unlinked)) {
  stop(
    "Supplementary figures are absent from manuscript legends: ",
    paste(figure_index$filename[seq_len(11L)][unlinked], collapse = ", ")
  )
}
if (!grepl("Supplementary Figures S12-S23", manuscript, fixed = TRUE)) {
  stop("Grouped legend for the frozen 12-event sashimi panel is missing")
}

cat("Supplementary figure package test passed\n")
