script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
if (length(script_file) != 1L) {
  stop("Run this audit with Rscript")
}

ROOT <- dirname(dirname(normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)))
manuscript_path <- file.path(ROOT, "manuscript", "manuscript.md")
lines <- readLines(manuscript_path, warn = FALSE)
text <- paste(lines, collapse = "\n")

word_count <- function(x) {
  tokens <- gregexpr("[[:alnum:]][[:alnum:]'-]*", x, perl = TRUE)[[1L]]
  if (identical(tokens, -1L)) 0L else length(tokens)
}

checks <- list()
add_check <- function(id, expected, observed, pass, details = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id,
    expected = as.character(expected),
    observed = as.character(observed),
    pass = isTRUE(pass),
    details = details,
    stringsAsFactors = FALSE
  )
}

abstract_match <- regexec(
  "(?s)## Abstract\\s+(.*?)\\s+## Introduction",
  text,
  perl = TRUE
)
abstract <- regmatches(text, abstract_match)[[1L]][2L]
abstract <- gsub("\\*\\*[^*]+:\\*\\*", "", abstract)
abstract <- gsub("`", "", abstract, fixed = TRUE)
abstract_words <- word_count(abstract)
add_check(
  "abstract_word_count", "150-250", abstract_words,
  abstract_words >= 150L && abstract_words <= 250L
)

keyword_match <- regexec(
  "(?s)\\*\\*Keywords:\\*\\*\\s*(.*?)\\s+## Introduction",
  text,
  perl = TRUE
)
keyword_text <- regmatches(text, keyword_match)[[1L]]
keyword_count <- if (length(keyword_text) == 2L) {
  length(strsplit(keyword_text[2L], ";", fixed = TRUE)[[1L]])
} else {
  NA_integer_
}
add_check(
  "keyword_count", "4-6", keyword_count,
  !is.na(keyword_count) && keyword_count >= 4L && keyword_count <= 6L
)

reference_heading <- grep("^## References$", lines)
if (length(reference_heading) != 1L) {
  stop("Expected exactly one References heading")
}
body_lines <- lines[seq_len(reference_heading - 1L)]
reference_lines <- lines[(reference_heading + 1L):length(lines)]
reference_starts <- grep("^[0-9]+\\. ", reference_lines)
reference_numbers <- as.integer(sub("^([0-9]+)\\..*$", "\\1", reference_lines[reference_starts]))
add_check(
  "reference_number_sequence", "1-19", paste(reference_numbers, collapse = ","),
  identical(reference_numbers, seq_len(19L))
)

reference_blocks <- character(length(reference_starts))
for (i in seq_along(reference_starts)) {
  start <- reference_starts[i]
  end <- if (i < length(reference_starts)) reference_starts[i + 1L] - 1L else length(reference_lines)
  reference_blocks[i] <- paste(reference_lines[start:end], collapse = " ")
}
doi_present <- grepl("https://doi.org/", reference_blocks, fixed = TRUE)
add_check(
  "reference_doi_links", "19/19", paste0(sum(doi_present), "/", length(doi_present)),
  length(doi_present) == 19L && all(doi_present)
)

body_text <- paste(body_lines, collapse = "\n")
body_normalized <- trimws(gsub("\\s+", " ", body_text))
citation_matches <- regmatches(
  body_normalized,
  gregexpr("\\[[0-9]+(?:\\s*[-,]\\s*[0-9]+)*\\]", body_normalized, perl = TRUE)
)[[1L]]
expand_citation <- function(x) {
  x <- gsub("[", "", x, fixed = TRUE)
  x <- gsub("]", "", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)
  pieces <- strsplit(x, ",", fixed = TRUE)[[1L]]
  unlist(lapply(pieces, function(piece) {
    if (grepl("-", piece, fixed = TRUE)) {
      bounds <- as.integer(strsplit(piece, "-", fixed = TRUE)[[1L]])
      seq.int(bounds[1L], bounds[2L])
    } else {
      as.integer(piece)
    }
  }), use.names = FALSE)
}
citations <- unlist(lapply(citation_matches, expand_citation), use.names = FALSE)
first_citation_order <- unique(citations)
add_check(
  "in_text_citation_order", "1-19", paste(first_citation_order, collapse = ","),
  identical(first_citation_order, seq_len(19L))
)
add_check(
  "all_references_cited", "1-19", paste(sort(unique(citations)), collapse = ","),
  identical(sort(unique(citations)), seq_len(19L))
)

main_figure_matches <- regmatches(
  body_normalized,
  gregexpr("Figure [0-9]+", body_normalized, perl = TRUE)
)[[1L]]
main_figure_order <- unique(as.integer(sub("Figure ", "", main_figure_matches, fixed = TRUE)))
add_check(
  "main_figure_citation_order", "1,2,3", paste(main_figure_order, collapse = ","),
  identical(main_figure_order, 1:3)
)

expand_resource_ranges <- function(text, resource) {
  pattern <- paste0("Supplementary ", resource, "s? S[0-9]+(?:-S?[0-9]+)?")
  matches <- regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
  if (length(matches) == 0L) {
    return(integer())
  }
  unlist(lapply(matches, function(match) {
    values <- as.integer(regmatches(match, gregexpr("[0-9]+", match))[[1L]])
    if (length(values) == 1L) values else seq.int(values[1L], values[2L])
  }), use.names = FALSE)
}

supp_figure_order <- unique(expand_resource_ranges(body_normalized, "Figure"))
supp_table_order <- unique(expand_resource_ranges(body_normalized, "Table"))
add_check(
  "supplementary_figure_citation_order", "1-25",
  paste(supp_figure_order, collapse = ","),
  identical(supp_figure_order, seq_len(25L))
)
add_check(
  "supplementary_table_citation_order", "1-15",
  paste(supp_table_order, collapse = ","),
  identical(supp_table_order, seq_len(15L))
)

manuscript_words <- word_count(body_text)
add_check(
  "original_communication_word_limit", "<=9000", manuscript_words,
  manuscript_words <= 9000L,
  "Approximate count excludes the reference list but includes abstract, legends, and declarations."
)

audit <- do.call(rbind, checks)
output_dir <- file.path(ROOT, "manuscript", "submission")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  audit,
  file.path(output_dir, "manuscript_structure_audit.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (!all(audit$pass)) {
  stop(
    "Manuscript structure audit failed: ",
    paste(audit$check_id[!audit$pass], collapse = ", ")
  )
}

marker <- data.frame(
  status = "COMPLETE",
  checked_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  checks = nrow(audit),
  failed = 0L,
  stringsAsFactors = FALSE
)
write.table(
  marker,
  file.path(output_dir, "MANUSCRIPT_STRUCTURE_AUDIT_COMPLETE.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Manuscript structure audit passed (", nrow(audit), " checks)\n", sep = "")
