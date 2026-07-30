script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_file),
    winslash = "/",
    mustWork = TRUE
  )
  ROOT <- dirname(dirname(script_path))
} else {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

project_library <- file.path(ROOT, ".Rlib")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

output_directories <- c(
  "results/r",
  "results/r/dataset_audit",
  "results/r/differential_expression",
  "results/r/splicing",
  "results/r/integration",
  "results/r/external_validation",
  "results/r/robustness",
  "results/r/cell_resolved",
  "results/r/figures",
  "results/r/raw_confirmation",
  "results/r/publication"
)
invisible(lapply(file.path(ROOT, output_directories), dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

read_count_matrix <- function(accession) {
  path <- file.path(
    ROOT,
    "data",
    "processed",
    paste0(accession, "_counts.csv.gz")
  )
  counts <- read.csv(
    gzfile(path),
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  counts
}

read_sample_metadata <- function(accession) {
  path <- file.path(
    ROOT,
    "data",
    "processed",
    paste0(accession, "_metadata.tsv")
  )
  metadata <- read.delim(
    path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata
}

assert_identical_samples <- function(counts, metadata, label) {
  if (!setequal(colnames(counts), rownames(metadata))) {
    stop(label, ": count and metadata sample IDs differ")
  }
  invisible(TRUE)
}

write_tsv <- function(frame, relative_path, row_names = FALSE) {
  write.table(
    frame,
    file.path(ROOT, relative_path),
    sep = "\t",
    quote = FALSE,
    row.names = row_names,
    col.names = TRUE,
    na = ""
  )
}

zscore <- function(values) {
  as.numeric(scale(as.numeric(values)))
}

normalize_excel_gene_symbol <- function(values) {
  text <- trimws(as.character(values))
  numeric_serial <- grepl("^[0-9]{5}(\\.0+)?$", text)
  if (any(numeric_serial)) {
    serial_dates <- as.Date(
      as.numeric(text[numeric_serial]),
      origin = "1899-12-30"
    )
    text[numeric_serial] <- format(serial_dates, "%Y-%m-%d")
  }
  for (index in seq_along(text)) {
    match <- regexec("^(20[0-9]{2})-([0-9]{2})-([0-9]{2})", text[index])
    pieces <- regmatches(text[index], match)[[1]]
    if (length(pieces) == 4L) {
      number <- as.integer(pieces[2]) - 2000L
      month <- as.integer(pieces[3])
      day <- as.integer(pieces[4])
      if (month == 9L && day == 1L) {
        text[index] <- paste0("SEPTIN", number)
      } else if (month == 3L && day == 1L) {
        text[index] <- paste0("MARCHF", number)
      }
    }

    day_month <- regexec(
      "^([0-9]{1,2})-(Mar|Sep)$",
      text[index],
      ignore.case = TRUE
    )
    day_month_pieces <- regmatches(text[index], day_month)[[1]]
    if (length(day_month_pieces) == 3L) {
      number <- as.integer(day_month_pieces[2])
      family <- tolower(day_month_pieces[3])
      text[index] <- if (family == "mar") {
        paste0("MARCHF", number)
      } else {
        paste0("SEPTIN", number)
      }
    }

    month_day <- regexec(
      "^(Mar|Sep)-([0-9]{1,2})$",
      text[index],
      ignore.case = TRUE
    )
    month_day_pieces <- regmatches(text[index], month_day)[[1]]
    if (length(month_day_pieces) == 3L) {
      family <- tolower(month_day_pieces[2])
      number <- as.integer(month_day_pieces[3])
      text[index] <- if (family == "mar") {
        paste0("MARCHF", number)
      } else {
        paste0("SEPTIN", number)
      }
    }
  }
  text
}

build_hgnc_resolver <- function(hgnc) {
  current <- unique(hgnc$symbol[!is.na(hgnc$symbol)])
  expand_aliases <- function(values) {
    values[is.na(values)] <- ""
    pieces <- strsplit(values, "|", fixed = TRUE)
    data.frame(
      alias = unlist(pieces, use.names = FALSE),
      current_symbol = rep(hgnc$symbol, lengths(pieces)),
      stringsAsFactors = FALSE
    )
  }
  mapping <- unique(rbind(
    expand_aliases(hgnc$prev_symbol),
    expand_aliases(hgnc$alias_symbol)
  ))
  mapping <- mapping[
    nzchar(mapping$alias) & !is.na(mapping$current_symbol),
    ,
    drop = FALSE
  ]
  target_counts <- table(mapping$alias)
  mapping <- mapping[
    as.integer(target_counts[mapping$alias]) == 1L,
    ,
    drop = FALSE
  ]
  list(
    current = current,
    aliases = setNames(mapping$current_symbol, mapping$alias)
  )
}

resolve_symbols <- function(symbols, resolver) {
  symbols <- as.character(symbols)
  resolved <- symbols
  current <- symbols %in% resolver$current
  alias_index <- match(symbols, names(resolver$aliases))
  use_alias <- !current & !is.na(alias_index)
  resolved[use_alias] <- unname(
    resolver$aliases[alias_index[use_alias]]
  )
  resolved
}

remap_unique <- function(frame, symbol_column, resolver, prefix) {
  frame$source_symbol <- frame[[symbol_column]]
  frame$gene_symbol <- resolve_symbols(frame$source_symbol, resolver)
  frame$current_symbol_preferred <- frame$source_symbol == frame$gene_symbol
  ordering <- order(
    frame$gene_symbol,
    -as.integer(frame$current_symbol_preferred)
  )
  frame <- frame[ordering, , drop = FALSE]
  frame <- frame[!duplicated(frame$gene_symbol), , drop = FALSE]
  selected <- setdiff(
    names(frame),
    c(symbol_column, "current_symbol_preferred")
  )
  names(frame)[match(
    setdiff(selected, "gene_symbol"),
    names(frame)
  )] <- paste0(prefix, setdiff(selected, "gene_symbol"))
  frame
}

check_full_rank <- function(design, label) {
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop(
      label,
      " design is not full rank (rank ",
      rank,
      " of ",
      ncol(design),
      ")"
    )
  }
  invisible(TRUE)
}

write_session_info <- function(relative_path) {
  lines <- capture.output(sessionInfo())
  writeLines(lines, file.path(ROOT, relative_path))
}
