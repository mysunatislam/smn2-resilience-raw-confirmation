script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
if (length(script_file) != 1L) {
  stop("This script must be run with Rscript")
}
script_path <- normalizePath(
  sub("^--file=", "", script_file),
  winslash = "/",
  mustWork = TRUE
)
ROOT <- dirname(dirname(script_path))

arguments <- commandArgs(trailingOnly = TRUE)
argument_value <- function(prefix, default = NULL) {
  matches <- arguments[startsWith(arguments, prefix)]
  if (length(matches) == 0L) {
    return(default)
  }
  if (length(matches) != 1L) {
    stop("Expected one argument beginning with ", prefix)
  }
  sub(prefix, "", matches, fixed = TRUE)
}

gtf_path <- argument_value("--gtf=")
fixed_event_directory <- argument_value(
  "--fixed-event-dir=",
  file.path(ROOT, "config", "rmats", "GSE290979", "fixed_events")
)
manifest_path <- argument_value(
  "--manifest=",
  file.path(
    ROOT,
    "config",
    "rmats",
    "GSE290979",
    "fixed_event_gene_id_mapping.tsv"
  )
)
check_only <- "--check-only" %in% arguments

event_specs <- list(
  ES = c(
    "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ),
  A5SS = c(
    "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
    "flankingES", "flankingEE"
  ),
  A3SS = c(
    "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
    "flankingES", "flankingEE"
  ),
  MXE = c(
    "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
    "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES",
    "downstreamEE"
  ),
  RI = c(
    "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  )
)
rmats_types <- c(ES = "SE", A5SS = "A5SS", A3SS = "A3SS", MXE = "MXE", RI = "RI")

read_fixed_events <- function() {
  frames <- lapply(names(event_specs), function(event_type) {
    path <- file.path(
      fixed_event_directory,
      paste0("fromGTF.", rmats_types[[event_type]], ".txt")
    )
    frame <- read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = ""
    )
    required <- c(
      "ID", "GeneID", "geneSymbol", "chr", "strand",
      event_specs[[event_type]]
    )
    missing <- setdiff(required, names(frame))
    if (length(missing) > 0L) {
      stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "))
    }
    coordinates <- as.matrix(frame[, event_specs[[event_type]], drop = FALSE])
    storage.mode(coordinates) <- "numeric"
    frame$event_type <- event_type
    frame$event_start_0base <- apply(coordinates, 1L, min)
    frame$event_end_1base <- apply(coordinates, 1L, max)
    frame$source_path <- path
    frame
  })
  names(frames) <- names(event_specs)
  frames
}

validate_outputs <- function(frames) {
  total <- sum(vapply(frames, nrow, integer(1)))
  if (total != 83L) {
    stop("Expected 83 fixed events, found ", total)
  }
  if (!file.exists(manifest_path)) {
    stop("Missing mapping manifest: ", manifest_path)
  }
  manifest <- read.delim(
    manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (
    nrow(manifest) != 83L ||
      anyDuplicated(paste(manifest$event_type, manifest$event_id, sep = "|"))
  ) {
    stop("Mapping manifest must contain 83 unique fixed events")
  }
  observed <- do.call(rbind, lapply(names(frames), function(event_type) {
    data.frame(
      event_type = event_type,
      event_id = frames[[event_type]]$ID,
      fixed_gene_id = frames[[event_type]]$GeneID,
      stringsAsFactors = FALSE
    )
  }))
  manifest_key <- paste(manifest$event_type, manifest$event_id, sep = "|")
  observed_key <- paste(observed$event_type, observed$event_id, sep = "|")
  matched_manifest <- manifest[match(observed_key, manifest_key), , drop = FALSE]
  valid_mapped <- grepl("^ENSG[0-9]+\\.[0-9]+$", observed$fixed_gene_id) &
    observed$fixed_gene_id == matched_manifest$gencode_gene_id
  valid_unmapped <- matched_manifest$match_method == "unmapped_in_gencode_v47" &
    observed$fixed_gene_id == matched_manifest$previous_gene_id &
    matched_manifest$gencode_gene_id == ""
  if (anyNA(valid_mapped) || anyNA(valid_unmapped) || !all(valid_mapped | valid_unmapped)) {
    stop("Fixed-event GeneID values do not agree with the mapping manifest")
  }
  invisible(TRUE)
}

frames <- read_fixed_events()
if (check_only) {
  validate_outputs(frames)
  cat("fixed-event GENCODE GeneID check passed\n")
  quit(save = "no", status = 0L)
}

if (is.null(gtf_path) || !file.exists(gtf_path)) {
  stop("Provide the locked GENCODE GTF with --gtf=")
}

extract_attribute <- function(attributes, key) {
  pattern <- paste0("(^|;[[:space:]]*)", key, "[[:space:]]+\"([^\"]+)\"")
  matches <- regexec(pattern, attributes, perl = TRUE)
  values <- regmatches(attributes, matches)
  vapply(
    values,
    function(value) if (length(value) >= 3L) value[[3L]] else NA_character_,
    character(1)
  )
}

read_gtf_genes <- function(path, chunk_size = 100000L) {
  connection <- file(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  chunks <- vector("list", 0L)
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }
    fields <- strsplit(
      lines[grepl("\tgene\t", lines, fixed = TRUE)],
      "\t",
      fixed = TRUE
    )
    fields <- fields[vapply(fields, length, integer(1)) == 9L]
    if (length(fields) == 0L) {
      next
    }
    matrix <- do.call(rbind, fields)
    attributes <- matrix[, 9L]
    chunks[[length(chunks) + 1L]] <- data.frame(
      chromosome = matrix[, 1L],
      start_1base = as.numeric(matrix[, 4L]),
      end_1base = as.numeric(matrix[, 5L]),
      strand = matrix[, 7L],
      gene_id = extract_attribute(attributes, "gene_id"),
      gene_name = extract_attribute(attributes, "gene_name"),
      stringsAsFactors = FALSE
    )
  }
  genes <- do.call(rbind, chunks)
  rownames(genes) <- NULL
  genes <- genes[
    !is.na(genes$gene_id) &
      !is.na(genes$gene_name) &
      grepl("^ENSG[0-9]+\\.[0-9]+$", genes$gene_id),
    ,
    drop = FALSE
  ]
  genes
}

genes <- read_gtf_genes(gtf_path)
if (nrow(genes) < 70000L) {
  stop("Unexpectedly few GENCODE gene records: ", nrow(genes))
}

mapping_records <- vector("list", 0L)
for (event_type in names(frames)) {
  frame <- frames[[event_type]]
  for (index in seq_len(nrow(frame))) {
    event <- frame[index, , drop = FALSE]
    contained <- genes[
      genes$chromosome == event$chr &
        genes$strand == event$strand &
        genes$start_1base <= event$event_start_0base + 1 &
        genes$end_1base >= event$event_end_1base,
      ,
      drop = FALSE
    ]
    overlapping <- genes[
      genes$chromosome == event$chr &
        genes$strand == event$strand &
        genes$start_1base <= event$event_end_1base &
        genes$end_1base >= event$event_start_0base + 1,
      ,
      drop = FALSE
    ]
    exact_symbol <- contained[contained$gene_name == event$geneSymbol, , drop = FALSE]
    exact_symbol_overlap <- overlapping[
      overlapping$gene_name == event$geneSymbol,
      ,
      drop = FALSE
    ]
    if (nrow(exact_symbol) == 1L) {
      match <- exact_symbol
      match_method <- "symbol_chromosome_strand_containment"
    } else if (nrow(exact_symbol) > 1L) {
      stop("Multiple exact GENCODE genes for ", event_type, " event ", event$ID)
    } else if (nrow(exact_symbol_overlap) == 1L) {
      match <- exact_symbol_overlap
      match_method <- "symbol_chromosome_strand_overlap"
    } else if (nrow(exact_symbol_overlap) > 1L) {
      stop("Multiple overlapping GENCODE genes for ", event_type, " event ", event$ID)
    } else if (nrow(contained) == 1L) {
      match <- contained
      match_method <- "chromosome_strand_containment_alias"
    } else if (nrow(overlapping) == 1L) {
      match <- overlapping
      match_method <- "chromosome_strand_overlap_alias"
    } else {
      match <- NULL
      match_method <- "unmapped_in_gencode_v47"
    }
    previous_gene_id <- frame$GeneID[[index]]
    if (!is.null(match)) {
      frame$GeneID[[index]] <- match$gene_id[[1L]]
    }
    mapping_records[[length(mapping_records) + 1L]] <- data.frame(
      event_type = event_type,
      event_id = event$ID,
      frozen_gene_symbol = event$geneSymbol,
      previous_gene_id = previous_gene_id,
      gencode_gene_id = if (is.null(match)) "" else match$gene_id[[1L]],
      gencode_gene_name = if (is.null(match)) "" else match$gene_name[[1L]],
      chromosome = event$chr,
      strand = event$strand,
      event_start_0base = event$event_start_0base,
      event_end_1base = event$event_end_1base,
      match_method = match_method,
      stringsAsFactors = FALSE
    )
  }
  output_columns <- setdiff(
    names(frame),
    c("event_type", "event_start_0base", "event_end_1base", "source_path")
  )
  write.table(
    frame[, output_columns, drop = FALSE],
    unique(frame$source_path),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  frames[[event_type]] <- frame
}

manifest <- do.call(rbind, mapping_records)
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
write.table(
  manifest,
  manifest_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

validate_outputs(read_fixed_events())
cat(
  "Mapped 83 frozen events to GENCODE GeneID values; alias matches: ",
  sum(grepl("alias$", manifest$match_method)),
  "; unmapped: ",
  sum(manifest$match_method == "unmapped_in_gencode_v47"),
  "\n",
  sep = ""
)
