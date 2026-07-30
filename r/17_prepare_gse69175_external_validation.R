source(file.path("r", "common.R"))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required; run r/00_install_packages.R")
}

accession <- "GSE69175"
raw_directory <- file.path(ROOT, "data", "raw", accession)
metadata_directory <- file.path(ROOT, "data", "metadata")
dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_directory, recursive = TRUE, showWarnings = FALSE)

resources <- data.frame(
  resource = c("gene_expression", "geo_soft"),
  url = c(
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69175/",
      "suppl/GSE69175_gene_exp.txt.gz"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69175/",
      "soft/GSE69175_family.soft.gz"
    )
  ),
  relative_path = c(
    "data/raw/GSE69175/GSE69175_gene_exp.txt.gz",
    "data/metadata/GSE69175_family.soft.gz"
  ),
  expected_sha256 = c(
    "3B8AE04FB61F52EFECACEEFA4A83EDC489BF64C9E6B2CCEACF6BDF44343ECDBC",
    "FDE8CA9192B563FDE51EEC0FFAF6324C783F46862A521B7957DD0C4945C109BD"
  ),
  stringsAsFactors = FALSE
)

for (index in seq_len(nrow(resources))) {
  path <- file.path(ROOT, resources$relative_path[index])
  if (!file.exists(path)) {
    download.file(resources$url[index], path, mode = "wb", quiet = FALSE)
  }
  observed <- toupper(digest::digest(path, algo = "sha256", file = TRUE))
  if (!identical(observed, resources$expected_sha256[index])) {
    stop(
      "SHA-256 mismatch for ", resources$resource[index],
      ": observed ", observed
    )
  }
  resources$observed_sha256[index] <- observed
  resources$bytes[index] <- file.info(path)$size
}
resources$sha256_pass <-
  resources$observed_sha256 == resources$expected_sha256

soft <- readLines(
  gzfile(file.path(metadata_directory, "GSE69175_family.soft.gz")),
  warn = FALSE
)
sample_starts <- grep("^\\^SAMPLE = ", soft)
sample_ends <- c(sample_starts[-1] - 1L, length(soft))
soft_value <- function(lines, prefix) {
  matches <- lines[startsWith(lines, prefix)]
  if (length(matches) == 0L) stop("Missing SOFT field: ", prefix)
  sub(prefix, "", matches[1], fixed = TRUE)
}

sample_records <- vector("list", length(sample_starts))
for (index in seq_along(sample_starts)) {
  block <- soft[sample_starts[index]:sample_ends[index]]
  title <- soft_value(block, "!Sample_title = ")
  characteristics <- sub(
    "!Sample_characteristics_ch1 = ",
    "",
    block[startsWith(block, "!Sample_characteristics_ch1 = ")],
    fixed = TRUE
  )
  line_value <- characteristics[startsWith(characteristics, "cell line: ")]
  cell_value <- characteristics[startsWith(
    characteristics,
    "derived cell typed: "
  )]
  stopifnot(length(line_value) == 1L, length(cell_value) == 1L)
  group <- if (startsWith(title, "control_")) {
    "Control"
  } else if (startsWith(title, "smn_")) {
    "SMA"
  } else {
    stop("Unexpected GSE69175 sample title: ", title)
  }
  sample_records[[index]] <- data.frame(
    sample_id = soft_value(block, "!Sample_geo_accession = "),
    title = title,
    group = group,
    line_id = sub("^cell line: ", "", line_value),
    replicate = as.integer(sub("^.*rep([0-9]+)_RNA-seq$", "\\1", title)),
    source = soft_value(block, "!Sample_source_name_ch1 = "),
    organism = soft_value(block, "!Sample_organism_ch1 = "),
    cell_type = sub("^derived cell typed: ", "", cell_value),
    platform = soft_value(block, "!Sample_platform_id = "),
    instrument = soft_value(block, "!Sample_instrument_model = "),
    stringsAsFactors = FALSE
  )
}
metadata <- do.call(rbind, sample_records)
rownames(metadata) <- NULL
line_groups <- unique(metadata[, c("line_id", "group")])
stopifnot(
  nrow(metadata) == 4L,
  identical(as.integer(table(metadata$group)), c(2L, 2L)),
  length(unique(metadata$line_id)) == 2L,
  nrow(line_groups) == 2L,
  all(table(metadata$line_id) == 2L),
  all(metadata$organism == "Homo sapiens"),
  all(metadata$cell_type == "motor neurons")
)
write_tsv(metadata, "data/processed/GSE69175_metadata.tsv")

expression <- read.delim(
  gzfile(file.path(raw_directory, "GSE69175_gene_exp.txt.gz")),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_columns <- c(
  "test_id", "gene_id", "gene", "locus", "sample_1", "sample_2",
  "status", "value_1", "value_2", "log2(fold_change)", "p_value",
  "q_value", "significant"
)
stopifnot(
  nrow(expression) == 19286L,
  all(required_columns %in% names(expression)),
  all(expression$sample_1 == "ctrl"),
  all(expression$sample_2 == "smn")
)
expression_row_count <- nrow(expression)
expression_ok_row_count <- sum(expression$status == "OK")
expression_fdr05_row_count <- sum(
  expression$status == "OK" & expression$q_value < 0.05,
  na.rm = TRUE
)

expression$published_log2_sma_vs_control <- suppressWarnings(as.numeric(
  expression$`log2(fold_change)`
))
orientation_rows <-
  expression$status == "OK" & expression$value_1 > 0 &
  expression$value_2 > 0 &
  is.finite(expression$published_log2_sma_vs_control)
calculated_sma_vs_control <- log2(
  expression$value_2[orientation_rows] /
    expression$value_1[orientation_rows]
)
orientation_error <- abs(
  calculated_sma_vs_control -
    expression$published_log2_sma_vs_control[orientation_rows]
)
orientation_correlation <- cor(
  calculated_sma_vs_control,
  expression$published_log2_sma_vs_control[orientation_rows],
  method = "spearman"
)
stopifnot(
  median(orientation_error) < 1e-4,
  orientation_correlation > 0.999
)

hgnc <- read.delim(
  file.path(ROOT, "data", "metadata", "hgnc_complete_set.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)
resolver <- build_hgnc_resolver(hgnc)
expression$source_gene_symbol <- expression$gene
expression$normalized_source_gene_symbol <- normalize_excel_gene_symbol(
  expression$source_gene_symbol
)
expression$gene_symbol <- resolve_symbols(
  expression$normalized_source_gene_symbol,
  resolver
)
approved <- unique(hgnc$symbol[hgnc$status == "Approved"])
expression$approved_hgnc_gene <- expression$gene_symbol %in% approved
expression$sma_vs_control_log2_effect <-
  expression$published_log2_sma_vs_control
expression$control_fpkm <- expression$value_1
expression$sma_fpkm <- expression$value_2
expression$official_tested <-
  expression$status == "OK" &
  is.finite(expression$sma_vs_control_log2_effect)
expression$external_direction_usable <-
  expression$approved_hgnc_gene & expression$official_tested &
  pmax(expression$control_fpkm, expression$sma_fpkm) >= 1

expression <- expression[
  expression$approved_hgnc_gene & nzchar(expression$gene_symbol),
  ,
  drop = FALSE
]
expression <- expression[order(
  expression$gene_symbol,
  -as.integer(expression$external_direction_usable),
  -pmax(expression$control_fpkm, expression$sma_fpkm),
  expression$q_value,
  expression$test_id
), , drop = FALSE]
expression <- expression[!duplicated(expression$gene_symbol), , drop = FALSE]

external <- expression[, c(
  "gene_symbol", "source_gene_symbol", "normalized_source_gene_symbol",
  "gene_id", "locus", "status", "control_fpkm", "sma_fpkm",
  "sma_vs_control_log2_effect", "p_value", "q_value", "significant",
  "official_tested", "external_direction_usable", "approved_hgnc_gene"
)]
names(external)[names(external) == "gene_id"] <- "source_gene_id"
stopifnot(
  !anyDuplicated(external$gene_symbol),
  all(c("HSPA5", "DDIT3", "XBP1", "ATF4", "CASP4") %in%
    external$gene_symbol)
)
write_tsv(
  external,
  "results/r/external_validation/GSE69175_expression_SMA_vs_control.tsv"
)

resampling <- data.frame(
  accession = accession,
  contrast = "SMA_minus_control",
  requested_scheme = "LOLO",
  biological_units = 2L,
  control_units = 1L,
  sma_units = 1L,
  estimable = FALSE,
  reason = paste(
    "One control line and one SMA line; holding out either line removes",
    "one genotype. The two libraries per line are not independent lines."
  ),
  stringsAsFactors = FALSE
)
summary <- data.frame(
  metric = c(
    "libraries", "independent_lines", "control_lines", "sma_lines",
    "official_expression_rows", "official_expression_OK_rows",
    "collapsed_approved_genes", "direction_usable_genes",
    "official_expression_FDR05_rows", "orientation_spearman",
    "orientation_median_absolute_error", "LOLO_estimable",
    "validation_role"
  ),
  value = c(
    nrow(metadata), length(unique(metadata$line_id)),
    sum(line_groups$group == "Control"), sum(line_groups$group == "SMA"),
    expression_row_count, expression_ok_row_count, nrow(external),
    sum(external$external_direction_usable),
    expression_fdr05_row_count,
    sprintf("%.8f", orientation_correlation),
    sprintf("%.8g", median(orientation_error)), FALSE,
    "external_directional_sensitivity_only"
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  summary,
  "results/r/external_validation/GSE69175_audit_summary.tsv"
)
write_tsv(
  resampling,
  "results/r/external_validation/GSE69175_resampling_feasibility.tsv"
)
write_tsv(
  resources,
  "results/r/external_validation/GSE69175_resource_manifest.tsv"
)
write_session_info(
  "results/r/external_validation/GSE69175_R_sessionInfo.txt"
)

cat(
  "GSE69175 external sensitivity validation prepared:",
  nrow(metadata), "libraries from", length(unique(metadata$line_id)),
  "lines;", sum(external$external_direction_usable),
  "direction-usable approved genes; LOLO is structurally not estimable\n"
)
