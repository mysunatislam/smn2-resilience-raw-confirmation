source(file.path("r", "common.R"))

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required; run r/00_install_packages.R")
}

accession <- "GSE108094"
raw_directory <- file.path(ROOT, "data", "raw", accession)
metadata_directory <- file.path(ROOT, "data", "metadata")
dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_directory, recursive = TRUE, showWarnings = FALSE)

resources <- data.frame(
  resource = c("gene_expression", "exon_skipping", "geo_soft"),
  url = c(
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108094/",
      "suppl/GSE108094_gene_exp.diff.gz"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108094/",
      "suppl/GSE108094_SE.MATS.ReadsOnTargetAndJunctionCounts.txt.gz"
    ),
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108094/",
      "soft/GSE108094_family.soft.gz"
    )
  ),
  relative_path = c(
    "data/raw/GSE108094/GSE108094_gene_exp.diff.gz",
    paste0(
      "data/raw/GSE108094/",
      "GSE108094_SE.MATS.ReadsOnTargetAndJunctionCounts.txt.gz"
    ),
    "data/metadata/GSE108094_family.soft.gz"
  ),
  expected_sha256 = c(
    "F2C1872E0FFF1D827F1558BD330AA7DB22D6B30B4B2B672ECE0A76667565EDC6",
    "9F709932198BEF2748885B1A95EF0F1C72F2A5B736E6AF2056023774DC661459",
    "5BAFB5F169DC287FE494AC5B7AA874DDB0E900B914AD989494C0AD40025C9151"
  ),
  stringsAsFactors = FALSE
)

for (index in seq_len(nrow(resources))) {
  path <- file.path(ROOT, resources$relative_path[index])
  if (!file.exists(path)) {
    download.file(
      resources$url[index],
      path,
      mode = "wb",
      quiet = FALSE
    )
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

soft_path <- file.path(metadata_directory, "GSE108094_family.soft.gz")
soft <- readLines(gzfile(soft_path), warn = FALSE)
sample_starts <- grep("^\\^SAMPLE = ", soft)
sample_ends <- c(sample_starts[-1] - 1L, length(soft))
soft_value <- function(lines, prefix, required = TRUE) {
  matches <- lines[startsWith(lines, prefix)]
  if (length(matches) == 0L) {
    if (required) stop("Missing SOFT field: ", prefix)
    return(NA_character_)
  }
  sub(prefix, "", matches[1], fixed = TRUE)
}

sample_records <- vector("list", length(sample_starts))
for (index in seq_along(sample_starts)) {
  block <- soft[sample_starts[index]:sample_ends[index]]
  title <- soft_value(block, "!Sample_title = ")
  group <- if (startsWith(title, "Control")) {
    "Control"
  } else if (startsWith(title, "SMA")) {
    "SMA"
  } else {
    stop("Unexpected GSE108094 title: ", title)
  }
  characteristics <- sub(
    "!Sample_characteristics_ch1 = ",
    "",
    block[startsWith(block, "!Sample_characteristics_ch1 = ")],
    fixed = TRUE
  )
  cell_type <- characteristics[startsWith(characteristics, "cell type: ")]
  if (length(cell_type) != 1L) {
    stop("Could not uniquely parse cell type for ", title)
  }
  line_number <- sub("^.* MN ([0-9]+) Replicate.*$", "\\1", title)
  replicate <- as.integer(sub("^.* Replicate ([0-9]+)$", "\\1", title))
  sample_records[[index]] <- data.frame(
    sample_id = soft_value(block, "!Sample_geo_accession = "),
    title = title,
    group = group,
    line_id = paste0(ifelse(group == "Control", "CTL", "SMA"), "_MN", line_number),
    replicate = replicate,
    organism = soft_value(block, "!Sample_organism_ch1 = "),
    cell_type = sub("^cell type: ", "", cell_type),
    platform = soft_value(block, "!Sample_platform_id = "),
    instrument = soft_value(block, "!Sample_instrument_model = "),
    stringsAsFactors = FALSE
  )
}
metadata <- do.call(rbind, sample_records)
rownames(metadata) <- NULL
stopifnot(
  nrow(metadata) == 8L,
  identical(as.integer(table(metadata$group)), c(4L, 4L)),
  length(unique(metadata$line_id)) == 4L,
  all(table(metadata$line_id) == 2L),
  all(metadata$organism == "Homo sapiens")
)
write_tsv(metadata, "data/processed/GSE108094_metadata.tsv")

hgnc <- read.delim(
  file.path(ROOT, "data", "metadata", "hgnc_complete_set.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)
resolver <- build_hgnc_resolver(hgnc)
approved_hgnc <- hgnc[
  hgnc$status == "Approved" &
    !is.na(hgnc$ensembl_gene_id) & nzchar(hgnc$ensembl_gene_id),
  c("ensembl_gene_id", "symbol")
]
approved_hgnc <- approved_hgnc[!duplicated(approved_hgnc$ensembl_gene_id), ]
ensembl_to_symbol <- setNames(
  approved_hgnc$symbol,
  approved_hgnc$ensembl_gene_id
)

expression_path <- file.path(
  raw_directory,
  "GSE108094_gene_exp.diff.gz"
)
expression <- read.delim(
  gzfile(expression_path),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expression_row_count <- nrow(expression)
expression_ok_row_count <- sum(expression$status == "OK")
expression_fdr05_row_count <- sum(
  expression$status == "OK" & expression$q_value < 0.05,
  na.rm = TRUE
)
required_expression_columns <- c(
  "test_id", "gene_id", "gene", "sample_1", "sample_2", "status",
  "value_1", "value_2", "log2(fold_change)", "p_value", "q_value",
  "significant"
)
stopifnot(
  nrow(expression) == 63657L,
  all(required_expression_columns %in% names(expression)),
  all(expression$sample_1 == "SMA"),
  all(expression$sample_2 == "CTL")
)

orientation_rows <-
  expression$status == "OK" &
  expression$value_1 > 0 & expression$value_2 > 0 &
  is.finite(expression$`log2(fold_change)`)
calculated_ctl_minus_sma <- log2(
  expression$value_2[orientation_rows] /
    expression$value_1[orientation_rows]
)
orientation_difference <- abs(
  calculated_ctl_minus_sma -
    expression$`log2(fold_change)`[orientation_rows]
)
orientation_correlation <- cor(
  calculated_ctl_minus_sma,
  expression$`log2(fold_change)`[orientation_rows],
  method = "spearman",
  use = "complete.obs"
)
stopifnot(
  median(orientation_difference, na.rm = TRUE) < 1e-4,
  orientation_correlation > 0.999
)

expression$source_gene_symbol <- expression$gene
ensembl_match <- match(expression$gene_id, names(ensembl_to_symbol))
expression$gene_symbol <- expression$source_gene_symbol
has_ensembl_mapping <- !is.na(ensembl_match)
expression$gene_symbol[has_ensembl_mapping] <- unname(
  ensembl_to_symbol[ensembl_match[has_ensembl_mapping]]
)
expression$gene_symbol <- resolve_symbols(expression$gene_symbol, resolver)
expression$hgnc_ensembl_mapping <- has_ensembl_mapping
expression$sma_vs_control_log2_effect <- -expression$`log2(fold_change)`
expression$published_cuffdiff_log2_control_vs_sma <-
  expression$`log2(fold_change)`
expression$sma_fpkm <- expression$value_1
expression$control_fpkm <- expression$value_2
expression$official_tested <-
  expression$status == "OK" &
  is.finite(expression$sma_vs_control_log2_effect)
expression$external_direction_usable <-
  expression$official_tested &
  pmax(expression$sma_fpkm, expression$control_fpkm) >= 1

expression <- expression[
  !is.na(expression$gene_symbol) &
    nzchar(expression$gene_symbol) & expression$gene_symbol != "-",
  ,
  drop = FALSE
]
expression <- expression[order(
  expression$gene_symbol,
  -as.integer(expression$official_tested),
  -as.integer(expression$hgnc_ensembl_mapping),
  -pmax(expression$sma_fpkm, expression$control_fpkm),
  expression$q_value,
  expression$test_id
), , drop = FALSE]
expression <- expression[!duplicated(expression$gene_symbol), , drop = FALSE]

external_expression <- expression[, c(
  "gene_symbol", "source_gene_symbol", "gene_id", "locus", "status",
  "sma_fpkm", "control_fpkm", "sma_vs_control_log2_effect",
  "published_cuffdiff_log2_control_vs_sma", "p_value", "q_value",
  "significant", "official_tested", "external_direction_usable",
  "hgnc_ensembl_mapping"
)]
names(external_expression)[names(external_expression) == "gene_id"] <-
  "ensembl_gene_id"
stopifnot(
  !anyDuplicated(external_expression$gene_symbol),
  all(c("TSPOAP1", "SEPTIN11", "OGA") %in%
    external_expression$gene_symbol)
)
write_tsv(
  external_expression,
  "results/r/external_validation/GSE108094_expression_SMA_vs_control.tsv"
)

splicing_path <- file.path(
  raw_directory,
  "GSE108094_SE.MATS.ReadsOnTargetAndJunctionCounts.txt.gz"
)
splicing <- read.delim(
  gzfile(splicing_path),
  skip = 15L,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  all(c(
    "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base",
    "exonEnd", "upstreamES", "upstreamEE", "downstreamES",
    "downstreamEE", "PValue", "FDR", "IncLevelDifference"
  ) %in% names(splicing))
)
splicing$source_gene_symbol <- splicing$geneSymbol
splicing$gene_symbol <- resolve_symbols(splicing$geneSymbol, resolver)
splicing$sma_vs_control_delta_psi <- -splicing$IncLevelDifference
splicing$external_significant <-
  is.finite(splicing$FDR) & splicing$FDR < 0.05 &
  is.finite(splicing$sma_vs_control_delta_psi) &
  abs(splicing$sma_vs_control_delta_psi) > 0.05
external_splicing <- splicing[, c(
  "gene_symbol", "source_gene_symbol", "GeneID", "chr", "strand",
  "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
  "downstreamES", "downstreamEE", "PValue", "FDR",
  "sma_vs_control_delta_psi", "external_significant"
)]
names(external_splicing)[names(external_splicing) == "GeneID"] <-
  "ensembl_gene_id"
write_tsv(
  external_splicing,
  "results/r/external_validation/GSE108094_exon_skipping_hg19.tsv"
)

summary <- data.frame(
  metric = c(
    "libraries", "independent_lines", "control_lines", "sma_lines",
    "official_expression_rows", "official_expression_OK_rows",
    "collapsed_expression_genes", "direction_usable_expression_genes",
    "official_expression_FDR05_rows", "external_ES_rows",
    "external_ES_FDR05_abs_dPSI_gt_005", "orientation_spearman",
    "orientation_median_absolute_error"
  ),
  value = c(
    nrow(metadata), length(unique(metadata$line_id)),
    length(unique(metadata$line_id[metadata$group == "Control"])),
    length(unique(metadata$line_id[metadata$group == "SMA"])),
    expression_row_count, expression_ok_row_count,
    nrow(external_expression),
    sum(external_expression$external_direction_usable),
    expression_fdr05_row_count,
    nrow(external_splicing), sum(external_splicing$external_significant),
    sprintf("%.8f", orientation_correlation),
    sprintf("%.8g", median(orientation_difference, na.rm = TRUE))
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  summary,
  "results/r/external_validation/GSE108094_audit_summary.tsv"
)
write_tsv(
  resources,
  "results/r/external_validation/GSE108094_resource_manifest.tsv"
)
write_session_info(
  "results/r/external_validation/GSE108094_R_sessionInfo.txt"
)

cat(
  "GSE108094 external validation prepared:",
  nrow(metadata), "libraries from", length(unique(metadata$line_id)),
  "lines;", nrow(external_expression), "collapsed genes;",
  sum(external_splicing$external_significant),
  "official hg19 ES events pass FDR/effect filters\n"
)
