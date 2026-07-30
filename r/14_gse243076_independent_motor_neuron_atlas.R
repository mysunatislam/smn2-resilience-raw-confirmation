source(file.path("r", "common.R"))

for (package in c("readxl", "digest", "Matrix")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package is unavailable: ", package)
  }
}

raw_root <- file.path(ROOT, "data", "raw", "GSE243076")
soft_path <- file.path(ROOT, "data", "metadata", "GSE243076_family.soft.gz")
raw_tar_path <- file.path(raw_root, "GSE243076_RAW.tar")
atlas_path <- file.path(raw_root, "elife-92046-supp1-v1.xlsx")
resource_paths <- c(soft_path, raw_tar_path, atlas_path)
resource_names <- c("GEO SOFT metadata", "GEO raw 10x archive", "eLife cluster averages")
expected_sha256 <- c(
  "f7faabc68a18ce841106c96f504a2623eaabb2ec4f1a05ca40b1ea2c77f4fe36",
  "8d9d7a225beb4b63a14b64bdb737e81a47f9a22b4049a285d367b243f7635b10",
  "9677e085fa33e3b7b07bdee4d91e0a82044b8a77088237674c2f0f33cf7d31db"
)
if (!all(file.exists(resource_paths))) {
  stop(
    "Missing GSE243076 inputs: ",
    paste(resource_paths[!file.exists(resource_paths)], collapse = ", ")
  )
}
observed_sha256 <- vapply(resource_paths, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1))
resource_manifest <- data.frame(
  resource = resource_names,
  path = normalizePath(resource_paths, winslash = "/"),
  expected_sha256 = expected_sha256,
  observed_sha256 = observed_sha256,
  sha256_pass = tolower(observed_sha256) == expected_sha256,
  stringsAsFactors = FALSE
)
stopifnot(all(resource_manifest$sha256_pass))
write_tsv(
  resource_manifest,
  "results/r/cell_resolved/GSE243076_resource_manifest.tsv"
)

soft <- readLines(gzfile(soft_path), warn = FALSE)
sample_accessions <- sub("^\\^SAMPLE = ", "", soft[grepl("^\\^SAMPLE = ", soft)])
stopifnot(
  length(sample_accessions) == 9L,
  !anyDuplicated(sample_accessions)
)

inner_pattern <- "_filtered_feature_bc_matrix[.]tar[.]gz$"
inner_tars <- list.files(raw_root, pattern = inner_pattern, full.names = TRUE)
if (length(inner_tars) != 9L) {
  utils::untar(raw_tar_path, exdir = raw_root)
  inner_tars <- list.files(raw_root, pattern = inner_pattern, full.names = TRUE)
}
stopifnot(length(inner_tars) == 9L)
inner_accessions <- sub("_.*$", "", basename(inner_tars))
stopifnot(setequal(inner_accessions, sample_accessions))
inner_tars <- inner_tars[match(sample_accessions, inner_accessions)]

audit_10x_donor <- function(tar_path, sample_accession) {
  cat(sample_accession, ": extracting and reading raw 10x matrix\n")
  flush.console()
  extract_root <- file.path(tempdir(), paste0("GSE243076_", sample_accession))
  unlink(extract_root, recursive = TRUE, force = TRUE)
  dir.create(extract_root, recursive = TRUE, showWarnings = FALSE)
  utils::untar(tar_path, exdir = extract_root)
  matrix_path <- list.files(
    extract_root,
    pattern = "^matrix[.]mtx[.]gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  feature_path <- list.files(
    extract_root,
    pattern = "^features[.]tsv[.]gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  barcode_path <- list.files(
    extract_root,
    pattern = "^barcodes[.]tsv[.]gz$",
    recursive = TRUE,
    full.names = TRUE
  )
  stopifnot(
    length(matrix_path) == 1L,
    length(feature_path) == 1L,
    length(barcode_path) == 1L
  )
  matrix <- Matrix::readMM(gzfile(matrix_path))
  matrix <- methods::as(matrix, "CsparseMatrix")
  features <- read.delim(
    gzfile(feature_path),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  barcodes <- readLines(gzfile(barcode_path), warn = FALSE)
  stopifnot(nrow(features) == nrow(matrix), length(barcodes) == ncol(matrix))

  umi <- as.numeric(Matrix::colSums(matrix))
  genes_detected <- diff(matrix@p)
  mitochondrial <- grepl("^MT-", features[[2]])
  mitochondrial_umi <- if (any(mitochondrial)) {
    as.numeric(Matrix::colSums(matrix[mitochondrial, , drop = FALSE]))
  } else {
    rep(0, ncol(matrix))
  }
  mitochondrial_percent <- 100 * mitochondrial_umi / pmax(umi, 1)
  complexity <- log10(pmax(genes_detected, 1)) / log10(pmax(umi, 10))
  threshold_only_qc_proxy <-
    genes_detected >= 200 & genes_detected <= 10000 &
    umi >= 500 & complexity >= 0.8 & mitochondrial_percent <= 10

  result <- data.frame(
    sample_accession = sample_accession,
    raw_archive = basename(tar_path),
    features = nrow(matrix),
    filtered_nuclei = ncol(matrix),
    nonzero_entries = length(matrix@x),
    total_umi = sum(umi),
    median_umi = stats::median(umi),
    median_genes_detected = stats::median(genes_detected),
    median_mitochondrial_percent = stats::median(mitochondrial_percent),
    nuclei_passing_threshold_only_qc_proxy = sum(threshold_only_qc_proxy),
    fraction_passing_threshold_only_qc_proxy = mean(threshold_only_qc_proxy),
    stringsAsFactors = FALSE
  )
  rm(matrix)
  invisible(gc())
  unlink(extract_root, recursive = TRUE, force = TRUE)
  cat(sample_accession, ":", result$filtered_nuclei, "filtered nuclei audited\n")
  flush.console()
  result
}

qc_cache_root <- file.path(
  ROOT,
  "results",
  "r",
  "cell_resolved",
  "GSE243076_raw_qc_cache"
)
dir.create(qc_cache_root, recursive = TRUE, showWarnings = FALSE)
audit_10x_donor_cached <- function(tar_path, sample_accession) {
  cache_path <- file.path(qc_cache_root, paste0(sample_accession, ".tsv"))
  if (file.exists(cache_path)) {
    cached <- read.delim(cache_path, stringsAsFactors = FALSE)
    names(cached)[names(cached) == "nuclei_passing_published_qc_proxy"] <-
      "nuclei_passing_threshold_only_qc_proxy"
    names(cached)[names(cached) == "fraction_passing_published_qc_proxy"] <-
      "fraction_passing_threshold_only_qc_proxy"
    if (
      nrow(cached) == 1L &&
        identical(cached$sample_accession, sample_accession) &&
        identical(cached$raw_archive, basename(tar_path))
    ) {
      cat(sample_accession, ": using cached raw 10x QC\n")
      flush.console()
      return(cached)
    }
  }
  result <- audit_10x_donor(tar_path, sample_accession)
  write.table(
    result,
    cache_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  result
}

raw_qc <- do.call(rbind, Map(
  audit_10x_donor_cached,
  inner_tars,
  sample_accessions
))
rownames(raw_qc) <- NULL
stopifnot(
  nrow(raw_qc) == 9L,
  all(raw_qc$features == 36601L),
  sum(raw_qc$filtered_nuclei) == 74711L,
  sum(raw_qc$nuclei_passing_threshold_only_qc_proxy) == 72834L
)
write_tsv(
  raw_qc,
  "results/r/cell_resolved/GSE243076_raw_10x_donor_qc.tsv"
)

atlas <- as.data.frame(readxl::read_excel(atlas_path, sheet = "information"))
cluster_columns <- paste0("C", 0:20)
stopifnot(
  nrow(atlas) == 36601L,
  identical(names(atlas), c("gene", cluster_columns))
)
for (column in cluster_columns) {
  atlas[[column]] <- as.numeric(atlas[[column]])
}
atlas$source_gene_symbol <- normalize_excel_gene_symbol(atlas$gene)
hgnc <- read.delim(
  file.path(ROOT, "data", "metadata", "hgnc_complete_set.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
resolver <- build_hgnc_resolver(hgnc)
atlas$current_gene_symbol <- resolve_symbols(atlas$source_gene_symbol, resolver)
atlas$current_symbol_preferred <-
  atlas$source_gene_symbol == atlas$current_gene_symbol

other_neurons <- as.matrix(atlas[paste0("C", 0:19)])
c20 <- atlas$C20
atlas$c20_other_neuron_mean <- rowMeans(other_neurons)
atlas$c20_other_neuron_median <- apply(other_neurons, 1L, stats::median)
atlas$c20_other_neuron_max <- apply(other_neurons, 1L, max)
atlas$c20_rank_of_21 <- 1L + rowSums(other_neurons > c20)
atlas$c20_is_max_neuronal_cluster <- atlas$c20_rank_of_21 == 1L
for (pseudocount in c(0.001, 0.01, 0.1)) {
  label <- sub("[.]", "_", format(pseudocount, scientific = FALSE))
  atlas[[paste0("c20_log2_enrichment_pc_", label)]] <- log2(
    (c20 + pseudocount) /
      (atlas$c20_other_neuron_mean + pseudocount)
  )
}
atlas$c20_expression_percentile <-
  (rank(c20, ties.method = "average") - 0.5) / length(c20)
atlas$c20_robust_localization <-
  c20 > 0 & atlas$c20_is_max_neuronal_cluster &
  atlas$c20_log2_enrichment_pc_0_001 >= 1 &
  atlas$c20_log2_enrichment_pc_0_01 >= 1 &
  atlas$c20_log2_enrichment_pc_0_1 >= 1

write_tsv(
  atlas,
  "results/r/cell_resolved/GSE243076_neuronal_cluster_expression.tsv"
)

ordering <- order(
  atlas$current_gene_symbol,
  -as.integer(atlas$current_symbol_preferred),
  -atlas$C20
)
atlas_unique <- atlas[ordering, , drop = FALSE]
atlas_unique <- atlas_unique[
  nzchar(atlas_unique$current_gene_symbol) &
    !duplicated(atlas_unique$current_gene_symbol),
  ,
  drop = FALSE
]
write_tsv(
  atlas_unique,
  "results/r/cell_resolved/GSE243076_C20_motor_neuron_context.tsv"
)

canonical_markers <- c(
  "CHAT", "SLC5A7", "SLC18A3", "MNX1", "ISL1", "ISL2", "PRPH", "NEFH"
)
marker_validation <- atlas_unique[
  match(canonical_markers, atlas_unique$current_gene_symbol),
  c(
    "current_gene_symbol", "C20", "c20_other_neuron_mean",
    "c20_rank_of_21", "c20_log2_enrichment_pc_0_01",
    "c20_robust_localization"
  ),
  drop = FALSE
]
stopifnot(
  all(marker_validation$current_gene_symbol == canonical_markers),
  all(marker_validation$C20 > 0),
  all(marker_validation$c20_log2_enrichment_pc_0_01[1:6] > 0)
)
write_tsv(
  marker_validation,
  "results/r/cell_resolved/GSE243076_C20_marker_validation.tsv"
)

summary <- data.frame(
  metric = c(
    "independent_donors", "deposited_filtered_nuclei",
    "nuclei_passing_threshold_only_qc_proxy",
    "author_reported_post_qc_nuclei", "features",
    "neuronal_clusters", "author_motor_neuron_cluster",
    "canonical_markers_detected", "canonical_markers_c20_enriched",
    "inferential_role"
  ),
  value = c(
    9L, sum(raw_qc$filtered_nuclei),
    sum(raw_qc$nuclei_passing_threshold_only_qc_proxy), 64021L,
    nrow(atlas), 21L, "C20",
    sum(marker_validation$C20 > 0),
    sum(marker_validation$c20_log2_enrichment_pc_0_01 > 0),
    "independent adult motor-neuron localization; no SMA contrast"
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "results/r/cell_resolved/GSE243076_summary.tsv")

top <- atlas_unique[
  atlas_unique$C20 > 0 & is.finite(atlas_unique$c20_log2_enrichment_pc_0_01),
  ,
  drop = FALSE
]
top <- head(top[order(-top$c20_log2_enrichment_pc_0_01), ], 20L)
png(
  file.path(ROOT, "results", "r", "figures", "GSE243076_C20_top_markers.png"),
  width = 1700,
  height = 1400,
  res = 180
)
par(mar = c(5, 10, 4, 2))
barplot(
  rev(top$c20_log2_enrichment_pc_0_01),
  names.arg = rev(top$current_gene_symbol),
  horiz = TRUE,
  las = 1,
  col = "#167D8D",
  border = NA,
  xlab = "C20 log2 enrichment over other neuronal clusters",
  main = "Independent adult human spinal motor-neuron cluster"
)
dev.off()

write_session_info(
  "results/r/cell_resolved/GSE243076_R_sessionInfo.txt"
)
cat(
  "GSE243076 independent atlas complete:",
  nrow(raw_qc), "donors,", sum(raw_qc$filtered_nuclei),
  "deposited filtered nuclei,",
  sum(raw_qc$nuclei_passing_threshold_only_qc_proxy),
  "passing the threshold-only QC proxy (before unreproduced doublet removal), and",
  nrow(atlas),
  "genes across 21 neuronal clusters\n"
)
