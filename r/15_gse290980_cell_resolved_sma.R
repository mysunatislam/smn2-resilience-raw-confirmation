source(file.path("r", "common.R"))

for (package in c("readxl", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package is unavailable: ", package)
  }
}

raw_root <- file.path(ROOT, "data", "raw", "GSE290980")
soft_path <- file.path(ROOT, "data", "metadata", "GSE290980_family.soft.gz")
marker_path <- file.path(raw_root, "supplementary_data_4.xlsx")
deg_path <- file.path(raw_root, "supplementary_data_5.xlsx")
description_path <- file.path(raw_root, "supplementary_file_description.pdf")
resource_paths <- c(soft_path, marker_path, deg_path, description_path)
expected_sha256 <- c(
  "54a0c9aea94bbba9ebfa535b955155e5e495921c1d3dc764fdcd08bc3f3503f2",
  "0641848d6d332d50202e400f81ef19952ce024097045b57fa025b812463489d6",
  "83e04fe3bf9c35fd40a941d7f3d69e417e8174a2f3b8c0adae5d79109a887df3",
  "567a23de33c8a79adea4b2c2d2aee2bca888b652184aaa473aff8b655693f2fe"
)
if (!all(file.exists(resource_paths))) {
  stop(
    "Missing GSE290980 inputs: ",
    paste(resource_paths[!file.exists(resource_paths)], collapse = ", ")
  )
}
observed_sha256 <- vapply(resource_paths, function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}, character(1))
resource_manifest <- data.frame(
  resource = c(
    "GEO SOFT metadata", "Supplementary Data 4 markers",
    "Supplementary Data 5 pseudobulk DEGs", "supplementary file description"
  ),
  path = normalizePath(resource_paths, winslash = "/"),
  expected_sha256 = expected_sha256,
  observed_sha256 = observed_sha256,
  sha256_pass = tolower(observed_sha256) == expected_sha256,
  stringsAsFactors = FALSE
)
stopifnot(all(resource_manifest$sha256_pass))
write_tsv(
  resource_manifest,
  "results/r/cell_resolved/GSE290980_resource_manifest.tsv"
)

soft <- readLines(gzfile(soft_path), warn = FALSE)
sample_accessions <- sub("^\\^SAMPLE = ", "", soft[grepl("^\\^SAMPLE = ", soft)])
sample_lines <- sub(
  "^!Sample_characteristics_ch1 = cell line: ",
  "",
  soft[grepl("^!Sample_characteristics_ch1 = cell line: ", soft)]
)
sample_genotypes <- trimws(sub(
  "^!Sample_characteristics_ch1 = genotype: ",
  "",
  soft[grepl("^!Sample_characteristics_ch1 = genotype: ", soft)]
))
stopifnot(
  length(sample_accessions) == 8L,
  length(sample_lines) == 8L,
  length(sample_genotypes) == 8L,
  length(unique(sample_lines)) == 4L,
  identical(sort(unique(sample_genotypes)), c("CTRL", "SMA"))
)
sample_metadata <- data.frame(
  sample_accession = sample_accessions,
  donor_line = sample_lines,
  genotype = sample_genotypes,
  independent_donor_line = !duplicated(sample_lines),
  stringsAsFactors = FALSE
)
write_tsv(
  sample_metadata,
  "results/r/cell_resolved/GSE290980_sample_metadata.tsv"
)

hgnc <- read.delim(
  file.path(ROOT, "data", "metadata", "hgnc_complete_set.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
resolver <- build_hgnc_resolver(hgnc)

marker_columns <- c(
  "worksheet_cell_type", "cluster", "source_gene_symbol", "avg_log2fc",
  "p_value", "q_value", "blank", "summary_cluster", "summary_cell_type"
)
markers <- as.data.frame(readxl::read_excel(
  marker_path,
  sheet = "Suppl Data 4 ",
  skip = 3,
  col_names = marker_columns
))
cluster_map <- unique(markers[
  grepl("^[0-9]+$", as.character(markers$summary_cluster)),
  c("summary_cluster", "summary_cell_type"),
  drop = FALSE
])
cluster_map$summary_cluster <- as.integer(cluster_map$summary_cluster)
markers <- markers[
  grepl("^[0-9]+$", as.character(markers$cluster)) &
    !is.na(markers$source_gene_symbol),
  c("cluster", "source_gene_symbol", "avg_log2fc", "p_value", "q_value"),
  drop = FALSE
]
markers$cluster <- as.integer(markers$cluster)
markers$cell_type <- cluster_map$summary_cell_type[
  match(markers$cluster, cluster_map$summary_cluster)
]
markers$source_gene_symbol <- normalize_excel_gene_symbol(
  markers$source_gene_symbol
)
markers$current_gene_symbol <- resolve_symbols(
  markers$source_gene_symbol,
  resolver
)
markers$current_symbol_preferred <-
  markers$source_gene_symbol == markers$current_gene_symbol
for (column in c("avg_log2fc", "p_value", "q_value")) {
  markers[[column]] <- as.numeric(markers[[column]])
}
marker_counts <- table(factor(markers$cluster, levels = 0:18))
stopifnot(
  nrow(cluster_map) == 19L,
  nrow(markers) == 5700L,
  all(marker_counts == 300L),
  identical(sort(cluster_map$summary_cluster[cluster_map$summary_cell_type == "MN"]),
    c(6L, 7L, 10L, 18L))
)
write_tsv(markers, "results/r/cell_resolved/GSE290980_cluster_markers.tsv")
write_tsv(cluster_map, "results/r/cell_resolved/GSE290980_cluster_annotation.tsv")

mn_markers <- markers[markers$cell_type == "MN", , drop = FALSE]
mn_marker_split <- split(mn_markers, mn_markers$current_gene_symbol)
mn_marker_summary <- do.call(rbind, lapply(mn_marker_split, function(frame) {
  data.frame(
    current_gene_symbol = frame$current_gene_symbol[1],
    gse290980_mn_marker_cluster_count = length(unique(frame$cluster)),
    gse290980_mn_marker_clusters = paste(sort(unique(frame$cluster)), collapse = ","),
    gse290980_mn_marker_max_log2fc = max(frame$avg_log2fc, na.rm = TRUE),
    gse290980_mn_marker_min_q_value = min(frame$q_value, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(mn_marker_summary) <- NULL
write_tsv(
  mn_marker_summary,
  "results/r/cell_resolved/GSE290980_motor_neuron_marker_summary.tsv"
)

deg_wide <- as.data.frame(readxl::read_excel(
  deg_path,
  sheet = "Suppl Data 5 - DEGs",
  col_names = FALSE
))
block_starts <- seq.int(1L, ncol(deg_wide), by = 4L)
deg_frames <- lapply(block_starts, function(start) {
  cell_type <- as.character(deg_wide[[start]][1])
  frame <- data.frame(
    cell_type = cell_type,
    source_gene_symbol = as.character(deg_wide[[start]][-(1:2)]),
    p_value = suppressWarnings(as.numeric(deg_wide[[start + 1L]][-(1:2)])),
    q_value = suppressWarnings(as.numeric(deg_wide[[start + 2L]][-(1:2)])),
    published_log2fc_ctrl_vs_sma = suppressWarnings(
      as.numeric(deg_wide[[start + 3L]][-(1:2)])
    ),
    stringsAsFactors = FALSE
  )
  frame <- frame[
    !is.na(frame$source_gene_symbol) & nzchar(frame$source_gene_symbol),
    ,
    drop = FALSE
  ]
  frame$source_gene_symbol <- normalize_excel_gene_symbol(
    frame$source_gene_symbol
  )
  frame$current_gene_symbol <- resolve_symbols(
    frame$source_gene_symbol,
    resolver
  )
  frame$current_symbol_preferred <-
    frame$source_gene_symbol == frame$current_gene_symbol
  # The workbook legend defines positive published values as CTRL > SMA.
  frame$sma_vs_control_log2_effect <-
    -frame$published_log2fc_ctrl_vs_sma
  frame
})
deg <- do.call(rbind, deg_frames)
rownames(deg) <- NULL
expected_counts <- c(
  MN = 2031L, pFP = 834L, Progenitor = 385L, Astroglia = 303L,
  V2b = 141L, V2a = 128L, pMN = 121L
)
observed_counts <- table(factor(deg$cell_type, levels = names(expected_counts)))
stopifnot(
  identical(as.integer(observed_counts), as.integer(expected_counts)),
  all(deg$q_value < 0.05),
  all(abs(deg$published_log2fc_ctrl_vs_sma) > 0.5)
)
ordering <- order(
  deg$cell_type,
  deg$current_gene_symbol,
  deg$q_value,
  -as.integer(deg$current_symbol_preferred)
)
deg <- deg[ordering, , drop = FALSE]
deg <- deg[!duplicated(deg[c("cell_type", "current_gene_symbol")]), , drop = FALSE]
write_tsv(deg, "results/r/cell_resolved/GSE290980_cell_type_pseudobulk_DEGs.tsv")

mn_deg <- deg[deg$cell_type == "MN", , drop = FALSE]
mn_deg <- merge(
  mn_deg,
  mn_marker_summary,
  by = "current_gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
write_tsv(
  mn_deg,
  "results/r/cell_resolved/GSE290980_motor_neuron_pseudobulk_DEGs.tsv"
)

deg_summary <- data.frame(
  cell_type = names(expected_counts),
  significant_degs = as.integer(observed_counts),
  sma_higher = vapply(names(expected_counts), function(cell_type) {
    sum(deg$sma_vs_control_log2_effect[deg$cell_type == cell_type] > 0)
  }, integer(1)),
  sma_lower = vapply(names(expected_counts), function(cell_type) {
    sum(deg$sma_vs_control_log2_effect[deg$cell_type == cell_type] < 0)
  }, integer(1)),
  stringsAsFactors = FALSE
)
write_tsv(
  deg_summary,
  "results/r/cell_resolved/GSE290980_cell_type_DEG_summary.tsv"
)

canonical_markers <- c(
  "CHAT", "SLC5A7", "SLC18A3", "MNX1", "ISL1", "ISL2", "PRPH", "NEFH"
)
marker_validation <- data.frame(
  current_gene_symbol = canonical_markers,
  present_in_any_mn_top300 = canonical_markers %in%
    mn_marker_summary$current_gene_symbol,
  mn_marker_cluster_count = mn_marker_summary$gse290980_mn_marker_cluster_count[
    match(canonical_markers, mn_marker_summary$current_gene_symbol)
  ],
  stringsAsFactors = FALSE
)
write_tsv(
  marker_validation,
  "results/r/cell_resolved/GSE290980_motor_neuron_marker_validation.tsv"
)

summary <- data.frame(
  metric = c(
    "geo_libraries", "donor_lines", "control_donor_lines",
    "sma_donor_lines", "reported_cells", "clusters", "motor_neuron_clusters",
    "motor_neuron_pseudobulk_degs", "analysis_role", "independence_status"
  ),
  value = c(
    8L, 4L, 2L, 2L, 27114L, 19L, 4L, nrow(mn_deg),
    "cell-resolved SMA disease localization",
    "same study as GSE290979; not independent validation"
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "results/r/cell_resolved/GSE290980_summary.tsv")

top <- head(mn_deg[order(-abs(mn_deg$sma_vs_control_log2_effect)), ], 20L)
png(
  file.path(ROOT, "results", "r", "figures", "GSE290980_MN_top_disease_DEGs.png"),
  width = 1750,
  height = 1400,
  res = 180
)
par(mar = c(5, 11, 4, 2))
barplot(
  rev(top$sma_vs_control_log2_effect),
  names.arg = rev(top$current_gene_symbol),
  horiz = TRUE,
  las = 1,
  col = ifelse(
    rev(top$sma_vs_control_log2_effect) > 0,
    "#B44B4B",
    "#2878A5"
  ),
  border = NA,
  xlab = "SMA minus control log2 effect",
  main = "GSE290980 motor-neuron pseudobulk disease effects"
)
abline(v = 0, col = "grey30")
dev.off()

write_session_info(
  "results/r/cell_resolved/GSE290980_R_sessionInfo.txt"
)
cat(
  "GSE290980 cell-resolved context complete:",
  nrow(mn_deg), "motor-neuron DEGs across", length(sample_accessions),
  "libraries from", length(unique(sample_lines)), "donor lines\n"
)
