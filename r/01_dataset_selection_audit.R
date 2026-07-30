source(file.path("r", "common.R"))

selection <- read.delim(
  file.path(ROOT, "config", "dataset_selection.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

core <- selection[selection$decision == "CORE", , drop = FALSE]
external <- selection[
  selection$decision == "EXTERNAL_VALIDATION",
  ,
  drop = FALSE
]
sensitivity <- selection[
  selection$decision == "EXTERNAL_SENSITIVITY_VALIDATION",
  ,
  drop = FALSE
]
stopifnot(
  identical(sort(core$accession), sort(c("GSE93939", "GSE290979"))),
  all(core$organism == "Homo sapiens"),
  nrow(external) == 1L,
  identical(external$accession, "GSE108094"),
  identical(external$organism, "Homo sapiens"),
  nrow(sensitivity) == 1L,
  identical(sensitivity$accession, "GSE69175"),
  identical(sensitivity$organism, "Homo sapiens")
)

audit_rows <- list()
group_rows <- list()
for (accession in core$accession) {
  counts <- read_count_matrix(accession)
  metadata <- read_sample_metadata(accession)
  assert_identical_samples(counts, metadata, accession)
  metadata <- metadata[colnames(counts), , drop = FALSE]

  expected_libraries <- if (accession == "GSE93939") 39L else 31L
  stopifnot(ncol(counts) == expected_libraries)

  if (accession == "GSE93939") {
    primary <- metadata$group %in% c("OMN", "SC")
    independent_units <- length(unique(metadata$donor_id[primary]))
    primary_libraries <- sum(primary)
    grouped <- aggregate(
      donor_id ~ group,
      metadata,
      function(values) length(unique(values))
    )
    library_counts <- as.data.frame(table(metadata$group))
    names(library_counts) <- c("group", "libraries")
    grouped <- merge(library_counts, grouped, by = "group", all.x = TRUE)
    names(grouped)[3] <- "independent_units"
  } else {
    independent_units <- length(unique(metadata$`cell line`))
    primary_libraries <- nrow(metadata)
    library_counts <- as.data.frame(table(metadata$analysis_group))
    names(library_counts) <- c("group", "libraries")
    line_counts <- aggregate(
      metadata$`cell line`,
      list(group = metadata$analysis_group),
      function(values) length(unique(values))
    )
    names(line_counts)[2] <- "independent_units"
    grouped <- merge(library_counts, line_counts, by = "group", all.x = TRUE)
  }
  grouped$accession <- accession
  group_rows[[accession]] <- grouped[, c(
    "accession",
    "group",
    "libraries",
    "independent_units"
  )]

  audit_rows[[accession]] <- data.frame(
    accession = accession,
    organism = core$organism[core$accession == accession],
    matrix_features = nrow(counts),
    libraries = ncol(counts),
    primary_contrast_libraries = primary_libraries,
    primary_independent_units = independent_units,
    count_metadata_exact_match = TRUE,
    selected_role = core$role[core$accession == accession],
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, audit_rows)
groups <- do.call(rbind, group_rows)
rownames(audit) <- NULL
rownames(groups) <- NULL

write_tsv(
  selection,
  "results/r/dataset_audit/dataset_selection_decisions.tsv"
)
write_tsv(
  audit,
  "results/r/dataset_audit/core_dataset_audit.tsv"
)
write_tsv(
  groups,
  "results/r/dataset_audit/core_dataset_group_counts.tsv"
)
write_session_info("results/r/dataset_audit/R_sessionInfo.txt")

cat(
  "Dataset audit complete:",
  paste(audit$accession, audit$libraries, sep = "="),
  "libraries; independent units",
  paste(audit$accession, audit$primary_independent_units, sep = "="),
  "\n"
)
