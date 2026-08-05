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

read_table <- function(filename, directory = "manuscript") {
  path <- file.path(ROOT, directory, filename)
  if (!file.exists(path)) {
    stop("Missing stability output: ", path)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
assert_close <- function(actual, expected, tolerance = 1e-12) {
  if (
    length(actual) != length(expected) ||
      any(!is.finite(actual)) ||
      any(abs(actual - expected) > tolerance)
  ) {
    stop("Unexpected numeric result")
  }
}

permutation <- read_table(
  "supplementary_table_S10_permutation_summary.tsv"
)
stopifnot(
  identical(
    permutation$test,
    c(
      "frozen_37_complete_pattern",
      "raw_processed_robust_overlap_4"
    )
  ),
  identical(permutation$observed, c(37L, 4L)),
  identical(permutation$permutations, c(10000L, 10000L))
)
assert_close(
  permutation$empirical_p_value,
  c(0.835716428357164, 0.0440955904409559)
)

negative_controls <- read_table(
  "supplementary_table_S11_negative_controls.tsv"
)
stopifnot(
  identical(negative_controls$selected_genes, c(37L, 36L, 69L, 468L))
)

score_sensitivity <- read_table(
  "supplementary_table_S12_score_sensitivity.tsv"
)
stopifnot(
  nrow(score_sensitivity) == 9L,
  all(score_sensitivity$eligible_genes == 11326L),
  min(score_sensitivity$spearman_rank_vs_equal) > 0.81,
  max(score_sensitivity$spearman_rank_vs_equal) == 1
)

rank_stability <- read_table(
  "supplementary_table_S13_candidate_rank_stability.tsv"
)
stopifnot(
  nrow(rank_stability) == 37L,
  !anyDuplicated(rank_stability$gene_symbol),
  all(rank_stability$variants_tested == 9L)
)
tier1 <- c("LY6H", "HS3ST5", "ZNF853", "IL17D")
tier1_rank <- rank_stability[match(tier1, rank_stability$gene_symbol), ]
stopifnot(
  identical(tier1_rank$variants_in_top100, c(9L, 9L, 5L, 4L))
)

bootstrap <- read_table(
  "supplementary_table_S14_bootstrap_candidate_stability.tsv"
)
stopifnot(
  nrow(bootstrap) == 37L,
  !anyDuplicated(bootstrap$gene_symbol),
  all(bootstrap$iterations == 1000L)
)
tier1_bootstrap <- bootstrap[
  match(tier1, bootstrap$gene_symbol),
  ,
  drop = FALSE
]
assert_close(
  tier1_bootstrap$selection_frequency,
  c(0.898, 0.959, 0.949, 0.941)
)

bootstrap_summary <- read_table(
  "supplementary_table_S15_bootstrap_summary.tsv"
)
summary_values <- setNames(
  bootstrap_summary$value,
  bootstrap_summary$metric
)
stopifnot(
  summary_values[["iterations"]] == "1000",
  summary_values[["median_complete_pattern_genes"]] == "38",
  summary_values[["median_frozen_37_selected"]] == "24"
)

permutation_audit <- read_table(
  "cross_model_permutation_null.tsv",
  file.path("manuscript", "audit")
)
bootstrap_audit <- read_table(
  "biological_unit_bootstrap_iterations.tsv",
  file.path("manuscript", "audit")
)
provenance <- read_table(
  "biological_unit_bootstrap_provenance.tsv",
  file.path("manuscript", "audit")
)
stopifnot(
  nrow(permutation_audit) == 10000L,
  nrow(bootstrap_audit) == 1000L,
  provenance$value[provenance$metric == "random_library_split"] == "FALSE",
  provenance$value[provenance$metric == "external_model"] ==
    "GSE108094_effect_fixed_not_unit_resampled"
)

png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
for (filename in c(
  "supplementary_figure_S24_permutation_null.png",
  "supplementary_figure_S25_bootstrap_stability.png"
)) {
  path <- file.path(ROOT, "manuscript", "figures", filename)
  connection <- file(path, open = "rb")
  signature <- readBin(connection, what = "raw", n = 8L)
  close(connection)
  stopifnot(identical(signature, png_signature))
}
for (filename in c(
  "supplementary_figure_S24_permutation_null.pdf",
  "supplementary_figure_S25_bootstrap_stability.pdf"
)) {
  path <- file.path(ROOT, "manuscript", "figures", filename)
  connection <- file(path, open = "rb")
  signature <- rawToChar(readBin(connection, what = "raw", n = 5L))
  close(connection)
  stopifnot(identical(signature, "%PDF-"))
}

for (filename in c(
  "CROSS_MODEL_PERMUTATION_SENSITIVITY_COMPLETE.tsv",
  "BIOLOGICAL_UNIT_BOOTSTRAP_COMPLETE.tsv"
)) {
  path <- file.path(ROOT, "manuscript", "audit", filename)
  stopifnot(file.exists(path), file.info(path)$size > 0)
}

cat("Cross-model stability publication test passed\n")
