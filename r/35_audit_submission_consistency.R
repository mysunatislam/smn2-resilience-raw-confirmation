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

read_tsv <- function(path) {
  read.delim(
    file.path(ROOT, path),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

as_flag <- function(x) {
  !is.na(x) & toupper(as.character(x)) == "TRUE"
}

checks <- list()
add_check <- function(id, category, source, expected, observed, pass, details = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id,
    category = category,
    source = source,
    expected = as.character(expected),
    observed = as.character(observed),
    pass = isTRUE(pass),
    details = details,
    stringsAsFactors = FALSE
  )
}

summary_table <- read_tsv("manuscript/supplementary_table_S2_summary.tsv")
summary_values <- setNames(summary_table$value, summary_table$metric)
candidates <- read_tsv("manuscript/supplementary_table_S1_candidates.tsv")
permutation <- read_tsv("manuscript/supplementary_table_S10_permutation_summary.tsv")
bootstrap <- read_tsv("manuscript/supplementary_table_S14_bootstrap_candidate_stability.tsv")
splicing <- read_tsv("manuscript/supplementary_table_S4_raw_splice_confirmation.tsv")
sensitivity <- read_tsv("manuscript/supplementary_table_S12_score_sensitivity.tsv")
manuscript <- paste(
  readLines(file.path(ROOT, "manuscript", "manuscript.md"), warn = FALSE),
  collapse = "\n"
)
manuscript_normalized <- trimws(gsub("\\s+", " ", manuscript))

eligible <- unique(as.integer(sensitivity$eligible_genes))
add_check(
  "ranked_gene_universe", "numerical", "Supplementary Table S12",
  11326L, paste(eligible, collapse = ","),
  length(eligible) == 1L && identical(eligible, 11326L)
)

add_check(
  "frozen_candidate_count", "numerical", "Supplementary Tables S1 and S2",
  37L, nrow(candidates),
  nrow(candidates) == 37L && as.integer(summary_values[["frozen_candidate_count"]]) == 37L
)

raw_full <- sum(as_flag(candidates$local9_raw_full_directional_pattern))
add_check(
  "raw_full_direction_count", "numerical", "Supplementary Tables S1 and S2",
  21L, raw_full,
  raw_full == 21L && as.integer(summary_values[["raw_full_directional_pattern"]]) == 21L
)

raw_robust <- sum(as_flag(candidates$local9_raw_direction_and_unit_robust))
add_check(
  "raw_unit_robust_count", "numerical", "Supplementary Tables S1 and S2",
  9L, raw_robust,
  raw_robust == 9L && as.integer(summary_values[["raw_unit_robust_candidates"]]) == 9L
)

both_robust <- sum(
  as_flag(candidates$local9_raw_direction_and_unit_robust) &
    as_flag(candidates$all_estimable_biological_unit_checks_robust)
)
tier_1 <- candidates$gene_symbol[
  candidates$publication_evidence_tier == "tier_1_raw_and_processed_unit_robust"
]
expected_tier_1 <- c("LY6H", "HS3ST5", "ZNF853", "IL17D")
add_check(
  "tier_1_count_and_identity", "numerical", "Supplementary Tables S1 and S2",
  paste(expected_tier_1, collapse = ","), paste(tier_1, collapse = ","),
  both_robust == 4L && identical(tier_1, expected_tier_1) &&
    as.integer(summary_values[["raw_and_processed_unit_robust_candidates"]]) == 4L
)

candidate_p <- permutation$empirical_p_value[
  permutation$test == "frozen_37_complete_pattern"
]
overlap_p <- permutation$empirical_p_value[
  permutation$test == "raw_processed_robust_overlap_4"
]
add_check(
  "candidate_count_permutation_p", "numerical", "Supplementary Table S10",
  0.835716428357164, candidate_p,
  length(candidate_p) == 1L && isTRUE(all.equal(candidate_p, 0.835716428357164))
)
add_check(
  "robustness_overlap_permutation_p", "numerical", "Supplementary Table S10",
  0.0440955904409559, overlap_p,
  length(overlap_p) == 1L && isTRUE(all.equal(overlap_p, 0.0440955904409559))
)

tier_1_boot <- bootstrap$selection_frequency[
  match(expected_tier_1, bootstrap$gene_symbol)
]
add_check(
  "tier_1_bootstrap_range", "numerical", "Supplementary Table S14",
  "0.898-0.959", sprintf("%.3f-%.3f", min(tier_1_boot), max(tier_1_boot)),
  length(tier_1_boot) == 4L && !anyNA(tier_1_boot) &&
    min(tier_1_boot) == 0.898 && max(tier_1_boot) == 0.959
)

recovered <- sum(as_flag(splicing$structurally_recovered))
confirmed <- sum(as_flag(splicing$strong_raw_confirmation))
add_check(
  "raw_splice_exact_recovery", "numerical", "Supplementary Table S4",
  "6/83", paste0(recovered, "/", nrow(splicing)),
  nrow(splicing) == 83L && recovered == 6L
)
add_check(
  "raw_splice_full_confirmation", "numerical", "Supplementary Table S4",
  0L, confirmed,
  confirmed == 0L
)

required_text <- c(
  ranked_to_frozen = "11,326 ranked approved genes, 37",
  raw_direction_and_robustness = "21 recovered the full raw direction pattern, and nine",
  permutation_p = "p = 0.836",
  overlap_p = "p = 0.044",
  bootstrap_range = "89.8-95.9%",
  splice_recovery = "six of 83 frozen events",
  splice_failure = "No event met the complete"
)
for (name in names(required_text)) {
  phrase <- required_text[[name]]
  found <- grepl(phrase, manuscript_normalized, fixed = TRUE)
  add_check(
    paste0("manuscript_", name), "manuscript", "manuscript/manuscript.md",
    phrase, if (found) phrase else "NOT FOUND",
    found
  )
}

audit <- do.call(rbind, checks)
output_dir <- file.path(ROOT, "manuscript", "submission")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  audit,
  file.path(output_dir, "numerical_consistency_audit.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

if (!all(audit$pass)) {
  failed <- audit$check_id[!audit$pass]
  stop("Submission consistency audit failed: ", paste(failed, collapse = ", "))
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
  file.path(output_dir, "NUMERICAL_CONSISTENCY_AUDIT_COMPLETE.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Submission numerical consistency audit passed (", nrow(audit), " checks)\n", sep = "")
