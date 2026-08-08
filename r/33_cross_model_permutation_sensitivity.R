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

integration_path <- argument_value("--integration-table=")
candidate_path <- argument_value(
  "--candidate-table=",
  file.path(ROOT, "manuscript", "supplementary_table_S1_candidates.tsv")
)
output_root <- argument_value("--output-root=")
permutations <- as.integer(argument_value("--permutations=", "10000"))
seed <- as.integer(argument_value("--seed=", "290979"))
if (
  is.null(integration_path) ||
    is.null(output_root) ||
    !is.finite(permutations) ||
    permutations < 1000L
) {
  stop(
    "Required: --integration-table=PATH --output-root=PATH ",
    "[--permutations=10000] [--seed=290979]"
  )
}

integration_path <- normalizePath(
  integration_path,
  winslash = "/",
  mustWork = TRUE
)
candidate_path <- normalizePath(
  candidate_path,
  winslash = "/",
  mustWork = TRUE
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)

integration <- read.delim(
  integration_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
candidates <- read.delim(
  candidate_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_integration <- c(
  "gene_symbol",
  "cross_model_ranking_eligible",
  "exploratory_natural_resistance_support",
  "omn_log2_effect_omn_vs_sc",
  "disease_sma_vs_control_log2_effect",
  "treatment_r6_vs_scramble_log2_effect",
  "gse108094_sma_vs_control_log2_effect",
  "natural_resistance_percentile",
  "gse290979_sma_depletion_percentile",
  "gse290979_r6_alignment_percentile",
  "gse108094_sma_depletion_percentile"
)
required_candidates <- c(
  "gene_symbol",
  "all_estimable_biological_unit_checks_robust",
  "local9_raw_direction_and_unit_robust"
)
stopifnot(
  all(required_integration %in% names(integration)),
  all(required_candidates %in% names(candidates)),
  !anyDuplicated(integration$gene_symbol),
  !anyDuplicated(candidates$gene_symbol),
  nrow(candidates) == 37L
)

as_flag <- function(values) {
  !is.na(values) & as.logical(values)
}
eligible <- as_flag(integration$cross_model_ranking_eligible)
universe <- integration[eligible, , drop = FALSE]
stopifnot(nrow(universe) == 11326L)

natural <- universe$omn_log2_effect_omn_vs_sc
disease <- universe$disease_sma_vs_control_log2_effect
treatment <- universe$treatment_r6_vs_scramble_log2_effect
external <- universe$gse108094_sma_vs_control_log2_effect
natural_support <- as_flag(
  universe$exploratory_natural_resistance_support
)

selection_count <- function(disease_values, treatment_values, external_values) {
  sum(
    natural_support &
      natural > 0 &
      disease_values < 0 &
      treatment_values > 0 &
      external_values < 0
  )
}
observed_candidate_count <- selection_count(disease, treatment, external)
stopifnot(observed_candidate_count == 37L)

set.seed(seed)
permutation_counts <- integer(permutations)
for (iteration in seq_len(permutations)) {
  permutation_counts[[iteration]] <- selection_count(
    sample(disease, replace = FALSE),
    sample(treatment, replace = FALSE),
    sample(external, replace = FALSE)
  )
}
candidate_empirical_p <- (
  1 + sum(permutation_counts >= observed_candidate_count)
) / (permutations + 1)

processed_robust <- as_flag(
  candidates$all_estimable_biological_unit_checks_robust
)
raw_robust <- as_flag(candidates$local9_raw_direction_and_unit_robust)
observed_robust_overlap <- sum(processed_robust & raw_robust)
stopifnot(
  sum(processed_robust) == 7L,
  sum(raw_robust) == 9L,
  observed_robust_overlap == 4L
)
robust_overlap_null <- integer(permutations)
for (iteration in seq_len(permutations)) {
  robust_overlap_null[[iteration]] <- sum(
    processed_robust & sample(raw_robust, replace = FALSE)
  )
}
robust_overlap_empirical_p <- (
  1 + sum(robust_overlap_null >= observed_robust_overlap)
) / (permutations + 1)

permutation_summary <- data.frame(
  test = c(
    "frozen_37_complete_pattern",
    "raw_processed_robust_overlap_4"
  ),
  observed = c(observed_candidate_count, observed_robust_overlap),
  permutations = permutations,
  null_mean = c(mean(permutation_counts), mean(robust_overlap_null)),
  null_median = c(
    median(permutation_counts),
    median(robust_overlap_null)
  ),
  null_q025 = c(
    unname(quantile(permutation_counts, 0.025)),
    unname(quantile(robust_overlap_null, 0.025))
  ),
  null_q975 = c(
    unname(quantile(permutation_counts, 0.975)),
    unname(quantile(robust_overlap_null, 0.975))
  ),
  null_at_least_observed = c(
    sum(permutation_counts >= observed_candidate_count),
    sum(robust_overlap_null >= observed_robust_overlap)
  ),
  empirical_p_value = c(
    candidate_empirical_p,
    robust_overlap_empirical_p
  ),
  stringsAsFactors = FALSE
)

negative_controls <- data.frame(
  control = c(
    "observed_complete_pattern",
    "both_disease_directions_reversed",
    "treatment_direction_reversed",
    "full_inverse_pattern_without_natural_support_filter"
  ),
  selected_genes = c(
    observed_candidate_count,
    sum(
      natural_support & natural > 0 & disease > 0 &
        treatment > 0 & external > 0
    ),
    sum(
      natural_support & natural > 0 & disease < 0 &
        treatment < 0 & external < 0
    ),
    sum(natural < 0 & disease > 0 & treatment < 0 & external > 0)
  ),
  interpretation = c(
    "frozen_selection_rule",
    "negative_control",
    "negative_control",
    "negative_control"
  ),
  stringsAsFactors = FALSE
)

components <- cbind(
  natural = universe$natural_resistance_percentile,
  organoid_disease = universe$gse290979_sma_depletion_percentile,
  treatment = universe$gse290979_r6_alignment_percentile,
  external_disease = universe$gse108094_sma_depletion_percentile
)
stopifnot(all(is.finite(components)))

weighted_score <- function(weights) {
  as.numeric(components %*% weights / sum(weights))
}
variant_scores <- list(
  equal_weights = weighted_score(c(1, 1, 1, 1)),
  external_double_weight = weighted_score(c(1, 1, 1, 2)),
  natural_double_weight = weighted_score(c(2, 1, 1, 1)),
  treatment_double_weight = weighted_score(c(1, 1, 2, 1)),
  omit_natural = weighted_score(c(0, 1, 1, 1)),
  omit_organoid_disease = weighted_score(c(1, 0, 1, 1)),
  omit_treatment = weighted_score(c(1, 1, 0, 1)),
  omit_external = weighted_score(c(1, 1, 1, 0)),
  median_percentile = apply(components, 1L, median)
)
variant_ranks <- lapply(
  variant_scores,
  function(score) rank(-score, ties.method = "average")
)
baseline_rank <- variant_ranks$equal_weights
baseline_top100 <- order(baseline_rank)[seq_len(100L)]

score_summary <- do.call(rbind, lapply(
  names(variant_scores),
  function(variant) {
    current_rank <- variant_ranks[[variant]]
    current_top100 <- order(current_rank)[seq_len(100L)]
    data.frame(
      score_variant = variant,
      eligible_genes = nrow(universe),
      spearman_rank_vs_equal = suppressWarnings(cor(
        baseline_rank,
        current_rank,
        method = "spearman"
      )),
      top100_overlap_with_equal = length(intersect(
        baseline_top100,
        current_top100
      )),
      stringsAsFactors = FALSE
    )
  }
))

candidate_index <- match(candidates$gene_symbol, universe$gene_symbol)
stopifnot(!anyNA(candidate_index))
candidate_rank_long <- do.call(rbind, lapply(
  names(variant_scores),
  function(variant) {
    data.frame(
      gene_symbol = candidates$gene_symbol,
      score_variant = variant,
      score = variant_scores[[variant]][candidate_index],
      eligible_universe_rank = variant_ranks[[variant]][candidate_index],
      stringsAsFactors = FALSE
    )
  }
))
candidate_rank_stability <- do.call(rbind, lapply(
  split(candidate_rank_long, candidate_rank_long$gene_symbol),
  function(frame) {
    data.frame(
      gene_symbol = frame$gene_symbol[[1L]],
      median_rank = median(frame$eligible_universe_rank),
      minimum_rank = min(frame$eligible_universe_rank),
      maximum_rank = max(frame$eligible_universe_rank),
      rank_range = diff(range(frame$eligible_universe_rank)),
      variants_in_top100 = sum(frame$eligible_universe_rank <= 100L),
      variants_tested = nrow(frame),
      stringsAsFactors = FALSE
    )
  }
))
candidate_rank_stability <- candidate_rank_stability[match(
  candidates$gene_symbol,
  candidate_rank_stability$gene_symbol
), , drop = FALSE]

write_tsv <- function(frame, filename) {
  write.table(
    frame,
    file.path(output_root, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}
write_tsv(
  permutation_summary,
  "cross_model_permutation_summary.tsv"
)
write_tsv(
  data.frame(
    iteration = seq_len(permutations),
    complete_pattern_genes = permutation_counts,
    raw_processed_robust_overlap = robust_overlap_null
  ),
  "cross_model_permutation_null.tsv"
)
write_tsv(negative_controls, "cross_model_negative_controls.tsv")
write_tsv(score_summary, "cross_model_score_sensitivity_summary.tsv")
write_tsv(candidate_rank_long, "cross_model_candidate_rank_variants.tsv")
write_tsv(
  candidate_rank_stability,
  "cross_model_candidate_rank_stability.tsv"
)

plot_permutation <- function(device) {
  device()
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
  hist(
    permutation_counts,
    breaks = 35,
    col = "#8CBEB2",
    border = "white",
    xlab = "Complete-pattern genes under permutation",
    main = ""
  )
  abline(v = observed_candidate_count, col = "#B23A48", lwd = 2)
  hist(
    robust_overlap_null,
    breaks = seq(-0.5, 7.5, by = 1),
    col = "#E6B655",
    border = "white",
    xlab = "Raw and processed robust overlap",
    main = ""
  )
  abline(v = observed_robust_overlap, col = "#B23A48", lwd = 2)
}
plot_permutation(function() {
  grDevices::pdf(
    file.path(output_root, "cross_model_permutation_null.pdf"),
    width = 8.5,
    height = 4.25,
    useDingbats = FALSE
  )
})
plot_permutation(function() {
  grDevices::png(
    file.path(output_root, "cross_model_permutation_null.png"),
    width = 5100,
    height = 2550,
    res = 600
  )
})

mandatory <- c(
  "cross_model_permutation_summary.tsv",
  "cross_model_permutation_null.tsv",
  "cross_model_negative_controls.tsv",
  "cross_model_score_sensitivity_summary.tsv",
  "cross_model_candidate_rank_variants.tsv",
  "cross_model_candidate_rank_stability.tsv",
  "cross_model_permutation_null.pdf",
  "cross_model_permutation_null.png"
)
stopifnot(
  all(file.exists(file.path(output_root, mandatory))),
  all(file.info(file.path(output_root, mandatory))$size > 0L)
)
write_tsv(
  data.frame(
    status = "COMPLETE",
    completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    permutations = permutations,
    seed = seed,
    frozen_candidate_count = observed_candidate_count,
    candidate_empirical_p = candidate_empirical_p,
    robust_overlap = observed_robust_overlap,
    robust_overlap_empirical_p = robust_overlap_empirical_p,
    stringsAsFactors = FALSE
  ),
  "CROSS_MODEL_PERMUTATION_SENSITIVITY_COMPLETE.tsv"
)

cat(
  "Permutation and score sensitivity complete:",
  "candidate p =", format(candidate_empirical_p, digits = 4),
  "; robust-overlap p =", format(robust_overlap_empirical_p, digits = 4),
  "\n"
)
