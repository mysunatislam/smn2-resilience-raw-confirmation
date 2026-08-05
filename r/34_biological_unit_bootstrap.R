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

data_root <- argument_value("--data-root=")
integration_path <- argument_value("--integration-table=")
candidate_path <- argument_value(
  "--candidate-table=",
  file.path(ROOT, "manuscript", "supplementary_table_S1_candidates.tsv")
)
output_root <- argument_value("--output-root=")
iterations <- as.integer(argument_value("--iterations=", "500"))
seed <- as.integer(argument_value("--seed=", "939390"))
if (
  is.null(data_root) ||
    is.null(integration_path) ||
    is.null(output_root) ||
    !is.finite(iterations) ||
    iterations < 200L
) {
  stop(
    "Required: --data-root=PATH --integration-table=PATH ",
    "--output-root=PATH [--iterations=500]"
  )
}

suppressPackageStartupMessages(library(edgeR))

data_root <- normalizePath(data_root, winslash = "/", mustWork = TRUE)
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

read_counts <- function(accession) {
  path <- file.path(data_root, paste0(accession, "_counts.csv.gz"))
  frame <- read.csv(
    gzfile(path),
    row.names = 1L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  matrix <- as.matrix(frame)
  storage.mode(matrix) <- "integer"
  matrix
}
read_metadata <- function(accession) {
  path <- file.path(data_root, paste0(accession, "_metadata.tsv"))
  read.delim(
    path,
    row.names = 1L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
aggregate_counts <- function(counts, groups) {
  t(rowsum(t(counts), group = groups, reorder = FALSE))
}
normalized_log_cpm <- function(counts) {
  y <- DGEList(counts = counts)
  y <- normLibSizes(y, method = "TMM")
  cpm(y, log = TRUE, prior.count = 0.5)
}
finite <- function(values) !is.na(values) & is.finite(values)
as_flag <- function(values) !is.na(values) & as.logical(values)
zscore <- function(values) {
  result <- as.numeric(scale(values))
  result[!is.finite(result)] <- 0
  result
}
percentile_rank <- function(values, use) {
  output <- rep(NA_real_, length(values))
  eligible <- use & finite(values)
  output[eligible] <- (
    rank(values[eligible], ties.method = "average") - 0.5
  ) / sum(eligible)
  output
}

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
stopifnot(
  !anyDuplicated(integration$gene_symbol),
  !anyDuplicated(candidates$gene_symbol),
  nrow(candidates) == 37L
)
integration <- integration[
  as_flag(integration$cross_model_ranking_eligible),
  ,
  drop = FALSE
]
stopifnot(nrow(integration) == 11326L)

gse93939_counts <- read_counts("GSE93939")
gse93939_metadata <- read_metadata("GSE93939")
gse93939_metadata <- gse93939_metadata[
  gse93939_metadata$group %in% c("OMN", "SC"),
  ,
  drop = FALSE
]
stopifnot(
  all(rownames(gse93939_metadata) %in% colnames(gse93939_counts)),
  nrow(gse93939_metadata) == 32L,
  length(unique(gse93939_metadata$donor_id)) == 19L
)
gse93939_counts <- gse93939_counts[
  ,
  rownames(gse93939_metadata),
  drop = FALSE
]
gse93939_log_cpm <- normalized_log_cpm(gse93939_counts)
natural_gene_index <- match(
  integration$gene_symbol,
  rownames(gse93939_log_cpm)
)
natural_expression <- matrix(
  NA_real_,
  nrow = ncol(gse93939_log_cpm),
  ncol = nrow(integration)
)
natural_present <- !is.na(natural_gene_index)
natural_expression[, natural_present] <- t(
  gse93939_log_cpm[natural_gene_index[natural_present], , drop = FALSE]
)

natural_design <- cbind(
  `(Intercept)` = 1,
  OMN_vs_SC = as.integer(gse93939_metadata$group == "OMN"),
  sex_M = as.integer(gse93939_metadata$sex == "M"),
  age_z = zscore(gse93939_metadata$age_at_death),
  postmortem_delay_z = zscore(
    gse93939_metadata$postmortem_delay_hours
  ),
  source_NDRI = as.integer(gse93939_metadata$tissue_source == "NDRI"),
  source_NIH = as.integer(gse93939_metadata$tissue_source == "NIH"),
  platform_HiSeq2500 = as.integer(
    gse93939_metadata$platform == "HiSeq2500"
  )
)
variable_design_columns <- vapply(
  seq_len(ncol(natural_design)),
  function(index) {
    colnames(natural_design)[index] == "(Intercept)" ||
      stats::sd(natural_design[, index]) > 0
  },
  logical(1)
)
natural_design <- natural_design[
  ,
  variable_design_columns,
  drop = FALSE
]
stopifnot(
  qr(natural_design)$rank == ncol(natural_design),
  "OMN_vs_SC" %in% colnames(natural_design)
)
donors <- unique(gse93939_metadata$donor_id)
sample_donor_index <- match(gse93939_metadata$donor_id, donors)

gse290979_counts <- read_counts("GSE290979")
gse290979_metadata <- read_metadata("GSE290979")
gse290979_metadata <- gse290979_metadata[
  colnames(gse290979_counts),
  ,
  drop = FALSE
]

untreated_metadata <- gse290979_metadata[
  gse290979_metadata$treatment == "NT",
  ,
  drop = FALSE
]
disease_counts <- aggregate_counts(
  gse290979_counts[, rownames(untreated_metadata), drop = FALSE],
  untreated_metadata$`cell line`
)
disease_log_cpm <- normalized_log_cpm(disease_counts)
disease_line_genotype <- vapply(
  colnames(disease_counts),
  function(cell_line) {
    unique(untreated_metadata$genotype[
      untreated_metadata$`cell line` == cell_line
    ])
  },
  character(1)
)
control_lines <- which(disease_line_genotype == "CTRL")
sma_lines <- which(disease_line_genotype == "SMA")
stopifnot(length(control_lines) == 3L, length(sma_lines) == 2L)

treatment_metadata <- gse290979_metadata[
  gse290979_metadata$treatment %in% c("Scramble", "R6-Mo"),
  ,
  drop = FALSE
]
treatment_groups <- paste(
  treatment_metadata$`cell line`,
  treatment_metadata$treatment,
  sep = "__"
)
treatment_counts <- aggregate_counts(
  gse290979_counts[, rownames(treatment_metadata), drop = FALSE],
  treatment_groups
)
treatment_log_cpm <- normalized_log_cpm(treatment_counts)
treatment_lines <- sort(unique(treatment_metadata$`cell line`))
stopifnot(length(treatment_lines) == 2L)
treatment_line_effect <- vapply(
  treatment_lines,
  function(cell_line) {
    treatment_log_cpm[, paste0(cell_line, "__R6-Mo")] -
      treatment_log_cpm[, paste0(cell_line, "__Scramble")]
  },
  numeric(nrow(treatment_log_cpm))
)
rownames(treatment_line_effect) <- rownames(treatment_log_cpm)

disease_gene_index <- match(
  integration$gene_symbol,
  rownames(disease_log_cpm)
)
treatment_gene_index <- match(
  integration$gene_symbol,
  rownames(treatment_line_effect)
)
external_effect <- integration$gse108094_sma_vs_control_log2_effect
natural_support <- as_flag(
  integration$exploratory_natural_resistance_support
)
candidate_index <- match(candidates$gene_symbol, integration$gene_symbol)
stopifnot(!anyNA(candidate_index))

candidate_count <- nrow(candidates)
selection <- matrix(FALSE, nrow = iterations, ncol = candidate_count)
natural_direction <- matrix(FALSE, nrow = iterations, ncol = candidate_count)
disease_direction <- matrix(FALSE, nrow = iterations, ncol = candidate_count)
treatment_direction <- matrix(FALSE, nrow = iterations, ncol = candidate_count)
candidate_rank <- matrix(
  NA_real_,
  nrow = iterations,
  ncol = candidate_count
)
selected_universe_count <- integer(iterations)
bootstrap_eligible_count <- integer(iterations)
natural_resamples_attempted <- integer(iterations)

set.seed(seed)
for (iteration in seq_len(iterations)) {
  natural_attempts <- 0L
  repeat {
    natural_attempts <- natural_attempts + 1L
    sampled_donors <- sample(
      seq_along(donors),
      length(donors),
      replace = TRUE
    )
    donor_multiplicity <- tabulate(
      sampled_donors,
      nbins = length(donors)
    )
    sample_weights <- donor_multiplicity[sample_donor_index]
    weighted_design <- natural_design * sqrt(sample_weights)
    if (qr(weighted_design)$rank == ncol(natural_design)) {
      break
    }
    if (natural_attempts >= 100L) {
      stop("Could not obtain a full-rank donor bootstrap design")
    }
  }
  weighted_expression <- natural_expression * sqrt(sample_weights)
  natural_coefficients <- qr.coef(
    qr(weighted_design),
    weighted_expression
  )
  natural_effect <- natural_coefficients[
    match("OMN_vs_SC", rownames(natural_coefficients)),
    ,
    drop = TRUE
  ]

  sampled_controls <- sample(
    control_lines,
    length(control_lines),
    replace = TRUE
  )
  sampled_sma <- sample(sma_lines, length(sma_lines), replace = TRUE)
  disease_effect_source <- rowMeans(
    disease_log_cpm[, sampled_sma, drop = FALSE]
  ) - rowMeans(
    disease_log_cpm[, sampled_controls, drop = FALSE]
  )
  disease_effect <- rep(NA_real_, nrow(integration))
  disease_present <- !is.na(disease_gene_index)
  disease_effect[disease_present] <- disease_effect_source[
    disease_gene_index[disease_present]
  ]

  sampled_treatment_lines <- sample(
    seq_along(treatment_lines),
    length(treatment_lines),
    replace = TRUE
  )
  treatment_effect_source <- rowMeans(
    treatment_line_effect[
      ,
      sampled_treatment_lines,
      drop = FALSE
    ]
  )
  treatment_effect <- rep(NA_real_, nrow(integration))
  treatment_present <- !is.na(treatment_gene_index)
  treatment_effect[treatment_present] <- treatment_effect_source[
    treatment_gene_index[treatment_present]
  ]

  bootstrap_eligible <- finite(natural_effect) &
    finite(disease_effect) &
    finite(treatment_effect) &
    finite(external_effect)
  bootstrap_eligible_count[[iteration]] <- sum(bootstrap_eligible)
  selected <- bootstrap_eligible &
    natural_support &
    natural_effect > 0 &
    disease_effect < 0 &
    treatment_effect > 0 &
    external_effect < 0
  selected_universe_count[[iteration]] <- sum(selected)

  component_matrix <- cbind(
    percentile_rank(natural_effect, bootstrap_eligible),
    percentile_rank(-disease_effect, bootstrap_eligible),
    percentile_rank(treatment_effect, bootstrap_eligible),
    percentile_rank(-external_effect, bootstrap_eligible)
  )
  score <- rowMeans(component_matrix, na.rm = FALSE)
  ranks <- rep(NA_real_, nrow(integration))
  ranks[bootstrap_eligible] <- rank(
    -score[bootstrap_eligible],
    ties.method = "average"
  )

  selection[iteration, ] <- selected[candidate_index]
  natural_direction[iteration, ] <-
    finite(natural_effect[candidate_index]) &
    natural_effect[candidate_index] > 0
  disease_direction[iteration, ] <-
    finite(disease_effect[candidate_index]) &
    disease_effect[candidate_index] < 0
  treatment_direction[iteration, ] <-
    finite(treatment_effect[candidate_index]) &
    treatment_effect[candidate_index] > 0
  candidate_rank[iteration, ] <- ranks[candidate_index]
  natural_resamples_attempted[[iteration]] <- natural_attempts
}

rank_quantile <- function(values, probability) {
  values <- values[finite(values)]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  unname(quantile(values, probability, type = 8))
}
candidate_stability <- data.frame(
  gene_symbol = candidates$gene_symbol,
  selection_frequency = colMeans(selection),
  natural_positive_frequency = colMeans(natural_direction),
  disease_negative_frequency = colMeans(disease_direction),
  treatment_positive_frequency = colMeans(treatment_direction),
  external_negative_fixed = external_effect[candidate_index] < 0,
  median_rank = apply(candidate_rank, 2L, median, na.rm = TRUE),
  rank_q025 = apply(candidate_rank, 2L, rank_quantile, probability = 0.025),
  rank_q975 = apply(candidate_rank, 2L, rank_quantile, probability = 0.975),
  rank_estimable_iterations = colSums(is.finite(candidate_rank)),
  iterations = iterations,
  stringsAsFactors = FALSE
)
candidate_stability$median_rank[
  !is.finite(candidate_stability$median_rank)
] <- NA_real_

iteration_summary <- data.frame(
  iteration = seq_len(iterations),
  bootstrap_eligible_genes = bootstrap_eligible_count,
  selected_complete_pattern_genes = selected_universe_count,
  frozen_37_selected = rowSums(selection),
  natural_design_attempts = natural_resamples_attempted,
  stringsAsFactors = FALSE
)
summary_table <- data.frame(
  metric = c(
    "iterations",
    "seed",
    "natural_resampled_donors",
    "disease_resampled_control_lines",
    "disease_resampled_sma_lines",
    "treatment_resampled_pairs",
    "external_effect_resampled",
    "median_bootstrap_eligible_genes",
    "median_complete_pattern_genes",
    "complete_pattern_genes_q025",
    "complete_pattern_genes_q975",
    "median_frozen_37_selected",
    "frozen_37_selected_q025",
    "frozen_37_selected_q975"
  ),
  value = c(
    iterations,
    seed,
    length(donors),
    length(control_lines),
    length(sma_lines),
    length(treatment_lines),
    FALSE,
    median(bootstrap_eligible_count),
    median(selected_universe_count),
    unname(quantile(selected_universe_count, 0.025)),
    unname(quantile(selected_universe_count, 0.975)),
    median(rowSums(selection)),
    unname(quantile(rowSums(selection), 0.025)),
    unname(quantile(rowSums(selection), 0.975))
  ),
  stringsAsFactors = FALSE
)
provenance <- data.frame(
  metric = c(
    "analysis_role",
    "natural_model",
    "disease_model",
    "treatment_model",
    "external_model",
    "resampling_unit",
    "random_library_split"
  ),
  value = c(
    "biological_unit_cluster_bootstrap_sensitivity",
    "covariate_adjusted_TMM_logCPM_weighted_OLS",
    "TMM_logCPM_donor_line_mean_difference",
    "TMM_logCPM_paired_line_difference",
    "GSE108094_effect_fixed_not_unit_resampled",
    "donor_or_donor_line",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

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
  candidate_stability,
  "biological_unit_bootstrap_candidate_stability.tsv"
)
write_tsv(
  iteration_summary,
  "biological_unit_bootstrap_iterations.tsv"
)
write_tsv(summary_table, "biological_unit_bootstrap_summary.tsv")
write_tsv(provenance, "biological_unit_bootstrap_provenance.tsv")

plot_stability <- function(device) {
  device()
  on.exit(grDevices::dev.off(), add = TRUE)
  order_index <- order(
    candidate_stability$selection_frequency,
    decreasing = FALSE
  )
  par(mar = c(4.5, 8.5, 2.5, 1))
  plot(
    candidate_stability$selection_frequency[order_index],
    seq_along(order_index),
    xlim = c(0, 1),
    ylim = c(0.5, length(order_index) + 0.5),
    pch = 19,
    col = "#176B87",
    axes = FALSE,
    xlab = "Biological-unit bootstrap selection frequency",
    ylab = "",
    main = "Frozen candidate stability"
  )
  axis(1)
  axis(
    2,
    at = seq_along(order_index),
    labels = candidate_stability$gene_symbol[order_index],
    las = 1,
    cex.axis = 0.65
  )
  abline(v = c(0.5, 0.8), lty = c(3, 2), col = "grey55")
  box()
}
plot_stability(function() {
  grDevices::pdf(
    file.path(output_root, "biological_unit_bootstrap_stability.pdf"),
    width = 7.5,
    height = 9,
    useDingbats = FALSE
  )
})
plot_stability(function() {
  grDevices::png(
    file.path(output_root, "biological_unit_bootstrap_stability.png"),
    width = 1500,
    height = 1800,
    res = 200
  )
})

mandatory <- c(
  "biological_unit_bootstrap_candidate_stability.tsv",
  "biological_unit_bootstrap_iterations.tsv",
  "biological_unit_bootstrap_summary.tsv",
  "biological_unit_bootstrap_provenance.tsv",
  "biological_unit_bootstrap_stability.pdf",
  "biological_unit_bootstrap_stability.png"
)
stopifnot(
  all(file.exists(file.path(output_root, mandatory))),
  all(file.info(file.path(output_root, mandatory))$size > 0L),
  nrow(candidate_stability) == 37L,
  nrow(iteration_summary) == iterations
)
write_tsv(
  data.frame(
    status = "COMPLETE",
    completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    iterations = iterations,
    seed = seed,
    random_library_split = FALSE,
    stringsAsFactors = FALSE
  ),
  "BIOLOGICAL_UNIT_BOOTSTRAP_COMPLETE.tsv"
)

cat(
  "Biological-unit bootstrap complete:",
  iterations, "iterations;",
  "median frozen candidates selected =",
  median(rowSums(selection)), "\n"
)
