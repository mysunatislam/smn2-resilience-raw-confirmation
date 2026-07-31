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

read_frozen <- function(name) {
  path <- file.path(ROOT, "config", name)
  if (!file.exists(path)) {
    stop("Missing frozen file: ", name)
  }
  read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

genes <- read_frozen("frozen_37_gene_shortlist_2026-08-01.tsv")
events <- read_frozen("frozen_83_splice_events_2026-08-01.tsv")
panel <- read_frozen("frozen_12_primary_splice_events_2026-08-01.tsv")

expected_genes <- c(
  "LY6H", "HS3ST5", "ZNF853", "IL17D", "LSAMP", "FAAH", "C12orf60",
  "KCNJ4", "CABLES2", "PNCK", "CLPTM1", "PDPR", "HEY1", "NRN1",
  "KCNIP2", "ROBO2", "UGT8", "CTBP1-DT", "RAB28", "SHQ1", "SPIN1",
  "SLC38A9", "ADCYAP1R1", "HLA-DPB1", "ZMAT4", "LINC00665", "FCSK",
  "LRRC4", "PFKFB2", "ADHFE1", "ZNF615", "TTTY15", "ATP9B",
  "ZNF204P", "EIF4ENIF1", "ZNF827", "RIOX1"
)
if (!identical(genes$gene_symbol, expected_genes)) {
  stop("Frozen 37-gene order changed")
}

expected_tier_1 <- c("LY6H", "HS3ST5", "ZNF853", "IL17D")
observed_tier_1 <- genes$gene_symbol[
  genes$publication_evidence_tier ==
    "tier_1_raw_and_processed_unit_robust"
]
if (!identical(observed_tier_1, expected_tier_1)) {
  stop("Frozen Tier 1 genes changed")
}

if (nrow(events) != 83L || anyDuplicated(events$event_key)) {
  stop("Frozen 83-event manifest changed")
}
expected_event_counts <- c(A3SS = 5L, A5SS = 6L, ES = 35L, MXE = 22L, RI = 15L)
observed_event_counts <- table(events$event_type)
if (!identical(
  as.integer(observed_event_counts[names(expected_event_counts)]),
  as.integer(expected_event_counts)
)) {
  stop("Frozen event-class counts changed")
}

expected_panel_genes <- c(
  "MRFAP1", "IFI27L1", "PDE9A", "HERPUD1", "ABCA1", "SAT2",
  "TNPO2", "SEPTIN11", "FAT3", "MGEA5", "COL5A2", "EIF5"
)
expected_panel_types <- c(
  "RI", "ES", "ES", "A5SS", "RI", "ES",
  "A3SS", "MXE", "MXE", "MXE", "RI", "A5SS"
)
if (
  nrow(panel) != 12L ||
  !identical(panel$source_gene_symbol, expected_panel_genes) ||
  !identical(panel$event_type, expected_panel_types) ||
  !identical(panel$panel_order, seq_len(12L))
) {
  stop("Frozen 12-event primary panel changed")
}
if (!all(panel$event_key %in% events$event_key)) {
  stop("Primary panel contains an event outside the frozen 83-event set")
}

declaration <- file.path(ROOT, "docs", "validation_freeze_2026-08-01.md")
checksums <- file.path(ROOT, "config", "validation_freeze_2026-08-01.sha256")
if (!file.exists(declaration) || !file.exists(checksums)) {
  stop("Freeze declaration or checksum list is missing")
}

cat("Validation freeze test passed\n")
