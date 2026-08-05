# Cross-model human motor-neuron resilience in SMA

Reviewer-facing code and outputs for a donor-aware integration of natural
motor-neuron resistance, human SMA disease, SMN-directed response,
cell-resolved context, and raw-read confirmation.

## Verified result

- Core datasets: GSE93939 (39 libraries; 19 donors in the primary contrast)
  and GSE290979 (31 libraries; 5 donor lines).
- External direction: GSE108094 and GSE69175.
- Cell-resolved context: GSE290980 and GSE243076.
- Frozen cross-model candidates: 37.
- Raw confirmation: 9 pre-specified GSE290979 libraries, 18 deposited
  checksum-verified FASTQ files, Salmon 1.10.2, GENCODE v47.
- Raw-versus-processed concordance:
  - SMA versus control: 14,636 genes, Spearman rho 0.778, direction agreement
    0.807.
  - R6-MO versus scramble: 14,385 genes, Spearman rho 0.588, direction
    agreement 0.732.
- Raw full-direction candidates: 21.
- Raw biological-unit robust candidates: 9.
- Robust in both raw and processed analyses: `LY6H`, `HS3ST5`, `ZNF853`,
  `IL17D`.
- Raw splice-junction confirmation: 6/83 exact structures recovered in both
  contrasts; 1/12 primary events recovered; 0 events met the full two-line
  correction criterion.

The original discovery ranks are unchanged. Raw confirmation is same-study
sensitivity evidence, not independent validation.

## Repository contents

```text
config/       Locked dataset roles, sample sheets, contrasts, and resampling rules
data/metadata GEO records and the locked ENA run manifest
docs/         Dataset rationale, statistical methods, and protocol details
local/        Windows/WSL all-nine Salmon and STAR-rMATS raw-read workflows
r/            Authoritative R analysis
scripts/      Source-data preparation
tests/        Deterministic regression tests
manuscript/   Manuscript, figures, and supplementary tables
```

Large public inputs and generated local results are excluded. GEO accessions,
ENA runs, checksums, sample selections, and expected paths are recorded in
`config/`, `data/metadata/`, and `docs/`.

## Main outputs

- [`manuscript/manuscript.md`](manuscript/manuscript.md)
- [`manuscript/supplementary_table_S1_candidates.tsv`](manuscript/supplementary_table_S1_candidates.tsv)
- [`manuscript/supplementary_table_S2_summary.tsv`](manuscript/supplementary_table_S2_summary.tsv)
- [`manuscript/supplementary_table_S3_concordance.tsv`](manuscript/supplementary_table_S3_concordance.tsv)
- [`manuscript/supplementary_figure_index.tsv`](manuscript/supplementary_figure_index.tsv)
- [`docs/validation_freeze_2026-08-01.md`](docs/validation_freeze_2026-08-01.md)
- [`manuscript/figures/figure_1_raw_vs_processed_concordance.pdf`](manuscript/figures/figure_1_raw_vs_processed_concordance.pdf)
- [`manuscript/figures/figure_2_candidate_evidence_matrix.pdf`](manuscript/figures/figure_2_candidate_evidence_matrix.pdf)
- [`manuscript/figures/figure_3_raw_splice_confirmation.pdf`](manuscript/figures/figure_3_raw_splice_confirmation.pdf)
- [`manuscript/figures/`](manuscript/figures/)

## Validation

With R 4.6.0 and the project packages installed:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_local9_profile.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_local9_count_analysis.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_local9_raw_concordance.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_local9_salmon_import.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_publication_integration.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_supplementary_figure_package.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_validation_freeze.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_local9_star_rmats_workflow.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_local9_star_rmats_analysis.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' tests\test_raw_splice_publication_package.R
```

The completed gene-level analysis passed its seven regression tests and 66
full-project output checks. The raw splice-junction workflow and tracked
publication package add three focused regression tests.

## Scope

The raw result includes all-nine transcriptome re-quantification and targeted
STAR-rMATS confirmation at the frozen splice loci. It is not independent
cohort validation, complete 31-library splice discovery, or paralog-specific
`SMN1`/`SMN2` quantification. Candidate genes and splice events are
experimental hypotheses, not therapeutic recommendations.
