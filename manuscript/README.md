# Publication Package

This directory contains the tracked, submission-facing outputs for the
cross-model human motor-neuron resilience analysis.

## Contents

- `manuscript.md`: complete working manuscript draft.
- `figures/figure_1_raw_vs_processed_concordance.*`: genome-wide raw versus
  deposited processed-matrix effect concordance.
- `figures/figure_2_candidate_evidence_matrix.*`: evidence map for the frozen
  37-candidate set and publication tiers.
- `figures/supplementary_figure_S1_*.png` through
  `figures/supplementary_figure_S11_*.png`: QC, cross-model integration,
  holdout robustness, cell-resolved context, and splicing-priority figures.
- `supplementary_figure_index.tsv`: figure titles, evidence roles, and
  generating R scripts.
- `supplementary_table_S1_candidates.tsv`: candidate-level effects,
  biological-unit checks, raw confirmation, and external annotations.
- `supplementary_table_S2_summary.tsv`: locked analysis scope and tier counts.
- `supplementary_table_S3_concordance.tsv`: raw-versus-processed concordance
  statistics.

The two main figures are supplied as PNG and vector PDF. Supplementary figures
are supplied as publication-resolution PNG files. The tracked figures and
tables are generated from validated local results and should not be edited
manually.

## Regenerate

From the repository root:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\28_integrate_local9_salmon_publication.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\06_validate_outputs.R
```

The first command requires the completed local9 Salmon outputs under
`results/r/raw_confirmation/local9_salmon/`. The second command verifies the
full project and locks the publication counts, candidate identities, and
figure presence.

## Claim Boundary

The completed raw analysis is an all-nine, checksum-verified, Salmon
transcriptome re-quantification of GSE290979. It is not independent cohort
validation, whole-genome alignment, splice-junction discovery, or
SMN1-versus-SMN2 allele-specific quantification.
