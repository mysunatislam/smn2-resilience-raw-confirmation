# Publication Package

This directory contains the tracked, submission-facing outputs for the
cross-model human motor-neuron resilience analysis.

## Contents

- `manuscript.md`: complete working manuscript draft.
- `figures/figure_1_raw_vs_processed_concordance.*`: genome-wide raw versus
  deposited processed-matrix effect concordance.
- `figures/figure_2_candidate_evidence_matrix.*`: evidence map for the frozen
  37-candidate set and publication tiers.
- `figures/figure_3_raw_splice_confirmation.*`: exact raw-versus-processed
  delta-PSI comparison for the six events recovered in both contrasts.
- `figures/supplementary_figure_S1_*.png` through
  `figures/sashimi/supplementary_figure_S23_*.png`: QC, cross-model
  integration, holdout robustness, cell-resolved context, splicing priority,
  and the 12-event raw sashimi panel.
- `figures/supplementary_figure_S24_*` and
  `figures/supplementary_figure_S25_*`: permutation nulls and
  biological-unit bootstrap stability.
- `supplementary_figure_index.tsv`: figure titles, evidence roles, and
  generating R scripts.
- `supplementary_table_index.tsv`: filenames and titles for Supplementary
  Tables S1-S15.
- `supplementary_table_S1_candidates.tsv`: candidate-level effects,
  biological-unit checks, raw confirmation, and external annotations.
- `supplementary_table_S2_summary.tsv`: locked analysis scope and tier counts.
- `supplementary_table_S3_concordance.tsv`: raw-versus-processed concordance
  statistics.
- `supplementary_table_S4_*.tsv` through `supplementary_table_S9_*.tsv`:
  targeted raw splice results, mismatch reasons, support, summary, and
  provenance.
- `supplementary_table_S10_*.tsv` through `supplementary_table_S15_*.tsv`:
  permutation, negative-control, score-sensitivity, and biological-unit
  bootstrap results.
- `audit/`: iteration-level permutation and bootstrap records, provenance,
  and completion markers.
- `figures/sashimi/`: PDF and PNG raw-junction plots for all 12 frozen primary
  events.

The three main figures are supplied as PNG and vector PDF. Supplementary figures
are supplied as publication-resolution PNG files. The tracked figures and
tables are generated from validated local results and should not be edited
manually.

## Regenerate

From the repository root:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\28_integrate_local9_salmon_publication.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\06_validate_outputs.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\33_cross_model_permutation_sensitivity.R --integration-table=PATH --output-root=PATH
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\34_biological_unit_bootstrap.R --data-root=PATH --integration-table=PATH --output-root=PATH
& 'C:\path\to\python.exe' scripts\audit_figure_assets.py
& 'C:\path\to\python.exe' scripts\build_submission_supplement.py --mode review
```

The first command requires the completed local9 Salmon outputs under
`results/r/raw_confirmation/local9_salmon/`. The second command verifies the
full project and locks the publication counts, candidate identities, and
figure presence.

## Claim Boundary

The completed raw analysis includes all-nine, checksum-verified Salmon
transcriptome re-quantification and targeted STAR-rMATS junction confirmation
of the frozen splice loci. It is not independent cohort validation, complete
31-library splice discovery, or SMN1-versus-SMN2 allele-specific
quantification. The 37-gene count is not permutation-enriched; the stronger
computational result is reproducibility of the four-gene raw/processed
unit-robust overlap and its biological-unit bootstrap stability.
