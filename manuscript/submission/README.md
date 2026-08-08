# Submission Readiness

This directory separates completed evidence from work that cannot be inferred
or fabricated.

- `numerical_consistency_audit.tsv`: executable audit of the frozen numerical
  claims.
- `manuscript_structure_audit.tsv`: abstract, reference, DOI, and citation-order
  checks.
- `experimental_validation_status.tsv`: explicit status of functional and
  RT-PCR validation.
- `figure_asset_audit.tsv`: dimensions, resolution, and vector-PDF coverage
  for all 28 main and supplementary figures.
- `metadata_required.tsv`: author-supplied information still required.
- `submission_readiness.tsv`: current pass and blocker summary.
- `cover_letter_draft.md`: scientific cover-letter body, not a signed final
  letter.
- `release_and_archive_checklist.md`: final GitHub and Zenodo freeze sequence.
- `supplementary_figure_captions.tsv`: frozen captions used by the PDF
  packager.

`scripts/build_submission_supplement.py` creates the review-only PDF under
`output/pdf/` and the machine-readable table archive under
`output/submission/`. Final mode fails while required metadata is incomplete.

The package is not ready for submission while experimental claims and author
metadata remain unresolved. Synthetic fixtures and protocol templates are not
experimental evidence.
