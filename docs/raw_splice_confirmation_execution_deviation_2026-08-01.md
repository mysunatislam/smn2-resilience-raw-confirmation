# Raw Splice Confirmation Execution Note

> Implementation update (2026-08-05): fixed-event rMATS inference returned
> empty MATS tables despite direct STAR junction support. The frozen target
> set and success criteria were unchanged, but inference was amended to
> target-locus `--novelSS` discovery followed by exact frozen-event matching.
> See `raw_splice_confirmation_targeted_discovery_amendment_2026-08-05.md`.

Date: `2026-08-01`

Status: `RECORDED_BEFORE_STAR_RMATS_RESULTS`

The frozen targets and success criteria in
`docs/validation_freeze_2026-08-01.md` remain unchanged. This note records one
compute-limited execution detail before any STAR-rMATS validation result is
examined.

## Targeted BAM Retention

STAR aligns every read against the complete GENCODE v47 GRCh38 primary
assembly. Its SAM stream is filtered to the 83 frozen event loci with 1 kb
padding and stored as coordinate-sorted BAM. STAR alignment summaries and the
genome-wide `SJ.out.tab` junction table are retained. A full-genome BAM is not
created because the local machine cannot preserve nine full BAMs without
deleting verified FASTQs or existing BAMs.

The retained BAMs therefore support prospective testing of the frozen 83
events, but they are not suitable for de novo genome-wide event discovery.

## FDR Scope

rMATS is run with the frozen 83-event fixed set. Its reported `FDR` is a
fixed-set quantity and must not be described as genome-wide. The confirmatory
multiplicity result is a separately computed Benjamini-Hochberg q-value across
the 83 frozen events, as specified in the freeze.

The originally requested genome-wide rMATS FDR is recorded as unavailable in
this local targeted execution. It can only be added after a full-BAM
genome-wide run on storage-capable infrastructure. Its absence does not change
the frozen panel denominator, direction thresholds, junction-support
threshold, line-specific correction criteria, or primary success decision.

## Sashimi Plots

The 12 frozen primary events are plotted with Xinglab
`rmats2sashimiplot` v4.0.0. Each plot contains grouped tracks for control,
untreated SMA, scramble, and R6-MO libraries. Junction-count differences
between rMATS and the plotting backend are treated as visualization-method
differences; inferential counts come only from rMATS.
