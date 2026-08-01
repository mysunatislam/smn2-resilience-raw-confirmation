# Raw STAR Benchmark Decision

Date: `2026-08-01`

Status: `LOCKED_BEFORE_RMATS_RESULTS`

## Observations

The first `BULK-SAM-148` attempt stopped before mapping because STAR requires
FIFO support and the E-drive NTFS mount does not provide it. The failed logs
were preserved under `star_rmats/failed_runs/20260801_fifo_ntfs_BULK-SAM-148`.
The reproducible fix places only STAR's ephemeral FIFO directory on the WSL
Linux filesystem. Permanent outputs remain on E:.

The corrected two-pass attempt loaded the complete sparse GENCODE v47 index
and reached 38,732,418 of 65,414,230 read pairs in first-pass mapping before
the command was externally interrupted. STAR reported approximately 7.3
million read pairs per hour and 71.2% unique mapping at that point. No rMATS
table or event-level result existed.

At that measured rate, one full first pass requires about nine hours for this
library. Two-pass mapping would require approximately 18 hours per library
and more than six days for the nine-library cohort.

## Locked Decision

The full-depth cohort uses annotation-guided one-pass STAR:

```text
twopassMode=None
readMapNumber=-1
```

Every read remains eligible for alignment. No FASTQ is downsampled. The STAR
index already includes splice junctions from the same frozen GENCODE v47 GTF,
and all 83 confirmatory targets are annotated fixed-set events. One-pass
alignment therefore answers the prospective annotated-junction confirmation
question without the runtime cost of de novo two-pass junction augmentation.

The sparse-index parameters, genome, annotation, sample set, target loci,
rMATS settings, event-matching rules, and success thresholds remain
unchanged.
