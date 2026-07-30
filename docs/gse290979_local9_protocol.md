# GSE290979 Local9 Workstation Protocol

**Purpose:** obtain limited workbook-independent raw-read evidence without HPC,
paid cloud compute, random library splitting, or partial FASTQ truncation.

**Execution status, 30 July 2026:** all nine complete FASTQ pairs passed
deposited ENA MD5 verification and completed the expedited Salmon workflow.
The 13-check integrity gate passed and wrote
`LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv`. The whole-genome completion marker
remains absent by design.

## Scientific role

The local9 profile preserves every independent biological unit available for
the two core GSE290979 questions:

| Contrast | Complete libraries | Independent units |
| --- | ---: | ---: |
| Untreated SMA versus control | 5 | 2 SMA and 3 control donor lines |
| R6-MO versus scramble | 4 | 2 paired SMA donor lines |

One technical preparation is selected per donor-line condition by minimum
distance from that group's median read-pair count, followed by run accession
as a deterministic tie-breaker. The nine selected libraries total 82.193 GiB.
No outcome, gene, splice event, random seed, or processed-matrix result enters
selection.

The full 31-library processed matrix remains the genome-wide discovery source.
Local9 is a sensitivity analysis of raw expression directions and
pre-specified candidate splice evidence. It cannot create either
`RAW_VALIDATION21_COMPLETE` or `RAW_REANALYSIS_COMPLETE`.

## Frozen libraries

| Role | Donor line | Treatment | Sample | Run |
| --- | --- | --- | --- | --- |
| Disease | C1 | NT | BULK-SAM-109 | SRR32559040 |
| Disease | C2 | NT | BULK-SAM-104 | SRR32559042 |
| Disease | C3 | NT | BULK-SAM-107 | SRR32559038 |
| Disease | S2 | NT | BULK-SAM-160 | SRR32559021 |
| Disease and benchmark | S3 | NT | BULK-SAM-148 | SRR32559034 |
| Treatment | S2 | Scramble | BULK-SAM-163 | SRR32559019 |
| Treatment | S3 | Scramble | BULK-SAM-151 | SRR32559030 |
| Treatment | S2 | R6-MO | BULK-SAM-166 | SRR32559016 |
| Treatment | S3 | R6-MO | BULK-SAM-157 | SRR32559025 |

`BULK-SAM-148` is benchmarked first because its complete FASTQ pair is the
smallest in local9 at 7.794 GiB. It is not analyzed as a standalone biological
result.

## Computational design

1. Install the project-local `Rsubread` binary.
2. Require at least 30 GiB free for the benchmark alignment phase.
3. Download and MD5-verify GENCODE v47 GRCh38 primary-assembly FASTA and GTF.
4. Build a gapped, split Rsubread index with a 2,000 MB index block.
5. Download one complete paired library by either deposited ENA FASTQ transfer
   or NCBI SRA-object transfer.
6. For ENA, verify both deposited FASTQ MD5 values and the locked total byte
   count. For SRA, validate the SRA object with `vdb-validate`, convert it to
   paired FASTQ with `fasterq-dump`, and record local FASTQ MD5 values.
7. Align with Subjunc using one thread and up to 16 equally best locations.
8. Count genes in unstranded, forward, and reverse modes and record all three.
   Select a directional mode only when its forward/reverse assignment ratio is
   at least 4 and it retains at least 80% of the unstranded assignments;
   otherwise use unstranded counts.
9. Require a nonempty BAM, a nonempty junction BED, and at least 10% assigned
   fragments relative to expected read pairs.
10. Write `ALIGNMENT_COMPLETE.tsv`; retain inputs until the marker is audited.

Reporting up to 16 best locations prevents arbitrary single-locus treatment of
ordinary multi-mappers. Gene-level featureCounts excludes multi-mapping reads.
Consequently, this workflow is not a definitive SMN1-versus-SMN2
allele-specific quantification method; highly homologous SMN reads require a
separate locus-aware analysis.

## Deadline-constrained all-nine route

If nine retained BAMs cannot fit locally or the measured Subjunc runtime
exceeds the study deadline, complete raw FASTQ coverage takes priority over
partial-library genomic alignment. The expedited route retains completed and
interrupted BAMs, keeps the genomic benchmark as orthogonal feasibility
evidence, and performs Salmon transcriptome quantification for all nine frozen
libraries.

The expedited computation is pre-specified as follows:

1. Download the official GENCODE v47 transcript FASTA and verify its official
   release MD5.
2. Build a Salmon 1.10.2 transcriptome-only index with `--gencode`.
3. Require both deposited ENA FASTQ MD5 values for every library.
4. Quantify complete paired FASTQs with automatic library-type detection,
   `--validateMappings`, `--seqBias`, and `--gcBias`.
5. Aggregate transcripts to GENCODE genes with Bioconductor `tximport` using
   `countsFromAbundance="lengthScaledTPM"`.
6. Apply the same donor-aware disease and treatment models, five disease LOLO
   folds, two linewise treatment effects, concordance analysis, and
   pre-specified candidate ranking.
7. Validate the frozen sample order, model units, outputs, and absence of
   random splits before writing
   `LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv`.

This marker records
`quantification_method=salmon_transcriptome_quasimapping` and
`whole_genome_alignment=FALSE`. It is deliberately distinct from
`LOCAL9_RAW_ANALYSIS_COMPLETE.tsv`. The approach supports an all-nine
raw-read re-quantification claim, but not an all-nine genome realignment,
splice-junction discovery, or allele-specific SMN1/SMN2 claim. Because the
deadline workflow uses a transcriptome-only rather than decoy-aware index,
possible genomic sequence misassignment remains a stated limitation and
raw-versus-processed concordance is reported as an empirical sensitivity
check.

## Feasibility decision

Proceed from one sample to all nine only when:

- the reference checksums and raw-input integrity checks pass;
- Rsubread finishes without exhausting memory;
- BAM and junction outputs are nonempty;
- the featureCounts assignment gate passes;
- peak disk use leaves enough room for the next sample; and
- elapsed time is practical for sequential processing.

Failure at any gate ends the raw workstation attempt. The manuscript then
continues as a processed-data, donor-aware, cross-model hypothesis-generating
study with the absent raw validation reported as a limitation.

## Statistical analysis after nine samples

The disease contrast remains three control versus two SMA donor lines. The
treatment contrast remains two paired SMA lines. Analysis must therefore
emphasize effect sizes and linewise consistency:

- disease: donor-line pseudobulk effect plus five leave-one-line-out checks;
- treatment: paired overall effect plus separate S2 and S3 effects;
- genome-wide raw-versus-processed effect correlations and direction
  agreement for both contrasts;
- raw candidate support only when the natural-resilience direction is
  recovered and the disease LOLO and treatment donor-line checks agree;
- no random train/test split;
- no claim that technical preparations increase biological sample size;
- no therapeutic efficacy claim.

The processed matrix is used only as a comparison target after raw alignment
and quantification. It never supplies, imputes, or replaces a raw count. Since
both analyses use the same GSE290979 biological material, agreement is
sensitivity confirmation and not independent validation.

After these analyses, a separate integrity gate checks the frozen library
order, biological-unit counts, complete LOLO and linewise outputs, concordance
tables, candidate tables, and absence of random splits. It writes
`LOCAL9_RAW_ANALYSIS_COMPLETE.tsv` only after every computational check passes.
The marker records completion regardless of effect direction or significance;
it is not a biological-success criterion.

## Manuscript wording

Allowed after all nine libraries pass:

> A pre-specified nine-library raw-read sensitivity analysis retained every
> available donor line while selecting one centrally sequenced preparation per
> line-condition for workstation-feasible Subjunc alignment.

Allowed after only the expedited marker passes:

> Complete checksum-verified FASTQ pairs from nine pre-specified libraries were
> re-quantified against the GENCODE v47 transcriptome with Salmon, followed by
> donor-aware differential analysis and leave-one-donor-line-out sensitivity
> checks.

Not allowed:

- "complete raw-read reanalysis";
- "all-nine genome realignment" when only the expedited marker passes;
- "independent cohort validation";
- "genome-wide raw rMATS replication";
- "proof of an SMA treatment target."

## Completed Results

The completed Salmon analysis quantified all nine frozen libraries and
retained the planned five disease donor-line units and two paired treatment
lines. Genome-wide raw-versus-processed comparisons found:

| Contrast | Matched genes | Spearman rho | Pearson r | Direction agreement |
| --- | ---: | ---: | ---: | ---: |
| SMA versus control | 14,636 | 0.778 | 0.784 | 0.807 |
| R6-MO versus scramble | 14,385 | 0.588 | 0.553 | 0.732 |

Among the 37 frozen candidates, 35 mapped in both raw contrasts, 21 retained
the full raw direction pattern, and nine also passed five disease LOLO folds
and both treatment-line checks. Four candidates (`LY6H`, `HS3ST5`, `ZNF853`,
`IL17D`) were unit-robust in both raw and processed analyses.

All verified FASTQs and existing BAMs remain preserved. These results support
gene-level raw re-quantification language only. They do not change the
requirements for genome-wide raw splicing, junction, or SMN1/SMN2
paralog-specific claims.
