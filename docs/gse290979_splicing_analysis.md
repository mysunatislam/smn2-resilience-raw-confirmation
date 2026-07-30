# GSE290979 SMN2 and Splicing Analysis

**Analysis date:** 21 July 2026
**Status:** Reproducible reanalysis of official sample-level supplementary
rMATS results; the raw-first reanalysis workflow is implemented but requires
HPC execution, and experimental validation is still required

## Executive Summary

This phase extends the human-first OculoRescue-SMA project from gene-level
expression to alternative splicing. It analyzes all 31 human libraries in
[GSE290979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290979)
through the official sample-level rMATS event tables and source data associated
with the [primary publication](https://www.nature.com/articles/s41467-025-67725-1).

The deposited run manifest contains 31 paired-end RNA-seq runs, 2.378 billion
read pairs, and 285.3 GiB of compressed FASTQ. The raw FASTQ files were not
downloaded because they exceed available local storage. Instead, this analysis
uses the publication's sample-level junction counts, inclusion levels,
significance statistics, and source-data SMN assays. The
[ENA project record](https://www.ebi.ac.uk/ena/browser/view/PRJNA1231113)
and MD5-verified supplementary files provide provenance.

Main results:

- 10,553 disease-associated and 2,614 treatment-associated significant
  splicing events were present in the official workbooks.
- Full structural event matching identified 182 events significant in both
  comparisons.
- After orienting the treatment effect as R6-MO minus scramble, 112 events
  reversed the disease direction.
- Of those, 108 also moved closer to the untreated control mean.
- A stringent requirement for reversal and improvement in both SMA lines
  retained 83 events across 74 genes.
- Seventy of the 74 strict genes are current HGNC protein-coding genes.
- No strict splicing gene has robust OMN-positive expression in the
  authoritative R analysis at FDR 0.05.
- Three strict genes have nominal, same-platform-supported OMN-positive R
  effects and remain exploratory: **KCNAB3, LINC00665**, and **TSPOAP1**.
- Source-data assays strongly support R6-MO target engagement: the pooled
  full-length-to-SMN-delta7 transcript ratio rises from 7.02 with scramble to
  103.63 with R6-MO.

A critical audit finding is that the treatment rMATS workbook stores
`IncLevelDifference` as scramble minus R6-MO, although the published rescue
workbook labels the copied treatment values as R6-MO minus scramble. Under the
sample orientation supported by the deposited data, all 111 published rescue
pairs have treatment changes in the same biological direction as disease,
rather than the opposite direction. In addition, only 67 of the 111 pairs have
identical full rMATS event structures; 44 match only the reduced coordinates
included in the rescue workbook.

This is a computational audit finding, not yet a conclusion that the article is
wrong. A raw FASTQ or BAM rerun and, ideally, clarification from the authors are
required before presenting it publicly as a corrected result.

## 1. Research Question

The analysis asks two linked questions:

1. Which splicing abnormalities distinguish human SMA spinal-cord organoids
   from control organoids?
2. Which abnormalities are directionally corrected by an SMN2-targeting
   R6-morpholino in both SMA donor lines?

The project then asks whether these splice-restoration genes overlap the
expression program of relatively resistant human oculomotor neurons.

This design separates two possible therapeutic concepts:

- **SMN and splice restoration:** correct SMN2 exon 7 and downstream
  SMN-sensitive RNA processing.
- **Motor-neuron resilience:** support pathways naturally enriched in
  relatively resistant oculomotor neurons.

Neither computational track demonstrates a treatment by itself.

## 2. Dataset and Provenance

### RNA-seq cohort

| Group | Libraries | Donor lines |
| --- | ---: | ---: |
| Control, untreated | 9 | C1, C2, C3 |
| SMA, untreated | 6 | S2, S3 |
| SMA, scramble | 8 | S2, S3 |
| SMA, R6-MO | 8 | S2, S3 |
| **Total** | **31** | **5** |

The 20+ sample request is met at the library level. It is not met at the
independent-patient level. Treatment consistency across S2 and S3 is useful,
but these remain only two SMA donor lines.

### Official resources used

- Supplementary Data 7: SMA versus control rMATS events.
- Supplementary Data 8: scramble and R6-MO rMATS events.
- Supplementary Data 9: published rescued-event pairs.
- Article source data: full-length and SMN-delta7 assays.
- GEO count matrix and sample metadata.
- ENA run-level read counts, FASTQ URLs, sizes, and MD5 hashes.

The preparation script verifies the supplementary-file MD5 values and writes
the complete run manifest before analysis.

## 3. Orientation Audit

Correct effect orientation is essential. rMATS reports
`IncLevelDifference = IncLevel1 - IncLevel2`; it does not infer which direction
is biologically called "treatment rescue."

### Disease workbook

The disease workbook has six values in `IncLevel1` and nine in `IncLevel2`.
Anonymous rMATS columns were matched to GEO libraries using genome-wide
junction-count profiles and expression profiles with a global Hungarian
assignment.

- Group 1 maps to all six untreated SMA libraries.
- Group 2 maps to all nine untreated control libraries.
- Disease delta is therefore **SMA minus control**.

### Treatment workbook

The treatment workbook has eight values in each group.

- Group 1 maps to all eight scramble libraries.
- Group 2 maps to all eight R6-MO libraries.
- The inferred 16-sample order exactly reproduces every explicitly named
  per-sample PSI column in all five event sheets.
- Deposited treatment `IncLevelDifference` is therefore **scramble minus
  R6-MO**.
- The biologically oriented treatment delta is **R6-MO minus scramble**, which
  requires multiplying the deposited value by -1.

The program asserts all these relationships and stops if they do not hold.

### Published rescue workbook audit

The rescue workbook contains 222 rows forming 111 disease-treatment pairs and
97 distinct normalized gene symbols. Its treatment rows are labeled
`SMAmo-SMAscr`, but their values exactly equal the scramble-minus-R6 values in
the treatment source workbook.

After biological reorientation:

- 111 of 111 published pairs change in the same direction as the disease
  difference.
- 0 of 111 reverse the disease direction.
- 67 pairs match the same full rMATS event structure.
- 44 pairs share the reduced displayed coordinates but differ elsewhere in the
  full event definition.
- No strict corrected event uses the same disease and treatment source-event
  ID pair as the published rescue set.
- The strict and published sets overlap at the gene-symbol level only for
  **IFI27L1** and **NCALD**, through different source-event pairs.

These observations motivate raw-read confirmation. They should not be framed
as misconduct or as a definitive correction without that confirmation.

## 4. Event Reconstruction and Filters

Events were matched using complete rMATS structural keys:

- ES: target exon plus upstream and downstream exon boundaries.
- A5SS and A3SS: long and short exon boundaries plus flanking exon.
- MXE: both mutually exclusive exons plus both flanking exons.
- RI: retained intron boundaries plus upstream and downstream exons.

Gene, chromosome, and strand are also included. Matching only the target exon
coordinates can pair different transcript structures and was therefore not
used for the corrected set.

All source events already satisfy the publication's significance criteria.
The following additional filters were applied sequentially:

1. **Exact common:** same full event is significant in disease and treatment.
2. **Direction reversal:** `(SMA - control) * (R6 - scramble) < 0`.
3. **Pooled correction:** direction reverses and pooled R6 PSI is closer to the
   untreated control mean than pooled scramble PSI.
4. **Strict correction:** the pooled criterion holds, and S2 and S3 each
   independently reverse direction and move closer to control.

### Event counts

| Event type | Exact common | Direction reversal | Pooled correction | Strict both lines |
| --- | ---: | ---: | ---: | ---: |
| Exon skipping | 80 | 48 | 46 | 35 |
| Alternative 5' splice site | 15 | 8 | 8 | 6 |
| Alternative 3' splice site | 15 | 7 | 7 | 5 |
| Mutually exclusive exon | 36 | 28 | 27 | 22 |
| Retained intron | 36 | 21 | 20 | 15 |
| **Total** | **182** | **112** | **108** | **83** |

The strict filter is a consistency filter, not a formal donor-level
replication test. Four libraries per treatment within each line do not turn two
SMA lines into eight independent patients.

## 5. SMN2 Target Engagement

The article source data provide an orthogonal full-length-to-SMN-delta7
transcript ratio.

| Source-data line | Untreated mean | Scramble mean | R6-MO mean | R6 / scramble |
| --- | ---: | ---: | ---: | ---: |
| S1 | 0.985 | 6.940 | 229.540 | 33.08 |
| S2 | 1.007 | 6.724 | 53.938 | 8.02 |
| S3 | 1.180 | 7.398 | 27.404 | 3.70 |
| **Pooled** | **1.058** | **7.021** | **103.627** | **14.76** |

This supports strong molecular activity of R6-MO at the SMN transcript level.
The downstream event-orientation concern does not negate this target-engagement
result.

Because SMN1 and SMN2 are highly homologous, the RNA-seq event tables alone are
not used here to claim an SMN2-specific exon 7 PSI. The source-data RT-PCR and
RFLP assays are the stronger evidence for target engagement until the raw reads
are reanalyzed with paralog-aware methods.

## 6. Two-Track Candidate Results

### Splicing-restoration track

The strict set contains 83 events across 74 genes:

- 70 current HGNC protein-coding genes.
- 1 long non-coding RNA, **LINC00665**.
- 1 readthrough transcript, **ARPC4-TTLL3**.
- 2 source symbols without a current unambiguous HGNC mapping.

Genes with multiple strict events include:

| Gene | Strict events | Event classes |
| --- | ---: | --- |
| POFUT2 | 5 | A5SS, MXE |
| SEPTIN10 | 2 | MXE |
| DENND3 | 2 | ES, MXE |
| ORC3 | 2 | ES |
| ANAPC10 | 2 | MXE |
| SREK1IP1 | 2 | ES |

Large event-level restoration scores also occur for **ABCA1, FAT3, PDE9A,
TNPO2, MRFAP1, OGA**, and **IFI27L1**. These scores rank statistical and
directional evidence. They do not establish that increasing or decreasing the
gene is beneficial.

### Expression-resilience bridge

The R analysis finds no strict splicing gene with robust OMN-positive
expression at FDR 0.05. Three genes meet an explicitly exploratory criterion:
primary OMN `p < 0.05`, positive same-platform effect, and same-platform
`p < 0.10`.

| Gene | Locus | R OMN effect | Primary p | Primary FDR | Same-platform p | Strict event |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| KCNAB3 | Protein coding | 6.058 | 0.0078 | 0.667 | 0.0353 | ES |
| LINC00665 | Long non-coding RNA | 3.732 | 0.0243 | 0.763 | 0.0943 | ES |
| TSPOAP1 | Protein coding | 3.530 | 0.0386 | 0.817 | 0.0635 | MXE |

All three FDR values are high. They are hypotheses for junction-level and
cell-type-specific validation, not statistically validated resilience genes.
KCNAB3 and TSPOAP1 also decrease at the R gene-expression effect level after
R6-MO, so splice correction and total-expression behavior must be treated as
different mechanisms.

The cross-model extension retains GSE93939 as the natural-resistance anchor
and adds the independent GSE108094 iPSC-SMA motor-neuron disease direction.
LINC00665 is externally depleted in SMA and is the only strict splice gene that
also meets the exploratory GSE93939 criterion and all four total-expression
directions. KCNAB3 and TSPOAP1 retain external disease support but not the
GSE290979 total-expression R6 alignment. Full results are in
[`cross_model_human_motor_neuron_resilience.md`](cross_model_human_motor_neuron_resilience.md).

The lack of a robust R expression bridge argues against compressing expression
and splicing into one apparently definitive target list.

### Pathway sensitivity analysis

g:Profiler enrichment was run with two custom backgrounds:

1. Protein-coding genes represented among the disease-significant splicing
   events.
2. All current HGNC protein-coding genes.

No pathway is significant with the disease-splicing background. Three terms
are significant only with the broader HGNC background: cilium assembly, cell
cycle, and spina bifida. They are background-sensitive and exploratory, not
robust pathway conclusions.

## 7. A Publishable Research Strategy

The strongest novel direction is not simply another candidate-gene list. It is
an **orientation-aware, full-structure reanalysis of SMN-restored human
organoid splicing integrated with motor-neuron resilience**.

Two manuscript paths are realistic:

### Path A: Reanalysis and methods/resource paper

Required additions:

1. Re-download or cloud-stream all 31 raw runs.
2. Reproduce the published STAR and rMATS pipeline with explicit sample order.
3. Repeat with a second splice engine such as MAJIQ or LeafCutter.
4. Confirm the sign and structural-pairing observations.
5. Release sample sheets, event keys, code, environment, and checksums.
6. Seek author clarification before describing the discrepancy as an error.

This path could be publishable as a focused data reanalysis if the finding
survives raw-read reproduction.

### Path B: Targeted validation plus cross-model paper

Required additions:

1. Run the locked 21-library profile containing all treatment libraries and
   one median-depth untreated preparation per donor line.
2. Use its raw STAR counts for selected-library expression sensitivity and its
   compact event-locus BAMs only for the frozen 83-event recount.
3. Keep the deposited processed matrix explicitly labeled as the genome-wide
   discovery source.
4. Report targeted rMATS FDR as selection-biased and interpret effect
   direction and junction support instead.
5. Add independent cell-resolved and disease-direction evidence, plus
   orthogonal junction assays for the frozen panel.

This path reduces the raw input from 285.3 to 193.116 GiB without pretending
that partial processing is complete genome-wide reanalysis.

### Path C: Mechanistic SMA biology paper

Required additions:

1. Validate the frozen 12-event panel by junction assays in independent SMA
   and control iPSC motor-neuron lines.
2. Measure whether SMN2 correction normalizes each event at RNA and protein
   levels.
3. Perturb selected splice isoforms without changing total gene abundance.
4. Test motor-neuron survival, axon integrity, electrophysiology, and
   neuromuscular-junction phenotypes.
5. Compare SMN2 correction alone, isoform perturbation alone, and combination
   treatment.
6. Include oculomotor-like and spinal motor-neuron identities when feasible.

The R workflow now freezes a transparent 12-event primary panel:
**MRFAP1, IFI27L1, PDE9A, HERPUD1, ABCA1, SAT2, TNPO2, SEPTIN11, FAT3, OGA,
COL5A2**, and **EIF5**. It contains 12 unique approved HGNC genes and represents
ES, A5SS, A3SS, MXE, and RI events. The exploratory cross-track bridge remains
a separate three-event family: **KCNAB3, LINC00665**, and **TSPOAP1**. NCALD is
retained only as an orientation-audit comparator using exact source-event
structures.

The scoring formula, exact panel metrics, assay-development steps, donor-level
statistics, and decision gates are in
[`splicing_validation_panel.md`](splicing_validation_panel.md). The panel is
for preclinical validation, not patient treatment or therapeutic ranking.

## 8. Raw-Read Confirmation Protocol

This is a claim-dependent publication gate rather than an optional
confirmation phase. `full31` rebuilds total-expression and discovery-splicing
inference without the processed workbooks. `validation21` rebuilds selected-
library expression from raw counts and performs only a pre-specified 83-event
recount. Execution status and both completion contracts are in
[`gse290979_local9_protocol.md`](gse290979_local9_protocol.md).

### Compute and storage

The ENA manifest reports 285.3 GiB of compressed FASTQ. A practical rerun
should use an HPC or cloud workspace with substantially more capacity for
FASTQ, reference indices, temporary files, coordinate-sorted BAMs, and splice
outputs. Every downloaded file should be checked against its ENA MD5.

The reduced profile processes 193.116 GiB and uses a 350 GiB free-space gate
with two concurrent streaming alignments. It deletes complete BAMs and FASTQs
only after full-alignment counts, QC, RSeQC, and compact event-locus BAMs pass.

The pre-specified full31 STAR/rMATS workflow remains a claim-dependent future
analysis. It has not been run locally because the required storage exceeds
this workspace.

### Alignment and quantification

1. Use the same genome build and annotation reported by the article for the
   primary reproduction.
2. Run STAR two-pass with fixed software and reference versions.
3. Preserve multimapping information because SMN1 and SMN2 are highly
   homologous.
4. Re-run rMATS with explicit sample lists:
   - disease: six SMA versus nine control;
   - treatment: eight R6-MO versus eight scramble, with the contrast name and
     effect direction written into the sample sheet.
5. Calculate biological effects directly from sample PSI:
   - disease = SMA minus control;
   - treatment = R6-MO minus scramble.
6. Validate full event coordinates across comparisons.
7. Repeat event discovery with an independent splice method.
8. Quantify SMN exon 7 junctions and paralog-discriminating sequence features,
   while reporting ambiguous reads rather than assigning them arbitrarily.

### Statistics

1. Use donor line as the biological unit.
2. Report sample-level PSI but avoid treating preparation replicates as
   independent patients.
3. Require effect-size direction consistency in S2 and S3.
4. Correct for multiple testing within the prespecified event family.
5. Report confidence intervals and raw values, not only FDR thresholds.
6. Treat the two-line result as discovery requiring independent donor lines.

### Experimental confirmation

For each selected event:

1. Design primers in constitutive flanking exons.
2. Confirm a single expected product structure by sequencing.
3. Quantify isoform ratios by capillary electrophoresis or digital PCR.
4. Measure total transcript abundance separately.
5. Test protein isoform or pathway output when an assay exists.
6. Repeat in at least three independent SMA donor lines and matched controls.
7. Blind image-based survival and neurite analyses to treatment group.

## 9. Limitations

- The downstream analysis uses official rMATS tables, not a fresh raw alignment.
- Only significant source events are available, so pathway backgrounds are
  imperfect.
- Two SMA donor lines limit biological generalization.
- Bulk organoids mix motor neurons, other neurons, glia, and progenitors.
- PSI correction does not prove restoration of protein function.
- SMN1/SMN2 paralog ambiguity requires dedicated raw-read handling.
- Candidate scores are prioritization tools, not causal or clinical evidence.
- One ECD source event has Excel-coerced IJC vectors; its validation score uses
  a labeled conservative SJC-only support lower bound pending raw-read recovery.
- The published-workbook sign finding requires independent reproduction and
  author clarification.

## 10. Reproducible Outputs

Key files:

- `results/splicing/GSE290979_splicing_summary.json`
- `results/splicing/GSE290979_rmats_sample_orientation.tsv`
- `results/splicing/GSE290979_common_splicing_events.tsv`
- `results/splicing/GSE290979_strict_corrected_rescue_events.tsv`
- `results/splicing/GSE290979_published_rescue_sign_audit.tsv`
- `results/candidate_ranking/human_strict_splicing_candidates.tsv`
- `results/candidate_ranking/human_two_track_summary.tsv`
- `results/pathway_enrichment/human_strict_splicing_gprofiler.tsv`
- `results/r/splicing/GSE290979_R_validation_event_ranking.tsv`
- `results/r/splicing/GSE290979_R_primary_validation_panel.tsv`
- `results/r/splicing/GSE290979_R_exploratory_OMN_bridge_panel.tsv`
- `results/r/splicing/GSE290979_R_validation_score_specification.tsv`
- `results/r/figures/GSE290979_R_validation_priority_top20.png`
- `results/r/figures/GSE290979_R_validation_panel_components.png`

Run:

```powershell
.\.venv\Scripts\python.exe scripts\06_prepare_gse290979_splicing.py
.\.venv\Scripts\python.exe scripts\07_analyze_gse290979_splicing.py
.\.venv\Scripts\python.exe scripts\08_integrate_splicing_resilience.py
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\11_rank_splicing_validation_events.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\06_validate_outputs.R
```

## Bottom Line

R6-MO strongly improves the full-length-to-SMN-delta7 transcript ratio, and a
careful reorientation of the deposited event tables identifies 83 downstream
splicing events that move toward control in both analyzed SMA lines. These
events provide a focused human validation set.

The most journal-relevant observation is the combination of sample-orientation
forensics, full-event structural matching, and two-track biological
integration. It is worth pursuing, but raw-read reproduction is the gate that
must be passed before making a strong publication claim.
