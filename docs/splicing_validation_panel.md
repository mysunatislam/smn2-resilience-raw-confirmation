# GSE290979 Splicing Validation Panel and Experimental Protocol

**Analysis-freeze date:** 21 July 2026
**Scope:** preclinical human SMA research
**Status:** computational panel selected; raw-read and wet-lab confirmation
remain outstanding

## Decision in Plain Language

This project is worth continuing, but the next claim must be narrow. The
current evidence supports a reproducible set of candidate splice events, not a
new treatment and not proof that manipulating any listed gene will protect a
motor neuron.

The immediate study has three modules:

1. **Raw-read confirmation:** first run the locked 21-library profile to test
   the exact 83 event structures; use the full 31-library profile before any
   genome-wide raw-discovery claim.
2. **Primary experimental panel:** validate 12 splice events selected with a
   transparent, fixed score and representation of all five rMATS event classes.
3. **Exploratory OMN bridge:** test KCNAB3, LINC00665, and TSPOAP1 separately as
   hypotheses connecting SMA splice restoration to oculomotor-neuron biology.

The 12-event panel is the confirmatory assay-development family. The three OMN
bridge events are a separate exploratory family and must not be added to the
primary panel after outcomes are known.

## 1. Evidence Being Frozen

The input is the R-defined set of 83 strict events across 74 genes. Every event:

- is significant in the official disease and treatment workbooks;
- has the same complete rMATS event structure in both workbooks;
- changes in SMA relative to control;
- changes in the opposite biological direction after R6-MO;
- moves closer to control after R6-MO in pooled samples; and
- reverses direction and moves toward control separately in S2 and S3.

This is still evidence from official processed tables. The full analysis and
orientation audit are described in
[`gse290979_splicing_analysis.md`](gse290979_splicing_analysis.md). The raw-read
workflow is in
[`local/gse290979_raw/README.md`](../local/gse290979_raw/README.md).

## 2. Transparent Validation-Priority Score

The score ranks assay candidates from 0 to 1. It was defined retrospectively
for this project and is now frozen before raw-read confirmation and wet-lab
testing. It is not described as a prospective preregistration.

| Component | Weight | Definition |
| --- | ---: | --- |
| Minimum effect | 0.30 | Smaller of absolute disease and treatment delta-PSI, scaled to 0.20 and capped at 1 |
| Significance | 0.25 | Mean capped `-log10(FDR)/5` from disease and treatment workbooks |
| Both-line consistency | 0.25 | Smaller S2/S3 improvement toward control, scaled to 0.10 and capped at 1 |
| Junction support | 0.15 | Group-median informative junction support plus fraction of samples with at least 10 reads |
| Assay feasibility | 0.05 | ES 1.00, A5SS/A3SS 0.90, RI 0.85, and MXE 0.75 |

Selection first takes the highest-scoring approved HGNC event in each event
class, then fills the panel by score while requiring a different gene for each
entry. This prevents a technically convenient event class or a multi-event
gene from occupying the whole panel.

The score ranks **validation feasibility**. It does not estimate therapeutic
benefit, causality, protein consequence, or safety.

### Source-data count caveat

Eighty-two events have exact comma-separated IJC and SJC vectors. Excel
coerced the four ECD IJC vectors into large scientific-notation numbers and
destroyed some per-sample digits. ECD is therefore scored with a conservative
SJC-only support lower bound and is explicitly flagged in the ranking table.
No IJC values were guessed. ECD ranks 80th and is not in either validation
panel; its support must be recovered from raw FASTQ before quantitative use.

## 3. Primary 12-Event Panel

`dPSI disease` is SMA minus control. `dPSI R6` is R6-MO minus scramble. Opposite
signs therefore indicate directional correction.

| Order | Gene | Class | Score | dPSI disease | dPSI R6 | Minimum group-median junction support | Support |
| ---: | --- | --- | ---: | ---: | ---: | ---: | --- |
| 1 | MRFAP1 | RI | 0.949 | -0.237 | 0.184 | 21.0 | Moderate |
| 2 | IFI27L1 | ES | 0.870 | -0.318 | 0.196 | 10.0 | Low |
| 3 | PDE9A | ES | 0.868 | -0.247 | 0.277 | 4.5 | Low |
| 4 | HERPUD1 | A5SS | 0.853 | 0.449 | -0.177 | 22.5 | Moderate |
| 5 | ABCA1 | RI | 0.852 | -0.375 | 0.309 | 5.0 | Low |
| 6 | SAT2 | ES | 0.848 | -0.502 | 0.323 | 8.5 | Low |
| 7 | TNPO2 | A3SS | 0.847 | -0.414 | 0.281 | 6.0 | Low |
| 8 | SEPTIN11 | MXE | 0.840 | -0.362 | 0.200 | 17.0 | Low |
| 9 | FAT3 | MXE | 0.840 | 0.341 | -0.384 | 7.0 | Low |
| 10 | OGA | MXE | 0.828 | 0.361 | -0.313 | 6.0 | Low |
| 11 | COL5A2 | RI | 0.817 | 0.258 | -0.162 | 24.5 | Moderate |
| 12 | EIF5 | A5SS | 0.816 | -0.114 | 0.300 | 16.5 | Low |

Low source-table support is a reason to perform a sensitive junction assay; it
is not evidence that the event is false. Exact coordinates, FDR values,
per-line consistency metrics, and selection reasons are in
`results/r/splicing/GSE290979_R_primary_validation_panel.tsv`.

## 4. Separate Exploratory OMN Bridge

No strict splice-restoration gene is OMN-enriched at FDR 0.05 in the adjusted
GSE93939 R model. These three genes pass only the frozen exploratory OMN
criterion. Their high OMN FDR values must accompany every presentation.

| Gene | Class | Splicing rank | Score | dPSI disease | dPSI R6 | OMN log2 effect | OMN p | OMN FDR |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LINC00665 | ES | 24 | 0.744 | 0.138 | -0.212 | 3.732 | 0.0243 | 0.763 |
| KCNAB3 | ES | 63 | 0.590 | 0.177 | -0.158 | 6.058 | 0.0078 | 0.667 |
| TSPOAP1 | MXE | 81 | 0.420 | -0.101 | 0.144 | 3.530 | 0.0386 | 0.817 |

The TSPOAP1 source workbook symbol is BZRAP1 and is normalized through HGNC.
These events test a biological bridge; they are not replacements for the
higher-scoring primary events.

All three have usable GSE108094 external iPSC-SMA motor-neuron effects.
LINC00665 is the only one with the complete total-expression direction across
GSE93939, GSE290979 disease, GSE290979 R6 response, and GSE108094. KCNAB3 and
TSPOAP1 remain splice-level hypotheses because their GSE290979 total-expression
R6 effects do not align with the natural-resistance direction. The cross-model
analysis is documented in
[`cross_model_human_motor_neuron_resilience.md`](cross_model_human_motor_neuron_resilience.md).

## 5. Stage 1: Raw-Read Confirmation Gate

Run the locked HPC workflow before ordering a large primer panel or beginning
causal perturbation work.

### Required inputs and provenance

1. For targeted panel validation, use all 21 runs in the locked `validation21`
   sheet: every treatment library and one median-depth untreated preparation
   per donor line. For genome-wide raw discovery, use all 31 runs.
2. Verify both FASTQ MD5 values for every run.
3. Record the exact hg38 FASTA and GTF labels and SHA-256 hashes.
4. Use STAR 2.7.10a and rMATS turbo 4.3.0 for the primary confirmation.
5. Record the 151-bp deposited read length and the article/ENA instrument
   discrepancy rather than silently harmonizing it.

### QC gate

1. Require a validated retained BAM, BAM index, STAR final log, full-alignment
   Samtools flagstat, and RSeQC output for every selected library.
2. Require every selected library to support the configured stranded
   orientation at the frozen 0.80 RSeQC fraction threshold and by at least a
   0.20 margin over the opposite orientation.
3. Flag, inspect, and explain any sample with less than 50% uniquely mapped
   reads. This is a review threshold, not an automatic biological exclusion.
4. Inspect group-level mapping, insert-size, duplication, and gene-body
   patterns before interpreting differential splicing.
5. Do not remove a sample because it weakens the desired result. Any exclusion
   rule must be technical, documented, and applied without seeing panel
   outcomes.

### Two analysis tracks

1. **Discovery track:** test all reference-derived events. Discovery FDR can
   support independent statistical replication.
2. **Fixed-event recount:** quantify the 83 locked structures. This tests
   recoverability, effect direction, junction depth, and S2/S3 consistency.
   Its FDR is selection-biased and is not independent evidence.

An official event receives full raw confirmation only when the same complete
structure is recovered, disease and treatment absolute delta-PSI are each
greater than 0.05, the effects have opposite signs, both S2 and S3 move toward
control, and both discovery-track FDR values are below 0.05. Events failing
full confirmation remain reportable as targeted directional support or not
reproduced; they must not be silently dropped.

### Independent splice engine

After the STAR/rMATS result is frozen, rerun the discovery question with one
independent method such as MAJIQ or LeafCutter. Coordinate conversion and event
matching rules must be written before comparing results. Agreement should be
reported at junction and gene levels because event definitions differ among
engines.

## 6. Stage 2: Junction-Assay Development

Do not fabricate primer sequences from a gene name. Lock the raw-confirmed
hg38/GTF coordinates and inspect the relevant transcripts first.

### Assay design by event class

| Class | Preferred first assay |
| --- | --- |
| ES | RT-PCR primers in constitutive exons flanking the skipped exon, followed by capillary sizing |
| A5SS/A3SS | Shared constitutive primer plus junction-specific primer for each splice-site choice |
| MXE | Two junction-specific assays, one for each mutually exclusive exon, with a common flanking primer where possible |
| RI | One exon-intron junction assay plus a fully spliced junction assay; include genomic-DNA controls |

For every assay:

1. Confirm that primer-binding regions are constitutive in the locked
   annotation and do not contain common variants in the donor lines.
2. Check genome-wide specificity and known pseudogenes/paralogs in silico.
3. Avoid assigning ambiguous SMN1/SMN2 reads; use a dedicated paralog-aware
   SMN assay for target engagement.
4. Prefer short, non-overlapping amplicons for qPCR/ddPCR and resolvable size
   differences for capillary RT-PCR.
5. Verify one representative product from each isoform by Sanger sequencing.
6. Establish linear range, limit of detection, and technical precision before
   unblinding biological groups.
7. Require no-template and no-reverse-transcriptase controls. RI assays also
   require a genomic-DNA contamination check.

Capillary electrophoresis is the most direct first choice for ES events with
resolvable products. Junction qPCR or ddPCR is preferable when products are
low abundance, differ by only a few nucleotides, or require isoform-specific
junctions. Use the same quantification platform across biological groups for a
given event.

## 7. Stage 3: Independent Human Motor-Neuron Validation

### Biological material

- At least three independent SMA donor lines and three independent control
  lines are recommended for the first validation cohort.
- Isogenic corrected pairs are valuable but do not replace independent genetic
  backgrounds.
- Use at least two independent differentiation batches per line. Batches are
  nested technical/experimental replication; the donor line remains the main
  biological unit.
- Record SMN1 genotype, SMN2 copy number, sex, passage, differentiation batch,
  motor-neuron purity, RNA integrity, and treatment exposure.

### Initial confirmation groups

1. Control motor neurons, untreated or matched vehicle.
2. SMA motor neurons, untreated.
3. SMA motor neurons plus scramble oligo.
4. SMA motor neurons plus the SMN2 splice-correcting research condition.

The initial purpose is to confirm disease association and SMN-responsive
normalization. Candidate-specific perturbation should begin only after the
event structure and assay are validated.

### Required molecular measurements

- PSI or isoform fraction for all 12 primary events.
- The three bridge events as a separately labeled exploratory family.
- Total transcript abundance for each tested gene.
- Full-length SMN and SMN-delta7 transcript ratio with a paralog-aware assay.
- Total SMN protein and, where feasible, nuclear gems as target-engagement
  measurements.
- Cell-identity and purity markers appropriate to the differentiation.

Include NCALD as an orientation-audit comparator using the exact source-event
structures. NCALD is not a member of the strict 83-event set and must not be
counted as a primary success.

### Blinding and allocation

1. Randomize wells across plates and processing order within each line and
   batch.
2. Blind sample labels during capillary peak calling, imaging, and primary QC.
3. Define technical replicate aggregation and failed-well rules before
   unblinding.
4. Keep all raw electropherograms, amplification curves, exclusion reasons,
   and sample metadata.

## 8. Statistical Analysis Plan

### Primary endpoint

The primary endpoint for each event is PSI. Report raw sample values, donor
means, delta-PSI, confidence intervals, and exact p values.

For an event to show directional experimental confirmation:

1. SMA versus control has the source-consistent sign.
2. SMN2 correction versus scramble has the opposite sign.
3. SMN2-corrected PSI is closer to the control mean than scramble PSI.
4. The direction is not driven by only one donor line.

Analyze donor line as the biological unit. A mixed model may include condition
as a fixed effect and donor line and differentiation batch as structured
effects when the design supports them. For bounded PSI near 0 or 1, use a
model appropriate for proportions or a justified transformation; do not apply
ordinary Gaussian tests automatically.

Control Benjamini-Hochberg FDR across the 12 primary events for each frozen
primary contrast. Analyze the three bridge events as a separate exploratory
family. Do not pool technical wells as independent observations.

A useful descriptive rescue fraction is:

```text
(PSI_SMN2_correction - PSI_scramble) /
(PSI_control - PSI_scramble)
```

Report it with the underlying PSI values and confidence interval. Do not cap
the value at 0 or 1 because overshoot and movement away from control are
scientifically informative.

Power calculations should use pilot donor-to-donor variance from the validated
assay, not the 31 source libraries as if they were independent patients.

## 9. Stage 4: Isoform Mechanism and Functional Rescue

Advance an event only after raw-read and junction-assay confirmation.

### Causal perturbation design

Compare:

1. SMA plus negative control.
2. SMA plus SMN2 correction.
3. SMA plus isoform-specific perturbation.
4. SMA plus SMN2 correction and isoform-specific perturbation.
5. Matched control neurons.

The perturbation should alter the implicated isoform or splice junction while
minimizing changes in total gene abundance. Verify both isoform fraction and
total transcript before interpreting phenotype. Where an isoform-specific
protein assay is unavailable, sequence confirmation and an orthogonal RNA
assay are required.

### Functional outcomes

- longitudinal motor-neuron survival;
- apoptosis or caspase activation;
- axon length, branching, degeneration, and growth-cone integrity;
- electrophysiological maturation;
- neuromuscular-junction formation or maintenance in a validated co-culture;
- event-relevant protein or pathway output; and
- SMN target engagement to distinguish adjunct rescue from failed SMN
  correction.

Functional image analysis should be blinded. The combination arm tests whether
an isoform effect adds to SMN restoration; it does not justify clinical use.

## 10. Decision Gates

| Gate | Continue when | Stop or redesign when |
| --- | --- | --- |
| Raw confirmation | Exact structures, orientation, and effects are recoverable with acceptable QC | Reference mismatch, ambiguous orientation, or systematic non-reproduction remains unexplained |
| Assay validation | Expected products are sequenced and quantification is precise and specific | Multiple unresolved products, genomic contamination, or poor dynamic range |
| Independent lines | Disease and correction directions generalize beyond one donor | Signal is confined to one line or tracks a technical batch |
| Isoform causality | Isoform perturbation changes a relevant phenotype without simply changing total gene abundance | Phenotype follows toxicity, total-expression collapse, or off-target effects |
| Combination | Isoform perturbation adds reproducible benefit beyond SMN correction | No incremental effect or unacceptable cellular toxicity |

Negative results are useful. A publishable resource should report how many of
the frozen 12 events pass each gate rather than replacing failed events with
new favorable ones.

## 11. Manuscript-Ready Outputs

The minimum defensible paper package is:

1. Complete sample-selection and reference provenance for all selected runs.
2. Raw QC and strandedness results for every selected library.
3. Fixed-event rMATS output with explicit group orientation; add discovery
   rMATS output before making any genome-wide raw-discovery claim.
4. Concordance of official and independently processed effects.
5. A second-engine junction-level sensitivity analysis.
6. Blinded orthogonal validation of the frozen 12-event panel.
7. Separate reporting of the three exploratory OMN bridge events.
8. Donor-level raw values, analysis code, failed assays, and exclusions.
9. Author clarification regarding the supplementary-workbook orientation
   before wording the discrepancy as an error.

## 12. Reproducible Files

- `r/11_rank_splicing_validation_events.R`
- `results/r/splicing/GSE290979_R_validation_event_ranking.tsv`
- `results/r/splicing/GSE290979_R_primary_validation_panel.tsv`
- `results/r/splicing/GSE290979_R_exploratory_OMN_bridge_panel.tsv`
- `results/r/splicing/GSE290979_R_validation_score_specification.tsv`
- `results/r/splicing/GSE290979_R_validation_ranking_summary.tsv`
- `results/r/figures/GSE290979_R_validation_priority_top20.png`
- `results/r/figures/GSE290979_R_validation_panel_components.png`

## Bottom Line

The next experiment is not to increase or inhibit all 12 genes. It is to prove
that the exact splice events are real, reproducible in independent human motor
neurons, and responsive to SMN restoration. Only then should isoform-specific
causal and functional experiments begin.
