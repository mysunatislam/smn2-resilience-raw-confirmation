# Cross-Model Human Motor-Neuron Resilience: R Analysis and Dataset Selection

**Analysis date:** 22 July 2026
**Primary implementation:** R 4.6.0, edgeR 4.10.1, limma 3.68.4
**Status:** Reproducible exploratory human analysis; suitable for hypothesis
generation, not clinical use

## Executive Decision

The project uses two core GEO datasets, one primary external directional-
validation dataset, and one constrained external sensitivity dataset:

1. **GSE93939** as the natural-resistance anchor based on human
   oculomotor-versus-spinal motor-neuron expression differences.
2. **GSE290979** as the human SMA disease, SMN-directed response, and
   splice-restoration model.
3. **GSE108094** as an independent iPSC-SMA spinal motor-neuron directional
   validation set.
4. **GSE69175** as purified patient-derived motor-neuron directional
   sensitivity evidence, never as discovery.

They were not selected merely because each has more than 20 GEO entries. They
were selected because their biological roles are complementary, raw counts and
sample metadata are available, and each directly addresses one side of the
research question.

The important correction is that library count is not donor count:

| Dataset | GEO libraries | Primary independent units | Use |
| --- | ---: | ---: | --- |
| GSE93939 | 39 | 19 donors in OMN-versus-spinal analysis | Natural-resistance anchor |
| GSE290979 | 31 | 5 iPSC lines overall; 2 SMA lines for treatment | SMA disease, SMN response, and splice restoration |
| GSE108094 | 8 | 4 iPSC motor-neuron lines | External directional validation only |
| GSE69175 | 4 | 2 lines total: 1 control and 1 SMA | External directional sensitivity only |

The 20+ library requirement applies to the two core datasets. Neither external
dataset is promoted to discovery status or used to inflate the core sample
count.

The R analysis is now authoritative for statistical inference. The earlier
Python analysis remains useful for data preparation, provenance, and an
independent event-level implementation, but its 451-gene OMN result and
Tier-A/Tier-B expression lists should not be cited as confirmed findings.

## 1. Dataset Selection Criteria

A core dataset had to satisfy all of the following:

1. **Human origin:** `Homo sapiens` for primary inference.
2. **Direct biological relevance:** motor-neuron subtype resistance or SMA with
   an SMN-restoring intervention.
3. **At least 20 relevant sequencing libraries:** not a total created by mixing
   species, unrelated assays, or cell models.
4. **Usable raw counts:** count matrices rather than only normalized fold
   changes.
5. **Traceable sample metadata:** group, donor or line, treatment, and technical
   variables must be recoverable.
6. **Biological replication:** donor and cell-line counts must be reported
   separately from replicate libraries.
7. **Complementarity:** discovery and validation cohorts should answer
   different parts of the hypothesis rather than duplicate the same weak
   contrast.

The full screening table is in `config/dataset_selection.tsv`.

## 2. Why GSE93939 Was Selected

[GSE93939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93939)
contains 39 laser-capture RNA-seq libraries from post-mortem human motor
neurons:

- 20 oculomotor-neuron libraries.
- 12 cervical or lumbar spinal motor-neuron libraries.
- 7 Onuf motor-neuron libraries.

### Strengths

- It directly compares a relatively resistant human motor-neuron population
  with vulnerable spinal motor neurons.
- The cells were anatomically identified and laser captured, giving more
  specific neuronal signal than whole spinal cord.
- Raw gene counts are available.
- The 32-library OMN-versus-spinal contrast contains 19 donors.
- Donor, sex, age, post-mortem delay, tissue source, and platform metadata can
  be modeled.

### Limitations

- It contains control tissue, not SMA tissue.
- The original biological context includes ALS resistance, so transfer to SMA
  is a hypothesis.
- Sex and tissue bank are strongly imbalanced between OMN and spinal groups.
- Two sequencing platforms were used.
- Several donors contributed repeated capture libraries.

These limitations are why the R model includes donor blocking and covariates,
and why a HiSeq 2000-only sensitivity analysis is mandatory.

## 3. Why GSE290979 Was Selected

[GSE290979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290979)
contains 31 human spinal-cord-organoid RNA-seq libraries associated with the
[primary organoid study](https://www.nature.com/articles/s41467-025-67725-1):

- 9 untreated control libraries from three lines.
- 6 untreated SMA libraries from two lines.
- 8 scramble-treated SMA libraries from two lines.
- 8 R6-MO-treated SMA libraries from the same two lines.

### Strengths

- It is a direct human SMA type I model.
- It includes both untreated disease and an SMN2-directed intervention.
- Treatment and scramble conditions are represented in both SMA lines.
- Raw counts, paired-end raw reads, sample-level rMATS values, and orthogonal
  SMN full-length-to-delta7 assays are available.
- The intervention provides a way to distinguish disease association from
  SMN-responsive change.

### Limitations

- There are only three control and two SMA donor lines.
- The treatment contrast has only two biological SMA pairs.
- Organoid libraries contain mixed cell populations.
- Individual organoid preparations are not independent patients.
- SMN1 and SMN2 are difficult to separate from ordinary short-read alignment.

For this reason, the R expression analysis sums replicate libraries within each
donor-line condition before primary inference.

## 4. Why GSE108094 Was Added

[GSE108094](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108094)
contains eight deeply sequenced human iPSC-derived spinal motor-neuron
libraries from two control and two SMA lines, each represented by two
libraries. It provides an independent disease direction in the cell type of
interest and official expression and exon-skipping outputs.

It remains external validation because it has only four lines. Its processed
Cuffdiff result reports four libraries per group but does not provide a sample-
level count matrix for a line-aware model; two libraries share each line.
Source p/FDR values are therefore context, while the integration uses effect
direction. R verifies that the deposited fold change is control minus SMA and
reverses it to the project-wide SMA-minus-control convention. Its hg19/Ensembl
75 splice coordinates are not directly matched to GSE290979 hg38 structures.

## 5. Additional Validation and Secondary Datasets

### GSE69175

[GSE69175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69175)
profiles purified patient-derived SMA motor neurons and is now included as a
checksum-locked external directional sensitivity analysis. It contains two
libraries from one `BJ-riPS` control line and two from one type I SMA `1-38G`
line. R verifies that the deposited Cuffdiff effect is SMA minus control and
recovers 5,093 approved direction-usable genes.

This is not an independent line-replicated contrast: genotype and line are
perfectly confounded, and omitting either line removes one genotype. GSE69175
therefore cannot support LOLO and does not enter the discovery score. Its
published ER-stress biology and gene directions are retained as a stringent,
separately labeled sensitivity tier.

### GSE87281

[GSE87281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87281)
supports SMN-dependent intron-retention biology. Although the complete series
has 101 samples, it mixes mouse spinal cord, SH-SY5Y cells, human iPSC motor
neurons, RNA-seq, and DRIP-seq. The relevant human iPSC motor-neuron RNA-seq arm
has only seven libraries. Counting all 101 as one human motor-neuron cohort
would be invalid.

### GSE136719

[GSE136719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136719)
contains more than 880,000 human spinal-cord cells and nuclei. It is valuable
for cell-type and developmental annotation, but individual cells are not
independent donors and the study is not an SMA or OMN comparison. It is a
reference atlas, not a disease-discovery cohort.

### Mouse datasets

GSE115706 and GSE208629 remain useful for conservation and mechanism after the
human hypothesis set is fixed. They are not used to satisfy the human 20+
requirement.

## 6. R Analysis of GSE93939

### Model

The primary R pipeline performs:

1. CPM filtering at CPM >= 1 in at least seven of the 32 primary libraries.
2. TMM library normalization with edgeR.
3. Limma-voom mean-variance modeling.
4. Sample quality weights.
5. Donor blocking with `duplicateCorrelation`.
6. Fixed adjustment for sex, age, post-mortem delay, tissue source, and
   sequencing platform.
7. Robust empirical-Bayes moderation.
8. Benjamini-Hochberg FDR correction.

The same model is repeated in the 25 HiSeq 2000 libraries from 13 donors. The
`source_NIH` term is absent from this sensitivity model because no NIH sample
remains in that subset, making the term non-estimable.

A second sensitivity analysis sums HiSeq 2000 libraries to 13 donor
pseudobulks and fits an edgeR quasi-likelihood model.

### Results

| Analysis | Libraries or pseudobulks | Donors | Genes tested | FDR < 0.05 |
| --- | ---: | ---: | ---: | ---: |
| Adjusted blocked limma-voom | 32 libraries | 19 | 16,443 | 2 |
| HiSeq 2000 blocked sensitivity | 25 libraries | 13 | 16,443 | 1 |
| HiSeq 2000 donor pseudobulk | 13 pseudobulks | 13 | 15,435 | 0 |

The two full-model FDR genes are **HOXC8** and **TMEM233**, both lower in OMN
than spinal motor neurons. Therefore:

- **Zero OMN-positive genes reach FDR 0.05 in the primary R model.**
- Primary and same-platform effects remain directionally similar, with Pearson
  `r = 0.854`.
- The same-platform blocked and donor-pseudobulk effects correlate at
  `r = 0.803`.

The result does not mean that OMN biology is absent. It means that this dataset
does not support a large, covariate-adjusted list of individual OMN-enriched
genes at controlled FDR.

## 7. Why the Earlier Python Hit Count Changed

R and Python effect estimates correlate strongly (`r = 0.949` across 16,443
genes), so they broadly agree about direction and magnitude. They disagree
about uncertainty.

The R workflow adds RNA-seq-specific and repeated-sample modeling:

- TMM composition normalization.
- Voom mean-variance weights.
- Sample quality weights.
- A donor correlation model.
- Limma empirical-Bayes variance moderation.

The earlier Python model used normalized log expression and custom
cluster-robust OLS. Its 451 FDR genes and downstream 35 Tier-A/44 Tier-B lists
are therefore superseded for inferential claims. They may still be inspected
as exploratory effect rankings, but not presented as statistically validated
human resistance genes.

## 8. R Analysis of GSE290979 Expression

### Disease model

The 15 untreated libraries are summed into five donor-line pseudobulks:

- C1, C2, and C3 controls.
- S2 and S3 SMA.

An edgeR quasi-likelihood model compares the three control lines with the two
SMA lines. Genotype and donor line are perfectly confounded, so the individual
libraries cannot create more independent disease donors.

### Treatment model

The 16 treatment libraries are summed within line and treatment, producing:

- S2 scramble and S2 R6-MO.
- S3 scramble and S3 R6-MO.

The model is paired by line. It has only one residual degree of freedom and is
therefore used primarily for effect direction. A 16-library fixed-line model is
retained as a preparation-level sensitivity analysis, not donor replication.

### Results

| Contrast | Modeled samples | Independent lines | FDR < 0.05 |
| --- | ---: | ---: | ---: |
| SMA minus control | 5 pseudobulks | 5 | 0 |
| R6-MO minus scramble | 4 paired pseudobulks | 2 | 0 |
| R6-MO library-level sensitivity | 16 libraries | 2 | 0 |

Pseudobulk and library-level treatment effects correlate at `r = 0.986`.
This supports stable direction estimates but not gene-level significance.

## 9. R Splicing Audit

The R implementation independently reads the official rMATS workbooks,
normalizes Excel-corrupted gene symbols, constructs full structural event
keys, verifies group means, reverses the deposited treatment sign to biological
R6-MO minus scramble, and applies the two-line consistency filter.

It exactly reproduces the independent Python event set:

| Event class | Exact common | Direction reversal | Pooled corrected | Strict in S2 and S3 |
| --- | ---: | ---: | ---: | ---: |
| Exon skipping | 80 | 48 | 46 | 35 |
| Alternative 5' site | 15 | 8 | 8 | 6 |
| Alternative 3' site | 15 | 7 | 7 | 5 |
| Mutually exclusive exon | 36 | 28 | 27 | 22 |
| Retained intron | 36 | 21 | 20 | 15 |
| **Total** | **182** | **112** | **108** | **83** |

The 83 strict events represent 74 genes. R and Python match all 83 event-type,
disease-event-ID, and treatment-event-ID triples exactly.

The R audit also confirms:

- The named treatment PSI columns place scramble in rMATS group 1 and R6-MO in
  group 2.
- Biological treatment effect is `R6-MO - scramble`, the negative of the
  deposited rMATS difference.
- All 111 published rescue pairs move in the same biological direction as the
  disease effect under this orientation.
- Only 67 published pairs have identical full event structures.

This remains an audit of official processed tables. Raw FASTQ re-alignment is
required before presenting the sign observation as a definitive correction.

## 10. SMN Target Engagement

R independently extracts the source-data full-length-to-SMN-delta7 ratios:

| Line | Scramble mean | R6-MO mean | R6 / scramble |
| --- | ---: | ---: | ---: |
| S1 | 6.940 | 229.540 | 33.08 |
| S2 | 6.724 | 53.938 | 8.02 |
| S3 | 7.398 | 27.404 | 3.70 |

This strongly supports molecular activity of R6-MO even though donor-level
gene-expression tests have low power.

## 11. R-Defined Candidate Interpretation

No strict splicing gene also has robust OMN-positive expression at FDR 0.05.
Three genes have nominal primary OMN `p < 0.05`, positive same-platform effects,
and same-platform `p < 0.10`:

| Gene | Strict event | R OMN effect | Primary p | Primary FDR | Same-platform p |
| --- | --- | ---: | ---: | ---: | ---: |
| KCNAB3 | ES | 6.058 | 0.0078 | 0.667 | 0.0353 |
| LINC00665 | ES | 3.732 | 0.0243 | 0.763 | 0.0943 |
| TSPOAP1 | MXE | 3.530 | 0.0386 | 0.817 | 0.0635 |

These are **exploratory bridge candidates**, not validated genes. Their high
FDR values must be shown whenever they are discussed.

The most defensible candidate resource is the 74-gene strict splicing set,
prioritized for junction-level experimental validation without assuming that
total gene activation or inhibition will be beneficial.

The frozen validation design ranks all 83 events with a transparent score for
effect, source FDR, S2/S3 consistency, junction support, and assay feasibility.
It selects 12 unique genes while representing all five event classes. The
primary panel, separate exploratory natural-resistance bridge, and donor-aware wet-lab
protocol are documented in
[`splicing_validation_panel.md`](splicing_validation_panel.md).

## 12. Cross-Model Human Motor-Neuron Resilience

The cross-model analysis does not discard GSE93939 or require it to produce a
large FDR-significant list. It uses GSE93939 as a directional natural-
resistance prior, then asks whether the opposite disease direction and aligned
SMN-response direction recur in GSE290979 and GSE108094.

Among 11,326 approved genes with usable effects in all model arms, 626 have the
complete directional-resilience pattern. Thirty-seven also meet the frozen
exploratory GSE93939 criterion; none is called robust because the GSE93939 FDR
values remain high. LINC00665 is the only member of this 37-gene set that also
has a strict GSE290979 splice-restoration event. FAT3 and TNPO2 are the two
members of the frozen 12-event assay panel with all four total-expression
directions.

Genome-wide effect correlations are weak, with absolute Spearman rho below
0.12 in every comparison. The models therefore provide focused intersection,
not a universal shared transcriptome. Full methods, correlations, and external
validation limitations are in
[`cross_model_human_motor_neuron_resilience.md`](cross_model_human_motor_neuron_resilience.md).

## 13. Cell-Resolved Motor-Neuron Context

The bulk-organoid limitation is addressed with two explicitly different human
cell-resolved roles.

`GSE290980` is the single-cell arm linked to the GSE290979 organoid study. It
contains eight libraries from two control and two SMA donor lines and 27,114
reported cells. The R audit reconstructs 19 clusters, four motor-neuron
clusters, and 2,031 author-pseudobulk motor-neuron DEGs. Because it is the same
study/model, it supplies cell-type disease localization but not independent
validation.

`GSE243076` is an independent adult human spinal-cord snRNA-seq atlas from nine
donors. All nine deposited 10x matrices were audited sequentially: they contain
74,711 filtered nuclei and 36,601 features. A threshold-only QC proxy retains
72,834 nuclei; the paper's final 64,021 also reflects doublet removal and
downstream exclusions. The author-defined C20 cluster is validated by strong
`CHAT`, `SLC5A7`, `MNX1`, and `ISL1` localization.

Of the 37 exploratory cross-model candidates, five occur in the GSE290980
motor-neuron DEG list and two (`PFKFB2`, `SHQ1`) have the expected SMA-opposed
direction. Thirty-five are detected in adult C20, but none is C20-selective
under the stringent localization rule. This refines the candidate set; it does
not validate a general motor-neuron-specific resilience program.

The complete rationale, rejected alternatives, controlled-access future
dataset, and power guardrails are in
[`cell_resolved_dataset_strategy.md`](cell_resolved_dataset_strategy.md).

## 14. Biological-Unit and Dataset Holdouts

No random train/test split is used anywhere in the analytical workflow.
Libraries from one donor or iPSC line are never divided across folds. The
deterministic robustness analysis performs:

- 19 GSE93939 primary leave-one-donor-out refits;
- 13 GSE93939 HiSeq 2000-only donor refits;
- five GSE290979 disease leave-one-line-out refits;
- separate S2 and S3 treatment effects because two paired lines cannot support
  a formal LOLO refit; and
- three leave-one-dataset-out ranking folds, omitting GSE93939, GSE290979, or
  GSE108094 as complete evidence sources.

The GSE93939 jackknife fixes donor correlation to the corresponding complete-
model consensus value while re-estimating TMM, voom weights, covariate effects,
and empirical Bayes in each fold. This is a deterministic coefficient-
stability analysis, not predictive accuracy.

Of the 37 exploratory candidates, 35 retain the GSE93939 direction in every
primary donor fold, 15 retain the GSE290979 disease direction in every line
fold, and 14 have the same treatment direction in both SMA lines. Seven pass
all three estimable checks: `LY6H`, `HS3ST5`, `ZNF853`, `PNCK`, `IL17D`,
`CLPTM1`, and `PDPR`. Six different candidates have GSE69175 directional
support, and no candidate passes both evidence sets. Full methods and outputs
are in
[`biological_unit_holdout_validation.md`](biological_unit_holdout_validation.md).

## 15. Recommended Publication Direction

The strongest current manuscript is an **orientation-aware cross-model human
motor-neuron resilience and SMA splicing analysis**, not a broad claim that
dozens of OMN resistance genes have been validated.

The all-nine local9 Salmon raw-read expression analysis is complete. Remaining
claim-dependent next steps are:

1. Re-align all 31 GSE290979 raw runs and run full rMATS discovery if the
   manuscript claims genome-wide workbook-independent splicing inference.
2. Confirm events with MAJIQ, LeafCutter, or another independent engine.
3. Quantify SMN1/SMN2 exon 7 with paralog-aware handling.
4. Validate selected events by junction RT-PCR in additional SMA donor lines.
5. Measure protein or functional consequences for the affected isoforms.
6. Seek clarification from the source authors before describing the workbook
   orientation as an article error.

The OMN comparison can still motivate experiments, but its gene-level signals
should remain exploratory until an independent human OMN-versus-spinal cohort
is available.

## 16. Reproduce the R Analysis

From the project root:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\00_install_packages.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\01_dataset_selection_audit.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\02_gse93939_limma_voom.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\02b_gse93939_donor_pseudobulk.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\03_gse290979_edger_pseudobulk.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\04_gse290979_splicing_audit.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\05_integrate_and_compare.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\07_prepare_gse290979_raw_confirmation.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\08_build_gse290979_fixed_event_set.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\11_rank_splicing_validation_events.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\12_prepare_gse108094_external_validation.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\13_cross_model_human_motor_neuron_resilience.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\14_gse243076_independent_motor_neuron_atlas.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\15_gse290980_cell_resolved_sma.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\16_integrate_cell_resolved_context.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\17_prepare_gse69175_external_validation.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\19_leave_one_unit_out_robustness.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\20_integrate_holdout_validation.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\27_import_gse290979_local9_salmon.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\28_integrate_local9_salmon_publication.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\06_validate_outputs.R
```

The final command performs hard output checks for every model arm, source hash,
holdout fold count, candidate identity, and the absence of random analytical
splits. Package versions and full R session information are written under
`results/r/`.

The claim-dependent raw gate is implemented under `hpc/gse290979_raw/`. Both
profiles bootstrap checksum-locked GENCODE v47 references, lock STAR 2.7.10a
and rMATS 4.3.0, validate ENA MD5 values, empirically gate stranded orientation,
and run raw STAR-count edgeR models. `validation21` processes 193.116 GiB and
performs only pre-specified targeted recount; `full31` processes 285.3 GiB and
adds reference-derived discovery rMATS. Neither Slurm profile has been
executed. The separate workstation `local9` profile has completed all-nine
Salmon transcriptome re-quantification and donor-aware gene-expression
analysis, but it does not substitute for either STAR/rMATS splice profile. See
[`gse290979_local9_protocol.md`](gse290979_local9_protocol.md).

## Bottom Line

The careful R analysis changes the interpretation in a healthy way:

- The large OMN expression list is not statistically robust after appropriate
  RNA-seq and donor modeling.
- GSE93939 remains informative as the natural-resistance anchor in a
  cross-model directional analysis.
- GSE290979 gene-level expression is underpowered at the donor-line level.
- The corrected splicing result is highly reproducible across two independent
  workbook-audit implementations, but raw junction-level reproduction is still
  required.
- The all-nine raw expression sensitivity analysis identifies four candidates
  robust in both raw and processed biological-unit analyses.
- Cell-resolved analysis supports `PFKFB2` and `SHQ1` directionally in the
  same-study motor-neuron arm; neither the 37-gene set nor those two genes show
  independent adult C20 selectivity.
- Biological-unit and whole-dataset omissions materially narrow candidate
  stability, and GSE69175 does not converge on the seven all-unit-stable genes.

That is a narrower conclusion, but it is much more defensible and gives the
project a clear next experiment.
