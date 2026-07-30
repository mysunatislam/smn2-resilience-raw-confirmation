# Biological-Unit and Dataset Holdout Validation

**Analysis date:** 22 July 2026
**Scope:** deterministic robustness analysis; no random train/test partitioning

## Decision

The analytical workflow contains no random train/test split. The only
`set.seed()` in the repository generates synthetic counts for a unit test.
Randomly dividing libraries would be invalid here because libraries from the
same donor or iPSC line could enter both partitions and create biological
leakage.

The project therefore uses three explicitly named schemes:

- **LODO-donor:** leave one biological donor and all of that donor's libraries
  out, then refit the model.
- **LOLO-line:** leave one iPSC donor line and all of that line's libraries out,
  then refit the model.
- **LODO-dataset:** omit one complete evidence dataset from candidate scoring
  and recompute ranks.

These are jackknife sensitivity analyses, not predictive-model accuracy
estimates. The complete policy, including non-estimable cases, is locked in
`config/resampling_policy.tsv`.

## GSE69175 Validation Role

[GSE69175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69175)
profiles purified human iPSC-derived motor neurons and is directly relevant to
SMA motor-neuron death and ER-stress biology. The public design contains:

- two control libraries from one `BJ-riPS` control line;
- two SMA libraries from one type I SMA `1-38G` line; and
- one deposited group-level Cuffdiff result rather than sample-level counts.

The checksum-locked audit recovers 19,286 expression rows, 5,394 rows marked
`OK`, 968 source FDR hits, and 5,093 approved genes usable for direction. The
project verifies from the FPKM values that the deposited effect is SMA minus
control.

GSE69175 is an **external directional sensitivity validation**, not an
independent discovery cohort. Genotype and line are perfectly confounded.
Formal LOLO is structurally impossible because omitting either line removes
one genotype, and the two libraries per line cannot be relabeled as two donor
lines. Source p-values and FDR values are therefore retained as study context,
not treated as line-aware replication.

## Biological-Unit Holdouts

### GSE93939 donor holdout

The primary 32-library OMN-versus-spinal model is refit 19 times, omitting all
libraries from one donor per fold. The 25-library HiSeq 2000 sensitivity model
is refit in 13 donor folds.

Each fold re-estimates TMM factors, voom mean-variance weights, covariate
coefficients, and robust empirical-Bayes moderation. The donor correlation is
fixed to the consensus value from the corresponding complete model. Fold-
specific `duplicateCorrelation()` and array-quality-weight optimization were
not used because they made the complete all-gene jackknife impractical on the
local system. This fixed-correlation sensitivity is labeled in every long
table. All primary folds retain the planned design; the HiSeq-only model drops
the absent `source_NIH` term, and one fold also drops `sex_M` after sex becomes
constant.

### GSE290979 line holdout

The five untreated donor-line pseudobulks are refit five times, omitting one
complete line per fold. All 17,014 modeled disease genes have five estimable
LOLO effects.

The R6-versus-scramble contrast has only two paired SMA lines. Omitting one
line leaves no biological replication for a formal edgeR refit. The workflow
therefore reports normalized S2 and S3 paired effects separately and requires
both line directions to agree with the complete paired model. This is called
**line-wise replication**, not formal LOLO.

### Other datasets

- GSE108094 needs raw sample-level reprocessing before four-line LOLO is
  possible; its public processed table contains group-level effects.
- GSE290980 needs a line-level cell-count matrix and rerun pseudobulk model.
- GSE243076 is a localization atlas without an SMA contrast, so disease LODO is
  not applicable.

## Leave-One-Dataset-Out Ranking

The 11,326-gene cross-model ranking is recomputed three times:

1. omit GSE93939;
2. omit GSE290979 as one source, removing both disease and treatment
   components; and
3. omit GSE108094.

GSE69175 remains a held-out sensitivity annotation and is not added as an
equally weighted score component. This prevents one control line and one SMA
line from influencing discovery ranks.

## Results for the 37 Exploratory Candidates

| Robustness result | Candidates |
| --- | ---: |
| GSE93939 primary donor direction stable in all 19 folds | 35 |
| GSE290979 disease direction stable in all 5 line folds | 15 |
| GSE290979 treatment direction shared by S2 and S3 | 14 |
| All three estimable biological-unit checks | 7 |
| Direction measurable in GSE69175 | 9 |
| GSE69175 disease direction opposes natural resistance | 6 |
| All biological-unit checks plus GSE69175 support | **0** |

The seven candidates passing all estimable donor/line direction checks are
`LY6H`, `HS3ST5`, `ZNF853`, `PNCK`, `IL17D`, `CLPTM1`, and `PDPR`.

The six candidates with GSE69175 directional support are `ROBO2`, `ATP9B`,
`HEY1`, `LRRC4`, `SPIN1`, and `EIF4ENIF1`.

The sets do not overlap. The correct interpretation is that the current public
data do not identify a fully convergent candidate. The holdouts narrow and
audit the hypotheses; they do not repair low donor-line power.

## Raw FASTQ Re-quantification Update

The pre-specified nine-library local9 profile subsequently completed
checksum-verified Salmon re-quantification for all five untreated donor lines
and both paired treatment lines. This is a technical re-analysis of the same
GSE290979 biological material, so it is not independent validation and is not
added to the frozen discovery score.

Among the 37 frozen candidates:

| Raw confirmation result | Candidates |
| --- | ---: |
| Mapped in both raw contrasts | 35 |
| Full raw disease-opposition and R6-alignment pattern | 21 |
| Full raw pattern plus five disease LOLO folds and both treatment lines | 9 |
| Robust in both raw and processed biological-unit analyses | 4 |
| Robust in either raw or processed biological-unit analysis | 12 |

The four candidates robust in both analyses are `LY6H`, `HS3ST5`, `ZNF853`,
and `IL17D`. Five additional candidates are raw-unit robust: `LSAMP`, `FAAH`,
`C12orf60`, `KCNJ4`, and `CABLES2`. Three retain processed-unit robustness
without meeting the complete raw unit criterion: `PNCK`, `CLPTM1`, and
`PDPR`.

These groups are publication evidence tiers, not a replacement ranking.
Direction stability is emphasized because none of the 37 candidates reaches
raw treatment FDR below 0.05, and only three have raw disease nominal
`p < 0.05`; none reaches raw disease FDR below 0.05.

## Reproduce

Run after the existing model, external-validation, and cell-resolved scripts:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\17_prepare_gse69175_external_validation.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\19_leave_one_unit_out_robustness.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\20_integrate_holdout_validation.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\27_import_gse290979_local9_salmon.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\28_integrate_local9_salmon_publication.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\06_validate_outputs.R
```

Primary outputs:

- `results/r/external_validation/GSE69175_expression_SMA_vs_control.tsv`
- `results/r/robustness/GSE93939_primary_LODO_summary.tsv`
- `results/r/robustness/GSE290979_disease_LOLO_summary.tsv`
- `results/r/robustness/GSE290979_treatment_linewise_summary.tsv`
- `results/r/robustness/human_cross_model_holdout_gene_table.tsv`
- `results/r/robustness/human_cross_model_holdout_candidates.tsv`
- `results/r/figures/human_holdout_robustness_top25.png`
- `results/r/publication/human_raw_confirmed_publication_candidates.tsv`
- `manuscript/supplementary_table_S1_candidates.tsv`
