# Prospective Validation Freeze

Freeze timestamp: `2026-08-01T02:17:48+06:00`

Status: `FROZEN_BEFORE_RAW_STAR_RMATS_OR_NEW_WET_LAB_RESULTS`

This declaration fixes the candidate genes, evidence tiers, splice-event
targets, validation order, and success criteria before raw splice-junction or
new experimental results are examined. The dated manifests must not be
overwritten. Any change requires a new dated version that explains the reason
without deleting this version.

## Locked Files

| File | Rows | SHA-256 |
|---|---:|---|
| `config/frozen_37_gene_shortlist_2026-08-01.tsv` | 37 | `9c67056b1d45ec0a519038359d6fdd284396ba48821b94708b9ba5c11434bcb2` |
| `config/frozen_83_splice_events_2026-08-01.tsv` | 83 | `1ab84f7719677532cb59c08318fb29650d7268ef00f81d0be6e6f2d36a6460ec` |
| `config/frozen_12_primary_splice_events_2026-08-01.tsv` | 12 | `4b0e27166020feccdc22e8f9e10290be81496f9a24e804cf6913b2507ae7a34c` |

The standard checksum list is
`config/validation_freeze_2026-08-01.sha256`.

## Frozen 37-Gene Shortlist

The order and publication evidence tiers are fixed.

| Tier | Genes |
|---|---|
| Tier 1: raw and processed unit robust | `LY6H`, `HS3ST5`, `ZNF853`, `IL17D` |
| Tier 2: raw unit robust | `LSAMP`, `FAAH`, `C12orf60`, `KCNJ4`, `CABLES2` |
| Tier 3: processed unit robust | `PNCK`, `CLPTM1`, `PDPR` |
| Tier 4: raw direction only | `HEY1`, `NRN1`, `KCNIP2`, `ROBO2`, `UGT8`, `CTBP1-DT`, `RAB28`, `SHQ1`, `SPIN1`, `SLC38A9` |
| Tier 5: other cross-model | `ADCYAP1R1`, `HLA-DPB1`, `ZMAT4`, `LINC00665`, `FCSK`, `LRRC4`, `PFKFB2`, `ADHFE1`, `ZNF615`, `TTTY15`, `ATP9B`, `ZNF204P`, `EIF4ENIF1`, `ZNF827`, `RIOX1` |

The four primary expression candidates are fixed as `LY6H`, `HS3ST5`,
`ZNF853`, and `IL17D`. Later splice-junction, permutation, bootstrap, or
experimental results cannot change their Tier 1 designation or the 37-gene
ranking.

## Frozen Splice Events

All 83 events in the dated manifest are the fixed secondary confirmation set.
No event may be added to or removed from that set after raw results are seen.
The following 12 events are the fixed primary panel.

| Order | Gene | Class | Frozen event key | Processed disease dPSI | Processed treatment dPSI |
|---:|---|---|---|---:|---:|
| 1 | MRFAP1 | RI | `RI|MRFAP1|chr4|+|6641231|6641610|6641231|6641258|6641445|6641610` | -0.237 | 0.184 |
| 2 | IFI27L1 | ES | `ES|IFI27L1|chr14|+|94101813|94101975|94096840|94096965|94102476|94102714` | -0.318 | 0.196 |
| 3 | PDE9A | ES | `ES|PDE9A|chr21|+|42698912|42699011|42653751|42653883|42733355|42733426` | -0.247 | 0.277 |
| 4 | HERPUD1 | A5SS | `A5SS|HERPUD1|chr16|+|56935234|56935497|56935234|56935312|56936686|56936817` | 0.449 | -0.177 |
| 5 | ABCA1 | RI | `RI|ABCA1|chr9|-|104799571|104799988|104799571|104799733|104799818|104799988` | -0.375 | 0.309 |
| 6 | SAT2 | ES | `ES|SAT2|chr17|-|7626761|7626793|7626237|7626614|7627142|7627226` | -0.502 | 0.323 |
| 7 | TNPO2 | A3SS | `A3SS|TNPO2|chr19|-|12705681|12706109|12705681|12705830|12706195|12706367` | -0.414 | 0.281 |
| 8 | SEPTIN11 | MXE | `MXE|SEPTIN11|chr4|+|76987808|76987870|76995780|76995938|76949713|76949930|76996424|76996539` | -0.362 | 0.200 |
| 9 | FAT3 | MXE | `MXE|FAT3|chr11|+|92559409|92559705|92696745|92696783|92524633|92524948|92697383|92697445` | 0.341 | -0.384 |
| 10 | MGEA5/OGA | MXE | `MXE|MGEA5|chr10|-|101803734|101804019|101807823|101807901|101798841|101799455|101817823|101818465` | 0.361 | -0.313 |
| 11 | COL5A2 | RI | `RI|COL5A2|chr2|-|189100106|189100418|189100106|189100139|189100333|189100418` | 0.258 | -0.162 |
| 12 | EIF5 | A5SS | `A5SS|EIF5|chr14|+|103336676|103337143|103336676|103336849|103338326|103338472` | -0.114 | 0.300 |

## Raw Splice-Junction Design

### Inputs

The only primary raw libraries are the nine checksum-verified paired FASTQ
runs already fixed in `config/GSE290979_local9_sample_sheet.tsv`:

`SRR32559040`, `SRR32559042`, `SRR32559038`, `SRR32559021`,
`SRR32559034`, `SRR32559019`, `SRR32559030`, `SRR32559016`, and
`SRR32559025`.

The disease comparison is the two untreated SMA donor lines S2 and S3 versus
the three untreated control donor lines C1, C2, and C3. The treatment
comparison is R6-MO versus scramble in S2 and S3. Line-specific treatment
effects are evaluated separately for S2 and S3. Libraries are biological
units; reads and junctions are not treated as independent replicates.

The locked reference is GENCODE v47, GRCh38 primary assembly:

- `GRCh38.primary_assembly.genome.fa`
- `gencode.v47.primary_assembly.annotation.gtf`
- read length 151 and STAR `sjdbOverhang=150`
- STAR 2.7.10a
- rMATS turbo 4.3.0
- stranded paired-end setting consistent with the audited library orientation

Reference files, software versions, parameters, and checksums must be written
to the final provenance table. FASTQ or BAM deletion is not permitted.

### Event Matching

Primary structural recovery requires an exact match on event class,
chromosome, strand, and every event-defining exon or splice-site coordinate
encoded in the frozen event key. There is no coordinate tolerance in the
primary analysis. A near-coordinate or gene-only match is exploratory and
must not be counted as recovered.

If multiple raw events match one frozen event, the exact unique match is used.
Ambiguous exact matches are reported as ambiguous and not counted as primary
recovery. Every unrecovered event receives exactly one primary reason:

1. no exact structural match;
2. annotation or chromosome mismatch;
3. ambiguous exact match;
4. exact structure with no finite PSI;
5. insufficient junction support;
6. sample-level alignment or QC failure.

### Event-Level Criteria

For each event, report the raw disease and treatment dPSI, rMATS nominal
p-value and genome-wide FDR, BH-adjusted q-value within the frozen 83-event
set, group junction counts, line-specific PSI, and the following fixed flags.

- `structurally_recovered`: exact coordinate match and finite PSI.
- `adequate_junction_support`: median inclusion-plus-skipping informative
  junction count is at least 10 in every pooled comparison group.
- `disease_direction_reproduced`: structurally recovered, raw disease dPSI
  has the same sign as the frozen disease dPSI, and absolute raw disease dPSI
  is at least 0.05.
- `treatment_reversal_reproduced`: raw treatment dPSI has the same sign as
  the frozen treatment dPSI, has the opposite sign to raw disease dPSI, and
  absolute raw treatment dPSI is at least 0.05.
- `both_lines_corrected`: in both S2 and S3, R6-MO changes PSI opposite to the
  raw disease direction and reduces distance to the untreated control mean by
  at least 0.05 relative to scramble.
- `strong_raw_confirmation`: all five flags above are true and the
  within-83 BH q-values are below 0.05 for both disease and treatment.

Direction-only and low-support results remain reportable but cannot be called
strong raw confirmation.

### Panel-Level Criteria

The primary 12-event raw validation aim is considered successful only if all
12 are attempted and:

1. at least 8 of 12 are structurally recovered;
2. at least 6 of 12 reproduce disease direction;
3. at least 4 of 12 reproduce treatment reversal and correction in both
   donor lines with adequate junction support; and
4. among at least 30 structurally recovered events from the full 83-event
   set, raw versus processed disease dPSI has Spearman rho at least 0.50.

All counts, the full raw-versus-processed correlation, junction-support
table, one sashimi plot per primary event, and a reason for every unrecovered
event are mandatory regardless of whether these thresholds are met.

## Wet-Lab Validation Freeze

### RT-PCR Panel

The first RT-PCR panel is fixed before primer testing:

1. `MRFAP1` RI
2. `HERPUD1` A5SS
3. `COL5A2` RI
4. `IFI27L1` ES
5. `PDE9A` ES
6. `TNPO2` A3SS

This set cannot be replaced based on STAR/rMATS results. A target may be
declared technically infeasible only after documenting primer sequences,
tested annealing conditions, expected product sizes, and failure evidence.
Such a target remains in the denominator.

An RT-PCR event is confirmed only when:

1. both expected isoforms are resolved and at least one representative
   amplicon identity is confirmed by Sanger sequencing;
2. at least three independent differentiations are analyzed per condition;
3. disease PSI or isoform-ratio change has the frozen direction and an
   absolute PSI change of at least 0.05;
4. R6-MO change is opposite to disease with an absolute PSI change of at
   least 0.05;
5. S2 and S3 each move in the correction direction; and
6. the biological-unit model gives two-sided p below 0.05 and BH q below
   0.10 across the six fixed events.

The RT-PCR panel-level aim succeeds if at least three of six events meet all
six criteria.

### Functional Candidates

The functional order is fixed as `LY6H` first and `IL17D` second. The planned
direction is restoration in SMA motor neurons, with matched vector or
non-targeting controls. The primary phenotype is motor-neuron survival; the
first secondary phenotype is neurite or axon length. A single assay day must
be fixed before the first experiment and used for all groups.

The biological unit is an independent differentiation or donor-line
experiment, not an image field, well subdivision, neurite, or individual
cell. At least three independent experiments are required.

A candidate-only functional rescue requires improvement over SMA control with
two-sided p below 0.05, absolute survival improvement of at least 10
percentage points or standardized effect size at least 0.5, and concordant
improvement in the prespecified secondary phenotype. A claim of benefit
beyond SMN restoration additionally requires candidate plus SMN restoration
to outperform SMN restoration alone with p below 0.05 and standardized effect
size at least 0.5. Failure to meet the additional criterion must be reported
as no evidence of benefit beyond SMN restoration.

## Confirmatory and Exploratory Analyses

Confirmatory:

- exact raw matching and fixed metrics for all 83 frozen splice events;
- the 12-event primary raw panel and its panel-level thresholds;
- the six-event RT-PCR panel;
- functional tests of LY6H followed by IL17D;
- permutation tests of the observed 37-gene and four-gene counts;
- biological-unit bootstrap stability and the prespecified score-sensitivity
  schemes.

Exploratory:

- any de novo rMATS event outside the frozen 83;
- near-coordinate, gene-only, or alternative-annotation event matches;
- any replacement RT-PCR event or functional candidate;
- pathway analyses, new composite scores, machine learning, or additional
  public datasets prompted by the validation results;
- alternative thresholds not listed in this declaration.

Exploratory findings may be reported, but they cannot replace or redefine a
failed confirmatory target.
