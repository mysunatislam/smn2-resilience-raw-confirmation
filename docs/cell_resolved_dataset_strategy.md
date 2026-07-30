# Cell-Resolved Human Motor-Neuron Strategy

## Dataset Decision

No openly downloadable, independent human SMA motor-neuron single-cell or
single-nucleus disease cohort was found that can presently serve every desired
role. The defensible design separates three roles instead of calling one
dataset more independent than it is.

| Dataset | Role | Independent units | What it answers | What it cannot answer |
| --- | --- | ---: | --- | --- |
| [GSE290980](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290980) | Same-study SMA cell-resolved context | 4 donor lines: 2 control, 2 SMA | Which author-defined organoid cell types, including motor neurons, show SMA-associated expression changes? | It overlaps the GSE290979 study/model and is not independent validation. |
| [GSE243076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243076) | Independent adult human motor-neuron localization | 9 adult spinal-cord donors | Is a candidate detected or selectively localized in the author-defined `CHAT`/`SLC5A7` C20 motor-neuron cluster? | It has no SMA contrast and cannot validate disease direction. |
| [EGAD00001011259](https://ega-archive.org/datasets/EGAD00001011259) | Future independent SMA disease validation | 8 isogenic differentiation samples across 3 times | Independent time-resolved human SMA/control single-cell validation. | EGA data are controlled access and approximately 954 GB; DAC approval is required. |

Peripheral-blood and CSF single-cell SMA datasets were rejected for the
motor-neuron question because they contain immune cells rather than spinal
motor neurons. GSE290980 is retained because it directly resolves the bulk
organoid cell-mixture limitation, but every output labels it as same-study.

## GSE243076 Raw Audit and Adult Localization

The official 798 MB raw archive was downloaded and SHA-256 locked. The R
workflow reads all nine 10x matrices sequentially with `Matrix`, avoiding an
all-donor object that would exceed available workstation memory.

Observed raw dimensions:

- 9 donor matrices;
- 36,601 features in every matrix;
- 74,711 deposited filtered nuclei;
- 72,834 nuclei pass a threshold-only proxy using the published gene, UMI,
  complexity, and mitochondrial cutoffs;
- the article reports 64,021 final nuclei after its complete QC and doublet
  removal workflow.

The threshold-only audit is not claimed as an exact re-clustering. Doublet
calls and all downstream author exclusions are not recoverable from simple
thresholding alone.

Primary cell localization uses the authors' complete 36,601-gene average-
expression matrix for neuronal clusters C0-C20. C20 identity is supported by
`CHAT`, `SLC5A7`, `MNX1`, and `ISL1`, each maximized and strongly enriched in
C20. Candidate localization is descriptive; no cell-level p-value is computed
and cells are never treated as independent donors.

## GSE290980 Same-Study Disease Context

Official GEO metadata confirms eight libraries from four donor lines. The
analysis imports Supplementary Data 4 and 5 from the linked article, verifies
their SHA-256 values, reconstructs all 19 cluster annotations, and parses the
authors' cell-type pseudobulk results.

Key source results reproduced structurally:

- 5,700 marker rows: top 300 genes for each of 19 clusters;
- four author-defined motor-neuron clusters: 6, 7, 10, and 18;
- 2,031 significant motor-neuron pseudobulk DEGs;
- additional pseudobulk DEG lists for pFP, progenitor, astroglia, V2b, V2a,
  and pMN populations.

The workbook defines positive log2 fold change as control greater than SMA.
The R output stores both that source value and a clearly named
`sma_vs_control_log2_effect` with the sign reversed.

## Cross-Model Result

Among the 37 exploratory full-direction cross-model candidates:

- 5 occur in the GSE290980 motor-neuron DEG list;
- `PFKFB2` and `SHQ1` show the expected SMA-opposed direction;
- 35 are detected in adult C20 motor neurons;
- none meets the stringent independent adult C20-selective localization rule.

This is a useful filter, not a failed analysis. The result does not support a
claim that the current 37-gene set is generally motor-neuron-specific.
`PFKFB2` and `SHQ1` receive same-study motor-neuron disease-direction support,
but neither is C20-selective in the adult atlas.

## Power and Inference Rules

- The unit of replication is donor line or donor, never cell count.
- GSE290980 has only two donor lines per genotype and is contextual.
- GSE243076 contributes nine independent donors for localization, but no SMA
  disease contrast.
- Thousands of cells improve cell-state resolution; they do not turn four
  donor lines into thousands of biological replicates.
- Candidate conclusions remain exploratory and require independent donor or
  experimental validation.

## Reproduce

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\14_gse243076_independent_motor_neuron_atlas.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\15_gse290980_cell_resolved_sma.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\16_integrate_cell_resolved_context.R
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\06_validate_outputs.R
```

Primary outputs are under `results/r/cell_resolved/`; figures are under
`results/r/figures/`. Exact source URLs, local paths, and SHA-256 values are
frozen in [`config/cell_resolved_resources.tsv`](../config/cell_resolved_resources.tsv).
