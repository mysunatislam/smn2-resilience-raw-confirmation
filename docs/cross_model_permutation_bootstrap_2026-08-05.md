# Cross-model permutation and bootstrap analysis

**Frozen analysis date:** 2026-08-05

## Confirmatory questions

1. Is the frozen count of 37 complete-pattern candidates larger than expected
   after breaking gene identity across the disease, treatment, and external
   evidence components?
2. Is the observed overlap of four processed- and raw-unit-robust candidates
   larger than expected when raw robustness labels are exchangeable among the
   37 frozen candidates?
3. Do the frozen candidates retain their direction pattern when biological
   donors or donor lines are resampled?

The candidate set, publication tiers, 12-event splicing panel, and success
criteria were not changed after these analyses.

## Permutation design

- Eligible universe: 11,326 approved genes.
- Iterations: 10,000; seed: 290979.
- Candidate-count null: independently permute the GSE290979 disease,
  GSE290979 treatment, and GSE108094 disease effects across genes. Keep the
  GSE93939 effect and frozen exploratory natural-resistance support fixed.
- Robust-overlap null: permute the nine raw robustness labels among the 37
  candidates while keeping the seven processed robustness labels fixed.
- Empirical probability: `(1 + exceedances) / (1 + iterations)`.

The 37-gene count was not enriched (`p = 0.8357164`). The overlap of four
raw- and processed-unit-robust candidates had empirical `p = 0.0440956`.

## Biological-unit bootstrap

- Iterations: 1,000; seed: 939390.
- GSE93939: resample 19 donors with replacement and fit a
  covariate-adjusted weighted model to TMM logCPM.
- GSE290979 disease: resample three control and two SMA donor lines with
  replacement.
- GSE290979 treatment: resample the two paired SMA treatment lines with
  replacement.
- GSE108094: retain the external effect as fixed because unit-level values
  were unavailable.
- Random technical-library splitting: not used.

Tier 1 selection frequencies were 0.898 (`LY6H`), 0.959 (`HS3ST5`), 0.949
(`ZNF853`), and 0.941 (`IL17D`).

## Score sensitivity

Nine variants tested equal weights, double weight for an individual component,
omission of each component, and median-percentile aggregation. Genome-wide
rank correlations with the equal-weight score ranged from 0.815 to 0.956.

## Claim boundary

The directional shortlist size is not permutation-enriched. The stronger
result is the reproducibility of the four-gene raw/processed robustness
overlap and the high biological-unit bootstrap frequency of those four genes.
Bootstrap results are sensitivity evidence because they use TMM-logCPM models
and hold GSE108094 fixed.
