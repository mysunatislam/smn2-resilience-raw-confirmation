# Cover Letter Draft

Status: `DRAFT_ONLY_NOT_FOR_SUBMISSION`

8 August 2026

Editors  
Journal of Neurology

Dear Editors,

We wish to submit the Original Communication, "Cross-model analysis
prioritizes reproducible human motor-neuron resilience candidates in spinal
muscular atrophy," for consideration in the Journal of Neurology.

This study addresses selective motor-neuron vulnerability in spinal muscular
atrophy by integrating several human transcriptomic systems with distinct and
predefined roles. Human post-mortem oculomotor and spinal motor neurons provide
a natural-resistance anchor; spinal-cord organoids provide disease and
SMN-restoration contrasts; an independent iPSC-derived motor-neuron cohort
provides external disease direction; and single-cell or single-nucleus datasets
provide motor-neuron-resolved context. Analyses are performed within datasets
before effect-direction integration, and inference is based on donors or donor
lines rather than technical libraries.

The study emphasizes reproducibility and calibrated negative evidence. No
individual oculomotor-positive gene survives donor-aware false-discovery
correction. Of 11,326 ranked genes, 37 meet the frozen cross-model direction
pattern, but this count is not enriched under gene-identity permutation
(`p = 0.836`). Complete raw-read re-quantification of nine pre-specified
GSE290979 libraries retains the full direction pattern for 21 candidates; nine
are raw biological-unit robust, and four genes (`LY6H`, `HS3ST5`, `ZNF853`, and
`IL17D`) are robust in both raw and processed analyses. Their overlap exceeds
the robustness-label permutation expectation (`p = 0.044`), and their
biological-unit bootstrap selection frequencies are 89.8-95.9%.

The manuscript also reports a stringent negative splice-junction result. STAR
and rMATS recover only six of 83 frozen structures in both contrasts, and none
meets the complete direction, reversal, and two-line correction rule. We
therefore label the processed splicing panel hypothesis-generating and avoid
therapeutic or neuroprotective claims. The accompanying repository contains
the frozen manifests, donor-aware code, raw-data provenance, regression tests,
negative controls, machine-readable tables, and figures needed to reproduce
the reported computational results.

We believe the work is relevant to readers interested in SMA biology,
selective neuronal vulnerability, human disease models, and reproducible
transcriptomic inference.

The corresponding author must add the author-approved signature, institutional
address, contact information, competing-interest and funding confirmations,
and confirmation that the manuscript is not under consideration elsewhere
before this letter is submitted.
