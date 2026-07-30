# GSE290979 local9 raw-read workflow

Windows/WSL workflow used for the completed nine-library raw-read sensitivity
analysis.

## Scope

The locked local9 profile retains:

- all five untreated donor lines for SMA versus control;
- both SMA donor lines for R6-MO versus scramble;
- one complete technical preparation per line-condition; and
- 82.193 GiB of paired compressed FASTQ input.

The selected libraries and deposited ENA checksums are in
`config/GSE290979_local9_sample_sheet.tsv`.

## Completed route

1. Download both complete ENA FASTQs for each selected library.
2. Require both deposited MD5 values.
3. Download and verify the GENCODE v47 transcriptome.
4. Build a Salmon 1.10.2 transcriptome index with `--gencode`.
5. Quantify paired reads with automatic library detection,
   `--validateMappings`, `--seqBias`, and `--gcBias`.
6. Import transcripts with tximport and
   `countsFromAbundance="lengthScaledTPM"`.
7. Fit the five-line disease and paired two-line treatment models.
8. Run five disease LOLO folds and separate S2/S3 treatment checks.
9. Compare raw-derived effects with the deposited processed matrix.
10. Run the 13-check integrity gate.

All nine quantifications passed. The completion marker is:

```text
results/r/raw_confirmation/local9_salmon/
LOCAL9_SALMON_RAW_REQUANT_COMPLETE.tsv
```

It records:

```text
quantification_method=salmon_transcriptome_quasimapping
whole_genome_alignment=FALSE
random_split_used=FALSE
```

## Execute

The worker uses `run_local9.ps1` for resumable checksum-verified ENA transfer
and WSL for Salmon:

```powershell
& .\local\gse290979_raw\run_local9_salmon.ps1 -Phase all
```

Run it in a hidden background process:

```powershell
& .\local\gse290979_raw\start_local9_salmon_background.ps1
```

Inspect status without modifying files:

```powershell
& .\local\gse290979_raw\get_local9_status.ps1
```

The analysis phase calls:

```text
r/27_import_gse290979_local9_salmon.R
r/24_analyze_gse290979_local9_counts.R
r/25_compare_gse290979_local9_to_processed.R
r/26_validate_gse290979_local9_outputs.R
```

Publication integration is performed separately by:

```text
r/28_integrate_local9_salmon_publication.R
```

## Retention and claim boundary

Verified FASTQs and existing BAMs were preserved. No FASTQ or BAM deletion is
part of this workflow.

The completed analysis supports all-nine raw-read transcriptome
re-quantification. It does not support whole-genome alignment, raw
splice-junction discovery, independent cohort validation, or definitive
`SMN1`/`SMN2` allele-specific quantification.
