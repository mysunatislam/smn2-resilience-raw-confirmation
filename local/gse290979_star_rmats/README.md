# GSE290979 Local9 STAR-rMATS Confirmation

This workflow performs prospective raw splice-junction confirmation for the
frozen 83-event set and 12-event primary panel. It uses the nine
checksum-verified local9 FASTQ pairs.

## Scope

- GENCODE v47 GRCh38 primary assembly.
- STAR 2.7.10a with a sparse local index.
- Full-depth, annotation-guided one-pass STAR (`twopassMode=None`,
  `readMapNumber=-1`).
- rMATS turbo 4.3.0 with the frozen fixed-event set.
- rmats2sashimiplot 4.0.0 for the frozen 12-event panel.
- Disease: untreated SMA S2/S3 versus control C1/C2/C3.
- Treatment: R6-MO versus scramble in S2/S3.
- STAR output is streamed directly into coordinate-sorted BAMs restricted to
  the 83 event loci with 1 kb padding.
- STAR FIFO and sorting scratch files use WSL's Linux `/tmp`; NTFS does not
  support the FIFO files required by STAR.
- A full BAM is never created. FASTQs and retained BAMs are never deleted.

The sparse STAR index is a local feasibility adaptation. It does not alter the
reference sequence, annotation, event coordinates, or read depth, but it is
slower than a standard human STAR index. One-pass alignment is locked because
the confirmatory set contains annotated events already embedded in the STAR
splice-junction database.

## Memory Gate

The host has 11.9 GiB RAM. Standard human STAR requires at least 16 GiB and is
not locally feasible. The locked benchmark uses:

```text
genomeSAsparseD=8
genomeSAindexNbases=12
limitGenomeGenerateRAM=9000000000
```

Do not start the index phase unless WSL reports at least 8 GiB RAM. Benchmark
one library before starting the remaining eight.

## Prepare

From the repository root in Windows:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' r\29_prepare_local9_star_rmats_profile.R
Copy-Item local\gse290979_star_rmats\config.env.example local\gse290979_star_rmats\config.env
```

In WSL:

```bash
cd /mnt/c/Users/Kotha/OneDrive/Documents/smn2/publication_repo
bash local/gse290979_star_rmats/bootstrap_wsl_tools.sh
```

The bootstrap creates a version-pinned micromamba environment under the
E-drive work root. It does not modify the system Python.

## Run Order

```bash
bash local/gse290979_star_rmats/run_local9_star_rmats.sh reference
bash local/gse290979_star_rmats/run_local9_star_rmats.sh preflight
bash local/gse290979_star_rmats/run_local9_star_rmats.sh index
bash local/gse290979_star_rmats/run_local9_star_rmats.sh fastqc
bash local/gse290979_star_rmats/run_local9_star_rmats.sh align BULK-SAM-148
```

Review the sparse-index and one-library benchmark before running:

```bash
bash local/gse290979_star_rmats/run_local9_star_rmats.sh align
bash local/gse290979_star_rmats/run_local9_star_rmats.sh rmats
bash local/gse290979_star_rmats/run_local9_star_rmats.sh sashimi
```

Then run the frozen R analysis:

```powershell
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' `
  r\30_analyze_local9_star_rmats.R `
  --rmats-root=E:\smn2_gse290979_local9\star_rmats\rmats `
  --output-root=E:\smn2_gse290979_local9\star_rmats\analysis
```

The analysis must use the dated freeze in
`docs/validation_freeze_2026-08-01.md`.

The fixed-set FDR scope and targeted BAM-retention deviation are documented
prospectively in
`docs/raw_splice_confirmation_execution_deviation_2026-08-01.md`.
