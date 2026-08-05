# Raw splice confirmation implementation amendment

Date: 2026-08-05

## Reason

The preregistered fixed-event rMATS 4.3.0 runs emitted zero MATS rows,
including after GENCODE v47 GeneID reconciliation and `--novelSS`.
This was not an absence of raw splice evidence: STAR's `SJ.out.tab`
contained a frozen MRFAP1 junction with uniquely mapped reads, and a
target-locus de novo rMATS diagnostic recovered supported events from the
same BAMs.

## Confirmatory implementation

The primary raw inference therefore uses rMATS 4.3.0 with `--novelSS`
on BAMs restricted to the 83 frozen event loci plus 1 kb flanks. Disease
and treatment contrasts are run independently. The resulting raw events
are matched back to the frozen 83-event manifest by event class, gene,
chromosome, strand, and every event-defining coordinate.

Only exact matches to the frozen 83 events are confirmatory. Other
events discovered inside the target loci are not interpreted as new
findings and are outside the confirmatory result set. The frozen genes,
event panel, effect directions, thresholds, and success criteria remain
unchanged.

## Provenance

This amendment was made before inspecting matched disease-treatment
concordance for the frozen panel. Failed fixed-event outputs are retained
under `failed_runs`; FASTQ and BAM files are preserved.
