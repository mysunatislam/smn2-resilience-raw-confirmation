#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_PATH="${CONFIG_PATH:-${SCRIPT_DIR}/config.env}"
PHASE="${1:-}"
SAMPLE_ID="${2:-}"

[[ -n "${PHASE}" ]] || {
  printf 'Usage: %s {preflight|reference|index|fastqc|align|rmats|sashimi} [sample_id]\n' \
    "$0" >&2
  exit 1
}
[[ -f "${CONFIG_PATH}" ]] || {
  printf 'Missing config: %s\n' "${CONFIG_PATH}" >&2
  exit 1
}
# shellcheck source=/dev/null
source "${CONFIG_PATH}"

if [[ "${DELETE_FASTQ:-0}" != "0" || "${DELETE_BAM:-0}" != "0" ]]; then
  printf 'FASTQ and BAM deletion must remain disabled\n' >&2
  exit 1
fi

ENV_ROOT="${TOOLS_ROOT}/env"
export PATH="${ENV_ROOT}/bin:${PATH}"
PROFILE_DIR="${REPO_ROOT}/config/rmats/GSE290979/local9"
FIXED_EVENT_DIR="${REPO_ROOT}/config/rmats/GSE290979/fixed_events"
SAMPLE_SHEET="${REPO_ROOT}/config/GSE290979_local9_sample_sheet.tsv"
TARGET_BED="${PROFILE_DIR}/fixed_events_padded_1kb_merged.bed"
REFERENCE_FASTA_GZ="${REFERENCE_SOURCE_ROOT}/GRCh38.primary_assembly.genome.fa.gz"
REFERENCE_GTF_GZ="${REFERENCE_SOURCE_ROOT}/gencode.v47.primary_assembly.annotation.gtf.gz"
REFERENCE_FASTA="${STAR_REFERENCE_ROOT}/GRCh38.primary_assembly.genome.fa"
REFERENCE_GTF="${STAR_REFERENCE_ROOT}/gencode.v47.primary_assembly.annotation.gtf"
STAR_INDEX="${STAR_REFERENCE_ROOT}/index_sparseD${STAR_GENOME_SA_SPARSE_D}_SA${STAR_GENOME_SA_INDEX_NBASES}"
BAM_ROOT="${WORK_ROOT}/bam"
QC_ROOT="${WORK_ROOT}/qc"
RMATS_ROOT="${WORK_ROOT}/rmats"
SASHIMI_ROOT="${WORK_ROOT}/sashimi"
PROVENANCE_ROOT="${WORK_ROOT}/provenance"
TMP_ROOT="${WORK_ROOT}/tmp"

mkdir -p \
  "${BAM_ROOT}" "${QC_ROOT}" "${RMATS_ROOT}" "${SASHIMI_ROOT}" \
  "${PROVENANCE_ROOT}" "${TMP_ROOT}"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

require_file() {
  [[ -s "$1" ]] || die "Missing or empty file: $1"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "Missing command: $1"
}

check_tools() {
  require_command STAR
  require_command samtools
  require_command fastqc
  require_command multiqc
  require_command rmats2sashimiplot
  require_command python
  require_file "${ENV_ROOT}/bin/rmats.py"
  [[ "$(STAR --version)" = "${EXPECTED_STAR_VERSION}" ]] ||
    die "Unexpected STAR version"
  python "${ENV_ROOT}/bin/rmats.py" --version 2>&1 |
    grep -F "${EXPECTED_RMATS_VERSION}" >/dev/null ||
    die "Unexpected rMATS version"
  python -c \
    'import importlib.metadata; print(importlib.metadata.version("rmats2sashimiplot"))' |
    grep -Fx "${EXPECTED_RMATS2SASHIMIPLOT_VERSION}" >/dev/null ||
    die "Unexpected rmats2sashimiplot version"
}

sample_row() {
  local sample_id="$1"
  awk -F '\t' -v sample="${sample_id}" '
    NR == 1 { next }
    $2 == sample {
      gsub(/\r$/, "", $0)
      print
      found++
    }
    END {
      if (found != 1) {
        exit 2
      }
    }
  ' "${SAMPLE_SHEET}"
}

fastq_paths() {
  local sample_id="$1"
  local row run_accession
  row="$(sample_row "${sample_id}")"
  run_accession="$(printf '%s\n' "${row}" | cut -f3)"
  printf '%s/%s/%s_1.fastq.gz\t%s/%s/%s_2.fastq.gz\n' \
    "${FASTQ_ROOT}" "${sample_id}" "${run_accession}" \
    "${FASTQ_ROOT}" "${sample_id}" "${run_accession}"
}

verify_fastq_pair() {
  local sample_id="$1"
  local row r1 r2 r1_md5 r2_md5
  row="$(sample_row "${sample_id}")"
  r1_md5="$(printf '%s\n' "${row}" | cut -f16)"
  r2_md5="$(printf '%s\n' "${row}" | cut -f17)"
  IFS=$'\t' read -r r1 r2 <<< "$(fastq_paths "${sample_id}")"
  require_file "${r1}"
  require_file "${r2}"
  printf '%s  %s\n' "${r1_md5}" "${r1}" | md5sum -c -
  printf '%s  %s\n' "${r2_md5}" "${r2}" | md5sum -c -
}

all_sample_ids() {
  awk -F '\t' 'NR > 1 { gsub(/\r$/, "", $2); print $2 }' "${SAMPLE_SHEET}"
}

phase_preflight() {
  require_file "${SAMPLE_SHEET}"
  require_file "${TARGET_BED}"
  require_file "${REFERENCE_FASTA_GZ}"
  require_file "${REFERENCE_GTF_GZ}"
  require_file "${REFERENCE_MD5_FILE}"
  check_tools

  local memory_kib
  memory_kib="$(awk '/^MemTotal:/ { print $2 }' /proc/meminfo)"
  {
    printf 'metric\tvalue\n'
    printf 'completed_utc\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    printf 'wsl_memory_gib\t%.3f\n' \
      "$(awk -v value="${memory_kib}" 'BEGIN { print value / 1024 / 1024 }')"
    printf 'cpu_threads\t%s\n' "$(nproc)"
    printf 'target_bed_sha256\t%s\n' \
      "$(sha256sum "${TARGET_BED}" | awk '{ print $1 }')"
    printf 'frozen_events_sha256\t%s\n' \
      "$(sha256sum "${REPO_ROOT}/config/frozen_83_splice_events_2026-08-01.tsv" | awk '{ print $1 }')"
    printf 'delete_fastq\tFALSE\n'
    printf 'delete_bam\tFALSE\n'
    df -BG "${WORK_ROOT}" |
      awk 'NR == 2 { printf "work_root_free_gib\t%s\n", $4 }'
  } > "${PROVENANCE_ROOT}/preflight.tsv"

  while IFS= read -r sample_id; do
    verify_fastq_pair "${sample_id}"
  done < <(all_sample_ids)
  printf 'Preflight passed for all nine checksum-verified libraries\n'
}

expected_md5() {
  local filename="$1"
  awk -v filename="${filename}" '$2 == filename { print $1 }' \
    "${REFERENCE_MD5_FILE}"
}

phase_reference() {
  require_file "${REFERENCE_FASTA_GZ}"
  require_file "${REFERENCE_GTF_GZ}"
  require_file "${REFERENCE_MD5_FILE}"
  mkdir -p "${STAR_REFERENCE_ROOT}"

  local fasta_name gtf_name fasta_md5 gtf_md5
  fasta_name="$(basename "${REFERENCE_FASTA_GZ}")"
  gtf_name="$(basename "${REFERENCE_GTF_GZ}")"
  fasta_md5="$(expected_md5 "${fasta_name}")"
  gtf_md5="$(expected_md5 "${gtf_name}")"
  [[ -n "${fasta_md5}" && -n "${gtf_md5}" ]] ||
    die "Official GENCODE MD5 entries are missing"
  printf '%s  %s\n' "${fasta_md5}" "${REFERENCE_FASTA_GZ}" | md5sum -c -
  printf '%s  %s\n' "${gtf_md5}" "${REFERENCE_GTF_GZ}" | md5sum -c -

  if [[ ! -s "${REFERENCE_FASTA}" ]]; then
    gzip -dc "${REFERENCE_FASTA_GZ}" > "${REFERENCE_FASTA}.partial"
    mv "${REFERENCE_FASTA}.partial" "${REFERENCE_FASTA}"
  fi
  if [[ ! -s "${REFERENCE_GTF}" ]]; then
    gzip -dc "${REFERENCE_GTF_GZ}" > "${REFERENCE_GTF}.partial"
    mv "${REFERENCE_GTF}.partial" "${REFERENCE_GTF}"
  fi

  {
    printf 'file\tsha256\n'
    printf '%s\t%s\n' "$(basename "${REFERENCE_FASTA}")" \
      "$(sha256sum "${REFERENCE_FASTA}" | awk '{ print $1 }')"
    printf '%s\t%s\n' "$(basename "${REFERENCE_GTF}")" \
      "$(sha256sum "${REFERENCE_GTF}" | awk '{ print $1 }')"
  } > "${PROVENANCE_ROOT}/reference_sha256.tsv"
  printf 'GENCODE v47 reference ready: %s\n' "${STAR_REFERENCE_ROOT}"
}

phase_index() {
  check_tools
  require_file "${REFERENCE_FASTA}"
  require_file "${REFERENCE_GTF}"
  local memory_kib minimum_kib
  memory_kib="$(awk '/^MemTotal:/ { print $2 }' /proc/meminfo)"
  minimum_kib=$((8 * 1024 * 1024))
  (( memory_kib >= minimum_kib )) ||
    die "Sparse STAR index benchmark requires at least 8 GiB WSL memory"

  if [[ -s "${STAR_INDEX}/Genome" && -s "${STAR_INDEX}/SA" ]]; then
    printf 'Sparse STAR index already complete: %s\n' "${STAR_INDEX}"
    return
  fi
  [[ ! -e "${STAR_INDEX}" ]] || die "Partial STAR index exists: ${STAR_INDEX}"
  local temporary_index
  temporary_index="${STAR_INDEX}.partial"
  mkdir -p "${temporary_index}"
  STAR \
    --runThreadN "${THREADS}" \
    --runMode genomeGenerate \
    --genomeDir "${temporary_index}" \
    --genomeFastaFiles "${REFERENCE_FASTA}" \
    --sjdbGTFfile "${REFERENCE_GTF}" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    --genomeSAsparseD "${STAR_GENOME_SA_SPARSE_D}" \
    --genomeSAindexNbases "${STAR_GENOME_SA_INDEX_NBASES}" \
    --limitGenomeGenerateRAM "${STAR_LIMIT_GENOME_GENERATE_RAM}"
  require_file "${temporary_index}/Genome"
  require_file "${temporary_index}/SA"
  mv "${temporary_index}" "${STAR_INDEX}"
  {
    printf 'metric\tvalue\n'
    printf 'STAR_version\t%s\n' "$(STAR --version)"
    printf 'genomeSAsparseD\t%s\n' "${STAR_GENOME_SA_SPARSE_D}"
    printf 'genomeSAindexNbases\t%s\n' "${STAR_GENOME_SA_INDEX_NBASES}"
    printf 'sjdbOverhang\t%s\n' "${SJDB_OVERHANG}"
    printf 'index_bytes\t%s\n' "$(du -sb "${STAR_INDEX}" | awk '{ print $1 }')"
  } > "${PROVENANCE_ROOT}/star_index.tsv"
  printf 'Sparse STAR index complete: %s\n' "${STAR_INDEX}"
}

phase_fastqc() {
  check_tools
  local fastqc_root
  fastqc_root="${QC_ROOT}/fastqc"
  mkdir -p "${fastqc_root}"
  while IFS= read -r sample_id; do
    local r1 r2 r1_name r2_name r1_output r2_output
    IFS=$'\t' read -r r1 r2 <<< "$(fastq_paths "${sample_id}")"
    require_file "${r1}"
    require_file "${r2}"
    r1_name="$(basename "${r1}")"
    r2_name="$(basename "${r2}")"
    r1_output="${fastqc_root}/${r1_name%.fastq.gz}_fastqc.zip"
    r2_output="${fastqc_root}/${r2_name%.fastq.gz}_fastqc.zip"
    if [[ ! -s "${r1_output}" || ! -s "${r2_output}" ]]; then
      fastqc --threads "${FASTQC_THREADS}" --outdir "${fastqc_root}" \
        "${r1}" "${r2}"
    fi
    require_file "${r1_output}"
    require_file "${r2_output}"
  done < <(all_sample_ids)
  multiqc --force --outdir "${QC_ROOT}/multiqc" "${fastqc_root}"
  require_file "${QC_ROOT}/multiqc/multiqc_report.html"
  printf 'FastQC and MultiQC complete\n'
}

align_one() {
  local sample_id="$1"
  check_tools
  require_file "${STAR_INDEX}/Genome"
  require_file "${STAR_INDEX}/SA"
  require_file "${TARGET_BED}"
  verify_fastq_pair "${sample_id}"

  local r1 r2 output_directory final_bam temporary_directory partial_bam
  local star_linux_tmp
  IFS=$'\t' read -r r1 r2 <<< "$(fastq_paths "${sample_id}")"
  output_directory="${BAM_ROOT}/${sample_id}"
  final_bam="${output_directory}/${sample_id}.star.event_loci.bam"
  if [[ -s "${final_bam}" && -s "${final_bam}.bai" ]]; then
    samtools quickcheck -v "${final_bam}"
    printf 'Alignment already complete: %s\n' "${sample_id}"
    return
  fi
  [[ ! -e "${output_directory}" ]] ||
    die "Partial alignment output exists: ${output_directory}"

  temporary_directory="${TMP_ROOT}/align_${sample_id}"
  [[ ! -e "${temporary_directory}" ]] ||
    die "Partial alignment temporary directory exists: ${temporary_directory}"
  mkdir -p "${temporary_directory}"
  partial_bam="${temporary_directory}/${sample_id}.star.event_loci.partial.bam"
  star_linux_tmp="${STAR_LINUX_TMP_ROOT}/${sample_id}"
  mkdir -p "${STAR_LINUX_TMP_ROOT}"
  [[ ! -e "${star_linux_tmp}" ]] ||
    die "Partial Linux STAR temporary directory exists: ${star_linux_tmp}"

  STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${r1}" "${r2}" \
    --readFilesCommand zcat \
    --outTmpDir "${star_linux_tmp}" \
    --twopassMode Basic \
    --outFileNamePrefix "${temporary_directory}/star." \
    --outStd SAM \
    --outSAMtype SAM \
    --outSAMattributes NH HI AS nM NM MD \
    --outSAMstrandField intronMotif \
    --outSAMunmapped None \
    --quantMode GeneCounts |
    samtools view -@ 1 -u -L "${TARGET_BED}" - |
    samtools sort \
      -@ "${SORT_THREADS}" \
      -m 768M \
      -T "${temporary_directory}/sort" \
      -o "${partial_bam}" -

  samtools quickcheck -v "${partial_bam}"
  samtools index -@ "${SORT_THREADS}" "${partial_bam}"
  samtools flagstat -@ "${SORT_THREADS}" "${partial_bam}" \
    > "${temporary_directory}/samtools.flagstat.txt"
  samtools stats -@ "${SORT_THREADS}" "${partial_bam}" \
    > "${temporary_directory}/samtools.stats.txt"
  local retained_alignments
  retained_alignments="$(samtools view -c "${partial_bam}")"
  (( retained_alignments > 0 )) ||
    die "No target-locus alignments retained for ${sample_id}"

  mkdir -p "${output_directory}"
  mv "${partial_bam}" "${final_bam}"
  mv "${partial_bam}.bai" "${final_bam}.bai"
  find "${temporary_directory}" -maxdepth 1 -type f -exec \
    mv -t "${output_directory}" -- {} +
  rmdir "${temporary_directory}"
  {
    printf 'metric\tvalue\n'
    printf 'status\tCOMPLETE\n'
    printf 'sample_id\t%s\n' "${sample_id}"
    printf 'completed_utc\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    printf 'STAR_version\t%s\n' "$(STAR --version)"
    printf 'retention\tevent_loci_only_streamed_from_STAR_SAM\n'
    printf 'full_BAM_created\tFALSE\n'
    printf 'FASTQ_deleted\tFALSE\n'
    printf 'BAM_deleted\tFALSE\n'
    printf 'retained_alignments\t%s\n' "${retained_alignments}"
    printf 'BAM_bytes\t%s\n' "$(stat -c '%s' "${final_bam}")"
    printf 'BAM_sha256\t%s\n' \
      "$(sha256sum "${final_bam}" | awk '{ print $1 }')"
    printf 'target_BED_sha256\t%s\n' \
      "$(sha256sum "${TARGET_BED}" | awk '{ print $1 }')"
  } > "${output_directory}/ALIGNMENT_COMPLETE.tsv"
  printf 'Target-locus STAR alignment complete: %s\n' "${sample_id}"
}

sample_list_to_bams() {
  local list_path="$1"
  local bam_paths=()
  while IFS= read -r sample_id; do
    local bam
    sample_id="${sample_id%$'\r'}"
    bam="${BAM_ROOT}/${sample_id}/${sample_id}.star.event_loci.bam"
    require_file "${bam}"
    require_file "${bam}.bai"
    samtools quickcheck -v "${bam}"
    bam_paths+=("${bam}")
  done < "${list_path}"
  local joined
  joined="$(IFS=,; printf '%s' "${bam_paths[*]}")"
  printf '%s\n' "${joined}"
}

run_rmats_contrast() {
  local contrast="$1"
  local group1="$2"
  local group2="$3"
  local output_directory tmp_directory
  output_directory="${RMATS_ROOT}/${contrast}"
  tmp_directory="${RMATS_ROOT}/tmp_${contrast}"
  if [[ -s "${output_directory}/SUCCESS.tsv" ]]; then
    printf 'rMATS already complete: %s\n' "${contrast}"
    return
  fi
  [[ ! -e "${output_directory}" ]] ||
    die "Partial rMATS output exists: ${output_directory}"
  [[ ! -e "${tmp_directory}" ]] ||
    die "Partial rMATS temporary directory exists: ${tmp_directory}"
  mkdir -p "${output_directory}" "${tmp_directory}"

  local b1_file b2_file
  b1_file="${output_directory}/b1.txt"
  b2_file="${output_directory}/b2.txt"
  sample_list_to_bams "${PROFILE_DIR}/${group1}.txt" > "${b1_file}"
  sample_list_to_bams "${PROFILE_DIR}/${group2}.txt" > "${b2_file}"

  python "${ENV_ROOT}/bin/rmats.py" \
    --b1 "${b1_file}" \
    --b2 "${b2_file}" \
    --gtf "${REFERENCE_GTF}" \
    -t paired \
    --readLength "${READ_LENGTH}" \
    --libType "${LIB_TYPE}" \
    --nthread "${THREADS}" \
    --tstat "${THREADS}" \
    --od "${output_directory}" \
    --tmp "${tmp_directory}" \
    --task both \
    --individual-counts \
    --fixed-event-set "${FIXED_EVENT_DIR}"

  for event_type in SE A5SS A3SS MXE RI; do
    require_file "${output_directory}/${event_type}.MATS.JC.txt"
  done
  {
    printf 'metric\tvalue\n'
    printf 'status\tCOMPLETE\n'
    printf 'contrast\t%s\n' "${contrast}"
    printf 'completed_utc\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    printf 'rMATS_version\t%s\n' \
      "$(python "${ENV_ROOT}/bin/rmats.py" --version 2>&1 | head -n 1)"
    printf 'fixed_event_set\tTRUE\n'
    printf 'event_count\t83\n'
    printf 'group1\t%s\n' "${group1}"
    printf 'group2\t%s\n' "${group2}"
  } > "${output_directory}/SUCCESS.tsv"
  printf 'rMATS complete: %s\n' "${contrast}"
}

phase_rmats() {
  check_tools
  require_file "${REFERENCE_GTF}"
  run_rmats_contrast \
    disease disease_group1_sma disease_group2_control
  run_rmats_contrast \
    treatment treatment_group1_r6 treatment_group2_scramble
}

phase_sashimi() {
  check_tools
  require_file "${RMATS_ROOT}/disease/SUCCESS.tsv"
  require_file \
    "${REPO_ROOT}/config/frozen_12_primary_splice_events_2026-08-01.tsv"
  local complete_marker
  complete_marker="${SASHIMI_ROOT}/SASHIMI_12_COMPLETE.tsv"
  if [[ -s "${complete_marker}" ]]; then
    printf 'Twelve-event sashimi panel already complete\n'
    return
  fi

  local disease_bams treatment_bams group_file manifest
  disease_bams="$(
    sample_list_to_bams \
      "${PROFILE_DIR}/disease_group1_sma.txt"
  ),$(
    sample_list_to_bams \
      "${PROFILE_DIR}/disease_group2_control.txt"
  )"
  treatment_bams="$(
    sample_list_to_bams \
      "${PROFILE_DIR}/treatment_group1_r6.txt"
  ),$(
    sample_list_to_bams \
      "${PROFILE_DIR}/treatment_group2_scramble.txt"
  )"
  group_file="${SASHIMI_ROOT}/four_groups.gf"
  manifest="${SASHIMI_ROOT}/sashimi_manifest.tsv"
  printf '%s\n' "${disease_bams}" > "${SASHIMI_ROOT}/disease_bams.txt"
  printf '%s\n' "${treatment_bams}" > "${SASHIMI_ROOT}/treatment_bams.txt"
  {
    printf 'SMA untreated: 1-2\n'
    printf 'Control: 3-5\n'
    printf 'R6-MO: 6-7\n'
    printf 'Scramble: 8-9\n'
  } > "${group_file}"
  printf 'panel_order\tgene_symbol\tevent_type\tdisease_event_id\tpdf_path\tsha256\n' \
    > "${manifest}"

  while IFS=$'\t' read -r panel_order gene_symbol event_type event_id; do
    local rmats_type source_file event_file output_directory pdf_count pdf_path
    rmats_type="${event_type}"
    [[ "${event_type}" != "ES" ]] || rmats_type="SE"
    source_file="${RMATS_ROOT}/disease/${rmats_type}.MATS.JC.txt"
    require_file "${source_file}"
    output_directory="$(
      printf '%s/%02d_%s_%s\n' \
        "${SASHIMI_ROOT}" "${panel_order}" "${gene_symbol}" "${event_type}"
    )"
    event_file="${output_directory}/single_event.MATS.JC.txt"
    if [[ ! -e "${output_directory}" ]]; then
      mkdir -p "${output_directory}"
      awk -F '\t' -v id="${event_id}" \
        'NR == 1 || $1 == id { print }' \
        "${source_file}" > "${event_file}"
      [[ "$(wc -l < "${event_file}")" -eq 2 ]] ||
        die "Expected one raw event for ${gene_symbol} ${event_type}"
      rmats2sashimiplot \
        --b1 "${SASHIMI_ROOT}/disease_bams.txt" \
        --b2 "${SASHIMI_ROOT}/treatment_bams.txt" \
        --event-type "${rmats_type}" \
        -e "${event_file}" \
        --group-info "${group_file}" \
        --exon_s 1 \
        --intron_s 20 \
        --min-counts 1 \
        --font-size 8 \
        --fig-height 8 \
        --fig-width 10 \
        -o "${output_directory}/plot"
    fi
    pdf_count="$(
      find "${output_directory}" -type f -name '*.pdf' | wc -l
    )"
    [[ "${pdf_count}" -eq 1 ]] ||
      die "Expected one sashimi PDF for ${gene_symbol}, found ${pdf_count}"
    pdf_path="$(
      find "${output_directory}" -type f -name '*.pdf' -print -quit
    )"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${panel_order}" "${gene_symbol}" "${event_type}" "${event_id}" \
      "${pdf_path}" "$(sha256sum "${pdf_path}" | awk '{ print $1 }')" \
      >> "${manifest}"
  done < <(
    awk -F '\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          column[$i] = i
        }
        next
      }
      {
        print $column["panel_order"] "\t" \
          $column["source_gene_symbol"] "\t" \
          $column["event_type"] "\t" \
          $column["disease_event_id"]
      }
    ' "${REPO_ROOT}/config/frozen_12_primary_splice_events_2026-08-01.tsv"
  )

  [[ "$(($(wc -l < "${manifest}") - 1))" -eq 12 ]] ||
    die "Sashimi manifest does not contain 12 events"
  {
    printf 'metric\tvalue\n'
    printf 'status\tCOMPLETE\n'
    printf 'completed_utc\t%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    printf 'primary_events\t12\n'
    printf 'group_tracks\t4\n'
    printf 'rmats2sashimiplot_version\t%s\n' \
      "$(python -c 'import importlib.metadata; print(importlib.metadata.version("rmats2sashimiplot"))')"
    printf 'FASTQ_deleted\tFALSE\n'
    printf 'BAM_deleted\tFALSE\n'
  } > "${complete_marker}"
  printf 'Twelve-event sashimi panel complete: %s\n' "${manifest}"
}

case "${PHASE}" in
  preflight)
    phase_preflight
    ;;
  reference)
    phase_reference
    ;;
  index)
    phase_index
    ;;
  fastqc)
    phase_fastqc
    ;;
  align)
    if [[ -n "${SAMPLE_ID}" ]]; then
      align_one "${SAMPLE_ID}"
    else
      while IFS= read -r current_sample; do
        align_one "${current_sample}"
      done < <(all_sample_ids)
    fi
    ;;
  rmats)
    phase_rmats
    ;;
  sashimi)
    phase_sashimi
    ;;
  *)
    die "Unknown phase: ${PHASE}"
    ;;
esac
