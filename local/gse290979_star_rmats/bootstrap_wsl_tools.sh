#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_PATH="${1:-${SCRIPT_DIR}/config.env}"
[[ -f "${CONFIG_PATH}" ]] || {
  printf 'Missing config: %s\n' "${CONFIG_PATH}" >&2
  exit 1
}
# shellcheck source=/dev/null
source "${CONFIG_PATH}"

MICROMAMBA="${TOOLS_ROOT}/bin/micromamba"
ENV_ROOT="${TOOLS_ROOT}/env"
PROVENANCE_ROOT="${WORK_ROOT}/provenance"
mkdir -p "${TOOLS_ROOT}/bin" "${PROVENANCE_ROOT}"

if [[ ! -x "${MICROMAMBA}" ]]; then
  archive="${TOOLS_ROOT}/micromamba-linux-64.tar.bz2"
  curl -fL --retry 5 --retry-delay 3 \
    https://micro.mamba.pm/api/micromamba/linux-64/latest \
    -o "${archive}"
  tar -xjf "${archive}" -C "${TOOLS_ROOT}" bin/micromamba
fi

packages=(
  "python=3.10"
  "star=${EXPECTED_STAR_VERSION}"
  "rmats=${EXPECTED_RMATS_VERSION}"
  "rmats2sashimiplot=${EXPECTED_RMATS2SASHIMIPLOT_VERSION}"
  "fastqc=${EXPECTED_FASTQC_VERSION}"
  "samtools=${EXPECTED_SAMTOOLS_VERSION}"
  "multiqc=${EXPECTED_MULTIQC_VERSION}"
)
if [[ ! -d "${ENV_ROOT}/conda-meta" ]]; then
  "${MICROMAMBA}" create -y -p "${ENV_ROOT}" \
    -c conda-forge -c bioconda \
    "${packages[@]}"
else
  "${MICROMAMBA}" install -y -p "${ENV_ROOT}" \
    -c conda-forge -c bioconda \
    "${packages[@]}"
fi

export PATH="${ENV_ROOT}/bin:${PATH}"
[[ "$(STAR --version)" = "${EXPECTED_STAR_VERSION}" ]]
python "${ENV_ROOT}/bin/rmats.py" --version 2>&1 |
  grep -F "${EXPECTED_RMATS_VERSION}" >/dev/null
fastqc --version 2>&1 | grep -F "${EXPECTED_FASTQC_VERSION}" >/dev/null
samtools --version | head -n 1 |
  grep -F "${EXPECTED_SAMTOOLS_VERSION}" >/dev/null
multiqc --version 2>&1 | grep -F "${EXPECTED_MULTIQC_VERSION}" >/dev/null
python -c \
  'import importlib.metadata; print(importlib.metadata.version("rmats2sashimiplot"))' |
  grep -Fx "${EXPECTED_RMATS2SASHIMIPLOT_VERSION}" >/dev/null

{
  printf 'tool\tversion\n'
  printf 'STAR\t%s\n' "$(STAR --version)"
  printf 'rMATS\t%s\n' \
    "$(python "${ENV_ROOT}/bin/rmats.py" --version 2>&1 | head -n 1)"
  printf 'FastQC\t%s\n' "$(fastqc --version 2>&1 | head -n 1)"
  printf 'samtools\t%s\n' "$(samtools --version | head -n 1)"
  printf 'MultiQC\t%s\n' "$(multiqc --version 2>&1 | head -n 1)"
  printf 'rmats2sashimiplot\t%s\n' \
    "$(python -c 'import importlib.metadata; print(importlib.metadata.version("rmats2sashimiplot"))')"
  printf 'micromamba\t%s\n' "$("${MICROMAMBA}" --version)"
} > "${PROVENANCE_ROOT}/tool_versions.tsv"

"${MICROMAMBA}" list -p "${ENV_ROOT}" --explicit \
  > "${PROVENANCE_ROOT}/conda_environment_explicit.txt"
printf 'Tool environment ready: %s\n' "${ENV_ROOT}"
