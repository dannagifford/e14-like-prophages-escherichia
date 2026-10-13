#!/usr/bin/env bash
set -euo pipefail

# 00_setup.sh
#
# Idempotent setup for the e14-related prophage comparison pipeline.
#
# Behaviour:
#   * reuses existing conda environments
#   * ensures required packages are installed in existing environments
#   * creates missing environments
#   * reuses existing Pharokka/Phold databases by default
#   * optionally refreshes databases with REFRESH_DBS=1
#   * creates the project directory structure
#   * verifies required executables/imports
#   * records versions for reproducibility
#
# Usage:
#   ./00_setup.sh
#
# Optional:
#   THREADS=16 ./00_setup.sh
#   REFRESH_DBS=1 ./00_setup.sh
#
# Note:
#   Updating Pharokka can occasionally require a matching newer database.
#   If Pharokka reports a database-version incompatibility later, rerun:
#
#       REFRESH_DBS=1 ./00_setup.sh

THREADS="${THREADS:-8}"
REFRESH_DBS="${REFRESH_DBS:-0}"

ROOT="$(pwd)"
DB_ROOT="${ROOT}/databases"
PHAROKKA_DB="${PHAROKKA_DB:-${DB_ROOT}/pharokka}"
PHOLD_DB="${PHOLD_DB:-${DB_ROOT}/phold}"

mkdir -p \
    "${DB_ROOT}" \
    "${ROOT}/genomes" \
    "${ROOT}/context" \
    "${ROOT}/boundaries" \
    "${ROOT}/prophages" \
    "${ROOT}/annotation/pharokka" \
    "${ROOT}/annotation/phold" \
    "${ROOT}/annotation/final_gbk" \
    "${ROOT}/annotation/mash" \
    "${ROOT}/comparisons/blast" \
    "${ROOT}/comparisons/clinker" \
    "${ROOT}/logs" \
    "${ROOT}/versions"

echo "Project root: ${ROOT}"
echo "Threads:      ${THREADS}"
echo "Refresh DBs:  ${REFRESH_DBS}"
echo

# ---------------------------------------------------------------------------
# Basic checks
# ---------------------------------------------------------------------------

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available in PATH." >&2
    exit 1
fi

if [[ ! -f "${ROOT}/manifest.tsv" ]]; then
    echo "WARNING: manifest.tsv is not present in ${ROOT}."
    echo "         Setup can continue, but later stages require it."
    echo
fi

env_exists () {
    conda env list \
        | awk 'NF && $1 !~ /^#/ {print $1}' \
        | grep -Fxq "$1"
}

# Ensure packages are present whether or not the environment already exists.
ensure_env () {
    local env_name="$1"
    shift
    local packages=("$@")

    if env_exists "${env_name}"; then
        echo "[update] ensuring packages are present in: ${env_name}"
        conda install -y \
            -n "${env_name}" \
            -c conda-forge \
            -c bioconda \
            "${packages[@]}"
    else
        echo "[create] conda environment: ${env_name}"
        conda create -y \
            -n "${env_name}" \
            -c conda-forge \
            -c bioconda \
            "${packages[@]}"
    fi
}

# ---------------------------------------------------------------------------
# Conda environments
# ---------------------------------------------------------------------------

echo "=== Conda environments ==="

# General sequence handling, NCBI retrieval, boundary inference and tabular work.
# pandas and Biopython are explicitly required by 02_find_boundaries.py.
ensure_env prophage_compare \
    blast \
    samtools \
    seqkit \
    entrez-direct \
    biopython \
    pandas \
    mafft

# Standardised phage annotation.
ensure_env pharokka \
    pharokka

# Structure-informed annotation.
# Existing GPU-enabled environments are preserved as environments; this command
# ensures phold itself is installed but does not intentionally replace CUDA
# configuration.
ensure_env phold \
    phold

# Comparative gene-cluster visualisation.
ensure_env clinker_py \
    clinker-py

echo

# ---------------------------------------------------------------------------
# Verify executables and Python imports
# ---------------------------------------------------------------------------

echo "=== Verifying tools ==="

conda run -n prophage_compare python - <<'PY'
import pandas
import Bio
print("pandas", pandas.__version__)
print("Biopython", Bio.__version__)
PY

conda run -n prophage_compare blastn -version | head -n 1
conda run -n prophage_compare samtools --version | head -n 1
conda run -n prophage_compare seqkit version

# edirect installs efetch as a shell/perl executable.
conda run -n prophage_compare bash -c 'command -v efetch >/dev/null && echo "efetch OK"'

conda run -n pharokka pharokka --version
conda run -n phold phold --version
conda run -n clinker_py clinker --version

echo

# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

dir_has_content () {
    local dir="$1"
    [[ -d "${dir}" ]] &&
    [[ -n "$(find "${dir}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]
}

# ---------------------------------------------------------------------------
# Pharokka database
# ---------------------------------------------------------------------------

echo "=== Pharokka database ==="

if [[ "${REFRESH_DBS}" == "1" ]]; then
    echo "[refresh] removing existing Pharokka DB: ${PHAROKKA_DB}"
    rm -rf "${PHAROKKA_DB}"
fi

if dir_has_content "${PHAROKKA_DB}"; then
    echo "[reuse] Pharokka database: ${PHAROKKA_DB}"
else
    echo "[install] Pharokka database: ${PHAROKKA_DB}"
    mkdir -p "$(dirname "${PHAROKKA_DB}")"

    conda run --no-capture-output -n pharokka \
        pharokka install \
        -o "${PHAROKKA_DB}"
fi

echo

# ---------------------------------------------------------------------------
# Phold database
# ---------------------------------------------------------------------------

echo "=== Phold database ==="

if [[ "${REFRESH_DBS}" == "1" ]]; then
    echo "[refresh] removing existing Phold DB: ${PHOLD_DB}"
    rm -rf "${PHOLD_DB}"
fi

if dir_has_content "${PHOLD_DB}"; then
    echo "[reuse] Phold database: ${PHOLD_DB}"
else
    echo "[install] Phold database: ${PHOLD_DB}"
    mkdir -p "$(dirname "${PHOLD_DB}")"

    conda run --no-capture-output -n phold \
        phold install \
        -d "${PHOLD_DB}" \
        -t "${THREADS}"
fi

echo

# ---------------------------------------------------------------------------
# Reproducibility record
# ---------------------------------------------------------------------------

VERSIONS="${ROOT}/versions/setup_versions.txt"

{
    echo "Generated: $(date --iso-8601=seconds)"
    echo "Project root: ${ROOT}"
    echo

    echo "== Conda =="
    conda --version
    echo

    echo "== prophage_compare =="
    conda run -n prophage_compare python - <<'PY'
import sys
import pandas
import Bio
print("Python", sys.version.replace("\n", " "))
print("pandas", pandas.__version__)
print("Biopython", Bio.__version__)
PY
    conda run -n prophage_compare blastn -version | head -n 1
    conda run -n prophage_compare samtools --version | head -n 1
    conda run -n prophage_compare seqkit version
    echo

    echo "== Pharokka =="
    conda run -n pharokka pharokka --version
    echo "database=${PHAROKKA_DB}"
    echo

    echo "== Phold =="
    conda run -n phold phold --version
    echo "database=${PHOLD_DB}"
    echo

    echo "== Clinker =="
    conda run -n clinker_py clinker --version
} > "${VERSIONS}"

echo "=== Setup complete ==="
echo
echo "Pharokka DB: ${PHAROKKA_DB}"
echo "Phold DB:    ${PHOLD_DB}"
echo "Versions:    ${VERSIONS}"
echo
echo "If a later Pharokka/Phold run reports a database-version mismatch:"
echo "  REFRESH_DBS=1 ./00_setup.sh"
echo
echo "Next step:"
echo "  ./01_download.sh"
