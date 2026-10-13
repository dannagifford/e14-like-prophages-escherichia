#!/usr/bin/env bash
set -euo pipefail

# Wrapper for the icd-aware boundary finder.
#
# Usage:
#   ./02_find_boundaries.sh
#
# Optional thresholds can be passed directly to the Python script, e.g.
#   ./02_find_boundaries.sh --min-icd-identity 80

MANIFEST="${MANIFEST:-manifest.tsv}"
TOOLS_ENV="prophage_compare"

if ! conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "${TOOLS_ENV}"; then
    echo "ERROR: conda environment '${TOOLS_ENV}' does not exist." >&2
    echo "Run 00_setup.sh first." >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2
    exit 1
fi

if [[ ! -f "02_find_boundaries.py" ]]; then
    echo "ERROR: 02_find_boundaries.py not found in $(pwd)" >&2
    exit 1
fi

conda run --no-capture-output -n "${TOOLS_ENV}" \
    python 02_find_boundaries.py \
    --manifest "${MANIFEST}" \
    "$@"

echo "Next step:"
echo "  ./03_annotate.sh"
