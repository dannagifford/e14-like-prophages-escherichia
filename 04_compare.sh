#!/usr/bin/env bash
set -euo pipefail

# 04_compare.sh
#
# Final comparison stage:
#
#   1. all-vs-all nucleotide BLAST summary of final extracted sequences
#   2. multi-way Clinker using consistently annotated Phold GenBanks
#   3. both auto-ordered and manifest-ordered Clinker HTML plots
#   4. MG1655-curated vs PHANOTATE e14 annotation comparison
#   5. manifest-order Clinker with the original MG1655 e14 annotation
#
# Inputs:
#   manifest.tsv
#   prophages/<label>/<label>.fna
#   annotation/final_gbk/<label>.gbk
#
# Outputs:
#   comparisons/blast/
#   comparisons/clinker/
#
# Usage:
#   ./04_compare.sh
#
# Optional:
#   THREADS=16 ./04_compare.sh
#   CLINKER_IDENTITY=0.5 ./04_compare.sh
#   FORCE=1 ./04_compare.sh

THREADS="${THREADS:-8}"
CLINKER_IDENTITY="${CLINKER_IDENTITY:-0.3}"
FORCE="${FORCE:-0}"

ROOT="$(pwd)"
MANIFEST="${MANIFEST:-${ROOT}/manifest.tsv}"

COMPARE_ENV="prophage_compare"
CLINKER_ENV="clinker_py"

FINAL_GBK_ROOT="${ROOT}/annotation/final_gbk"
BLAST_OUT="${ROOT}/comparisons/blast"
CLINKER_OUT="${ROOT}/comparisons/clinker"
ANNOTATION_OUT="${ROOT}/comparisons/annotation_check"

E14_LABEL="${E14_LABEL:-Ecoli_K12_e14}"
E14_ACCESSION="${E14_ACCESSION:-NC_000913.3}"
E14_START="${E14_START:-1196209}"
E14_END="${E14_END:-1211412}"
REFERENCE_GBK="${REFERENCE_GBK:-${ROOT}/genomes/${E14_LABEL}/${E14_LABEL}.gbk}"

mkdir -p "${BLAST_OUT}" "${CLINKER_OUT}" "${ANNOTATION_OUT}" "${ROOT}/logs"

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2
    exit 1
fi

if [[ ! -f "${ROOT}/04_pairwise_similarity.py" ]]; then
    echo "ERROR: 04_pairwise_similarity.py not found in ${ROOT}" >&2
    exit 1
fi

if [[ ! -f "${ROOT}/compare_e14_annotations.py" ]]; then
    echo "ERROR: compare_e14_annotations.py not found in ${ROOT}" >&2
    exit 1
fi

env_exists () {
    conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$1"
}

for ENV in "${COMPARE_ENV}" "${CLINKER_ENV}"; do
    if ! env_exists "${ENV}"; then
        echo "ERROR: conda environment '${ENV}' does not exist." >&2
        echo "Run 00_setup.sh first." >&2
        exit 1
    fi
done

# Build ordered arrays from manifest.
LABELS=()
GBKS=()
E14_FOUND=0

while IFS=$'\t' read -r label accession type; do
    [[ -n "${label}" ]] || continue

    GBK="${FINAL_GBK_ROOT}/${label}.gbk"
    FASTA="${ROOT}/prophages/${label}/${label}.fna"

    if [[ ! -s "${GBK}" ]]; then
        echo "ERROR: missing final annotated GenBank: ${GBK}" >&2
        echo "Run 03_annotate.sh first." >&2
        exit 1
    fi

    if [[ ! -s "${FASTA}" ]]; then
        echo "ERROR: missing final FASTA: ${FASTA}" >&2
        exit 1
    fi

    LABELS+=("${label}")
    GBKS+=("${GBK}")
    if [[ "${label}" == "${E14_LABEL}" ]]; then
        E14_FOUND=1
    fi
done < <(tail -n +2 "${MANIFEST}" | cut -f1,3,4)

echo "Sequences in comparison:"
printf '  %s\n' "${LABELS[@]}"
echo

if [[ "${E14_FOUND}" != "1" ]]; then
    echo "ERROR: e14 label '${E14_LABEL}' was not found in manifest.tsv." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 0. K-12 e14 annotation-control comparison
# ---------------------------------------------------------------------------
#
# The curated MG1655 track is the SAME e14 nucleotide sequence with the
# original NC_000913.3 CDS models. It is deliberately excluded from
# manifest.tsv, pairwise nucleotide comparisons and the standard Clinker
# session consumed by stage 05.
# ---------------------------------------------------------------------------

PIPELINE_E14_GBK="${FINAL_GBK_ROOT}/${E14_LABEL}.gbk"
CURATED_E14_GBK="${ANNOTATION_OUT}/MG1655_curated_e14.gbk"
ANNOTATION_REPORT="${ANNOTATION_OUT}/report.txt"

if [[ "${FORCE}" == "1" || ! -s "${CURATED_E14_GBK}" || ! -s "${ANNOTATION_REPORT}" ]]; then
    echo "=== e14 annotation comparison: MG1655 curated vs PHANOTATE ==="

    ANNOTATION_ARGS=(
        python "${ROOT}/compare_e14_annotations.py"
        --pharokka-gbk "${PIPELINE_E14_GBK}"
        --accession "${E14_ACCESSION}"
        --region-start "${E14_START}"
        --region-end "${E14_END}"
        --outdir "${ANNOTATION_OUT}"
    )

    if [[ -s "${REFERENCE_GBK}" ]]; then
        ANNOTATION_ARGS+=(--reference-gbk "${REFERENCE_GBK}")
    fi

    conda run --no-capture-output -n "${COMPARE_ENV}" \
        "${ANNOTATION_ARGS[@]}" \
        2>&1 | tee "${ROOT}/logs/04_e14_annotation_comparison.log"
else
    echo "[reuse] e14 annotation comparison: ${ANNOTATION_REPORT}"
fi

if [[ ! -s "${CURATED_E14_GBK}" ]]; then
    echo "ERROR: annotation comparison did not create ${CURATED_E14_GBK}" >&2
    exit 1
fi

GBKS_WITH_CURATED=()
for i in "${!LABELS[@]}"; do
    GBKS_WITH_CURATED+=("${GBKS[$i]}")
    if [[ "${LABELS[$i]}" == "${E14_LABEL}" ]]; then
        GBKS_WITH_CURATED+=("${CURATED_E14_GBK}")
    fi
done

echo

# ---------------------------------------------------------------------------
# 1. Pairwise nucleotide similarity
# ---------------------------------------------------------------------------

PAIRWISE_TSV="${BLAST_OUT}/pairwise_similarity.tsv"

if [[ "${FORCE}" == "1" || ! -s "${PAIRWISE_TSV}" ]]; then
    echo "=== All-vs-all nucleotide BLAST ==="

    conda run --no-capture-output -n "${COMPARE_ENV}" \
        python "${ROOT}/04_pairwise_similarity.py" \
        --manifest "${MANIFEST}" \
        --root "${ROOT}" \
        --outdir "comparisons/blast" \
        2>&1 | tee "${ROOT}/logs/04_pairwise_similarity.log"
else
    echo "[reuse] pairwise similarity: ${PAIRWISE_TSV}"
fi

echo

# ---------------------------------------------------------------------------
# 2. Clinker: automatic ordering
# ---------------------------------------------------------------------------

AUTO_HTML="${CLINKER_OUT}/e14_relatives_auto_order.html"
AUTO_ALIGN="${CLINKER_OUT}/e14_relatives_auto_order_alignments.csv"
AUTO_MATRIX="${CLINKER_OUT}/e14_relatives_auto_order_matrix.csv"
AUTO_SESSION="${CLINKER_OUT}/e14_relatives_auto_order_session.json"

if [[ "${FORCE}" == "1" || ! -s "${AUTO_HTML}" ]]; then
    echo "=== Clinker: automatic similarity-based order ==="

    # Always start a fresh Clinker session.
    rm -f "${AUTO_SESSION}"

    echo "=== Clinker: automatic similarity-based order ==="

    conda run --no-capture-output -n "${CLINKER_ENV}" \
        clinker \
        "${GBKS[@]}" \
        --identity "${CLINKER_IDENTITY}" \
        --jobs "${THREADS}" \
        --session "${AUTO_SESSION}" \
        --output "${AUTO_ALIGN}" \
        --delimiter "," \
        --decimals 4 \
        --matrix_out "${AUTO_MATRIX}" \
        --plot "${AUTO_HTML}" \
        --force \
        2>&1 | tee "${ROOT}/logs/04_clinker_auto.log"
else
    echo "[reuse] auto-order Clinker: ${AUTO_HTML}"
fi

echo

# ---------------------------------------------------------------------------
# 3. Clinker: manifest order
# ---------------------------------------------------------------------------

MANIFEST_HTML="${CLINKER_OUT}/e14_relatives_manifest_order.html"
MANIFEST_ALIGN="${CLINKER_OUT}/e14_relatives_manifest_order_alignments.csv"
MANIFEST_MATRIX="${CLINKER_OUT}/e14_relatives_manifest_order_matrix.csv"
MANIFEST_SESSION="${CLINKER_OUT}/e14_relatives_manifest_order_session.json"

if [[ "${FORCE}" == "1" || ! -s "${MANIFEST_HTML}" ]]; then
    echo "=== Clinker: manifest order ==="

    # Always start a fresh Clinker session.
    rm -f "${MANIFEST_SESSION}"

    conda run --no-capture-output -n "${CLINKER_ENV}" \
        clinker \
        "${GBKS[@]}" \
        --identity "${CLINKER_IDENTITY}" \
        --jobs "${THREADS}" \
        --session "${MANIFEST_SESSION}" \
        --output "${MANIFEST_ALIGN}" \
        --delimiter "," \
        --decimals 4 \
        --matrix_out "${MANIFEST_MATRIX}" \
        --use_file_order \
        --plot "${MANIFEST_HTML}" \
        --force \
        2>&1 | tee "${ROOT}/logs/04_clinker_manifest.log"
else
    echo "[reuse] manifest-order Clinker: ${MANIFEST_HTML}"
fi

echo

# ---------------------------------------------------------------------------
# 4. Clinker: manifest order plus original MG1655 e14 annotation
# ---------------------------------------------------------------------------
#
# This uses a separate session because stage 05 requires the standard session
# to contain manifest taxa only.
# ---------------------------------------------------------------------------

CURATED_HTML="${CLINKER_OUT}/e14_relatives_manifest_order_with_MG1655.html"
CURATED_ALIGN="${CLINKER_OUT}/e14_relatives_manifest_order_alignments_with_MG1655.csv"
CURATED_MATRIX="${CLINKER_OUT}/e14_relatives_manifest_order_matrix_with_MG1655.csv"
CURATED_SESSION="${CLINKER_OUT}/e14_relatives_manifest_order_with_MG1655_session.json"

if [[ "${FORCE}" == "1" || ! -s "${CURATED_HTML}" ]]; then
    echo "=== Clinker: manifest order + original MG1655 e14 annotation ==="

    rm -f "${CURATED_SESSION}"

    conda run --no-capture-output -n "${CLINKER_ENV}" \
        clinker \
        "${GBKS_WITH_CURATED[@]}" \
        --identity "${CLINKER_IDENTITY}" \
        --jobs "${THREADS}" \
        --session "${CURATED_SESSION}" \
        --output "${CURATED_ALIGN}" \
        --delimiter "," \
        --decimals 4 \
        --matrix_out "${CURATED_MATRIX}" \
        --use_file_order \
        --plot "${CURATED_HTML}" \
        --force \
        2>&1 | tee "${ROOT}/logs/04_clinker_manifest_with_MG1655.log"
else
    echo "[reuse] manifest-order Clinker: ${CURATED_HTML}"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

echo
echo "======================================================================"
echo "Stage 04 complete."
echo
echo "Pairwise nucleotide comparison:"
echo "  ${BLAST_OUT}/pairwise_similarity.tsv"
echo "  ${BLAST_OUT}/mean_reciprocal_coverage_matrix.tsv"
echo "  ${BLAST_OUT}/mean_weighted_identity_matrix.tsv"
echo
echo "Clinker:"
echo "  ${AUTO_HTML}"
echo "  ${MANIFEST_HTML}"
echo "  ${CURATED_HTML}  (includes original MG1655 annotation control)"
echo
echo "e14 annotation comparison:"
echo "  ${ANNOTATION_OUT}/mg1655_vs_phanotate_summary.tsv"
echo "  ${ANNOTATION_OUT}/overlapping_cds_pairs.tsv"
echo "  ${ANNOTATION_OUT}/annotation_tracks.svg"
echo "  ${ANNOTATION_OUT}/report.txt"
echo "  ${CURATED_E14_GBK}"
echo
echo "Clinker identity threshold: ${CLINKER_IDENTITY}"
echo
echo "I recommend inspecting the auto-order plot first, then using the"
echo "manifest-order version when you want a fixed biological/taxonomic layout."
echo "Next step:"
echo "  ./05_finalize.sh"
