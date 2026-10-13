#!/usr/bin/env bash
set -euo pipefail

# 03_annotate.sh
#
# Standardised annotation:
#   final FASTA -> sanitised FASTA -> Pharokka -> Phold
#
# Important:
# Phold rejects ':' in GenBank DEFINITION lines. Pharokka can propagate the
# FASTA description into DEFINITION, so this script first writes a clean FASTA
# whose header is simply:
#
#   >LABEL
#
# This prevents genomic coordinate strings such as
# NC_000913.3:1196209-1211412 from reaching the Pharokka GenBank DEFINITION.
#
# Existing valid outputs are reused unless FORCE=1. Existing Pharokka GenBanks
# containing ':' in DEFINITION are treated as invalid and regenerated.
#
# Usage:
#   ./03_annotate.sh
#
# Optional:
#   THREADS=16 ./03_annotate.sh
#   FORCE=1 ./03_annotate.sh
#   PHOLD_GPU=1 ./03_annotate.sh

THREADS="${THREADS:-8}"
FORCE="${FORCE:-0}"
PHOLD_GPU="${PHOLD_GPU:-0}"

ROOT="$(pwd)"
MANIFEST="${MANIFEST:-${ROOT}/manifest.tsv}"

PHAROKKA_ENV="pharokka"
PHOLD_ENV="phold"

PHAROKKA_DB="${PHAROKKA_DB:-${ROOT}/databases/pharokka}"
PHOLD_DB="${PHOLD_DB:-${ROOT}/databases/phold}"

PROPHAGE_ROOT="${ROOT}/prophages"
CLEAN_INPUT_ROOT="${ROOT}/annotation/clean_input"
PHAROKKA_ROOT="${ROOT}/annotation/pharokka"
PHOLD_ROOT="${ROOT}/annotation/phold"
FINAL_GBK_ROOT="${ROOT}/annotation/final_gbk"
MASH_ROOT="${ROOT}/annotation/mash"
SUMMARY="${ROOT}/annotation/annotation_summary.tsv"

mkdir -p \
    "${CLEAN_INPUT_ROOT}" \
    "${PHAROKKA_ROOT}" \
    "${PHOLD_ROOT}" \
    "${FINAL_GBK_ROOT}" \
    "${MASH_ROOT}" \
    "${ROOT}/logs"

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2
    exit 1
fi

env_exists () {
    conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$1"
}

for ENV in "${PHAROKKA_ENV}" "${PHOLD_ENV}"; do
    if ! env_exists "${ENV}"; then
        echo "ERROR: conda environment '${ENV}' does not exist." >&2
        echo "Run 00_setup.sh first." >&2
        exit 1
    fi
done

if [[ ! -d "${PHAROKKA_DB}" ]]; then
    echo "ERROR: Pharokka database not found: ${PHAROKKA_DB}" >&2
    exit 1
fi

if [[ ! -d "${PHOLD_DB}" ]]; then
    echo "ERROR: Phold database not found: ${PHOLD_DB}" >&2
    exit 1
fi

valid_fasta () {
    local file="$1"
    [[ -s "${file}" ]] || return 1
    head -n 1 "${file}" | grep -q '^>'
}

valid_genbank () {
    local file="$1"
    [[ -s "${file}" ]] || return 1
    head -n 1 "${file}" | grep -q '^LOCUS' || return 1
    tail -n 20 "${file}" | grep -q '^//$'
}

# Phold specifically rejects ':' in DEFINITION.
definition_has_colon () {
    local gbk="$1"
    grep '^DEFINITION' "${gbk}" | grep -q ':'
}

find_single_gbk () {
    local dir="$1"

    mapfile -t gbks < <(
        find "${dir}" -maxdepth 1 -type f \
            \( -name '*.gbk' -o -name '*.gb' -o -name '*.gbff' \) \
            -print | sort
    )

    if (( ${#gbks[@]} == 0 )); then
        return 1
    fi

    if (( ${#gbks[@]} > 1 )); then
        echo "ERROR: expected one top-level GenBank file in ${dir}, found:" >&2
        printf '  %s\n' "${gbks[@]}" >&2
        return 2
    fi

    printf '%s\n' "${gbks[0]}"
}

count_cds () {
    local gbk="$1"
    grep -c '^     CDS' "${gbk}" || true
}

fasta_length () {
    local fasta="$1"
    awk '
        /^>/ {next}
        {
            gsub(/[[:space:]]/, "", $0)
            n += length($0)
        }
        END {print n+0}
    ' "${fasta}"
}

make_clean_fasta () {
    local input="$1"
    local output="$2"
    local label="$3"

    # These inputs should all be single-contig sequences.
    local nheaders
    nheaders="$(grep -c '^>' "${input}" || true)"

    if [[ "${nheaders}" -ne 1 ]]; then
        echo "ERROR: ${input} contains ${nheaders} FASTA records; expected 1." >&2
        return 1
    fi

    {
        printf ">%s\n" "${label}"
        grep -v '^>' "${input}" | tr -d '\r'
    } > "${output}"

    valid_fasta "${output}"
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

printf "label\taccession\ttype\tlength_bp\tclean_fasta\tpharokka_gbk\tphold_gbk\tfinal_gbk\tcds_count\tmash_top_hits\tstatus\n" \
    > "${SUMMARY}"

# ---------------------------------------------------------------------------
# Annotate
# ---------------------------------------------------------------------------

tail -n +2 "${MANIFEST}" | cut -f1,3,4 | \
while IFS=$'\t' read -r label accession type; do

    [[ -n "${label}" ]] || continue

    SOURCE_FASTA="${PROPHAGE_ROOT}/${label}/${label}.fna"
    CLEAN_FASTA="${CLEAN_INPUT_ROOT}/${label}.fna"

    PHAROKKA_OUT="${PHAROKKA_ROOT}/${label}"
    PHOLD_OUT="${PHOLD_ROOT}/${label}"
    FINAL_GBK="${FINAL_GBK_ROOT}/${label}.gbk"
    LOG="${ROOT}/logs/03_annotate_${label}.log"

    echo "======================================================================"
    echo "${label}"
    echo "  accession: ${accession}"
    echo "  type:      ${type}"

    if ! valid_fasta "${SOURCE_FASTA}"; then
        echo "ERROR: missing/invalid input FASTA: ${SOURCE_FASTA}" >&2
        echo "Run and review stage 02 first." >&2
        exit 1
    fi

    LENGTH="$(fasta_length "${SOURCE_FASTA}")"

    # Always recreate the clean header copy: tiny, deterministic, and avoids
    # stale descriptions from an earlier pipeline version.
    make_clean_fasta "${SOURCE_FASTA}" "${CLEAN_FASTA}" "${label}"

    echo "  source:      ${SOURCE_FASTA}"
    echo "  clean input: ${CLEAN_FASTA}"
    echo "  header:      $(head -n 1 "${CLEAN_FASTA}")"
    echo "  length:      ${LENGTH} bp"

    # -----------------------------------------------------------------------
    # Reuse final output if it is valid and Phold-compatible.
    # -----------------------------------------------------------------------

    if [[ "${FORCE}" != "1" ]] && valid_genbank "${FINAL_GBK}" \
        && ! definition_has_colon "${FINAL_GBK}"; then

        echo "  [reuse] final Phold GenBank: ${FINAL_GBK}"

        PHAROKKA_GBK=""
        if [[ -d "${PHAROKKA_OUT}" ]]; then
            PHAROKKA_GBK="$(find_single_gbk "${PHAROKKA_OUT}" 2>/dev/null || true)"
        fi

        PHOLD_GBK=""
        if [[ -d "${PHOLD_OUT}" ]]; then
            PHOLD_GBK="$(find_single_gbk "${PHOLD_OUT}" 2>/dev/null || true)"
        fi

        CDS_COUNT="$(count_cds "${FINAL_GBK}")"
        MASH_FILE="$(find "${PHAROKKA_OUT}" -maxdepth 1 -type f \
            -name '*top_hits_mash_inphared.tsv' -print -quit 2>/dev/null || true)"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${label}" \
            "${accession}" \
            "${type}" \
            "${LENGTH}" \
            "${CLEAN_FASTA#${ROOT}/}" \
            "${PHAROKKA_GBK#${ROOT}/}" \
            "${PHOLD_GBK#${ROOT}/}" \
            "${FINAL_GBK#${ROOT}/}" \
            "${CDS_COUNT}" \
            "${MASH_FILE#${ROOT}/}" \
            "reused" \
            >> "${SUMMARY}"

        continue
    fi

    # -----------------------------------------------------------------------
    # Pharokka
    # -----------------------------------------------------------------------

    PHAROKKA_GBK=""
    NEED_PHAROKKA=1

    if [[ "${FORCE}" != "1" ]] && [[ -d "${PHAROKKA_OUT}" ]]; then
        PHAROKKA_GBK="$(find_single_gbk "${PHAROKKA_OUT}" 2>/dev/null || true)"

        if [[ -n "${PHAROKKA_GBK}" ]] \
            && valid_genbank "${PHAROKKA_GBK}" \
            && ! definition_has_colon "${PHAROKKA_GBK}"; then
            NEED_PHAROKKA=0
            echo "  [reuse] Pharokka: ${PHAROKKA_GBK}"
        elif [[ -n "${PHAROKKA_GBK}" ]] && definition_has_colon "${PHAROKKA_GBK}"; then
            echo "  [rerun] existing Pharokka GenBank has ':' in DEFINITION"
        fi
    fi

    if [[ "${NEED_PHAROKKA}" == "1" ]]; then
        echo "  [run] Pharokka from sanitised FASTA"

        # Remove the old directory entirely so no stale GenBank is accidentally
        # picked up after the rerun.
        rm -rf "${PHAROKKA_OUT}"

        conda run --no-capture-output -n "${PHAROKKA_ENV}" \
            pharokka run \
            -i "${CLEAN_FASTA}" \
            -o "${PHAROKKA_OUT}" \
            -d "${PHAROKKA_DB}" \
            -t "${THREADS}" \
            -g phanotate \
            -p "${label}" \
            -l "${label}" \
            -f \
            2>&1 | tee "${LOG}"

        PHAROKKA_GBK="$(find_single_gbk "${PHAROKKA_OUT}")"

        if ! valid_genbank "${PHAROKKA_GBK}"; then
            echo "ERROR: Pharokka GenBank failed validation: ${PHAROKKA_GBK}" >&2
            exit 1
        fi

        if definition_has_colon "${PHAROKKA_GBK}"; then
            echo "ERROR: Pharokka GenBank still contains ':' in DEFINITION:" >&2
            grep '^DEFINITION' "${PHAROKKA_GBK}" >&2
            exit 1
        fi
    fi

    echo "  Pharokka DEFINITION:"
    grep '^DEFINITION' "${PHAROKKA_GBK}" || true

    # -----------------------------------------------------------------------
    # Mash / INPHARED
    # -----------------------------------------------------------------------

    MASH_FILE="$(find "${PHAROKKA_OUT}" -maxdepth 1 -type f \
        -name '*top_hits_mash_inphared.tsv' -print -quit || true)"

    if [[ -n "${MASH_FILE}" ]] && [[ -s "${MASH_FILE}" ]]; then
        MASH_COPY="${MASH_ROOT}/${label}_top_hits_mash_inphared.tsv"
        cp "${MASH_FILE}" "${MASH_COPY}"
        echo "  Mash hits: ${MASH_COPY}"
    else
        MASH_COPY=""
        echo "WARNING: no Pharokka Mash/INPHARED table found for ${label}" >&2
    fi

    # -----------------------------------------------------------------------
    # Phold
    # -----------------------------------------------------------------------

    echo "  [run] Phold"

    # Remove any partial output left by a previous failed Phold attempt.
    rm -rf "${PHOLD_OUT}"

    PHOLD_ARGS=(
        phold run
        -i "${PHAROKKA_GBK}"
        -o "${PHOLD_OUT}"
        -d "${PHOLD_DB}"
        -t "${THREADS}"
        -p "${label}"
        -f
    )

    # Do not use --hyps: annotate all CDS consistently.
    if [[ "${PHOLD_GPU}" == "1" ]]; then
        PHOLD_ARGS+=(--foldseek_gpu)
        echo "  Phold mode: GPU + Foldseek-GPU"
    else
        PHOLD_ARGS+=(--cpu)
        echo "  Phold mode: CPU"
    fi

    conda run --no-capture-output -n "${PHOLD_ENV}" \
        "${PHOLD_ARGS[@]}" \
        2>&1 | tee -a "${LOG}"

    PHOLD_GBK="$(find_single_gbk "${PHOLD_OUT}")"

    if ! valid_genbank "${PHOLD_GBK}"; then
        echo "ERROR: Phold GenBank failed validation: ${PHOLD_GBK}" >&2
        exit 1
    fi

    # -----------------------------------------------------------------------
    # Stable downstream copy
    # -----------------------------------------------------------------------

    cp "${PHOLD_GBK}" "${FINAL_GBK}"

    CDS_COUNT="$(count_cds "${FINAL_GBK}")"

    echo "  final GBK: ${FINAL_GBK}"
    echo "  CDS:       ${CDS_COUNT}"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${label}" \
        "${accession}" \
        "${type}" \
        "${LENGTH}" \
        "${CLEAN_FASTA#${ROOT}/}" \
        "${PHAROKKA_GBK#${ROOT}/}" \
        "${PHOLD_GBK#${ROOT}/}" \
        "${FINAL_GBK#${ROOT}/}" \
        "${CDS_COUNT}" \
        "${MASH_COPY#${ROOT}/}" \
        "completed" \
        >> "${SUMMARY}"

done

echo
echo "======================================================================"
echo "Annotation stage complete."
echo
echo "Summary:"
if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "${SUMMARY}"
else
    cat "${SUMMARY}"
fi

echo
echo "Final GenBanks:"
find "${FINAL_GBK_ROOT}" -maxdepth 1 -type f -name '*.gbk' -print | sort
echo "Next step:"
echo "  ./04_compare.sh"
