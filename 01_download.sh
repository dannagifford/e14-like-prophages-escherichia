#!/usr/bin/env bash
set -euo pipefail

# 01_download.sh
#
# Download the complete NCBI nucleotide record for every accession in
# manifest.tsv. For each entry, save:
#
#   genomes/<label>/<label>.fna
#   genomes/<label>/<label>.gbk
#
# The script is idempotent:
#   - valid existing files are reused
#   - missing or invalid files are re-downloaded
#
# It also writes:
#   genomes/downloaded_genomes.tsv
#   genomes/checksums.sha256
#
# Usage:
#   bash 01_download.sh
#
# Optional:
#   MANIFEST=path/to/manifest.tsv bash 01_download.sh

MANIFEST="${MANIFEST:-manifest.tsv}"
TOOLS_ENV="prophage_compare"

ROOT="$(pwd)"
GENOME_ROOT="${ROOT}/genomes"
SUMMARY="${GENOME_ROOT}/downloaded_genomes.tsv"
CHECKSUMS="${GENOME_ROOT}/checksums.sha256"

mkdir -p "${GENOME_ROOT}" "${ROOT}/logs"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2
    exit 1
fi

if ! conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "${TOOLS_ENV}"; then
    echo "ERROR: conda environment '${TOOLS_ENV}' does not exist." >&2
    echo "Run 00_setup.sh first." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

valid_fasta () {
    local file="$1"

    [[ -s "${file}" ]] || return 1
    head -n 1 "${file}" | grep -q '^>' || return 1

    # Require at least one nucleotide after removing FASTA headers/whitespace.
    awk '
        /^>/ {next}
        {
            gsub(/[[:space:]]/, "", $0)
            n += length($0)
        }
        END {exit !(n > 0)}
    ' "${file}"
}

valid_genbank () {
    local file="$1"

    [[ -s "${file}" ]] || return 1
    head -n 1 "${file}" | grep -q '^LOCUS' || return 1

    # A complete GenBank flatfile record should terminate with //.
    tail -n 20 "${file}" | grep -q '^//$'
}

retry () {
    local attempts="$1"
    shift

    local n=1
    until "$@"; do
        if (( n >= attempts )); then
            echo "ERROR: command failed after ${attempts} attempts:" >&2
            printf '  %q' "$@" >&2
            echo >&2
            return 1
        fi

        echo "Download failed (attempt ${n}/${attempts}); retrying in 5 s..." >&2
        sleep 5
        ((n++))
    done
}

fetch_and_validate () {
    local accession="$1"
    local format="$2"
    local validator="$3"
    local output="$4"

    conda run --no-capture-output -n "${TOOLS_ENV}" \
        efetch \
        -db nuccore \
        -id "${accession}" \
        -format "${format}" \
        > "${output}" \
        && "${validator}" "${output}"
}

download_fasta () {
    local accession="$1"
    local output="$2"
    local tmp="${output}.tmp"

    rm -f "${tmp}"

    retry 3 fetch_and_validate \
        "${accession}" fasta valid_fasta "${tmp}"

    if ! valid_fasta "${tmp}"; then
        echo "ERROR: downloaded FASTA failed validation for ${accession}" >&2
        rm -f "${tmp}"
        return 1
    fi

    mv "${tmp}" "${output}"
}

download_genbank () {
    local accession="$1"
    local output="$2"
    local tmp="${output}.tmp"

    rm -f "${tmp}"

    retry 3 fetch_and_validate \
        "${accession}" gbwithparts valid_genbank "${tmp}"

    if ! valid_genbank "${tmp}"; then
        echo "ERROR: downloaded GenBank failed validation for ${accession}" >&2
        rm -f "${tmp}"
        return 1
    fi

    mv "${tmp}" "${output}"
}

fasta_length () {
    local file="$1"

    awk '
        /^>/ {next}
        {
            gsub(/[[:space:]]/, "", $0)
            n += length($0)
        }
        END {print n}
    ' "${file}"
}

fasta_header () {
    local file="$1"
    head -n 1 "${file}" | sed 's/^>//'
}

# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

printf "label\torganism\taccession\ttype\tfasta\tgenbank\tlength_bp\tfasta_header\n" \
    > "${SUMMARY}"

# We only need columns 1-4 here, so using cut avoids problems with the empty
# seed_start/seed_end fields in the free-phage row.
tail -n +2 "${MANIFEST}" | cut -f1-4 | \
while IFS=$'\t' read -r label organism accession type; do

    [[ -n "${label}" ]] || continue
    [[ -n "${accession}" ]] || {
        echo "ERROR: missing accession for ${label}" >&2
        exit 1
    }

    OUTDIR="${GENOME_ROOT}/${label}"
    FASTA="${OUTDIR}/${label}.fna"
    GBK="${OUTDIR}/${label}.gbk"

    mkdir -p "${OUTDIR}"

    echo "============================================================"
    echo "${label}"
    echo "  organism:  ${organism}"
    echo "  accession: ${accession}"
    echo "  type:      ${type}"

    if valid_fasta "${FASTA}"; then
        echo "  [reuse] FASTA   ${FASTA}"
    else
        echo "  [fetch] FASTA   ${accession}"
        rm -f "${FASTA}"
        download_fasta "${accession}" "${FASTA}"
    fi

    if valid_genbank "${GBK}"; then
        echo "  [reuse] GenBank ${GBK}"
    else
        echo "  [fetch] GenBank ${accession}"
        rm -f "${GBK}"
        download_genbank "${accession}" "${GBK}"
    fi

    LENGTH="$(fasta_length "${FASTA}")"
    HEADER="$(fasta_header "${FASTA}" | tr '\t' ' ')"

    echo "  length:    ${LENGTH} bp"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${label}" \
        "${organism}" \
        "${accession}" \
        "${type}" \
        "${FASTA#${ROOT}/}" \
        "${GBK#${ROOT}/}" \
        "${LENGTH}" \
        "${HEADER}" \
        >> "${SUMMARY}"

done

# ---------------------------------------------------------------------------
# Checksums
# ---------------------------------------------------------------------------

find "${GENOME_ROOT}" -type f \( -name '*.fna' -o -name '*.gbk' \) \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKSUMS}"

echo
echo "============================================================"
echo "Download stage complete."
echo
echo "Summary:"
column -t -s $'\t' "${SUMMARY}" || cat "${SUMMARY}"
echo
echo "Checksums:"
echo "  ${CHECKSUMS}"
echo
echo "Next step:"
echo "  02_find_boundaries / attL-attR inference"
