#!/usr/bin/env bash
set -euo pipefail

# 05_finalize.sh
#
# This is the original Stage 05 logic, but run on the Clinker
# comparison that also contains:
#
#   MG1655_curated_e14
#
# i.e. the same K-12 e14 nucleotide sequence with the original NC_000913.3
# CDS annotation.
#
# It deliberately keeps the original group-relabelling method:
#   - Clinker defines the sequence-similarity groups
#   - 05_finalize.py reviews/harmonises their product annotations
#   - the saved Clinker session JSON is relabelled afterwards
#
# It DOES NOT use --gene_functions, so ungrouped/singleton genes are not
# artificially assigned colours.
#
# First run:
#   ./05_finalize.sh
#
# This creates:
#   results/final/group_label_review.tsv
#
# Inspect rows where review_required=1. For each, choose final_product and set:
#
#   approved=yes
#
# Then:
#   MODE=apply ./05_finalize.sh
#
# Apply mode:
#   1. writes harmonised GenBanks while preserving original products in /note
#   2. rebuilds Clinker from a clean session
#   3. writes manifest-order and auto-order Clinker plots
#   4. relabels ONLY existing Clinker homology groups using reviewed products
#   5. writes product-labelled plots
#
# Requires the Stage-04 outputs:
#   comparisons/annotation_check/MG1655_curated_e14.gbk
#   comparisons/clinker/
#     e14_relatives_manifest_order_session.json

MODE="${MODE:-review}"
THREADS="${THREADS:-8}"
CLINKER_IDENTITY="${CLINKER_IDENTITY:-0.3}"

ROOT="$(pwd)"
ENV="clinker_py"

BASE_MANIFEST="${MANIFEST:-${ROOT}/manifest.tsv}"
BASE_GBK_DIR="${ROOT}/annotation/final_gbk"

EXTRA_LABEL="MG1655_curated_e14"
EXTRA_GBK="${ROOT}/comparisons/annotation_check/MG1655_curated_e14.gbk"
E14_LABEL="${E14_LABEL:-Ecoli_K12_e14}"

SESSION="${ROOT}/comparisons/clinker/e14_relatives_manifest_order_with_MG1655_session.json"

OUTDIR="${ROOT}/results/final"
AUG_MANIFEST="${OUTDIR}/manifest_with_MG1655_annotation.tsv"
INPUT_GBK_DIR="${OUTDIR}/input_gbk"

REVIEW="${OUTDIR}/group_label_review.tsv"
GROUPS_TSV="${OUTDIR}/homologous_gene_groups.tsv"
HARMONISED="${OUTDIR}/harmonised_gbk"
CLINKER_OUT="${OUTDIR}/clinker"

MANIFEST_SESSION="${CLINKER_OUT}/final_manifest_order_session.json"
AUTO_SESSION="${CLINKER_OUT}/final_auto_order_session.json"

MANIFEST_PRODUCT_SESSION="${CLINKER_OUT}/final_manifest_order_product_labels.json"
AUTO_PRODUCT_SESSION="${CLINKER_OUT}/final_auto_order_product_labels.json"

mkdir -p \
    "${OUTDIR}" \
    "${INPUT_GBK_DIR}" \
    "${CLINKER_OUT}" \
    "${ROOT}/logs"

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

if [[ ! -f "${BASE_MANIFEST}" ]]; then
    echo "ERROR: manifest.tsv not found." >&2
    exit 1
fi

if [[ ! -f "${ROOT}/05_finalize.py" ]]; then
    echo "ERROR: 05_finalize.py not found." >&2
    exit 1
fi

if [[ ! -s "${EXTRA_GBK}" ]]; then
    echo "ERROR: curated MG1655 e14 GenBank not found:" >&2
    echo "  ${EXTRA_GBK}" >&2
    echo "Run the stage 04 first." >&2
    exit 1
fi

if [[ ! -s "${SESSION}" ]]; then
    echo "ERROR: Clinker session not found:" >&2
    echo "  ${SESSION}" >&2
    echo "Run the stage 04 first." >&2
    exit 1
fi

if ! conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "${ENV}"; then
    echo "ERROR: conda environment '${ENV}' does not exist." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build a manifest including e14
#
# MG1655_curated_e14 is inserted immediately after Ecoli_K12_e14 so that:
#
#   Ecoli_K12_e14
#   MG1655_curated_e14
#   ...
#
# are adjacent in the fixed-order Clinker plot.
#
# This generated manifest is used ONLY for this augmented analysis. The main
# manifest.tsv is not changed.
# ---------------------------------------------------------------------------

export BASE_MANIFEST AUG_MANIFEST EXTRA_LABEL E14_LABEL

conda run --no-capture-output -n "${ENV}" python 05_augment_manifest_MG1655.py

# ---------------------------------------------------------------------------
# Build a directory containing exactly the GenBanks represented by the
# augmented Clinker session.
#
# The seven pipeline GenBanks are symlinked from annotation/final_gbk.
# The original MG1655 annotation is copied in as MG1655_curated_e14.gbk.
# ---------------------------------------------------------------------------

rm -rf "${INPUT_GBK_DIR}"
mkdir -p "${INPUT_GBK_DIR}"

while IFS=$'\t' read -r label rest; do
    [[ -n "${label}" ]] || continue

    if [[ "${label}" == "${EXTRA_LABEL}" ]]; then
        cp -f "${EXTRA_GBK}" "${INPUT_GBK_DIR}/${label}.gbk"
    else
        SRC="${BASE_GBK_DIR}/${label}.gbk"

        if [[ ! -s "${SRC}" ]]; then
            echo "ERROR: missing final annotated GenBank: ${SRC}" >&2
            exit 1
        fi

        ln -s "$(realpath "${SRC}")" "${INPUT_GBK_DIR}/${label}.gbk"
    fi
done < <(tail -n +2 "${AUG_MANIFEST}")

# ---------------------------------------------------------------------------
# Review / apply
# ---------------------------------------------------------------------------

case "${MODE}" in
    review)
        echo "=== Stage 05 + MG1655: homology-group review ==="

        conda run --no-capture-output -n "${ENV}" \
            python "${ROOT}/05_finalize.py" review \
            --root "${ROOT}" \
            --manifest "${OUTDIR}/manifest_with_MG1655_annotation.tsv" \
            --session "${SESSION}" \
            --outdir "${OUTDIR}"

        echo
        echo "Review table:"
        echo "  ${REVIEW}"
        echo
        echo "Inspect rows with review_required=1."
        echo "Set final_product and approved=yes for those rows."
        echo
        echo "Then run:"
        echo "  MODE=apply ./05_finalize.sh"
        ;;

    apply)
        echo "=== Stage 05 + MG1655: apply reviewed labels ==="

        if [[ ! -s "${REVIEW}" ]]; then
            echo "ERROR: review table not found: ${REVIEW}" >&2
            echo "Run ./05_finalize.sh first." >&2
            exit 1
        fi

        conda run --no-capture-output -n "${ENV}" \
            python "${ROOT}/05_finalize.py" apply \
            --root "${ROOT}" \
            --manifest "${OUTDIR}/manifest_with_MG1655_annotation.tsv" \
            --session "${SESSION}" \
            --outdir "${OUTDIR}" \
            --review-table "${OUTDIR}/group_label_review.tsv" \
            --input-gbk-dir "${OUTDIR}/input_gbk" \
            --output-gbk-dir "${OUTDIR}/harmonised_gbk"

        echo
        echo "=== Rebuilding final Clinker from harmonised GenBanks ==="

        GBKS=()
        while IFS=$'\t' read -r label rest; do
            [[ -n "${label}" ]] || continue
            GBK="${HARMONISED}/${label}.gbk"

            if [[ ! -s "${GBK}" ]]; then
                echo "ERROR: missing harmonised GenBank: ${GBK}" >&2
                exit 1
            fi

            GBKS+=("${GBK}")
        done < <(tail -n +2 "${AUG_MANIFEST}")

        echo
        echo "Sequences in final Clinker:"
        printf '  %s\n' "${GBKS[@]##*/}"
        echo

        # ---------------------------------------------------------------
        # Manifest-order final plot
        # ---------------------------------------------------------------

        rm -f "${MANIFEST_SESSION}"

        conda run --no-capture-output -n "${ENV}" \
            clinker \
            "${GBKS[@]}" \
            --identity "${CLINKER_IDENTITY}" \
            --jobs "${THREADS}" \
            --session "${MANIFEST_SESSION}" \
            --output "${CLINKER_OUT}/final_manifest_order_alignments.csv" \
            --delimiter "," \
            --decimals 4 \
            --matrix_out "${CLINKER_OUT}/final_manifest_order_matrix.csv" \
            --use_file_order \
            --plot "${CLINKER_OUT}/final_manifest_order.html" \
            --force \
            2>&1 | tee "${ROOT}/logs/05_MG1655_clinker_manifest.log"

        # ---------------------------------------------------------------
        # Auto-order final plot
        # ---------------------------------------------------------------

        rm -f "${AUTO_SESSION}"

        conda run --no-capture-output -n "${ENV}" \
            clinker \
            "${GBKS[@]}" \
            --identity "${CLINKER_IDENTITY}" \
            --jobs "${THREADS}" \
            --session "${AUTO_SESSION}" \
            --output "${CLINKER_OUT}/final_auto_order_alignments.csv" \
            --delimiter "," \
            --decimals 4 \
            --matrix_out "${CLINKER_OUT}/final_auto_order_matrix.csv" \
            --plot "${CLINKER_OUT}/final_auto_order.html" \
            --force \
            2>&1 | tee "${ROOT}/logs/05_MG1655_clinker_auto.log"

        # ---------------------------------------------------------------
        # Relabel EXISTING Clinker homology groups using reviewed products.
        #
        # This is copied from the working original 05_finalize.sh logic.
        # It changes group["label"] in the saved session JSON only.
        #
        # Crucially:
        #   - group membership is NOT changed
        #   - no new functional groups are created
        #   - singleton/unlinked genes remain singleton/unlinked
        #   - colours therefore remain restricted to Clinker's actual
        #     sequence-defined homology groups
        # ---------------------------------------------------------------

        echo
        echo "=== Relabelling Clinker groups using reviewed products ==="

        export ROOT OUTDIR REVIEW GROUPS_TSV CLINKER_OUT
        export MANIFEST_SESSION AUTO_SESSION
        export MANIFEST_PRODUCT_SESSION AUTO_PRODUCT_SESSION

        conda run --no-capture-output -n "${ENV}" python 05_relabel_clinker_reviewed.py

        # ---------------------------------------------------------------
        # Render the relabelled sessions.
        #
        # No GenBanks are supplied here, so these commands only render the
        # existing sessions and cannot add or recolour genes.
        # ---------------------------------------------------------------

        echo
        echo "=== Rendering product-labelled Clinker plots ==="

        rm -f \
            "${CLINKER_OUT}/final_manifest_order_product_labels.html" \
            "${CLINKER_OUT}/final_auto_order_product_labels.html"

        conda run --no-capture-output -n "${ENV}" \
            clinker \
            --session "${MANIFEST_PRODUCT_SESSION}" \
            --plot "${CLINKER_OUT}/final_manifest_order_product_labels.html" \
            --force \
            2>&1 | tee "${ROOT}/logs/05_MG1655_clinker_manifest_product_labels.log"

        conda run --no-capture-output -n "${ENV}" \
            clinker \
            --session "${AUTO_PRODUCT_SESSION}" \
            --plot "${CLINKER_OUT}/final_auto_order_product_labels.html" \
            --force \
            2>&1 | tee "${ROOT}/logs/05_MG1655_clinker_auto_product_labels.log"

        echo
        echo "======================================================================"
        echo "Stage 05 + MG1655 complete."
        echo
        echo "Final tables:"
        echo "  ${OUTDIR}/homologous_gene_groups.tsv"
        echo "  ${OUTDIR}/gene_membership.tsv"
        echo "  ${OUTDIR}/gene_presence_absence.tsv"
        echo "  ${OUTDIR}/gene_copy_number.tsv"
        echo "  ${OUTDIR}/harmonisation_changes.tsv"
        echo
        echo "Final harmonised GenBanks:"
        echo "  ${HARMONISED}/"
        echo
        echo "Final Clinker:"
        echo "  ${CLINKER_OUT}/final_manifest_order.html"
        echo "  ${CLINKER_OUT}/final_auto_order.html"
        echo
        echo "Product-labelled Clinker:"
        echo "  ${CLINKER_OUT}/final_manifest_order_product_labels.html"
        echo "  ${CLINKER_OUT}/final_auto_order_product_labels.html"
        echo "Next step:"
	echo "  ./06_phylogeny.sh"
        ;;

    *)
        echo "ERROR: MODE must be 'review' or 'apply'." >&2
        exit 1
        ;;
esac
