#!/usr/bin/env bash
set -euo pipefail

# 06_phylogeny.sh
#
# 1. Search the installed INPHARED Mash database using manifest discovery seeds.
# 2. Download references passing the distance threshold; remove near-duplicates
#    with pairwise Mash distances and cap the retained reference set.
# 3. Predict proteins with PHANOTATE and cluster them with MMseqs2.
# 4. Select single-copy protein groups containing e14, align with MAFFT,
#    trim with ClipKIT, and infer one partitioned ML tree with IQ-TREE.
#
# Usage (from the project root, with both Stage 06 files installed there):
#   chmod +x 06_phylogeny.sh
#   ./06_phylogeny.sh
#
# Useful overrides:
#   THREADS=16 ./06_phylogeny.sh
#   MAX_MASH_DISTANCE=0.10 MAX_TOTAL=100 ./06_phylogeny.sh
#   E14_CORE_FRACTION=0.30 E14_MIN_TAXON_GROUPS=4 ./06_phylogeny.sh
#   REDUNDANCY_MASH=0.005 ./06_phylogeny.sh
#   FORCE=1 ./06_phylogeny.sh
#
# Defaults: distance <= 0.15; near-duplicate distance < 0.002;
# maximum 150 focal + reference genomes (all focal genomes are retained);
# protein-group occupancy >= 0.25; >= 4 groups; >= 3 groups per tree taxon.
# Occupancy is measured across all retained genomes before taxon filtering.
#
# Requirements: existing 'pharokka' conda environment and its INPHARED database.
# A separate 'prophage_phylogeny' environment is created if absent.
# Cached Mash searches and protein predictions are reused unless FORCE=1.
# Clusters, alignments and the protein tree are rebuilt on every run.
# Previous Stage 06 outputs for removed analyses are left untouched.

THREADS="${THREADS:-8}"
MAX_MASH_DISTANCE="${MAX_MASH_DISTANCE:-0.15}"

# Only these manifest type(s) are used to discover additional INPHARED
# neighbours. All manifest genomes remain predefined downstream taxa.
# Supply multiple types comma-separated if ever needed.
DISCOVERY_TYPES="${DISCOVERY_TYPES:-prophage}"
MAX_TOTAL="${MAX_TOTAL:-150}"
REDUNDANCY_MASH="${REDUNDANCY_MASH:-0.002}"

# Space-separated accessions known not to be retrievable from NCBI under the
# identifier used in the installed INPHARED metadata. They remain visible in
# the threshold-hit/download tables but are explicitly logged and skipped.
SKIP_DOWNLOAD_ACCESSIONS="${SKIP_DOWNLOAD_ACCESSIONS:-SFU82619}"
MASH_K="${MASH_K:-21}"
MASH_SKETCH="${MASH_SKETCH:-1000}"

PROTEIN_MIN_ID="${PROTEIN_MIN_ID:-0.30}"
PROTEIN_COV="${PROTEIN_COV:-0.50}"
# e14-inclusive partial backbone: e14 must be present, but lower occupancy is
# expected because e14 is strongly reduced; the group count depends on the
# current retained taxa and protein clusters.
E14_CORE_FRACTION="${E14_CORE_FRACTION:-0.25}"
# Require each taxon in the e14-inclusive concatenated tree to contribute to
# at least this many selected protein groups after trimming. This limits
# sparse taxon representation; inspect coverage before interpreting the tree.
E14_MIN_TAXON_GROUPS="${E14_MIN_TAXON_GROUPS:-3}"

MIN_PROTEIN_AA="${MIN_PROTEIN_AA:-50}"
MIN_CORE_GROUPS="${MIN_CORE_GROUPS:-4}"
ANCHOR_LABEL="${ANCHOR_LABEL:-Ecoli_K12_e14}"


FORCE="${FORCE:-0}"

ROOT="$(pwd)"
MANIFEST="${MANIFEST:-${ROOT}/manifest.tsv}"
HELPER="${ROOT}/06_phylogeny.py"

PHYLO_ENV="${PHYLO_ENV:-prophage_phylogeny}"
PHAROKKA_ENV="${PHAROKKA_ENV:-pharokka}"

PHAROKKA_DB="${PHAROKKA_DB:-${ROOT}/databases/pharokka}"

WORK="${ROOT}/06_phylogeny"
RESULTS="${ROOT}/results/phylogeny"
LOGS="${ROOT}/logs"

RAW_HITS="${WORK}/inphared/raw_mash"
DOWNLOADS="${WORK}/inphared/downloads"
GENOMES_ALL="${WORK}/genomes_all"
GENOMES_KEEP="${WORK}/genomes_retained"
PROTEINS="${WORK}/proteins"
MMSEQS="${WORK}/mmseqs"
BACKBONE_ROOT="${RESULTS}/backbone"
BACKBONE_E14="${BACKBONE_ROOT}/e14_inclusive"

mkdir -p \
    "${RAW_HITS}" "${DOWNLOADS}" "${GENOMES_ALL}" "${GENOMES_KEEP}" \
    "${PROTEINS}" "${MMSEQS}" "${BACKBONE_E14}" \
    "${LOGS}"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2
    exit 1
fi
if [[ ! -f "${HELPER}" ]]; then
    echo "ERROR: helper not found: ${HELPER}" >&2
    exit 1
fi

env_exists () {
    conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$1"
}

if ! env_exists "${PHAROKKA_ENV}"; then
    echo "ERROR: existing Pharokka environment '${PHAROKKA_ENV}' not found." >&2
    exit 1
fi

if ! env_exists "${PHYLO_ENV}"; then
    echo "=== Creating ${PHYLO_ENV} ==="
    conda create -y -n "${PHYLO_ENV}" \
        -c conda-forge -c bioconda \
        mash mmseqs2 mafft clipkit iqtree entrez-direct \
        biopython pandas
fi

# Verify central tools.
for cmd in mash mmseqs mafft clipkit efetch; do
    if ! conda run -n "${PHYLO_ENV}" bash -c "command -v ${cmd} >/dev/null"; then
        echo "ERROR: '${cmd}' missing from ${PHYLO_ENV}" >&2
        exit 1
    fi
done

if conda run -n "${PHYLO_ENV}" iqtree2 -version >/dev/null 2>&1; then
    IQTREE="iqtree2"
elif conda run -n "${PHYLO_ENV}" iqtree -version >/dev/null 2>&1; then
    IQTREE="iqtree"
else
    echo "ERROR: IQ-TREE executable not found in ${PHYLO_ENV}" >&2
    exit 1
fi

if ! conda run -n "${PHAROKKA_ENV}" bash -c "command -v phanotate.py >/dev/null"; then
    echo "ERROR: phanotate.py not found in ${PHAROKKA_ENV}" >&2
    exit 1
fi

# Locate the exact INPHARED files installed with Pharokka.
mapfile -t MSH_FILES < <(
    find "${PHAROKKA_DB}" -maxdepth 1 -type f         -name '*genomes.fa.msh' ! -name '._*' | sort
)
mapfile -t DATA_FILES < <(
    find "${PHAROKKA_DB}" -maxdepth 1 -type f         -name '*data.tsv' ! -name '._*' | sort
)

if (( ${#MSH_FILES[@]} != 1 )); then
    echo "ERROR: expected exactly one INPHARED Mash sketch in ${PHAROKKA_DB}" >&2
    printf '  %s\n' "${MSH_FILES[@]:-}" >&2
    exit 1
fi
if (( ${#DATA_FILES[@]} != 1 )); then
    echo "ERROR: expected exactly one INPHARED metadata TSV in ${PHAROKKA_DB}" >&2
    printf '  %s\n' "${DATA_FILES[@]:-}" >&2
    exit 1
fi

INPHARED_MSH="${MSH_FILES[0]}"
INPHARED_META="${DATA_FILES[0]}"

echo "INPHARED Mash sketch: ${INPHARED_MSH}"
echo "INPHARED metadata:    ${INPHARED_META}"
echo

# ---------------------------------------------------------------------------
# 1. Query discovery-seed sequences against INPHARED
# ---------------------------------------------------------------------------

echo "=== 1. INPHARED Mash searches from discovery seeds only ==="
mkdir -p "${WORK}/queries"

DISCOVERY_SEEDS="${RESULTS}/discovery_seeds.tsv"
{
    printf 'label\torganism\taccession\ttype\n'
    awk -F'\t' -v types="${DISCOVERY_TYPES}" '
        BEGIN {
            n=split(types, a, ",")
            for (i=1; i<=n; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", a[i])
                wanted[a[i]]=1
            }
        }
        NR>1 && ($4 in wanted) {
            print $1 "\t" $2 "\t" $3 "\t" $4
        }
    ' "${MANIFEST}"
} > "${DISCOVERY_SEEDS}"

N_SEEDS="$(tail -n +2 "${DISCOVERY_SEEDS}" | wc -l | tr -d ' ')"
if (( N_SEEDS < 1 )); then
    echo "ERROR: no manifest rows match DISCOVERY_TYPES=${DISCOVERY_TYPES}" >&2
    exit 1
fi

echo "Discovery seed types: ${DISCOVERY_TYPES}"
echo "Discovery seeds (${N_SEEDS}):"
tail -n +2 "${DISCOVERY_SEEDS}" | cut -f1 | sed 's/^/  /'

while IFS=$'\t' read -r label organism accession type; do
    [[ -n "${label}" ]] || continue

    SRC="${ROOT}/prophages/${label}/${label}.fna"
    QUERY="${WORK}/queries/${label}.fna"
    SKETCH="${WORK}/queries/${label}.msh"
    RAW="${RAW_HITS}/${label}.mash.tsv"

    if [[ ! -s "${SRC}" ]]; then
        echo "ERROR: missing focal FASTA: ${SRC}" >&2
        echo "Run stages 01-03 first for every manifest entry." >&2
        exit 1
    fi

    conda run -n "${PHYLO_ENV}" python "${HELPER}" normalise-fasta \
        --input "${SRC}" --output "${QUERY}" --id "${label}"

    if [[ "${FORCE}" == "1" || ! -s "${RAW}" ]]; then
        rm -f "${SKETCH}"
        conda run -n "${PHYLO_ENV}" mash sketch \
            -i -k "${MASH_K}" -s "${MASH_SKETCH}" \
            -o "${WORK}/queries/${label}" "${QUERY}"

        conda run -n "${PHYLO_ENV}" mash dist \
            "${INPHARED_MSH}" "${SKETCH}" > "${RAW}"
    else
        echo "  [reuse] ${RAW}"
    fi
done < <(tail -n +2 "${DISCOVERY_SEEDS}")

INPHARED_HITS="${RESULTS}/inphared_hits_within_threshold.tsv"
DOWNLOAD_LIST="${WORK}/inphared/download_list.tsv"

conda run -n "${PHYLO_ENV}" python "${HELPER}" parse-hits \
    --manifest "${MANIFEST}" \
    --metadata "${INPHARED_META}" \
    --raw-dir "${RAW_HITS}" \
    --out "${INPHARED_HITS}" \
    --download-list "${DOWNLOAD_LIST}" \
    --max-distance "${MAX_MASH_DISTANCE}" \
    --discovery-types "${DISCOVERY_TYPES}"

echo "All INPHARED hits within Mash distance <= ${MAX_MASH_DISTANCE}:"
echo "  ${INPHARED_HITS}"
echo "Unique non-focal accessions queued for download:"
echo "  ${DOWNLOAD_LIST}"
echo

# ---------------------------------------------------------------------------
# 2. Download unique INPHARED references
# ---------------------------------------------------------------------------

echo "=== 2. Downloading INPHARED reference genomes ==="

tail -n +2 "${DOWNLOAD_LIST}" | cut -f1 | \
while IFS=$'\t' read -r accession; do
    [[ -n "${accession}" ]] || continue

    skip=0
    for blocked in ${SKIP_DOWNLOAD_ACCESSIONS}; do
        if [[ "${accession}" == "${blocked}" ]]; then
            echo "  [skip] ${accession} (configured as unretrievable)"
            skip=1
            break
        fi
    done
    [[ "${skip}" == "1" ]] && continue

    RAW="${DOWNLOADS}/${accession}.fna"

    if [[ "${FORCE}" != "1" && -s "${RAW}" ]]; then
        echo "  [reuse] ${accession}"
        continue
    fi

    echo "  [download] ${accession}"
    tmp="${RAW}.tmp"
    rm -f "${tmp}"

    ok=0
    for attempt in 1 2 3; do
        if conda run -n "${PHYLO_ENV}" efetch \
            -db nuccore -id "${accession}" -format fasta > "${tmp}" \
            && grep -q '^>' "${tmp}"; then
            ok=1
            break
        fi
        sleep 2
    done

    if [[ "${ok}" == "1" ]]; then
        mv "${tmp}" "${RAW}"
    else
        echo "WARNING: failed to download ${accession}; skipping." >&2
        rm -f "${tmp}"
    fi
done

# ---------------------------------------------------------------------------
# 3. Build candidate genome set and lightly de-replicate
# ---------------------------------------------------------------------------

echo
echo "=== 3. Downloaded-set de-replication and final selection ==="
echo "Selection order:"
echo "  1. all unique INPHARED hits within Mash distance <= ${MAX_MASH_DISTANCE}"
echo "     of discovery seeds only (${DISCOVERY_TYPES})"
echo "  2. all-vs-all Mash comparison"
echo "  3. remove near-duplicates at Mash distance < ${REDUNDANCY_MASH}"
echo "  4. rank surviving references by distance to the discovery-seed set"
echo "  5. apply final MAX_TOTAL=${MAX_TOTAL} cap"
rm -f "${GENOMES_ALL}"/*.fna

# Focal genomes always retained.
while IFS=$'\t' read -r label accession type; do
    [[ -n "${label}" ]] || continue
    SRC="${ROOT}/prophages/${label}/${label}.fna"
    conda run -n "${PHYLO_ENV}" python "${HELPER}" normalise-fasta \
        --input "${SRC}" \
        --output "${GENOMES_ALL}/${label}.fna" \
        --id "${label}"
done < <(tail -n +2 "${MANIFEST}" | cut -f1,3,4)

# Downloaded INPHARED references.
tail -n +2 "${DOWNLOAD_LIST}" | cut -f1 | \
while IFS=$'\t' read -r accession; do
    [[ -n "${accession}" ]] || continue
    SRC="${DOWNLOADS}/${accession}.fna"
    [[ -s "${SRC}" ]] || continue
    taxon="INPHARED_${accession%%.*}"
    conda run -n "${PHYLO_ENV}" python "${HELPER}" normalise-fasta \
        --input "${SRC}" \
        --output "${GENOMES_ALL}/${taxon}.fna" \
        --id "${taxon}"
done

ALL_TAXA="${WORK}/all_taxa.tsv"
conda run -n "${PHYLO_ENV}" python "${HELPER}" make-taxa \
    --manifest "${MANIFEST}" \
    --download-list "${DOWNLOAD_LIST}" \
    --genome-dir "${GENOMES_ALL}" \
    --out "${ALL_TAXA}"

ALL_FASTA="${WORK}/all_candidates.fna"
cat "${GENOMES_ALL}"/*.fna > "${ALL_FASTA}"

ALL_SKETCH_PREFIX="${WORK}/all_candidates"
ALL_SKETCH="${ALL_SKETCH_PREFIX}.msh"
ALL_PAIRWISE="${WORK}/all_candidates_pairwise.tsv"

# Candidate composition may change whenever the threshold-hit set changes, so
# rebuild this inexpensive all-vs-all Mash comparison rather than reusing a
# potentially stale matrix.
rm -f "${ALL_SKETCH}" "${ALL_PAIRWISE}"
conda run -n "${PHYLO_ENV}" mash sketch \
    -i -k "${MASH_K}" -s "${MASH_SKETCH}" \
    -o "${ALL_SKETCH_PREFIX}" "${ALL_FASTA}"

conda run -n "${PHYLO_ENV}" mash dist \
    "${ALL_SKETCH}" "${ALL_SKETCH}" > "${ALL_PAIRWISE}"

RETAINED="${RESULTS}/retained_taxa.tsv"
DROPPED="${RESULTS}/dropped_near_duplicates.tsv"

conda run -n "${PHYLO_ENV}" python "${HELPER}" curate-genomes \
    --pairwise "${ALL_PAIRWISE}" \
    --taxa "${ALL_TAXA}" \
    --out "${RETAINED}" \
    --dropped "${DROPPED}" \
    --redundancy-distance "${REDUNDANCY_MASH}" \
    --max-total "${MAX_TOTAL}"

rm -f "${GENOMES_KEEP}"/*.fna
tail -n +2 "${RETAINED}" | cut -f1 | while read -r taxon; do
    cp "${GENOMES_ALL}/${taxon}.fna" "${GENOMES_KEEP}/${taxon}.fna"
done

N_RETAINED="$(tail -n +2 "${RETAINED}" | wc -l | tr -d ' ')"
echo "Retained genomes: ${N_RETAINED}"
echo

# ---------------------------------------------------------------------------
# 4. Consistent PHANOTATE calls for the retained genomes
# ---------------------------------------------------------------------------

echo
echo "=== 4. PHANOTATE protein prediction ==="

rm -f "${PROTEINS}/all_proteins.faa"

tail -n +2 "${RETAINED}" | cut -f1 | while read -r taxon; do
    FNA="${GENOMES_KEEP}/${taxon}.fna"
    RAWFAA="${PROTEINS}/${taxon}.raw.faa"
    FAA="${PROTEINS}/${taxon}.faa"

    if [[ "${FORCE}" == "1" || ! -s "${FAA}" ]]; then
        echo "  [PHANOTATE] ${taxon}"
        conda run --no-capture-output -n "${PHAROKKA_ENV}" \
            phanotate.py "${FNA}" -f faa > "${RAWFAA}"

        conda run -n "${PHYLO_ENV}" python "${HELPER}" sanitise-proteins \
            --input "${RAWFAA}" \
            --output "${FAA}" \
            --taxon "${taxon}"
    fi
done

# Build the combined FASTA once for human inspection/provenance.
{
    while IFS=$'\t' read -r taxon; do
        [[ -n "${taxon}" ]] || continue
        cat "${PROTEINS}/${taxon}.faa"
    done < <(tail -n +2 "${RETAINED}" | cut -f1)
} > "${PROTEINS}/all_proteins.faa"

# ---------------------------------------------------------------------------
# 5. MMseqs2 protein clustering + e14-anchored backbone
# ---------------------------------------------------------------------------

echo
echo "=== 5. e14-inclusive protein clustering ==="

# Rebuild a canonical FASTA on every run from the CURRENT retained taxa.
# This prevents stale clustering and avoids parser sensitivity to upstream
# FASTA wrapping/formatting.
MMSEQS_INPUT="${PROTEINS}/all_proteins_mmseqs.faa"
MMSEQS_INPUT_SUMMARY="${PROTEINS}/all_proteins_mmseqs_summary.tsv"

conda run --no-capture-output -n "${PHYLO_ENV}" \
    python "${HELPER}" prepare-mmseqs-fasta \
    --taxa "${RETAINED}" \
    --protein-dir "${PROTEINS}" \
    --output "${MMSEQS_INPUT}" \
    --summary "${MMSEQS_INPUT_SUMMARY}"

EXPECTED_PROTEINS="$(grep -c '^>' "${MMSEQS_INPUT}" || true)"
if (( EXPECTED_PROTEINS < 1 )); then
    echo "ERROR: no proteins in canonical MMseqs FASTA." >&2
    exit 1
fi

echo "Canonical protein count: ${EXPECTED_PROTEINS}"

# Always rebuild MMseqs from the canonical input. This step is cheap for this
# dataset and is safer than reusing a cluster file after the taxon set changes.
rm -rf "${MMSEQS}"
mkdir -p "${MMSEQS}"

MMSEQS_DB="${MMSEQS}/proteins_db"
MMSEQS_CHECK="${MMSEQS}/proteins_db_check.faa"
MMSEQS_CLU="${MMSEQS}/proteins_clu"
MMSEQS_TMP="${MMSEQS}/tmp"
MMSEQS_CLUSTER="${MMSEQS}/proteins_cluster.tsv"

echo "  [MMseqs] createdb"
conda run --no-capture-output -n "${PHYLO_ENV}" \
    mmseqs createdb \
    "${MMSEQS_INPUT}" \
    "${MMSEQS_DB}" \
    --dbtype 1 \
    --shuffle 0 \
    2>&1 | tee "${LOGS}/06_mmseqs.log"

# Hard validation: MMseqs must parse every supplied protein.
conda run --no-capture-output -n "${PHYLO_ENV}" \
    mmseqs convert2fasta \
    "${MMSEQS_DB}" \
    "${MMSEQS_CHECK}" \
    >> "${LOGS}/06_mmseqs.log" 2>&1

PARSED_PROTEINS="$(grep -c '^>' "${MMSEQS_CHECK}" || true)"
if [[ "${PARSED_PROTEINS}" -ne "${EXPECTED_PROTEINS}" ]]; then
    echo "ERROR: MMseqs parsed ${PARSED_PROTEINS}/${EXPECTED_PROTEINS} proteins." >&2
    echo "Input:  ${MMSEQS_INPUT}" >&2
    echo "Check:  ${MMSEQS_CHECK}" >&2
    exit 1
fi
echo "  parser check: ${PARSED_PROTEINS}/${EXPECTED_PROTEINS}"

echo "  [MMseqs] cluster"
conda run --no-capture-output -n "${PHYLO_ENV}" \
    mmseqs cluster \
    "${MMSEQS_DB}" \
    "${MMSEQS_CLU}" \
    "${MMSEQS_TMP}" \
    --min-seq-id "${PROTEIN_MIN_ID}" \
    -c "${PROTEIN_COV}" \
    --cov-mode 0 \
    --cluster-mode 1 \
    --threads "${THREADS}" \
    2>&1 | tee -a "${LOGS}/06_mmseqs.log"

conda run --no-capture-output -n "${PHYLO_ENV}" \
    mmseqs createtsv \
    "${MMSEQS_DB}" \
    "${MMSEQS_DB}" \
    "${MMSEQS_CLU}" \
    "${MMSEQS_CLUSTER}" \
    >> "${LOGS}/06_mmseqs.log" 2>&1

CLUSTERED_PROTEINS="$(wc -l < "${MMSEQS_CLUSTER}" | tr -d ' ')"
if [[ "${CLUSTERED_PROTEINS}" -ne "${EXPECTED_PROTEINS}" ]]; then
    echo "ERROR: MMseqs cluster table has ${CLUSTERED_PROTEINS} memberships" >&2
    echo "for ${EXPECTED_PROTEINS} input proteins." >&2
    exit 1
fi
echo "  cluster membership check: ${CLUSTERED_PROTEINS}/${EXPECTED_PROTEINS}"

# ---------------------------------------------------------------------------
# 6. e14-inclusive protein phylogeny
# ---------------------------------------------------------------------------

run_e14_analysis () {
    local OUTDIR="${BACKBONE_E14}"
    local TAG="e14_backbone"
    local IQ_PREFIX="${OUTDIR}/backbone_iqtree"

    # Invalidate previous tree outputs before rebuilding its inputs.
    rm -f "${IQ_PREFIX}".* "${OUTDIR}/backbone_supermatrix.faa" \
        "${OUTDIR}/backbone_partitions.nex" "${OUTDIR}/backbone_taxon_coverage.tsv"

    rm -rf "${OUTDIR}/groups" "${OUTDIR}/alignments"
    mkdir -p "${OUTDIR}/groups" "${OUTDIR}/alignments"

    local SELECT_ARGS=(
        python "${HELPER}" select-backbone
        --cluster-tsv "${MMSEQS_CLUSTER}"
        --proteins "${MMSEQS_INPUT}"
        --taxa "${RETAINED}"
        --min-fraction "${E14_CORE_FRACTION}"
        --min-aa "${MIN_PROTEIN_AA}"
        --min-groups "${MIN_CORE_GROUPS}"
        --outdir "${OUTDIR}"
        --anchor "${ANCHOR_LABEL}"
    )

    conda run --no-capture-output -n "${PHYLO_ENV}" "${SELECT_ARGS[@]}"

    local N_GROUPS
    N_GROUPS="$(tail -n +2 "${OUTDIR}/backbone_groups.tsv" | wc -l | tr -d ' ')"
    echo "  ${TAG}: ${N_GROUPS} protein groups"

    for faa in "${OUTDIR}"/groups/*.faa; do
        local gid
        gid="$(basename "${faa}" .faa)"
        local aln="${OUTDIR}/alignments/${gid}.mafft.faa"
        local trim="${OUTDIR}/alignments/${gid}.clipkit.faa"

        conda run --no-capture-output -n "${PHYLO_ENV}" \
            mafft --auto --thread 1 "${faa}" \
            > "${aln}" 2> "${LOGS}/06_${TAG}_${gid}_mafft.log"

        conda run --no-capture-output -n "${PHYLO_ENV}" \
            clipkit "${aln}" -m smart-gap -o "${trim}" \
            > "${LOGS}/06_${TAG}_${gid}_clipkit.log" 2>&1
    done

    local SUPERMATRIX="${OUTDIR}/backbone_supermatrix.faa"
    local PARTITIONS="${OUTDIR}/backbone_partitions.nex"
    local COVERAGE="${OUTDIR}/backbone_taxon_coverage.tsv"

    conda run --no-capture-output -n "${PHYLO_ENV}" \
        python "${HELPER}" concatenate \
        --align-dir "${OUTDIR}/alignments" \
        --taxa "${RETAINED}" \
        --groups "${OUTDIR}/backbone_groups.tsv" \
        --output "${SUPERMATRIX}" \
        --partitions "${PARTITIONS}" \
        --coverage "${COVERAGE}" \
        --min-taxon-groups "${E14_MIN_TAXON_GROUPS}" \
        --anchor "${ANCHOR_LABEL}" \
        --min-groups "${MIN_CORE_GROUPS}"

    conda run --no-capture-output -n "${PHYLO_ENV}" \
        "${IQTREE}" \
        -s "${SUPERMATRIX}" \
        -p "${PARTITIONS}" \
        -m MFP+MERGE \
        -B 1000 \
        --alrt 1000 \
        -T "${THREADS}" \
        --prefix "${IQ_PREFIX}" \
        2>&1 | tee "${LOGS}/06_${TAG}_iqtree.log"
}

echo
echo "=== 6. e14-inclusive protein phylogeny ==="
echo "Single-copy groups must contain ${ANCHOR_LABEL} and occur in"
echo ">= ${E14_CORE_FRACTION} of retained taxa."
run_e14_analysis

# ---------------------------------------------------------------------------
# 7. Provenance
# ---------------------------------------------------------------------------

VERSIONS="${RESULTS}/software_versions.txt"
{
    echo "Generated: $(date --iso-8601=seconds)"
    echo
    echo "mash:"
    conda run -n "${PHYLO_ENV}" mash --version 2>&1 || true
    echo
    echo "mmseqs:"
    conda run -n "${PHYLO_ENV}" mmseqs version 2>&1 || true
    echo
    echo "mafft:"
    conda run -n "${PHYLO_ENV}" mafft --version 2>&1 || true
    echo
    echo "clipkit:"
    conda run -n "${PHYLO_ENV}" clipkit --version 2>&1 || true
    echo
    echo "IQ-TREE:"
    conda run -n "${PHYLO_ENV}" "${IQTREE}" -version 2>&1 || true
    echo
    echo "PHANOTATE:"
    conda run -n "${PHAROKKA_ENV}" phanotate.py --version 2>&1 || true
    echo
    echo "INPHARED sketch:"
    basename "${INPHARED_MSH}"
    echo "INPHARED metadata:"
    basename "${INPHARED_META}"
    echo
    echo "Parameters:"
    echo "MAX_MASH_DISTANCE=${MAX_MASH_DISTANCE}"
    echo "DISCOVERY_TYPES=${DISCOVERY_TYPES}"
    echo "MAX_TOTAL=${MAX_TOTAL}"
    echo "REDUNDANCY_MASH=${REDUNDANCY_MASH}"
    echo "SKIP_DOWNLOAD_ACCESSIONS=${SKIP_DOWNLOAD_ACCESSIONS}"
    echo "MASH_K=${MASH_K}"
    echo "MASH_SKETCH=${MASH_SKETCH}"
    echo "PROTEIN_MIN_ID=${PROTEIN_MIN_ID}"
    echo "PROTEIN_COV=${PROTEIN_COV}"
    echo "E14_CORE_FRACTION=${E14_CORE_FRACTION}"
    echo "E14_MIN_TAXON_GROUPS=${E14_MIN_TAXON_GROUPS}"
    echo "MIN_PROTEIN_AA=${MIN_PROTEIN_AA}"
    echo "MIN_CORE_GROUPS=${MIN_CORE_GROUPS}"
    echo "ANCHOR_LABEL=${ANCHOR_LABEL}"
} > "${VERSIONS}"

echo
echo "======================================================================"
echo "Stage 06 complete."
echo
echo "INPHARED sampling:"
echo "  ${INPHARED_HITS}"
echo "  ${RETAINED}"
echo "  ${DROPPED}"
echo
echo "e14-inclusive partial-backbone ML phylogeny:"
echo "  ${BACKBONE_E14}/backbone_groups.tsv"
echo "  ${BACKBONE_E14}/backbone_taxon_coverage.tsv"
echo "  ${BACKBONE_E14}/backbone_supermatrix.faa"
echo "  ${BACKBONE_E14}/backbone_iqtree.treefile"
echo
echo "Provenance:"
echo "  ${VERSIONS}"
echo
echo "Interpretation:"
echo "  - Mash is used for reference discovery and near-duplicate removal only."
echo "  - IQ-TREE estimates relationships from e14-containing protein groups."
echo "  - Inspect backbone_groups.tsv and backbone_taxon_coverage.tsv."
echo "  - Sequence clusters are candidate orthologue groups; review their evidence."
echo "  - The tree is unrooted; e14 is a selection anchor, not an outgroup."
echo "Analysis finished."
echo "======================================================================"

