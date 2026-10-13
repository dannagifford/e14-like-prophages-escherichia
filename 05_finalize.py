#!/usr/bin/env python3
"""
05_finalize.py

Review and optionally harmonise Clinker homology-group annotations.

This stage deliberately separates *sequence-defined homology groups* from
*functional labels*. Clinker groups are single-linkage groups derived from
protein-sequence links above the identity threshold used in stage 04; they are
not treated here as formal orthogroups.

Modes
-----
review
    Read the manifest-order Clinker session and produce:
      - homologous_gene_groups.tsv
      - gene_membership.tsv
      - gene_presence_absence.tsv
      - gene_copy_number.tsv
      - group_label_review.tsv

    No annotation is modified.

apply
    Read an edited group_label_review.tsv, require approval for any groups with
    conflicting product labels, and write harmonised GenBank files in which
    CDS /product qualifiers are made consistent within each Clinker homology
    group. Original products are retained in /note qualifiers.

Typical workflow
----------------
    python 05_finalize.py review
    # inspect/edit results/final/group_label_review.tsv
    python 05_finalize.py apply
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from clinker.align import Globaligner


UNINFORMATIVE_PATTERNS = (
    r"^hypothetical protein$",
    r"^unknown protein$",
    r"^unknown function$",
    r"^uncharacteri[sz]ed protein$",
    r"^conserved hypothetical protein$",
    r"^putative protein$",
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["review", "apply"])
    p.add_argument("--manifest", default="manifest.tsv")
    p.add_argument(
        "--session",
        default="comparisons/clinker/e14_relatives_manifest_order_session.json",
    )
    p.add_argument("--root", default=".")
    p.add_argument("--outdir", default="results/final")
    p.add_argument(
        "--review-table",
        default="results/final/group_label_review.tsv",
    )
    p.add_argument(
        "--input-gbk-dir",
        default="annotation/final_gbk",
    )
    p.add_argument(
        "--output-gbk-dir",
        default="results/final/harmonised_gbk",
    )
    return p.parse_args()


def read_manifest(path: Path):
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    labels = [r["label"] for r in rows if r.get("label")]
    return rows, labels


def load_globaligner(path: Path) -> Globaligner:
    with path.open() as fh:
        return Globaligner.from_json(fh)


def normalise_product(s: str) -> str:
    s = re.sub(r"\s+", " ", (s or "").strip())
    s = s.rstrip(".;")
    return s.lower()


def is_uninformative(product: str) -> bool:
    p = normalise_product(product)
    if not p:
        return True
    return any(re.fullmatch(pattern, p) for pattern in UNINFORMATIVE_PATTERNS)


def first_name(gene, key: str) -> str:
    value = gene.names.get(key, "")
    if isinstance(value, list):
        return value[0] if value else ""
    return str(value) if value is not None else ""


def strand_int(value):
    try:
        return int(value)
    except Exception:
        return 0


def build_gene_catalogue(ga: Globaligner, manifest_labels: list[str]):
    """
    Return gene metadata keyed by Clinker gene UID.

    Cluster names should correspond to final GenBank stems (= manifest labels).
    """
    order_index = {label: i for i, label in enumerate(manifest_labels)}
    
    session_labels = {cluster.name for cluster in ga.clusters.values()}
    expected_labels = set(manifest_labels)

    missing = expected_labels - session_labels
    unexpected = session_labels - expected_labels

    if missing or unexpected:
        problems = []
        if missing:
            problems.append(
                "Manifest labels missing from Clinker session: "
                + ", ".join(sorted(missing))
            )
        if unexpected:
            problems.append(
                "Clinker clusters absent from manifest: "
                + ", ".join(sorted(unexpected))
            )

        raise ValueError(
            "\n".join(problems)
            + "\nRebuild the correct stage-04 Clinker session before stage 05."
        )
    
    genes = {}

    unknown_clusters = []

    for cluster in ga.clusters.values():
        cluster_name = cluster.name
        if cluster_name not in order_index:
            unknown_clusters.append(cluster_name)

        for locus in cluster.loci:
            for gene in locus.genes:
                genes[gene.uid] = {
                    "uid": gene.uid,
                    "cluster": cluster_name,
                    "cluster_order": order_index.get(cluster_name, 10**9),
                    "locus": locus.name,
                    "start": int(gene.start),
                    "end": int(gene.end),
                    "strand": strand_int(gene.strand),
                    "label": gene.label or "",
                    "locus_tag": first_name(gene, "locus_tag"),
                    "protein_id": first_name(gene, "protein_id"),
                    "gene_name": first_name(gene, "gene"),
                    "product": first_name(gene, "product") or "hypothetical protein",
                }

    if unknown_clusters:
        raise ValueError(
            "Clinker session contains cluster names absent from manifest: "
            + ", ".join(sorted(set(unknown_clusters)))
        )

    return genes


def build_groups(ga: Globaligner, genes: dict):
    """
    Build all Clinker homology groups plus singleton genes.

    Returns:
      groups: list[dict]
      uid_to_group_index: dict[gene_uid -> raw group index]
    """
    raw_groups = []
    assigned = set()

    for gr in ga.groups:
        members = [uid for uid in gr.genes if uid in genes]
        if not members:
            continue
        raw_groups.append(
            {
                "clinker_label": gr.label or "",
                "clinker_uid": gr.uid or "",
                "members": members,
                "singleton": False,
            }
        )
        assigned.update(members)

    for uid in genes:
        if uid not in assigned:
            raw_groups.append(
                {
                    "clinker_label": "singleton",
                    "clinker_uid": "",
                    "members": [uid],
                    "singleton": True,
                }
            )

    def group_sort_key(group):
        member_infos = [genes[uid] for uid in group["members"]]
        return min(
            (g["cluster_order"], g["start"], g["end"], g["label"])
            for g in member_infos
        )

    raw_groups.sort(key=group_sort_key)

    uid_to_group = {}
    for i, gr in enumerate(raw_groups, start=1):
        gr["group_id"] = f"G{i:03d}"
        for uid in gr["members"]:
            if uid in uid_to_group:
                raise ValueError(
                    f"Gene {uid} appears in more than one Clinker group"
                )
            uid_to_group[uid] = gr["group_id"]

    return raw_groups, uid_to_group


def link_identities_by_group(ga: Globaligner, uid_to_group: dict):
    values = defaultdict(list)

    for link in ga._links.values():
        q = link.query.uid
        t = link.target.uid
        gq = uid_to_group.get(q)
        gt = uid_to_group.get(t)

        if gq and gq == gt:
            values[gq].append(100.0 * float(link.identity))

    return values


def choose_product(products: list[str], singleton: bool):
    cleaned = [re.sub(r"\s+", " ", (p or "").strip()) for p in products]
    cleaned = [p if p else "hypothetical protein" for p in cleaned]

    if singleton:
        return cleaned[0], "singleton", 0

    informative = [p for p in cleaned if not is_uninformative(p)]
    informative_norm = [normalise_product(p) for p in informative]
    distinct_informative = sorted(set(informative_norm))

    if not informative:
        # All members are uninformative; use the most frequent original wording.
        count = Counter(cleaned)
        proposed = count.most_common(1)[0][0]
        return proposed, "uninformative", 0

    if len(distinct_informative) == 1:
        norm = distinct_informative[0]
        variants = [p for p in informative if normalise_product(p) == norm]
        proposed = Counter(variants).most_common(1)[0][0]

        if len(informative) == len(cleaned):
            return proposed, "consistent", 0
        return proposed, "single_informative_with_hypotheticals", 0

    # Do not silently resolve genuinely different functional labels.
    counts = Counter(informative_norm)
    winning_norm, _ = counts.most_common(1)[0]
    variants = [p for p in informative if normalise_product(p) == winning_norm]
    proposed = Counter(variants).most_common(1)[0][0]
    return proposed, "conflicting_products", 1


def summarise_groups(ga, groups, genes, uid_to_group, manifest_labels):
    identities = link_identities_by_group(ga, uid_to_group)

    summaries = []

    for gr in groups:
        gid = gr["group_id"]
        member_infos = [genes[uid] for uid in gr["members"]]

        products = [g["product"] for g in member_infos]
        proposed, status, review_required = choose_product(
            products, gr["singleton"]
        )

        present_genomes = [
            label
            for label in manifest_labels
            if any(g["cluster"] == label for g in member_infos)
        ]

        gene_names = sorted(
            {g["gene_name"] for g in member_infos if g["gene_name"]}
        )

        distinct_products = []
        seen = set()
        for p in products:
            n = normalise_product(p)
            if n not in seen:
                distinct_products.append(p)
                seen.add(n)

        vals = identities.get(gid, [])
        if vals:
            med_ident = statistics.median(vals)
            min_ident = min(vals)
            max_ident = max(vals)
        else:
            med_ident = min_ident = max_ident = None

        summaries.append(
            {
                "group_id": gid,
                "clinker_group_label": gr["clinker_label"],
                "n_genomes": len(set(g["cluster"] for g in member_infos)),
                "n_genes": len(member_infos),
                "genomes": ";".join(present_genomes),
                "proposed_product": proposed,
                "final_product": proposed,
                "status": status,
                "review_required": review_required,
                "approved": "" if review_required else "yes",
                "distinct_products": " | ".join(distinct_products),
                "gene_names": ";".join(gene_names),
                "median_link_identity_pct": (
                    f"{med_ident:.3f}" if med_ident is not None else ""
                ),
                "min_link_identity_pct": (
                    f"{min_ident:.3f}" if min_ident is not None else ""
                ),
                "max_link_identity_pct": (
                    f"{max_ident:.3f}" if max_ident is not None else ""
                ),
                "notes": (
                    "Manual product review required"
                    if review_required
                    else ""
                ),
            }
        )

    return summaries


def write_tsv(path: Path, fieldnames: list[str], rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)


def review_mode(args):
    root = Path(args.root).resolve()
    manifest_path = (root / args.manifest).resolve()
    session_path = (root / args.session).resolve()
    outdir = (root / args.outdir).resolve()

    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    if not session_path.exists():
        raise FileNotFoundError(
            f"{session_path}\nRun stage 04 before stage 05."
        )

    manifest_labels = read_manifest(manifest_path)
    ga = load_globaligner(session_path)

    genes = build_gene_catalogue(ga, manifest_labels)
    groups, uid_to_group = build_groups(ga, genes)
    summaries = summarise_groups(
        ga, groups, genes, uid_to_group, manifest_labels
    )
    summary_by_gid = {r["group_id"]: r for r in summaries}

    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Group summary
    # ------------------------------------------------------------------
    group_fields = [
        "group_id",
        "clinker_group_label",
        "n_genomes",
        "n_genes",
        "genomes",
        "proposed_product",
        "status",
        "review_required",
        "distinct_products",
        "gene_names",
        "median_link_identity_pct",
        "min_link_identity_pct",
        "max_link_identity_pct",
    ]
    write_tsv(
        outdir / "homologous_gene_groups.tsv",
        group_fields,
        summaries,
    )

    # ------------------------------------------------------------------
    # Manual review table
    # ------------------------------------------------------------------
    review_fields = [
        "group_id",
        "n_genomes",
        "n_genes",
        "genomes",
        "proposed_product",
        "final_product",
        "status",
        "review_required",
        "approved",
        "distinct_products",
        "gene_names",
        "median_link_identity_pct",
        "min_link_identity_pct",
        "max_link_identity_pct",
        "notes",
    ]
    write_tsv(
        outdir / "group_label_review.tsv",
        review_fields,
        summaries,
    )

    # ------------------------------------------------------------------
    # One row per CDS
    # ------------------------------------------------------------------
    membership_rows = []
    for gr in groups:
        gid = gr["group_id"]
        consensus = summary_by_gid[gid]["proposed_product"]

        member_infos = sorted(
            (genes[uid] for uid in gr["members"]),
            key=lambda g: (
                g["cluster_order"],
                g["start"],
                g["end"],
                g["label"],
            ),
        )

        for g in member_infos:
            membership_rows.append(
                {
                    "group_id": gid,
                    "genome": g["cluster"],
                    "locus": g["locus"],
                    "start_0based": g["start"],
                    "end_0based_exclusive": g["end"],
                    "strand": g["strand"],
                    "gene_label": g["label"],
                    "locus_tag": g["locus_tag"],
                    "protein_id": g["protein_id"],
                    "gene_name": g["gene_name"],
                    "product": g["product"],
                    "proposed_group_product": consensus,
                }
            )

    membership_fields = [
        "group_id",
        "genome",
        "locus",
        "start_0based",
        "end_0based_exclusive",
        "strand",
        "gene_label",
        "locus_tag",
        "protein_id",
        "gene_name",
        "product",
        "proposed_group_product",
    ]
    write_tsv(
        outdir / "gene_membership.tsv",
        membership_fields,
        membership_rows,
    )

    # ------------------------------------------------------------------
    # Presence/absence and copy-number matrices
    # ------------------------------------------------------------------
    p_fields = [
        "group_id",
        "proposed_product",
        "n_genomes",
        "n_genes",
        *manifest_labels,
    ]

    presence_rows = []
    copy_rows = []

    for gr in groups:
        gid = gr["group_id"]
        s = summary_by_gid[gid]
        by_genome = defaultdict(list)

        for uid in gr["members"]:
            g = genes[uid]
            by_genome[g["cluster"]].append(g["label"])

        presence = {
            "group_id": gid,
            "proposed_product": s["proposed_product"],
            "n_genomes": s["n_genomes"],
            "n_genes": s["n_genes"],
        }
        copies = dict(presence)

        for label in manifest_labels:
            labels = sorted(x for x in by_genome.get(label, []) if x)
            presence[label] = ";".join(labels)
            copies[label] = len(by_genome.get(label, []))

        presence_rows.append(presence)
        copy_rows.append(copies)

    write_tsv(
        outdir / "gene_presence_absence.tsv",
        p_fields,
        presence_rows,
    )
    write_tsv(
        outdir / "gene_copy_number.tsv",
        p_fields,
        copy_rows,
    )

    conflicting = sum(int(r["review_required"]) for r in summaries)

    print(f"Wrote {outdir / 'homologous_gene_groups.tsv'}")
    print(f"Wrote {outdir / 'gene_membership.tsv'}")
    print(f"Wrote {outdir / 'gene_presence_absence.tsv'}")
    print(f"Wrote {outdir / 'gene_copy_number.tsv'}")
    print(f"Wrote {outdir / 'group_label_review.tsv'}")
    print()
    print(f"Total CDS: {len(genes)}")
    print(f"Total Clinker groups including singletons: {len(groups)}")
    print(f"Groups requiring manual product review: {conflicting}")

    if conflicting:
        print()
        print(
            "Edit group_label_review.tsv: inspect rows with "
            "review_required=1, set final_product as desired, and set "
            "approved=yes before apply mode."
        )


def approved(value: str) -> bool:
    return str(value).strip().lower() in {"yes", "y", "true", "1", "approved"}


def read_review_table(path: Path):
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    required = {
        "group_id",
        "final_product",
        "review_required",
        "approved",
    }
    missing = required - set(rows[0].keys() if rows else [])
    if missing:
        raise ValueError(
            f"Review table missing columns: {', '.join(sorted(missing))}"
        )

    return rows


def add_note(feature, text: str):
    notes = feature.qualifiers.setdefault("note", [])
    if text not in notes:
        notes.append(text)


def apply_mode(args):
    root = Path(args.root).resolve()
    manifest_path = (root / args.manifest).resolve()
    session_path = (root / args.session).resolve()
    review_path = (root / args.review_table).resolve()
    input_gbk_dir = (root / args.input_gbk_dir).resolve()
    output_gbk_dir = (root / args.output_gbk_dir).resolve()

    for p in (manifest_path, session_path, review_path):
        if not p.exists():
            raise FileNotFoundError(p)

    _, manifest_labels = read_manifest(manifest_path)
    ga = load_globaligner(session_path)
    genes = build_gene_catalogue(ga, manifest_labels)
    groups, uid_to_group = build_groups(ga, genes)

    review_rows = read_review_table(review_path)
    review = {r["group_id"]: r for r in review_rows}

    expected_gids = {gr["group_id"] for gr in groups}
    missing_review = sorted(expected_gids - set(review))
    extra_review = sorted(set(review) - expected_gids)

    if missing_review:
        raise ValueError(
            "Review table is missing groups: " + ", ".join(missing_review)
        )
    if extra_review:
        raise ValueError(
            "Review table contains unknown groups: " + ", ".join(extra_review)
        )

    unresolved = []
    for gid in sorted(expected_gids):
        row = review[gid]
        final_product = (row.get("final_product") or "").strip()
        if not final_product:
            unresolved.append(f"{gid}: blank final_product")
        if str(row.get("review_required", "")).strip() in {"1", "true", "True"}:
            if not approved(row.get("approved", "")):
                unresolved.append(f"{gid}: manual approval required")

    if unresolved:
        msg = "\n".join(f"  - {x}" for x in unresolved)
        raise ValueError(
            "Cannot apply harmonised labels until review is complete:\n" + msg
        )

    final_product = {
        gid: review[gid]["final_product"].strip()
        for gid in expected_gids
    }

    # Coordinate-based mapping is robust to different display labels.
    coord_to_group = defaultdict(dict)

    for uid, g in genes.items():
        gid = uid_to_group[uid]
        key = (g["start"], g["end"], g["strand"])
        if key in coord_to_group[g["cluster"]]:
            raise ValueError(
                f"Duplicate CDS coordinate key in {g['cluster']}: {key}"
            )
        coord_to_group[g["cluster"]][key] = gid

    output_gbk_dir.mkdir(parents=True, exist_ok=True)

    change_rows = []

    for label in manifest_labels:
        src = input_gbk_dir / f"{label}.gbk"
        dst = output_gbk_dir / f"{label}.gbk"

        if not src.exists():
            raise FileNotFoundError(src)

        records = list(SeqIO.parse(src, "genbank"))
        if len(records) != 1:
            raise ValueError(
                f"{src}: expected one GenBank record, found {len(records)}"
            )

        rec = records[0]
        matched = 0
        cds_count = 0

        for feature in rec.features:
            if feature.type != "CDS":
                continue

            cds_count += 1
            key = (
                int(feature.location.start),
                int(feature.location.end),
                int(feature.location.strand or 0),
            )

            gid = coord_to_group[label].get(key)
            if gid is None:
                continue

            matched += 1
            new_product = final_product[gid]
            old_products = feature.qualifiers.get("product", [])
            old_product = old_products[0] if old_products else ""

            if old_product != new_product:
                add_note(
                    feature,
                    f"Original Pharokka/Phold product: {old_product or 'not specified'}",
                )
                feature.qualifiers["product"] = [new_product]

            add_note(feature, f"Clinker homology group: {gid}")

            change_rows.append(
                {
                    "genome": label,
                    "group_id": gid,
                    "start_0based": key[0],
                    "end_0based_exclusive": key[1],
                    "old_product": old_product,
                    "final_product": new_product,
                    "changed": int(old_product != new_product),
                }
            )

        if matched != len(coord_to_group[label]):
            raise RuntimeError(
                f"{label}: matched {matched} session CDS to GenBank features, "
                f"but expected {len(coord_to_group[label])}. "
                "Do not use partially harmonised output."
            )

        SeqIO.write(rec, dst, "genbank")
        print(
            f"{label}: wrote {dst} "
            f"({matched}/{cds_count} CDS mapped to Clinker groups)"
        )

    change_fields = [
        "genome",
        "group_id",
        "start_0based",
        "end_0based_exclusive",
        "old_product",
        "final_product",
        "changed",
    ]
    write_tsv(
        Path(args.root).resolve() / args.outdir / "harmonisation_changes.tsv",
        change_fields,
        change_rows,
    )

    changed = sum(int(r["changed"]) for r in change_rows)
    print()
    print(f"Product qualifiers changed: {changed}")
    print(
        "Original products were retained in CDS /note qualifiers where changed."
    )


def main():
    args = parse_args()
    try:
        if args.mode == "review":
            review_mode(args)
        else:
            apply_mode(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
