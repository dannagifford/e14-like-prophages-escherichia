#!/usr/bin/env python3
import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

from clinker.align import Globaligner


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def membership_key(genome, start, end, strand):
    return (
        str(genome),
        int(start),
        int(end),
        int(strand or 0),
    )


groups_path = Path(os.environ["GROUPS_TSV"])
review_path = Path(os.environ["REVIEW"])
membership_path = groups_path.parent / "gene_membership.tsv"

groups = read_tsv(groups_path)
review = {
    row["group_id"]: row
    for row in read_tsv(review_path)
}

# Original reviewed group membership.
members_by_gid = defaultdict(set)

for row in read_tsv(membership_path):
    members_by_gid[row["group_id"]].add(
        membership_key(
            row["genome"],
            row["start_0based"],
            row["end_0based_exclusive"],
            row["strand"],
        )
    )

records = []

for row in groups:
    gid = row["group_id"]

    # Singletons remain outside Clinker's coloured homology groups.
    if row["clinker_group_label"] == "singleton":
        continue

    product = (review[gid].get("final_product") or "").strip()
    if not product:
        raise SystemExit(f"ERROR: blank reviewed product for {gid}")

    members = frozenset(members_by_gid[gid])
    if not members:
        raise SystemExit(f"ERROR: no membership recorded for {gid}")

    records.append((gid, product, members))

product_counts = Counter(product for _, product, _ in records)
labels_by_members = {}

for gid, product, members in records:
    if members in labels_by_members:
        raise SystemExit("ERROR: duplicate reviewed group membership")

    labels_by_members[members] = (
        f"{product} [{gid}]"
        if product_counts[product] > 1
        else product
    )

sessions = [
    (
        Path(os.environ["MANIFEST_SESSION"]),
        Path(os.environ["MANIFEST_PRODUCT_SESSION"]),
    ),
    (
        Path(os.environ["AUTO_SESSION"]),
        Path(os.environ["AUTO_PRODUCT_SESSION"]),
    ),
]

# Validate both sessions before writing either output.
outputs = []

for src, dst in sessions:
    with src.open() as fh:
        ga = Globaligner.from_json(fh)

    with src.open() as fh:
        data = json.load(fh)

    gene_keys = {}

    for cluster in ga.clusters.values():
        for locus in cluster.loci:
            for gene in locus.genes:
                gene_keys[gene.uid] = membership_key(
                    cluster.name,
                    gene.start,
                    gene.end,
                    gene.strand,
                )

    labels_by_uid = {}
    seen_members = set()

    for group in ga.groups:
        members = frozenset(gene_keys[uid] for uid in group.genes)

        if members not in labels_by_members:
            raise SystemExit(
                f"ERROR: group membership changed in {src.name}. "
                "Regenerate and review stage 05 using the current inputs "
                "and the same Clinker identity threshold."
            )

        if members in seen_members:
            raise SystemExit(f"ERROR: duplicate group membership in {src}")

        seen_members.add(members)
        labels_by_uid[group.uid] = labels_by_members[members]

    if seen_members != set(labels_by_members):
        raise SystemExit(
            f"ERROR: reviewed groups are missing from {src.name}"
        )

    raw_groups = data.get("groups")
    if not isinstance(raw_groups, list):
        raise SystemExit(f"ERROR: missing JSON groups list in {src}")

    raw_uids = [group.get("uid") for group in raw_groups]
    if (
        len(raw_uids) != len(set(raw_uids))
        or set(raw_uids) != set(labels_by_uid)
    ):
        raise SystemExit(f"ERROR: inconsistent JSON group UIDs in {src}")

    for group in raw_groups:
        group["label"] = labels_by_uid[group["uid"]]

    outputs.append((dst, data))

for dst, data in outputs:
    with dst.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"Wrote: {dst}")
