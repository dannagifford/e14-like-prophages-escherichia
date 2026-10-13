#!/usr/bin/env python3
"""
Helper utilities for 06_phylogeny.sh.

Search INPHARED with Mash, curate reference genomes, and prepare one
concatenated protein phylogeny requiring e14 in every selected group.
Mash pairwise distances are used only to remove near-duplicate references.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ACC_RE = re.compile(r"(?<![A-Z0-9])([A-Z]{1,5}_?\d{5,10}(?:\.\d+)?)(?![A-Z0-9])")


def accession_base(x: str) -> str:
    return str(x).strip().split(".")[0]


def read_manifest(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fieldnames, delimiter="\t",
            extrasaction="ignore", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)


def normalise_fasta(args):
    records = list(SeqIO.parse(args.input, "fasta"))
    if len(records) != 1:
        raise ValueError(f"{args.input}: expected one FASTA record, found {len(records)}")
    rec = records[0]
    rec.id = args.id
    rec.name = args.id
    rec.description = ""
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(rec, args.output, "fasta")


def load_inphared_metadata(path: Path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "Accession" not in df.columns:
        raise ValueError(f"{path}: expected an 'Accession' column")
    by_base = {}
    for _, r in df.iterrows():
        acc = r["Accession"].strip()
        if acc:
            by_base[accession_base(acc)] = r.to_dict()
    return df, by_base


def accession_from_ref(ref: str, metadata_by_base: dict) -> str | None:
    # First try accession-like strings embedded in the Mash reference name.
    for m in ACC_RE.findall(ref):
        b = accession_base(m)
        if b in metadata_by_base:
            return b

    # Common fallback: reference itself or its basename is the accession.
    candidates = {
        accession_base(ref),
        accession_base(Path(ref).name),
        accession_base(Path(ref).stem),
        accession_base(ref.split()[0]),
    }
    for c in candidates:
        if c in metadata_by_base:
            return c
    return None


def parse_hits(args):
    manifest = read_manifest(Path(args.manifest))

    # All manifest genomes are predefined downstream taxa, so none should be
    # downloaded again from INPHARED.
    focal_accessions = {
        accession_base(r["accession"])
        for r in manifest
        if r.get("accession")
    }

    discovery_types = {
        x.strip()
        for x in args.discovery_types.split(",")
        if x.strip()
    }
    if not discovery_types:
        raise ValueError("At least one discovery type is required.")

    # Only selected biological seed types spawn INPHARED neighbourhood
    # searches. Free-phage comparators remain in downstream analyses but do
    # not recruit a second ring of reference genomes.
    labels = [
        r["label"]
        for r in manifest
        if r.get("label") and r.get("type") in discovery_types
    ]

    if not labels:
        raise ValueError(
            "No manifest entries matched discovery types: "
            + ",".join(sorted(discovery_types))
        )

    _, meta = load_inphared_metadata(Path(args.metadata))

    all_rows = []
    for label in labels:
        p = Path(args.raw_dir) / f"{label}.mash.tsv"
        if not p.exists():
            raise FileNotFoundError(p)

        rows = []
        with p.open() as fh:
            for line in fh:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                ref, query, dist, pval, hashes = f[:5]
                try:
                    distf = float(dist)
                except ValueError:
                    continue
                if distf > args.max_distance:
                    continue

                acc = accession_from_ref(ref, meta)
                if not acc:
                    continue
                m = meta[acc]
                rows.append({
                    "query_label": label,
                    "reference_name": ref,
                    "accession": acc,
                    "mash_distance": distf,
                    "mash_pval": pval,
                    "matching_hashes": hashes,
                    "description": m.get("Description", ""),
                    "host": m.get("Host", ""),
                    "genome_length_bp": m.get("Genome_Length_(bp)", ""),
                })

        rows.sort(key=lambda x: (x["mash_distance"], x["accession"]))
        # Keep every hit that passes the Mash-distance threshold.
        # De-replication and the global cap are applied later, after
        # all-vs-all Mash distances have been calculated.
        all_rows.extend(rows)

    if not all_rows:
        raise ValueError("No INPHARED hits could be parsed from the Mash outputs.")

    hit_fields = [
        "query_label", "reference_name", "accession", "mash_distance",
        "mash_pval", "matching_hashes", "description", "host",
        "genome_length_bp",
    ]
    write_tsv(Path(args.out), all_rows, hit_fields)

    # Rank unique references by their closest focal-genome Mash distance.
    by_acc = {}
    for r in all_rows:
        acc = r["accession"]
        if acc in focal_accessions:
            continue
        cur = by_acc.get(acc)
        if cur is None or r["mash_distance"] < cur["min_mash_distance"]:
            by_acc[acc] = {
                "accession": acc,
                "description": r["description"],
                "host": r["host"],
                "genome_length_bp": r["genome_length_bp"],
                "min_mash_distance": r["mash_distance"],
                "best_query_label": r["query_label"],
            }

    refs = sorted(
        by_acc.values(),
        key=lambda x: (x["min_mash_distance"], x["accession"])
    )

    print("Mash discovery seeds: " + ", ".join(labels))
    print(f"Mash hits within threshold: {len(all_rows)}")
    print(f"Unique non-focal INPHARED accessions to download: {len(refs)}")

    write_tsv(
        Path(args.download_list),
        refs,
        [
            "accession", "description", "host", "genome_length_bp",
            "min_mash_distance", "best_query_label",
        ],
    )


def make_taxa(args):
    manifest = read_manifest(Path(args.manifest))
    rows = []

    for r in manifest:
        label = r["label"]
        fna = Path(args.genome_dir) / f"{label}.fna"
        if not fna.exists():
            raise FileNotFoundError(fna)
        rows.append({
            "taxon": label,
            "source": "focal",
            "is_focal": 1,
            "accession": accession_base(r.get("accession", "")),
            "description": r.get("organism", ""),
            "min_query_distance": 0.0,
            "path": str(fna),
        })

    with Path(args.download_list).open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            taxon = f"INPHARED_{accession_base(r['accession'])}"
            fna = Path(args.genome_dir) / f"{taxon}.fna"
            if not fna.exists():
                continue
            rows.append({
                "taxon": taxon,
                "source": "INPHARED",
                "is_focal": 0,
                "accession": accession_base(r["accession"]),
                "description": r.get("description", ""),
                "min_query_distance": float(r["min_mash_distance"]),
                "path": str(fna),
            })

    write_tsv(
        Path(args.out), rows,
        ["taxon", "source", "is_focal", "accession", "description",
         "min_query_distance", "path"],
    )


def read_pairwise(path: Path):
    d = {}
    names = set()
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            a, b, dist = f[0], f[1], f[2]
            try:
                val = float(dist)
            except ValueError:
                continue
            names.update([a, b])
            d[(a, b)] = val
    return names, d


def pair_distance(a, b, d):
    if a == b:
        return 0.0
    vals = []
    if (a, b) in d:
        vals.append(d[(a, b)])
    if (b, a) in d:
        vals.append(d[(b, a)])
    if not vals:
        return math.nan
    return sum(vals) / len(vals)


def curate_genomes(args):
    taxa = pd.read_csv(args.taxa, sep="\t", dtype={"taxon": str})
    _, d = read_pairwise(Path(args.pairwise))

    focal = taxa[taxa["is_focal"].astype(int) == 1].copy()
    refs = taxa[taxa["is_focal"].astype(int) == 0].copy()
    refs = refs.sort_values(["min_query_distance", "accession"])

    retained = list(focal["taxon"])
    retained_refs = []
    dropped = []

    max_refs = max(0, args.max_total - len(retained))

    for _, row in refs.iterrows():
        taxon = row["taxon"]
        if len(retained_refs) >= max_refs:
            dropped.append({**row.to_dict(), "drop_reason": "max_total"})
            continue

        near = None
        near_dist = None
        for kept in retained:
            x = pair_distance(taxon, kept, d)
            if math.isnan(x):
                continue
            if x < args.redundancy_distance:
                near = kept
                near_dist = x
                break

        if near is not None:
            dropped.append({
                **row.to_dict(),
                "drop_reason": f"near_duplicate_of:{near}",
                "nearest_retained_distance": near_dist,
            })
        else:
            retained.append(taxon)
            retained_refs.append(taxon)

    keep_df = taxa[taxa["taxon"].isin(retained)].copy()
    order = {x: i for i, x in enumerate(retained)}
    keep_df["_order"] = keep_df["taxon"].map(order)
    keep_df = keep_df.sort_values("_order").drop(columns="_order")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    keep_df.to_csv(args.out, sep="\t", index=False)

    drop_fields = list(taxa.columns) + ["drop_reason", "nearest_retained_distance"]
    write_tsv(Path(args.dropped), dropped, drop_fields)


def sanitise_proteins(args):
    out_records = []
    for i, rec in enumerate(SeqIO.parse(args.input, "fasta"), start=1):
        seq = str(rec.seq).strip().rstrip("*").replace("*", "X")
        if not seq:
            continue
        rid = f"{args.taxon}|CDS{i:05d}"
        out_records.append(SeqRecord(rec.seq.__class__(seq), id=rid, name=rid, description=""))
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_records, args.output, "fasta")


def prepare_mmseqs_fasta(args):
    """Write a strict, deterministic one-line protein FASTA for MMseqs."""
    taxa_df = pd.read_csv(args.taxa, sep="\t")
    taxa = list(taxa_df["taxon"])

    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    rows = []
    proteins = []
    seen = set()

    for taxon in taxa:
        faa = Path(args.protein_dir) / f"{taxon}.faa"
        if not faa.exists():
            raise FileNotFoundError(faa)

        recs = list(SeqIO.parse(faa, "fasta"))
        if not recs:
            raise ValueError(f"{faa}: no protein records")

        total_aa = 0
        for i, rec in enumerate(recs, start=1):
            rid = f"{taxon}|CDS{i:05d}"
            if rid in seen:
                raise ValueError(f"Duplicate protein ID: {rid}")
            seen.add(rid)

            raw = str(rec.seq).upper().replace("*", "X")
            seq = "".join(c if c in allowed else "X" for c in raw)
            if not seq:
                raise ValueError(f"{faa}: empty sequence for {rid}")

            total_aa += len(seq)
            proteins.append((rid, seq))

        rows.append({
            "taxon": taxon,
            "n_proteins": len(recs),
            "total_aa": total_aa,
        })

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.output).open("w") as fh:
        for rid, seq in proteins:
            fh.write(f">{rid}\n{seq}\n")

    write_tsv(
        Path(args.summary),
        rows,
        ["taxon", "n_proteins", "total_aa"],
    )

    print(f"Taxa: {len(taxa)}")
    print(f"Proteins: {len(proteins)}")
    print(f"Wrote: {args.output}")

def select_backbone(args):
    taxa_df = pd.read_csv(args.taxa, sep="\t")
    taxa = list(taxa_df["taxon"])
    taxa_set = set(taxa)

    anchor = args.anchor.strip()
    if not anchor or anchor not in taxa_set:
        raise ValueError(f"Anchor taxon '{anchor}' is not in retained taxa.")

    seqs = {rec.id: rec for rec in SeqIO.parse(args.proteins, "fasta")}
    seq_to_taxon = {sid: sid.split("|", 1)[0] for sid in seqs}

    clusters = defaultdict(list)
    with Path(args.cluster_tsv).open() as fh:
        for line in fh:
            if not line.strip():
                continue
            rep, member = line.rstrip("\n").split("\t")[:2]
            clusters[rep].append(member)

    selected = []
    for rep, members in clusters.items():
        # Some MMseqs outputs may omit the representative as an explicit member.
        members = list(dict.fromkeys(members + ([rep] if rep in seqs else [])))
        members = [m for m in members if m in seqs and seq_to_taxon[m] in taxa_set]
        if not members:
            continue

        by_taxon = defaultdict(list)
        for m in members:
            by_taxon[seq_to_taxon[m]].append(m)

        if anchor not in by_taxon:
            continue
        if any(len(v) != 1 for v in by_taxon.values()):
            continue

        occupancy = len(by_taxon) / len(taxa)
        if occupancy < args.min_fraction:
            continue

        lengths = [len(seqs[m].seq) for m in members]
        if statistics.median(lengths) < args.min_aa:
            continue

        selected.append({
            "rep": rep,
            "members": members,
            "by_taxon": by_taxon,
            "occupancy": occupancy,
            "median_aa": statistics.median(lengths),
        })

    selected.sort(key=lambda x: (-x["occupancy"], -x["median_aa"], x["rep"]))

    outdir = Path(args.outdir)
    groupdir = outdir / "groups"
    groupdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for i, c in enumerate(selected, start=1):
        gid = f"BACKBONE_{i:03d}"
        recs = []
        for taxon in taxa:
            if taxon not in c["by_taxon"]:
                continue
            sid = c["by_taxon"][taxon][0]
            rec = seqs[sid]
            recs.append(SeqRecord(rec.seq, id=taxon, name=taxon, description=""))
        SeqIO.write(recs, groupdir / f"{gid}.faa", "fasta")
        rows.append({
            "group_id": gid,
            "mmseqs_representative": c["rep"],
            "n_taxa": len(c["by_taxon"]),
            "total_taxa": len(taxa),
            "occupancy": f"{c['occupancy']:.4f}",
            "median_length_aa": c["median_aa"],
            "anchor": anchor,
        })

    write_tsv(
        outdir / "backbone_groups.tsv",
        rows,
        ["group_id", "mmseqs_representative", "n_taxa", "total_taxa",
         "occupancy", "median_length_aa", "anchor"],
    )

    if len(rows) < args.min_groups:
        raise ValueError(
            f"Only {len(rows)} single-copy groups anchored on {anchor} passed the "
            f"occupancy threshold ({args.min_fraction}); "
            f"{args.min_groups} are required."
        )


def concatenate(args):
    taxa_df = pd.read_csv(args.taxa, sep="\t")
    taxa = list(taxa_df["taxon"])
    groups = pd.read_csv(args.groups, sep="\t")
    anchor = args.anchor.strip()
    if not anchor or anchor not in taxa:
        raise ValueError(f"Anchor taxon '{anchor}' is not in retained taxa.")
    if groups.empty or not (groups["anchor"] == anchor).all():
        raise ValueError("All selected groups must require the requested e14 anchor.")

    concatenated = {t: [] for t in taxa}
    partitions = []
    start = 1

    for gid in groups["group_id"]:
        p = Path(args.align_dir) / f"{gid}.clipkit.faa"
        if not p.exists():
            raise FileNotFoundError(p)

        recs = list(SeqIO.parse(p, "fasta"))
        if not recs:
            raise ValueError(f"{p}: no alignment records after trimming")
        lengths = {len(r.seq) for r in recs}
        if len(lengths) != 1:
            raise ValueError(f"{p}: alignment sequences have unequal lengths")
        L = lengths.pop()
        if L == 0:
            raise ValueError(f"{p}: zero alignment columns after trimming")
        by_taxon = {r.id: str(r.seq) for r in recs}
        if len(by_taxon) != len(recs):
            raise ValueError(f"{p}: duplicate taxon IDs")
        if set(by_taxon) - set(taxa):
            raise ValueError(f"{p}: unknown taxa in alignment")
        if anchor not in by_taxon:
            raise ValueError(f"{p}: missing anchor taxon '{anchor}'")

        for t in taxa:
            concatenated[t].append(by_taxon.get(t, "-" * L))

        end = start + L - 1
        partitions.append((gid, start, end))
        start = end + 1

    out_records = []
    coverage_rows = []
    from Bio.Seq import Seq

    for t in taxa:
        chunks = concatenated[t]
        seq = "".join(chunks)
        nongap = sum(1 for x in seq.upper() if x not in "-?X")
        coverage = nongap / len(seq) if seq else 0
        n_groups_present = sum(
            1 for chunk in chunks if any(x not in "-?X" for x in chunk.upper())
        )
        fraction_groups_present = (
            n_groups_present / len(chunks) if chunks else 0
        )
        retain = n_groups_present >= args.min_taxon_groups

        coverage_rows.append({
            "taxon": t,
            "aligned_sites": len(seq),
            "non_gap_sites": nongap,
            "fraction_non_gap": f"{coverage:.4f}",
            "groups_present": n_groups_present,
            "total_groups": len(chunks),
            "fraction_groups_present": f"{fraction_groups_present:.4f}",
            "retained_in_tree": int(retain),
        })

        if retain:
            out_records.append(
                SeqRecord(Seq(seq), id=t, name=t, description="")
            )

    # Keep diagnostics even if the tree cannot be built.
    write_tsv(
        Path(args.coverage), coverage_rows,
        ["taxon", "aligned_sites", "non_gap_sites", "fraction_non_gap",
         "groups_present", "total_groups", "fraction_groups_present",
         "retained_in_tree"],
    )
    if anchor not in {r.id for r in out_records}:
        raise ValueError(
            f"Anchor '{anchor}' has fewer than {args.min_taxon_groups} groups "
            "with residues after trimming; refusing to build a tree without e14."
        )
    if len(out_records) < 4:
        raise ValueError(
            f"Only {len(out_records)} taxa have at least "
            f"{args.min_taxon_groups} backbone groups; at least 4 are required."
        )

    represented = [
        any(any(x not in "-?X" for x in concatenated[r.id][i].upper())
            for r in out_records)
        for i in range(len(partitions))
    ]
    if not all(represented):
        raise ValueError("A selected group has no residues in the retained tree taxa.")
    if len(partitions) < args.min_groups:
        raise ValueError(
            f"Only {len(partitions)} trimmed groups; {args.min_groups} are required."
        )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_records, args.output, "fasta")
    with Path(args.partitions).open("w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for gid, a, b in partitions:
            fh.write(f"    charset {gid} = {a}-{b};\n")
        fh.write("end;\n")


def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers(dest="cmd", required=True)

    q = sp.add_parser("normalise-fasta")
    q.add_argument("--input", required=True)
    q.add_argument("--output", required=True)
    q.add_argument("--id", required=True)
    q.set_defaults(func=normalise_fasta)

    q = sp.add_parser("parse-hits")
    q.add_argument("--manifest", required=True)
    q.add_argument("--metadata", required=True)
    q.add_argument("--raw-dir", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--download-list", required=True)
    q.add_argument("--max-distance", type=float, default=0.15)
    q.add_argument("--discovery-types", default="prophage")
    q.set_defaults(func=parse_hits)

    q = sp.add_parser("make-taxa")
    q.add_argument("--manifest", required=True)
    q.add_argument("--download-list", required=True)
    q.add_argument("--genome-dir", required=True)
    q.add_argument("--out", required=True)
    q.set_defaults(func=make_taxa)

    q = sp.add_parser("curate-genomes")
    q.add_argument("--pairwise", required=True)
    q.add_argument("--taxa", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--dropped", required=True)
    q.add_argument("--redundancy-distance", type=float, default=0.002)
    q.add_argument("--max-total", type=int, default=150)
    q.set_defaults(func=curate_genomes)

    q = sp.add_parser("sanitise-proteins")
    q.add_argument("--input", required=True)
    q.add_argument("--output", required=True)
    q.add_argument("--taxon", required=True)
    q.set_defaults(func=sanitise_proteins)

    q = sp.add_parser("prepare-mmseqs-fasta")
    q.add_argument("--taxa", required=True)
    q.add_argument("--protein-dir", required=True)
    q.add_argument("--output", required=True)
    q.add_argument("--summary", required=True)
    q.set_defaults(func=prepare_mmseqs_fasta)

    q = sp.add_parser("select-backbone")
    q.add_argument("--cluster-tsv", required=True)
    q.add_argument("--proteins", required=True)
    q.add_argument("--taxa", required=True)
    q.add_argument("--anchor", required=True)
    q.add_argument("--min-fraction", type=float, default=0.25)
    q.add_argument("--min-aa", type=int, default=50)
    q.add_argument("--min-groups", type=int, default=4)
    q.add_argument("--outdir", required=True)
    q.set_defaults(func=select_backbone)

    q = sp.add_parser("concatenate")
    q.add_argument("--align-dir", required=True)
    q.add_argument("--taxa", required=True)
    q.add_argument("--groups", required=True)
    q.add_argument("--output", required=True)
    q.add_argument("--partitions", required=True)
    q.add_argument("--coverage", required=True)
    q.add_argument("--min-taxon-groups", type=int, default=3)
    q.add_argument("--anchor", required=True)
    q.add_argument("--min-groups", type=int, default=4)
    q.set_defaults(func=concatenate)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
