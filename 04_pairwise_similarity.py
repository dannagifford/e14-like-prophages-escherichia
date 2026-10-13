#!/usr/bin/env python3
"""
04_pairwise_similarity.py

All-vs-all nucleotide comparison of the final prophage/phage FASTAs produced by
02_find_boundaries.py.

For each unordered pair A/B, BLAST is run in both directions:
    A -> B
    B -> A

For each direction, HSPs are greedily filtered by bitscore so that query
intervals do not overlap. This avoids double-counting overlapping HSPs when
calculating:

    * unique query coverage
    * aligned query bp
    * length-weighted nucleotide identity

Outputs:
    comparisons/blast/pairwise_similarity.tsv
    comparisons/blast/mean_reciprocal_coverage_matrix.tsv
    comparisons/blast/mean_weighted_identity_matrix.tsv

These summaries are descriptive. They are not intended to infer phylogeny.
"""

from __future__ import annotations

import argparse
import itertools
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="manifest.tsv")
    p.add_argument("--root", default=".")
    p.add_argument("--outdir", default="comparisons/blast")
    p.add_argument("--task", default="blastn")
    p.add_argument("--evalue", default="1e-10")
    return p.parse_args()


def load_one_fasta(path: Path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one record, found {len(records)}")
    return records[0]


def run_blast(query: Path, subject: Path, task: str, evalue: str):
    with tempfile.TemporaryDirectory(prefix="pairwise_blast_") as td:
        out = Path(td) / "hits.tsv"
        cmd = [
            "blastn",
            "-task", task,
            "-query", str(query),
            "-subject", str(subject),
            "-dust", "no",
            "-soft_masking", "false",
            "-evalue", evalue,
            "-outfmt",
            "6 qstart qend sstart send pident length nident bitscore evalue",
            "-out", str(out),
        ]
        subprocess.run(cmd, check=True)

        hits = []
        if out.exists():
            with out.open() as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    f = line.rstrip("\n").split("\t")
                    hits.append({
                        "qstart": int(f[0]),
                        "qend": int(f[1]),
                        "sstart": int(f[2]),
                        "send": int(f[3]),
                        "pident": float(f[4]),
                        "length": int(f[5]),
                        "nident": int(f[6]),
                        "bitscore": float(f[7]),
                        "evalue": float(f[8]),
                    })
        return hits


def overlap(a0, a1, b0, b1):
    a0, a1 = sorted((a0, a1))
    b0, b1 = sorted((b0, b1))
    return max(a0, b0) <= min(a1, b1)


def select_nonoverlapping_query_hsps(hits):
    """
    Greedily retain high-scoring HSPs without overlap on the query.
    This is deliberately conservative for query-coverage calculation.
    """
    chosen = []
    for h in sorted(hits, key=lambda x: (x["bitscore"], x["length"]), reverse=True):
        q0, q1 = sorted((h["qstart"], h["qend"]))
        if any(overlap(q0, q1, c["qstart"], c["qend"]) for c in chosen):
            continue
        chosen.append({
            **h,
            "qstart": q0,
            "qend": q1,
        })
    return sorted(chosen, key=lambda x: x["qstart"])


def summarize_direction(query_path: Path, subject_path: Path, task: str, evalue: str):
    qrec = load_one_fasta(query_path)
    qlen = len(qrec.seq)

    hits = run_blast(query_path, subject_path, task, evalue)
    chosen = select_nonoverlapping_query_hsps(hits)

    aligned_bp = sum(h["qend"] - h["qstart"] + 1 for h in chosen)
    coverage = 100.0 * aligned_bp / qlen if qlen else 0.0

    total_alignment_len = sum(h["length"] for h in chosen)
    weighted_identity = (
        100.0 * sum(h["nident"] for h in chosen) / total_alignment_len
        if total_alignment_len
        else 0.0
    )

    return {
        "query_length": qlen,
        "aligned_bp": aligned_bp,
        "coverage": coverage,
        "weighted_identity": weighted_identity,
        "hsp_count": len(chosen),
    }


def main():
    args = parse_args()

    root = Path(args.root).resolve()
    manifest = root / args.manifest
    outdir = root / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(manifest, sep="\t", dtype=str).fillna("")
    labels = df["label"].tolist()

    fasta_paths = {}
    lengths = {}

    for label in labels:
        path = root / "prophages" / label / f"{label}.fna"
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Run/review stage 02 before stage 04."
            )
        fasta_paths[label] = path
        lengths[label] = len(load_one_fasta(path).seq)

    rows = []

    for a, b in itertools.combinations(labels, 2):
        print(f"{a} <-> {b}")

        ab = summarize_direction(
            fasta_paths[a], fasta_paths[b], args.task, args.evalue
        )
        ba = summarize_direction(
            fasta_paths[b], fasta_paths[a], args.task, args.evalue
        )

        mean_cov = (ab["coverage"] + ba["coverage"]) / 2.0

        identities = [
            x["weighted_identity"]
            for x in (ab, ba)
            if x["aligned_bp"] > 0
        ]
        mean_identity = sum(identities) / len(identities) if identities else 0.0

        rows.append({
            "genome_A": a,
            "genome_B": b,
            "length_A": lengths[a],
            "length_B": lengths[b],
            "A_to_B_aligned_bp": ab["aligned_bp"],
            "A_to_B_query_coverage_pct": ab["coverage"],
            "A_to_B_weighted_identity_pct": ab["weighted_identity"],
            "A_to_B_nonoverlap_hsps": ab["hsp_count"],
            "B_to_A_aligned_bp": ba["aligned_bp"],
            "B_to_A_query_coverage_pct": ba["coverage"],
            "B_to_A_weighted_identity_pct": ba["weighted_identity"],
            "B_to_A_nonoverlap_hsps": ba["hsp_count"],
            "mean_reciprocal_query_coverage_pct": mean_cov,
            "mean_weighted_identity_pct": mean_identity,
        })

    result = pd.DataFrame(rows)
    detail_path = outdir / "pairwise_similarity.tsv"
    result.to_csv(detail_path, sep="\t", index=False, float_format="%.3f")

    cov_matrix = pd.DataFrame(100.0, index=labels, columns=labels)
    ident_matrix = pd.DataFrame(100.0, index=labels, columns=labels)

    for row in rows:
        a = row["genome_A"]
        b = row["genome_B"]
        cov = row["mean_reciprocal_query_coverage_pct"]
        ident = row["mean_weighted_identity_pct"]

        cov_matrix.loc[a, b] = cov
        cov_matrix.loc[b, a] = cov

        ident_matrix.loc[a, b] = ident
        ident_matrix.loc[b, a] = ident

    cov_path = outdir / "mean_reciprocal_coverage_matrix.tsv"
    ident_path = outdir / "mean_weighted_identity_matrix.tsv"

    cov_matrix.to_csv(cov_path, sep="\t", float_format="%.3f")
    ident_matrix.to_csv(ident_path, sep="\t", float_format="%.3f")

    print()
    print(f"Wrote {detail_path}")
    print(f"Wrote {cov_path}")
    print(f"Wrote {ident_path}")


if __name__ == "__main__":
    main()
