#!/usr/bin/env python3
"""
02_find_boundaries.py

Infer conservative attL/attR-associated boundaries for the e14-related
prophage candidates listed in manifest.tsv.

The key biological assumption is that e14 integrates at the unusual icd
locus. Rather than searching first for arbitrary short direct repeats, this
script:

1. Uses canonical K-12 e14 as a positive control.
2. Derives the longer icd-associated homology shared across the two e14
   boundaries.
3. Searches each candidate host-context region for two homologous copies of
   that e14/icd boundary sequence bracketing the seed prophage locus.
4. Projects the K-12 e14 boundary positions through the alignments to obtain
   attL/attR point estimates where possible.
5. Uses the outer edges of the homologous blocks as conservative extraction
   boundaries, so the attachment-associated sequence is retained.
6. Falls back to direct left-vs-right boundary homology, then to the seed
   interval, if the icd-reference search is inconclusive.

This stage intentionally produces *best guesses*, not claims of exact
crossover nucleotides. The exact ~11-bp crossover core should be confirmed
manually before publication.

Inputs expected from 01_download.sh:
    genomes/<label>/<label>.fna
    genomes/<label>/<label>.gbk

Outputs:
    context/<label>/<label>_context.fna
    context/<label>/<label>_context.gbk
    boundaries/e14_icd_att_reference.fna
    boundaries/e14_reference_hsp.tsv
    boundaries/<label>_boundary_candidates.tsv
    boundaries/boundaries.tsv
    prophages/<label>/<label>.fna

Usage:
    conda run --no-capture-output -n prophage_compare \
        python 02_find_boundaries.py --manifest manifest.tsv
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class Hit:
    qstart: int
    qend: int
    sstart: int
    send: int
    pident: float
    length: int
    bitscore: float
    evalue: float
    qseq: str
    sseq: str

    @property
    def subject_strand(self) -> str:
        return "+" if self.send >= self.sstart else "-"

    @property
    def subject_min(self) -> int:
        return min(self.sstart, self.send)

    @property
    def subject_max(self) -> int:
        return max(self.sstart, self.send)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="manifest.tsv")
    p.add_argument("--root", default=".")
    p.add_argument(
        "--reference-window",
        type=int,
        default=2000,
        help="bp either side of known K-12 e14 boundaries used to derive icd homology",
    )
    p.add_argument(
        "--min-icd-hit-length",
        type=int,
        default=70,
        help="minimum aligned length for e14/icd-reference hits",
    )
    p.add_argument(
        "--min-icd-identity",
        type=float,
        default=75.0,
        help="minimum nucleotide identity for e14/icd-reference hits",
    )
    p.add_argument(
        "--fallback-max-repeat",
        type=int,
        default=1000,
        help="maximum HSP length considered in generic boundary-repeat fallback",
    )
    return p.parse_args()


def as_int(value) -> Optional[int]:
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return None
    return int(float(s))


def load_single_fasta(path: Path) -> SeqRecord:
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one FASTA record, found {len(records)}")
    return records[0]


def load_single_genbank(path: Path) -> SeqRecord:
    records = list(SeqIO.parse(path, "genbank"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one GenBank record, found {len(records)}")
    return records[0]


def write_fasta(seq: Seq, path: Path, record_id: str, description: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(SeqRecord(seq, id=record_id, description=description), path, "fasta")


def slice_record(record: SeqRecord, start_1: int, end_1: int, new_id: str) -> SeqRecord:
    """Slice a SeqRecord using 1-based inclusive genomic coordinates."""
    sub = record[start_1 - 1 : end_1]
    sub.id = new_id
    sub.name = new_id[:16]
    sub.description = f"{record.id}:{start_1}-{end_1}"
    sub.annotations["molecule_type"] = "DNA"
    return sub


def run_blast(query: Seq, subject: Seq, tmpdir: Path) -> list[Hit]:
    qfile = tmpdir / "query.fna"
    sfile = tmpdir / "subject.fna"
    outfile = tmpdir / "blast.tsv"

    write_fasta(query, qfile, "query")
    write_fasta(subject, sfile, "subject")

    outfmt = "6 qstart qend sstart send pident length bitscore evalue qseq sseq"
    cmd = [
        "blastn",
        "-task",
        "blastn",
        "-query",
        str(qfile),
        "-subject",
        str(sfile),
        "-dust",
        "no",
        "-soft_masking",
        "false",
        "-evalue",
        "1e-5",
        "-outfmt",
        outfmt,
        "-out",
        str(outfile),
    ]
    subprocess.run(cmd, check=True)

    hits: list[Hit] = []
    if not outfile.exists():
        return hits

    with outfile.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            hits.append(
                Hit(
                    qstart=int(f[0]),
                    qend=int(f[1]),
                    sstart=int(f[2]),
                    send=int(f[3]),
                    pident=float(f[4]),
                    length=int(f[5]),
                    bitscore=float(f[6]),
                    evalue=float(f[7]),
                    qseq=f[8],
                    sseq=f[9],
                )
            )
    return hits


def distance_to_interval(x: int, a: int, b: int) -> int:
    lo, hi = min(a, b), max(a, b)
    if lo <= x <= hi:
        return 0
    return min(abs(x - lo), abs(x - hi))


def map_query_position_to_subject(hit: Hit, query_pos: int) -> Optional[int]:
    """
    Map a 1-based position in the *full query sequence used for BLAST*
    to the corresponding 1-based subject position in this HSP.
    """
    qstep = 1 if hit.qend >= hit.qstart else -1
    sstep = 1 if hit.send >= hit.sstart else -1

    qpos = hit.qstart
    spos = hit.sstart

    for qchar, schar in zip(hit.qseq, hit.sseq):
        current_q = qpos if qchar != "-" else None
        current_s = spos if schar != "-" else None

        if current_q == query_pos:
            return current_s

        if qchar != "-":
            qpos += qstep
        if schar != "-":
            spos += sstep

    return None


def map_subject_position_to_query(hit: Hit, subject_pos: int) -> Optional[int]:
    """
    Map a 1-based position in the *full subject sequence used for BLAST*
    to the corresponding 1-based query position in this HSP.
    """
    qstep = 1 if hit.qend >= hit.qstart else -1
    sstep = 1 if hit.send >= hit.sstart else -1

    qpos = hit.qstart
    spos = hit.sstart

    for qchar, schar in zip(hit.qseq, hit.sseq):
        current_q = qpos if qchar != "-" else None
        current_s = spos if schar != "-" else None

        if current_s == subject_pos:
            return current_q

        if qchar != "-":
            qpos += qstep
        if schar != "-":
            spos += sstep

    return None


def find_icd_features(record: SeqRecord, context_start: int, context_end: int) -> list[str]:
    """Return simple descriptions of annotated icd/isocitrate-dehydrogenase features."""
    found: list[str] = []
    for feat in record.features:
        if feat.type not in {"CDS", "gene"}:
            continue

        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        if end < context_start or start > context_end:
            continue

        gene = ";".join(feat.qualifiers.get("gene", []))
        product = ";".join(feat.qualifiers.get("product", []))
        label = f"{gene} {product}".lower()

        if gene.lower() == "icd" or "isocitrate dehydrogenase" in label:
            strand = feat.location.strand
            found.append(f"{start}-{end}({strand:+d})" if strand else f"{start}-{end}")

    return found


def derive_k12_reference(
    known_row: pd.Series,
    fasta_record: SeqRecord,
    boundaries_dir: Path,
    reference_window: int,
) -> tuple[Seq, int, int, dict]:
    """
    Derive the longer icd-associated homology from the two canonical e14
    boundaries and locate the known e14 start/end within that homology.

    Returns:
        reference_seq
        left_boundary_pos_in_reference (1-based)
        right_boundary_pos_in_reference (1-based)
        metadata
    """
    seed_start = as_int(known_row["seed_start"])
    seed_end = as_int(known_row["seed_end"])
    if seed_start is None or seed_end is None:
        raise ValueError("known_icd row requires seed_start and seed_end")

    genome_len = len(fasta_record.seq)
    l0 = max(1, seed_start - reference_window)
    l1 = min(genome_len, seed_start + reference_window)
    r0 = max(1, seed_end - reference_window)
    r1 = min(genome_len, seed_end + reference_window)

    left = fasta_record.seq[l0 - 1 : l1]
    right = fasta_record.seq[r0 - 1 : r1]

    with tempfile.TemporaryDirectory(prefix="e14_ref_") as td:
        hits = run_blast(left, right, Path(td))

    candidates = []
    for h in hits:
        if h.length < 60:
            continue
        if h.pident < 70:
            continue

        qg0 = l0 + min(h.qstart, h.qend) - 1
        qg1 = l0 + max(h.qstart, h.qend) - 1
        sg0 = r0 + min(h.sstart, h.send) - 1
        sg1 = r0 + max(h.sstart, h.send) - 1

        dleft = distance_to_interval(seed_start, qg0, qg1)
        dright = distance_to_interval(seed_end, sg0, sg1)

        # Strongly favour an HSP that actually spans the known boundaries.
        span_bonus = 5000 if dleft == 0 else 0
        span_bonus += 5000 if dright == 0 else 0

        # Prefer a plausible icd-scale homology rather than a huge unrelated repeat.
        length_penalty = abs(h.length - 216) * 0.5
        score = h.bitscore + span_bonus - (dleft + dright) - length_penalty

        candidates.append((score, h, qg0, qg1, sg0, sg1, dleft, dright))

    if not candidates:
        raise RuntimeError(
            "Could not derive an icd-associated boundary homology from canonical K-12 e14."
        )

    candidates.sort(key=lambda x: x[0], reverse=True)
    score, best, qg0, qg1, sg0, sg1, dleft, dright = candidates[0]

    # Reference is the raw left-side genomic sequence corresponding to the HSP.
    qlo, qhi = min(best.qstart, best.qend), max(best.qstart, best.qend)
    reference = left[qlo - 1 : qhi]

    left_ref_pos = seed_start - qg0 + 1
    if not (1 <= left_ref_pos <= len(reference)):
        # If the published/EcoCyc boundary lies just outside the best HSP,
        # use the closest edge and flag this in metadata.
        left_ref_pos = min(max(left_ref_pos, 1), len(reference))

    # Convert the known right boundary to the corresponding position in the
    # left HSP/query, then into reference coordinates.
    right_local_subject = seed_end - r0 + 1
    qpos_for_right = map_subject_position_to_query(best, right_local_subject)

    if qpos_for_right is None:
        # Conservative fallback: use whichever end of the reference is closer
        # in the K-12 HSP geometry.
        right_ref_pos = left_ref_pos
        right_projection = "fallback_same_as_left"
    else:
        right_ref_pos = qpos_for_right - qlo + 1
        right_ref_pos = min(max(right_ref_pos, 1), len(reference))
        right_projection = "mapped"

    write_fasta(
        reference,
        boundaries_dir / "e14_icd_att_reference.fna",
        "e14_icd_att_reference",
        f"derived from K-12 e14 boundary homology; {qg0}-{qg1}",
    )

    with (boundaries_dir / "e14_reference_hsp.tsv").open("w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "left_homology_start",
                "left_homology_end",
                "right_homology_start",
                "right_homology_end",
                "identity",
                "alignment_length",
                "bitscore",
                "known_attL",
                "known_attR",
                "attL_ref_pos",
                "attR_ref_pos",
                "attR_projection",
            ]
        )
        writer.writerow(
            [
                qg0,
                qg1,
                sg0,
                sg1,
                f"{best.pident:.3f}",
                best.length,
                f"{best.bitscore:.1f}",
                seed_start,
                seed_end,
                left_ref_pos,
                right_ref_pos,
                right_projection,
            ]
        )

    meta = {
        "left_block_start": qg0,
        "left_block_end": qg1,
        "right_block_start": sg0,
        "right_block_end": sg1,
        "identity": best.pident,
        "length": best.length,
        "score": score,
        "left_distance": dleft,
        "right_distance": dright,
        "right_projection": right_projection,
    }
    return reference, left_ref_pos, right_ref_pos, meta


def hit_to_genome(hit: Hit, context_start: int) -> tuple[int, int]:
    return (
        context_start + hit.subject_min - 1,
        context_start + hit.subject_max - 1,
    )


def choose_icd_pair(
    hits: list[Hit],
    context_start: int,
    seed_start: int,
    seed_end: int,
    attL_ref_pos: int,
    attR_ref_pos: int,
) -> Optional[dict]:
    """
    Choose two e14/icd-reference hits that best bracket the seed locus.
    """
    usable = []
    for h in hits:
        gs, ge = hit_to_genome(h, context_start)
        centre = (gs + ge) / 2

        attL_local = map_query_position_to_subject(h, attL_ref_pos)
        attR_local = map_query_position_to_subject(h, attR_ref_pos)

        attL_proj = context_start + attL_local - 1 if attL_local is not None else None
        attR_proj = context_start + attR_local - 1 if attR_local is not None else None

        usable.append(
            {
                "hit": h,
                "start": gs,
                "end": ge,
                "centre": centre,
                "strand": h.subject_strand,
                "attL_proj": attL_proj,
                "attR_proj": attR_proj,
            }
        )

    if len(usable) < 2:
        return None

    midpoint = (seed_start + seed_end) / 2
    pairs = []

    for left in usable:
        for right in usable:
            if left is right:
                continue
            if left["centre"] >= right["centre"]:
                continue

            # Prefer bracketing the candidate rather than two nearby paralogues.
            bracket_penalty = 0
            if left["centre"] > midpoint:
                bracket_penalty += abs(left["centre"] - seed_start)
            if right["centre"] < midpoint:
                bracket_penalty += abs(seed_end - right["centre"])

            edge_distance = abs(left["centre"] - seed_start) + abs(right["centre"] - seed_end)

            strand_bonus = 500 if left["strand"] == right["strand"] else 0
            projection_bonus = 250 if left["attL_proj"] is not None else 0
            projection_bonus += 250 if right["attR_proj"] is not None else 0

            score = (
                left["hit"].bitscore
                + right["hit"].bitscore
                + strand_bonus
                + projection_bonus
                - 0.05 * edge_distance
                - 0.2 * bracket_penalty
            )

            pairs.append((score, left, right))

    if not pairs:
        return None

    pairs.sort(key=lambda x: x[0], reverse=True)
    score, left, right = pairs[0]

    return {"score": score, "left": left, "right": right}


def choose_fallback_pair(
    left_seq: Seq,
    right_seq: Seq,
    left_window_start: int,
    right_window_start: int,
    seed_start: int,
    seed_end: int,
    max_repeat: int,
) -> Optional[dict]:
    """
    Generic direct-repeat fallback: align the two boundary search windows
    directly and choose a reasonably short homologous block near both seed
    edges.
    """
    with tempfile.TemporaryDirectory(prefix="att_fallback_") as td:
        hits = run_blast(left_seq, right_seq, Path(td))

    candidates = []
    for h in hits:
        if h.length < 20 or h.length > max_repeat:
            continue
        if h.pident < 70:
            continue

        lg0 = left_window_start + min(h.qstart, h.qend) - 1
        lg1 = left_window_start + max(h.qstart, h.qend) - 1
        rg0 = right_window_start + min(h.sstart, h.send) - 1
        rg1 = right_window_start + max(h.sstart, h.send) - 1

        lcentre = (lg0 + lg1) / 2
        rcentre = (rg0 + rg1) / 2

        edge_distance = abs(lcentre - seed_start) + abs(rcentre - seed_end)
        length_penalty = abs(h.length - 216) * 0.2
        orientation_bonus = 250 if h.subject_strand == "+" else 0

        score = h.bitscore + orientation_bonus - 0.05 * edge_distance - length_penalty
        candidates.append(
            (
                score,
                h,
                {
                    "start": lg0,
                    "end": lg1,
                    "centre": lcentre,
                    "strand": "+",
                },
                {
                    "start": rg0,
                    "end": rg1,
                    "centre": rcentre,
                    "strand": h.subject_strand,
                },
            )
        )

    if not candidates:
        return None

    candidates.sort(key=lambda x: x[0], reverse=True)
    score, hit, left, right = candidates[0]
    return {"score": score, "hit": hit, "left": left, "right": right}


def confidence_for_icd_pair(pair: dict, seed_start: int, seed_end: int) -> str:
    left = pair["left"]
    right = pair["right"]

    min_ident = min(left["hit"].pident, right["hit"].pident)
    min_len = min(left["hit"].length, right["hit"].length)
    edge_dist = abs(left["centre"] - seed_start) + abs(right["centre"] - seed_end)

    if (
        min_ident >= 85
        and min_len >= 100
        and edge_dist <= 10000
        and left["strand"] == right["strand"]
    ):
        return "high"
    if min_ident >= 75 and min_len >= 70:
        return "medium"
    return "low"


def write_candidate_hits(
    path: Path,
    hits: list[Hit],
    context_start: int,
    attL_ref_pos: int,
    attR_ref_pos: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "subject_start",
                "subject_end",
                "strand",
                "identity",
                "alignment_length",
                "bitscore",
                "evalue",
                "query_start",
                "query_end",
                "projected_attL",
                "projected_attR",
            ]
        )

        for h in sorted(hits, key=lambda x: x.bitscore, reverse=True):
            gs, ge = hit_to_genome(h, context_start)
            pL = map_query_position_to_subject(h, attL_ref_pos)
            pR = map_query_position_to_subject(h, attR_ref_pos)
            pL = context_start + pL - 1 if pL is not None else ""
            pR = context_start + pR - 1 if pR is not None else ""

            w.writerow(
                [
                    gs,
                    ge,
                    h.subject_strand,
                    f"{h.pident:.3f}",
                    h.length,
                    f"{h.bitscore:.1f}",
                    f"{h.evalue:.3g}",
                    h.qstart,
                    h.qend,
                    pL,
                    pR,
                ]
            )


def main() -> int:
    args = parse_args()

    root = Path(args.root).resolve()
    manifest = (root / args.manifest).resolve()
    genomes_dir = root / "genomes"
    context_dir = root / "context"
    boundaries_dir = root / "boundaries"
    prophages_dir = root / "prophages"

    context_dir.mkdir(parents=True, exist_ok=True)
    boundaries_dir.mkdir(parents=True, exist_ok=True)
    prophages_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(manifest, sep="\t", dtype=str).fillna("")

    required = {
        "label",
        "organism",
        "accession",
        "type",
        "seed_start",
        "seed_end",
        "flank_bp",
        "boundary_mode",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Manifest is missing columns: {', '.join(sorted(missing))}")

    known = df[df["boundary_mode"] == "known_icd"]
    if len(known) != 1:
        raise ValueError(
            "Expected exactly one boundary_mode=known_icd row to define the e14 positive control"
        )

    known_row = known.iloc[0]
    known_label = known_row["label"]

    known_fasta = genomes_dir / known_label / f"{known_label}.fna"
    if not known_fasta.exists():
        raise FileNotFoundError(f"Missing {known_fasta}; run 01_download.sh first")

    k12_fasta_record = load_single_fasta(known_fasta)

    print("Deriving K-12 e14 icd-associated attachment reference...")
    (
        icd_reference,
        attL_ref_pos,
        attR_ref_pos,
        ref_meta,
    ) = derive_k12_reference(
        known_row,
        k12_fasta_record,
        boundaries_dir,
        args.reference_window,
    )

    print(
        f"  reference length: {len(icd_reference)} bp; "
        f"attL ref position={attL_ref_pos}; attR ref position={attR_ref_pos}"
    )
    print(
        "  K-12 paired blocks: "
        f"{ref_meta['left_block_start']}-{ref_meta['left_block_end']} and "
        f"{ref_meta['right_block_start']}-{ref_meta['right_block_end']}"
    )
    print()

    output_rows = []

    for _, row in df.iterrows():
        label = row["label"]
        organism = row["organism"]
        accession = row["accession"]
        seq_type = row["type"]
        mode = row["boundary_mode"]

        fasta_path = genomes_dir / label / f"{label}.fna"
        gbk_path = genomes_dir / label / f"{label}.gbk"

        if not fasta_path.exists() or not gbk_path.exists():
            raise FileNotFoundError(
                f"Missing downloaded files for {label}; run 01_download.sh first"
            )

        fasta_rec = load_single_fasta(fasta_path)
        gbk_rec = load_single_genbank(gbk_path)
        genome_len = len(fasta_rec.seq)

        print("=" * 72)
        print(label)
        print(f"  {organism}")
        print(f"  {accession} | {genome_len:,} bp | mode={mode}")

        # ---------------------------------------------------------------
        # Free phage: no attachment inference.
        # ---------------------------------------------------------------
        if seq_type == "free_phage" or mode == "none":
            extract_start = 1
            extract_end = genome_len

            outdir = prophages_dir / label
            outdir.mkdir(parents=True, exist_ok=True)
            write_fasta(
                fasta_rec.seq,
                outdir / f"{label}.fna",
                label,
                f"{accession} complete free-phage genome",
            )

            output_rows.append(
                {
                    "label": label,
                    "organism": organism,
                    "accession": accession,
                    "type": seq_type,
                    "genome_length": genome_len,
                    "seed_start": "",
                    "seed_end": "",
                    "context_start": 1,
                    "context_end": genome_len,
                    "method": "whole_free_phage",
                    "confidence": "not_applicable",
                    "icd_features_in_context": "",
                    "attL_block_start": "",
                    "attL_block_end": "",
                    "attL_identity": "",
                    "attL_length": "",
                    "attL_strand": "",
                    "attL_guess": "",
                    "attR_block_start": "",
                    "attR_block_end": "",
                    "attR_identity": "",
                    "attR_length": "",
                    "attR_strand": "",
                    "attR_guess": "",
                    "extract_start": extract_start,
                    "extract_end": extract_end,
                    "extract_length": extract_end - extract_start + 1,
                    "notes": "Complete isolated phage genome; no att-site inference",
                }
            )
            print("  complete free phage; no boundary inference")
            continue

        seed_start = as_int(row["seed_start"])
        seed_end = as_int(row["seed_end"])
        flank = as_int(row["flank_bp"]) or 10000

        if seed_start is None or seed_end is None:
            raise ValueError(f"{label}: prophage row requires seed_start and seed_end")

        if seed_start > seed_end:
            seed_start, seed_end = seed_end, seed_start

        context_start = max(1, seed_start - flank)
        context_end = min(genome_len, seed_end + flank)

        # Save broad host context for inspection.
        cdir = context_dir / label
        cdir.mkdir(parents=True, exist_ok=True)

        context_seq = fasta_rec.seq[context_start - 1 : context_end]
        write_fasta(
            context_seq,
            cdir / f"{label}_context.fna",
            f"{label}_context",
            f"{accession}:{context_start}-{context_end}",
        )

        context_gbk = slice_record(
            gbk_rec,
            context_start,
            context_end,
            f"{label}_context",
        )
        SeqIO.write(context_gbk, cdir / f"{label}_context.gbk", "genbank")

        icd_features = find_icd_features(gbk_rec, context_start, context_end)
        print(
            "  annotated icd/isocitrate-dehydrogenase features in context: "
            + (", ".join(icd_features) if icd_features else "none")
        )

        # ---------------------------------------------------------------
        # Canonical e14 positive control.
        # ---------------------------------------------------------------
        if mode == "known_icd":
            extract_start = seed_start
            extract_end = seed_end

            outdir = prophages_dir / label
            outdir.mkdir(parents=True, exist_ok=True)
            write_fasta(
                fasta_rec.seq[extract_start - 1 : extract_end],
                outdir / f"{label}.fna",
                label,
                f"{accession}:{extract_start}-{extract_end}; canonical e14",
            )

            output_rows.append(
                {
                    "label": label,
                    "organism": organism,
                    "accession": accession,
                    "type": seq_type,
                    "genome_length": genome_len,
                    "seed_start": seed_start,
                    "seed_end": seed_end,
                    "context_start": context_start,
                    "context_end": context_end,
                    "method": "known_icd_positive_control",
                    "confidence": "known",
                    "icd_features_in_context": ";".join(icd_features),
                    "attL_block_start": ref_meta["left_block_start"],
                    "attL_block_end": ref_meta["left_block_end"],
                    "attL_identity": f"{ref_meta['identity']:.3f}",
                    "attL_length": ref_meta["length"],
                    "attL_strand": "+",
                    "attL_guess": seed_start,
                    "attR_block_start": ref_meta["right_block_start"],
                    "attR_block_end": ref_meta["right_block_end"],
                    "attR_identity": f"{ref_meta['identity']:.3f}",
                    "attR_length": ref_meta["length"],
                    "attR_strand": "+",
                    "attR_guess": seed_end,
                    "extract_start": extract_start,
                    "extract_end": extract_end,
                    "extract_length": extract_end - extract_start + 1,
                    "notes": "Canonical e14 coordinates retained exactly",
                }
            )
            print(f"  known e14: {extract_start}-{extract_end}")
            continue

        # ---------------------------------------------------------------
        # Primary strategy: find two copies of the K-12 e14/icd boundary
        # homology in this host context.
        # ---------------------------------------------------------------
        with tempfile.TemporaryDirectory(prefix=f"{label}_icd_") as td:
            hits = run_blast(icd_reference, context_seq, Path(td))

        hits = [
            h
            for h in hits
            if h.length >= args.min_icd_hit_length
            and h.pident >= args.min_icd_identity
        ]

        write_candidate_hits(
            boundaries_dir / f"{label}_boundary_candidates.tsv",
            hits,
            context_start,
            attL_ref_pos,
            attR_ref_pos,
        )

        pair = choose_icd_pair(
            hits,
            context_start,
            seed_start,
            seed_end,
            attL_ref_pos,
            attR_ref_pos,
        )

        if pair is not None:
            left = pair["left"]
            right = pair["right"]

            attL_guess = left["attL_proj"]
            attR_guess = right["attR_proj"]

            # Conservative extraction includes the entire icd-associated
            # homology block at each end.
            extract_start = min(left["start"], left["end"])
            extract_end = max(right["start"], right["end"])

            confidence = confidence_for_icd_pair(pair, seed_start, seed_end)
            method = "e14_icd_reference_pair"

            notes = (
                "Conservative extraction uses outer edges of paired "
                "e14/icd-homologous blocks; projected att positions are best guesses"
            )

            out_left_identity = f"{left['hit'].pident:.3f}"
            out_right_identity = f"{right['hit'].pident:.3f}"
            out_left_len = left["hit"].length
            out_right_len = right["hit"].length
            out_left_strand = left["strand"]
            out_right_strand = right["strand"]

            print(
                f"  icd pair: {left['start']}-{left['end']} "
                f"({left['hit'].pident:.1f}%, {left['hit'].length} bp) | "
                f"{right['start']}-{right['end']} "
                f"({right['hit'].pident:.1f}%, {right['hit'].length} bp)"
            )
            print(f"  projected attL/attR: {attL_guess} / {attR_guess}")
            print(
                f"  conservative extraction: {extract_start}-{extract_end} "
                f"({extract_end - extract_start + 1:,} bp), confidence={confidence}"
            )

        else:
            # -----------------------------------------------------------
            # Fallback: directly compare broad left/right boundary windows.
            # -----------------------------------------------------------
            left_window_start = max(1, seed_start - flank)
            left_window_end = min(genome_len, seed_start + flank)
            right_window_start = max(1, seed_end - flank)
            right_window_end = min(genome_len, seed_end + flank)

            left_seq = fasta_rec.seq[left_window_start - 1 : left_window_end]
            right_seq = fasta_rec.seq[right_window_start - 1 : right_window_end]

            fallback = choose_fallback_pair(
                left_seq,
                right_seq,
                left_window_start,
                right_window_start,
                seed_start,
                seed_end,
                args.fallback_max_repeat,
            )

            if fallback is not None:
                left = fallback["left"]
                right = fallback["right"]
                h = fallback["hit"]

                extract_start = min(left["start"], left["end"])
                extract_end = max(right["start"], right["end"])

                attL_guess = round(left["centre"])
                attR_guess = round(right["centre"])

                method = "direct_boundary_homology_fallback"
                confidence = "medium" if h.length >= 50 and h.pident >= 80 else "low"
                notes = (
                    "No convincing paired e14/icd-reference hits; used direct "
                    "left-vs-right homology. att point estimates are block midpoints."
                )

                out_left_identity = f"{h.pident:.3f}"
                out_right_identity = f"{h.pident:.3f}"
                out_left_len = h.length
                out_right_len = h.length
                out_left_strand = "+"
                out_right_strand = h.subject_strand

                print(
                    f"  fallback homology: {left['start']}-{left['end']} | "
                    f"{right['start']}-{right['end']} "
                    f"({h.pident:.1f}%, {h.length} bp)"
                )
                print(
                    f"  conservative extraction: {extract_start}-{extract_end} "
                    f"({extract_end - extract_start + 1:,} bp), confidence={confidence}"
                )
            else:
                # Last resort: preserve seed interval and demand manual review.
                extract_start = seed_start
                extract_end = seed_end
                attL_guess = seed_start
                attR_guess = seed_end
                method = "seed_interval_fallback"
                confidence = "low"
                notes = (
                    "No convincing icd-reference pair or boundary repeat found; "
                    "seed interval retained and manual boundary review required"
                )

                left = {
                    "start": "",
                    "end": "",
                    "strand": "",
                }
                right = {
                    "start": "",
                    "end": "",
                    "strand": "",
                }
                out_left_identity = ""
                out_right_identity = ""
                out_left_len = ""
                out_right_len = ""
                out_left_strand = ""
                out_right_strand = ""

                print("  WARNING: no convincing boundary pair; retaining seed interval")

        # Sanity check.
        extract_start = int(extract_start)
        extract_end = int(extract_end)
        if extract_start < 1 or extract_end > genome_len or extract_start >= extract_end:
            raise RuntimeError(
                f"{label}: invalid inferred extraction {extract_start}-{extract_end}"
            )

        # Final candidate FASTA for the annotation stage.
        outdir = prophages_dir / label
        outdir.mkdir(parents=True, exist_ok=True)
        write_fasta(
            fasta_rec.seq[extract_start - 1 : extract_end],
            outdir / f"{label}.fna",
            label,
            f"{accession}:{extract_start}-{extract_end}; method={method}",
        )

        output_rows.append(
            {
                "label": label,
                "organism": organism,
                "accession": accession,
                "type": seq_type,
                "genome_length": genome_len,
                "seed_start": seed_start,
                "seed_end": seed_end,
                "context_start": context_start,
                "context_end": context_end,
                "method": method,
                "confidence": confidence,
                "icd_features_in_context": ";".join(icd_features),
                "attL_block_start": left["start"],
                "attL_block_end": left["end"],
                "attL_identity": out_left_identity,
                "attL_length": out_left_len,
                "attL_strand": out_left_strand,
                "attL_guess": attL_guess if attL_guess is not None else "",
                "attR_block_start": right["start"],
                "attR_block_end": right["end"],
                "attR_identity": out_right_identity,
                "attR_length": out_right_len,
                "attR_strand": out_right_strand,
                "attR_guess": attR_guess if attR_guess is not None else "",
                "extract_start": extract_start,
                "extract_end": extract_end,
                "extract_length": extract_end - extract_start + 1,
                "notes": notes,
            }
        )

    out_df = pd.DataFrame(output_rows)
    out_path = boundaries_dir / "boundaries.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)

    print()
    print("=" * 72)
    print(f"Wrote: {out_path}")
    print(f"Wrote final candidate FASTAs under: {prophages_dir}")
    print()
    print(
        "IMPORTANT: these are conservative best-guess boundaries. "
        "Review boundaries.tsv and the candidate-hit tables before annotation."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
