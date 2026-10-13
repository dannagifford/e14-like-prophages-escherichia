#!/usr/bin/env python3
"""
compare_e14_annotations.py

Line up the curated E. coli K-12 MG1655 CDS annotations against the CDS calls
made by Pharokka/PHANOTATE on the same e14 nucleotide sequence.

Default e14 interval:
    NC_000913.3:1196209-1211412  (1-based, inclusive)

The script:
  1. reads the original MG1655 GenBank annotation;
  2. extracts CDSs in the e14 interval and converts them to e14-relative coords;
  3. reads the Pharokka GenBank annotation for the extracted e14 sequence;
  4. matches CDSs by genomic overlap and strand;
  5. classifies exact matches, alternative starts/stops, frame changes,
     apparent splits and apparent merges;
  6. writes tabular outputs plus a two-track SVG with both annotations on the
     same nucleotide coordinate system;
  7. writes MG1655_curated_e14.gbk for use as an annotation-control track in Clinker.

If --reference-gbk is omitted, the script uses NCBI Entrez Direct `efetch`
to download NC_000913.3 in GenBank format.

Designed to run in the existing `prophage_compare` conda environment.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from copy import deepcopy
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


@dataclass
class CDS:
    source: str
    ident: str
    gene: str
    locus_tag: str
    product: str
    protein_id: str
    start: int          # 1-based inclusive, relative to e14
    end: int            # 1-based inclusive, relative to e14
    strand: int
    nt_length: int
    translation: str
    chrom_start: Optional[int] = None
    chrom_end: Optional[int] = None

    @property
    def aa_length(self) -> int:
        return len(self.translation.rstrip("*"))

    @property
    def coding_start(self) -> int:
        """Coordinate of the first translated nucleotide on the e14 axis."""
        return self.start if self.strand == 1 else self.end

    @property
    def coding_end(self) -> int:
        """Coordinate of the last translated nucleotide on the e14 axis."""
        return self.end if self.strand == 1 else self.start

    @property
    def frame_key(self) -> int:
        # Sufficient for comparing calls on the identical underlying sequence.
        if self.strand == 1:
            return (self.start - 1) % 3
        return (self.end - 1) % 3


def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument(
        "--pharokka-gbk",
        required=True,
        help="Pipeline GenBank for e14, with PHANOTATE CDS calls."
    )
    p.add_argument(
        "--reference-gbk",
        default=None,
        help="Full original MG1655 GenBank. If omitted, download with efetch."
    )
    p.add_argument("--accession", default="NC_000913.3")
    p.add_argument("--region-start", type=int, default=1196209)
    p.add_argument("--region-end", type=int, default=1211412)
    p.add_argument(
        "--outdir",
        default="results/e14_annotation_comparison"
    )
    p.add_argument(
        "--include-partial-reference-cds",
        action="store_true",
        help="Include MG1655 CDSs crossing an e14 boundary. By default only CDSs "
             "fully contained in the e14 interval are compared."
    )
    p.add_argument(
        "--allow-sequence-mismatch",
        action="store_true",
        help="Continue even if the Pharokka sequence is not exactly the same as "
             "the MG1655 e14 interval. Coordinate comparisons are unsafe if used."
    )
    return p.parse_args()


def q(feat, key: str) -> str:
    vals = feat.qualifiers.get(key, [])
    return ";".join(str(x) for x in vals)


def feature_translation(feat, record) -> str:
    tr = q(feat, "translation").replace(" ", "").replace("\n", "")
    if tr:
        return tr.rstrip("*")

    seq = feat.extract(record.seq)
    codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0])
    if codon_start > 1:
        seq = seq[codon_start - 1:]

    # Avoid Biopython warnings/errors if an annotation is not divisible by 3.
    usable = len(seq) - (len(seq) % 3)
    if usable <= 0:
        return ""
    return str(seq[:usable].translate(table=11, to_stop=False)).rstrip("*")


def pick_id(feat, fallback: str) -> str:
    for key in ("locus_tag", "gene", "protein_id", "ID"):
        val = q(feat, key)
        if val:
            return val
    return fallback


def download_reference(accession: str, out: Path) -> Path:
    efetch = shutil.which("efetch")
    if efetch is None:
        raise RuntimeError(
            "No --reference-gbk was supplied and `efetch` is not in PATH.\n"
            "Either supply the original NC_000913.3 GenBank file or run this "
            "inside the `prophage_compare` environment, where Entrez Direct "
            "was installed by the project setup."
        )

    out.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {accession} GenBank annotation -> {out}", file=sys.stderr)
    with out.open("w") as fh:
        subprocess.run(
            [efetch, "-db", "nuccore", "-id", accession, "-format", "gbwithparts"],
            stdout=fh,
            check=True
        )
    return out


def one_record(path: Path):
    records = list(SeqIO.parse(path, "genbank"))
    if len(records) != 1:
        raise ValueError(f"{path}: expected one GenBank record, found {len(records)}")
    return records[0]


def reference_cds(record, region_start: int, region_end: int,
                  include_partial: bool) -> list[CDS]:
    rows = []
    n = 0

    for feat in record.features:
        if feat.type != "CDS":
            continue

        chrom_start = int(feat.location.start) + 1
        chrom_end = int(feat.location.end)
        strand = int(feat.location.strand or 0)

        overlaps = max(chrom_start, region_start) <= min(chrom_end, region_end)
        contained = chrom_start >= region_start and chrom_end <= region_end

        if not overlaps:
            continue
        if not include_partial and not contained:
            continue

        n += 1
        rel_start = chrom_start - region_start + 1
        rel_end = chrom_end - region_start + 1

        # Clip only for display/overlap if a boundary-crossing CDS was requested.
        display_start = max(1, rel_start)
        display_end = min(region_end - region_start + 1, rel_end)

        tr = feature_translation(feat, record)

        rows.append(CDS(
            source="MG1655",
            ident=pick_id(feat, f"MG1655_CDS_{n:04d}"),
            gene=q(feat, "gene"),
            locus_tag=q(feat, "locus_tag"),
            product=q(feat, "product"),
            protein_id=q(feat, "protein_id"),
            start=display_start,
            end=display_end,
            strand=strand,
            nt_length=chrom_end - chrom_start + 1,
            translation=tr,
            chrom_start=chrom_start,
            chrom_end=chrom_end,
        ))

    return sorted(rows, key=lambda x: (x.start, x.end, x.strand))


def pharokka_cds(record) -> list[CDS]:
    rows = []
    n = 0

    for feat in record.features:
        if feat.type != "CDS":
            continue
        n += 1
        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        strand = int(feat.location.strand or 0)
        tr = feature_translation(feat, record)

        rows.append(CDS(
            source="PHANOTATE",
            ident=pick_id(feat, f"PHANOTATE_CDS_{n:04d}"),
            gene=q(feat, "gene"),
            locus_tag=q(feat, "locus_tag"),
            product=q(feat, "product"),
            protein_id=q(feat, "protein_id"),
            start=start,
            end=end,
            strand=strand,
            nt_length=end - start + 1,
            translation=tr,
        ))

    return sorted(rows, key=lambda x: (x.start, x.end, x.strand))



def write_curated_e14_genbank(record, region_start: int, region_end: int, out_path: Path) -> int:
    """Write e14 with the original NC_000913.3 CDS annotations retained."""
    start0 = region_start - 1
    end0 = region_end

    out = SeqRecord(
        record.seq[start0:end0],
        id="MG1655_curated_e14",
        name="MG1655_curated_e14",
        description=(
            f"E. coli K-12 MG1655 e14 ({record.id}:{region_start}-{region_end}); "
            "original reference CDS annotation"
        ),
    )
    out.annotations["molecule_type"] = "DNA"
    out.annotations["topology"] = "linear"

    n = 0
    for feat in record.features:
        if feat.type != "CDS":
            continue
        fstart0 = int(feat.location.start)
        fend0 = int(feat.location.end)

        # Retain only complete reference CDSs. Do not clip boundary-crossing genes.
        if fstart0 < start0 or fend0 > end0:
            continue

        copied_feature = deepcopy(feat)
        copied_feature.location = feat.location + (-start0)
        out.features.append(copied_feature)
        n += 1

    SeqIO.write(out, out_path, "genbank")
    return n


def overlap_bp(a: CDS, b: CDS) -> int:
    if a.strand != b.strand:
        return 0
    return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)


def translation_relation(a: CDS, b: CDS) -> str:
    x = a.translation.rstrip("*")
    y = b.translation.rstrip("*")
    if not x or not y:
        return "translation_unavailable"
    if x == y:
        return "exact_translation"
    if x in y:
        return "MG1655_translation_is_subsequence"
    if y in x:
        return "PHANOTATE_translation_is_subsequence"
    return "different_translation"


def classify_one_to_one(a: CDS, b: CDS) -> str:
    if a.start == b.start and a.end == b.end and a.strand == b.strand:
        return "exact_gene_model"

    same_frame = a.frame_key == b.frame_key
    same_coding_start = a.coding_start == b.coding_start
    same_coding_end = a.coding_end == b.coding_end

    if same_frame and same_coding_end and not same_coding_start:
        return "alternative_start_same_frame"
    if same_frame and same_coding_start and not same_coding_end:
        return "alternative_stop_same_frame"
    if same_frame:
        return "boundary_difference_same_frame"
    return "different_frame_overlap"


def interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    s, e = intervals[0]
    for ns, ne in intervals[1:]:
        if ns <= e + 1:
            e = max(e, ne)
        else:
            total += e - s + 1
            s, e = ns, ne
    total += e - s + 1
    return total


def overlaps_for(query: CDS, targets: list[CDS]) -> list[tuple[CDS, int]]:
    hits = []
    for t in targets:
        ov = overlap_bp(query, t)
        if ov > 0:
            hits.append((t, ov))
    return sorted(hits, key=lambda x: (-x[1], x[0].start))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def esc(s: str) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
                  .replace(">", "&gt;").replace('"', "&quot;"))


def gene_label(c: CDS, pharokka=False) -> str:
    if pharokka:
        short = c.ident.replace("Ecoli_K12_e14_", "")
        if c.product:
            return f"{short} | {c.product}"
        return short
    return c.gene or c.locus_tag or c.protein_id or c.ident


def arrow_points(c: CDS, x0: float, scale: float, y: float, h: float) -> str:
    left = x0 + (c.start - 1) * scale
    right = x0 + c.end * scale
    head = min(10.0, max(4.0, (right - left) * 0.25))
    mid = y + h / 2

    if c.strand == -1:
        pts = [
            (left, mid),
            (left + head, y),
            (right, y),
            (right, y + h),
            (left + head, y + h),
        ]
    else:
        pts = [
            (left, y),
            (right - head, y),
            (right, mid),
            (right - head, y + h),
            (left, y + h),
        ]
    return " ".join(f"{x:.1f},{yy:.1f}" for x, yy in pts)


def make_svg(path: Path, ref: list[CDS], ph: list[CDS], region_len: int):
    width = 1800
    left_margin = 110
    right_margin = 40
    plot_w = width - left_margin - right_margin
    scale = plot_w / region_len
    height = 520

    # Alternating label rows reduce collisions without trying to solve a full
    # label-placement problem.
    out = []
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    )
    out.append('<rect width="100%" height="100%" fill="white"/>')
    out.append(
        '<style>'
        'text{font-family:Arial,sans-serif;fill:#111}'
        '.track{font-size:16px;font-weight:bold}'
        '.gene{fill:none;stroke:#111;stroke-width:1.4}'
        '.label{font-size:10px}'
        '.axis{stroke:#444;stroke-width:1}'
        '.tick{stroke:#777;stroke-width:1}'
        '.small{font-size:10px}'
        '</style>'
    )

    # Axis
    axis_y = 55
    out.append(
        f'<line class="axis" x1="{left_margin}" y1="{axis_y}" '
        f'x2="{left_margin + plot_w}" y2="{axis_y}"/>'
    )
    tick = 1000
    for pos in range(0, region_len + 1, tick):
        x = left_margin + pos * scale
        out.append(f'<line class="tick" x1="{x:.1f}" y1="{axis_y-4}" x2="{x:.1f}" y2="{axis_y+4}"/>')
        out.append(f'<text class="small" x="{x:.1f}" y="{axis_y-9}" text-anchor="middle">{pos}</text>')
    out.append(f'<text class="small" x="{left_margin + plot_w/2:.1f}" y="22" text-anchor="middle">e14-relative nucleotide coordinate</text>')

    tracks = [
        ("MG1655 curated", ref, 130, False),
        ("Pharokka / PHANOTATE", ph, 335, True),
    ]

    for title, genes, y, is_ph in tracks:
        out.append(f'<text class="track" x="12" y="{y+15}">{esc(title)}</text>')
        out.append(
            f'<line class="axis" x1="{left_margin}" y1="{y+14}" '
            f'x2="{left_margin + plot_w}" y2="{y+14}"/>'
        )

        for i, c in enumerate(genes):
            gy = y
            out.append(f'<polygon class="gene" points="{arrow_points(c, left_margin, scale, gy, 28)}"/>')

            centre = left_margin + (((c.start + c.end) / 2) - 1) * scale
            # alternate above/below labels
            if i % 2 == 0:
                ly = gy - 8
                anchor = "end"
                rot = -38
            else:
                ly = gy + 52
                anchor = "start"
                rot = 38

            label = gene_label(c, pharokka=is_ph)
            out.append(
                f'<text class="label" x="{centre:.1f}" y="{ly}" '
                f'text-anchor="{anchor}" transform="rotate({rot} {centre:.1f} {ly})">'
                f'{esc(label)}</text>'
            )

    out.append('</svg>')
    path.write_text("\n".join(out))


def main():
    args = parse_args()

    ph_path = Path(args.pharokka_gbk)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not ph_path.exists():
        raise FileNotFoundError(ph_path)

    if args.reference_gbk:
        ref_path = Path(args.reference_gbk)
        if not ref_path.exists():
            raise FileNotFoundError(ref_path)
    else:
        ref_path = outdir / f"{args.accession}.gbk"
        if not ref_path.exists():
            download_reference(args.accession, ref_path)

    ref_record = one_record(ref_path)
    ph_record = one_record(ph_path)

    region_len = args.region_end - args.region_start + 1
    if region_len <= 0:
        raise ValueError("--region-end must be >= --region-start")

    ref_region = ref_record.seq[args.region_start - 1:args.region_end]
    ph_seq = ph_record.seq

    seq_exact = str(ref_region).upper() == str(ph_seq).upper()
    seq_revcomp = str(ref_region.reverse_complement()).upper() == str(ph_seq).upper()

    if not seq_exact:
        msg = (
            f"Sequence check failed: MG1655 interval length={len(ref_region)}, "
            f"Pipeline record length={len(ph_seq)}, exact_match={seq_exact}, "
            f"reverse_complement_match={seq_revcomp}."
        )
        if not args.allow_sequence_mismatch:
            raise RuntimeError(
                msg + "\nThe annotations cannot safely be placed on the same "
                "coordinate axis. Use --allow-sequence-mismatch only if you "
                "understand why the sequences differ."
            )
        print("WARNING:", msg, file=sys.stderr)

    curated_gbk = outdir / "MG1655_curated_e14.gbk"
    n_curated_gbk = write_curated_e14_genbank(
        ref_record,
        args.region_start,
        args.region_end,
        curated_gbk,
    )

    ref = reference_cds(
        ref_record,
        args.region_start,
        args.region_end,
        args.include_partial_reference_cds,
    )
    ph = pharokka_cds(ph_record)

    pair_rows = []
    ref_summary = []

    # Precompute reverse mapping so we can identify apparent merges.
    ph_to_ref = {x.ident: overlaps_for(x, ref) for x in ph}

    for a in ref:
        hits = overlaps_for(a, ph)
        hit_intervals = [
            (max(a.start, b.start), min(a.end, b.end))
            for b, ov in hits if ov > 0
        ]
        covered = interval_union_length(hit_intervals)
        ref_cov = covered / (a.end - a.start + 1) if a.end >= a.start else 0

        if not hits:
            overall = "MG1655_only"
        elif len(hits) > 1 and ref_cov >= 0.80:
            overall = "MG1655_CDS_split_by_PHANOTATE"
        elif len(hits) > 1:
            overall = "multiple_PHANOTATE_overlaps"
        else:
            overall = classify_one_to_one(a, hits[0][0])

        ref_summary.append({
            "mg1655_id": a.ident,
            "gene": a.gene,
            "locus_tag": a.locus_tag,
            "protein_id": a.protein_id,
            "mg1655_product": a.product,
            "strand": a.strand,
            "mg1655_rel_start": a.start,
            "mg1655_rel_end": a.end,
            "mg1655_chr_start": a.chrom_start,
            "mg1655_chr_end": a.chrom_end,
            "mg1655_nt_length": a.nt_length,
            "mg1655_aa_length": a.aa_length,
            "n_overlapping_phanotate_cds": len(hits),
            "phanotate_ids": ";".join(b.ident for b, _ in hits),
            "phanotate_products": ";".join(b.product for b, _ in hits),
            "reference_fraction_covered": f"{ref_cov:.4f}",
            "classification": overall,
        })

        for b, ov in hits:
            b_hits = ph_to_ref[b.ident]
            merge = len(b_hits) > 1
            pair_rows.append({
                "mg1655_id": a.ident,
                "gene": a.gene,
                "locus_tag": a.locus_tag,
                "protein_id": a.protein_id,
                "mg1655_product": a.product,
                "phanotate_id": b.ident,
                "phanotate_product": b.product,
                "strand": a.strand,
                "mg1655_rel_start": a.start,
                "mg1655_rel_end": a.end,
                "phanotate_rel_start": b.start,
                "phanotate_rel_end": b.end,
                "mg1655_aa_length": a.aa_length,
                "phanotate_aa_length": b.aa_length,
                "overlap_bp": ov,
                "mg1655_fraction_overlapped": f"{ov / (a.end-a.start+1):.4f}",
                "phanotate_fraction_overlapped": f"{ov / (b.end-b.start+1):.4f}",
                "same_frame": int(a.frame_key == b.frame_key),
                "same_coding_start": int(a.coding_start == b.coding_start),
                "same_coding_end": int(a.coding_end == b.coding_end),
                "pair_classification": classify_one_to_one(a, b),
                "translation_relation": translation_relation(a, b),
                "phanotate_overlaps_multiple_mg1655_cds": int(merge),
            })

    unmatched_ph = []
    for b in ph:
        hits = ph_to_ref[b.ident]
        if not hits:
            unmatched_ph.append({
                "phanotate_id": b.ident,
                "product": b.product,
                "strand": b.strand,
                "rel_start": b.start,
                "rel_end": b.end,
                "nt_length": b.nt_length,
                "aa_length": b.aa_length,
            })

    write_tsv(
        outdir / "mg1655_vs_phanotate_summary.tsv",
        [
            "mg1655_id", "gene", "locus_tag", "protein_id", "mg1655_product",
            "strand", "mg1655_rel_start", "mg1655_rel_end",
            "mg1655_chr_start", "mg1655_chr_end",
            "mg1655_nt_length", "mg1655_aa_length",
            "n_overlapping_phanotate_cds", "phanotate_ids",
            "phanotate_products", "reference_fraction_covered",
            "classification"
        ],
        ref_summary
    )

    write_tsv(
        outdir / "overlapping_cds_pairs.tsv",
        [
            "mg1655_id", "gene", "locus_tag", "protein_id", "mg1655_product",
            "phanotate_id", "phanotate_product", "strand",
            "mg1655_rel_start", "mg1655_rel_end",
            "phanotate_rel_start", "phanotate_rel_end",
            "mg1655_aa_length", "phanotate_aa_length",
            "overlap_bp", "mg1655_fraction_overlapped",
            "phanotate_fraction_overlapped", "same_frame",
            "same_coding_start", "same_coding_end",
            "pair_classification", "translation_relation",
            "phanotate_overlaps_multiple_mg1655_cds"
        ],
        pair_rows
    )

    write_tsv(
        outdir / "phanotate_only.tsv",
        [
            "phanotate_id", "product", "strand", "rel_start", "rel_end",
            "nt_length", "aa_length"
        ],
        unmatched_ph
    )

    make_svg(
        outdir / "annotation_tracks.svg",
        ref,
        ph,
        region_len
    )

    counts = {}
    for r in ref_summary:
        counts[r["classification"]] = counts.get(r["classification"], 0) + 1

    report = [
        "e14 annotation comparison",
        "========================",
        f"Reference: {ref_path}",
        f"Pipeline annotation: {ph_path}",
        f"Region:    {args.accession}:{args.region_start}-{args.region_end}",
        f"Region bp: {region_len}",
        f"Reference sequence == pipeline sequence: {seq_exact}",
        f"Reference sequence reverse-complement == pipeline sequence: {seq_revcomp}",
        "",
        f"MG1655 CDSs compared:   {len(ref)}",
        f"PHANOTATE CDSs:         {len(ph)}",
        f"PHANOTATE-only CDSs:    {len(unmatched_ph)}",
        f"Curated Clinker-control CDSs: {n_curated_gbk}",
        f"Curated Clinker control: {curated_gbk}",
        "",
        "MG1655-centred classifications:",
    ]
    for k in sorted(counts):
        report.append(f"  {k}: {counts[k]}")

    split = [r for r in ref_summary if r["classification"] == "MG1655_CDS_split_by_PHANOTATE"]
    report.extend(["", "Potential MG1655 CDSs split into multiple PHANOTATE calls:"])
    if split:
        for r in split:
            report.append(
                f"  {r['mg1655_id']} ({r['gene'] or r['mg1655_product']}) -> "
                f"{r['phanotate_ids']}"
            )
    else:
        report.append("  none")

    merges = {}
    for b in ph:
        hits = ph_to_ref[b.ident]
        if len(hits) > 1:
            merges[b.ident] = ";".join(a.ident for a, _ in hits)

    report.extend(["", "Potential PHANOTATE CDSs spanning multiple MG1655 CDSs:"])
    if merges:
        for pid, refs in sorted(merges.items()):
            report.append(f"  {pid} -> {refs}")
    else:
        report.append("  none")

    (outdir / "report.txt").write_text("\n".join(report) + "\n")

    print("\n".join(report))
    print()
    print("Outputs:")
    for name in [
        "mg1655_vs_phanotate_summary.tsv",
        "overlapping_cds_pairs.tsv",
        "phanotate_only.tsv",
        "annotation_tracks.svg",
        "MG1655_curated_e14.gbk",
        "report.txt",
    ]:
        print(f"  {outdir / name}")


if __name__ == "__main__":
    main()
