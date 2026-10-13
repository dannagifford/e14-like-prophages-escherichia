#!/usr/bin/env python3
import csv
import os
from pathlib import Path

src = Path(os.environ["BASE_MANIFEST"])
dst = Path(os.environ["AUG_MANIFEST"])
extra_label = os.environ["EXTRA_LABEL"]
e14_label = os.environ.get("E14_LABEL", "Ecoli_K12_e14")

with src.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)

if not fields or "label" not in fields:
    raise SystemExit("ERROR: manifest.tsv must contain a 'label' column")

if any((r.get("label") or "").strip() == extra_label for r in rows):
    raise SystemExit(
        f"ERROR: {extra_label} is already present in manifest.tsv. "
        "It should remain an annotation-control track, not a main manifest entry."
    )

new_rows = []
inserted = False

for row in rows:
    new_rows.append(row)

    if (row.get("label") or "").strip() == e14_label:
        extra = {field: "" for field in fields}
        extra["label"] = extra_label

        if "accession" in extra:
            extra["accession"] = "NC_000913.3"
        if "type" in extra:
            extra["type"] = "annotation_control"
        if "organism" in extra:
            extra["organism"] = "Escherichia coli K-12 MG1655"
        if "description" in extra:
            extra["description"] = (
                "Original NC_000913.3 CDS annotation of the same e14 sequence"
            )

        new_rows.append(extra)
        inserted = True

if not inserted:
    raise SystemExit(f"ERROR: {e14_label} was not found in manifest.")

dst.parent.mkdir(parents=True, exist_ok=True)

with dst.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(new_rows)

print(f"Wrote e14-augmented manifest: {dst}")
