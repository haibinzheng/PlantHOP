from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "01_投稿材料" / "TEX" / "supplementary"


def portable_text(value: str) -> str:
    value = value.replace("\\", "/")
    project_prefix = str(ROOT).replace("\\", "/") + "/"
    value = value.replace(project_prefix, "project_record:")
    value = value.replace(
        "/workspace/projects/PhyloOpenCell/",
        "internal_source_record:not_publicly_portable:workspace/projects/PhyloOpenCell/",
    )
    value = value.replace(
        "/data/",
        "internal_source_record:not_publicly_portable:data/",
    )
    value = value.replace("audit_20260910", "resource_inventory_20260910")
    value = value.replace("metadata_and_audit_products_only", "metadata_and_derived_records_only")
    value = value.replace(
        "P1 authority package",
        "available annotation and mapping evidence",
    )
    value = value.replace(
        "Frozen authority version recorded by manuscript_evidence_v1.",
        "Provenance record version used for the submitted evidence table.",
    )
    return value


def rewrite_tsv(path: Path) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames is None:
            return
        fields = ["provenance_version" if f == "authority_version" else f for f in reader.fieldnames]

    normalized: list[dict[str, str]] = []
    for row in rows:
        out: dict[str, str] = {}
        for key, value in row.items():
            new_key = "provenance_version" if key == "authority_version" else key
            new_value = portable_text(value or "")
            if new_key == "field" and new_value == "authority_version":
                new_value = "provenance_version"
            out[new_key] = new_value
        normalized.append(out)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(normalized)


def main() -> None:
    for path in sorted(SUPP.glob("*.tsv")):
        rewrite_tsv(path)


if __name__ == "__main__":
    main()
