#!/usr/bin/env python3
"""Freeze two metadata-selected datasets before validation expression is read."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(".")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"

SELECTIONS = [
    {
        "dataset_id": "Populus trichocarpa_PRJCA014789",
        "panel": "vascular",
        "validation_label": "Phloem",
        "validation_role": "unseen_expression_species_validation",
        "selection_reason": "new expression species; official source resolved; 92.4% orthogroup coverage; 580 target cells; no discovery dataset or study overlap",
    },
    {
        "dataset_id": "Arabidopsis thaliana_PRJNA796288",
        "panel": "leaf",
        "validation_label": "Leaf guard cell",
        "validation_role": "independent_study_validation",
        "selection_reason": "independent stomatal-lineage study; official source resolved; 85.2% orthogroup coverage; 6254 target cells; no discovery dataset or study overlap",
    },
]


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    candidates = read_tsv(META / "independent_validation_candidates_v1.tsv")
    candidate_lookup = {(r["dataset_id"], r["panel"]): r for r in candidates}
    panels = read_tsv(META / "panel_dataset_freeze_v1.tsv")
    manifest = {r["dataset_id"]: r for r in read_tsv(META / "dataset_manifest.tsv")}
    source_rows = read_tsv(META / "dataset_source_metadata.tsv")
    source = {r["dataset_id"]: r for r in source_rows if r["status"] == "found"}
    accepted_datasets = {r["dataset_id"] for r in panels if r["accepted"] == "True"}
    accepted_studies = {r["study_group_id"] for r in panels if r["accepted"] == "True"}
    panel_lookup = {(r["dataset_id"], r["panel"]): r for r in panels}

    rows = []
    for selected in SELECTIONS:
        key = (selected["dataset_id"], selected["panel"])
        candidate = candidate_lookup[key]
        panel = panel_lookup[key]
        man = manifest[selected["dataset_id"]]
        src = source[selected["dataset_id"]]
        if selected["dataset_id"] in accepted_datasets or panel["study_group_id"] in accepted_studies:
            raise RuntimeError(f"discovery overlap detected for {selected['dataset_id']}")
        if src["species_matches_official"] != "True" or man["readable"] != "True":
            raise RuntimeError(f"source/readability gate failed for {selected['dataset_id']}")
        if candidate["orthogroup_coverage_fraction"] == "" or float(candidate["orthogroup_coverage_fraction"]) < 0.60:
            raise RuntimeError(f"orthogroup coverage gate failed for {selected['dataset_id']}")
        rows.append({
            **selected,
            "study_group_id": panel["study_group_id"],
            "species": panel["species"],
            "h5ad_path": man["h5ad_path"],
            "h5ad_file_bytes": man["file_bytes"],
            "source_expression_sha256": man["source_expression_sha256"],
            "source_annotation_sha256": man["source_annotation_sha256"],
            "expression_class": man["expression_class"],
            "target_cells_precount": panel["candidate_cells"],
            "orthogroup_coverage_fraction": candidate["orthogroup_coverage_fraction"],
            "official_repository": src["repository"],
            "official_accession": src["official_accession"],
            "official_title": src["title"],
            "discovery_dataset_overlap": "False",
            "discovery_study_overlap": "False",
            "expression_inspected_for_selection": "False",
            "marker_results_inspected_for_selection": "False",
            "freeze_status": "frozen_for_independent_validation_v1",
        })

    output = META / "independent_validation_freeze_v1.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "independent_validation_v1_frozen_before_expression_access",
        "datasets": len(rows),
        "roles": {row["dataset_id"]: row["validation_role"] for row in rows},
        "output": str(output),
        "output_sha256": sha256(output),
        "frozen_claims": [
            "Populus phloem tests vascular/phloem conserved programs in an expression species not used for discovery.",
            "Arabidopsis leaf guard cells test epidermal conserved programs in an independent study not used for discovery.",
        ],
        "primary_comparison": "complete orthogroup-family score versus single-copy orthogroup restriction, with matched random-orthogroup negative controls",
        "safety": {
            "source_h5ad_accessed_by_this_freeze": False,
            "selection_used_expression_values": False,
            "selection_used_marker_results": False,
            "discovery_freeze_modified": False,
            "other_task_written": False,
        },
    }
    audit_path = REPORTS / "independent_validation_freeze_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
