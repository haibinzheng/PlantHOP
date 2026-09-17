#!/usr/bin/env python3
"""Select metadata-only candidates for a frozen independent validation cohort.

Selection never inspects expression values or marker results. A dataset is
eligible only if neither it nor its study group contributed to any accepted
discovery panel.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(".")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    panels = read_tsv(META / "panel_dataset_freeze_v1.tsv")
    manifest = {r["dataset_id"]: r for r in read_tsv(META / "dataset_manifest.tsv")}
    coverage = {r["dataset_id"]: r for r in read_tsv(META / "dataset_orthogroup_coverage_v1.tsv")}
    source_rows = read_tsv(META / "dataset_source_metadata.tsv")
    source = {}
    for row in source_rows:
        # Prefer a resolved official record when duplicate dataset rows exist.
        current = source.get(row["dataset_id"])
        if current is None or (current["status"] != "found" and row["status"] == "found"):
            source[row["dataset_id"]] = row

    accepted_datasets = {r["dataset_id"] for r in panels if r["accepted"] == "True"}
    accepted_studies = {r["study_group_id"] for r in panels if r["accepted"] == "True"}
    candidates = []
    for row in panels:
        dataset = row["dataset_id"]
        if dataset in accepted_datasets or row["study_group_id"] in accepted_studies:
            continue
        labels = int(row["labels_ge_100_cells"])
        cells = int(row["candidate_cells"])
        cov = coverage.get(dataset)
        src = source.get(dataset)
        man = manifest.get(dataset)
        reasons = []
        if labels < 2:
            reasons.append("fewer_than_two_labels_ge_100")
        if cells < 500:
            reasons.append("fewer_than_500_panel_cells")
        if cov is None or cov["passes_provisional_60pct_gate"] != "True":
            reasons.append("orthogroup_coverage_below_60pct_or_missing")
        if src is None or src["status"] != "found":
            reasons.append("official_source_not_resolved")
        if src is None or src["species_matches_official"] != "True":
            reasons.append("official_species_not_confirmed")
        if man is None or man["readable"] != "True":
            reasons.append("h5ad_not_readable")
        eligible = not reasons
        candidates.append({
            "dataset_id": dataset,
            "study_group_id": row["study_group_id"],
            "species": row["species"],
            "panel": row["panel"],
            "labels_ge_100_cells": str(labels),
            "candidate_cells": str(cells),
            "orthogroup_coverage_fraction": cov["orthogroup_coverage_fraction"] if cov else "",
            "expression_class": man["expression_class"] if man else "",
            "official_repository": src["repository"] if src else "",
            "official_accession": src["official_accession"] if src else "",
            "original_exclusion_reason": row["reason"],
            "metadata_only_eligible": str(eligible),
            "ineligibility_reasons": ";".join(reasons),
            "discovery_dataset_overlap": "False",
            "discovery_study_overlap": "False",
        })

    candidates.sort(key=lambda r: (
        r["metadata_only_eligible"] != "True",
        -int(r["labels_ge_100_cells"]),
        -int(r["candidate_cells"]),
        r["dataset_id"], r["panel"],
    ))
    out = META / "independent_validation_candidates_v1.tsv"
    write_tsv(out, candidates)
    eligible = [r for r in candidates if r["metadata_only_eligible"] == "True"]
    by_panel = defaultdict(int)
    by_species = defaultdict(int)
    for row in eligible:
        by_panel[row["panel"]] += 1
        by_species[row["species"]] += 1
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "metadata_only_validation_candidate_screen_complete",
        "selection_used_expression_values": False,
        "selection_used_marker_results": False,
        "discovery_datasets": len(accepted_datasets),
        "discovery_study_groups": len(accepted_studies),
        "candidate_panel_rows": len(candidates),
        "eligible_panel_rows": len(eligible),
        "eligible_by_panel": dict(by_panel),
        "eligible_by_species": dict(by_species),
        "output": str(out),
        "output_sha256": sha256(out),
        "next_gate": "Freeze one or two validation datasets using only metadata, label support, provenance, and mapping coverage before reading their expression values.",
        "safety": {
            "source_h5ad_accessed": False,
            "model_training_started": False,
            "external_queries": False,
            "discovery_freeze_modified": False,
        },
    }
    audit_path = REPORTS / "independent_validation_candidate_screen_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    print("\nTOP ELIGIBLE")
    for row in eligible[:20]:
        print("\t".join([row["dataset_id"], row["panel"], row["species"], row["labels_ge_100_cells"], row["candidate_cells"], row["orthogroup_coverage_fraction"], row["expression_class"], row["official_repository"]]))


if __name__ == "__main__":
    main()
