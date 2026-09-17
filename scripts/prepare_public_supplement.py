#!/usr/bin/env python3
"""Prepare reviewer-facing supplementary tables without machine-specific paths."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
REPO_SUPP = REPO / "results" / "supplementary"


def read_table(path: Path, delimiter: str = "\t") -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or []), list(reader)


def write_table(path: Path, fields: list[str], rows: list[dict[str, str]], delimiter: str = "\t") -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def basename_record(value: str) -> str:
    records: list[str] = []
    for item in filter(None, value.split(";")):
        name = item.replace("\\", "/").rstrip("/").split("/")[-1]
        if name in {"dataset_manifest.tsv", "audit_datasets.csv", "coarse_label_hierarchy_v1.tsv", "panel_label_freeze_v1.tsv", "conserved_program_candidates_v2.tsv"}:
            records.append(f"metadata/{name}")
        elif name.startswith("Table_S") or name in {"field_dictionary.tsv", "Methods_S1_detailed_methods.tex"}:
            records.append(f"results/supplementary/{name}")
        else:
            records.append(name)
    return ";".join(records)


def clean_metadata_files() -> tuple[str, str]:
    manifest = REPO / "metadata" / "dataset_manifest.tsv"
    fields, rows = read_table(manifest)
    for row in rows:
        for key in ("h5ad_path", "source_expression_path", "source_annotation_path"):
            row[key] = basename_record(row.get(key, ""))
    write_table(manifest, fields, rows)

    audit = REPO / "metadata" / "audit_datasets.csv"
    fields, rows = read_table(audit, delimiter=",")
    for row in rows:
        row["path"] = basename_record(row.get("path", ""))
    write_table(audit, fields, rows, delimiter=",")
    return sha256(manifest), sha256(audit)


def clean_supplement(directory: Path, manifest_sha: str, audit_sha: str) -> None:
    checks_path = directory / "Table_S1b_discovery_object_checksums.tsv"
    _, checks = read_table(checks_path)
    checksum_by_id = {row["dataset_id"]: row for row in checks}

    s1 = directory / "Table_S1_dataset_inventory.tsv"
    fields, rows = read_table(s1)
    if "analysis_role" not in fields:
        fields.insert(fields.index("dataset_role") + 1, "analysis_role")
    if "object_size_bytes" not in fields:
        fields.insert(fields.index("downloaded_object_filename") + 1, "object_size_bytes")

    for row in rows:
        original_role = row["dataset_role"]
        if original_role == "discovery_resource":
            row["dataset_role"] = "audited_resource"
            panels = [panel for panel in ("root", "leaf", "vascular") if row.get(f"{panel}_accepted") == "True"]
            if row["dataset_id"] == "Populus trichocarpa_PRJCA014789":
                row["analysis_role"] = "held_out_validation_input;excluded_from_discovery"
            elif panels:
                row["analysis_role"] = "discovery_input:" + ",".join(panels)
            else:
                row["analysis_role"] = "audited_not_used_in_discovery"
        else:
            row["analysis_role"] = "external_validation_input"

        check = checksum_by_id.get(row["dataset_id"])
        if check:
            row["downloaded_object_filename"] = check["standardized_object_filename"]
            row["object_size_bytes"] = check["object_size_bytes"]
            row["direct_source_object_path"] = check["standardized_object_filename"]
            row["direct_source_object_sha256"] = check["direct_source_object_sha256"]
            row["checksum_scope"] = "standardized_h5ad_object"
            gaps = [g for g in row.get("gap_flags", "").split(";") if g and g != "direct_source_object_path_and_sha256_not_available_in_manuscript_evidence_v1"]
            if not row.get("actual_download_url") and "public_raw_download_or_conversion_chain_not_resolved" not in gaps:
                gaps.append("public_raw_download_or_conversion_chain_not_resolved")
            row["gap_flags"] = ";".join(gaps)

        row["source_path"] = "metadata/audit_datasets.csv;metadata/dataset_manifest.tsv"
        row["source_sha256"] = f"{audit_sha};{manifest_sha}"
    write_table(s1, fields, rows)

    path_columns = {
        "Table_S2_label_hierarchy_dictionary.tsv": ["source_path"],
        "Table_S2_source_label_mapping.tsv": ["source_path"],
        "Table_S3_dataset_orthogroup_mapping_coverage.tsv": ["source_h5ad_path", "source_annotation_path", "mapping_source_path"],
        "Table_S4_xylem_orthogroup_program.tsv": ["source_path"],
        "Table_S5_replicate_metrics_controls_exclusions.tsv": ["source_path"],
    }
    for filename, columns in path_columns.items():
        path = directory / filename
        fields, rows = read_table(path)
        for row in rows:
            for column in columns:
                row[column] = basename_record(row.get(column, ""))
        write_table(path, fields, rows)


def sync_submission_to_repo(submission: Path) -> None:
    for path in submission.glob("*.tsv"):
        (REPO_SUPP / path.name).write_bytes(path.read_bytes())
    for name in ("Methods_S1_detailed_methods.tex", "README.md"):
        source = submission / name
        if source.exists():
            (REPO_SUPP / name).write_bytes(source.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submission-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest_sha, audit_sha = clean_metadata_files()
    clean_supplement(args.submission_dir, manifest_sha, audit_sha)
    sync_submission_to_repo(args.submission_dir)


if __name__ == "__main__":
    main()
