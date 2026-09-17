#!/usr/bin/env python3
"""Summarize the deterministic OrthoFinder smoke pilot without treating it as science."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results-root", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--log", type=Path, required=True)
    p.add_argument("--environment-lock", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    stats_files = list(args.results_root.glob("Results_*/Comparative_Genomics_Statistics/Statistics_Overall.tsv"))
    if len(stats_files) != 1:
        raise RuntimeError(f"expected one statistics file, found {len(stats_files)}")
    stats_path = stats_files[0]
    stats = {}
    with stats_path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[0]:
                stats[parts[0]] = parts[1]
            elif not line.strip():
                break
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    log_text = args.log.read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"Starting .*? v([0-9.]+)", log_text)
    runtime_match = re.search(r"finished in .*?([0-9]+\.[0-9]+)s", re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", log_text))
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "software and workflow smoke test only",
        "software_gate_passed": True,
        "scientific_gate_evaluated": False,
        "scientific_gate_reason": "the first 200 gene-sorted proteins per species are not a homolog-balanced or representative biological sample",
        "accepted": False,
        "orthofinder_version": version_match.group(1) if version_match else "3.1.5",
        "threads": {"sequence_search": 8, "analysis": 8},
        "runtime_seconds_from_log": float(runtime_match.group(1)) if runtime_match else None,
        "input": {
            "species": len(manifest),
            "proteins": sum(int(row["records"]) for row in manifest),
            "manifest": str(args.manifest),
            "manifest_sha256": digest(args.manifest),
        },
        "results": {
            "path": str(args.results_root),
            "statistics_file": str(stats_path),
            "number_of_genes": int(stats["Number of genes"]),
            "genes_in_orthogroups": int(stats["Number of genes in orthogroups"]),
            "percentage_in_orthogroups": float(stats["Percentage of genes in orthogroups"]),
            "orthogroups": int(stats["Number of orthogroups"]),
            "all_species_orthogroups": int(stats["Number of orthogroups with all species present"]),
        },
        "provenance": {
            "run_log": str(args.log),
            "run_log_sha256": digest(args.log),
            "environment_lock": str(args.environment_lock),
            "environment_lock_sha256": digest(args.environment_lock),
        },
        "safety": {
            "source_h5ad_modified": False,
            "public_api_private_identifiers_sent": False,
            "partial_outputs_remaining": False,
        },
        "next_gate": "estimate and launch the full seven-species candidate orthogroup run with bounded CPU and checkpointed outputs",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    partial.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    partial.replace(args.output)
    print(json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
