#!/usr/bin/env python3
"""Convert completed OrthoFinder output into auditable gene and dataset candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


HEADER_TO_SPECIES = {
    "arabidopsis_thaliana.representative_proteins": "Arabidopsis thaliana",
    "brassica_rapa.representative_proteins": "Brassica rapa",
    "catharanthus_roseus.representative_proteins": "Catharanthus roseus",
    "glycine_max.representative_proteins": "Glycine max",
    "oryza_sativa.representative_proteins": "Oryza sativa",
    "populus_trichocarpa.representative_proteins": "Populus trichocarpa",
    "sorghum_bicolor.representative_proteins": "Sorghum bicolor",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_groups(path: Path, status: str, records: list[dict], seen: set[tuple[str, str]]) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        species_headers = [h for h in headers if h != "Orthogroup"]
        if set(species_headers) != set(HEADER_TO_SPECIES):
            raise RuntimeError(f"unexpected species columns in {path}")
        for row in reader:
            members = []
            present_species = 0
            for header in species_headers:
                tokens = [x.strip() for x in row[header].split(",") if x.strip()]
                if tokens:
                    present_species += 1
                for token in tokens:
                    if "|" not in token:
                        raise RuntimeError(f"missing gene/protein separator: {token}")
                    gene, protein = token.split("|", 1)
                    key = (HEADER_TO_SPECIES[header], gene)
                    if key in seen:
                        raise RuntimeError(f"duplicate gene across orthogroups: {key}")
                    seen.add(key)
                    members.append((HEADER_TO_SPECIES[header], gene, protein))
            group_size = len(members)
            for species, gene, protein in members:
                records.append({
                    "species": species,
                    "candidate_reference_gene_id": gene,
                    "representative_protein_accession": protein,
                    "orthogroup_id": row["Orthogroup"],
                    "orthogroup_status": status,
                    "orthogroup_species_count": present_species,
                    "orthogroup_gene_count": group_size,
                    "accepted": "False",
                })


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--gene-mapping", type=Path, required=True)
    p.add_argument("--run-log", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--dataset-summary", type=Path, required=True)
    p.add_argument("--audit", type=Path, required=True)
    args = p.parse_args()

    assigned_path = args.results_dir / "Orthogroups" / "Orthogroups.tsv"
    unassigned_path = args.results_dir / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
    stats_path = args.results_dir / "Comparative_Genomics_Statistics" / "Statistics_Overall.tsv"
    tree_path = args.results_dir / "Species_Tree" / "SpeciesTree_rooted.txt"
    records: list[dict] = []
    seen: set[tuple[str, str]] = set()
    read_groups(assigned_path, "assigned_candidate", records, seen)
    read_groups(unassigned_path, "unassigned_singleton", records, seen)
    records.sort(key=lambda r: (r["species"], r["candidate_reference_gene_id"]))
    assigned_lookup = {
        (r["species"], r["candidate_reference_gene_id"]): r["orthogroup_id"]
        for r in records if r["orthogroup_status"] == "assigned_candidate"
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = args.output.with_suffix(args.output.suffix + ".partial")
    with gzip.open(output_tmp, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    os.replace(output_tmp, args.output)

    dataset_counts = defaultdict(lambda: {"features": 0, "gene_mapped": 0, "orthogroup_assigned": 0})
    dataset_species = {}
    with gzip.open(args.gene_mapping, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row["dataset_id"]
            dataset_species[key] = row["species"]
            counts = dataset_counts[key]
            counts["features"] += 1
            if row["candidate_reference_gene_id"]:
                counts["gene_mapped"] += 1
            if (row["species"], row["candidate_reference_gene_id"]) in assigned_lookup:
                counts["orthogroup_assigned"] += 1
    summary_rows = []
    for dataset_id in sorted(dataset_counts):
        c = dataset_counts[dataset_id]
        total = c["features"]
        summary_rows.append({
            "dataset_id": dataset_id,
            "species": dataset_species[dataset_id],
            "features": total,
            "gene_mapped_candidates": c["gene_mapped"],
            "orthogroup_assigned_candidates": c["orthogroup_assigned"],
            "orthogroup_coverage_fraction": c["orthogroup_assigned"] / total if total else 0.0,
            "passes_provisional_60pct_gate": str(c["orthogroup_assigned"] / total >= 0.60 if total else False),
            "accepted": "False",
        })
    args.dataset_summary.parent.mkdir(parents=True, exist_ok=True)
    summary_tmp = args.dataset_summary.with_suffix(args.dataset_summary.suffix + ".partial")
    with summary_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    os.replace(summary_tmp, args.dataset_summary)

    stats = {}
    with stats_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 2 and fields[0]:
                stats[fields[0]] = fields[1]
            elif not line.strip():
                break
    log_text = args.run_log.read_text(encoding="utf-8", errors="replace")
    clean_log = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", log_text)
    elapsed = re.search(r"Elapsed \(wall clock\) time.*\):\s*(\S+)", clean_log)
    peak = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", clean_log)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "state": "complete_candidate",
        "accepted": False,
        "software": {"orthofinder": "3.1.5", "sequence_search": "diamond", "tree_method": "msa/famsa/fasttree"},
        "inputs": {"species": 7, "proteins": len(records), "gene_mapping_sha256": sha256(args.gene_mapping)},
        "results": {
            "genes_in_orthogroups": int(stats["Number of genes in orthogroups"]),
            "unassigned_genes": int(stats["Number of unassigned genes"]),
            "percentage_in_orthogroups": float(stats["Percentage of genes in orthogroups"]),
            "orthogroups": int(stats["Number of orthogroups"]),
            "species_specific_orthogroups": int(stats["Number of species-specific orthogroups"]),
            "all_species_orthogroups": int(stats["Number of orthogroups with all species present"]),
            "single_copy_orthogroups": int(stats["Number of single-copy orthogroups"]),
            "datasets_evaluated": len(summary_rows),
            "datasets_passing_60pct_gate": sum(r["passes_provisional_60pct_gate"] == "True" for r in summary_rows),
        },
        "runtime": {"elapsed_wall_clock": elapsed.group(1) if elapsed else None, "peak_rss_kbytes": int(peak.group(1)) if peak else None, "threads": 12},
        "artifacts": {
            "gene_orthogroup_candidates": str(args.output),
            "gene_orthogroup_candidates_sha256": sha256(args.output),
            "dataset_summary": str(args.dataset_summary),
            "dataset_summary_sha256": sha256(args.dataset_summary),
            "orthogroups_sha256": sha256(assigned_path),
            "unassigned_sha256": sha256(unassigned_path),
            "species_tree": tree_path.read_text(encoding="utf-8").strip(),
            "species_tree_sha256": sha256(tree_path),
            "run_log_sha256": sha256(args.run_log),
        },
        "limitations": [
            "Zea mays is excluded pending a source-matched Zm00001d.1 protein FASTA.",
            "All gene, protein, and orthogroup mappings remain candidates pending acceptance review.",
            "Orthogroup coverage does not by itself validate cell labels or model performance.",
        ],
        "safety": {"source_h5ad_modified": False, "external_private_data_sent": False, "partial_outputs_remaining": False},
        "next_gate": "freeze curator-reviewed root, leaf, and vascular label panels and leakage-proof splits",
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    audit_tmp = args.audit.with_suffix(args.audit.suffix + ".partial")
    audit_tmp.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(audit_tmp, args.audit)
    print(json.dumps(audit["results"], ensure_ascii=False))


if __name__ == "__main__":
    main()
