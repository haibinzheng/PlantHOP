#!/usr/bin/env python3
"""Attach locally available FASTA-header annotations to conserved OG candidates."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def atomic_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        frame.to_csv(temp, sep="\t", index=False)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        temp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def fasta_headers(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def parse_ensembl(path: Path, species: str) -> dict[str, dict[str, str]]:
    result = {}
    for header in fasta_headers(path):
        gene = re.search(r"(?:^|\s)gene:([^\s]+)", header)
        if not gene:
            continue
        gene_id = gene.group(1)
        symbol = re.search(r"(?:^|\s)gene_symbol:([^\s]+)", header)
        description = re.search(r"(?:^|\s)description:(.*?)(?:\s\[Source:|$)", header)
        record = result.setdefault(gene_id, {"species": species, "gene_id": gene_id, "gene_symbol": "", "description": "", "source": str(path)})
        if symbol and not record["gene_symbol"]:
            record["gene_symbol"] = symbol.group(1)
        if description and not record["description"]:
            record["description"] = description.group(1).strip()
    return result


def parse_uniprot(path: Path, species: str) -> dict[str, dict[str, str]]:
    result = {}
    for header in fasta_headers(path):
        gene = re.search(r"(?:^|\s)GN=([^\s]+)", header)
        if not gene:
            continue
        gene_id = gene.group(1)
        desc = header.split(" ", 1)[1] if " " in header else ""
        desc = desc.split(" OS=", 1)[0].strip()
        result[gene_id] = {"species": species, "gene_id": gene_id, "gene_symbol": "", "description": desc, "source": str(path)}
    return result


def parse_gene_field(text: str):
    for block in str(text).split(" | "):
        if ": " not in block:
            continue
        species, genes = block.split(": ", 1)
        for gene in genes.split(","):
            gene = gene.strip()
            if gene and gene != "NA":
                yield species, gene


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-dir", type=Path, required=True)
    p.add_argument("--references-dir", type=Path, required=True)
    p.add_argument("--uniprot-dir", type=Path, required=True)
    p.add_argument("--version", default="v2")
    args = p.parse_args()
    candidate_path = args.project_dir / "metadata" / f"conserved_program_candidates_{args.version}.tsv"
    candidates = pd.read_csv(candidate_path, sep="\t")

    files = {
        "Arabidopsis thaliana": args.references_dir / "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz",
        "Brassica rapa": args.references_dir / "Brassica_rapa.Brapa_1.0.pep.all.fa.gz",
        "Glycine max": args.references_dir / "Glycine_max.Glycine_max_v2.1.pep.all.fa.gz",
        "Oryza sativa": args.references_dir / "Oryza_sativa.IRGSP-1.0.pep.all.fa.gz",
        "Sorghum bicolor": args.references_dir / "Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz",
    }
    annotations: dict[tuple[str, str], dict[str, str]] = {}
    source_counts = {}
    for species, path in files.items():
        parsed = parse_ensembl(path, species)
        source_counts[species] = len(parsed)
        annotations.update({(species, gene): row for gene, row in parsed.items()})
    cath_path = args.uniprot_dir / "UP001060085.uniprot.fasta.gz"
    cath = parse_uniprot(cath_path, "Catharanthus roseus")
    source_counts["Catharanthus roseus"] = len(cath)
    annotations.update({("Catharanthus roseus", gene): row for gene, row in cath.items()})

    long_rows = []
    for row in candidates.itertuples(index=False):
        for species, gene in parse_gene_field(row.representative_genes_by_species):
            ann = annotations.get((species, gene), {})
            long_rows.append({
                "panel": row.panel,
                "coarse_label": row.coarse_label,
                "orthogroup_id": row.orthogroup_id,
                "candidate_species_count": int(row.species_count),
                "species": species,
                "gene_id": gene,
                "gene_symbol": ann.get("gene_symbol", ""),
                "local_description": ann.get("description", ""),
                "annotation_found": bool(ann),
                "annotation_source": ann.get("source", ""),
            })
    long = pd.DataFrame(long_rows)

    summary_rows = []
    for row in candidates.itertuples(index=False):
        subset = long[long["orthogroup_id"].eq(row.orthogroup_id) & long["panel"].eq(row.panel) & long["coarse_label"].eq(row.coarse_label)]
        descriptions = sorted({x for x in subset["local_description"] if x})
        symbols = sorted({x for x in subset["gene_symbol"] if x})
        arabidopsis = subset[subset["species"].eq("Arabidopsis thaliana")]
        summary_rows.append({
            "panel": row.panel,
            "coarse_label": row.coarse_label,
            "orthogroup_id": row.orthogroup_id,
            "species_count": int(row.species_count),
            "median_within_species_rank": float(row.median_within_species_rank),
            "mean_rank_score_difference": float(row.mean_rank_score_difference),
            "species_with_multi_gene_family": int(row.species_with_multi_gene_family),
            "low_coarse_specificity_flag": bool(row.low_coarse_specificity_flag),
            "representative_genes_checked": len(subset),
            "genes_with_local_annotation": int(subset["annotation_found"].sum()),
            "annotation_coverage_fraction": float(subset["annotation_found"].mean()) if len(subset) else 0.0,
            "arabidopsis_genes": ";".join(arabidopsis["gene_id"]),
            "arabidopsis_symbols": ";".join(sorted({x for x in arabidopsis["gene_symbol"] if x})),
            "local_symbols_all_species": ";".join(symbols),
            "local_descriptions": " | ".join(descriptions),
            "priority_tier": "A" if int(row.species_count) >= 4 and not bool(row.low_coarse_specificity_flag) else "B",
            "paralog_interpretation": "family_expansion_hint_gene_level_test_required" if int(row.species_with_multi_gene_family) > 0 else "single_representative_per_species_in_checked_set",
        })
    summary = pd.DataFrame(summary_rows).sort_values(
        ["priority_tier", "species_count", "median_within_species_rank", "mean_rank_score_difference"],
        ascending=[True, False, True, False],
    )
    summary_path = args.project_dir / "metadata" / f"conserved_program_function_evidence_{args.version}.tsv"
    long_path = args.project_dir / "metadata" / f"conserved_program_gene_annotations_{args.version}.tsv"
    audit_path = args.project_dir / "reports" / f"conserved_program_function_evidence_{args.version}_audit.json"
    atomic_tsv(summary_path, summary)
    atomic_tsv(long_path, long)
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "local_header_annotation_complete",
        "candidate_rows": len(summary),
        "priority_A_rows": int(summary["priority_tier"].eq("A").sum()),
        "candidate_rows_with_any_description": int(summary["local_descriptions"].ne("").sum()),
        "representative_gene_rows": len(long),
        "representative_gene_rows_annotated": int(long["annotation_found"].sum()),
        "source_gene_counts": source_counts,
        "outputs": {"candidate_summary": str(summary_path), "gene_annotations": str(long_path)},
        "limitations": [
            "Annotations come only from local FASTA headers and are incomplete for several species.",
            "A matching description does not validate cell-type specificity or conserved function.",
            "Orthogroup-level rank features cannot identify which paralog drives expression; gene-level matrices are required.",
            "No gene identifiers or expression data were sent to external services.",
        ],
    }
    atomic_json(audit_path, audit)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
