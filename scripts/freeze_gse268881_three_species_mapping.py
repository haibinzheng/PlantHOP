#!/usr/bin/env python3
"""Freeze expression-blind Ath-to-three-species mappings for GSE268881."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def main(audit_path: Path, mapping_path: Path, summary_path: Path) -> None:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    frozen = audit["frozen_orthogroups"]
    rank_by_og = {og: rank for rank, og in enumerate(frozen, start=1)}
    gene_to_og = {gene: og for og, genes in frozen.items() for gene in genes}

    author_markers = set()
    mappings = {}
    for hit in audit["hits"]:
        values = hit["values"]
        if hit["workbook"] == "Supplementary_Data_2.xlsx" and hit["sheet"] == "TableS2":
            author_markers.update(hit["target_genes"])
        if hit["workbook"] == "Supplementary_Data_8.xlsx" and hit["sheet"] == "ABA":
            for gene in hit["target_genes"]:
                mappings.setdefault(gene, {
                    "author_orthogroup": values[0],
                    "ath_gene": values[1],
                    "esa_gene": values[2],
                    "sir_gene": values[3],
                    "spa_gene": values[4],
                })

    representative_by_og = {}
    for og, genes in frozen.items():
        mapped = sorted(gene for gene in genes if gene in mappings)
        representative_by_og[og] = mapped[0] if mapped else ""

    rows = []
    for og, genes in frozen.items():
        for gene in genes:
            source = mappings.get(gene, {})
            rows.append({
                "frozen_rank": rank_by_og[og],
                "project_orthogroup": og,
                "ath_gene": gene,
                "author_orthogroup": source.get("author_orthogroup", ""),
                "esa_gene": source.get("esa_gene", ""),
                "sir_gene": source.get("sir_gene", ""),
                "spa_gene": source.get("spa_gene", ""),
                "mapped_in_author_one_to_one_table": str(gene in mappings).lower(),
                "single_gene_representative": str(gene == representative_by_og[og]).lower(),
                "in_author_celltype_marker_table": str(gene in author_markers).lower(),
                "mapping_source": "Supplementary_Data_8:ABA:first_five_identifier_columns" if gene in mappings else "",
                "marker_audit_source": "Supplementary_Data_2:TableS2",
            })

    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    mapped_ogs = {
        og for og, genes in frozen.items()
        if any(gene in mappings for gene in genes)
    }
    overlap_ogs = {
        og for og, genes in frozen.items()
        if any(gene in author_markers for gene in genes)
    }
    sensitivity_ogs = set(frozen) - overlap_ogs
    sensitivity_mapped = sensitivity_ogs & mapped_ogs
    summary = {
        "status": "expression_blind_mapping_freeze_complete",
        "candidate_orthogroups": len(frozen),
        "candidate_genes": len(gene_to_og),
        "mapped_candidate_orthogroups": len(mapped_ogs),
        "mapped_candidate_genes": len(mappings),
        "mapping_coverage": len(mapped_ogs) / len(frozen),
        "unmapped_orthogroups": sorted(set(frozen) - mapped_ogs),
        "author_marker_overlap_genes": sorted(author_markers),
        "author_marker_overlap_orthogroups": sorted(overlap_ogs),
        "anti_circularity_candidate_orthogroups": len(sensitivity_ogs),
        "anti_circularity_mapped_orthogroups": len(sensitivity_mapped),
        "anti_circularity_mapping_coverage": len(sensitivity_mapped) / len(sensitivity_ogs),
        "single_gene_rule": "lexicographically_smallest_mapped_ath_gene_per_project_orthogroup",
        "expression_values_read": False,
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
