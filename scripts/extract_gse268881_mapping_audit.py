#!/usr/bin/env python3
"""Extract unique mapping and marker-overlap evidence from the blind audit JSON."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def main(audit_path: Path, mapping_path: Path, overlap_path: Path) -> None:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    frozen = audit["frozen_orthogroups"]
    gene_to_frozen = {gene: og for og, genes in frozen.items() for gene in genes}

    mapping = {}
    marker_genes = set()
    for hit in audit["hits"]:
        vals = hit["values"]
        if hit["workbook"] == "Supplementary_Data_8.xlsx" and hit["sheet"] == "ABA":
            for gene in hit["target_genes"]:
                mapping.setdefault(gene, {
                    "frozen_og": gene_to_frozen[gene],
                    "ath_gene": gene,
                    "author_og": vals[0],
                    "esa_gene": vals[2],
                    "source_workbook": hit["workbook"],
                    "source_sheet": hit["sheet"],
                })
        if hit["workbook"] == "Supplementary_Data_2.xlsx" and hit["sheet"] == "TableS2":
            marker_genes.update(hit["target_genes"])

    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "frozen_og", "ath_gene", "author_og", "esa_gene",
            "source_workbook", "source_sheet", "in_author_marker_table",
        ], delimiter="\t")
        writer.writeheader()
        for gene in sorted(gene_to_frozen, key=lambda g: (gene_to_frozen[g], g)):
            row = mapping.get(gene, {
                "frozen_og": gene_to_frozen[gene], "ath_gene": gene,
                "author_og": "", "esa_gene": "", "source_workbook": "", "source_sheet": "",
            })
            row["in_author_marker_table"] = str(gene in marker_genes).lower()
            writer.writerow(row)

    with overlap_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["frozen_og", "ath_gene", "in_author_marker_table"])
        for gene in sorted(gene_to_frozen, key=lambda g: (gene_to_frozen[g], g)):
            writer.writerow([gene_to_frozen[gene], gene, str(gene in marker_genes).lower()])

    mapped_ogs = {gene_to_frozen[g] for g, row in mapping.items() if row["esa_gene"]}
    overlap_ogs = {gene_to_frozen[g] for g in marker_genes}
    print(json.dumps({
        "candidate_genes": len(gene_to_frozen),
        "mapped_genes": len(mapping),
        "candidate_ogs": len(frozen),
        "mapped_ogs": len(mapped_ogs),
        "mapping_coverage": len(mapped_ogs) / len(frozen),
        "author_marker_overlap_genes": len(marker_genes),
        "author_marker_overlap_ogs": len(overlap_ogs),
    }))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
