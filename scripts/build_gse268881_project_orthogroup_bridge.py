#!/usr/bin/env python3
"""Bridge GSE268881 public one-to-one orthologs to the frozen project OG namespace."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(orthogroups: Path, public_table: Path, output: Path, audit_path: Path) -> None:
    if output.exists() or audit_path.exists():
        raise RuntimeError("refusing to overwrite bridge outputs")
    gene_ogs = defaultdict(set)
    with orthogroups.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        ath_col = "arabidopsis_thaliana.representative_proteins"
        for row in reader:
            for token in row[ath_col].split(","):
                token = token.strip()
                if token:
                    gene_ogs[token.split("|", 1)[0]].add(row["Orthogroup"])
    unique = {gene: next(iter(ogs)) for gene, ogs in gene_ogs.items() if len(ogs) == 1}

    rows = []
    with public_table.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            project_og = unique.get(row["ath_gene"], "")
            if project_og:
                rows.append({**row, "project_orthogroup": project_og})
    if not rows:
        raise RuntimeError("bridge is empty")
    fields = [
        "author_orthogroup", "ath_gene", "project_orthogroup",
        "esa_gene", "sir_gene", "spa_gene",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "status": "gse268881_project_orthogroup_bridge_complete",
        "orthogroups_source": str(orthogroups),
        "orthogroups_sha256": sha256(orthogroups),
        "public_table_source": str(public_table),
        "public_table_sha256": sha256(public_table),
        "public_rows": 15198,
        "bridged_rows": len(rows),
        "bridged_project_orthogroups": len({row["project_orthogroup"] for row in rows}),
        "ambiguous_arabidopsis_gene_ids": sum(len(ogs) > 1 for ogs in gene_ogs.values()),
        "expression_values_read": False,
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main(*(Path(arg) for arg in sys.argv[1:5]))
