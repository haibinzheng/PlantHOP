#!/usr/bin/env python3
"""Measure Zea mays Zm00001d coverage against the pinned MaizeGDB v4 list."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import pandas as pd


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--reference", required=True)
    args = parser.parse_args()
    project = Path(args.project)
    reference = Path(args.reference)
    manifest = pd.read_csv(project / "metadata" / "dataset_manifest.tsv", sep="\t", keep_default_na=False)
    with gzip.open(reference, "rt", encoding="utf-8") as handle:
        if reference.name.endswith(".gff3.gz"):
            ids = set()
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) == 9 and fields[2] == "gene":
                    attrs = dict(
                        item.split("=", 1) for item in fields[8].split(";") if "=" in item
                    )
                    if attrs.get("gene_id"):
                        ids.add(attrs["gene_id"])
        else:
            ids = {line.split()[0] for line in handle if line.strip() and not line.lstrip().startswith("gene_name")}
    rows = []
    for _, dataset in manifest[manifest.proposed_species.eq("Zea mays")].iterrows():
        adata = ad.read_h5ad(dataset.h5ad_path, backed="r")
        try:
            values = [str(x) for x in adata.var_names]
        finally:
            adata.file.close()
        matches = sum(x in ids for x in values)
        rows.append({
            "dataset_id": dataset.dataset_id,
            "species": "Zea mays",
            "n_features": len(values),
            "exact_matches": matches,
            "exact_fraction": matches / len(values) if values else 0,
            "reference_name": "Zm-B73-REFERENCE-GRAMENE-4.0 Zm00001d source-matched annotation",
            "reference_path": str(reference),
            "reference_sha256": file_sha256(reference),
            "source_url": "https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.genemodel_locus.txt.gz",
            "accepted": False,
        })
    out = pd.DataFrame(rows)
    out.to_csv(project / "metadata" / "maize_v4_exact_id_coverage.tsv", sep="\t", index=False)
    summary = {
        "datasets": len(out),
        "reference_gene_ids": len(ids),
        "mean_exact_fraction": float(out.exact_fraction.mean()),
        "minimum_exact_fraction": float(out.exact_fraction.min()),
        "datasets_ge_0_9": int(out.exact_fraction.ge(0.9).sum()),
        "private_identifiers_sent_externally": False,
        "accepted": False,
    }
    (project / "reports" / "maize_v4_coverage_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
