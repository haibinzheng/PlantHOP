#!/usr/bin/env python3
"""Download a public UniProt proteome and compare its gene names locally."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import pandas as pd
import requests


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--proteome", required=True)
    parser.add_argument("--reference-dir", required=True)
    args = parser.parse_args()
    project = Path(args.project)
    ref_dir = Path(args.reference_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": f"(proteome:{args.proteome})",
        "format": "tsv",
        "fields": "accession,gene_names,organism_id",
    }
    reference = ref_dir / f"{args.proteome}.uniprot_gene_names.tsv.gz"
    if not reference.exists():
        response = requests.get(url, params=params, timeout=180, headers={"User-Agent": "PhyloOpenCell-local-reference-audit/0.1"})
        response.raise_for_status()
        partial = reference.with_suffix(reference.suffix + ".partial")
        with gzip.open(partial, "wb") as handle:
            handle.write(response.content)
        partial.replace(reference)
        resolved_url = response.url
    else:
        resolved_url = requests.Request("GET", url, params=params).prepare().url

    table = pd.read_csv(reference, sep="\t", compression="gzip", keep_default_na=False)
    gene_col = next(column for column in table.columns if column.casefold().startswith("gene names"))
    ids = set()
    for value in table[gene_col]:
        ids.update(str(value).split())

    manifest = pd.read_csv(project / "metadata" / "dataset_manifest.tsv", sep="\t", keep_default_na=False)
    rows = []
    for _, dataset in manifest[manifest.proposed_species.eq(args.species)].iterrows():
        adata = ad.read_h5ad(dataset.h5ad_path, backed="r")
        try:
            values = [str(x) for x in adata.var_names]
        finally:
            adata.file.close()
        matches = sum(value in ids for value in values)
        transform_matches = sum(
            value not in ids and (value.replace("-", "_") in ids or value.replace("_", "-") in ids)
            for value in values
        )
        rows.append({
            "dataset_id": dataset.dataset_id,
            "species": args.species,
            "n_features": len(values),
            "exact_gene_name_matches": matches,
            "exact_gene_name_fraction": matches / len(values) if values else 0,
            "punctuation_transform_candidates": transform_matches,
            "punctuation_transform_fraction": transform_matches / len(values) if values else 0,
            "uniprot_proteome": args.proteome,
            "reference_path": str(reference),
            "reference_sha256": sha256(reference),
            "source_url": resolved_url,
            "mapping_semantics": "exact match to a UniProt gene-name token; gene-to-protein cardinality review required",
            "accepted": False,
        })
    output = pd.DataFrame(rows)
    safe_species = args.species.replace(" ", "_")
    output.to_csv(project / "metadata" / f"uniprot_gene_name_coverage_{safe_species}.tsv", sep="\t", index=False)
    summary = {
        "species": args.species,
        "proteome": args.proteome,
        "proteins": len(table),
        "unique_gene_name_tokens": len(ids),
        "datasets": len(output),
        "mean_exact_gene_name_fraction": float(output.exact_gene_name_fraction.mean()),
        "minimum_exact_gene_name_fraction": float(output.exact_gene_name_fraction.min()),
        "datasets_ge_0_9": int(output.exact_gene_name_fraction.ge(0.9).sum()),
        "mean_exact_plus_transform_fraction": float(
            (output.exact_gene_name_fraction + output.punctuation_transform_fraction).mean()
        ),
        "minimum_exact_plus_transform_fraction": float(
            (output.exact_gene_name_fraction + output.punctuation_transform_fraction).min()
        ),
        "private_identifiers_sent_externally": False,
        "accepted": False,
    }
    (project / "reports" / f"uniprot_coverage_{safe_species}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
