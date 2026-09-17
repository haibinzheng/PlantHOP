#!/usr/bin/env python3
"""Download public Ensembl Plants GTFs and compare H5AD gene IDs locally."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urljoin

import anndata as ad
import pandas as pd
import requests


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(session: requests.Session, url: str, destination: Path) -> None:
    partial = destination.with_suffix(destination.suffix + ".partial")
    with session.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with partial.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    partial.replace(destination)


def gtf_gene_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    pattern = re.compile(r'gene_id "([^"]+)"')
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            match = pattern.search(line)
            if match:
                ids.add(match.group(1))
    return ids


def transformed(value: str) -> set[str]:
    candidates = {value}
    if value.startswith("gene:"):
        candidates.add(value[5:])
    poplar = re.match(r"^POPTR-(\d+G\d+)v4$", value, flags=re.I)
    if poplar:
        candidates.add(f"Potri.{poplar.group(1)}.v4.1")
    for _ in range(2):
        for item in list(candidates):
            candidates.add(re.sub(r"\.\d+$", "", item))
            candidates.add(re.sub(r"-T\d+$", "", item))
            candidates.add(re.sub(r"v\d+$", "", item))
            candidates.add(item.replace("-", "."))
            candidates.add(item.replace("-", "_"))
    return candidates


def choose_ensembl_name(row: pd.Series) -> str:
    names = [x for x in str(row.get("ensembl_name", "")).split("|") if x]
    target = re.sub(r"[^a-z0-9]+", "_", str(row["query_species"]).casefold()).strip("_")
    if target in names:
        return target
    return names[0] if len(names) == 1 else ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--reference-registry", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--release", type=int, required=True)
    parser.add_argument("--species", required=True, help="Comma-separated source species")
    args = parser.parse_args()

    wanted = {x.strip() for x in args.species.split(",") if x.strip()}
    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    registry = pd.read_csv(args.reference_registry, sep="\t", keep_default_na=False)
    ref_dir = Path(args.reference_dir)
    out_dir = Path(args.output_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "PhyloOpenCell-local-reference-audit/0.1"})

    coverage_rows = []
    reference_rows = []
    for _, row in registry[registry.source_species.isin(wanted)].iterrows():
        species = row["source_species"]
        ensembl_name = choose_ensembl_name(row)
        if not ensembl_name:
            reference_rows.append({"species": species, "status": "ambiguous_ensembl_name"})
            continue
        directory_url = (
            f"https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-{args.release}/plants/gtf/{ensembl_name}/"
        )
        listing = session.get(directory_url, timeout=90)
        if listing.status_code != 200:
            reference_rows.append({"species": species, "ensembl_name": ensembl_name, "status": f"listing_http_{listing.status_code}", "source_url": directory_url})
            continue
        names = re.findall(r'href="([^"]+\.gtf\.gz)"', listing.text, flags=re.I)
        names = [
            x for x in names
            if "abinitio" not in x.casefold() and ".chr.gtf.gz" not in x.casefold()
        ]
        if len(names) != 1:
            reference_rows.append({"species": species, "ensembl_name": ensembl_name, "status": f"gtf_candidates_{len(names)}", "source_url": directory_url})
            continue
        filename = names[0]
        url = urljoin(directory_url, filename)
        destination = ref_dir / filename
        if not destination.exists():
            print("download", species, filename, flush=True)
            download(session, url, destination)
        ids = gtf_gene_ids(destination)
        reference_rows.append({
            "species": species,
            "ensembl_name": ensembl_name,
            "status": "downloaded_and_parsed",
            "source_url": url,
            "local_path": str(destination),
            "file_bytes": destination.stat().st_size,
            "sha256": sha256(destination),
            "reference_gene_ids": len(ids),
            "release": args.release,
            "accepted": False,
        })

        for _, dataset in manifest[manifest.proposed_species.eq(species)].iterrows():
            adata = ad.read_h5ad(dataset.h5ad_path, backed="r")
            try:
                values = [str(x) for x in adata.var_names]
            finally:
                if getattr(adata, "file", None) is not None:
                    adata.file.close()
            exact = sum(value in ids for value in values)
            normalized = sum(value not in ids and bool(transformed(value) & ids) for value in values)
            coverage_rows.append({
                "dataset_id": dataset.dataset_id,
                "species": species,
                "n_features": len(values),
                "exact_reference_matches": exact,
                "exact_reference_fraction": exact / len(values) if values else 0,
                "transform_only_candidates": normalized,
                "transform_candidate_fraction": normalized / len(values) if values else 0,
                "reference_release": args.release,
                "reference_file": filename,
                "mapping_status": "exact_match_measured_locally",
                "note": "Transform-only matches are candidates, not accepted mappings.",
            })
            print("coverage", dataset.dataset_id, exact, len(values), flush=True)

    references = pd.DataFrame(reference_rows)
    coverage = pd.DataFrame(coverage_rows)
    references.to_csv(out_dir / "ensembl_reference_downloads.tsv", sep="\t", index=False)
    coverage.to_csv(out_dir / "ensembl_exact_id_coverage.tsv", sep="\t", index=False)
    summary = {
        "requested_species": len(wanted),
        "references_parsed": int(references.status.eq("downloaded_and_parsed").sum()) if "status" in references else 0,
        "datasets_measured": len(coverage),
        "datasets_exact_coverage_ge_0_6": int(coverage.exact_reference_fraction.ge(0.6).sum()) if len(coverage) else 0,
        "datasets_exact_coverage_ge_0_9": int(coverage.exact_reference_fraction.ge(0.9).sum()) if len(coverage) else 0,
        "private_identifiers_sent_externally": False,
        "comparison_location": "GPU server local filesystem",
        "release": args.release,
    }
    (out_dir.parent / "reports" / "ensembl_coverage_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
