#!/usr/bin/env python3
"""Fetch versioned public proteomes with atomic writes and provenance hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ENSEMBL_BASE = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/fasta"
SOURCES = [
    ("Arabidopsis thaliana", "ensembl_plants_63", "Arabidopsis_thaliana.TAIR10.pep.all.fa.gz", f"{ENSEMBL_BASE}/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"),
    ("Brassica rapa", "ensembl_plants_63", "Brassica_rapa.Brapa_1.0.pep.all.fa.gz", f"{ENSEMBL_BASE}/brassica_rapa/pep/Brassica_rapa.Brapa_1.0.pep.all.fa.gz"),
    ("Glycine max", "ensembl_plants_63", "Glycine_max.Glycine_max_v2.1.pep.all.fa.gz", f"{ENSEMBL_BASE}/glycine_max/pep/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz"),
    ("Oryza sativa", "ensembl_plants_63", "Oryza_sativa.IRGSP-1.0.pep.all.fa.gz", f"{ENSEMBL_BASE}/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz"),
    ("Populus trichocarpa", "ensembl_plants_63", "Populus_trichocarpa.Pop_tri_v4.pep.all.fa.gz", f"{ENSEMBL_BASE}/populus_trichocarpa/pep/Populus_trichocarpa.Pop_tri_v4.pep.all.fa.gz"),
    ("Sorghum bicolor", "ensembl_plants_63", "Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz", f"{ENSEMBL_BASE}/sorghum_bicolor/pep/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz"),
    ("Catharanthus roseus", "uniprot_UP001060085", "UP001060085.uniprot.fasta.gz", "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3AUP001060085%29"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, target: Path) -> str:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".partial")
    if partial.exists():
        partial.unlink()
    request = urllib.request.Request(url, headers={"User-Agent": "PhyloOpenCell/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as out:
            shutil.copyfileobj(response, out, length=1024 * 1024)
        if partial.stat().st_size == 0:
            raise RuntimeError(f"empty download: {url}")
        os.replace(partial, target)
    finally:
        if partial.exists():
            partial.unlink()
    return sha256(target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for species, source, filename, url in SOURCES:
        subdir = "uniprot" if source.startswith("uniprot") else "ensembl_plants_release_63"
        target = args.reference_root / subdir / filename
        status = "existing"
        if not target.exists():
            digest = download(url, target)
            status = "downloaded"
        else:
            digest = sha256(target)
        rows.append({
            "species": species,
            "source": source,
            "url": url,
            "local_path": str(target),
            "bytes": target.stat().st_size,
            "sha256": digest,
            "status": status,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        })
        print(f"{status}\t{species}\t{target.stat().st_size}\t{digest}", flush=True)
    args.registry.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.registry.with_suffix(args.registry.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, args.registry)


if __name__ == "__main__":
    main()
