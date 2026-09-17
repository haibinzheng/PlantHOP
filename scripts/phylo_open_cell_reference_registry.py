#!/usr/bin/env python3
"""Build reference-resource candidates from official APIs without selecting them."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests


QUERY_OVERRIDES = {
    "Manihot esculenta Crantz": "Manihot esculenta",
    "Hylocereus undatus": "Selenicereus undatus",
}


def clean_name(value: object) -> str:
    text = re.sub(r"\s*\([^)]*\)\s*$", "", str(value)).casefold()
    return re.sub(r"[^a-z0-9]+", " ", text).strip()


def get_json(session: requests.Session, url: str, params: dict | None = None, tries: int = 5) -> tuple[dict, requests.Response]:
    last: Exception | None = None
    for attempt in range(tries):
        try:
            response = session.get(url, params=params, timeout=90)
            if response.status_code == 200:
                return response.json(), response
            if response.status_code not in {429, 500, 502, 503, 504}:
                response.raise_for_status()
        except Exception as exc:
            last = exc
        time.sleep(min(20, 2 ** attempt))
    raise RuntimeError(f"official API failed: {url}: {last}")


def busco(record: dict) -> tuple[object, object]:
    report = record.get("proteomeCompletenessReport", {}).get("buscoReport", {})
    return report.get("score", ""), report.get("lineageDb", "")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-summary", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir)
    cache = output / "reference_cache"
    cache.mkdir(parents=True, exist_ok=True)
    species = pd.read_csv(args.species_summary, sep="\t", keep_default_na=False)["species"].tolist()

    session = requests.Session()
    session.headers.update({"User-Agent": "PhyloOpenCell-reference-audit/0.1"})
    retrieved = dt.datetime.now(dt.timezone.utc).isoformat()
    ensembl_raw, ensembl_response = get_json(
        session,
        "https://rest.ensembl.org/info/species",
        {"division": "EnsemblPlants", "content-type": "application/json"},
    )
    (cache / "ensembl_species.json").write_text(json.dumps(ensembl_raw, indent=2), encoding="utf-8")
    ensembl_records = ensembl_raw.get("species", [])

    rows = []
    for source_species in species:
        query_species = QUERY_OVERRIDES.get(source_species, source_species)
        query = f'(organism_name:"{query_species}")'
        error = ""
        raw = {"results": []}
        response_url = ""
        try:
            raw, response = get_json(
                session,
                "https://rest.uniprot.org/proteomes/search",
                {"query": query, "format": "json", "size": 100},
            )
            response_url = response.url
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
        (cache / f"uniprot_{re.sub(r'[^A-Za-z0-9]+', '_', source_species)}.json").write_text(
            json.dumps(raw, indent=2, ensure_ascii=False), encoding="utf-8"
        )

        exact = [
            item for item in raw.get("results", [])
            if clean_name(item.get("taxonomy", {}).get("scientificName", "")) == clean_name(query_species)
        ]
        used_legacy_taxon_name = False
        if not exact and query_species != source_species:
            legacy_raw, legacy_response = get_json(
                session,
                "https://rest.uniprot.org/proteomes/search",
                {"query": f'(organism_name:"{source_species}")', "format": "json", "size": 100},
            )
            exact = [
                item for item in legacy_raw.get("results", [])
                if clean_name(item.get("taxonomy", {}).get("scientificName", "")) == clean_name(source_species)
            ]
            if exact:
                used_legacy_taxon_name = True
                response_url = legacy_response.url
        exact.sort(key=lambda item: (
            item.get("proteomeType") == "Reference proteome",
            float(busco(item)[0] or -1),
            float(item.get("annotationScore") or -1),
            str(item.get("modified", "")),
        ), reverse=True)
        ens = []
        qkey = clean_name(query_species)
        for item in ensembl_records:
            names = [item.get("display_name", ""), item.get("name", "")] + item.get("aliases", [])
            internal_name = str(item.get("name", ""))
            if qkey in {clean_name(x) for x in names} or internal_name.startswith(qkey.replace(" ", "_") + "_gca"):
                ens.append(item)

        if " x " in source_species.casefold() or "euramericana" in source_species.casefold():
            status = "hybrid_requires_parent_policy"
        elif not exact:
            status = "no_exact_uniprot_candidate"
        elif used_legacy_taxon_name:
            status = "candidate_legacy_taxon_name"
        elif source_species != query_species:
            status = "candidate_via_taxonomy_override"
        else:
            status = "candidate_available_not_selected"

        chosen = exact[0] if exact else {}
        assembly = chosen.get("genomeAssembly", {})
        score, lineage = busco(chosen)
        rows.append({
            "source_species": source_species,
            "query_species": query_species,
            "candidate_status": status,
            "candidate_count_exact_name": len(exact),
            "uniprot_candidate_id": chosen.get("id", ""),
            "uniprot_candidate_organism": chosen.get("taxonomy", {}).get("scientificName", ""),
            "uniprot_taxon_id": chosen.get("taxonomy", {}).get("taxonId", ""),
            "uniprot_proteome_type": chosen.get("proteomeType", ""),
            "uniprot_modified": chosen.get("modified", ""),
            "uniprot_gene_count": chosen.get("geneCount", ""),
            "uniprot_protein_count": chosen.get("proteinCount", ""),
            "uniprot_annotation_score": chosen.get("annotationScore", ""),
            "uniprot_busco_score": score,
            "uniprot_busco_lineage": lineage,
            "assembly_accession": assembly.get("assemblyId", ""),
            "assembly_source": assembly.get("source", ""),
            "assembly_url": assembly.get("genomeAssemblyUrl", ""),
            "ensembl_match_count": len(ens),
            "ensembl_name": "|".join(str(x.get("name", "")) for x in ens),
            "ensembl_assembly": "|".join(str(x.get("assembly", "")) for x in ens),
            "ensembl_accession": "|".join(str(x.get("accession", "")) for x in ens),
            "ensembl_release": "|".join(str(x.get("release", "")) for x in ens),
            "ensembl_division": "|".join(str(x.get("division", "")) for x in ens),
            "uniprot_query_url": response_url,
            "retrieved_at": retrieved,
            "accepted": False,
            "review_note": error,
        })
        print(source_species, status, chosen.get("id", ""), flush=True)
        time.sleep(0.15)

    registry = pd.DataFrame(rows)
    registry.to_csv(output / "species_reference_candidates.tsv", sep="\t", index=False)
    summary = {
        "species": len(registry),
        "uniprot_exact_candidate_available": int(registry.uniprot_candidate_id.ne("").sum()),
        "uniprot_reference_candidate": int(registry.uniprot_proteome_type.eq("Reference proteome").sum()),
        "ensembl_match_available": int(registry.ensembl_match_count.gt(0).sum()),
        "no_exact_uniprot_candidate": int(registry.candidate_status.eq("no_exact_uniprot_candidate").sum()),
        "selected": 0,
        "source_endpoints": ["https://rest.uniprot.org/proteomes/search", "https://rest.ensembl.org/info/species"],
        "ensembl_payload_sha256": hashlib.sha256(json.dumps(ensembl_raw, sort_keys=True).encode()).hexdigest(),
        "retrieved_at": retrieved,
        "note": "Candidates are not accepted references; assembly and gene-ID compatibility require review.",
    }
    (output.parent / "reports" / "reference_registry_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
