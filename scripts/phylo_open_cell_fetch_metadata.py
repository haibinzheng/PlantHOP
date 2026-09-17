#!/usr/bin/env python3
"""Fetch authoritative project metadata without modifying source H5AD files."""

from __future__ import annotations

import argparse
import json
import re
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests


ACCESSION_RE = re.compile(r"(PRJ(?:NA|EB|DB|CA)\d+|GSE\d+|E-?MTAB-?\d+|CRA\d+|CNP\d+)", re.I)


def accession_for(row: pd.Series) -> str:
    recorded = str(row.get("accession", "")).strip()
    if recorded and recorded.lower() != "not explicitly recorded":
        return recorded
    match = ACCESSION_RE.search(str(row["dataset_id"]))
    return match.group(1).upper() if match else ""


def get_json(session: requests.Session, url: str, params: dict | None = None, tries: int = 4) -> dict:
    for attempt in range(tries):
        response = session.get(url, params=params, timeout=45)
        if response.status_code == 200:
            return response.json()
        if response.status_code in {429, 500, 502, 503, 504}:
            time.sleep(2 ** attempt)
            continue
        response.raise_for_status()
    raise RuntimeError(f"Failed after {tries} attempts: {url}")


def fetch_ncbi_bioproject(session: requests.Session, accession: str) -> tuple[dict, dict]:
    search = get_json(session, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", {
        "db": "bioproject", "term": f"{accession}[Project Accession]", "retmode": "json"
    })
    ids = search.get("esearchresult", {}).get("idlist", [])
    time.sleep(0.36)
    if not ids:
        return {}, {"search": search}
    summary = get_json(session, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", {
        "db": "bioproject", "id": ids[0], "retmode": "json"
    })
    time.sleep(0.36)
    record = summary.get("result", {}).get(ids[0], {})
    return {
        "status": "found",
        "repository": "NCBI BioProject",
        "official_accession": record.get("project_acc", accession),
        "organism": record.get("organism_name", ""),
        "title": record.get("project_title", ""),
        "description": record.get("project_description", ""),
        "registration_date": record.get("registration_date", ""),
        "submitter": record.get("submitter_organization", ""),
        "url": f"https://www.ncbi.nlm.nih.gov/bioproject/{accession}",
    }, {"search": search, "summary": summary}


def fetch_geo(session: requests.Session, accession: str) -> tuple[dict, dict]:
    search = get_json(session, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", {
        "db": "gds", "term": f"{accession}[Accession]", "retmode": "json"
    })
    ids = search.get("esearchresult", {}).get("idlist", [])
    time.sleep(0.36)
    if not ids:
        return {}, {"search": search}
    summary = get_json(session, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", {
        "db": "gds", "id": ids[0], "retmode": "json"
    })
    time.sleep(0.36)
    record = summary.get("result", {}).get(ids[0], {})
    taxon = record.get("taxon", "")
    if isinstance(taxon, list):
        taxon = "|".join(str(item) for item in taxon)
    else:
        taxon = str(taxon)
    return {
        "status": "found",
        "repository": "NCBI GEO",
        "official_accession": accession,
        "organism": taxon,
        "title": record.get("title", ""),
        "description": record.get("summary", ""),
        "registration_date": record.get("pdat", ""),
        "submitter": "",
        "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
    }, {"search": search, "summary": summary}


def fetch_ngdc_bioproject(session: requests.Session, accession: str) -> tuple[dict, dict]:
    raw = get_json(session, f"https://ngdc.cncb.ac.cn/gwh/api/public/bioProject/{accession}")
    organism = raw.get("organism") or raw.get("organismName") or raw.get("species") or raw.get("taxons") or ""
    if isinstance(organism, list):
        organism = "|".join(
            str(item.get("name", "")) if isinstance(item, dict) else str(item)
            for item in organism
        )
    if isinstance(organism, dict):
        organism = organism.get("organismName") or organism.get("name") or ""
    data_types = raw.get("dataTypes", [])
    submitter = raw.get("organization", raw.get("submitterOrganization", raw.get("submitter", "")))
    if isinstance(submitter, dict):
        submitter = submitter.get("organization") or " ".join(
            x for x in [submitter.get("firstName", ""), submitter.get("lastName", "")] if x
        )
    return {
        "status": "found",
        "repository": "CNCB-NGDC BioProject",
        "official_accession": raw.get("accession", raw.get("prjAccession", accession)),
        "organism": str(organism),
        "title": raw.get("title", raw.get("projectTitle", "")),
        "description": raw.get("description", ""),
        "registration_date": raw.get("releaseDate", raw.get("releaseTime", raw.get("submissionDate", ""))),
        "submitter": submitter,
        "url": f"https://ngdc.cncb.ac.cn/bioproject/browse/{accession}",
        "data_types": "|".join(str(x.get("dataTypeName", "")) for x in data_types if isinstance(x, dict)),
    }, raw


def fetch_biostudies(session: requests.Session, accession: str) -> tuple[dict, dict]:
    normalized = accession.replace("EMTAB", "E-MTAB-").replace("E--", "E-")
    raw = get_json(session, f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{normalized}")
    attrs = {}
    for item in raw.get("section", {}).get("attributes", []):
        attrs[str(item.get("name", "")).lower()] = str(item.get("value", ""))
    return {
        "status": "found",
        "repository": "EMBL-EBI BioStudies",
        "official_accession": normalized,
        "organism": attrs.get("organism", ""),
        "title": attrs.get("title", attrs.get("study title", "")),
        "description": attrs.get("description", attrs.get("study description", "")),
        "registration_date": attrs.get("release date", ""),
        "submitter": "",
        "url": f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{normalized}",
    }, raw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--refresh-pattern",
        default="",
        help="Regex of accessions to refetch; other rows are reused from an existing registry.",
    )
    args = parser.parse_args()
    output = Path(args.output_dir)
    cache = output / "source_cache"
    cache.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    manifest["query_accession"] = manifest.apply(accession_for, axis=1)
    datasets = defaultdict(list)
    for _, row in manifest.iterrows():
        datasets[row["query_accession"]].append(row["dataset_id"])

    registry_path = output / "source_project_registry.tsv"
    existing = {}
    if registry_path.exists():
        previous = pd.read_csv(registry_path, sep="\t", keep_default_na=False)
        existing = {str(row["query_accession"]): row.to_dict() for _, row in previous.iterrows()}
    refresh_re = re.compile(args.refresh_pattern, re.I) if args.refresh_pattern else None

    session = requests.Session()
    session.headers.update({"User-Agent": "PhyloOpenCell-governance/0.1 (research metadata audit)"})
    rows = []
    for accession in sorted(x for x in datasets if x):
        if refresh_re is not None and not refresh_re.search(accession) and accession in existing:
            old = existing[accession]
            old["dataset_count"] = len(datasets[accession])
            old["dataset_ids"] = "|".join(datasets[accession])
            rows.append(old)
            print(accession, "reused", flush=True)
            continue
        cache_path = cache / f"{accession}.json"
        record = {}
        raw = {}
        error = ""
        try:
            if accession.upper().startswith("GSE"):
                record, raw = fetch_geo(session, accession)
            elif "MTAB" in accession.upper():
                record, raw = fetch_biostudies(session, accession)
            elif accession.upper().startswith("PRJCA"):
                record, raw = fetch_ngdc_bioproject(session, accession)
            elif accession.upper().startswith(("CNP", "CRA")):
                record = {
                    "status": "pending_repository_adapter",
                    "repository": "CNGBdb" if accession.upper().startswith("CNP") else "CNCB-NGDC GSA",
                    "official_accession": accession,
                    "url": (
                        f"https://db.cngb.org/data_resources/project/{accession}/"
                        if accession.upper().startswith("CNP")
                        else f"https://ngdc.cncb.ac.cn/gsa/browse/{accession}"
                    ),
                }
            else:
                record, raw = fetch_ncbi_bioproject(session, accession)
            if not record:
                record = {"status": "not_found", "repository": "NCBI BioProject"}
        except Exception as exc:
            record = {"status": "error", "repository": ""}
            error = f"{type(exc).__name__}: {exc}"
        cache_path.write_text(json.dumps(raw, ensure_ascii=False, indent=2), encoding="utf-8")
        rows.append({
            "query_accession": accession,
            "dataset_count": len(datasets[accession]),
            "dataset_ids": "|".join(datasets[accession]),
            "status": record.get("status", ""),
            "repository": record.get("repository", ""),
            "official_accession": record.get("official_accession", ""),
            "organism": record.get("organism", ""),
            "title": record.get("title", ""),
            "description": record.get("description", ""),
            "registration_date": record.get("registration_date", ""),
            "submitter": record.get("submitter", ""),
            "url": record.get("url", ""),
            "error": error,
        })
        print(accession, rows[-1]["status"], flush=True)

    registry = pd.DataFrame(rows)
    registry.to_csv(output / "source_project_registry.tsv", sep="\t", index=False)
    joined = manifest.merge(registry, on="query_accession", how="left", suffixes=("", "_official"))
    joined["species_matches_official"] = joined.apply(
        lambda r: str(r["proposed_species"]).lower() in str(r["organism"]).lower()
        or str(r["organism"]).lower() in str(r["proposed_species"]).lower(), axis=1
    )
    joined.to_csv(output / "dataset_source_metadata.tsv", sep="\t", index=False)
    summary = {
        "accessions": len(registry),
        "found": int((registry.status == "found").sum()),
        "not_found": int((registry.status == "not_found").sum()),
        "errors": int((registry.status == "error").sum()),
        "datasets_with_official_source": int(joined["status"].eq("found").sum()),
        "species_match": int((joined["status"].eq("found") & joined["species_matches_official"]).sum()),
        "species_mismatch_or_complex": int((joined["status"].eq("found") & ~joined["species_matches_official"]).sum()),
    }
    (output.parent / "reports" / "source_metadata_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
