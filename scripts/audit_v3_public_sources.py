#!/usr/bin/env python3

import csv
import json
import os
from pathlib import Path

import requests


PROJECT = Path(".")
REPORT = PROJECT / "reports" / "external_validation_v3_source_audit_v1"
REPORT.mkdir(parents=True, exist_ok=False)

session = requests.Session()
session.headers.update({"User-Agent": "PhyloOpenCell-source-audit/1.0"})

gitlab_project = "hhu-plant-biochemistry%2Ftriesch2025_moricandia_snrna_seq"
tree_url = f"https://git.nfdi4plants.org/api/v4/projects/{gitlab_project}/repository/tree"
entries = []
page = 1
while True:
    response = session.get(
        tree_url,
        params={"recursive": "true", "per_page": 100, "page": page},
        timeout=60,
    )
    response.raise_for_status()
    batch = response.json()
    if not batch:
        break
    entries.extend(batch)
    next_page = response.headers.get("X-Next-Page")
    if not next_page:
        break
    page = int(next_page)

interesting_suffixes = (
    ".h5ad", ".h5", ".mtx", ".mtx.gz", ".csv", ".csv.gz",
    ".tsv", ".tsv.gz", ".rds", ".qs", ".loom", ".zip", ".tar.gz",
)
interesting = []
for entry in entries:
    if entry.get("type") != "blob":
        continue
    path = entry["path"]
    lower = path.lower()
    if (
        "assays/02_final_experiment/dataset" in lower
        or "metadata" in lower
        or lower.endswith(interesting_suffixes)
        or lower.endswith(".arc")
    ):
        blob = session.get(
            f"https://git.nfdi4plants.org/api/v4/projects/{gitlab_project}/repository/blobs/{entry['id']}",
            timeout=60,
        )
        blob.raise_for_status()
        info = blob.json()
        interesting.append(
            {
                "path": path,
                "mode": entry.get("mode"),
                "blob_id": entry["id"],
                "repository_blob_size": info.get("size"),
            }
        )

(REPORT / "moricandia_repository_tree_interesting.json").write_text(
    json.dumps(interesting, indent=2), encoding="utf-8"
)

runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
runinfo_response = session.get(runinfo_url, params={"acc": "PRJNA865791"}, timeout=60)
runinfo_response.raise_for_status()
(REPORT / "PRJNA865791_runinfo.csv").write_bytes(runinfo_response.content)
run_rows = list(csv.DictReader(runinfo_response.text.splitlines()))

size_fields = ["size_MB", "Size_MB", "size", "bases"]
summary = {
    "moricandia": {
        "accession": "PRJNA1186371",
        "repository_entries": len(entries),
        "interesting_blobs": len(interesting),
        "candidate_files": interesting,
    },
    "maize": {
        "accession": "PRJNA865791",
        "sra_runs": len(run_rows),
        "run_accessions": [row.get("Run") for row in run_rows],
        "sample_names": [row.get("SampleName") or row.get("LibraryName") for row in run_rows],
        "reported_size_fields": {
            field: [row.get(field) for row in run_rows] for field in size_fields if field in (run_rows[0] if run_rows else {})
        },
    },
    "audit": {
        "candidate_expression_inspected": False,
        "files_downloaded": False,
        "raw_data_modified": False,
        "https_proxy_present": bool(os.environ.get("HTTPS_PROXY")),
    },
}
(REPORT / "source_audit_summary.json").write_text(
    json.dumps(summary, indent=2), encoding="utf-8"
)
print(json.dumps(summary, indent=2))
