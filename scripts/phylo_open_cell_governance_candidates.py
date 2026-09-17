#!/usr/bin/env python3
"""Build reviewable governance candidates; never modifies source H5AD files."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


def norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def parse_obo(path: Path) -> list[dict]:
    terms: list[dict] = []
    current: dict | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "[Term]":
            if current and current.get("id") and current.get("name") and not current.get("obsolete"):
                terms.append(current)
            current = {"synonyms": [], "obsolete": False}
        elif line.startswith("["):
            if current and current.get("id") and current.get("name") and not current.get("obsolete"):
                terms.append(current)
            current = None
        elif current is not None and line.startswith("id: "):
            current["id"] = line[4:]
        elif current is not None and line.startswith("name: "):
            current["name"] = line[6:]
        elif current is not None and line.startswith("synonym: "):
            match = re.match(r'synonym: "(.*?)"', line)
            if match:
                current["synonyms"].append(match.group(1))
        elif current is not None and line == "is_obsolete: true":
            current["obsolete"] = True
    if current and current.get("id") and current.get("name") and not current.get("obsolete"):
        terms.append(current)
    return terms


def ontology_candidates(vocabulary: pd.DataFrame, obo_path: Path) -> pd.DataFrame:
    exact_names: dict[str, list[dict]] = defaultdict(list)
    exact_synonyms: dict[str, list[dict]] = defaultdict(list)
    for term in parse_obo(obo_path):
        exact_names[norm(term["name"])].append(term)
        for synonym in term["synonyms"]:
            exact_synonyms[norm(synonym)].append(term)

    rows = []
    for _, source in vocabulary.iterrows():
        key = norm(source["source_label"])
        candidates = exact_names.get(key, [])
        relation = "exact_name"
        if not candidates:
            candidates = exact_synonyms.get(key, [])
            relation = "exact_synonym"
        if bool(source.get("is_unknown_label", False)) or not key:
            status = "excluded_unknown"
            candidates = []
            relation = ""
        elif len(candidates) == 1:
            status = "candidate_auto_exact"
        elif len(candidates) > 1:
            status = "candidate_ambiguous"
        else:
            status = "pending_manual_curation"
        rows.append({
            "source_label": source["source_label"],
            "normalized_label": key,
            "cell_count": source["cell_count"],
            "dataset_count": source["dataset_count"],
            "species_count": source["species_count"],
            "candidate_status": status,
            "mapping_relation": relation,
            "candidate_ids": "|".join(x["id"] for x in candidates),
            "candidate_names": "|".join(x["name"] for x in candidates),
            "evidence": "Plant Ontology exact lexical match; requires curator acceptance" if candidates else "",
            "accepted": False,
        })
    return pd.DataFrame(rows)


def source_reconciliation(source: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, item in source.iterrows():
        proposed = str(item.get("proposed_species", ""))
        official = str(item.get("organism", ""))
        status = str(item.get("status", ""))
        relation = "exact_or_contained" if bool(item.get("species_matches_official", False)) else ""
        action = "retain_proposed"
        note = ""
        if status != "found":
            relation = "source_adapter_pending"
            action = "hold_for_source_review"
            note = "Repository page recorded, but structured metadata adapter is pending."
        elif proposed == "Hylocereus undatus" and official == "Selenicereus undatus":
            relation = "taxonomy_synonym_candidate"
            action = "review_taxonomic_synonym"
            note = "Do not rename until an authoritative taxonomic source is pinned."
        elif " x " in proposed.casefold() and set(norm(proposed).split()) == set(norm(official).split()):
            relation = "hybrid_parent_order"
            action = "retain_source_order_and_add_normalized_hybrid"
            note = "Parent order differs; preserve both source and normalized names."
        elif official.endswith(" sp.") and proposed.startswith(official[:-4]):
            relation = "official_broad_taxon"
            action = "retain_proposed_with_warning"
            note = "Official project metadata is less specific than the dataset label."
        elif not bool(item.get("species_matches_official", False)):
            relation = "mismatch_requires_review"
            action = "hold_for_manual_review"
        rows.append({
            "dataset_id": item.get("dataset_id", ""),
            "query_accession": item.get("query_accession", ""),
            "proposed_species": proposed,
            "official_organism": official,
            "source_status": status,
            "repository": item.get("repository", ""),
            "source_url": item.get("url", ""),
            "reconciliation_relation": relation,
            "proposed_action": action,
            "review_note": note,
            "accepted": False,
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    args = parser.parse_args()
    project = Path(args.project)
    metadata = project / "metadata"
    reports = project / "reports"

    vocabulary = pd.read_csv(metadata / "cell_type_vocabulary.tsv", sep="\t", keep_default_na=False)
    ontology = ontology_candidates(vocabulary, project / "resources" / "ontology" / "po.obo")
    ontology.to_csv(metadata / "ontology_mapping_candidates.tsv", sep="\t", index=False)

    source = pd.read_csv(metadata / "dataset_source_metadata.tsv", sep="\t", keep_default_na=False)
    reconciled = source_reconciliation(source)
    reconciled.to_csv(metadata / "source_reconciliation_candidates.tsv", sep="\t", index=False)

    summary = {
        "ontology_labels": len(ontology),
        "ontology_auto_exact_candidates": int(ontology.candidate_status.eq("candidate_auto_exact").sum()),
        "ontology_ambiguous_candidates": int(ontology.candidate_status.eq("candidate_ambiguous").sum()),
        "ontology_pending_manual": int(ontology.candidate_status.eq("pending_manual_curation").sum()),
        "ontology_excluded_unknown": int(ontology.candidate_status.eq("excluded_unknown").sum()),
        "source_dataset_rows": len(reconciled),
        "source_exact_or_contained": int(reconciled.reconciliation_relation.eq("exact_or_contained").sum()),
        "source_adapter_pending": int(reconciled.reconciliation_relation.eq("source_adapter_pending").sum()),
        "source_special_review": int((~reconciled.reconciliation_relation.isin(["exact_or_contained", "source_adapter_pending"])).sum()),
        "note": "All mappings are candidates. accepted=false is intentional.",
    }
    (reports / "governance_candidate_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
