#!/usr/bin/env python3
"""Turn measured local reference coverage into a review queue and phase report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def triage(row: pd.Series) -> str:
    exact = float(row["preferred_exact_fraction"])
    transformed = float(row["transform_candidate_fraction"])
    potential = exact + transformed
    if exact >= 0.90:
        return "direct_exact_candidate"
    if exact >= 0.60:
        return "partial_exact_review"
    if potential >= 0.90 and transformed >= 0.20:
        return "deterministic_transform_candidate"
    if potential >= 0.60:
        return "partial_transform_review"
    if potential < 0.10:
        return "reference_version_or_namespace_mismatch"
    return "low_coverage_manual_review"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    args = parser.parse_args()
    project = Path(args.project)
    metadata = project / "metadata"
    reports = project / "reports"
    coverage = pd.read_csv(metadata / "ensembl_exact_id_coverage.tsv", sep="\t")
    catharanthus_path = metadata / "uniprot_gene_name_coverage_Catharanthus_roseus.tsv"
    if catharanthus_path.exists():
        uniprot = pd.read_csv(catharanthus_path, sep="\t")
        extra = pd.DataFrame({
            "dataset_id": uniprot.dataset_id,
            "species": uniprot.species,
            "n_features": uniprot.n_features,
            "exact_reference_matches": uniprot.exact_gene_name_matches,
            "exact_reference_fraction": uniprot.exact_gene_name_fraction,
            "transform_only_candidates": uniprot.punctuation_transform_candidates,
            "transform_candidate_fraction": uniprot.punctuation_transform_fraction,
            "reference_release": uniprot.uniprot_proteome,
            "reference_file": uniprot.reference_path,
            "mapping_status": "uniprot_gene_name_match_measured_locally",
            "note": uniprot.mapping_semantics,
        })
        coverage = pd.concat([coverage, extra], ignore_index=True)
    coverage["preferred_exact_fraction"] = coverage["exact_reference_fraction"]
    coverage["preferred_reference"] = "Ensembl Plants release 63"
    maize_path = metadata / "maize_v4_exact_id_coverage.tsv"
    if maize_path.exists():
        maize = pd.read_csv(maize_path, sep="\t")[["dataset_id", "exact_fraction", "reference_name", "reference_sha256", "source_url"]]
        maize = maize.rename(columns={
            "exact_fraction": "source_matched_exact_fraction",
            "reference_name": "source_matched_reference",
            "reference_sha256": "source_matched_reference_sha256",
            "source_url": "source_matched_source_url",
        })
        coverage = coverage.merge(maize, on="dataset_id", how="left")
        mask = coverage.source_matched_exact_fraction.notna()
        coverage.loc[mask, "preferred_exact_fraction"] = coverage.loc[mask, "source_matched_exact_fraction"]
        coverage.loc[mask, "preferred_reference"] = coverage.loc[mask, "source_matched_reference"]
    manifest = pd.read_csv(metadata / "dataset_manifest.tsv", sep="\t", keep_default_na=False)
    coverage["potential_reference_fraction"] = (
        coverage.preferred_exact_fraction + coverage.transform_candidate_fraction
    )
    coverage["triage_status"] = coverage.apply(triage, axis=1)
    coverage = coverage.merge(manifest[["dataset_id", "n_cells", "expression_class"]], on="dataset_id", how="left")
    coverage["accepted"] = False
    coverage.to_csv(metadata / "reference_mapping_action_queue.tsv", sep="\t", index=False)

    species = coverage.groupby("species").agg(
        datasets=("dataset_id", "count"),
        cells=("n_cells", "sum"),
        mean_exact=("preferred_exact_fraction", "mean"),
        minimum_exact=("preferred_exact_fraction", "min"),
        mean_potential=("potential_reference_fraction", "mean"),
        minimum_potential=("potential_reference_fraction", "min"),
    ).reset_index()
    statuses = coverage.groupby(["species", "triage_status"]).size().unstack(fill_value=0)
    species = species.merge(statuses, left_on="species", right_index=True, how="left")
    species.to_csv(metadata / "species_reference_coverage_summary.tsv", sep="\t", index=False)

    counts = coverage.triage_status.value_counts().to_dict()
    payload = {
        "species_measured": int(coverage.species.nunique()),
        "datasets_measured": len(coverage),
        "cells_represented": int(coverage.n_cells.sum()),
        "reference_release": "Ensembl Plants 63 plus source-matched MaizeGDB and UniProt candidates",
        "triage_counts": {str(k): int(v) for k, v in counts.items()},
        "accepted_mappings": 0,
        "private_identifiers_sent_externally": False,
    }
    (reports / "reference_coverage_phase1.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "# Reference coverage phase 1",
        "",
        "## Scope",
        "",
        f"Public Ensembl Plants release 63, source-matched MaizeGDB, and UniProt reference files were compared locally against {payload['datasets_measured']} H5AD files from {payload['species_measured']} species ({payload['cells_represented']:,} cells). No H5AD identifier was sent to an external service. No mapping is accepted automatically.",
        "",
        "## Species-level results",
        "",
        "| Species | Datasets | Mean exact | Minimum exact | Mean exact+transform candidate | Minimum potential |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in species.iterrows():
        lines.append(
            f"| {row.species} | {int(row.datasets)} | {row.mean_exact:.1%} | {row.minimum_exact:.1%} | {row.mean_potential:.1%} | {row.minimum_potential:.1%} |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "- Arabidopsis, Brassica rapa, and Oryza sativa already support high exact-ID coverage for most or all measured datasets.",
        "- Glycine max, Sorghum bicolor, and Catharanthus roseus have near-complete deterministic punctuation-format candidates; Populus trichocarpa contains one direct-v4 dataset and one explicit POPTR-v4 conversion candidate. These rules require provenance review before acceptance.",
        "- Manihot esculenta, Medicago truncatula, Nicotiana attenuata, and Pisum sativum require source-version or namespace-specific mappings. Current-release IDs must not be substituted by string similarity.",
        "- Zea mays data use Zm00001d identifiers. The source-matched MaizeGDB B73 v4 Zm00001d.1 annotation restores >90% exact coverage in all ten datasets; the current v5 reference is retained only to document the version mismatch.",
        "",
        "## Next gate",
        "",
        "Review deterministic mappings against source annotations, then acquire the source-matched maize B73 v4, Medicago A17, cassava, Nicotiana, and pea annotations. Orthogroup construction begins only after each included species has an accepted protein-to-gene mapping and at least 60% measured coverage.",
    ]
    (reports / "reference_coverage_phase1.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
