#!/usr/bin/env python3
"""Build the GSE232863 validation-v2 freeze without reading expression values."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(".")
META = PROJECT / "metadata"
REPORTS = PROJECT / "reports"
AUDIT_DIR = REPORTS / "gse232863_structure_audit_v1"
ORTHOGROUPS = Path(
    "data/runs/PhyloOpenCell/orthofinder_full7_v1/results/Results_Sep10/Orthogroups/Orthogroups.tsv"
)
SOURCE_RDS = Path(
    "data/external/GSE232863/GSE232863_scRNA_omics.Rds"
)
SOURCE_SHA256 = "7a308241e07575770f227057b1cd716b4cb9c255200c4ff4189c567bc2fbc958"
EXPECTED_BYTES = 5_521_719_712


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]):
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_values(text: str):
    return [value for value in text.split(";") if value]


def main():
    if SOURCE_RDS.stat().st_size != EXPECTED_BYTES:
        raise RuntimeError("source RDS size does not match the frozen official response size")

    features = [row["feature_id"] for row in read_tsv(AUDIT_DIR / "feature_ids.tsv")]
    feature_set = set(features)
    candidates = read_tsv(META / "conserved_program_candidates_v2.tsv")
    sample_counts = read_tsv(AUDIT_DIR / "sample_cluster_counts.tsv")

    with ORTHOGROUPS.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rice_column = "oryza_sativa.representative_proteins"
        rice_gene_to_og: dict[str, str] = {}
        rice_genes_by_og: dict[str, list[str]] = defaultdict(list)
        for row in reader:
            orthogroup = row["Orthogroup"]
            for entry in row[rice_column].split(","):
                entry = entry.strip()
                if not entry:
                    continue
                gene_id = entry.split("|", 1)[0]
                rice_gene_to_og[gene_id] = orthogroup
                rice_genes_by_og[orthogroup].append(gene_id)

    mapping_rows = []
    mapped_features_by_og: dict[str, list[str]] = defaultdict(list)
    for feature in features:
        orthogroup = rice_gene_to_og.get(feature, "")
        status = "exact_rice_gene_to_orthogroup" if orthogroup else "unmapped"
        mapping_rows.append(
            {
                "source_feature_id": feature,
                "normalized_feature_id": feature,
                "orthogroup_id": orthogroup,
                "mapping_status": status,
            }
        )
        if orthogroup:
            mapped_features_by_og[orthogroup].append(feature)

    grouped_candidates: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in candidates:
        grouped_candidates[(row["panel"], row["coarse_label"])].append(row)

    coverage_rows = []
    coverage_by_key: dict[tuple[str, str], dict] = {}
    for (panel, coarse_label), rows in sorted(grouped_candidates.items()):
        target_ogs = sorted({row["orthogroup_id"] for row in rows})
        mapped_ogs = [og for og in target_ogs if mapped_features_by_og.get(og)]
        single_gene_ogs = [og for og in mapped_ogs if len(mapped_features_by_og[og]) == 1]
        expanded_family_ogs = [og for og in mapped_ogs if len(mapped_features_by_og[og]) > 1]
        record = {
            "panel": panel,
            "coarse_label": coarse_label,
            "candidate_orthogroups": len(target_ogs),
            "orthogroups_with_rice_reference_members": sum(bool(rice_genes_by_og.get(og)) for og in target_ogs),
            "orthogroups_mapped_to_rds_features": len(mapped_ogs),
            "mapping_coverage_fraction": len(mapped_ogs) / len(target_ogs),
            "single_gene_orthogroups_in_rds": len(single_gene_ogs),
            "expanded_family_orthogroups_in_rds": len(expanded_family_ogs),
            "mapped_rds_features": sum(len(mapped_features_by_og[og]) for og in mapped_ogs),
        }
        coverage_rows.append(record)
        coverage_by_key[(panel, coarse_label)] = record

    definitions = [
        {
            "contrast_id": "rice_leaf_photosynthetic_v1",
            "panel": "leaf",
            "coarse_label": "photosynthetic_ground_tissue",
            "tissues": "leaf;Flag",
            "samples": "leaf_1;leaf_2;Flag_1;Flag_2",
            "positive_labels": "Mesophyll",
            "negative_labels": "Epidermis;Guard_cell;Vascular_cylinder",
            "role": "primary",
        },
        {
            "contrast_id": "rice_leaf_epidermal_v1",
            "panel": "leaf",
            "coarse_label": "epidermal_system",
            "tissues": "leaf;Flag",
            "samples": "leaf_1;leaf_2;Flag_1;Flag_2",
            "positive_labels": "Epidermis;Guard_cell",
            "negative_labels": "Mesophyll;Vascular_cylinder",
            "role": "primary",
        },
        {
            "contrast_id": "rice_root_epidermal_v1",
            "panel": "root",
            "coarse_label": "epidermal_system",
            "tissues": "Root",
            "samples": "Root_1;Root_2",
            "positive_labels": "Epidermis;epidermis_near_root_hair",
            "negative_labels": "Cortex;Endodermis;Vascular_cylinder;Root_cap",
            "role": "primary",
        },
        {
            "contrast_id": "rice_root_ground_v1",
            "panel": "root",
            "coarse_label": "ground_tissue",
            "tissues": "Root",
            "samples": "Root_1;Root_2",
            "positive_labels": "Cortex;Endodermis",
            "negative_labels": "Epidermis;epidermis_near_root_hair;Vascular_cylinder;Root_cap",
            "role": "primary",
        },
        {
            "contrast_id": "rice_root_stele_v1",
            "panel": "root",
            "coarse_label": "stele_lineage",
            "tissues": "Root",
            "samples": "Root_1;Root_2",
            "positive_labels": "Vascular_cylinder",
            "negative_labels": "Epidermis;epidermis_near_root_hair;Cortex;Endodermis;Root_cap",
            "role": "primary",
        },
        {
            "contrast_id": "rice_root_cap_v1",
            "panel": "root",
            "coarse_label": "root_cap",
            "tissues": "Root",
            "samples": "Root_1;Root_2",
            "positive_labels": "Root_cap",
            "negative_labels": "Epidermis;epidermis_near_root_hair;Cortex;Endodermis;Vascular_cylinder",
            "role": "sensitivity_only_one_candidate_orthogroup",
        },
    ]

    count_lookup: dict[tuple[str, str, str], int] = defaultdict(int)
    for row in sample_counts:
        count_lookup[(row["sample"], row["tissue"], row["cluster_names"])] += int(row["cell_count"])

    freeze_rows = []
    for definition in definitions:
        tissues = set(split_values(definition["tissues"]))
        samples = split_values(definition["samples"])
        positives = set(split_values(definition["positive_labels"]))
        negatives = set(split_values(definition["negative_labels"]))
        sample_positive_counts = {}
        sample_negative_counts = {}
        for sample in samples:
            sample_rows = [key for key in count_lookup if key[0] == sample and key[1] in tissues]
            sample_positive_counts[sample] = sum(count_lookup[key] for key in sample_rows if key[2] in positives)
            sample_negative_counts[sample] = sum(count_lookup[key] for key in sample_rows if key[2] in negatives)
        coverage = coverage_by_key[(definition["panel"], definition["coarse_label"])]
        min_positive = min(sample_positive_counts.values())
        min_negative = min(sample_negative_counts.values())
        eligible = (
            min_positive >= 100
            and min_negative >= 100
            and coverage["mapping_coverage_fraction"] >= 0.70
            and (coverage["candidate_orthogroups"] >= 5 or definition["role"] != "primary")
        )
        freeze_rows.append(
            {
                **definition,
                "positive_cells_total": sum(sample_positive_counts.values()),
                "negative_cells_total": sum(sample_negative_counts.values()),
                "minimum_positive_cells_per_sample": min_positive,
                "minimum_negative_cells_per_sample": min_negative,
                "positive_counts_by_sample": json.dumps(sample_positive_counts, sort_keys=True),
                "negative_counts_by_sample": json.dumps(sample_negative_counts, sort_keys=True),
                "candidate_orthogroups": coverage["candidate_orthogroups"],
                "mapped_candidate_orthogroups": coverage["orthogroups_mapped_to_rds_features"],
                "mapping_coverage_fraction": coverage["mapping_coverage_fraction"],
                "eligibility_status": "frozen_eligible" if eligible else "frozen_ineligible",
                "expression_inspected_for_selection": False,
            }
        )

    label_rows = []
    for definition in definitions:
        for source_label in split_values(definition["positive_labels"]):
            label_rows.append(
                {
                    "contrast_id": definition["contrast_id"],
                    "panel": definition["panel"],
                    "source_label_column": "cluster_names",
                    "source_label": source_label,
                    "coarse_label": definition["coarse_label"],
                    "contrast_role": "positive",
                    "mapping_basis": "pre_expression_match_to_frozen_coarse_hierarchy_v1",
                }
            )
        for source_label in split_values(definition["negative_labels"]):
            label_rows.append(
                {
                    "contrast_id": definition["contrast_id"],
                    "panel": definition["panel"],
                    "source_label_column": "cluster_names",
                    "source_label": source_label,
                    "coarse_label": "defined_other_system",
                    "contrast_role": "negative",
                    "mapping_basis": "pre_expression_cross_system_negative_control",
                }
            )

    mapping_path = META / "external_validation_v2_rice_feature_orthogroup_map_v1.tsv"
    coverage_path = META / "external_validation_v2_program_coverage_v1.tsv"
    freeze_path = META / "external_validation_v2_freeze_v1.tsv"
    labels_path = META / "external_validation_v2_label_map_v1.tsv"
    write_tsv(mapping_path, mapping_rows)
    write_tsv(coverage_path, coverage_rows)
    write_tsv(freeze_path, freeze_rows)
    write_tsv(labels_path, label_rows)

    protocol = "# External validation v2 frozen protocol\n\n"
    protocol += "Dataset: GSE232863, Oryza sativa, 115,395 cells across 16 samples and eight tissues. "
    protocol += "The object structure, labels, and gene identifiers were inspected before freeze; expression values and candidate-program outcomes were not inspected.\n\n"
    protocol += "## Primary endpoint\n\n"
    protocol += "Use the unchanged top-128 rank/orthogroup representation from validation v1. For each eligible contrast, report pooled cell-level AUROC/AP/effect as descriptive measures and sample-stratified AUROC/effect as the biological-replicate check. A contrast passes only if pooled AUROC > 0.60, standardized mean difference > 0.25, empirical random-program p <= 0.05, and the direction is positive (AUROC > 0.50) in at least 3 of 4 leaf samples or both root samples.\n\n"
    protocol += "## Overall interpretation\n\n"
    protocol += "Positive external replication requires at least two passing primary coarse systems spanning both leaf and root panels, with no eligible primary system showing a majority-replicate direction reversal. One passing system or incomplete threshold support is mixed evidence. Zero passing systems is external nonreplication under this frozen method.\n\n"
    protocol += "## Method comparison\n\n"
    protocol += "Compare complete orthogroup families, target-species single-gene orthogroups, target-species expanded families, strict single-copy families when non-empty, 500 size-matched random programs, and the pre-specified off-target program. A family-method advantage may be claimed only when its pooled AUROC exceeds the single-gene subset by at least 0.02 and its median sample-level AUROC is also higher. Biological replication does not depend on demonstrating this method advantage.\n\n"
    protocol += "## Boundaries\n\n"
    protocol += "Unknown, meristematic, exodermal, sclerenchymal, reproductive, and otherwise ambiguous labels are excluded rather than reassigned. Cell-level metrics do not substitute for sample-level consistency. No result establishes causal function or paralog replacement. No post-hoc relabeling, endpoint substitution, or v1 overwrite is allowed.\n"
    protocol_path = REPORTS / "external_validation_v2_frozen_protocol_v1.md"
    protocol_path.write_text(protocol, encoding="utf-8")

    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "external_validation_v2_frozen_before_expression_scoring",
        "dataset": "GSE232863",
        "source_rds": str(SOURCE_RDS),
        "source_bytes": EXPECTED_BYTES,
        "source_sha256": SOURCE_SHA256,
        "cells": 115395,
        "features": len(features),
        "mapped_features": sum(bool(row["orthogroup_id"]) for row in mapping_rows),
        "primary_contrasts": sum(row["role"] == "primary" for row in freeze_rows),
        "eligible_primary_contrasts": sum(
            row["role"] == "primary" and row["eligibility_status"] == "frozen_eligible"
            for row in freeze_rows
        ),
        "outputs": {
            str(path): sha256(path)
            for path in (mapping_path, coverage_path, freeze_path, labels_path, protocol_path)
        },
        "safety": {
            "source_rds_modified": False,
            "expression_values_inspected": False,
            "candidate_program_outcomes_inspected": False,
            "model_training_started": False,
            "v1_outputs_modified": False,
            "other_task_written": False,
        },
    }
    audit_path = REPORTS / "external_validation_v2_freeze_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
