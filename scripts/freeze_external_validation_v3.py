#!/usr/bin/env python3
"""Freeze third-round external validation without reading expression values."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path("/workspace/projects/PhyloOpenCell")
META = PROJECT / "metadata"
REPORTS = PROJECT / "reports"
MORICANDIA = Path(
    "/data/datasets/Bioinformatics/external/PRJNA1186371_Moricandia_paper_object/m_arvensis_bbknn.h5ad"
)
WHEAT = Path(
    "/data/datasets/Bioinformatics/external/GSE270342_Triticum_aestivum/GSE270342_seuratObj_for_publication.rds.gz"
)
MORICANDIA_SHA256 = "be5f9be6570b78029fc138dd9c81ad908532776bbe665b9356b2055f554fd802"
WHEAT_SHA256 = "7c0054b23b782e538ead563fb1cbe1688db89f185753cf5ff14ccd2cb72fd58f"
RANDOM_SEED = 20260915
RANDOM_PROGRAMS = 500


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
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


def json_load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def split_values(text: str):
    return [value for value in text.split(";") if value]


def sum_labels(counts, labels):
    return sum(int(counts.get(label, 0)) for label in labels)


def main():
    source_hashes = {
        "moricandia": sha256(MORICANDIA),
        "wheat": sha256(WHEAT),
    }
    if source_hashes["moricandia"] != MORICANDIA_SHA256:
        raise RuntimeError("Moricandia source hash changed")
    if source_hashes["wheat"] != WHEAT_SHA256:
        raise RuntimeError("wheat source hash changed")

    candidates = read_tsv(META / "conserved_program_candidates_v2.tsv")
    programs = defaultdict(set)
    for row in candidates:
        programs[(row["panel"], row["coarse_label"])].add(row["orthogroup_id"])

    mor_map = read_tsv(META / "external_validation_v3_moricandia_feature_orthogroup_map_v1.tsv")
    wheat_map = read_tsv(META / "external_validation_v3_wheat_feature_orthogroup_map_v1.tsv")
    mapped_ogs = {
        "moricandia": {row["project_orthogroup_id"] for row in mor_map if row["project_orthogroup_id"]},
        "wheat": {
            row["project_orthogroup_id"]
            for row in wheat_map
            if row["mapping_status"].startswith("unanimous_") and row["project_orthogroup_id"]
        },
    }

    mor_counts = json_load(REPORTS / "moricandia_paper_object_count_linkage_audit_v1.json")
    wheat_audit = json_load(REPORTS / "wheat_gse270342_seurat_structure_audit_v1.json")
    count_tables = {
        "moricandia": mor_counts["cell_type_by_replicate"],
        "wheat": wheat_audit["annotation_by_replicate_named"],
    }
    sample_names = {
        "moricandia": ["0", "1", "2"],
        "wheat": ["wheat-root-1", "wheat-root-2", "wheat-root-3"],
    }

    definitions = [
        {
            "contrast_id": "moricandia_leaf_photosynthetic_v1",
            "dataset": "moricandia",
            "panel": "leaf",
            "coarse_label": "photosynthetic_ground_tissue",
            "positive_labels": "Spongy mesophyll cell;Palisade and spongy mesophyll",
            "negative_labels": "Adaxial pavement cell;Phloem parenchyma;Companion cell",
            "excluded_labels": "Bundle sheath;Mesophyll cell/hydathode cell",
            "role": "primary",
            "off_target_program": "leaf:epidermal_system",
        },
        {
            "contrast_id": "moricandia_leaf_epidermal_v1",
            "dataset": "moricandia",
            "panel": "leaf",
            "coarse_label": "epidermal_system",
            "positive_labels": "Adaxial pavement cell",
            "negative_labels": "Spongy mesophyll cell;Palisade and spongy mesophyll;Bundle sheath;Phloem parenchyma;Companion cell",
            "excluded_labels": "Mesophyll cell/hydathode cell",
            "role": "primary",
            "off_target_program": "leaf:photosynthetic_ground_tissue",
        },
        {
            "contrast_id": "moricandia_leaf_phloem_v1",
            "dataset": "moricandia",
            "panel": "vascular",
            "coarse_label": "phloem_lineage",
            "positive_labels": "Phloem parenchyma;Companion cell",
            "negative_labels": "Spongy mesophyll cell;Palisade and spongy mesophyll;Adaxial pavement cell;Bundle sheath",
            "excluded_labels": "Mesophyll cell/hydathode cell",
            "role": "sensitivity_sparse_one_replicate",
            "off_target_program": "vascular:xylem_lineage",
        },
        {
            "contrast_id": "wheat_root_epidermal_v1",
            "dataset": "wheat",
            "panel": "root",
            "coarse_label": "epidermal_system",
            "positive_labels": "Epidermis;Root Hair",
            "negative_labels": "Cortex;Endodermis;Provascular cells;Pericycle;Xylem;Root Cap",
            "excluded_labels": "Unknown;Meristems;Dividing Cells;Endodermis/Phloem;Phloem",
            "role": "primary_candidate",
            "off_target_program": "root:ground_tissue",
        },
        {
            "contrast_id": "wheat_root_ground_v1",
            "dataset": "wheat",
            "panel": "root",
            "coarse_label": "ground_tissue",
            "positive_labels": "Cortex;Endodermis",
            "negative_labels": "Epidermis;Root Hair;Provascular cells;Pericycle;Xylem;Root Cap",
            "excluded_labels": "Unknown;Meristems;Dividing Cells;Endodermis/Phloem;Phloem",
            "role": "primary_candidate",
            "off_target_program": "root:epidermal_system",
        },
        {
            "contrast_id": "wheat_root_stele_v1",
            "dataset": "wheat",
            "panel": "root",
            "coarse_label": "stele_lineage",
            "positive_labels": "Provascular cells;Pericycle;Xylem",
            "negative_labels": "Epidermis;Root Hair;Cortex;Endodermis;Root Cap",
            "excluded_labels": "Unknown;Meristems;Dividing Cells;Endodermis/Phloem;Phloem",
            "role": "primary",
            "off_target_program": "root:epidermal_system",
        },
        {
            "contrast_id": "wheat_root_xylem_v1",
            "dataset": "wheat",
            "panel": "vascular",
            "coarse_label": "xylem_lineage",
            "positive_labels": "Xylem",
            "negative_labels": "Epidermis;Root Hair;Cortex;Endodermis;Root Cap",
            "excluded_labels": "Unknown;Meristems;Dividing Cells;Endodermis/Phloem;Phloem;Provascular cells;Pericycle",
            "role": "sensitivity_sparse_replicates",
            "off_target_program": "vascular:phloem_lineage",
        },
        {
            "contrast_id": "wheat_root_phloem_v1",
            "dataset": "wheat",
            "panel": "vascular",
            "coarse_label": "phloem_lineage",
            "positive_labels": "Phloem",
            "negative_labels": "Epidermis;Root Hair;Cortex;Endodermis;Root Cap",
            "excluded_labels": "Unknown;Meristems;Dividing Cells;Endodermis/Phloem;Provascular cells;Pericycle;Xylem",
            "role": "sensitivity_sparse_replicates",
            "off_target_program": "vascular:xylem_lineage",
        },
        {
            "contrast_id": "wheat_root_cap_v1",
            "dataset": "wheat",
            "panel": "root",
            "coarse_label": "root_cap",
            "positive_labels": "Root Cap",
            "negative_labels": "Epidermis;Root Hair;Cortex;Endodermis;Provascular cells;Pericycle;Xylem",
            "excluded_labels": "Unknown;Meristems;Dividing Cells;Endodermis/Phloem;Phloem",
            "role": "sensitivity_sparse_replicates_single_candidate_orthogroup",
            "off_target_program": "root:ground_tissue",
        },
    ]

    coverage_rows = []
    for dataset in ("moricandia", "wheat"):
        for (panel, label), orthogroups in sorted(programs.items()):
            covered = orthogroups & mapped_ogs[dataset]
            coverage_rows.append(
                {
                    "dataset": dataset,
                    "panel": panel,
                    "coarse_label": label,
                    "candidate_orthogroups": len(orthogroups),
                    "mapped_candidate_orthogroups": len(covered),
                    "mapping_coverage_fraction": len(covered) / len(orthogroups),
                    "expression_inspected": False,
                }
            )
    coverage_lookup = {
        (row["dataset"], row["panel"], row["coarse_label"]): row for row in coverage_rows
    }

    freeze_rows = []
    label_rows = []
    for definition in definitions:
        dataset = definition["dataset"]
        positives = split_values(definition["positive_labels"])
        negatives = split_values(definition["negative_labels"])
        exclusions = split_values(definition["excluded_labels"])
        positive_by_sample = {}
        negative_by_sample = {}
        excluded_by_sample = {}
        for sample in sample_names[dataset]:
            sample_counts = {
                label: int(values.get(sample, 0))
                for label, values in count_tables[dataset].items()
            }
            positive_by_sample[sample] = sum_labels(sample_counts, positives)
            negative_by_sample[sample] = sum_labels(sample_counts, negatives)
            excluded_by_sample[sample] = sum_labels(sample_counts, exclusions)

        coverage = coverage_lookup[(dataset, definition["panel"], definition["coarse_label"])]
        basic_eligible = (
            min(positive_by_sample.values()) >= 100
            and min(negative_by_sample.values()) >= 100
            and int(coverage["mapped_candidate_orthogroups"]) >= 5
            and float(coverage["mapping_coverage_fraction"]) >= 0.70
        )
        eligible = basic_eligible and definition["role"] in {"primary", "primary_candidate"}
        final_role = "primary" if eligible else definition["role"]
        if definition["role"] == "primary_candidate" and not eligible:
            final_role = "frozen_ineligible_mapping_coverage"
        freeze_rows.append(
            {
                **definition,
                "role": final_role,
                "samples": ";".join(sample_names[dataset]),
                "positive_cells_total": sum(positive_by_sample.values()),
                "negative_cells_total": sum(negative_by_sample.values()),
                "excluded_cells_total": sum(excluded_by_sample.values()),
                "minimum_positive_cells_per_sample": min(positive_by_sample.values()),
                "minimum_negative_cells_per_sample": min(negative_by_sample.values()),
                "positive_counts_by_sample": json.dumps(positive_by_sample, sort_keys=True),
                "negative_counts_by_sample": json.dumps(negative_by_sample, sort_keys=True),
                "excluded_counts_by_sample": json.dumps(excluded_by_sample, sort_keys=True),
                "candidate_orthogroups": coverage["candidate_orthogroups"],
                "mapped_candidate_orthogroups": coverage["mapped_candidate_orthogroups"],
                "mapping_coverage_fraction": coverage["mapping_coverage_fraction"],
                "eligibility_status": "frozen_eligible" if eligible else "frozen_ineligible_or_sensitivity",
                "expression_inspected_for_selection": False,
            }
        )
        for contrast_role, labels in (
            ("positive", positives), ("negative", negatives), ("excluded", exclusions)
        ):
            for label in labels:
                label_rows.append(
                    {
                        "contrast_id": definition["contrast_id"],
                        "dataset": dataset,
                        "source_label_column": "cell_type" if dataset == "moricandia" else "annotation",
                        "source_label": label,
                        "coarse_label": definition["coarse_label"] if contrast_role == "positive" else ("defined_other_system" if contrast_role == "negative" else "unresolved_or_out_of_scope"),
                        "contrast_role": contrast_role,
                        "mapping_basis": "author_label_to_predeclared_coarse_hierarchy_before_expression_scoring",
                    }
                )

    coverage_path = META / "external_validation_v3_program_coverage_v1.tsv"
    freeze_path = META / "external_validation_v3_freeze_v1.tsv"
    labels_path = META / "external_validation_v3_label_map_v1.tsv"
    write_tsv(coverage_path, coverage_rows)
    write_tsv(freeze_path, freeze_rows)
    write_tsv(labels_path, label_rows)

    protocol = f"""# External validation v3 frozen protocol

Frozen before candidate-program expression scoring. This round uses two author-labelled, previously unseen queues: Moricandia arvensis leaf nuclei (10,129 retained cells, three replicates) and Triticum aestivum root cells (7,388 cells, three replicates). Source expression objects are read-only.

## Primary contrasts

Three contrasts are eligible and frozen as primary: Moricandia photosynthetic/ground tissue, Moricandia epidermal system, and wheat stele lineage. Wheat epidermal and ground contrasts fail the pre-score mapping gate and remain ineligible; sparse phloem, xylem, root-cap, and Moricandia phloem contrasts are sensitivity-only. Ambiguous or out-of-scope labels are explicitly excluded and may not be reassigned after scoring.

## Mapping

Moricandia features use the author-provided Arabidopsis reference gene stored in the H5AD variable index and exact membership in the frozen project OrthoFinder namespace; duplicate Moricandia contigs are retained as family members. Wheat uses the paper author's orthogroup table (SHA-256 63ede47fdec74f7f7fedad4939ebba1ec26e405e88ff3c282b1e0a8e1641effc), bridged to the frozen project namespace only when all exact shared Arabidopsis/rice genes unanimously indicate one project orthogroup. Ambiguous and unsupported families remain unresolved.

## Eligibility gate

A primary contrast requires at least 100 positive and 100 negative cells in every biological replicate, at least five mapped candidate orthogroups, and candidate-program mapping coverage of at least 0.70. This gate uses labels and identifiers only.

## Frozen scoring endpoint

Use the unchanged top-128 within-cell rank/orthogroup scoring representation from validation v2. Pooled cell-level AUROC, average precision, and standardized mean difference are descriptive. A primary contrast passes only when pooled AUROC is greater than 0.60, standardized mean difference is greater than 0.25, a {RANDOM_PROGRAMS}-draw size-matched empirical random-program p-value is at most 0.05, and AUROC is greater than 0.50 in all three biological replicates. Random seed: {RANDOM_SEED}.

Overall third-round replication requires at least two of the three eligible primary contrasts to pass, the passing set must span both leaf and root, and no eligible primary contrast may show a majority-replicate direction reversal. Because only wheat stele is root-eligible, a positive overall result necessarily requires wheat stele plus at least one Moricandia leaf contrast.

## Method and negative controls

Report complete mapped orthogroup families, target-species single-gene orthogroups, expanded-family orthogroups, and a strict wheat two-reference-species crosswalk sensitivity. Compare against {RANDOM_PROGRAMS} size-matched random programs, the contrast-specific frozen off-target program, and an OMG-style unweighted orthogroup-set baseline in which each orthogroup contributes once regardless of family size. Family advantage is secondary and requires pooled AUROC at least 0.02 above the single-gene subset plus higher median replicate AUROC; biological replication does not require family advantage.

## Boundaries

No post-hoc relabelling, endpoint substitution, mapping relaxation, candidate replacement, or overwrite of v1/v2 outputs is allowed. Cell-level significance cannot substitute for replicate consistency. Results establish external association and conservation only, not causal function, paralog replacement, or field-level performance. Deep learning, scVI/scANVI, and strong-model benchmarking remain outside this round.
"""
    protocol_path = REPORTS / "external_validation_v3_frozen_protocol_v1.md"
    protocol_path.write_text(protocol, encoding="utf-8")

    output_paths = [
        META / "external_validation_v3_moricandia_feature_orthogroup_map_v1.tsv",
        META / "external_validation_v3_wheat_feature_ids_v1.tsv",
        META / "external_validation_v3_wheat_feature_orthogroup_map_v1.tsv",
        coverage_path,
        freeze_path,
        labels_path,
        protocol_path,
    ]
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "external_validation_v3_frozen_before_expression_scoring",
        "sources": {
            "moricandia": {"path": str(MORICANDIA), "sha256": source_hashes["moricandia"]},
            "wheat": {"path": str(WHEAT), "sha256": source_hashes["wheat"]},
        },
        "random_seed": RANDOM_SEED,
        "random_programs": RANDOM_PROGRAMS,
        "eligible_primary_contrasts": [
            row["contrast_id"] for row in freeze_rows if row["eligibility_status"] == "frozen_eligible"
        ],
        "outputs": {str(path): sha256(path) for path in output_paths},
        "safety": {
            "source_objects_modified": False,
            "expression_values_inspected": False,
            "candidate_program_outcomes_inspected": False,
            "model_training_started": False,
            "v1_or_v2_outputs_modified": False,
            "other_task_written": False,
        },
    }
    audit_path = REPORTS / "external_validation_v3_freeze_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
