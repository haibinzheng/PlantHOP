#!/usr/bin/env python3
"""Lock the interpretation of independent validation v1 without rerunning it."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/workspace/projects/PhyloOpenCell")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"
RUN = Path("/data/runs/PhyloOpenCell/independent_program_validation_v1")


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    summaries = {r["dataset_id"]: r for r in read_tsv(META / "independent_program_validation_summary_v1.tsv")}
    methods = read_tsv(META / "independent_program_method_comparison_v1.tsv")
    cells = read_tsv(RUN / "sampled_cell_scores.tsv")
    labels = defaultdict(Counter)
    for row in cells:
        labels[row["dataset_id"]][row["source_label"]] += 1

    interpretation = [
        {
            "dataset_id": "Populus trichocarpa_PRJCA014789",
            "validation_question": "Do frozen phloem-lineage orthogroups enrich Populus cells labelled Phloem relative to xylem/cambial cells?",
            "endpoint_identifiability": "identifiable_coarse_contrast",
            "result": "external_nonreplication",
            "evidence": "Complete program AUROC 0.3208, standardized difference -0.6032, and empirical random-program p=0.9741; the frozen program scores lower in labelled Phloem than in xylem/cambial negatives.",
            "negative_label_composition": "; ".join(f"{k}:{v}" for k, v in labels["Populus trichocarpa_PRJCA014789"].most_common() if k != "Phloem"),
            "allowed_conclusion": "The current phloem program does not generalize to this Populus wood-forming dataset under the frozen rank-orthogroup method.",
            "prohibited_conclusion": "Do not infer absence of phloem biology or invalidate all phloem markers; source tissue, developmental state, and label granularity remain possible moderators.",
        },
        {
            "dataset_id": "Arabidopsis thaliana_PRJNA796288",
            "validation_question": "Do frozen epidermal-system orthogroups distinguish mature leaf guard cells from other cells in a stomatal-lineage-enriched experiment?",
            "endpoint_identifiability": "not_identifiable_for_coarse_epidermal_program",
            "result": "technical_endpoint_mismatch_not_biological_nonreplication",
            "evidence": "The 3,000 negative cells are stomatal-lineage ground cells, young guard cells, meristemoids, guard mother cells, and proliferating cells; these largely share the epidermal/stomatal lineage with the positive class.",
            "negative_label_composition": "; ".join(f"{k}:{v}" for k, v in labels["Arabidopsis thaliana_PRJNA796288"].most_common() if k != "Leaf guard cell"),
            "allowed_conclusion": "This dataset can test fine stomatal-state markers but cannot test a coarse epidermal-versus-non-epidermal program.",
            "prohibited_conclusion": "Do not count AUROC 0.2854 as external failure of the coarse epidermal program.",
        },
    ]
    interpretation_path = META / "independent_validation_v1_interpretation.tsv"
    write_tsv(interpretation_path, interpretation)

    full = {r["dataset_id"]: r for r in methods if r["method"] == "complete_orthogroup_family_program"}
    restricted = {r["dataset_id"]: r for r in methods if r["method"] == "target_species_single_gene_orthogroups"}
    method_conclusion = {
        "strict_single_copy_coverage": {
            dataset: int(summaries[dataset]["strict_single_copy_orthogroups"]) for dataset in summaries
        },
        "complete_vs_target_species_single_gene_auc": {
            dataset: {
                "complete": float(full[dataset]["roc_auc"]),
                "single_gene_orthogroup_subset": float(restricted[dataset]["roc_auc"]),
            } for dataset in summaries
        },
        "conclusion": "Strict single-copy restriction retained zero frozen target orthogroups in both tests. The complete-family program did not outperform the target-species single-gene subset in either dataset; therefore no method advantage is supported by validation v1.",
    }

    matrix_v4 = read_tsv(META / "jeb_biology_evidence_matrix_v4.tsv")
    matrix_v5 = []
    for row in matrix_v4:
        updated = dict(row)
        if row["claim_id"] == "C3":
            updated["current_status"] = "discovery_supported_external_validation_not_yet_passed"
            updated["supporting_evidence"] += " In frozen external validation v1, the identifiable Populus phloem test did not replicate; the Arabidopsis guard-cell dataset was not identifiable for the coarse epidermal endpoint."
            updated["claim_boundary"] += " No positive independent external replication is currently available."
            updated["primary_source"] += "; metadata/independent_validation_v1_interpretation.tsv; reports/independent_validation_v1_interpretation_audit.json"
            updated["recommended_placement"] = "Discovery result only; external validation must be shown as a negative/boundary result unless a new pre-frozen cohort replicates."
            updated["next_required_evidence"] = "Acquire or curate at least one dataset containing two or more mapped coarse systems, preferably in a species absent from expression discovery."
        matrix_v5.append(updated)
    matrix_path = META / "jeb_biology_evidence_matrix_v5.tsv"
    write_tsv(matrix_path, matrix_v5)

    report_text = """# Independent validation v1 interpretation\n\n## Outcome\n\nThe first frozen external round did not provide positive independent replication. The Populus phloem contrast was biologically identifiable and failed under the frozen method. The Arabidopsis guard-cell contrast cannot answer the coarse epidermal question because its controls are almost entirely cells from the same stomatal/epidermal lineage.\n\n## Method comparison\n\nStrict single-copy filtering retained zero target orthogroups in both datasets, illustrating severe coverage loss. However, complete-family orthogroup scoring also failed to outperform the target-species single-gene subset, so this round does not demonstrate a benefit from complete gene families.\n\n## Consequence for JEB\n\nThe resource/governance and internal recurrence results remain valid, but “externally conserved programs” is not yet an accepted claim. The negative Populus result should remain visible as a boundary condition. A second validation round requires a pre-frozen dataset with at least two mapped coarse cell systems; single-lineage enrichment studies are unsuitable for this endpoint.\n\n## Data requirement\n\nThe existing local inventory contains no fully eligible unused dataset meeting all strict gates (resolved official provenance, adequate orthogroup coverage, at least two coarse labels with sufficient cells, and no discovery-study overlap). A new public dataset or a provenance-resolved existing candidate is therefore required for a decisive second round.\n"""
    report_path = REPORTS / "independent_validation_v1_interpretation.md"
    report_path.write_text(report_text, encoding="utf-8")
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "independent_validation_v1_interpreted_without_posthoc_relabeling",
        "identifiable_tests": 1,
        "positive_external_replications": 0,
        "external_nonreplications": 1,
        "endpoint_mismatches": 1,
        "method_comparison": method_conclusion,
        "outputs": {
            "interpretation": str(interpretation_path),
            "report": str(report_path),
            "updated_jeb_evidence_matrix": str(matrix_path),
        },
        "next_gate": "Freeze a multi-lineage, provenance-resolved independent dataset before expression inspection.",
        "safety": {"source_h5ad_accessed": False, "validation_rerun": False, "posthoc_label_change": False, "model_training_started": False},
    }
    audit_path = REPORTS / "independent_validation_v1_interpretation_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
