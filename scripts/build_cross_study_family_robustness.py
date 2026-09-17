#!/usr/bin/env python3
"""Leave-one-study-out robustness checks for priority orthogroup families.

For each orthogroup and species with at least two datasets, the dominant gene
is selected using all datasets except one. Its standardized mean difference is
then evaluated in the omitted dataset. This is an internal study-robustness
check, not external validation. Existing source and run data are read only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(".")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def median(values):
    return statistics.median(values) if values else math.nan


def main():
    effects = read_tsv(META / "priority_paralog_dataset_effects_v2.tsv")
    tree_rows = read_tsv(META / "priority_paralog_gene_tree_evidence_v2.tsv")
    matrix_v1 = read_tsv(META / "jeb_biology_evidence_matrix_v1.tsv")
    tree_by_og = {row["orthogroup_id"]: row for row in tree_rows}

    grouped = defaultdict(list)
    for row in effects:
        grouped[(row["panel"], row["coarse_label"], row["orthogroup_id"], row["species"])].append(row)

    detail = []
    species_summary = []
    for (panel, label, og, species), rows in sorted(grouped.items()):
        datasets = sorted({row["dataset_id"] for row in rows})
        if len(datasets) < 2:
            species_summary.append({
                "panel": panel,
                "coarse_label": label,
                "orthogroup_id": og,
                "species": species,
                "datasets_available": str(len(datasets)),
                "heldout_tests": "0",
                "heldout_positive_tests": "0",
                "heldout_positive_fraction": "",
                "median_heldout_standardized_difference": "",
                "dominant_gene_stability": "",
                "study_robustness_status": "insufficient_multi_study_coverage",
            })
            continue

        by_gene_dataset = {(r["reference_gene_id"], r["dataset_id"]): float(r["standardized_mean_difference"]) for r in rows}
        genes = sorted({row["reference_gene_id"] for row in rows})
        selected_genes = []
        heldout_effects = []
        for heldout in datasets:
            scored = []
            for gene in genes:
                train_values = [by_gene_dataset[(gene, d)] for d in datasets if d != heldout and (gene, d) in by_gene_dataset]
                if not train_values:
                    continue
                score = median(train_values)
                positive_fraction = sum(v > 0 for v in train_values) / len(train_values)
                scored.append((score, positive_fraction, len(train_values), gene))
            if not scored:
                continue
            _, _, training_n, selected_gene = max(scored)
            heldout_value = by_gene_dataset.get((selected_gene, heldout))
            if heldout_value is None:
                detail.append({
                    "panel": panel, "coarse_label": label, "orthogroup_id": og, "species": species,
                    "heldout_dataset": heldout, "training_datasets_with_gene": str(training_n),
                    "selected_gene_without_heldout": selected_gene, "heldout_standardized_mean_difference": "",
                    "heldout_positive": "", "evaluation_status": "selected_gene_missing_in_heldout",
                })
                continue
            selected_genes.append(selected_gene)
            heldout_effects.append(heldout_value)
            detail.append({
                "panel": panel, "coarse_label": label, "orthogroup_id": og, "species": species,
                "heldout_dataset": heldout, "training_datasets_with_gene": str(training_n),
                "selected_gene_without_heldout": selected_gene,
                "heldout_standardized_mean_difference": f"{heldout_value:.9g}",
                "heldout_positive": str(heldout_value > 0), "evaluation_status": "evaluated",
            })
        n = len(heldout_effects)
        positive_n = sum(v > 0 for v in heldout_effects)
        positive_fraction = positive_n / n if n else math.nan
        top_gene_count = Counter(selected_genes).most_common(1)[0][1] if selected_genes else 0
        stability = top_gene_count / len(selected_genes) if selected_genes else math.nan
        if n < 2:
            status = "insufficient_evaluable_holdouts"
        elif positive_fraction >= 0.75 and median(heldout_effects) > 0:
            status = "study_robust"
        elif positive_fraction >= 0.5 and median(heldout_effects) > 0:
            status = "mixed_positive"
        else:
            status = "not_study_robust"
        species_summary.append({
            "panel": panel, "coarse_label": label, "orthogroup_id": og, "species": species,
            "datasets_available": str(len(datasets)), "heldout_tests": str(n),
            "heldout_positive_tests": str(positive_n),
            "heldout_positive_fraction": f"{positive_fraction:.6g}" if n else "",
            "median_heldout_standardized_difference": f"{median(heldout_effects):.9g}" if n else "",
            "dominant_gene_stability": f"{stability:.6g}" if selected_genes else "",
            "study_robustness_status": status,
        })

    family_summary = []
    for og in sorted({r["orthogroup_id"] for r in species_summary}):
        rows = [r for r in species_summary if r["orthogroup_id"] == og]
        evaluable = [r for r in rows if int(r["heldout_tests"]) >= 2]
        robust = [r for r in evaluable if r["study_robustness_status"] == "study_robust"]
        mixed = [r for r in evaluable if r["study_robustness_status"] == "mixed_positive"]
        heldout = [d for d in detail if d["orthogroup_id"] == og and d["evaluation_status"] == "evaluated"]
        values = [float(d["heldout_standardized_mean_difference"]) for d in heldout]
        if len(robust) >= 3:
            status = "cross_species_study_robust"
        elif len(robust) + len(mixed) >= 3:
            status = "cross_species_mixed_support"
        elif len(evaluable) < 3:
            status = "insufficient_multispecies_study_coverage"
        else:
            status = "not_cross_species_study_robust"
        tr = tree_by_og.get(og, {})
        family_summary.append({
            "panel": rows[0]["panel"], "coarse_label": rows[0]["coarse_label"], "orthogroup_id": og,
            "species_in_priority_set": str(len(rows)), "species_with_multistudy_evaluation": str(len(evaluable)),
            "study_robust_species": str(len(robust)), "mixed_positive_species": str(len(mixed)),
            "heldout_tests": str(len(values)),
            "heldout_positive_fraction": f"{sum(v > 0 for v in values) / len(values):.6g}" if values else "",
            "median_heldout_standardized_difference": f"{median(values):.9g}" if values else "",
            "cross_study_family_status": status,
            "gene_tree_pattern": tr.get("gene_tree_pattern", ""),
            "paralog_substitution_claim": "not_established",
        })

    detail_path = META / "priority_family_leave_one_study_out_v1.tsv"
    species_path = META / "priority_family_study_robustness_by_species_v1.tsv"
    family_path = META / "priority_family_study_robustness_summary_v1.tsv"
    write_tsv(detail_path, detail, list(detail[0]))
    write_tsv(species_path, species_summary, list(species_summary[0]))
    write_tsv(family_path, family_summary, list(family_summary[0]))

    # Preserve v1 and publish a versioned evidence matrix with the new internal
    # robustness evidence. C3 remains discovery evidence, not external proof.
    family_counts = Counter(r["cross_study_family_status"] for r in family_summary)
    matrix_v2 = []
    for row in matrix_v1:
        new = dict(row)
        if row["claim_id"] == "C3":
            new["supporting_evidence"] += (
                f" Internal leave-one-study-out analysis classified {family_counts.get('cross_species_study_robust', 0)} of 19 priority families as cross-species study-robust and "
                f"{family_counts.get('cross_species_mixed_support', 0)} as mixed support."
            )
            new["claim_boundary"] += " This is internal cross-study robustness within the governed resource, not external validation."
            new["primary_source"] += "; metadata/priority_family_study_robustness_summary_v1.tsv"
        if row["claim_id"] in {"C6", "C7"}:
            og = "OG0003912" if row["claim_id"] == "C7" else "OG0000378/OG0000114"
            new["next_required_evidence"] += f" Consult {og} leave-one-study-out rows before promoting the claim."
        matrix_v2.append(new)
    matrix_path = META / "jeb_biology_evidence_matrix_v2.tsv"
    write_tsv(matrix_path, matrix_v2, list(matrix_v2[0]))

    selected = [r for r in family_summary if r["orthogroup_id"] in {"OG0000378", "OG0000114", "OG0003912"}]
    evaluation_status_counts = Counter(r["evaluation_status"] for r in detail)
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "internal_leave_one_study_out_complete",
        "method": "Within each orthogroup and species, select the dominant gene on all but one dataset and evaluate its standardized class-versus-other effect in the omitted dataset.",
        "family_status_counts": dict(family_counts),
        "orthogroups": len(family_summary),
        "species_orthogroup_rows": len(species_summary),
        "heldout_detail_rows": len(detail),
        "evaluation_status_counts": dict(evaluation_status_counts),
        "selected_family_results": selected,
        "interpretation_boundary": "Internal study-level robustness reduces same-study overfitting risk but is not an external cohort validation and does not establish function or paralog substitution.",
        "outputs": {str(p.relative_to(ROOT)): {"sha256": sha256(p)} for p in [detail_path, species_path, family_path, matrix_path]},
        "safety": {
            "source_h5ad_accessed": False,
            "model_training_started": False,
            "external_queries": False,
            "other_task_written": False,
            "existing_v1_run_overwritten": False,
        },
    }
    audit_path = REPORTS / "priority_family_study_robustness_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
