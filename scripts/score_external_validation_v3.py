#!/usr/bin/env python3
"""Score the frozen third-round external validation from encoded matrices."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.io import mmread


PROJECT = Path(".")
META = PROJECT / "metadata"
REPORTS = PROJECT / "reports"
RUN_DIR = Path("data/runs/PhyloOpenCell/independent_program_validation_v3")
PARTIAL_DIR = Path(str(RUN_DIR) + ".partial")
ENCODED_DIRS = {
    "moricandia": Path("data/derived/PhyloOpenCell/external_validation_v3/moricandia_top128_formal"),
    "wheat": Path("data/derived/PhyloOpenCell/external_validation_v3/wheat_top128_formal"),
}
RANDOM_SEED = 20260915
RANDOM_PROGRAMS = 500


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_values(value):
    return [part for part in value.split(";") if part]


def roc_auc(y, score):
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=float)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if not n_pos or not n_neg:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    sorted_scores = score[order]
    ranks = np.empty(len(score), dtype=float)
    start = 0
    while start < len(score):
        end = start + 1
        while end < len(score) and sorted_scores[end] == sorted_scores[start]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2
        start = end
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(y, score):
    y = np.asarray(y, dtype=np.int8)
    if not int((y == 1).sum()):
        return float("nan")
    order = np.lexsort((np.arange(len(score)), -np.asarray(score)))
    sorted_y = y[order]
    precision = np.cumsum(sorted_y == 1) / np.arange(1, len(y) + 1)
    return float(precision[sorted_y == 1].mean())


def standardized_effect(y, score):
    positive = np.asarray(score)[np.asarray(y) == 1]
    negative = np.asarray(score)[np.asarray(y) == 0]
    pooled = np.sqrt((positive.var(ddof=1) + negative.var(ddof=1)) / 2)
    if not np.isfinite(pooled) or pooled == 0:
        return 0.0
    return float((positive.mean() - negative.mean()) / pooled)


def metrics(y, score):
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=float)
    return {
        "roc_auc": roc_auc(y, score),
        "average_precision": average_precision(y, score),
        "standardized_mean_difference": standardized_effect(y, score),
        "mean_positive": float(score[y == 1].mean()),
        "mean_negative": float(score[y == 0].mean()),
    }


def load_encoded(dataset):
    directory = ENCODED_DIRS[dataset]
    matrix = mmread(directory / "encoded.mtx").tocsr()
    orthogroups = (directory / "orthogroups.txt").read_text(encoding="utf-8").splitlines()
    metadata = read_tsv(directory / "cell_metadata.tsv")
    if matrix.shape != (len(metadata), len(orthogroups)):
        raise RuntimeError(f"{dataset}: encoded dimensions do not match metadata/vocabulary")
    if not np.allclose(np.sqrt(matrix.multiply(matrix).sum(axis=1)).A1, 1.0):
        raise RuntimeError(f"{dataset}: encoded rows are not L2-normalized")
    return matrix, orthogroups, metadata


def main():
    if RUN_DIR.exists() or PARTIAL_DIR.exists():
        raise RuntimeError("refusing to overwrite existing v3 scoring run")
    PARTIAL_DIR.mkdir(parents=True)
    started = datetime.now(timezone.utc).isoformat()

    freeze = read_tsv(META / "external_validation_v3_freeze_v1.tsv")
    candidates = read_tsv(META / "conserved_program_candidates_v2.tsv")
    candidate_sets = defaultdict(set)
    for row in candidates:
        candidate_sets[(row["panel"], row["coarse_label"])].add(row["orthogroup_id"])

    feature_counts = {}
    strict_wheat_ogs = set()
    mor_rows = read_tsv(META / "external_validation_v3_moricandia_feature_orthogroup_map_v1.tsv")
    mor_counter = Counter(row["project_orthogroup_id"] for row in mor_rows if row["project_orthogroup_id"])
    feature_counts["moricandia"] = mor_counter
    wheat_rows = read_tsv(META / "external_validation_v3_wheat_feature_orthogroup_map_v1.tsv")
    accepted_wheat = [row for row in wheat_rows if row["mapping_status"].startswith("unanimous_") and row["project_orthogroup_id"]]
    wheat_counter = Counter(row["project_orthogroup_id"] for row in accepted_wheat)
    feature_counts["wheat"] = wheat_counter
    statuses_by_og = defaultdict(set)
    for row in accepted_wheat:
        statuses_by_og[row["project_orthogroup_id"]].add(row["mapping_status"])
    strict_wheat_ogs = {
        orthogroup for orthogroup, statuses in statuses_by_og.items()
        if statuses == {"unanimous_two_reference_species"}
    }

    rng = np.random.default_rng(RANDOM_SEED)
    pooled_rows = []
    replicate_rows = []
    summary_rows = []
    for dataset in ("moricandia", "wheat"):
        matrix, orthogroups, metadata = load_encoded(dataset)
        og_to_col = {orthogroup: index for index, orthogroup in enumerate(orthogroups)}
        labels = np.asarray([row["source_label"] for row in metadata])
        replicates = np.asarray([row["replicate"] for row in metadata])
        available_ogs = set(orthogroups)

        for frozen in [row for row in freeze if row["dataset"] == dataset]:
            positives = set(split_values(frozen["positive_labels"]))
            negatives = set(split_values(frozen["negative_labels"]))
            selected = np.flatnonzero(np.isin(labels, list(positives | negatives)))
            contrast_matrix = matrix[selected]
            contrast_labels = labels[selected]
            contrast_replicates = replicates[selected]
            y = np.isin(contrast_labels, list(positives)).astype(np.int8)
            target = sorted(candidate_sets[(frozen["panel"], frozen["coarse_label"])] & available_ogs)
            off_panel, off_label = frozen["off_target_program"].split(":", 1)
            off_target = sorted(candidate_sets[(off_panel, off_label)] & available_ogs)
            single = [og for og in target if feature_counts[dataset][og] == 1]
            expanded = [og for og in target if feature_counts[dataset][og] > 1]
            strict = [og for og in target if dataset == "wheat" and og in strict_wheat_ogs]

            method_sets = {
                "complete_orthogroup_family_program": target,
                "target_species_single_gene_orthogroups": single,
                "target_species_expanded_family_orthogroups": expanded,
                "off_target_program_negative_control": off_target,
            }
            if dataset == "wheat":
                method_sets["strict_two_reference_species_crosswalk"] = strict

            complete_score = None
            complete_replicate_aucs = []
            for method, method_ogs in method_sets.items():
                columns = [og_to_col[og] for og in method_ogs]
                score = np.asarray(contrast_matrix[:, columns].sum(axis=1)).ravel() if columns else np.zeros(len(y))
                if method == "complete_orthogroup_family_program":
                    complete_score = score
                record = metrics(y, score)
                pooled_rows.append({
                    "contrast_id": frozen["contrast_id"], "dataset": dataset,
                    "panel": frozen["panel"], "coarse_label": frozen["coarse_label"],
                    "role": frozen["role"], "method": method,
                    "candidate_orthogroups": len(method_ogs),
                    "positive_cells": int((y == 1).sum()), "negative_cells": int((y == 0).sum()),
                    **record,
                })
                for replicate in split_values(frozen["samples"]):
                    mask = contrast_replicates == replicate
                    replicate_record = metrics(y[mask], score[mask])
                    replicate_rows.append({
                        "contrast_id": frozen["contrast_id"], "dataset": dataset,
                        "replicate": replicate, "method": method,
                        "positive_cells": int((y[mask] == 1).sum()),
                        "negative_cells": int((y[mask] == 0).sum()),
                        **replicate_record,
                    })
                    if method == "complete_orthogroup_family_program":
                        complete_replicate_aucs.append(replicate_record["roc_auc"])

            binary = contrast_matrix[:, [og_to_col[og] for og in target]].copy() if target else sp.csr_matrix((len(y), 0))
            if binary.nnz:
                binary.data[:] = 1
            omg_score = np.asarray(binary.sum(axis=1)).ravel() / max(len(target), 1)
            pooled_rows.append({
                "contrast_id": frozen["contrast_id"], "dataset": dataset,
                "panel": frozen["panel"], "coarse_label": frozen["coarse_label"],
                "role": frozen["role"], "method": "omg_style_unweighted_orthogroup_presence",
                "candidate_orthogroups": len(target),
                "positive_cells": int((y == 1).sum()), "negative_cells": int((y == 0).sum()),
                **metrics(y, omg_score),
            })
            for replicate in split_values(frozen["samples"]):
                mask = contrast_replicates == replicate
                replicate_rows.append({
                    "contrast_id": frozen["contrast_id"], "dataset": dataset,
                    "replicate": replicate, "method": "omg_style_unweighted_orthogroup_presence",
                    "positive_cells": int((y[mask] == 1).sum()),
                    "negative_cells": int((y[mask] == 0).sum()),
                    **metrics(y[mask], omg_score[mask]),
                })

            random_pool = sorted(available_ogs - set(target) - set(off_target))
            random_size = len(target)
            random_aucs = np.empty(RANDOM_PROGRAMS)
            random_columns = []
            for _ in range(RANDOM_PROGRAMS):
                chosen = rng.choice(random_pool, size=random_size, replace=False)
                random_columns.append([og_to_col[og] for og in chosen])
            indicator_i = np.concatenate([np.asarray(columns) for columns in random_columns])
            indicator_j = np.repeat(np.arange(RANDOM_PROGRAMS), random_size)
            indicator = sp.csc_matrix((np.ones(len(indicator_i)), (indicator_i, indicator_j)), shape=(len(orthogroups), RANDOM_PROGRAMS))
            random_scores = (contrast_matrix @ indicator).toarray()
            for index in range(RANDOM_PROGRAMS):
                random_aucs[index] = roc_auc(y, random_scores[:, index])

            complete_metrics = metrics(y, complete_score)
            empirical_p = float((1 + np.sum(random_aucs >= complete_metrics["roc_auc"])) / (RANDOM_PROGRAMS + 1))
            positive_directions = int(np.sum(np.asarray(complete_replicate_aucs) > 0.5))
            passed = bool(
                frozen["role"] == "primary"
                and complete_metrics["roc_auc"] > 0.60
                and complete_metrics["standardized_mean_difference"] > 0.25
                and empirical_p <= 0.05
                and positive_directions == 3
            )
            summary_rows.append({
                "contrast_id": frozen["contrast_id"], "dataset": dataset,
                "panel": frozen["panel"], "coarse_label": frozen["coarse_label"],
                "role": frozen["role"], "positive_cells": int((y == 1).sum()),
                "negative_cells": int((y == 0).sum()), "mapped_target_orthogroups": len(target),
                "complete_program_roc_auc": complete_metrics["roc_auc"],
                "complete_program_average_precision": complete_metrics["average_precision"],
                "complete_program_standardized_mean_difference": complete_metrics["standardized_mean_difference"],
                "random_program_auc_median": float(np.median(random_aucs)),
                "random_program_auc_q95": float(np.quantile(random_aucs, 0.95)),
                "empirical_p_random_auc_ge_complete": empirical_p,
                "positive_replicate_directions": positive_directions,
                "required_positive_replicate_directions": 3,
                "frozen_contrast_pass": passed,
            })

    write_tsv(PARTIAL_DIR / "method_comparison_pooled.tsv", pooled_rows)
    write_tsv(PARTIAL_DIR / "method_comparison_by_replicate.tsv", replicate_rows)
    write_tsv(PARTIAL_DIR / "validation_summary.tsv", summary_rows)

    pooled_lookup = {(row["contrast_id"], row["method"]): row for row in pooled_rows}
    replicate_groups = defaultdict(list)
    for row in replicate_rows:
        replicate_groups[(row["contrast_id"], row["method"])].append(float(row["roc_auc"]))
    advantage_rows = []
    for summary in summary_rows:
        contrast = summary["contrast_id"]
        complete = pooled_lookup[(contrast, "complete_orthogroup_family_program")]
        single = pooled_lookup[(contrast, "target_species_single_gene_orthogroups")]
        delta = float(complete["roc_auc"]) - float(single["roc_auc"])
        complete_median = float(np.nanmedian(replicate_groups[(contrast, "complete_orthogroup_family_program")]))
        single_median = float(np.nanmedian(replicate_groups[(contrast, "target_species_single_gene_orthogroups")]))
        advantage_rows.append({
            "contrast_id": contrast,
            "complete_family_roc_auc": complete["roc_auc"],
            "single_gene_roc_auc": single["roc_auc"],
            "pooled_auc_delta_complete_minus_single": delta,
            "complete_family_replicate_auc_median": complete_median,
            "single_gene_replicate_auc_median": single_median,
            "frozen_method_advantage_supported": bool(
                int(single["candidate_orthogroups"]) > 0
                and delta >= 0.02 and complete_median > single_median
            ),
        })
    write_tsv(PARTIAL_DIR / "method_advantage.tsv", advantage_rows)

    primary = [row for row in summary_rows if row["role"] == "primary"]
    passing = [row for row in primary if row["frozen_contrast_pass"]]
    majority_reversal = [row for row in primary if int(row["positive_replicate_directions"]) < 2]
    panels = {row["panel"] for row in passing}
    if len(passing) >= 2 and {"leaf", "root"}.issubset(panels) and not majority_reversal:
        overall = "positive_external_replication"
    elif passing:
        overall = "mixed_external_support"
    else:
        overall = "external_nonreplication"

    audit = {
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "status": "independent_program_validation_v3_complete",
        "overall_interpretation": overall,
        "primary_contrasts": len(primary),
        "passing_primary_contrasts": len(passing),
        "passing_contrast_ids": [row["contrast_id"] for row in passing],
        "majority_reversal_contrast_ids": [row["contrast_id"] for row in majority_reversal],
        "top_k": 128,
        "random_programs": RANDOM_PROGRAMS,
        "random_seed": RANDOM_SEED,
        "freeze_sha256": sha256(META / "external_validation_v3_freeze_v1.tsv"),
        "protocol_sha256": sha256(REPORTS / "external_validation_v3_frozen_protocol_v1.md"),
        "execution_addendum_sha256": sha256(REPORTS / "external_validation_v3_execution_addendum_v1.json"),
        "safety": {
            "source_objects_modified": False,
            "freeze_modified_after_scoring": False,
            "model_training_started": False,
            "v1_or_v2_outputs_modified": False,
            "other_task_written": False,
        },
    }
    (PARTIAL_DIR / "audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    PARTIAL_DIR.rename(RUN_DIR)

    copies = {
        RUN_DIR / "validation_summary.tsv": META / "independent_program_validation_summary_v3.tsv",
        RUN_DIR / "method_comparison_pooled.tsv": META / "independent_program_method_comparison_pooled_v3.tsv",
        RUN_DIR / "method_comparison_by_replicate.tsv": META / "independent_program_method_comparison_by_replicate_v3.tsv",
        RUN_DIR / "method_advantage.tsv": META / "independent_program_method_advantage_v3.tsv",
        RUN_DIR / "audit.json": REPORTS / "independent_program_validation_v3_audit.json",
    }
    for source, target in copies.items():
        if target.exists():
            raise RuntimeError(f"refusing to overwrite {target}")
        target.write_bytes(source.read_bytes())
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
