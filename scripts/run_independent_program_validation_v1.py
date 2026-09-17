#!/usr/bin/env python3
"""Independent validation of frozen conserved programs without model training.

The two validation datasets and target labels are fixed in
independent_validation_freeze_v1.tsv. This script reads their H5AD files in
backed/read-only mode, builds the same top-128 rank/orthogroup representation,
and compares complete family programs with restricted subsets and random
orthogroup controls. It refuses to overwrite an existing run.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import random
import shutil
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import normalize


PROJECT = Path("/workspace/projects/PhyloOpenCell")
META = PROJECT / "metadata"
REPORTS = PROJECT / "reports"
RUN = Path("/data/r/PhyloOpenCell")
RUN = Path("/data/runs/PhyloOpenCell/independent_program_validation_v1")
PARTIAL = RUN.with_name(RUN.name + ".partial")
GENE_MAP = Path("/data/datasets/PhyloOpenCell/mappings/gene_mapping_candidates_v1.tsv.gz")
OG_MAP = Path("/data/datasets/PhyloOpenCell/mappings/gene_orthogroup_candidates_v1.tsv.gz")
ORTHOGROUPS = Path("/data/runs/PhyloOpenCell/orthofinder_full7_v1/results/Results_Sep10/Orthogroups/Orthogroups.tsv")
ALLOWED = {"exact_candidate", "normalized_candidate", "deterministic_transform_candidate"}
SEED = 20260911
TOP_K = 128
MAX_PER_CLASS = 3000
RANDOM_SETS = 500
TARGETS = {
    ("vascular", "Phloem"): ("phloem_lineage", "xylem_lineage"),
    ("leaf", "Leaf guard cell"): ("epidermal_system", "photosynthetic_ground_tissue"),
}


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stable_seed(text):
    return int.from_bytes(hashlib.sha256(f"{SEED}|{text}".encode()).digest()[:8], "little") % (2**32)


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_maps(dataset_ids):
    gene_maps = defaultdict(dict)
    status_counts = defaultdict(Counter)
    with gzip.open(GENE_MAP, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dataset = row["dataset_id"]
            if dataset not in dataset_ids:
                continue
            status_counts[dataset][row["mapping_status"]] += 1
            if row["mapping_status"] in ALLOWED and row["candidate_reference_gene_id"]:
                gene_maps[dataset][row["source_gene_id"]] = row["candidate_reference_gene_id"]
    ref_to_og = {}
    with gzip.open(OG_MAP, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["orthogroup_status"] == "assigned_candidate" and row["orthogroup_id"]:
                ref_to_og[(row["species"], row["candidate_reference_gene_id"])] = row["orthogroup_id"]
    return gene_maps, ref_to_og, status_counts


def orthogroup_gene_counts():
    result = {}
    with ORTHOGROUPS.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = [c for c in reader.fieldnames if c != "Orthogroup"]
        for row in reader:
            counts = {}
            for column in columns:
                species = column.split(".representative_proteins", 1)[0].replace("_", " ").capitalize()
                counts[species] = len([x for x in row[column].split(",") if x.strip()])
            result[row["Orthogroup"]] = counts
    return result


def encode(matrix, feature_ogs, og_to_col, nonnegative):
    mat = matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix)
    rr, cc, vv = [], [], []
    for i in range(mat.shape[0]):
        if sparse.issparse(mat):
            start, end = mat.indptr[i], mat.indptr[i + 1]
            idx, vals = mat.indices[start:end], mat.data[start:end]
            keep = np.isfinite(vals) & (vals > 0)
            idx, vals = idx[keep], vals[keep]
        else:
            vals_all = mat[i]
            keep = np.isfinite(vals_all) & ((vals_all > 0) if nonnegative else np.ones(vals_all.shape, dtype=bool))
            idx = np.flatnonzero(keep)
            vals = vals_all[idx]
        if not len(idx):
            continue
        take = min(TOP_K, len(idx))
        if take < len(idx):
            chosen = np.argpartition(vals, -take)[-take:]
            idx, vals = idx[chosen], vals[chosen]
        order = np.argsort(-vals, kind="stable")
        best = {}
        for rank, feature_index in enumerate(idx[order], start=1):
            col = og_to_col[str(feature_ogs[feature_index])]
            best[col] = max(best.get(col, 0.0), 1.0 / math.log2(rank + 1.0))
        for col, value in best.items():
            rr.append(i); cc.append(col); vv.append(value)
    return sparse.csr_matrix((vv, (rr, cc)), shape=(mat.shape[0], len(og_to_col)), dtype=np.float32)


def score(x, columns):
    if not columns:
        return np.zeros(x.shape[0], dtype=float)
    return np.asarray(x[:, columns].sum(axis=1)).ravel()


def effect(y, values):
    pos, neg = values[y == 1], values[y == 0]
    pooled = math.sqrt((float(np.var(pos, ddof=1)) + float(np.var(neg, ddof=1))) / 2) if len(pos) > 1 and len(neg) > 1 else 0
    return (float(np.mean(pos)) - float(np.mean(neg))) / pooled if pooled > 0 else 0.0


def metrics(y, values):
    return {
        "roc_auc": float(roc_auc_score(y, values)),
        "average_precision": float(average_precision_score(y, values)),
        "standardized_mean_difference": effect(y, values),
        "mean_positive": float(np.mean(values[y == 1])),
        "mean_negative": float(np.mean(values[y == 0])),
    }


def main():
    if RUN.exists() or PARTIAL.exists():
        raise RuntimeError(f"refusing to overwrite {RUN} or {PARTIAL}")
    PARTIAL.mkdir(parents=True)
    freeze = read_tsv(META / "independent_validation_freeze_v1.tsv")
    candidates = read_tsv(META / "conserved_program_candidates_v2.tsv")
    dataset_ids = {r["dataset_id"] for r in freeze}
    gene_maps, ref_to_og, mapping_counts = load_maps(dataset_ids)
    counts_by_og = orthogroup_gene_counts()
    all_ogs = sorted(set(ref_to_og.values()))
    og_to_col = {og: i for i, og in enumerate(all_ogs)}
    summary_rows, cell_rows, method_rows = [], [], []

    for frozen in freeze:
        dataset = frozen["dataset_id"]
        species = frozen["species"]
        panel = frozen["panel"]
        target_label = frozen["validation_label"]
        target_coarse, off_coarse = TARGETS[(panel, target_label)]
        target_candidates = [r for r in candidates if r["panel"] == panel and r["coarse_label"] == target_coarse]
        off_candidates = [r for r in candidates if r["panel"] == panel and r["coarse_label"] == off_coarse]
        target_ogs = [r["orthogroup_id"] for r in target_candidates]
        off_ogs = [r["orthogroup_id"] for r in off_candidates]
        single_target = [og for og in target_ogs if counts_by_og.get(og, {}).get(species, 0) == 1]
        expanded_target = [og for og in target_ogs if counts_by_og.get(og, {}).get(species, 0) > 1]
        strict_single = []
        for candidate in target_candidates:
            supporting = [x.strip() for x in candidate["species"].split(";")]
            og = candidate["orthogroup_id"]
            if counts_by_og.get(og, {}).get(species, 0) == 1 and all(counts_by_og.get(og, {}).get(s, 0) == 1 for s in supporting):
                strict_single.append(og)

        obj = ad.read_h5ad(Path(frozen["h5ad_path"]), backed="r")
        labels = obj.obs["cell_type_original"].astype(str).to_numpy()
        positive = np.flatnonzero(labels == target_label)
        unknown_terms = {"unknown", "nan", "unannotated", "unassigned", "none", "not available"}
        negative = np.flatnonzero(np.array([(x != target_label and x.strip().lower() not in unknown_terms) for x in labels]))
        rng = np.random.default_rng(stable_seed(dataset))
        if len(positive) > MAX_PER_CLASS:
            positive = np.sort(rng.choice(positive, MAX_PER_CLASS, replace=False))
        if len(negative) > MAX_PER_CLASS:
            negative = np.sort(rng.choice(negative, MAX_PER_CLASS, replace=False))
        selected = np.concatenate([positive, negative])
        y = np.concatenate([np.ones(len(positive), dtype=int), np.zeros(len(negative), dtype=int)])
        var_names = np.asarray(obj.var_names.astype(str))
        mapped_indices, mapped_ogs = [], []
        for index, source_gene in enumerate(var_names):
            reference = gene_maps[dataset].get(source_gene)
            og = ref_to_og.get((species, reference)) if reference else None
            if og:
                mapped_indices.append(index); mapped_ogs.append(og)
        mapped_indices = np.asarray(mapped_indices, dtype=np.int64)
        mapped_ogs = np.asarray(mapped_ogs, dtype=object)
        blocks = []
        for offset in range(0, len(selected), 256):
            batch = selected[offset:offset + 256]
            blocks.append(encode(obj[batch, mapped_indices].X, mapped_ogs, og_to_col, frozen["expression_class"] != "scaled_continuous"))
        x = normalize(sparse.vstack(blocks, format="csr"), norm="l2", axis=1, copy=False)
        if getattr(obj, "file", None) is not None:
            obj.file.close()

        methods = {
            "complete_orthogroup_family_program": target_ogs,
            "target_species_single_gene_orthogroups": single_target,
            "target_species_expanded_family_orthogroups": expanded_target,
            "strict_single_copy_across_supporting_species": strict_single,
            "off_target_program_negative_control": off_ogs,
        }
        method_scores = {}
        for method, ogs in methods.items():
            columns = [og_to_col[og] for og in ogs if og in og_to_col]
            values = score(x, columns)
            method_scores[method] = values
            result = metrics(y, values) if columns else {k: None for k in ["roc_auc", "average_precision", "standardized_mean_difference", "mean_positive", "mean_negative"]}
            method_rows.append({
                "dataset_id": dataset, "species": species, "panel": panel, "validation_label": target_label,
                "target_coarse_label": target_coarse, "method": method,
                "candidate_orthogroups": len(ogs), "mapped_candidate_orthogroups": len(columns),
                **result,
            })

        full_values = method_scores["complete_orthogroup_family_program"]
        random_pool = [og for og in all_ogs if og not in set(target_ogs) | set(off_ogs) and og in og_to_col]
        random_aucs = []
        random_rng = random.Random(stable_seed(dataset + "|random_programs"))
        for _ in range(RANDOM_SETS):
            picked = random_rng.sample(random_pool, min(len(target_ogs), len(random_pool)))
            random_aucs.append(float(roc_auc_score(y, score(x, [og_to_col[og] for og in picked]))))
        full_metrics = metrics(y, full_values)
        empirical_p = (1 + sum(v >= full_metrics["roc_auc"] for v in random_aucs)) / (RANDOM_SETS + 1)
        for row_index, cell_index in enumerate(selected):
            cell_rows.append({
                "dataset_id": dataset, "source_cell_index": int(cell_index), "source_label": labels[cell_index],
                "target": int(y[row_index]), "complete_program_score": float(full_values[row_index]),
                "single_gene_og_program_score": float(method_scores["target_species_single_gene_orthogroups"][row_index]),
                "expanded_family_program_score": float(method_scores["target_species_expanded_family_orthogroups"][row_index]),
                "off_target_program_score": float(method_scores["off_target_program_negative_control"][row_index]),
            })
        summary_rows.append({
            "dataset_id": dataset, "species": species, "panel": panel, "validation_label": target_label,
            "target_coarse_label": target_coarse, "positive_cells_sampled": len(positive), "negative_cells_sampled": len(negative),
            "mapped_gene_features": len(mapped_indices), "nonzero_encoded_cells": int(np.count_nonzero(x.getnnz(axis=1))),
            "target_candidate_orthogroups": len(target_ogs), "target_single_gene_orthogroups": len(single_target),
            "target_expanded_family_orthogroups": len(expanded_target), "strict_single_copy_orthogroups": len(strict_single),
            "complete_program_roc_auc": full_metrics["roc_auc"], "complete_program_average_precision": full_metrics["average_precision"],
            "complete_program_standardized_mean_difference": full_metrics["standardized_mean_difference"],
            "random_program_auc_median": statistics.median(random_aucs),
            "random_program_auc_q95": float(np.quantile(random_aucs, .95)),
            "empirical_p_random_program_auc_ge_complete": empirical_p,
            "mapping_status_counts": json.dumps(mapping_counts[dataset], sort_keys=True),
        })

    write_tsv(PARTIAL / "validation_summary.tsv", summary_rows)
    write_tsv(PARTIAL / "method_comparison.tsv", method_rows)
    write_tsv(PARTIAL / "sampled_cell_scores.tsv", cell_rows)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "independent_program_validation_complete",
        "datasets": len(summary_rows), "top_k": TOP_K, "max_cells_per_class": MAX_PER_CLASS,
        "random_control_sets": RANDOM_SETS, "summary": summary_rows,
        "interpretation_boundary": "Datasets and target labels were frozen before expression access. Cell-level metrics are descriptive because cells are not independent biological replicates. This validates program transfer, not causal function or paralog substitution.",
        "safety": {"source_h5ad_mode": "read_only_backed", "source_h5ad_modified": False, "model_training_started": False, "other_task_written": False, "discovery_outputs_modified": False},
    }
    (PARTIAL / "audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(PARTIAL, RUN)
    report = REPORTS / "independent_program_validation_v1_audit.json"
    report.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    shutil.copy2(RUN / "validation_summary.tsv", META / "independent_program_validation_summary_v1.tsv")
    shutil.copy2(RUN / "method_comparison.tsv", META / "independent_program_method_comparison_v1.tsv")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
