#!/usr/bin/env python3
"""Build gene-level dominance evidence for priority conserved orthogroups.

The analysis reads frozen sampled cell indices and source H5AD files in backed mode.
It does not train a model and does not treat orthogroup recurrence as proof of
paralog substitution.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ALLOWED = {"exact_candidate", "normalized_candidate", "deterministic_transform_candidate"}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def atomic_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        frame.to_csv(temp, sep="\t", index=False)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        temp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def moments(matrix):
    if sparse.issparse(matrix):
        m = np.asarray(matrix.mean(axis=0)).ravel()
        q = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        nz = np.asarray((matrix != 0).mean(axis=0)).ravel()
    else:
        a = np.asarray(matrix, dtype=np.float64)
        m = np.nanmean(a, axis=0)
        q = np.nanmean(a * a, axis=0)
        nz = np.nanmean(a != 0, axis=0)
    return m, np.maximum(q - m * m, 0.0), nz


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-dir", type=Path, required=True)
    p.add_argument("--feature-dir", type=Path, required=True)
    p.add_argument("--gene-map", type=Path, required=True)
    p.add_argument("--orthogroup-map", type=Path, required=True)
    p.add_argument("--min-cells-per-side", type=int, default=20)
    p.add_argument("--version", default="v2")
    args = p.parse_args()
    metadata = args.project_dir / "metadata"

    evidence = read_tsv(metadata / f"conserved_program_function_evidence_{args.version}.tsv")
    priority = evidence[evidence["priority_tier"].eq("A")].copy()
    targets = {(r.panel, r.coarse_label, r.orthogroup_id) for r in priority.itertuples(index=False)}
    target_ogs = {r.orthogroup_id for r in priority.itertuples(index=False)}
    manifest = read_tsv(metadata / "dataset_manifest.tsv").set_index("dataset_id")
    hierarchy = read_tsv(metadata / "coarse_label_hierarchy_v1.tsv")
    hierarchy = hierarchy[hierarchy["accepted"].astype(str).eq("True")]
    label_to_coarse = {(r.panel, r.fine_source_label): r.coarse_label for r in hierarchy.itertuples(index=False)}

    obs_by_panel = {pnl: read_tsv(args.feature_dir / f"{pnl}_sampled_cells.tsv") for pnl in ("leaf", "root", "vascular")}
    selected_datasets = set(pd.concat(obs_by_panel.values())["dataset_id"])
    gene_maps: dict[str, dict[str, str]] = defaultdict(dict)
    with gzip.open(args.gene_map, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["dataset_id"] in selected_datasets and row["mapping_status"] in ALLOWED and row["candidate_reference_gene_id"]:
                gene_maps[row["dataset_id"]][row["source_gene_id"]] = row["candidate_reference_gene_id"]
    ref_to_og = {}
    family_size = defaultdict(int)
    with gzip.open(args.orthogroup_map, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["orthogroup_status"] == "assigned_candidate" and row["orthogroup_id"] in target_ogs:
                key = (row["species"], row["candidate_reference_gene_id"])
                ref_to_og[key] = row["orthogroup_id"]
                family_size[(row["species"], row["orthogroup_id"])] += 1

    result_rows = []
    dataset_checks = []
    for panel, obs in obs_by_panel.items():
        panel_targets = [(coarse, og) for pnl, coarse, og in targets if pnl == panel]
        for dataset_id, group in obs.groupby("dataset_id", sort=True):
            species = group["species"].iloc[0]
            path = Path(manifest.loc[dataset_id, "h5ad_path"])
            obj = ad.read_h5ad(path, backed="r")
            var_names = np.asarray(obj.var_names.astype(str))
            columns, refs, ogs = [], [], []
            source_map = gene_maps.get(dataset_id, {})
            for idx, source_gene in enumerate(var_names):
                ref = source_map.get(source_gene)
                og = ref_to_og.get((species, ref)) if ref else None
                if og in target_ogs:
                    columns.append(idx); refs.append(ref); ogs.append(og)
            row_indices = group["source_cell_index"].astype(int).to_numpy()
            coarse_labels = np.asarray([label_to_coarse[(panel, x)] for x in group["source_label"]])
            if columns:
                matrix = obj[row_indices, np.asarray(columns, dtype=int)].X
                matrix = matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix)
            else:
                matrix = sparse.csr_matrix((len(row_indices), 0))
            tests = 0
            for coarse, og in panel_targets:
                in_mask = coarse_labels == coarse
                out_mask = coarse_labels != coarse
                if in_mask.sum() < args.min_cells_per_side or out_mask.sum() < args.min_cells_per_side:
                    continue
                og_cols = [i for i, value in enumerate(ogs) if value == og]
                if not og_cols:
                    continue
                mi, vi, di = moments(matrix[in_mask][:, og_cols])
                mo, vo, do = moments(matrix[out_mask][:, og_cols])
                denom = np.sqrt((vi + vo) / 2.0 + 1e-12)
                smd = np.divide(mi - mo, denom, out=np.zeros_like(mi, dtype=float), where=denom > 0)
                for local_idx, col in enumerate(og_cols):
                    result_rows.append({
                        "panel": panel, "coarse_label": coarse, "orthogroup_id": og,
                        "species": species, "dataset_id": dataset_id, "reference_gene_id": refs[col],
                        "class_cells": int(in_mask.sum()), "other_cells": int(out_mask.sum()),
                        "mean_expression_class": float(mi[local_idx]), "mean_expression_other": float(mo[local_idx]),
                        "standardized_mean_difference": float(smd[local_idx]),
                        "nonzero_fraction_class": float(di[local_idx]), "nonzero_fraction_other": float(do[local_idx]),
                        "nonzero_fraction_difference": float(di[local_idx] - do[local_idx]),
                    })
                tests += 1
            dataset_checks.append({
                "panel": panel, "species": species, "dataset_id": dataset_id,
                "sampled_cells": len(row_indices), "priority_orthogroup_gene_columns": len(columns),
                "eligible_orthogroup_label_tests": tests,
            })
            if getattr(obj, "file", None) is not None:
                obj.file.close()

    detail = pd.DataFrame(result_rows)
    gene_rows = []
    for keys, group in detail.groupby(["panel", "coarse_label", "orthogroup_id", "species", "reference_gene_id"], sort=True):
        panel, coarse, og, species, gene = keys
        effects = group["standardized_mean_difference"].to_numpy()
        gene_rows.append({
            "panel": panel, "coarse_label": coarse, "orthogroup_id": og, "species": species,
            "reference_gene_id": gene, "species_orthogroup_gene_count": family_size[(species, og)],
            "datasets_tested": group["dataset_id"].nunique(),
            "median_standardized_mean_difference": float(np.median(effects)),
            "minimum_standardized_mean_difference": float(np.min(effects)),
            "positive_dataset_fraction": float(np.mean(effects > 0)),
            "median_nonzero_fraction_difference": float(np.median(group["nonzero_fraction_difference"])),
        })
    genes = pd.DataFrame(gene_rows)
    dominance_rows = []
    for keys, group in genes.groupby(["panel", "coarse_label", "orthogroup_id", "species"], sort=True):
        panel, coarse, og, species = keys
        ranked = group.sort_values(["median_standardized_mean_difference", "positive_dataset_fraction", "reference_gene_id"], ascending=[False, False, True])
        positive = ranked[ranked["median_standardized_mean_difference"] > 0]
        dominant = ranked.iloc[0]
        positive_sum = float(positive["median_standardized_mean_difference"].sum())
        dominance_rows.append({
            "panel": panel, "coarse_label": coarse, "orthogroup_id": og, "species": species,
            "species_orthogroup_gene_count": int(dominant["species_orthogroup_gene_count"]),
            "genes_tested": len(group), "positive_marker_genes": len(positive),
            "dominant_reference_gene_id": dominant["reference_gene_id"],
            "dominant_median_standardized_mean_difference": float(dominant["median_standardized_mean_difference"]),
            "dominant_positive_dataset_fraction": float(dominant["positive_dataset_fraction"]),
            "dominant_effect_share_among_positive_genes": float(max(dominant["median_standardized_mean_difference"], 0) / positive_sum) if positive_sum > 0 else 0.0,
            "within_species_pattern": "single_dominant" if len(positive) == 1 or (positive_sum > 0 and max(dominant["median_standardized_mean_difference"], 0) / positive_sum >= 0.60) else "distributed_across_paralogs",
        })
    dominance = pd.DataFrame(dominance_rows)
    og_rows = []
    for keys, group in dominance.groupby(["panel", "coarse_label", "orthogroup_id"], sort=True):
        panel, coarse, og = keys
        represented = group[group["dominant_median_standardized_mean_difference"] > 0]
        multigene = represented[represented["species_orthogroup_gene_count"] > 1]
        og_rows.append({
            "panel": panel, "coarse_label": coarse, "orthogroup_id": og,
            "species_with_positive_gene_level_marker": represented["species"].nunique(),
            "species_with_multi_gene_family_and_positive_driver": multigene["species"].nunique(),
            "species_with_distributed_paralog_signal": int((represented["within_species_pattern"] == "distributed_across_paralogs").sum()),
            "dominant_genes_by_species": " | ".join(f"{r.species}: {r.dominant_reference_gene_id}" for r in represented.itertuples(index=False)),
            "gene_level_followup_status": "ready_for_gene_tree_subclade_mapping" if multigene["species"].nunique() >= 2 else "limited_multi_gene_support",
            "paralog_substitution_claim": "not_established",
        })
    og_summary = pd.DataFrame(og_rows)

    detail_path = metadata / f"priority_paralog_dataset_effects_{args.version}.tsv"
    gene_path = metadata / f"priority_paralog_gene_effects_{args.version}.tsv"
    dominance_path = metadata / f"priority_paralog_species_dominance_{args.version}.tsv"
    summary_path = metadata / f"priority_paralog_orthogroup_summary_{args.version}.tsv"
    checks_path = metadata / f"priority_paralog_dataset_checks_{args.version}.tsv"
    atomic_tsv(detail_path, detail)
    atomic_tsv(gene_path, genes)
    atomic_tsv(dominance_path, dominance)
    atomic_tsv(summary_path, og_summary)
    atomic_tsv(checks_path, pd.DataFrame(dataset_checks))
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "gene_level_dominance_complete_substitution_not_established",
        "priority_orthogroups": len(priority), "dataset_gene_effect_rows": len(detail),
        "aggregated_gene_rows": len(genes), "species_dominance_rows": len(dominance),
        "orthogroups_ready_for_gene_tree_mapping": int(og_summary["gene_level_followup_status"].eq("ready_for_gene_tree_subclade_mapping").sum()),
        "safety": {"source_h5ad_mode": "read_only_backed", "source_h5ad_modified": False,
                   "model_training_started": False, "external_gene_queries": False},
        "outputs": {"dataset_effects": str(detail_path), "gene_effects": str(gene_path),
                    "species_dominance": str(dominance_path), "orthogroup_summary": str(summary_path),
                    "dataset_checks": str(checks_path)},
        "limitations": [
            "Standardized effects are calculated within each dataset and summarized across datasets; they are not formal differential-expression statistics.",
            "Dominant-gene differences across species do not establish paralog substitution without orthogroup gene-tree subclade mapping.",
            "Detection fractions are descriptive because source expression transformations vary among datasets.",
        ],
    }
    atomic_json(args.project_dir / "reports" / f"priority_paralog_evidence_{args.version}_audit.json", audit)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
