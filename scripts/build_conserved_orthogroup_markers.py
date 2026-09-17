#!/usr/bin/env python3
"""Build exploratory coarse-label orthogroup marker candidates from frozen rank features."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(temp_name)
    try:
        frame.to_csv(temp, sep="\t", index=False)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(temp_name)
    try:
        temp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument("--feature-dir", type=Path, required=True)
    parser.add_argument("--orthogroup-map", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--min-class-cells", type=int, default=100)
    parser.add_argument("--min-other-cells", type=int, default=100)
    parser.add_argument("--min-species", type=int, default=3)
    parser.add_argument("--version", default="v1")
    args = parser.parse_args()

    hierarchy_path = args.project_dir / "metadata" / "coarse_label_hierarchy_v1.tsv"
    hierarchy = pd.read_csv(hierarchy_path, sep="\t")
    hierarchy = hierarchy[hierarchy["accepted"].astype(str).eq("True")]
    label_to_coarse = {
        (row.panel, row.fine_source_label): row.coarse_label
        for row in hierarchy.itertuples(index=False)
    }

    with gzip.open(args.orthogroup_map, "rt", encoding="utf-8", newline="") as handle:
        mapping = pd.DataFrame(csv.DictReader(handle, delimiter="\t"))
    mapping = mapping[
        mapping["orthogroup_status"].eq("assigned_candidate") & mapping["orthogroup_id"].ne("")
    ].copy()
    all_ogs = sorted(mapping["orthogroup_id"].unique())
    gene_lookup: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in mapping.itertuples(index=False):
        key = (row.species, row.orthogroup_id)
        if row.candidate_reference_gene_id not in gene_lookup[key]:
            gene_lookup[key].append(row.candidate_reference_gene_id)

    detail_rows: list[dict] = []
    eligible_tests: list[dict] = []
    representation_rows: list[dict] = []
    for panel in ("leaf", "root", "vascular"):
        feature_path = args.feature_dir / f"{panel}_rank_orthogroup_features.npz"
        obs_path = args.feature_dir / f"{panel}_sampled_cells.tsv"
        x = sparse.load_npz(feature_path).tocsr()
        obs = pd.read_csv(obs_path, sep="\t")
        if x.shape[0] != len(obs) or x.shape[1] != len(all_ogs):
            raise RuntimeError(f"shape mismatch for {panel}: {x.shape}, obs={len(obs)}, ogs={len(all_ogs)}")
        obs["coarse_label"] = [label_to_coarse.get((panel, label), "") for label in obs["source_label"]]
        if (obs["coarse_label"] == "").any():
            missing = sorted(obs.loc[obs["coarse_label"] == "", "source_label"].unique())
            raise RuntimeError(f"unmapped accepted labels in {panel}: {missing}")

        for species in sorted(obs["species"].unique()):
            species_mask = obs["species"].eq(species).to_numpy()
            species_nnz = int(x[species_mask].nnz)
            nonzero_cells = int(np.count_nonzero(x[species_mask].getnnz(axis=1)))
            representation_rows.append({
                "panel": panel,
                "species": species,
                "cells": int(species_mask.sum()),
                "nonzero_feature_cells": nonzero_cells,
                "nonzero_cell_fraction": nonzero_cells / int(species_mask.sum()),
                "matrix_nnz": species_nnz,
                "usable_for_marker_analysis": nonzero_cells == int(species_mask.sum()) and species_nnz > 0,
            })
            if species_nnz == 0:
                continue
            for coarse in sorted(obs.loc[species_mask, "coarse_label"].unique()):
                in_mask = species_mask & obs["coarse_label"].eq(coarse).to_numpy()
                out_mask = species_mask & ~obs["coarse_label"].eq(coarse).to_numpy()
                n_in, n_out = int(in_mask.sum()), int(out_mask.sum())
                if n_in < args.min_class_cells or n_out < args.min_other_cells:
                    continue
                effect = np.asarray(x[in_mask].mean(axis=0)).ravel() - np.asarray(x[out_mask].mean(axis=0)).ravel()
                positive = np.flatnonzero(effect > 0)
                order = positive[np.argsort(-effect[positive], kind="stable")][: args.top_n]
                eligible_tests.append({"panel": panel, "species": species, "coarse_label": coarse,
                                       "n_in": n_in, "n_out": n_out, "positive_features": int(len(positive))})
                for rank, idx in enumerate(order, start=1):
                    genes = sorted(gene_lookup.get((species, all_ogs[idx]), []))
                    detail_rows.append({
                        "panel": panel,
                        "coarse_label": coarse,
                        "species": species,
                        "orthogroup_id": all_ogs[idx],
                        "within_species_rank": rank,
                        "mean_rank_score_difference": float(effect[idx]),
                        "class_cells": n_in,
                        "other_cells": n_out,
                        "representative_genes": ";".join(genes[:3]),
                        "genes_in_species_orthogroup": len(genes),
                        "multi_gene_family_in_species": len(genes) > 1,
                    })

    detail = pd.DataFrame(detail_rows)
    if detail.empty:
        raise RuntimeError("no eligible species/coarse-label marker comparisons")

    conserved_rows: list[dict] = []
    grouped = detail.groupby(["panel", "coarse_label", "orthogroup_id"], sort=True)
    for (panel, coarse, og), group in grouped:
        if group["species"].nunique() < args.min_species:
            continue
        species = sorted(group["species"].unique())
        gene_text = []
        for sp in species:
            genes = sorted(gene_lookup.get((sp, og), []))[:3]
            gene_text.append(f"{sp}: {','.join(genes) if genes else 'NA'}")
        conserved_rows.append({
            "panel": panel,
            "coarse_label": coarse,
            "orthogroup_id": og,
            "species_count": len(species),
            "species": "; ".join(species),
            "median_within_species_rank": float(group["within_species_rank"].median()),
            "mean_rank_score_difference": float(group["mean_rank_score_difference"].mean()),
            "minimum_rank_score_difference": float(group["mean_rank_score_difference"].min()),
            "species_with_multi_gene_family": int(group["multi_gene_family_in_species"].sum()),
            "representative_genes_by_species": " | ".join(gene_text),
        })
    conserved = pd.DataFrame(conserved_rows)
    if not conserved.empty:
        specificity = conserved.groupby(["panel", "orthogroup_id"])["coarse_label"].transform("nunique")
        conserved["coarse_labels_in_panel"] = specificity.astype(int)
        conserved["low_coarse_specificity_flag"] = specificity.gt(1)
        conserved = conserved.sort_values(
            ["panel", "coarse_label", "species_count", "median_within_species_rank", "orthogroup_id"],
            ascending=[True, True, False, True, True],
        )
    else:
        conserved = pd.DataFrame(columns=[
            "panel", "coarse_label", "orthogroup_id", "species_count", "species",
            "median_within_species_rank", "mean_rank_score_difference", "minimum_rank_score_difference",
            "species_with_multi_gene_family", "representative_genes_by_species",
            "coarse_labels_in_panel", "low_coarse_specificity_flag",
        ])

    detail_path = args.project_dir / "metadata" / f"coarse_marker_orthogroups_{args.version}.tsv"
    conserved_path = args.project_dir / "metadata" / f"conserved_program_candidates_{args.version}.tsv"
    representation_path = args.project_dir / "metadata" / f"marker_feature_representation_audit_{args.version}.tsv"
    audit_path = args.project_dir / "reports" / f"conserved_program_candidates_{args.version}_audit.json"
    atomic_tsv(detail_path, detail)
    atomic_tsv(conserved_path, conserved)
    representation = pd.DataFrame(representation_rows)
    atomic_tsv(representation_path, representation)
    unusable = representation.loc[~representation["usable_for_marker_analysis"], ["panel", "species"]]
    blocked_species = [f"{r.panel}|{r.species}" for r in unusable.itertuples(index=False)]
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "blocked_by_zero_feature_species" if blocked_species else "exploratory_candidates_not_independently_validated",
        "method": {
            "comparison": "coarse class versus all other coarse classes within the same species and panel",
            "feature": "L2-normalized rank-encoded orthogroup score from frozen baseline",
            "ranking": "positive mean score difference, descending",
            "top_n_per_species_class": args.top_n,
            "minimum_class_cells": args.min_class_cells,
            "minimum_other_cells": args.min_other_cells,
            "minimum_species_recurrence": args.min_species,
        },
        "eligible_species_class_tests": len(eligible_tests),
        "species_marker_rows": len(detail),
        "conserved_candidate_rows": len(conserved),
        "candidate_counts_by_panel": conserved.groupby("panel").size().astype(int).to_dict() if len(conserved) else {},
        "candidate_counts_by_panel_and_class": {
            f"{p}|{c}": int(n) for (p, c), n in conserved.groupby(["panel", "coarse_label"]).size().items()
        } if len(conserved) else {},
        "inputs": {
            "hierarchy": {"path": str(hierarchy_path), "sha256": sha256(hierarchy_path)},
            "orthogroup_map": {"path": str(args.orthogroup_map), "sha256": sha256(args.orthogroup_map)},
            "feature_dir": str(args.feature_dir),
        },
        "feature_representation": {
            "panel_species_combinations": len(representation),
            "usable_panel_species_combinations": int(representation["usable_for_marker_analysis"].sum()),
            "zero_feature_panel_species": blocked_species,
        },
        "outputs": {"detail": str(detail_path), "conserved": str(conserved_path),
                    "representation_audit": str(representation_path)},
        "caveats": [
            "Candidates are derived from the same frozen observational datasets and are not independent validation.",
            "Orthogroup recurrence does not establish conserved function or causal mechanism.",
            "Multi-gene orthogroups may reflect paralog substitution and require gene-level follow-up.",
            "Rank features were designed for robust cross-dataset comparison, not differential-expression inference.",
            "Zero-feature species are excluded from marker comparisons; a zero candidate count is not biological evidence of non-conservation.",
        ],
        "eligible_tests": eligible_tests,
    }
    atomic_json(audit_path, audit)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
