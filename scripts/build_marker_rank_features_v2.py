#!/usr/bin/env python3
"""Build corrected rank-orthogroup features for biological marker analysis only.

This script reads source H5AD files in backed/read-only mode and does not train or
evaluate a model. V2 corrects the V1 mapping-status filter by including curated
deterministic transform candidates.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import tempfile
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse
from sklearn.preprocessing import normalize


ALLOWED_MAPPING_STATUSES = {
    "exact_candidate",
    "normalized_candidate",
    "deterministic_transform_candidate",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".partial", dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        temp.write_text(text, encoding="utf-8")
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def write_tsv(path: Path, rows: list[dict]) -> None:
    from io import StringIO
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, buffer.getvalue())


def stable_seed(text: str, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}|{text}".encode()).digest()
    return int.from_bytes(digest[:8], "little") % (2**32)


def load_maps(gene_path: Path, og_path: Path, dataset_ids: set[str]):
    gene_maps: dict[str, dict[str, str]] = defaultdict(dict)
    status_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with gzip.open(gene_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["dataset_id"] not in dataset_ids:
                continue
            status_counts[row["dataset_id"]][row["mapping_status"]] += 1
            if row["mapping_status"] in ALLOWED_MAPPING_STATUSES and row["candidate_reference_gene_id"]:
                gene_maps[row["dataset_id"]][row["source_gene_id"]] = row["candidate_reference_gene_id"]
    ref_to_og: dict[tuple[str, str], str] = {}
    with gzip.open(og_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["orthogroup_status"] == "assigned_candidate" and row["orthogroup_id"]:
                ref_to_og[(row["species"], row["candidate_reference_gene_id"])] = row["orthogroup_id"]
    return gene_maps, ref_to_og, status_counts


def rank_encode_batch(matrix, feature_ogs: np.ndarray, og_to_col: dict[str, int], top_k: int, nonnegative: bool):
    mat = matrix.tocsr() if sparse.issparse(matrix) else np.asarray(matrix)
    out_r, out_c, out_v = [], [], []
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
        take = min(top_k, len(idx))
        if take < len(idx):
            local = np.argpartition(vals, -take)[-take:]
            idx, vals = idx[local], vals[local]
        order = np.argsort(-vals, kind="stable")
        best: dict[int, float] = {}
        for rank, feature_idx in enumerate(idx[order], start=1):
            col = og_to_col[str(feature_ogs[feature_idx])]
            score = 1.0 / math.log2(rank + 1.0)
            best[col] = max(score, best.get(col, 0.0))
        for col, score in best.items():
            out_r.append(i); out_c.append(col); out_v.append(score)
    return sparse.csr_matrix((out_v, (out_r, out_c)), shape=(mat.shape[0], len(og_to_col)), dtype=np.float32)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-dir", type=Path, required=True)
    p.add_argument("--gene-map", type=Path, required=True)
    p.add_argument("--orthogroup-map", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--max-cells-per-dataset-label", type=int, default=750)
    p.add_argument("--top-k", type=int, default=128)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--seed", type=int, default=20260911)
    args = p.parse_args()
    started = time.time()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata = args.project_dir / "metadata"

    accepted_labels: dict[str, set[str]] = defaultdict(set)
    for row in read_tsv(metadata / "panel_label_freeze_v1.tsv"):
        if row["accepted"] == "True":
            accepted_labels[row["panel"]].add(row["source_label"])
    dataset_rows = [r for r in read_tsv(metadata / "panel_dataset_freeze_v1.tsv") if r["accepted"] == "True"]
    manifest = {r["dataset_id"]: r for r in read_tsv(metadata / "dataset_manifest.tsv")}
    dataset_ids = {r["dataset_id"] for r in dataset_rows}
    gene_maps, ref_to_og, status_counts = load_maps(args.gene_map, args.orthogroup_map, dataset_ids)
    all_ogs = sorted(set(ref_to_og.values()))
    og_to_col = {og: i for i, og in enumerate(all_ogs)}
    dataset_panels: dict[str, set[str]] = defaultdict(set)
    study_lookup = {}
    for row in dataset_rows:
        dataset_panels[row["dataset_id"]].add(row["panel"])
        study_lookup[(row["dataset_id"], row["panel"])] = row["study_group_id"]

    blocks: dict[str, list[sparse.csr_matrix]] = defaultdict(list)
    obs_rows: dict[str, list[dict]] = defaultdict(list)
    audit_rows = []
    for dataset_id in sorted(dataset_ids):
        info = manifest[dataset_id]
        species = info["proposed_species"]
        obj = ad.read_h5ad(Path(info["h5ad_path"]), backed="r")
        var_names = np.asarray(obj.var_names.astype(str))
        mapped_indices, mapped_ogs = [], []
        for idx, source_gene in enumerate(var_names):
            ref = gene_maps.get(dataset_id, {}).get(source_gene)
            og = ref_to_og.get((species, ref)) if ref else None
            if og:
                mapped_indices.append(idx); mapped_ogs.append(og)
        mapped_indices = np.asarray(mapped_indices, dtype=np.int64)
        feature_ogs = np.asarray(mapped_ogs, dtype=object)
        cell_labels = obj.obs["cell_type_original"].astype(str).to_numpy()
        for panel in sorted(dataset_panels[dataset_id]):
            selected_parts = []
            for label in sorted(accepted_labels[panel]):
                indices = np.flatnonzero(cell_labels == label)
                if len(indices) > args.max_cells_per_dataset_label:
                    rng = np.random.default_rng(stable_seed(f"{dataset_id}|{panel}|{label}", args.seed))
                    indices = np.sort(rng.choice(indices, args.max_cells_per_dataset_label, replace=False))
                if len(indices):
                    selected_parts.append((label, indices))
            if not selected_parts:
                continue
            selected = np.sort(np.concatenate([part[1] for part in selected_parts]))
            label_lookup = {int(i): label for label, indices in selected_parts for i in indices}
            encoded = []
            for offset in range(0, len(selected), args.batch_size):
                batch = selected[offset:offset + args.batch_size]
                encoded.append(rank_encode_batch(obj[batch, mapped_indices].X, feature_ogs, og_to_col,
                                                 args.top_k, info.get("expression_class") != "scaled_continuous"))
            feature_block = sparse.vstack(encoded, format="csr")
            blocks[panel].append(feature_block)
            for index in selected:
                obs_rows[panel].append({
                    "dataset_id": dataset_id,
                    "study_group_id": study_lookup[(dataset_id, panel)],
                    "species": species,
                    "panel": panel,
                    "source_label": label_lookup[int(index)],
                    "source_cell_index": int(index),
                })
            audit_rows.append({
                "dataset_id": dataset_id, "species": species, "panel": panel,
                "mapped_features": len(mapped_indices), "sampled_cells": len(selected),
                "nonzero_cells": int(np.count_nonzero(feature_block.getnnz(axis=1))),
                "rank_feature_nnz": int(feature_block.nnz),
                "exact_candidate_rows": status_counts[dataset_id].get("exact_candidate", 0),
                "deterministic_transform_candidate_rows": status_counts[dataset_id].get("deterministic_transform_candidate", 0),
            })
        if getattr(obj, "file", None) is not None:
            obj.file.close()

    panel_summary = {}
    for panel in sorted(blocks):
        x = normalize(sparse.vstack(blocks[panel], format="csr"), norm="l2", axis=1, copy=False)
        obs = obs_rows[panel]
        if x.shape[0] != len(obs):
            raise RuntimeError(f"row mismatch for {panel}")
        feature_path = args.output_dir / f"{panel}_rank_orthogroup_features.npz"
        sparse.save_npz(feature_path.with_suffix(".npz.partial"), x)
        os.replace(feature_path.with_suffix(".npz.partial.npz"), feature_path)
        obs_path = args.output_dir / f"{panel}_sampled_cells.tsv"
        write_tsv(obs_path, obs)
        panel_summary[panel] = {
            "cells": x.shape[0], "orthogroup_features": x.shape[1], "nnz": int(x.nnz),
            "nonzero_cells": int(np.count_nonzero(x.getnnz(axis=1))),
            "feature_path": str(feature_path), "feature_sha256": sha256(feature_path),
            "obs_path": str(obs_path), "obs_sha256": sha256(obs_path),
        }
    audit_path = args.output_dir / "dataset_feature_audit.tsv"
    write_tsv(audit_path, audit_rows)
    zero_datasets = [r["dataset_id"] for r in audit_rows if r["rank_feature_nnz"] == 0]
    payload = {
        "schema_version": "2.0.0", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete" if not zero_datasets else "failed_zero_feature_dataset",
        "purpose": "marker_and_conserved_program_analysis_only_no_model_training",
        "mapping_statuses_included": sorted(ALLOWED_MAPPING_STATUSES),
        "panels": panel_summary, "datasets": len(audit_rows), "zero_feature_datasets": zero_datasets,
        "dataset_audit": {"path": str(audit_path), "sha256": sha256(audit_path)},
        "runtime_seconds": time.time() - started,
        "safety": {"source_h5ad_mode": "read_only_backed", "source_h5ad_files_modified": 0,
                   "model_training_started": False, "v1_outputs_overwritten": False},
    }
    atomic_text(args.output_dir / "summary.json", json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
