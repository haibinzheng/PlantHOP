#!/usr/bin/env python3
"""Encode Moricandia raw counts under the frozen v3 top-128 rule."""

import argparse
import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite


PROJECT = Path("/workspace/projects/PhyloOpenCell")
H5AD = Path("/data/datasets/Bioinformatics/external/PRJNA1186371_Moricandia_paper_object/m_arvensis_bbknn.h5ad")
MAPPING = PROJECT / "metadata/external_validation_v3_moricandia_feature_orthogroup_map_v1.tsv"
COUNT_PATHS = {
    "0": Path("/data/datasets/Bioinformatics/external/PRJNA1186371_Moricandia/assays/rnaseq/dataset/902/902-6_GE_cellranger_count/outs/filtered_feature_bc_matrix.h5"),
    "1": Path("/data/datasets/Bioinformatics/external/PRJNA1186371_Moricandia/assays/rnaseq/dataset/902/902-7_GE_cellranger_count/outs/filtered_feature_bc_matrix.h5"),
    "2": Path("/data/datasets/Bioinformatics/external/PRJNA1186371_Moricandia/assays/rnaseq/dataset/902/902-8_GE_cellranger_count/outs/filtered_feature_bc_matrix.h5"),
}
TOP_K = 128


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("smoke", "formal"))
    args = parser.parse_args()
    output_dir = Path("/data/derived/PhyloOpenCell/external_validation_v3") / f"moricandia_top128_{args.mode}"
    partial_dir = Path(str(output_dir) + ".partial")
    if output_dir.exists() or partial_dir.exists():
        raise RuntimeError("refusing to overwrite Moricandia encoding output")
    partial_dir.mkdir(parents=True)

    source_gene_to_og = {}
    with MAPPING.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["project_orthogroup_id"]:
                source_gene_to_og[row["moricandia_source_gene_id"]] = row["project_orthogroup_id"]
    og_levels = sorted(set(source_gene_to_og.values()))
    og_to_index = {orthogroup: index for index, orthogroup in enumerate(og_levels)}

    object_view = ad.read_h5ad(H5AD, backed="r")
    obs = object_view.obs[["replicate", "cell_type"]].copy()
    obs.index = obs.index.astype(str)

    output_i = []
    output_j = []
    output_x = []
    metadata_rows = []
    output_row = 0
    replicate_audit = {}
    for replicate, count_path in COUNT_PATHS.items():
        counts = sc.read_10x_h5(count_path)
        counts.X = counts.X.tocsr()
        transformed_names = np.asarray([f"{name}-{replicate}" for name in counts.obs_names.astype(str)])
        expected_names = set(obs.index[obs["replicate"].astype(str) == replicate])
        selected = np.flatnonzero(np.isin(transformed_names, list(expected_names)))
        if len(selected) != len(expected_names):
            raise RuntimeError(f"replicate {replicate}: expected {len(expected_names)} linked cells, got {len(selected)}")
        if args.mode == "smoke":
            selected = selected[:80]

        feature_og_indices = np.asarray(
            [og_to_index.get(source_gene_to_og.get(str(gene), ""), -1) for gene in counts.var_names],
            dtype=np.int64,
        )
        selected_counts = counts.X[selected]
        for local_row, source_row in enumerate(selected):
            start, end = selected_counts.indptr[local_row : local_row + 2]
            feature_indices = selected_counts.indices[start:end]
            values = selected_counts.data[start:end]
            mapped_indices = feature_og_indices[feature_indices]
            keep = (mapped_indices >= 0) & np.isfinite(values) & (values > 0)
            feature_indices = feature_indices[keep]
            mapped_indices = mapped_indices[keep]
            values = values[keep]
            ordering = np.lexsort((feature_indices, -values))[:TOP_K]
            best = {}
            for rank, mapped_index in enumerate(mapped_indices[ordering], start=1):
                weight = 1.0 / np.log2(rank + 1)
                best[mapped_index] = max(best.get(mapped_index, 0.0), weight)
            norm = np.sqrt(sum(value * value for value in best.values()))
            if norm > 0:
                for mapped_index, value in best.items():
                    output_i.append(output_row)
                    output_j.append(mapped_index)
                    output_x.append(value / norm)
            full_name = transformed_names[source_row]
            metadata_rows.append(
                {
                    "row_index": output_row + 1,
                    "cell_id": full_name,
                    "replicate": replicate,
                    "source_label": str(obs.loc[full_name, "cell_type"]),
                }
            )
            output_row += 1
        replicate_audit[replicate] = {
            "expected_linked_cells": len(expected_names),
            "encoded_cells": len(selected),
        }

    encoded = sp.coo_matrix(
        (output_x, (output_i, output_j)),
        shape=(output_row, len(og_levels)),
        dtype=np.float64,
    ).tocsr()
    mmwrite(partial_dir / "encoded.mtx", encoded)
    (partial_dir / "orthogroups.txt").write_text("\n".join(og_levels) + "\n", encoding="utf-8")
    with (partial_dir / "cell_metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metadata_rows)
    manifest = {
        "status": f"moricandia_external_validation_v3_{args.mode}_encoded",
        "mode": args.mode,
        "source_h5ad": str(H5AD),
        "source_h5ad_sha256": "be5f9be6570b78029fc138dd9c81ad908532776bbe665b9356b2055f554fd802",
        "mapping_path": str(MAPPING),
        "mapping_sha256": "dd858f7361988b29b539755913aed2ee28e99693a5abeafe1342b4b0b835a659",
        "cells": encoded.shape[0],
        "mapped_orthogroups": encoded.shape[1],
        "nonzero": encoded.nnz,
        "top_k": TOP_K,
        "replicates": replicate_audit,
        "expression_values_used_for_mapping_or_selection": False,
        "source_modified": False,
    }
    (partial_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    partial_dir.rename(output_dir)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
