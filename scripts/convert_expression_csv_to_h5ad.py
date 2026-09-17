#!/usr/bin/env python3
"""Convert one expression/annotation CSV pair to the project's H5AD schema."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def log(message: str, log_path: Path) -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = f"{stamp}\t{message}"
    print(line, flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def categorical(series: pd.Series) -> pd.Series:
    return series.astype("category")


def convert(args: argparse.Namespace) -> dict[str, object]:
    expression_path = Path(args.expression)
    annotation_path = Path(args.annotation)
    output_path = Path(args.output)
    partial_path = output_path.with_suffix(output_path.suffix + ".partial")
    log_path = Path(args.log)

    for path in (expression_path, annotation_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {output_path}")
    if partial_path.exists():
        partial_path.unlink()

    started = time.monotonic()
    log(f"START\t{args.dataset_id}", log_path)
    source_hashes = {
        "expression": sha256(expression_path),
        "annotation": sha256(annotation_path),
    }
    log(f"HASHED\t{args.dataset_id}\t{json.dumps(source_hashes, sort_keys=True)}", log_path)

    expression_columns = pd.read_csv(expression_path, nrows=0).columns
    numeric_dtypes = {name: np.float32 for name in expression_columns[1:]}
    expression = pd.read_csv(expression_path, index_col=0, dtype=numeric_dtypes)
    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)
    if not expression.index.is_unique:
        raise ValueError("Expression cell identifiers are not unique")
    if not expression.columns.is_unique:
        raise ValueError("Expression gene identifiers are not unique")
    if expression.empty:
        raise ValueError("Expression matrix is empty")

    annotation_raw = pd.read_csv(annotation_path, dtype=str, keep_default_na=False)
    if annotation_raw.shape[1] < 2:
        raise ValueError("Annotation file has fewer than two columns")
    source_names = list(annotation_raw.columns)
    source_names[0] = ""
    stored_names = list(annotation_raw.columns)
    stored_names[0] = "source__unnamed"
    stored_names = ["source__index" if name == "_index" else name for name in stored_names]
    if len(set(stored_names)) != len(stored_names):
        raise ValueError("Annotation column names collide after normalization")
    annotation_raw.columns = stored_names
    annotation_raw["source__unnamed"] = annotation_raw["source__unnamed"].astype(str)
    annotation_raw.index = pd.Index(annotation_raw["source__unnamed"].astype(str), dtype=str)

    if len(expression) != len(annotation_raw):
        raise ValueError(
            f"Row-count mismatch: expression={len(expression)}, annotation={len(annotation_raw)}"
        )
    if not expression.index.equals(annotation_raw.index):
        mismatch = np.flatnonzero(expression.index.to_numpy() != annotation_raw.index.to_numpy())
        position = int(mismatch[0]) if len(mismatch) else -1
        raise ValueError(f"Cell identifier/order mismatch at row {position}")

    matrix = expression.to_numpy(dtype=np.float32, copy=False)
    finite_values = 0
    negative_values = 0
    zero_values = 0
    all_values_integer = True
    for start in range(0, matrix.shape[0], 1024):
        block = matrix[start : start + 1024]
        finite = np.isfinite(block)
        if not finite.all():
            bad = int(block.size - finite.sum())
            raise ValueError(f"Expression matrix contains {bad} non-finite values")
        finite_values += int(finite.sum())
        negative_values += int(np.count_nonzero(block < 0))
        zero_values += int(np.count_nonzero(block == 0))
        if all_values_integer and not np.equal(block, np.floor(block)).all():
            all_values_integer = False

    obs = annotation_raw.copy()
    if "celltype_after" not in obs.columns:
        obs["celltype_after"] = ""
    for column in obs.columns:
        if column not in {"source__unnamed", "source__index"}:
            obs[column] = categorical(obs[column])
    if "source__index" in obs.columns:
        obs["source__index"] = obs["source__index"].astype(str)

    original_type = obs["celltype_after"].astype(str)
    obs["dataset_id"] = categorical(pd.Series(args.dataset_id, index=obs.index))
    obs["species"] = categorical(pd.Series(args.species, index=obs.index))
    obs["accession"] = categorical(pd.Series(args.accession, index=obs.index))
    obs["tissue"] = categorical(pd.Series(args.tissue, index=obs.index))
    obs["cell_type_original"] = categorical(original_type)
    obs["cell_type_canonical"] = categorical(pd.Series("", index=obs.index))
    obs["cell_type_mapping_status"] = categorical(pd.Series("unmapped", index=obs.index))
    unknown_labels = {"", "unknown", "unassigned", "na", "n/a", "none", "nan"}
    obs["cell_type_is_unknown"] = original_type.str.strip().str.lower().isin(unknown_labels).to_numpy()

    genes = expression.columns
    var = pd.DataFrame(index=pd.Index(genes, dtype=str))
    var["gene_id_original"] = pd.Series(genes, index=var.index, dtype=str)
    var["gene_id_canonical"] = categorical(pd.Series("", index=var.index))
    var["gene_mapping_status"] = categorical(pd.Series("unmapped", index=var.index))
    var["species"] = categorical(pd.Series(args.species, index=var.index))
    var["dataset_id"] = categorical(pd.Series(args.dataset_id, index=var.index))

    zero_fraction = zero_values / matrix.size
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    result.uns["conversion_provenance"] = {
        "compression": "lzf",
        "converted_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "converter": "convert_expression_csv_to_h5ad.py",
        "zero_fraction_for_csr": 0.7,
    }
    result.uns["expression_representation"] = {
        "class": "scaled_continuous" if negative_values or not all_values_integer else "unverified_nonnegative",
        "dtype": "float32",
        "evidence": {
            "all_values_integer": bool(all_values_integer),
            "finite_values": finite_values,
            "negative_values": negative_values,
            "rule": "raw_counts_verified is never inferred from integer-valued CSV values",
            "zero_fraction": zero_fraction,
            "zero_values": zero_values,
        },
        "storage_kind": "dense",
    }
    result.uns["schema_version"] = "1.0.0"
    result.uns["source"] = {
        "annotation_path": str(annotation_path),
        "expression_path": str(expression_path),
        "sha256": source_hashes,
    }
    result.uns["source_annotation_column_name_map"] = {
        "source_names": np.asarray(source_names, dtype=object),
        "stored_names": np.asarray(stored_names, dtype=object),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(partial_path, compression="lzf")
    del result, expression, annotation_raw, matrix

    check = ad.read_h5ad(partial_path, backed="r")
    try:
        if check.shape != (len(obs), len(var)):
            raise ValueError(f"Written shape mismatch: {check.shape}")
        if list(check.obs_names) != list(obs.index):
            raise ValueError("Written observation identifiers differ")
        if list(check.var_names) != list(var.index):
            raise ValueError("Written variable identifiers differ")
        required_obs = {
            "dataset_id", "species", "accession", "tissue", "cell_type_original",
            "cell_type_canonical", "cell_type_mapping_status", "cell_type_is_unknown",
        }
        if not required_obs.issubset(check.obs.columns):
            raise ValueError("Written file is missing required observation columns")
    finally:
        check.file.close()

    os.replace(partial_path, output_path)
    summary = {
        "dataset_id": args.dataset_id,
        "output": str(output_path),
        "shape": [len(obs), len(var)],
        "bytes": output_path.stat().st_size,
        "sha256": sha256(output_path),
        "seconds": round(time.monotonic() - started, 2),
        "negative_values": negative_values,
        "zero_fraction": zero_fraction,
    }
    log(f"COMPLETE\t{json.dumps(summary, ensure_ascii=False, sort_keys=True)}", log_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--tissue", default="not explicitly recorded")
    parser.add_argument("--log", required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    try:
        print(json.dumps(convert(args), ensure_ascii=False, sort_keys=True), flush=True)
        return 0
    except Exception as exc:
        partial = Path(args.output).with_suffix(Path(args.output).suffix + ".partial")
        if partial.exists():
            partial.unlink()
        log(f"SKIPPED\t{args.dataset_id}\t{type(exc).__name__}: {exc}", Path(args.log))
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
