#!/usr/bin/env python3
"""Read-only, collection-level audit for plant single-cell H5AD files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np


UNKNOWN_TOKENS = {
    "", "na", "n/a", "nan", "none", "null", "unknown", "unassigned",
    "unannotated", "unmapped", "other", "others", "not available",
}
LABEL_NAMES = (
    "celltype_after", "cell_type", "celltype", "cell type", "annotation",
    "cell_type_annotation", "label", "labels", "predicted_cell_type",
)
META_GROUPS = {
    "species": ("species", "organism", "scientific_name"),
    "tissue": ("tissue", "organ", "organ_part"),
    "stage": ("stage", "development_stage", "developmental_stage", "time", "timepoint"),
    "condition": ("condition", "treatment", "stimulus", "genotype"),
    "sample": ("sample", "sample_id", "orig.ident", "batch", "replicate", "donor"),
}


def scalar(value):
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
        return value if math.isfinite(value) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def sha256(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def infer_species(dataset_id: str) -> str:
    marker = re.split(r"_(?:PRJ|EMTAB|GSE|CRA|SRP|ERP|DRP|CNP)", dataset_id, maxsplit=1)
    return marker[0].replace("_", " ").strip()


def choose_column(columns, preferred):
    lookup = {str(c).strip().lower(): str(c) for c in columns}
    for name in preferred:
        if name.lower() in lookup:
            return lookup[name.lower()]
    return None


def matrix_sample(matrix, n_obs: int, n_vars: int, max_rows: int, max_cols: int):
    if n_obs == 0 or n_vars == 0:
        return np.empty((0, 0), dtype=float)
    row_count = min(max_rows, n_obs)
    col_count = min(max_cols, n_vars)
    starts = sorted(set([0, max(0, (n_obs - row_count) // 2), max(0, n_obs - row_count)]))
    blocks = []
    for start in starts:
        block = matrix[start:start + row_count, :col_count]
        if hasattr(block, "toarray"):
            block = block.toarray()
        blocks.append(np.asarray(block).reshape(-1))
    return np.concatenate(blocks) if blocks else np.empty(0, dtype=float)


def analyze_values(values: np.ndarray):
    if values.size == 0:
        return {
            "sampled_values": 0, "finite_fraction": None, "zero_fraction": None,
            "negative_fraction": None, "integer_like_fraction": None,
            "sample_min": None, "sample_max": None, "expression_semantics": "empty",
        }
    try:
        arr = values.astype(np.float64, copy=False)
    except (TypeError, ValueError):
        return {"sampled_values": int(values.size), "expression_semantics": "non_numeric"}
    finite = np.isfinite(arr)
    f = arr[finite]
    if f.size == 0:
        return {
            "sampled_values": int(arr.size), "finite_fraction": 0.0,
            "zero_fraction": None, "negative_fraction": None,
            "integer_like_fraction": None, "sample_min": None, "sample_max": None,
            "expression_semantics": "non_finite",
        }
    zero_fraction = float(np.mean(f == 0))
    negative_fraction = float(np.mean(f < 0))
    integer_like = float(np.mean(np.isclose(f, np.round(f), atol=1e-6)))
    if negative_fraction > 0:
        semantics = "scaled_continuous_or_centered"
    elif integer_like >= 0.999:
        semantics = "count_like_nonnegative_integer"
    else:
        semantics = "nonnegative_continuous_or_log_normalized"
    return {
        "sampled_values": int(arr.size),
        "finite_fraction": float(np.mean(finite)),
        "zero_fraction": zero_fraction,
        "negative_fraction": negative_fraction,
        "integer_like_fraction": integer_like,
        "sample_min": float(np.min(f)),
        "sample_max": float(np.max(f)),
        "expression_semantics": semantics,
    }


def gene_id_profile(index_values):
    values = [str(x) for x in list(index_values[: min(len(index_values), 5000)])]
    patterns = {
        "arabidopsis_locus": re.compile(r"^AT[1-5CM]G\d+", re.I),
        "ensembl": re.compile(r"^ENS[A-Z]*G\d+", re.I),
        "rice_locus": re.compile(r"^(LOC_)?Os\d+g\d+", re.I),
        "maize_locus": re.compile(r"^(Zm\d+d\d+|GRMZM|Zm00001)", re.I),
        "symbol_like": re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{1,30}$"),
    }
    counts = Counter()
    for value in values:
        matched = False
        for name, pattern in patterns.items():
            if pattern.search(value):
                counts[name] += 1
                matched = True
                break
        if not matched:
            counts["other"] += 1
    return dict(counts)


def audit_one(path: Path, args):
    result = {
        "dataset_id": path.stem,
        "path": str(path),
        "file_size_bytes": path.stat().st_size,
        "species_inferred": infer_species(path.stem),
        "status": "ok",
        "warnings": [],
    }
    if args.sha256:
        result["sha256"] = sha256(path)
    adata = None
    try:
        adata = ad.read_h5ad(path, backed="r")
        result.update({
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "obs_names_unique": bool(adata.obs_names.is_unique),
            "var_names_unique": bool(adata.var_names.is_unique),
            "obs_columns": [str(c) for c in adata.obs.columns],
            "var_columns": [str(c) for c in adata.var.columns],
            "layers": [str(k) for k in adata.layers.keys()],
            "obsm": [str(k) for k in adata.obsm.keys()],
            "varm": [str(k) for k in adata.varm.keys()],
            "obsp": [str(k) for k in adata.obsp.keys()],
            "uns_keys": [str(k) for k in adata.uns.keys()],
            "x_class": type(adata.X).__name__,
            "x_dtype": str(getattr(adata.X, "dtype", "unknown")),
            "gene_id_profile": gene_id_profile(adata.var_names),
        })
        x_sample = matrix_sample(
            adata.X, adata.n_obs, adata.n_vars, args.sample_rows, args.sample_cols
        )
        result.update(analyze_values(x_sample))

        label_col = choose_column(adata.obs.columns, LABEL_NAMES)
        result["label_column"] = label_col
        if label_col:
            series = adata.obs[label_col]
            strings = series.astype("string").fillna("").str.strip()
            lower = strings.str.lower()
            unknown = lower.isin(UNKNOWN_TOKENS)
            known = strings[~unknown]
            counts = known.value_counts(dropna=False)
            result.update({
                "labeled_cells": int((~unknown).sum()),
                "unknown_or_missing_cells": int(unknown.sum()),
                "label_coverage": float((~unknown).mean()) if len(strings) else None,
                "n_known_cell_types": int(len(counts)),
                "min_known_type_cells": int(counts.min()) if len(counts) else None,
                "median_known_type_cells": float(counts.median()) if len(counts) else None,
                "rare_types_lt50": int((counts < 50).sum()) if len(counts) else 0,
                "rare_types_lt1pct": int((counts < max(1, len(strings) * 0.01)).sum()) if len(counts) else 0,
                "top_cell_types": {str(k): int(v) for k, v in counts.head(20).items()},
            })
            if result["label_coverage"] == 0:
                result["warnings"].append("cell_type_column_present_but_all_unknown_or_missing")
            elif result["label_coverage"] < 0.5:
                result["warnings"].append("cell_type_coverage_below_50_percent")
        else:
            result.update({
                "labeled_cells": 0, "unknown_or_missing_cells": int(adata.n_obs),
                "label_coverage": 0.0, "n_known_cell_types": 0,
                "min_known_type_cells": None, "median_known_type_cells": None,
                "rare_types_lt50": 0, "rare_types_lt1pct": 0, "top_cell_types": {},
            })
            result["warnings"].append("no_recognized_cell_type_column")

        metadata = {}
        for group, candidates in META_GROUPS.items():
            column = choose_column(adata.obs.columns, candidates)
            metadata[group] = column
            if column:
                values = adata.obs[column].astype("string").fillna("").str.strip()
                metadata[group + "_nonempty_fraction"] = float((values != "").mean()) if len(values) else None
                metadata[group + "_n_unique"] = int(values[values != ""].nunique())
        result["metadata_columns"] = metadata

        if "cell_type_canonical" in adata.obs.columns:
            canonical = adata.obs["cell_type_canonical"].astype("string").fillna("").str.strip()
            result["canonical_cell_type_nonempty_cells"] = int((canonical != "").sum())
            result["canonical_cell_type_coverage"] = float((canonical != "").mean()) if len(canonical) else None
        else:
            result["canonical_cell_type_nonempty_cells"] = 0
            result["canonical_cell_type_coverage"] = 0.0
        if "cell_type_mapping_status" in adata.obs.columns:
            values = adata.obs["cell_type_mapping_status"].astype("string").fillna("")
            result["cell_type_mapping_status_counts"] = {
                str(k): int(v) for k, v in values.value_counts(dropna=False).items()
            }
        if "gene_id_canonical" in adata.var.columns:
            canonical = adata.var["gene_id_canonical"].astype("string").fillna("").str.strip()
            result["canonical_gene_nonempty"] = int((canonical != "").sum())
            result["canonical_gene_coverage"] = float((canonical != "").mean()) if len(canonical) else None
        else:
            result["canonical_gene_nonempty"] = 0
            result["canonical_gene_coverage"] = 0.0
        if "gene_mapping_status" in adata.var.columns:
            values = adata.var["gene_mapping_status"].astype("string").fillna("")
            result["gene_mapping_status_counts"] = {
                str(k): int(v) for k, v in values.value_counts(dropna=False).items()
            }
        if result.get("canonical_cell_type_coverage", 0) == 0:
            result["warnings"].append("canonical_cell_type_mapping_empty")
        if result.get("canonical_gene_coverage", 0) == 0:
            result["warnings"].append("canonical_gene_mapping_empty")

        expression_meta = adata.uns.get("expression_representation", {})
        result["declared_expression_representation"] = scalar(expression_meta.get("class")) if hasattr(expression_meta, "get") else None
        result["declared_storage_kind"] = scalar(expression_meta.get("storage_kind")) if hasattr(expression_meta, "get") else None
        result["schema_version"] = scalar(adata.uns.get("schema_version"))

        for layer_name in adata.layers.keys():
            layer = adata.layers[layer_name]
            layer_values = matrix_sample(layer, adata.n_obs, adata.n_vars, args.sample_rows, args.sample_cols)
            result.setdefault("layer_profiles", {})[str(layer_name)] = {
                "class": type(layer).__name__,
                "dtype": str(getattr(layer, "dtype", "unknown")),
                "sample_equal_to_x": bool(
                    layer_values.shape == x_sample.shape
                    and np.array_equal(layer_values, x_sample, equal_nan=True)
                ),
                **analyze_values(layer_values),
            }
            if str(layer_name).strip().lower() in {"", "none", "null"}:
                result["warnings"].append("nonsemantic_expression_layer_name")

        if not result["obs_names_unique"]:
            result["warnings"].append("duplicate_cell_ids")
        if not result["var_names_unique"]:
            result["warnings"].append("duplicate_gene_ids")
        if result.get("finite_fraction") is not None and result["finite_fraction"] < 1:
            result["warnings"].append("non_finite_expression_values_in_sample")
        if not result["layers"]:
            result["warnings"].append("no_expression_layers")
        if not result["metadata_columns"].get("tissue"):
            result["warnings"].append("no_tissue_column")
        if not result["metadata_columns"].get("sample"):
            result["warnings"].append("no_sample_or_batch_column")
    except Exception as exc:
        result["status"] = "error"
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc(limit=5)
    finally:
        if adata is not None and getattr(adata, "file", None) is not None:
            adata.file.close()
    return result


def csv_value(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return scalar(value)


def build_summary(records, started_at, completed_at):
    ok = [r for r in records if r["status"] == "ok"]
    semantics = Counter(r.get("expression_semantics", "unknown") for r in ok)
    species = Counter(r.get("species_inferred", "") for r in ok)
    warning_counts = Counter(w for r in ok for w in r.get("warnings", []))
    return {
        "audit_version": "1.1",
        "started_at_utc": started_at,
        "completed_at_utc": completed_at,
        "files_total": len(records),
        "files_ok": len(ok),
        "files_error": len(records) - len(ok),
        "cells_total": sum(r.get("n_cells", 0) for r in ok),
        "genes_min": min((r.get("n_genes", 0) for r in ok), default=None),
        "genes_max": max((r.get("n_genes", 0) for r in ok), default=None),
        "file_size_bytes_total": sum(r.get("file_size_bytes", 0) for r in records),
        "species_inferred_count": len([s for s in species if s]),
        "species_dataset_counts": dict(sorted(species.items())),
        "expression_semantics_counts": dict(semantics),
        "datasets_with_recognized_labels": sum(bool(r.get("label_column")) for r in ok),
        "datasets_with_sample_metadata": sum(bool(r.get("metadata_columns", {}).get("sample")) for r in ok),
        "datasets_with_tissue_metadata": sum(bool(r.get("metadata_columns", {}).get("tissue")) for r in ok),
        "datasets_with_layers": sum(bool(r.get("layers")) for r in ok),
        "labeled_cells_total": sum(r.get("labeled_cells", 0) for r in ok),
        "unknown_or_missing_cells_total": sum(r.get("unknown_or_missing_cells", 0) for r in ok),
        "datasets_with_zero_label_coverage": sum(r.get("label_coverage", 0) == 0 for r in ok),
        "datasets_with_label_coverage_below_50_percent": sum(r.get("label_coverage", 0) < 0.5 for r in ok),
        "canonical_cell_type_nonempty_cells_total": sum(r.get("canonical_cell_type_nonempty_cells", 0) for r in ok),
        "canonical_gene_nonempty_total": sum(r.get("canonical_gene_nonempty", 0) for r in ok),
        "warning_counts": dict(warning_counts),
    }


def write_markdown(path: Path, summary, records):
    lines = [
        "# Plant single-cell H5AD audit",
        "",
        f"Completed: {summary['completed_at_utc']}",
        "",
        "## Collection summary",
        "",
        f"- Files: {summary['files_ok']}/{summary['files_total']} opened successfully",
        f"- Cells: {summary['cells_total']:,}",
        f"- Inferred species: {summary['species_inferred_count']}",
        f"- Total H5AD bytes: {summary['file_size_bytes_total']:,}",
        f"- Recognized label columns: {summary['datasets_with_recognized_labels']}/{summary['files_ok']}",
        f"- Sample or batch metadata: {summary['datasets_with_sample_metadata']}/{summary['files_ok']}",
        f"- Tissue metadata: {summary['datasets_with_tissue_metadata']}/{summary['files_ok']}",
        f"- Expression layers present: {summary['datasets_with_layers']}/{summary['files_ok']}",
        f"- Labeled cells: {summary['labeled_cells_total']:,}",
        f"- Unknown or missing labels: {summary['unknown_or_missing_cells_total']:,}",
        f"- Datasets below 50% label coverage: {summary['datasets_with_label_coverage_below_50_percent']}",
        f"- Canonical cell-type assignments populated: {summary['canonical_cell_type_nonempty_cells_total']:,}",
        f"- Canonical gene IDs populated: {summary['canonical_gene_nonempty_total']:,}",
        "",
        "## Sampled expression semantics",
        "",
    ]
    for key, value in sorted(summary["expression_semantics_counts"].items()):
        lines.append(f"- {key}: {value}")
    lines += ["", "## Warning counts", ""]
    for key, value in sorted(summary["warning_counts"].items()):
        lines.append(f"- {key}: {value}")
    lines += [
        "", "## Dataset table", "",
        "| Dataset | Species | Cells | Genes | Labels | Coverage | Matrix semantics | Warnings |",
        "|---|---|---:|---:|---:|---:|---|---|",
    ]
    for r in records:
        coverage = r.get("label_coverage")
        coverage_text = "" if coverage is None else f"{coverage:.1%}"
        lines.append(
            f"| {r['dataset_id']} | {r.get('species_inferred','')} | "
            f"{r.get('n_cells','')} | {r.get('n_genes','')} | {r.get('n_known_cell_types','')} | "
            f"{coverage_text} | {r.get('expression_semantics','')} | "
            f"{', '.join(r.get('warnings', []))} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample-rows", type=int, default=256)
    parser.add_argument("--sample-cols", type=int, default=1024)
    parser.add_argument("--sha256", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = sorted(args.input_dir.glob("*.h5ad"))
    started = datetime.now(timezone.utc).isoformat()
    records = []
    for index, path in enumerate(paths, start=1):
        print(f"[{index}/{len(paths)}] {path.name}", flush=True)
        records.append(audit_one(path, args))
    completed = datetime.now(timezone.utc).isoformat()
    summary = build_summary(records, started, completed)
    payload = {"summary": summary, "datasets": records}
    (args.output_dir / "audit_full.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=scalar) + "\n", encoding="utf-8"
    )
    fields = sorted({key for record in records for key in record.keys()} - {"traceback"})
    with (args.output_dir / "audit_datasets.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({key: csv_value(record.get(key)) for key in fields})
    write_markdown(args.output_dir / "audit_report.md", summary, records)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    raise SystemExit(1 if summary["files_error"] else 0)


if __name__ == "__main__":
    main()
