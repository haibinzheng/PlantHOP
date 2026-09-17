#!/usr/bin/env python3
"""Create non-destructive governance and feasibility artifacts for PhyloOpenCell."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd


UNKNOWN = {"", "unknown", "unassigned", "na", "n/a", "none", "nan"}
TISSUE_HINTS = {
    "root": "root",
    "shoot": "shoot",
    "leaf": "leaf",
    "flower": "flower",
    "inflorescence": "inflorescence",
    "spikelet": "spikelet",
    "stem": "stem",
    "nodule": "nodule",
    "seed": "seed",
    "cotyledon": "cotyledon",
    "anther": "anther",
}


def text_scalar(frame: pd.DataFrame, column: str, default: str = "") -> str:
    if column not in frame.columns or frame.empty:
        return default
    value = str(frame[column].iloc[0]).strip()
    return value if value.lower() not in {"nan", "none"} else default


def infer_species(dataset_id: str) -> str:
    return dataset_id.split("_", 1)[0].strip()


def infer_tissue(dataset_id: str) -> str:
    tail = dataset_id.rsplit("_", 1)[-1].lower()
    for hint, normalized in TISSUE_HINTS.items():
        if hint in tail:
            return normalized
    return ""


def normalize_label(label: str) -> str:
    value = re.sub(r"[_\-/]+", " ", label.strip().lower())
    return re.sub(r"\s+", " ", value)


def gene_id_style(names: list[str]) -> str:
    sample = names[: min(len(names), 5000)]
    rules = [
        ("arabidopsis_locus", re.compile(r"^AT[1-5MC]G\d+", re.I)),
        ("rice_locus", re.compile(r"^Os\d+g\d+", re.I)),
        ("maize_locus", re.compile(r"^(Zm\d+|GRMZM|AC\d+)", re.I)),
        ("ensembl_plant", re.compile(r"^[A-Z]{2,8}\d{4,}")),
    ]
    scores = [(name, sum(bool(rx.search(x)) for x in sample)) for name, rx in rules]
    best, count = max(scores, key=lambda item: item[1], default=("mixed_or_unknown", 0))
    return best if sample and count / len(sample) >= 0.5 else "mixed_or_unknown"


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad-dir", required=True)
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args()

    h5ad_dir = Path(args.h5ad_dir)
    project = Path(args.project_dir)
    metadata = project / "metadata"
    reports = project / "reports"
    metadata.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    dataset_rows: list[dict] = []
    schema_rows: list[dict] = []
    gene_rows: list[dict] = []
    repair_rows: list[dict] = []
    label_counts: Counter[str] = Counter()
    label_datasets: dict[str, set[str]] = defaultdict(set)
    label_species: dict[str, set[str]] = defaultdict(set)
    species_genes: dict[str, list[set[str]]] = defaultdict(list)
    errors: list[dict] = []

    files = sorted(h5ad_dir.glob("*.h5ad"))
    for path in files:
        try:
            obj = ad.read_h5ad(path, backed="r")
            dataset_id = text_scalar(obj.obs, "dataset_id", path.stem)
            observed_species = text_scalar(obj.obs, "species", "not explicitly recorded")
            proposed_species = infer_species(dataset_id)
            species_status = "recorded"
            if observed_species.lower() in UNKNOWN | {"not explicitly recorded"}:
                species_status = "candidate_from_dataset_id"
                repair_rows.append({
                    "dataset_id": dataset_id,
                    "field": "species",
                    "observed_value": observed_species,
                    "proposed_value": proposed_species,
                    "basis": "dataset_id prefix; requires source-publication verification",
                    "status": "candidate_not_applied",
                })

            observed_tissue = text_scalar(obj.obs, "tissue", "not explicitly recorded")
            proposed_tissue = infer_tissue(dataset_id)
            tissue_status = "recorded"
            if observed_tissue.lower() in UNKNOWN | {"not explicitly recorded"}:
                tissue_status = "missing"
                if proposed_tissue:
                    tissue_status = "candidate_from_dataset_id"
                    repair_rows.append({
                        "dataset_id": dataset_id,
                        "field": "tissue",
                        "observed_value": observed_tissue,
                        "proposed_value": proposed_tissue,
                        "basis": "dataset_id suffix; requires source-publication verification",
                        "status": "candidate_not_applied",
                    })

            unknown_count = int(obj.obs["cell_type_is_unknown"].sum()) if "cell_type_is_unknown" in obj.obs else obj.n_obs
            labels = obj.obs["cell_type_original"].astype(str) if "cell_type_original" in obj.obs else pd.Series([], dtype=str)
            for label, count in labels.value_counts(dropna=False).items():
                label = str(label)
                label_counts[label] += int(count)
                label_datasets[label].add(dataset_id)
                label_species[label].add(proposed_species)

            representation = obj.uns.get("expression_representation", {})
            names = [str(x) for x in obj.var_names]
            species_genes[proposed_species].append(set(names))
            mapped = 0
            if "gene_id_canonical" in obj.var:
                mapped = int((obj.var["gene_id_canonical"].astype(str).str.strip() != "").sum())
            gene_rows.append({
                "dataset_id": dataset_id,
                "species": proposed_species,
                "n_genes": obj.n_vars,
                "unique_gene_ids": len(set(names)),
                "canonical_gene_ids_present": mapped,
                "canonical_mapping_fraction": mapped / obj.n_vars if obj.n_vars else 0,
                "gene_id_style": gene_id_style(names),
                "example_gene_ids": "|".join(names[:10]),
            })

            for column in obj.obs.columns:
                schema_rows.append({
                    "dataset_id": dataset_id,
                    "axis": "obs",
                    "column": column,
                    "dtype": str(obj.obs[column].dtype),
                })
            for column in obj.var.columns:
                schema_rows.append({
                    "dataset_id": dataset_id,
                    "axis": "var",
                    "column": column,
                    "dtype": str(obj.var[column].dtype),
                })

            source = obj.uns.get("source", {})
            dataset_rows.append({
                "dataset_id": dataset_id,
                "h5ad_path": str(path),
                "readable": True,
                "file_bytes": path.stat().st_size,
                "n_cells": obj.n_obs,
                "n_genes": obj.n_vars,
                "observed_species": observed_species,
                "proposed_species": proposed_species,
                "species_status": species_status,
                "accession": text_scalar(obj.obs, "accession"),
                "observed_tissue": observed_tissue,
                "proposed_tissue": proposed_tissue,
                "tissue_status": tissue_status,
                "unknown_cells": unknown_count,
                "unknown_fraction": unknown_count / obj.n_obs if obj.n_obs else 1,
                "expression_class": representation.get("class", "missing"),
                "storage_kind": representation.get("storage_kind", type(obj.X).__name__),
                "expression_dtype": representation.get("dtype", str(obj.X.dtype)),
                "source_expression_path": source.get("expression_path", ""),
                "source_annotation_path": source.get("annotation_path", ""),
                "source_expression_sha256": source.get("sha256", {}).get("expression", ""),
                "source_annotation_sha256": source.get("sha256", {}).get("annotation", ""),
                "schema_version": obj.uns.get("schema_version", ""),
            })
            obj.file.close()
        except Exception as exc:
            errors.append({"path": str(path), "error": f"{type(exc).__name__}: {exc}"})

    label_rows = []
    for label, count in label_counts.most_common():
        normalized = normalize_label(label)
        label_rows.append({
            "source_label": label,
            "normalized_lexical_label": normalized,
            "cell_count": count,
            "dataset_count": len(label_datasets[label]),
            "species_count": len(label_species[label]),
            "is_unknown_label": normalized in UNKNOWN,
            "ontology_id": "",
            "ontology_label": "",
            "ontology_level": "",
            "mapping_confidence": "",
            "mapping_status": "pending_manual_curation" if normalized not in UNKNOWN else "unknown",
        })

    species_rows = []
    frame = pd.DataFrame(dataset_rows)
    for species, group in frame.groupby("proposed_species", sort=True):
        gene_sets = species_genes[species]
        union = set().union(*gene_sets) if gene_sets else set()
        intersection = set.intersection(*gene_sets) if gene_sets else set()
        species_rows.append({
            "species": species,
            "dataset_count": len(group),
            "cell_count": int(group["n_cells"].sum()),
            "unknown_cells": int(group["unknown_cells"].sum()),
            "unknown_fraction": float(group["unknown_cells"].sum() / group["n_cells"].sum()),
            "gene_union_size": len(union),
            "gene_intersection_size": len(intersection),
            "expression_classes": "|".join(sorted(set(group["expression_class"]))),
        })

    eligible_rows = [
        {
            "source_label": row["source_label"],
            "cell_count": row["cell_count"],
            "dataset_count": row["dataset_count"],
            "species_count": row["species_count"],
            "eligibility_rule": ">=3 species, >=3 datasets, >=500 cells, known label",
        }
        for row in label_rows
        if not row["is_unknown_label"]
        and row["species_count"] >= 3
        and row["dataset_count"] >= 3
        and row["cell_count"] >= 500
    ]

    write_tsv(dataset_rows, metadata / "dataset_manifest.tsv")
    write_tsv(species_rows, metadata / "species_summary.tsv")
    write_tsv(label_rows, metadata / "cell_type_vocabulary.tsv")
    write_tsv(eligible_rows, metadata / "cross_species_label_candidates.tsv")
    write_tsv(schema_rows, metadata / "schema_inventory.tsv")
    write_tsv(gene_rows, metadata / "gene_namespace_audit.tsv")
    write_tsv(repair_rows, metadata / "metadata_repair_candidates.tsv")

    total_cells = int(frame["n_cells"].sum()) if not frame.empty else 0
    unknown_cells = int(frame["unknown_cells"].sum()) if not frame.empty else 0
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_directory": str(h5ad_dir),
        "h5ad_files": len(files),
        "readable_files": len(dataset_rows),
        "errors": errors,
        "cells": total_cells,
        "unknown_cells": unknown_cells,
        "unknown_fraction": unknown_cells / total_cells if total_cells else 1,
        "proposed_species": len(species_rows),
        "source_cell_type_labels": len(label_rows),
        "cross_species_label_candidates": len(eligible_rows),
        "metadata_repair_candidates": len(repair_rows),
        "expression_class_counts": frame["expression_class"].value_counts().to_dict() if not frame.empty else {},
        "governance_policy": "Source H5AD files were opened read-only; proposed repairs were not applied.",
    }
    with (reports / "governance_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)

    feasibility = f"""# PhyloOpenCell initial feasibility assessment

Generated: {summary['generated_at_utc']}

## Evidence snapshot

- H5AD files: {summary['h5ad_files']} ({summary['readable_files']} readable; {len(errors)} errors)
- Cells: {total_cells:,}
- Unknown/unassigned cells: {unknown_cells:,} ({summary['unknown_fraction']:.1%})
- Proposed species identities: {len(species_rows)}
- Source cell-type labels: {len(label_rows)}
- Preliminary cross-species label candidates: {len(eligible_rows)}
- Metadata repair candidates: {len(repair_rows)}

## Decision

**Conditionally feasible.** The corpus is large enough for a compact cross-species encoder and strict leave-one-dataset/species-out experiments. Training must not begin until label ontology mapping and orthogroup coverage are quantified. Unknown labels must remain unlabeled rather than be merged into one biological class.

## Required gates before model training

1. Verify all proposed species and tissue repairs against source publications.
2. Map source labels to a hierarchical plant cell ontology with confidence and provenance.
3. Build gene-to-orthogroup/protein mappings and report coverage per dataset.
4. Freeze dataset-level and species-level train/validation/test splits before optimization.
5. Define rank-based expression inputs for heterogeneous matrix representations.
6. Select external datasets that are absent from the training corpus.

## Primary evaluation tasks

- Cross-dataset transfer within species.
- Leave-one-species-out transfer.
- Monocot-to-eudicot and eudicot-to-monocot transfer.
- Held-out-cell-type open-set recognition.
- Low-shot adaptation and probability calibration.
"""
    (reports / "feasibility_report.md").write_text(feasibility, encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
