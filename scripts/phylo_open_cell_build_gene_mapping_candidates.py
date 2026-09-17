#!/usr/bin/env python3
"""Build reversible gene mapping candidates for high-coverage datasets only."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gtf_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    pattern = re.compile(r'gene_id "([^"]+)"')
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                ids.add(match.group(1))
    return ids


def gff3_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            if attrs.get("gene_id"):
                ids.add(attrs["gene_id"])
    return ids


def uniprot_gene_index(path: Path) -> tuple[set[str], dict[str, set[str]]]:
    table = pd.read_csv(path, sep="\t", compression="gzip", keep_default_na=False)
    gene_col = next(column for column in table.columns if column.casefold().startswith("gene names"))
    accession_col = next(column for column in table.columns if column.casefold() == "entry")
    index: dict[str, set[str]] = defaultdict(set)
    for _, row in table.iterrows():
        for token in str(row[gene_col]).split():
            index[token].add(str(row[accession_col]))
    return set(index), index


def transformed(species: str, value: str) -> list[tuple[str, str]]:
    candidates: list[tuple[str, str]] = []
    if species in {"Glycine max", "Sorghum bicolor"}:
        candidates += [(value.replace("-", "."), "hyphen_to_dot"), (value.replace("-", "_"), "hyphen_to_underscore")]
    if species == "Catharanthus roseus":
        candidates.append((value.replace("-", "_"), "hyphen_to_underscore"))
    if species == "Populus trichocarpa":
        match = re.fullmatch(r"POPTR-(\d+G\d+)v4", value, flags=re.I)
        if match:
            candidates.append((f"Potri.{match.group(1)}.v4.1", "poptr_v4_to_potri_v4_1"))
    return candidates


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    project = Path(args.project)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".partial")

    action = pd.read_csv(project / "metadata" / "reference_mapping_action_queue.tsv", sep="\t", keep_default_na=False)
    allowed = action[action.triage_status.isin(["direct_exact_candidate", "deterministic_transform_candidate"])]
    manifest = pd.read_csv(project / "metadata" / "dataset_manifest.tsv", sep="\t", keep_default_na=False).set_index("dataset_id")
    references = pd.read_csv(project / "metadata" / "ensembl_reference_downloads.tsv", sep="\t", keep_default_na=False)
    reference_by_species = {
        row.species: Path(row.local_path)
        for _, row in references[references.status.eq("downloaded_and_parsed")].iterrows()
    }
    maize = Path("/data/datasets/PhyloOpenCell/references/ensembl_plants_release_63/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.gff3.gz")
    catharanthus = Path("/data/datasets/PhyloOpenCell/references/uniprot/UP001060085.uniprot_gene_names.tsv.gz")

    id_sets: dict[str, set[str]] = {}
    protein_index: dict[str, dict[str, set[str]]] = {}
    for species in allowed.species.unique():
        if species == "Zea mays":
            id_sets[species] = gff3_ids(maize)
        elif species == "Catharanthus roseus":
            id_sets[species], protein_index[species] = uniprot_gene_index(catharanthus)
        else:
            id_sets[species] = gtf_ids(reference_by_species[species])

    header = [
        "dataset_id", "species", "source_gene_id", "candidate_reference_gene_id",
        "mapping_status", "mapping_rule", "candidate_protein_accessions",
        "protein_cardinality", "reference_sha256", "accepted",
    ]
    reference_hashes = {}
    for species in id_sets:
        path = catharanthus if species == "Catharanthus roseus" else maize if species == "Zea mays" else reference_by_species[species]
        reference_hashes[species] = sha256(path)

    summaries = []
    with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for _, item in allowed.iterrows():
            dataset_id = item.dataset_id
            species = item.species
            adata = ad.read_h5ad(manifest.loc[dataset_id, "h5ad_path"], backed="r")
            try:
                values = [str(x) for x in adata.var_names]
            finally:
                adata.file.close()
            counts = Counter()
            targets = Counter()
            protein_cardinalities = Counter()
            for value in values:
                matches: list[tuple[str, str]] = []
                if value in id_sets[species]:
                    matches.append((value, "exact"))
                for candidate, rule in transformed(species, value):
                    if candidate in id_sets[species] and candidate != value:
                        matches.append((candidate, rule))
                unique = sorted(set(matches))
                if not unique:
                    status, target, rule = "unmapped", "", ""
                elif len({x[0] for x in unique}) > 1:
                    status, target, rule = "ambiguous_candidate", "|".join(x[0] for x in unique), "|".join(x[1] for x in unique)
                else:
                    target = unique[0][0]
                    rules = {x[1] for x in unique}
                    rule = "exact" if "exact" in rules else sorted(rules)[0]
                    status = "exact_candidate" if rule == "exact" else "deterministic_transform_candidate"
                    targets[target] += 1
                proteins = sorted(protein_index.get(species, {}).get(target, set())) if target and "|" not in target else []
                cardinality = len(proteins)
                protein_cardinalities[cardinality] += 1
                counts[status] += 1
                writer.writerow({
                    "dataset_id": dataset_id,
                    "species": species,
                    "source_gene_id": value,
                    "candidate_reference_gene_id": target,
                    "mapping_status": status,
                    "mapping_rule": rule,
                    "candidate_protein_accessions": "|".join(proteins),
                    "protein_cardinality": cardinality,
                    "reference_sha256": reference_hashes[species],
                    "accepted": False,
                })
            collision_targets = sum(count > 1 for count in targets.values())
            mapped = counts["exact_candidate"] + counts["deterministic_transform_candidate"]
            summaries.append({
                "dataset_id": dataset_id,
                "species": species,
                "features": len(values),
                "mapped_candidates": mapped,
                "mapped_fraction": mapped / len(values) if values else 0,
                "exact_candidates": counts["exact_candidate"],
                "transform_candidates": counts["deterministic_transform_candidate"],
                "ambiguous_candidates": counts["ambiguous_candidate"],
                "unmapped": counts["unmapped"],
                "target_collision_count": collision_targets,
                "gene_rows_with_multiple_proteins": sum(n for card, n in protein_cardinalities.items() if card > 1),
                "accepted": False,
            })
            print(dataset_id, mapped, len(values), flush=True)
    partial.replace(output)

    summary = pd.DataFrame(summaries)
    summary.to_csv(project / "metadata" / "gene_mapping_candidate_summary.tsv", sep="\t", index=False)
    payload = {
        "species": int(summary.species.nunique()),
        "datasets": len(summary),
        "mapping_rows": int(summary.features.sum()),
        "mapped_candidate_rows": int(summary.mapped_candidates.sum()),
        "weighted_mapped_fraction": float(summary.mapped_candidates.sum() / summary.features.sum()),
        "ambiguous_candidate_rows": int(summary.ambiguous_candidates.sum()),
        "target_collisions_across_source_ids": int(summary.target_collision_count.sum()),
        "gene_rows_with_multiple_proteins": int(summary.gene_rows_with_multiple_proteins.sum()),
        "mapping_file": str(output),
        "mapping_file_sha256": sha256(output),
        "accepted": False,
        "source_h5ad_modified": False,
    }
    (project / "reports" / "gene_mapping_candidate_audit.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
