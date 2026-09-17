#!/usr/bin/env python3
"""Map dominant marker genes to OrthoFinder resolved trees using a random-choice compactness null."""

from __future__ import annotations

import argparse
import json
import os
import random
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    leaves: list["Node"] = field(default_factory=list)


def parse_newick(text: str) -> Node:
    i = 0

    def token() -> str:
        nonlocal i
        start = i
        while i < len(text) and text[i] not in ",():;":
            i += 1
        return text[start:i].strip()

    def branch() -> float:
        nonlocal i
        if i < len(text) and text[i] == ":":
            i += 1
            start = i
            while i < len(text) and text[i] not in ",();":
                i += 1
            try:
                return float(text[start:i])
            except ValueError:
                return 0.0
        return 0.0

    def subtree() -> Node:
        nonlocal i
        if text[i] == "(":
            i += 1
            children = [subtree()]
            while text[i] == ",":
                i += 1
                children.append(subtree())
            if text[i] != ")":
                raise ValueError(f"expected ')' at {i}")
            i += 1
            node = Node(name=token(), children=children)
            for child in children:
                child.parent = node
            node.length = branch()
            return node
        node = Node(name=token())
        node.length = branch()
        return node

    root = subtree()
    return root


def populate_leaves(node: Node) -> list[Node]:
    if not node.children:
        node.leaves = [node]
    else:
        node.leaves = [leaf for child in node.children for leaf in populate_leaves(child)]
    return node.leaves


def ancestors(node: Node) -> list[Node]:
    result = []
    while node is not None:
        result.append(node)
        node = node.parent
    return result


def mrca(nodes: list[Node]) -> Node:
    common = set(ancestors(nodes[0]))
    for node in nodes[1:]:
        common &= set(ancestors(node))
    for node in ancestors(nodes[0]):
        if node in common:
            return node
    raise RuntimeError("no common ancestor")


def distance(a: Node, b: Node) -> float:
    aa = ancestors(a)
    bb = set(ancestors(b))
    common = next(node for node in aa if node in bb)
    total = 0.0
    node = a
    while node is not common:
        total += node.length
        node = node.parent
    node = b
    while node is not common:
        total += node.length
        node = node.parent
    return total


def mean_pairwise(nodes: list[Node]) -> float:
    values = [distance(nodes[i], nodes[j]) for i in range(len(nodes)) for j in range(i + 1, len(nodes))]
    return float(np.mean(values)) if values else 0.0


def leaf_identity(name: str) -> tuple[str, str]:
    marker = "_representative_proteins_"
    if marker not in name:
        return "", ""
    species, remainder = name.split(marker, 1)
    gene = remainder.split("|", 1)[0]
    return species.replace("_", " ").capitalize(), gene


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


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--project-dir", type=Path, required=True)
    p.add_argument("--resolved-trees", type=Path, required=True)
    p.add_argument("--version", default="v2")
    p.add_argument("--random-draws", type=int, default=5000)
    p.add_argument("--seed", type=int, default=20260911)
    args = p.parse_args()
    summary_path = args.project_dir / "metadata" / f"priority_paralog_orthogroup_summary_{args.version}.tsv"
    summary = pd.read_csv(summary_path, sep="\t", keep_default_na=False)
    target_ogs = set(summary["orthogroup_id"])
    tree_text = {}
    with args.resolved_trees.open("r", encoding="utf-8") as handle:
        for line in handle:
            og, newick = line.rstrip("\n").split(": ", 1)
            if og in target_ogs:
                tree_text[og] = newick

    rows = []
    for row in summary.itertuples(index=False):
        root = parse_newick(tree_text[row.orthogroup_id])
        leaves = populate_leaves(root)
        by_species: dict[str, list[Node]] = {}
        by_gene: dict[tuple[str, str], Node] = {}
        for leaf in leaves:
            species, gene = leaf_identity(leaf.name)
            if species and gene:
                by_species.setdefault(species, []).append(leaf)
                by_gene[(species, gene)] = leaf
        dominant = {}
        for block in row.dominant_genes_by_species.split(" | "):
            species, gene = block.split(": ", 1)
            dominant[species] = gene
        missing = [f"{species}:{gene}" for species, gene in dominant.items() if (species, gene) not in by_gene]
        dominant_nodes = [by_gene[(species, gene)] for species, gene in dominant.items() if (species, gene) in by_gene]
        if len(dominant_nodes) < 2:
            observed = len(leaves)
            null_sizes = [len(leaves)]
            observed_distance = 0.0
            null_distances = [0.0]
        else:
            observed = len(mrca(dominant_nodes).leaves)
            observed_distance = mean_pairwise(dominant_nodes)
            rng = random.Random(f"{args.seed}|{row.orthogroup_id}|{row.panel}|{row.coarse_label}")
            species_list = [species for species in dominant if species in by_species]
            null_sizes = []
            null_distances = []
            for _ in range(args.random_draws):
                sampled = [rng.choice(by_species[species]) for species in species_list]
                null_sizes.append(len(mrca(sampled).leaves))
                null_distances.append(mean_pairwise(sampled))
        percentile = float(np.mean(np.asarray(null_distances) <= observed_distance))
        q10, median, q90 = (float(x) for x in np.quantile(null_sizes, [0.10, 0.50, 0.90]))
        d10, dmedian, d90 = (float(x) for x in np.quantile(null_distances, [0.10, 0.50, 0.90]))
        multigene_species = int(row.species_with_multi_gene_family_and_positive_driver)
        if not missing and percentile <= 0.10:
            pattern = "conserved_subclade_signal"
        elif not missing and percentile >= 0.90 and multigene_species >= 2:
            pattern = "dispersed_drivers_possible_paralog_switching"
        else:
            pattern = "indeterminate_tree_pattern"
        rows.append({
            "panel": row.panel, "coarse_label": row.coarse_label, "orthogroup_id": row.orthogroup_id,
            "dominant_species_count": len(dominant), "tree_leaf_count": len(leaves),
            "species_with_multi_gene_family_and_positive_driver": multigene_species,
            "observed_dominant_gene_mrca_leaf_count": observed,
            "random_mrca_leaf_count_q10": q10, "random_mrca_leaf_count_median": median,
            "random_mrca_leaf_count_q90": q90,
            "observed_mean_pairwise_tree_distance": observed_distance,
            "random_mean_pairwise_distance_q10": d10,
            "random_mean_pairwise_distance_median": dmedian,
            "random_mean_pairwise_distance_q90": d90,
            "distance_percentile_fraction_random_le_observed": percentile,
            "gene_tree_pattern": pattern, "missing_dominant_tree_leaves": ";".join(missing),
            "dominant_genes_by_species": row.dominant_genes_by_species,
            "paralog_substitution_claim": "not_established",
        })
    result = pd.DataFrame(rows).sort_values(
        ["gene_tree_pattern", "distance_percentile_fraction_random_le_observed", "panel", "orthogroup_id"]
    )
    output = args.project_dir / "metadata" / f"priority_paralog_gene_tree_evidence_{args.version}.tsv"
    atomic_tsv(output, result)
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "exploratory_gene_tree_mapping_complete",
        "orthogroups": len(result), "random_draws_per_orthogroup": args.random_draws,
        "patterns": result["gene_tree_pattern"].value_counts().to_dict(),
        "missing_dominant_tree_leaves": int(result["missing_dominant_tree_leaves"].ne("").sum()),
        "output": str(output),
        "interpretation": {
            "conserved_subclade_signal": "Dominant genes have a lower mean pairwise tree distance than at least 90% of random one-gene-per-species selections.",
            "dispersed_drivers_possible_paralog_switching": "Dominant genes have a higher mean pairwise tree distance than at least 90% of random selections; this prioritizes, but does not prove, paralog switching.",
            "indeterminate_tree_pattern": "Current dominant genes do not meet either exploratory compactness threshold.",
        },
        "limitations": [
            "The random-choice compactness screen is exploratory and depends on the OrthoFinder resolved tree.",
            "Gene-tree dispersion can arise from tree uncertainty, annotation errors, lineage-specific duplication or incomplete sampling.",
            "No paralog substitution claim is accepted without sequence-tree review and independent expression evidence.",
        ],
        "safety": {"source_h5ad_accessed": False, "model_training_started": False, "external_queries": False},
    }
    atomic_json(args.project_dir / "reports" / f"priority_paralog_gene_tree_evidence_{args.version}_audit.json", audit)
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
