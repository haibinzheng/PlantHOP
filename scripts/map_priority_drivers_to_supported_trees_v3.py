#!/usr/bin/env python3
"""Map dominant marker genes to clean IQ-TREE priority-family trees.

This repeats the exploratory random one-gene-per-species distance screen on the
versioned MAFFT/IQ-TREE rebuild. It records topology changes relative to the
superseded OrthoFinder resolved-tree screen and never establishes a paralog
substitution claim.
"""

from __future__ import annotations

import csv
import json
import random
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(".")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"
TREE_ROOT = Path("data/runs/PhyloOpenCell/priority_gene_trees_v2")
ORTHOGROUPS = Path("data/runs/PhyloOpenCell/orthofinder_full7_v1/results/Results_Sep10/Orthogroups/Orthogroups.tsv")


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    leaves: list["Node"] = field(default_factory=list)


def parse_newick(text: str):
    i = 0

    def token():
        nonlocal i
        start = i
        while i < len(text) and text[i] not in ",():;":
            i += 1
        return text[start:i].strip()

    def branch():
        nonlocal i
        if i < len(text) and text[i] == ":":
            i += 1
            start = i
            while i < len(text) and text[i] not in ",();":
                i += 1
            return float(text[start:i])
        return 0.0

    def subtree():
        nonlocal i
        if text[i] == "(":
            i += 1
            children = [subtree()]
            while text[i] == ",":
                i += 1
                children.append(subtree())
            if text[i] != ")":
                raise ValueError(f"expected ')' at character {i}")
            i += 1
            node = Node(name=token(), children=children)
            for child in children:
                child.parent = node
            node.length = branch()
            return node
        node = Node(name=token())
        node.length = branch()
        return node

    return subtree()


def populate(node):
    node.leaves = [node] if not node.children else [leaf for child in node.children for leaf in populate(child)]
    return node.leaves


def ancestors(node):
    result = []
    while node is not None:
        result.append(node)
        node = node.parent
    return result


def mrca(nodes):
    common = set(ancestors(nodes[0]))
    for node in nodes[1:]:
        common &= set(ancestors(node))
    return next(node for node in ancestors(nodes[0]) if node in common)


def distance(a, b):
    path_a = ancestors(a)
    common = next(node for node in path_a if node in set(ancestors(b)))
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


def mean_pairwise(nodes):
    values = [distance(nodes[i], nodes[j]) for i in range(len(nodes)) for j in range(i + 1, len(nodes))]
    return statistics.mean(values) if values else 0.0


def quantile(values, q):
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = (len(ordered) - 1) * q
    lo = int(position)
    hi = min(lo + 1, len(ordered) - 1)
    fraction = position - lo
    return ordered[lo] * (1 - fraction) + ordered[hi] * fraction


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def orthogroup_gene_species(targets):
    result = {}
    with ORTHOGROUPS.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        species_columns = [column for column in reader.fieldnames if column != "Orthogroup"]
        for row in reader:
            og = row["Orthogroup"]
            if og not in targets:
                continue
            mapping = {}
            for column in species_columns:
                species = column.split(".representative_proteins", 1)[0].replace("_", " ").capitalize()
                for member in row[column].split(","):
                    member = member.strip()
                    if member:
                        mapping[member.split("|", 1)[0]] = species
            result[og] = mapping
    return result


def parse_support(label):
    if "/" not in label:
        return "", ""
    left, right = label.split("/", 1)
    try:
        return f"{float(left):.6g}", f"{float(right):.6g}"
    except ValueError:
        return "", ""


def main():
    summary = read_tsv(META / "priority_paralog_orthogroup_summary_v2.tsv")
    old = {r["orthogroup_id"]: r for r in read_tsv(META / "priority_paralog_gene_tree_evidence_v2.tsv")}
    species_by_og = orthogroup_gene_species({r["orthogroup_id"] for r in summary})
    rows = []
    for record in summary:
        og = record["orthogroup_id"]
        root = parse_newick((TREE_ROOT / f"trees/{og}.treefile").read_text(encoding="utf-8").strip())
        leaves = populate(root)
        species_map = species_by_og[og]
        by_species = defaultdict(list)
        by_gene = {}
        for leaf in leaves:
            gene = leaf.name.split("|", 1)[0]
            species = species_map[gene]
            by_species[species].append(leaf)
            by_gene[(species, gene)] = leaf
        dominant = {}
        for block in record["dominant_genes_by_species"].split(" | "):
            species, gene = block.split(": ", 1)
            dominant[species] = gene
        missing = [f"{species}:{gene}" for species, gene in dominant.items() if (species, gene) not in by_gene]
        dominant_nodes = [by_gene[(species, gene)] for species, gene in dominant.items() if (species, gene) in by_gene]
        observed_node = mrca(dominant_nodes)
        observed_distance = mean_pairwise(dominant_nodes)
        rng = random.Random(f"20260911|supported_v3|{og}|{record['panel']}|{record['coarse_label']}")
        species_list = [species for species in dominant if species in by_species]
        null_distances = []
        null_sizes = []
        for _ in range(5000):
            sampled = [rng.choice(by_species[species]) for species in species_list]
            null_distances.append(mean_pairwise(sampled))
            null_sizes.append(len(mrca(sampled).leaves))
        percentile = sum(v <= observed_distance for v in null_distances) / len(null_distances)
        multigene_species = int(record["species_with_multi_gene_family_and_positive_driver"])
        sh, ufboot = parse_support(observed_node.name)
        sh_value = float(sh) if sh else None
        ufboot_value = float(ufboot) if ufboot else None
        if (
            not missing and percentile <= 0.10 and len(observed_node.leaves) < len(leaves)
            and sh_value is not None and ufboot_value is not None
            and sh_value >= 80 and ufboot_value >= 80
        ):
            pattern = "well_supported_dominant_subclade"
        elif not missing and percentile <= 0.10:
            pattern = "compact_driver_distance_only"
        elif not missing and percentile >= 0.90 and multigene_species >= 2:
            pattern = "dispersed_driver_distance_signal"
        else:
            pattern = "indeterminate_tree_pattern"
        old_pattern = old[og]["gene_tree_pattern"]
        old_equivalent = {
            "conserved_subclade_signal": "compact",
            "dispersed_drivers_possible_paralog_switching": "dispersed",
            "indeterminate_tree_pattern": "indeterminate",
        }[old_pattern]
        new_equivalent = (
            "compact" if pattern in {"well_supported_dominant_subclade", "compact_driver_distance_only"}
            else "dispersed" if pattern == "dispersed_driver_distance_signal"
            else "indeterminate"
        )
        rows.append({
            "panel": record["panel"], "coarse_label": record["coarse_label"], "orthogroup_id": og,
            "dominant_species_count": str(len(dominant)), "tree_leaf_count": str(len(leaves)),
            "species_with_multi_gene_family_and_positive_driver": str(multigene_species),
            "observed_dominant_gene_mrca_leaf_count": str(len(observed_node.leaves)),
            "observed_mrca_sh_alrt": sh, "observed_mrca_ultrafast_bootstrap": ufboot,
            "random_mrca_leaf_count_q10": f"{quantile(null_sizes, .1):.9g}",
            "random_mrca_leaf_count_median": f"{quantile(null_sizes, .5):.9g}",
            "random_mrca_leaf_count_q90": f"{quantile(null_sizes, .9):.9g}",
            "observed_mean_pairwise_tree_distance": f"{observed_distance:.12g}",
            "random_mean_pairwise_distance_q10": f"{quantile(null_distances, .1):.12g}",
            "random_mean_pairwise_distance_median": f"{quantile(null_distances, .5):.12g}",
            "random_mean_pairwise_distance_q90": f"{quantile(null_distances, .9):.12g}",
            "distance_percentile_fraction_random_le_observed": f"{percentile:.6g}",
            "supported_tree_pattern": pattern, "superseded_v2_pattern": old_pattern,
            "pattern_stable_after_clean_rebuild": str(new_equivalent == old_equivalent),
            "missing_dominant_tree_leaves": ";".join(missing),
            "dominant_genes_by_species": record["dominant_genes_by_species"],
            "paralog_substitution_claim": "not_established",
        })

    rows.sort(key=lambda x: (x["supported_tree_pattern"], float(x["distance_percentile_fraction_random_le_observed"]), x["orthogroup_id"]))
    output = META / "priority_paralog_supported_tree_evidence_v3.tsv"
    write_tsv(output, rows)
    patterns = Counter(r["supported_tree_pattern"] for r in rows)
    transitions = Counter(f"{r['superseded_v2_pattern']}->{r['supported_tree_pattern']}" for r in rows)
    stable = sum(r["pattern_stable_after_clean_rebuild"] == "True" for r in rows)
    selected = [r for r in rows if r["orthogroup_id"] in {"OG0000114", "OG0000378", "OG0003912"}]
    robustness = {r["orthogroup_id"]: r for r in read_tsv(META / "priority_family_study_robustness_summary_v1.tsv")}
    functions = {r["orthogroup_id"]: r for r in read_tsv(META / "conserved_program_function_evidence_v2.tsv")}
    joint_rows = []
    for tree_row in rows:
        og = tree_row["orthogroup_id"]
        study = robustness[og]
        annotation = functions[og]
        tree_pattern = tree_row["supported_tree_pattern"]
        study_status = study["cross_study_family_status"]
        if tree_pattern == "well_supported_dominant_subclade" and study_status == "cross_species_study_robust":
            joint_status = "joint_tree_and_expression_priority"
        elif tree_pattern == "well_supported_dominant_subclade":
            joint_status = "supported_subclade_expression_not_robust"
        elif tree_pattern == "dispersed_driver_distance_signal" and study_status == "cross_species_study_robust":
            joint_status = "expression_robust_dispersed_hypothesis"
        elif tree_pattern == "compact_driver_distance_only" and study_status in {"cross_species_study_robust", "cross_species_mixed_support"}:
            joint_status = "expression_supported_compact_distance_only"
        elif study_status == "cross_species_study_robust":
            joint_status = "expression_robust_tree_indeterminate"
        else:
            joint_status = "preliminary_or_insufficient_joint_evidence"
        joint_rows.append({
            "panel": tree_row["panel"], "coarse_label": tree_row["coarse_label"], "orthogroup_id": og,
            "local_symbols": annotation["local_symbols_all_species"],
            "local_descriptions": annotation["local_descriptions"],
            "study_robustness_status": study_status,
            "heldout_tests": study["heldout_tests"],
            "heldout_positive_fraction": study["heldout_positive_fraction"],
            "supported_tree_pattern": tree_pattern,
            "distance_percentile": tree_row["distance_percentile_fraction_random_le_observed"],
            "mrca_sh_alrt": tree_row["observed_mrca_sh_alrt"],
            "mrca_ultrafast_bootstrap": tree_row["observed_mrca_ultrafast_bootstrap"],
            "joint_evidence_status": joint_status,
            "paralog_substitution_claim": "not_established",
        })
    joint_rows.sort(key=lambda r: (r["joint_evidence_status"], r["panel"], r["orthogroup_id"]))
    joint_output = META / "priority_family_joint_expression_tree_evidence_v1.tsv"
    write_tsv(joint_output, joint_rows)
    joint_counts = Counter(r["joint_evidence_status"] for r in joint_rows)
    matrix_v3 = read_tsv(META / "jeb_biology_evidence_matrix_v3.tsv")
    matrix_v4 = []
    for row in matrix_v3:
        updated = dict(row)
        if row["claim_id"] == "C6":
            updated["proposed_claim"] = "One family has a well-supported dominant-gene subclade, but no priority family yet passes both the tree and cross-study expression gates."
            updated["current_status"] = "exploratory_not_main_claim"
            updated["supporting_evidence"] = "OG0000299 (BCA3/BCA4/CA1 carbonic anhydrase family) forms a compact subclade with SH-aLRT/UFBoot 98.2/98, but its leave-one-study positive fraction is 0.5625 and it is not cross-species study-robust. AAP OG0000378 is compact by distance only and is not a monophyletic dominant subclade."
            updated["claim_boundary"] = "No family currently combines a well-supported dominant subclade with cross-species study-robust expression. LTP OG0000114 became indeterminate after clean reconstruction."
            updated["primary_source"] = "metadata/priority_family_joint_expression_tree_evidence_v1.tsv; metadata/priority_paralog_supported_tree_evidence_v3.tsv; metadata/priority_family_study_robustness_summary_v1.tsv"
            updated["recommended_placement"] = "Supplementary exploratory analysis; do not use as a main mechanistic claim."
            updated["next_required_evidence"] = "Independent expression replication for OG0000299 and greater multi-study species coverage for AAP are required before promotion."
        elif row["claim_id"] == "C7":
            updated["proposed_claim"] = "PSBO (OG0003912) and OG0002909 show dispersed dominant-gene distances and are candidates for evolutionary follow-up."
            updated["current_status"] = "hypothesis_only_after_clean_tree_rebuild"
            updated["supporting_evidence"] = "OG0003912 retained a high-distance-tail signal (percentile 0.9102) after clean reconstruction and is internally study-robust; OG0002909 also entered the high-distance tail (0.9358)."
            updated["claim_boundary"] = "Neither result establishes paralog switching: both dominant sets span the family tree, and dispersion may reflect duplication history, annotation uncertainty, or incomplete species/study coverage."
            updated["primary_source"] = "metadata/priority_paralog_supported_tree_evidence_v3.tsv; metadata/priority_family_study_robustness_summary_v1.tsv"
            updated["recommended_placement"] = "Supplementary hypothesis dossier; main text only as a future-testable observation."
            updated["next_required_evidence"] = "Curated duplication-node interpretation, independent expression replication, and reciprocal paralog patterns are required before any switching terminology."
        matrix_v4.append(updated)
    matrix_v4_path = META / "jeb_biology_evidence_matrix_v4.tsv"
    write_tsv(matrix_v4_path, matrix_v4)
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "clean_supported_tree_mapping_complete",
        "orthogroups": len(rows),
        "random_draws_per_orthogroup": 5000,
        "patterns": dict(patterns),
        "patterns_stable_vs_superseded_v2": stable,
        "pattern_transitions": dict(transitions),
        "missing_dominant_tree_leaves": sum(bool(r["missing_dominant_tree_leaves"]) for r in rows),
        "selected_family_results": selected,
        "joint_evidence_status_counts": dict(joint_counts),
        "joint_tree_and_expression_priority_families": [r["orthogroup_id"] for r in joint_rows if r["joint_evidence_status"] == "joint_tree_and_expression_priority"],
        "interpretation_boundary": "Clean supported trees restore auditable topology provenance. The compactness screen remains exploratory; neither dispersion nor different dominant genes establishes paralog substitution.",
        "input_run": str(TREE_ROOT),
        "output": str(output),
        "joint_evidence_output": str(joint_output),
        "updated_evidence_matrix": str(matrix_v4_path),
        "safety": {"source_h5ad_accessed": False, "model_training_started": False, "external_queries": False, "original_v1_modified": False},
    }
    audit_path = REPORTS / "priority_paralog_supported_tree_evidence_v3_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
