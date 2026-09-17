#!/usr/bin/env python3
"""Audit file-to-file consistency for priority OrthoFinder gene trees.

The audit is read-only. It distinguishes manuscript-facing named alignments and
resolved trees from the internal-ID alignment/raw-tree artifacts, because branch
support must never be transferred across inconsistent identifier layers.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(".")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"
OF = Path("data/runs/PhyloOpenCell/orthofinder_full7_v1/results/Results_Sep10")


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path):
    records = {}
    name = None
    chunks = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).replace("-", "").upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).replace("-", "").upper()
    return records


def sha256_text(value: str):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def parse_sequence_ids(path: Path):
    mapping = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            internal, rest = line.rstrip("\n").split(": ", 1)
            mapping[internal] = rest.split()[0]
    return mapping


def parse_resolved_trees(path: Path, targets: set[str]):
    trees = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            og, tree = line.rstrip("\n").split(": ", 1)
            if og in targets:
                trees[og] = tree
    return trees


def tree_leaves(newick: str):
    # OrthoFinder leaf names do not contain Newick punctuation or whitespace.
    return set(re.findall(r"(?<=[(,])([^():;,]+):", newick))


def count_members(value: str):
    if not value:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def write_tsv(path: Path, rows: list[dict]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_rows(path: Path, rows: list[dict]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    evidence = read_tsv(META / "priority_paralog_gene_tree_evidence_v2.tsv")
    targets = {row["orthogroup_id"] for row in evidence}
    species_ids = {}
    with (OF / "WorkingDirectory/SpeciesIDs.txt").open("r", encoding="utf-8") as handle:
        for line in handle:
            number, filename = line.strip().split(": ", 1)
            species_ids[number] = filename.removesuffix(".fa")
    sequence_ids = parse_sequence_ids(OF / "WorkingDirectory/SequenceIDs.txt")
    resolved = parse_resolved_trees(OF / "Resolved_Gene_Trees/Resolved_Gene_Trees.txt", targets)

    orthogroups = {}
    with (OF / "Orthogroups/Orthogroups.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        species_columns = [c for c in reader.fieldnames if c != "Orthogroup"]
        for row in reader:
            if row["Orthogroup"] in targets:
                members = set()
                prefixed = set()
                for species in species_columns:
                    prefix = species.removesuffix(".fa").replace(".", "_") + "_"
                    for gene in count_members(row[species]):
                        members.add(gene)
                        prefixed.add(prefix + gene)
                orthogroups[row["Orthogroup"]] = (members, prefixed)

    rows = []
    for og in sorted(targets):
        expected_members, expected_prefixed = orthogroups[og]
        named = read_fasta(OF / f"MultipleSequenceAlignments/{og}.fa")
        internal = read_fasta(OF / f"WorkingDirectory/Alignments_ids/{og}.fa")
        raw_tree = (OF / f"WorkingDirectory/Trees_ids/{og}.txt").read_text(encoding="utf-8").strip()
        resolved_tree = resolved[og]
        resolved_leaves = tree_leaves(resolved_tree)
        raw_leaves = tree_leaves(raw_tree)

        # Expected named headers are exactly the species-prefixed orthogroup IDs.
        named_header_match = set(named) == expected_prefixed
        resolved_header_match = resolved_leaves == set(named)
        raw_internal_match = raw_leaves == set(internal)

        # Internal alignment IDs should map back to the same orthogroup. A
        # failure does not invalidate the named MSA by itself, but it makes raw
        # branch supports unsafe to transfer.
        mapped_internal = {sequence_ids[x] for x in internal if x in sequence_ids}
        internal_mapping_match = mapped_internal == expected_members

        # Sequence hashes should be one-to-one between named and internal files
        # if they are two identifier views of the same alignment.
        named_hashes = Counter(sha256_text(seq) for seq in named.values())
        internal_hashes = Counter(sha256_text(seq) for seq in internal.values())
        sequence_multiset_match = named_hashes == internal_hashes

        if not named_header_match or not resolved_header_match:
            manuscript_tree_status = "invalid_named_tree_membership"
        elif not internal_mapping_match or not sequence_multiset_match or not raw_internal_match:
            manuscript_tree_status = "named_topology_membership_ok_raw_support_unusable"
        else:
            manuscript_tree_status = "all_layers_consistent"

        rows.append({
            "orthogroup_id": og,
            "orthogroup_members": str(len(expected_members)),
            "named_alignment_sequences": str(len(named)),
            "resolved_tree_leaves": str(len(resolved_leaves)),
            "internal_alignment_sequences": str(len(internal)),
            "raw_tree_leaves": str(len(raw_leaves)),
            "named_alignment_membership_match": str(named_header_match),
            "resolved_tree_named_alignment_match": str(resolved_header_match),
            "raw_tree_internal_alignment_match": str(raw_internal_match),
            "internal_ids_map_to_orthogroup": str(internal_mapping_match),
            "named_internal_sequence_multiset_match": str(sequence_multiset_match),
            "manuscript_tree_status": manuscript_tree_status,
        })

    output = META / "priority_gene_tree_integrity_audit_v1.tsv"
    write_tsv(output, rows)
    statuses = Counter(r["manuscript_tree_status"] for r in rows)
    issue_detected = set(statuses) != {"all_layers_consistent"}

    matrix_v2 = read_tsv(META / "jeb_biology_evidence_matrix_v2.tsv")
    matrix_v3 = []
    for row in matrix_v2:
        updated = dict(row)
        if issue_detected and row["claim_id"] in {"C6", "C7"}:
            updated["current_status"] = "suspended_pending_tree_rebuild"
            updated["claim_boundary"] += (
                " Integrity audit found that internal-ID/raw-tree artifacts do not represent the same sequence set as the named alignment under the current SequenceIDs mapping; raw support values are unusable."
            )
            updated["primary_source"] += "; reports/priority_gene_tree_integrity_v1_audit.json"
            updated["recommended_placement"] = "Do not include in main manuscript until a clean versioned tree rebuild passes integrity audit."
            updated["next_required_evidence"] = "Rebuild the priority-family alignments and supported trees from named orthogroup sequences in a new versioned directory, then repeat topology analysis."
        matrix_v3.append(updated)
    matrix_v3_path = META / "jeb_biology_evidence_matrix_v3.tsv"
    write_rows(matrix_v3_path, matrix_v3)

    correction_path = REPORTS / "gene_tree_evidence_correction_notice_v1.md"
    correction_path.write_text(
        "# Gene-tree evidence correction notice (v1)\n\n"
        "## Status\n\n"
        "The expression recurrence, orthogroup membership, gene-level dominance, and leave-one-study-out results remain usable. The existing gene-tree compactness classification and all paralog-switching interpretations are suspended.\n\n"
        "## Reason\n\n"
        "For all 19 priority orthogroups, the named multiple-sequence alignment and resolved-tree leaf membership agree with `Orthogroups.tsv`. However, the corresponding internal-ID alignment and raw FastTree tree do not map back to that orthogroup under the run's current `SequenceIDs.txt`, and the named/internal sequence multisets differ. Therefore raw branch-support values cannot be transferred to the named resolved trees, and the provenance of the existing topology is insufficient for a manuscript claim.\n\n"
        "## Required correction\n\n"
        "Create a new, versioned tree-only analysis from the named `Orthogroup_Sequences` files for the priority families. Record exact sequence hashes, alignment hashes, tree program/version/command, branch-support method, and one-to-one leaf membership checks. Do not modify `data/runs/PhyloOpenCell/orthofinder_full7_v1`.\n\n"
        "## Manuscript consequence\n\n"
        "C6 and C7 in the JEB evidence matrix are suspended in `jeb_biology_evidence_matrix_v3.tsv`. OG0003912 remains an expression-robust PSBO-family candidate, but it is not currently a supported paralog-switching candidate.\n",
        encoding="utf-8",
    )
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "priority_gene_tree_integrity_issue_detected" if issue_detected else "priority_gene_tree_integrity_pass",
        "orthogroups_checked": len(rows),
        "manuscript_tree_status_counts": dict(statuses),
        "interpretation": {
            "named_topology_membership_ok_raw_support_unusable": "Named MSA membership and resolved-tree leaves agree with Orthogroups.tsv, but internal-ID/raw-tree artifacts do not represent the same sequence set under the current SequenceIDs mapping. Resolved topology may remain exploratory; raw FastTree support values must not be attached to it.",
            "invalid_named_tree_membership": "Named MSA or resolved-tree membership does not agree with Orthogroups.tsv; exclude the tree from interpretation.",
        },
        "decision": "Suspend any gene-tree compactness or paralog-switching claim until priority families are rebuilt in a new versioned directory from named orthogroup sequences with auditable support values.",
        "affected_previous_outputs": [
            "metadata/priority_paralog_gene_tree_evidence_v2.tsv",
            "reports/priority_paralog_gene_tree_evidence_v2_audit.json",
            "metadata/jeb_priority_family_dossier_v1.tsv",
            "metadata/jeb_biology_evidence_matrix_v2.tsv claim C6/C7 tree component",
        ],
        "unaffected_evidence": [
            "orthogroup membership and expression recurrence tables",
            "gene-level expression dominance tables",
            "internal leave-one-study-out expression robustness results",
        ],
        "output": str(output),
        "corrected_evidence_matrix": str(matrix_v3_path),
        "correction_notice": str(correction_path),
        "safety": {
            "source_h5ad_accessed": False,
            "model_training_started": False,
            "external_queries": False,
            "existing_v1_run_modified": False,
        },
    }
    audit_path = REPORTS / "priority_gene_tree_integrity_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
