#!/usr/bin/env python3
"""Package the completed P1 sensitivity and evidence-gap results as Supplementary Table S7."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path("/workspace/projects/PhyloOpenCell")
SOURCE = PROJECT / "reports/discovery_threshold_ortholog_robustness_v1"
OUTPUT = PROJECT / "reports/manuscript_supplement_table_s7_v1"
TABLES = {
    "Table_S7a_discovery_threshold_sensitivity.tsv": SOURCE / "discovery_threshold_sensitivity.tsv",
    "Table_S7b_clade_balance_eligibility.tsv": SOURCE / "class_clade_eligibility.tsv",
    "Table_S7c_xylem_program_rediscovery.tsv": SOURCE / "xylem_program_discovery_ortholog_gap.tsv",
    "Table_S7d_single_gene_ortholog_evidence_gaps.tsv": SOURCE / "priority_family_single_gene_ortholog_evidence_gap.tsv",
}
AUTHORITY = SOURCE / "audit.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if not rows:
        raise RuntimeError(f"empty output: {path}")
    names = fields or list(rows[0])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


DESCRIPTIONS = {
    "top_n": "Within-species positive-effect rank cutoff.",
    "minimum_species_recurrence": "Minimum number of recurrent species required for candidate selection.",
    "scope": "All panels or one named anatomical panel.",
    "candidate_count": "Candidates meeting the top-N and species-recurrence condition.",
    "cross_clade_presence_count": "Selected candidates present in at least one eudicot and one monocot.",
    "strict_2plus2_clade_balance_count": "Selected candidates present in at least two eudicots and two monocots.",
    "frozen_102_retained_count": "Members of the frozen top-100/minimum-3 candidate set retained.",
    "frozen_102_retained_fraction": "Fraction of the complete frozen 102-candidate reference set retained.",
    "panel": "Frozen anatomical panel.", "coarse_label": "Frozen coarse cell-label group.",
    "eligible_species_count": "Species represented for this panel and label.", "eligible_species": "Semicolon-separated represented species.",
    "eudicot_count": "Number of represented eudicot species.", "monocot_count": "Number of represented monocot species.",
    "cross_clade_evaluable": "Whether both eudicots and monocots are represented.",
    "strict_2plus2_evaluable": "Whether at least two eudicots and two monocots are represented.",
    "frozen_rank": "Rank in the submitted frozen 18-orthogroup xylem program.", "orthogroup_id": "Frozen orthogroup identifier.",
    "mapped_in_all_three_validation_species": "Whether the program member maps in all three GSE268881 species assets.",
    "top25_species_count": "Species recurring within their top 25 vascular-xylem markers.",
    "top50_species_count": "Species recurring within their top 50 vascular-xylem markers.",
    "top100_species_count": "Species recurring within their top 100 vascular-xylem markers.",
    "top100_species": "Species supporting top-100 vascular-xylem rediscovery.",
    "top100_eudicot_count": "Rediscovery-supporting eudicot count.", "top100_monocot_count": "Rediscovery-supporting monocot count.",
    "selected_top100_min3": "Whether independently rediscovered under the frozen top-100/minimum-3 rule.",
    "cross_clade_presence_top100": "Whether top-100 recurrence includes a eudicot and a monocot.",
    "single_gene_representatives": "Deterministic representatives where available; not proof of one-to-one orthology.",
    "annotation_class": "Availability class of functional annotation.", "function_symbols": "Available gene or protein symbols.",
    "evidence_scope": "Maximum scope of the mapped/annotation evidence.", "gap_flags": "Semicolon-separated unresolved evidence gaps.",
    "allowed_claim": "Maximum claim supported by the row.",
    "discovery_species_count_top100": "Species count supporting the priority family in frozen top-100 discovery.",
    "cross_study_family_status": "Internal leave-one-study-out expression robustness category.",
    "heldout_tests": "Number of evaluable held-out study tests.", "heldout_positive_fraction": "Fraction of held-out tests with positive direction.",
    "species_with_positive_gene_level_marker": "Species containing a positive gene-level family member.",
    "species_with_single_dominant_pattern": "Species classified as having a single-dominant expression pattern.",
    "species_with_dominance_rows": "Species with gene-dominance evidence rows.",
    "dominant_genes_by_species": "Descriptive dominant genes by species; not proof of orthology or substitution.",
    "supported_tree_pattern": "Conservative pattern from clean supported-tree analysis.",
    "observed_mrca_sh_alrt": "SH-aLRT support at the observed dominant-gene MRCA, when available.",
    "observed_mrca_ultrafast_bootstrap": "Ultrafast-bootstrap support at the observed dominant-gene MRCA, when available.",
    "joint_evidence_status_existing": "Existing integrated expression/tree evidence category.",
    "conservative_evidence_tier": "Conservative evidence tier assigned from the available annotation and mapping evidence.",
    "one_to_one_ortholog_claim": "Whether one-to-one orthology is established; expected to remain not_established.",
    "paralog_substitution_claim": "Whether paralog substitution is established; expected to remain not_established.",
    "main_gap": "Principal unresolved evidence limitation.",
    "interpretation": "Interpretation boundary for the threshold-sensitivity row.",
    "design_note": "Design limitation relevant to clade-balance interpretation.",
}


def main() -> int:
    if OUTPUT.exists():
        raise RuntimeError(f"refusing to overwrite {OUTPUT}")
    if not AUTHORITY.is_file():
        raise FileNotFoundError(AUTHORITY)
    for path in TABLES.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    authority = json.loads(AUTHORITY.read_text(encoding="utf-8"))
    if authority.get("status") != "COMPLETE_WITH_GAPS" or not all(authority.get("checks", {}).values()):
        raise RuntimeError("Priority-family evidence record is not complete")

    expected = {
        "Table_S7a_discovery_threshold_sensitivity.tsv": 48,
        "Table_S7b_clade_balance_eligibility.tsv": 11,
        "Table_S7c_xylem_program_rediscovery.tsv": 18,
        "Table_S7d_single_gene_ortholog_evidence_gaps.tsv": 19,
    }
    tables = {name: read_tsv(path) for name, path in TABLES.items()}
    actual = {name: len(rows) for name, rows in tables.items()}
    if actual != expected:
        raise RuntimeError(f"unexpected source row counts: {actual}")
    if any(r["one_to_one_ortholog_claim"] != "not_established" for r in tables["Table_S7d_single_gene_ortholog_evidence_gaps.tsv"]):
        raise RuntimeError("unexpected one-to-one ortholog claim")
    if any(r["paralog_substitution_claim"] != "not_established" for r in tables["Table_S7d_single_gene_ortholog_evidence_gaps.tsv"]):
        raise RuntimeError("unexpected paralog-substitution claim")

    OUTPUT.mkdir(parents=False)
    (OUTPUT / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    for name, rows in tables.items():
        write_tsv(OUTPUT / name, rows)

    dictionary_rows = []
    for table_name, rows in tables.items():
        for field in rows[0]:
            dictionary_rows.append({
                "table": table_name, "field": field,
                "description": DESCRIPTIONS.get(field, field.replace("_", " ").capitalize() + "."),
                "missing_value_policy": "Blank means unavailable or not applicable; never interpret as zero.",
            })
    write_tsv(OUTPUT / "field_dictionary.tsv", dictionary_rows, ["table", "field", "description", "missing_value_policy"])

    sources = [{"source_id": "p1_authority_audit", "path": str(AUTHORITY), "sha256": sha256(AUTHORITY), "bytes": AUTHORITY.stat().st_size}]
    for name, path in TABLES.items():
        sources.append({"source_id": name, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size})
    write_tsv(OUTPUT / "input_source_manifest.tsv", sources)

    readme = [
        "# Supplementary Table S7", "",
        "This submission package has four linked TSV parts:",
        "- S7a: discovery sensitivity over top-N and species-recurrence settings (48 rows).",
        "- S7b: clade-balance eligibility by panel and label (11 rows).",
        "- S7c: independent rediscovery and evidence gaps for the submitted 18-OG xylem program (18 rows).",
        "- S7d: integrated single-gene, study-robustness, tree, and ortholog-evidence gaps for 19 priority families (19 rows).", "",
        "The top-100/minimum-3 reference exactly reproduces the frozen 102 candidates. Cross-clade and strict 2+2 fields are descriptive sensitivity fields, not replacement discovery rules.",
        "Only 6/18 submitted xylem-program orthogroups were rediscovered by the frozen vascular-xylem top-100/minimum-3 rule. This is a discovery-overlap result, not evidence against the remaining program members.",
        "No priority family satisfies the strict joint expression-plus-supported-subclade criterion. One-to-one orthology, conserved function, causality, and paralog substitution remain unestablished.",
    ]
    (OUTPUT / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    audit = {
        "status": "COMPLETE_WITH_GAPS", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "version": "manuscript_supplement_table_s7_v1", "table_rows": actual,
        "total_rows": sum(actual.values()),
        "checks": {
            "p1_authority_complete": True, "expected_row_counts": True,
            "one_to_one_ortholog_claims_remain_not_established": True,
            "paralog_substitution_claims_remain_not_established": True,
            "frozen_102_reproduction_inherited": authority["checks"]["frozen_102_exactly_reproduced"],
        },
        "key_counts": authority["counts"],
        "interpretation_boundary": authority["interpretation_boundary"],
        "source_p1_audit_sha256": sha256(AUTHORITY),
        "safety": {"p1_authority_modified": False, "source_data_modified": False, "existing_results_overwritten": False},
    }
    (OUTPUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUTPUT / "STATUS").write_text("COMPLETE_WITH_GAPS\n", encoding="utf-8")
    targets = [p for p in sorted(OUTPUT.iterdir()) if p.is_file() and p.name != "output_sha256.tsv"]
    write_tsv(OUTPUT / "output_sha256.tsv", [{"output_file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in targets])
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
