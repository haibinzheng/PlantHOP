#!/usr/bin/env python3
"""Build a frozen sensitivity and ortholog-evidence gap package from existing TSVs."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(".")
REPORT = PROJECT / "reports/discovery_threshold_ortholog_robustness_v1"
SUPPLEMENT = PROJECT / "reports/manuscript_supplement_tables_v1_1"
STAGED_FREEZE = Path(__file__).with_name("discovery_threshold_ortholog_robustness_v1_freeze_protocol.md")
INPUTS = {
    "detail": PROJECT / "metadata/coarse_marker_orthogroups_v2.tsv",
    "candidates": PROJECT / "metadata/conserved_program_candidates_v2.tsv",
    "study": PROJECT / "metadata/priority_family_study_robustness_summary_v1.tsv",
    "dominance": PROJECT / "metadata/priority_paralog_species_dominance_v2.tsv",
    "orthogroup": PROJECT / "metadata/priority_paralog_orthogroup_summary_v2.tsv",
    "tree": PROJECT / "metadata/priority_paralog_supported_tree_evidence_v3.tsv",
    "joint": PROJECT / "metadata/priority_family_joint_expression_tree_evidence_v1.tsv",
    "xylem_program": SUPPLEMENT / "Table_S4_xylem_orthogroup_program.tsv",
}
TOP_NS = (25, 50, 100)
MIN_SPECIES = (2, 3, 4, 5)
EUDICOTS = {"Arabidopsis thaliana", "Brassica rapa", "Catharanthus roseus", "Glycine max", "Populus trichocarpa"}
MONOCOTS = {"Oryza sativa", "Sorghum bicolor"}


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
        raise RuntimeError(f"refusing empty output: {path}")
    names = fields or list(rows[0])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def key(row: dict[str, str]) -> tuple[str, str, str]:
    return row["panel"], row["coarse_label"], row["orthogroup_id"]


def main() -> int:
    if REPORT.exists():
        raise RuntimeError(f"refusing to overwrite {REPORT}")
    if not STAGED_FREEZE.is_file():
        raise FileNotFoundError(STAGED_FREEZE)
    for path in INPUTS.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    REPORT.mkdir(parents=False)
    (REPORT / "freeze_protocol_v1.md").write_bytes(STAGED_FREEZE.read_bytes())
    (REPORT / Path(__file__).name).write_bytes(Path(__file__).read_bytes())

    detail = read_tsv(INPUTS["detail"])
    frozen = read_tsv(INPUTS["candidates"])
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    class_species: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in detail:
        grouped[key(row)].append(row)
        class_species[(row["panel"], row["coarse_label"])].add(row["species"])

    selections: dict[tuple[int, int], set[tuple[str, str, str]]] = {}
    group_stats: dict[tuple[int, tuple[str, str, str]], dict[str, object]] = {}
    for top_n in TOP_NS:
        for group_key, rows in grouped.items():
            kept = [r for r in rows if int(r["within_species_rank"]) <= top_n]
            species = sorted({r["species"] for r in kept})
            eudicots = sorted(set(species) & EUDICOTS)
            monocots = sorted(set(species) & MONOCOTS)
            group_stats[(top_n, group_key)] = {
                "species": species, "eudicots": eudicots, "monocots": monocots,
                "median_rank": statistics.median(int(r["within_species_rank"]) for r in kept) if kept else "",
            }
        for minimum in MIN_SPECIES:
            selections[(top_n, minimum)] = {
                group_key for group_key in grouped
                if len(group_stats[(top_n, group_key)]["species"]) >= minimum
            }

    frozen_keys = {key(row) for row in frozen}
    reconstructed = selections[(100, 3)]
    if frozen_keys != reconstructed:
        raise RuntimeError(f"baseline reconstruction mismatch: frozen={len(frozen_keys)} reconstructed={len(reconstructed)}")

    sensitivity_rows = []
    for top_n in TOP_NS:
        for minimum in MIN_SPECIES:
            selected = selections[(top_n, minimum)]
            for scope in ("all", "leaf", "root", "vascular"):
                scoped = selected if scope == "all" else {k for k in selected if k[0] == scope}
                cross = [k for k in scoped if group_stats[(top_n, k)]["eudicots"] and group_stats[(top_n, k)]["monocots"]]
                balanced = [k for k in scoped if len(group_stats[(top_n, k)]["eudicots"]) >= 2 and len(group_stats[(top_n, k)]["monocots"]) >= 2]
                sensitivity_rows.append({
                    "top_n": top_n, "minimum_species_recurrence": minimum, "scope": scope,
                    "candidate_count": len(scoped), "cross_clade_presence_count": len(cross),
                    "strict_2plus2_clade_balance_count": len(balanced),
                    "frozen_102_retained_count": len(scoped & frozen_keys),
                    "frozen_102_retained_fraction": len(scoped & frozen_keys) / len(frozen_keys) if frozen_keys else "",
                    "interpretation": "descriptive sensitivity; not a replacement discovery threshold",
                })
    write_tsv(REPORT / "discovery_threshold_sensitivity.tsv", sensitivity_rows)

    stability_rows = []
    union = selections[(100, 2)]
    for group_key in sorted(union):
        s100 = group_stats[(100, group_key)]
        selected_conditions = sum(group_key in selections[c] for c in selections)
        stability_rows.append({
            "panel": group_key[0], "coarse_label": group_key[1], "orthogroup_id": group_key[2],
            "species_top25": "; ".join(group_stats[(25, group_key)]["species"]),
            "species_count_top25": len(group_stats[(25, group_key)]["species"]),
            "species_top50": "; ".join(group_stats[(50, group_key)]["species"]),
            "species_count_top50": len(group_stats[(50, group_key)]["species"]),
            "species_top100": "; ".join(s100["species"]), "species_count_top100": len(s100["species"]),
            "eudicot_count_top100": len(s100["eudicots"]), "monocot_count_top100": len(s100["monocots"]),
            "cross_clade_presence_top100": str(bool(s100["eudicots"] and s100["monocots"])).lower(),
            "strict_2plus2_clade_balance_top100": str(len(s100["eudicots"]) >= 2 and len(s100["monocots"]) >= 2).lower(),
            "selected_in_frozen_top100_min3": str(group_key in frozen_keys).lower(),
            "selected_grid_conditions": selected_conditions, "grid_conditions_total": len(selections),
            "selection_fraction": selected_conditions / len(selections),
        })
    write_tsv(REPORT / "discovery_candidate_stability.tsv", stability_rows)

    class_rows = []
    for (panel, label), species_set in sorted(class_species.items()):
        eudicots, monocots = sorted(species_set & EUDICOTS), sorted(species_set & MONOCOTS)
        class_rows.append({
            "panel": panel, "coarse_label": label, "eligible_species_count": len(species_set),
            "eligible_species": "; ".join(sorted(species_set)), "eudicot_count": len(eudicots),
            "monocot_count": len(monocots), "cross_clade_evaluable": str(bool(eudicots and monocots)).lower(),
            "strict_2plus2_evaluable": str(len(eudicots) >= 2 and len(monocots) >= 2).lower(),
            "design_note": "absence of clade eligibility is not evidence of non-conservation",
        })
    write_tsv(REPORT / "class_clade_eligibility.tsv", class_rows)

    xylem_rows = []
    for row in read_tsv(INPUTS["xylem_program"]):
        group_key = ("vascular", "xylem_lineage", row["orthogroup_id"])
        s100 = group_stats.get((100, group_key), {"species": [], "eudicots": [], "monocots": []})
        flags = [x for x in row["gap_flags"].split(";") if x]
        if len(s100["species"]) < 3:
            flags.append("not_recurrent_in_frozen_v2_top100_min3_xylem_discovery")
        if not (s100["eudicots"] and s100["monocots"]):
            flags.append("no_cross_clade_recurrence_in_frozen_v2_xylem_discovery")
        xylem_rows.append({
            "frozen_rank": row["frozen_rank"], "orthogroup_id": row["orthogroup_id"],
            "mapped_in_all_three_validation_species": row["mapped_in_all_three_validation_species"],
            "top25_species_count": len(group_stats.get((25, group_key), {"species": []})["species"]),
            "top50_species_count": len(group_stats.get((50, group_key), {"species": []})["species"]),
            "top100_species_count": len(s100["species"]), "top100_species": "; ".join(s100["species"]),
            "top100_eudicot_count": len(s100["eudicots"]), "top100_monocot_count": len(s100["monocots"]),
            "selected_top100_min3": str(group_key in frozen_keys).lower(),
            "cross_clade_presence_top100": str(bool(s100["eudicots"] and s100["monocots"])).lower(),
            "single_gene_representatives": row["single_gene_representatives"],
            "annotation_class": row["annotation_class"], "function_symbols": row["function_symbols"],
            "evidence_scope": row["annotation_scope"], "gap_flags": ";".join(sorted(set(flags))),
            "allowed_claim": "program membership/mapping and exploratory recurrence only; not one-to-one orthology or causality",
        })
    if len(xylem_rows) != 18:
        raise RuntimeError("expected exactly 18 xylem program rows")
    write_tsv(REPORT / "xylem_program_discovery_ortholog_gap.tsv", xylem_rows)

    study = {r["orthogroup_id"]: r for r in read_tsv(INPUTS["study"])}
    orthogroup = {r["orthogroup_id"]: r for r in read_tsv(INPUTS["orthogroup"])}
    tree = {r["orthogroup_id"]: r for r in read_tsv(INPUTS["tree"])}
    joint = {r["orthogroup_id"]: r for r in read_tsv(INPUTS["joint"])}
    dominance: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(INPUTS["dominance"]):
        dominance[row["orthogroup_id"]].append(row)
    priority_rows = []
    for og in sorted(study):
        s, o, t = study[og], orthogroup[og], tree[og]
        d = dominance[og]
        j = joint.get(og, {})
        robust = s["cross_study_family_status"] == "cross_species_study_robust"
        supported = t["supported_tree_pattern"] == "well_supported_dominant_subclade"
        if robust and supported:
            tier = "joint_expression_and_supported_subclade_candidate"
        elif robust:
            tier = "expression_robust_tree_not_confirmatory"
        elif supported:
            tier = "supported_subclade_expression_not_robust"
        else:
            tier = "preliminary_or_insufficient_joint_evidence"
        priority_rows.append({
            "panel": s["panel"], "coarse_label": s["coarse_label"], "orthogroup_id": og,
            "discovery_species_count_top100": len(group_stats.get((100, key(s)), {"species": []})["species"]),
            "cross_study_family_status": s["cross_study_family_status"],
            "heldout_tests": s["heldout_tests"], "heldout_positive_fraction": s["heldout_positive_fraction"],
            "species_with_positive_gene_level_marker": o["species_with_positive_gene_level_marker"],
            "species_with_single_dominant_pattern": sum(r["within_species_pattern"] == "single_dominant" for r in d),
            "species_with_dominance_rows": len(d), "dominant_genes_by_species": o["dominant_genes_by_species"],
            "supported_tree_pattern": t["supported_tree_pattern"],
            "observed_mrca_sh_alrt": t["observed_mrca_sh_alrt"],
            "observed_mrca_ultrafast_bootstrap": t["observed_mrca_ultrafast_bootstrap"],
            "joint_evidence_status_existing": j.get("joint_evidence_status", ""),
            "conservative_evidence_tier": tier,
            "one_to_one_ortholog_claim": "not_established", "paralog_substitution_claim": "not_established",
            "main_gap": "no independent functional validation and no family-wide one-to-one orthology proof",
        })
    if len(priority_rows) != 19:
        raise RuntimeError("expected exactly 19 priority-family rows")
    write_tsv(REPORT / "priority_family_single_gene_ortholog_evidence_gap.tsv", priority_rows)

    strict_joint = sum(r["conservative_evidence_tier"] == "joint_expression_and_supported_subclade_candidate" for r in priority_rows)
    summary = {
        "status": "COMPLETE_WITH_GAPS", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "version": "discovery_threshold_ortholog_robustness_v1",
        "checks": {"frozen_102_exactly_reproduced": True, "xylem_rows_18": True, "priority_rows_19": True},
        "counts": {
            "frozen_candidates": len(frozen_keys), "top25_min3_candidates": len(selections[(25, 3)]),
            "top50_min3_candidates": len(selections[(50, 3)]), "top100_min3_candidates": len(selections[(100, 3)]),
            "top100_min3_cross_clade_candidates": sum(bool(group_stats[(100, k)]["eudicots"] and group_stats[(100, k)]["monocots"]) for k in frozen_keys),
            "top100_min3_strict_2plus2_candidates": sum(len(group_stats[(100, k)]["eudicots"]) >= 2 and len(group_stats[(100, k)]["monocots"]) >= 2 for k in frozen_keys),
            "submitted_xylem_ogs_recurrent_top100_min3": sum(r["selected_top100_min3"] == "true" for r in xylem_rows),
            "priority_families_joint_strict_candidates": strict_joint,
        },
        "interpretation_boundary": "threshold sensitivity and internal evidence integration only; no independent functional validation, one-to-one orthology proof, causal marker proof, or paralog-substitution proof",
        "inputs": {
            **{name: {"path": str(path), "sha256": sha256(path)} for name, path in INPUTS.items()},
            "freeze_protocol": {"path": str(STAGED_FREEZE), "sha256": sha256(STAGED_FREEZE)},
        },
        "safety": {"source_files_modified": False, "model_training_started": False, "new_data_downloaded": False, "existing_results_overwritten": False},
    }
    (REPORT / "audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = [
        "# Discovery-threshold and ortholog-evidence robustness v1", "",
        f"The frozen top-100/minimum-3 set was exactly reconstructed ({len(frozen_keys)} candidates).",
        f"Candidate counts at minimum three species were {len(selections[(25, 3)])}, {len(selections[(50, 3)])}, and {len(selections[(100, 3)])} for top 25, 50, and 100, respectively.",
        f"Among the frozen 102, {summary['counts']['top100_min3_cross_clade_candidates']} had both monocot and eudicot recurrence and {summary['counts']['top100_min3_strict_2plus2_candidates']} met the strict two-plus-two balance description.",
        f"Only {summary['counts']['submitted_xylem_ogs_recurrent_top100_min3']}/18 submitted xylem-program orthogroups were independently rediscovered by the frozen vascular-xylem top-100/minimum-3 rule.",
        f"Strict joint expression-plus-supported-subclade candidates among the 19 priority families: {strict_joint}.", "",
        "These results are sensitivity and gap evidence. They do not establish one-to-one orthology, conserved function, causality, or paralog substitution.",
    ]
    (REPORT / "report_en.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (REPORT / "STATUS").write_text("COMPLETE_WITH_GAPS\n", encoding="utf-8")
    targets = [p for p in sorted(REPORT.iterdir()) if p.is_file() and p.name != "output_sha256.tsv"]
    write_tsv(REPORT / "output_sha256.tsv", [{"output_file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in targets])
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
