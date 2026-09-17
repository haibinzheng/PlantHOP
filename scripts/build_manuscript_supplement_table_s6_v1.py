#!/usr/bin/env python3
"""Build submission-ready Supplementary Table S6 from the completed P0 package."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(".")
SOURCE = PROJECT / "reports/label_resolution_circularity_robustness_v1"
OUTPUT = PROJECT / "reports/manuscript_supplement_table_s6_v1"
FILES = {
    "baseline": SOURCE / "baseline_nested_label_metrics.tsv",
    "direction": SOURCE / "parameter_direction_concordance.tsv",
    "circularity": SOURCE / "circularity_partial_exclusion_comparison.tsv",
    "coverage": SOURCE / "circularity_coverage.tsv",
    "p0_audit": SOURCE / "audit.json",
    "p0_hash_manifest": SOURCE / "output_sha256.tsv",
}


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


FIELDS = [
    "section", "species", "endpoint_id", "scope", "program", "k", "permutation_id", "order_condition",
    "positive_cells", "negative_cells", "auroc", "average_precision", "smd", "mean_positive", "mean_negative",
    "conditions", "auroc_min", "auroc_max", "smd_min", "smd_max", "all_auroc_above_chance", "all_smd_positive",
    "joint_positive_direction_fraction", "full_auroc", "partial_exclusion_auroc", "partial_minus_full_auroc",
    "full_average_precision", "partial_exclusion_average_precision", "partial_minus_full_average_precision",
    "full_smd", "partial_exclusion_smd", "partial_minus_full_smd", "coverage", "allowed_claim", "residual_risk",
    "label_evidence", "scoring_layer", "interpretation",
]


def blank_row(section: str, row: dict[str, str]) -> dict[str, object]:
    out = {field: "" for field in FIELDS}
    out.update({"section": section, "species": row["species"], "endpoint_id": row["endpoint_id"], "scope": row["scope"]})
    return out


def main() -> int:
    if OUTPUT.exists():
        raise RuntimeError(f"refusing to overwrite {OUTPUT}")
    for path in FILES.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    source_audit = json.loads(FILES["p0_audit"].read_text(encoding="utf-8"))
    if source_audit.get("status") != "COMPLETE" or not all(source_audit.get("checks", {}).values()):
        raise RuntimeError("Nested-label analysis record is not fully complete")
    OUTPUT.mkdir(parents=False)
    (OUTPUT / Path(__file__).name).write_bytes(Path(__file__).read_bytes())

    table_rows: list[dict[str, object]] = []
    baseline = [r for r in read_tsv(FILES["baseline"]) if r["program"] == "complete_18og_family"]
    for row in baseline:
        out = blank_row("k128_original_order_metric", row)
        for field in ("program", "k", "permutation_id", "order_condition", "positive_cells", "negative_cells", "auroc", "average_precision", "smd", "mean_positive", "mean_negative", "label_evidence", "scoring_layer"):
            out[field] = row[field]
        out["interpretation"] = "Frozen K=128 immutable-source-order cell-level discrimination metric."
        table_rows.append(out)

    direction = [r for r in read_tsv(FILES["direction"]) if r["program"] == "complete_18og_family"]
    for row in direction:
        out = blank_row("k_and_tie_order_direction_robustness", row)
        for field in ("program", "conditions", "auroc_min", "auroc_max", "smd_min", "smd_max", "all_auroc_above_chance", "all_smd_positive", "joint_positive_direction_fraction", "interpretation"):
            out[field] = row[field]
        out["scoring_layer"] = "RNA_counts"
        table_rows.append(out)

    circularity = read_tsv(FILES["circularity"])
    for row in circularity:
        out = blank_row("partial_label_circularity_stress_test", row)
        for field in ("full_auroc", "partial_exclusion_auroc", "partial_minus_full_auroc", "full_average_precision", "partial_exclusion_average_precision", "partial_minus_full_average_precision", "full_smd", "partial_exclusion_smd", "partial_minus_full_smd", "coverage", "allowed_claim", "residual_risk"):
            out[field] = row[field]
        out["k"] = "128"
        out["permutation_id"] = "0"
        out["order_condition"] = "immutable_source_order"
        out["scoring_layer"] = "RNA_counts"
        out["interpretation"] = "Removal of the documented 5/18 published-marker-overlap orthogroups; incomplete circularity coverage."
        table_rows.append(out)

    section_counts = {section: sum(r["section"] == section for r in table_rows) for section in sorted({r["section"] for r in table_rows})}
    if len(table_rows) != 54 or set(section_counts.values()) != {18}:
        raise RuntimeError(f"unexpected S6 rows: total={len(table_rows)} sections={section_counts}")
    write_tsv(OUTPUT / "Table_S6_label_resolution_circularity_parameter_robustness.tsv", table_rows, FIELDS)

    descriptions = {
        "section": "Record type within Table S6.", "species": "Species code used by the frozen GSE268881 analysis (esa or sir).",
        "endpoint_id": "Frozen nested binary endpoint.", "scope": "Pooled cells or within-study replicate R1/R2.",
        "program": "Frozen orthogroup-program variant.", "k": "Frozen top-K rank-encoding value when applicable.",
        "permutation_id": "Feature tie-order condition identifier; 0 is immutable source order.", "order_condition": "Feature-order condition label.",
        "positive_cells": "Number of endpoint-positive cells.", "negative_cells": "Number of endpoint-negative cells.",
        "auroc": "Area under the receiver-operating-characteristic curve.", "average_precision": "Average precision under endpoint class prevalence.",
        "smd": "Standardized mean difference of program scores, positive minus negative.", "mean_positive": "Mean program score in endpoint-positive cells.",
        "mean_negative": "Mean program score in endpoint-negative cells.", "conditions": "Number of K and feature-order conditions summarized.",
        "auroc_min": "Minimum AUROC across summarized conditions.", "auroc_max": "Maximum AUROC across summarized conditions.",
        "smd_min": "Minimum SMD across summarized conditions.", "smd_max": "Maximum SMD across summarized conditions.",
        "all_auroc_above_chance": "Whether AUROC exceeded 0.5 in every summarized condition.",
        "all_smd_positive": "Whether SMD was positive in every summarized condition.",
        "joint_positive_direction_fraction": "Fraction of summarized conditions with AUROC > 0.5 and SMD > 0.",
        "full_auroc": "AUROC for the complete frozen 18-OG program.", "partial_exclusion_auroc": "AUROC after documented 5-OG overlap exclusion.",
        "partial_minus_full_auroc": "Partial-exclusion minus complete-program AUROC.", "full_average_precision": "Average precision for the complete program.",
        "partial_exclusion_average_precision": "Average precision after documented overlap exclusion.",
        "partial_minus_full_average_precision": "Partial-exclusion minus complete-program average precision.",
        "full_smd": "SMD for the complete program.", "partial_exclusion_smd": "SMD after documented overlap exclusion.",
        "partial_minus_full_smd": "Partial-exclusion minus complete-program SMD.", "coverage": "Coverage of the circularity stress test.",
        "allowed_claim": "Maximum supported circularity claim.", "residual_risk": "Unresolved label-circularity limitation.",
        "label_evidence": "Provenance of endpoint labels.", "scoring_layer": "Expression layer used for scoring.",
        "interpretation": "Section-specific interpretation boundary.",
    }
    write_tsv(OUTPUT / "field_dictionary.tsv", [{"field": f, "description": descriptions[f], "missing_value_policy": "Blank means not applicable or unavailable; never interpret as zero."} for f in FIELDS])
    write_tsv(OUTPUT / "input_source_manifest.tsv", [{"source_id": name, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size} for name, path in FILES.items()])

    coverage = read_tsv(FILES["coverage"])[0]
    readme = [
        "# Supplementary Table S6", "",
        "Submission-ready long-format table containing 54 records in three 18-row sections:",
        "1. K=128 immutable-source-order metrics for pooled, R1, and R2 scopes.",
        "2. Direction robustness across K=64/128/256 and 11 feature-order conditions.",
        "3. Complete-program versus documented 5/18 marker-overlap exclusion comparisons.", "",
        f"Circularity coverage is {coverage['excluded_frozen_orthogroups']}/{coverage['frozen_program_orthogroups']} frozen orthogroups. The complete label-transfer feature set is unavailable, so circularity is not completely excluded.",
        "All scores use RNA counts; integrated scale.data was excluded. Results are cell-level discrimination and within-study repeatability, not broad species-level generalization.",
    ]
    (OUTPUT / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    audit = {
        "status": "COMPLETE_WITH_LIMITATIONS", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "version": "manuscript_supplement_table_s6_v1", "table_rows": len(table_rows), "section_counts": section_counts,
        "checks": {"p0_authority_complete": True, "expected_total_rows_54": True, "each_section_rows_18": True, "coverage_5_of_18": coverage["excluded_frozen_orthogroups"] == "5" and coverage["frozen_program_orthogroups"] == "18"},
        "interpretation_boundary": "partial documented-marker overlap stress test only; residual label-circularity risk remains",
        "source_p0_audit_sha256": sha256(FILES["p0_audit"]),
        "safety": {"p0_directory_modified": False, "source_data_modified": False, "existing_results_overwritten": False},
    }
    (OUTPUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUTPUT / "STATUS").write_text("COMPLETE_WITH_LIMITATIONS\n", encoding="utf-8")
    targets = [p for p in sorted(OUTPUT.iterdir()) if p.is_file() and p.name != "output_sha256.tsv"]
    write_tsv(OUTPUT / "output_sha256.tsv", [{"output_file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in targets])
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
