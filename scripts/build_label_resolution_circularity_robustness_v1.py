#!/usr/bin/env python3
"""Summarize the frozen GSE268881 label-resolution robustness run."""

from __future__ import annotations

import csv
import hashlib
import json
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(".")
RUN = Path("data/runs/PhyloOpenCell/label_resolution_circularity_robustness_v1")
REPORT = PROJECT / "reports/label_resolution_circularity_robustness_v1"


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
        raise ValueError(f"empty output: {path}")
    names = fields or list(rows[0])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def main() -> int:
    if not REPORT.is_dir() or not RUN.is_dir():
        raise RuntimeError("freeze report directory and completed run are required")
    required_freeze = [REPORT / "freeze_protocol_v1.md", REPORT / "label_definition_freeze.tsv", REPORT / "input_manifest.tsv"]
    required_run = [RUN / "nested_label_parameter_metrics.tsv", RUN / "expression_layer_audit.tsv", RUN / "circularity_coverage.tsv", RUN / "audit.json"]
    for path in required_freeze + required_run:
        if not path.is_file():
            raise FileNotFoundError(path)
    run_audit = json.loads((RUN / "audit.json").read_text(encoding="utf-8"))
    if run_audit["status"] != "COMPLETE" or run_audit["metric_rows"] != run_audit["expected_metric_rows"]:
        raise RuntimeError("run audit is not complete")

    metrics = read_tsv(RUN / "nested_label_parameter_metrics.tsv")
    baseline = [r for r in metrics if r["k"] == "128" and r["permutation_id"] == "0"]
    write_tsv(REPORT / "baseline_nested_label_metrics.tsv", baseline)

    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        key = (row["species"], row["endpoint_id"], row["program"], row["scope"], row["k"])
        groups[key].append(row)
    parameter_rows = []
    for key in sorted(groups):
        rows = groups[key]
        aucs = [f(r, "auroc") for r in rows]
        aps = [f(r, "average_precision") for r in rows]
        smds = [f(r, "smd") for r in rows]
        parameter_rows.append({
            "species": key[0], "endpoint_id": key[1], "program": key[2], "scope": key[3], "k": key[4],
            "order_conditions": len(rows), "auroc_min": min(aucs), "auroc_median": statistics.median(aucs), "auroc_max": max(aucs),
            "average_precision_min": min(aps), "average_precision_median": statistics.median(aps), "average_precision_max": max(aps),
            "smd_min": min(smds), "smd_median": statistics.median(smds), "smd_max": max(smds),
            "positive_auroc_direction_fraction": sum(x > 0.5 for x in aucs) / len(aucs),
            "positive_smd_direction_fraction": sum(x > 0 for x in smds) / len(smds),
        })
    write_tsv(REPORT / "parameter_robustness_summary.tsv", parameter_rows)

    overall_groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in metrics:
        overall_groups[(row["species"], row["endpoint_id"], row["program"], row["scope"])].append(row)
    direction_rows = []
    for key in sorted(overall_groups):
        rows = overall_groups[key]
        aucs = [f(r, "auroc") for r in rows]
        smds = [f(r, "smd") for r in rows]
        direction_rows.append({
            "species": key[0], "endpoint_id": key[1], "program": key[2], "scope": key[3],
            "conditions": len(rows), "auroc_min": min(aucs), "auroc_max": max(aucs), "smd_min": min(smds), "smd_max": max(smds),
            "all_auroc_above_chance": str(all(x > 0.5 for x in aucs)).lower(),
            "all_smd_positive": str(all(x > 0 for x in smds)).lower(),
            "joint_positive_direction_fraction": sum(a > 0.5 and s > 0 for a, s in zip(aucs, smds)) / len(rows),
            "interpretation": "descriptive direction across K and tie-order conditions; not a new confirmatory endpoint",
        })
    write_tsv(REPORT / "parameter_direction_concordance.tsv", direction_rows)

    baseline_map = {(r["species"], r["endpoint_id"], r["scope"], r["program"]): r for r in baseline}
    circularity_rows = []
    for species, endpoint, scope, _ in sorted(k for k in baseline_map if k[3] == "complete_18og_family"):
        full = baseline_map[(species, endpoint, scope, "complete_18og_family")]
        partial = baseline_map[(species, endpoint, scope, "partial_published_marker_overlap_exclusion")]
        circularity_rows.append({
            "species": species, "endpoint_id": endpoint, "scope": scope,
            "full_auroc": full["auroc"], "partial_exclusion_auroc": partial["auroc"], "partial_minus_full_auroc": f(partial, "auroc") - f(full, "auroc"),
            "full_average_precision": full["average_precision"], "partial_exclusion_average_precision": partial["average_precision"],
            "partial_minus_full_average_precision": f(partial, "average_precision") - f(full, "average_precision"),
            "full_smd": full["smd"], "partial_exclusion_smd": partial["smd"], "partial_minus_full_smd": f(partial, "smd") - f(full, "smd"),
            "coverage": "5/18 frozen OGs; complete label-transfer feature set unavailable",
            "allowed_claim": "robustness to removal of documented published-marker overlap only",
            "residual_risk": "label circularity is not completely excluded",
        })
    write_tsv(REPORT / "circularity_partial_exclusion_comparison.tsv", circularity_rows)

    for source_name, target_name in (
        ("nested_label_parameter_metrics.tsv", "nested_label_parameter_metrics.tsv"),
        ("expression_layer_audit.tsv", "expression_layer_audit.tsv"),
        ("circularity_coverage.tsv", "circularity_coverage.tsv"),
        ("audit.json", "run_audit.json"),
    ):
        target = REPORT / target_name
        if target.exists():
            raise RuntimeError(f"refusing to overwrite {target}")
        target.write_bytes((RUN / source_name).read_bytes())

    fields = []
    for table_name, rows in {
        "baseline_nested_label_metrics": baseline,
        "parameter_robustness_summary": parameter_rows,
        "parameter_direction_concordance": direction_rows,
        "circularity_partial_exclusion_comparison": circularity_rows,
    }.items():
        for field in rows[0]:
            fields.append({
                "table": table_name, "field": field,
                "description": field.replace("_", " ").capitalize() + ".",
                "missing_value_policy": "Blank means unavailable or not applicable; never interpret as zero.",
            })
    write_tsv(REPORT / "field_dictionary.tsv", fields, ["table", "field", "description", "missing_value_policy"])

    main_rows = [r for r in baseline if r["endpoint_id"] == "xylem_vs_all_other_root" and r["program"] == "complete_18og_family" and r["scope"] == "pooled"]
    nested_rows = [r for r in baseline if r["program"] == "complete_18og_family" and r["scope"] == "pooled"]
    circ_rows = [r for r in circularity_rows if r["endpoint_id"] == "xylem_vs_all_other_root" and r["scope"] == "pooled"]
    all_full_direction = [r for r in direction_rows if r["program"] == "complete_18og_family"]
    report_en = [
        "# Label resolution, circularity and parameter robustness v1", "",
        "All metrics are cell-level discrimination of author-assigned/label-transfer-derived labels. R1/R2 are within-study biological replicates; two species do not establish broad species-level generalization.", "",
        "## Frozen K=128/source-order results", "",
    ]
    for row in sorted(nested_rows, key=lambda x: (x["species"], x["endpoint_id"])):
        report_en.append(f"- {row['species']} / {row['endpoint_id']}: AUROC {f(row,'auroc'):.3f}, AP {f(row,'average_precision'):.3f}, SMD {f(row,'smd'):.3f}; n+={row['positive_cells']}, n-={row['negative_cells']}.")
    report_en += ["", "## Circularity boundary", ""]
    for row in sorted(circ_rows, key=lambda x: x["species"]):
        report_en.append(f"- {row['species']}: removing the documented 5/18-OG published-marker overlap changed AUROC by {float(row['partial_minus_full_auroc']):+.3f}.")
    report_en += [
        "",
        "This is partial coverage only. The complete feature set used for label assignment/transfer is unavailable, so residual label-circularity risk remains and complete exclusion of circularity must not be claimed.",
        "",
        "## Parameter robustness", "",
        f"Across all K/order/scope combinations for the complete program, {sum(r['all_auroc_above_chance']=='true' and r['all_smd_positive']=='true' for r in all_full_direction)}/{len(all_full_direction)} species-endpoint-scope groups retained positive direction in every condition.",
        "The full condition-level summaries are in parameter_robustness_summary.tsv and parameter_direction_concordance.tsv. The integrated scaled matrix was excluded; all scores used RNA counts.",
    ]
    report_zh = [
        "# 标签分辨率、循环性与参数稳健性 v1", "",
        "所有指标均为对作者赋值/标签转移所得标签的细胞层面区分度。R1/R2仅代表研究内重复；两个物种不能支持广泛物种层面的泛化结论。", "",
        "## 冻结 K=128/原始顺序结果", "",
    ]
    for row in sorted(nested_rows, key=lambda x: (x["species"], x["endpoint_id"])):
        report_zh.append(f"- {row['species']} / {row['endpoint_id']}：AUROC {f(row,'auroc'):.3f}，AP {f(row,'average_precision'):.3f}，SMD {f(row,'smd'):.3f}；阳性={row['positive_cells']}，阴性={row['negative_cells']}。")
    report_zh += ["", "## 循环性边界", ""]
    for row in sorted(circ_rows, key=lambda x: x["species"]):
        report_zh.append(f"- {row['species']}：排除已记录的5/18个公开marker重叠OG后，AUROC变化 {float(row['partial_minus_full_auroc']):+.3f}。")
    report_zh += [
        "", "这只是部分覆盖。完整标签赋值/转移特征集不可得，因此仍存在残余标签循环性风险，不得声称已经完全排除循环性。", "",
        "## 参数稳健性", "",
        f"完整程序的全部K/顺序/范围组合中，{sum(r['all_auroc_above_chance']=='true' and r['all_smd_positive']=='true' for r in all_full_direction)}/{len(all_full_direction)}个“物种-端点-范围”组在所有条件下均保持正方向。",
        "完整明细见 parameter_robustness_summary.tsv 与 parameter_direction_concordance.tsv。所有评分均使用RNA counts，未使用integrated scale.data。",
    ]
    (REPORT / "report_en.md").write_text("\n".join(report_en) + "\n", encoding="utf-8")
    (REPORT / "report_zh.md").write_text("\n".join(report_zh) + "\n", encoding="utf-8")

    summary_audit = {
        "status": "COMPLETE", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "version": "label_resolution_circularity_robustness_v1", "run_audit_sha256": sha256(RUN / "audit.json"),
        "checks": {
            "run_complete": True, "metric_rows_expected": len(metrics) == 1188,
            "baseline_rows_expected": len(baseline) == 36, "parameter_groups_expected": len(parameter_rows) == 108,
            "direction_groups_expected": len(direction_rows) == 36, "circularity_comparison_rows_expected": len(circularity_rows) == 18,
            "authority_reproduction_passed": all(v["pooled_exact_1e12"] and v["replicate_auroc_exact_1e12"] for v in run_audit["authority_checks"].values()),
            "scaled_matrix_excluded": not run_audit["safety"]["scaled_matrix_used"],
        },
        "safety": run_audit["safety"],
        "interpretation_boundary": "cell-level discrimination and within-study repeatability; not broad replicate/species-level generalization; partial circularity coverage only",
    }
    (REPORT / "audit.json").write_text(json.dumps(summary_audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (REPORT / "STATUS").write_text("COMPLETE\n", encoding="utf-8")

    hash_targets = [p for p in sorted(REPORT.iterdir()) if p.is_file() and p.name != "output_sha256.tsv"]
    write_tsv(REPORT / "output_sha256.tsv", [{"output_file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in hash_targets], ["output_file", "sha256", "bytes"])
    print(json.dumps(summary_audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
