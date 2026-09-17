#!/usr/bin/env python3
"""Render submission-scale PlantHOP Figures 1-5 and Supplementary Figures S1-S2."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
import textwrap
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch
import numpy as np
from PIL import Image


REPOSITORY = Path(__file__).resolve().parents[1]
EVIDENCE = REPOSITORY / "results" / "figure_data"
OUT = REPOSITORY / "results" / "figures"
DESIGN = EVIDENCE / "figure_design_freeze_v1_1.json"
SCRIPT = Path(__file__).resolve()

SOURCE = {
    "figure1": EVIDENCE / "figure1_resource_flow.tsv",
    "figure2": EVIDENCE / "figure2_independent_validation.tsv",
    "figure3": EVIDENCE / "figure3_prospective_xylem.tsv",
    "figure4": EVIDENCE / "figure4_controls_comparators.tsv",
    "figure5": EVIDENCE / "figure5_xylem_program.tsv",
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


def fl(value: str) -> float:
    return float(value)


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=FS["panel_label"], fontweight="bold", va="top")


def clean_axes(ax, grid_axis: str | None = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=C["grid"], linewidth=0.7, zorder=0)
    ax.tick_params(labelsize=FS["tick"], length=3)


def add_chance(ax) -> None:
    ax.axhline(0.5, color="#777777", linewidth=0.9, linestyle="--", zorder=1)


def save(fig, stem: str, source_keys: list[str] | None = None, layout_key: str | None = None, manifest_key: str | None = None) -> dict[str, object]:
    fig.align_labels()
    png = OUT / f"{stem}.png"
    pdf = OUT / f"{stem}.pdf"
    fig.savefig(png, dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white", metadata={"Title": stem, "Creator": "PlantHOP manuscript_figures_v1_1"})
    plt.close(fig)
    with Image.open(png) as image:
        pixels = image.size
    return {
        "stem": stem,
        "png": png,
        "pdf": pdf,
        "png_width_px": pixels[0],
        "png_height_px": pixels[1],
        "source_keys": source_keys or [stem.split("_")[0]],
        "layout_key": layout_key or stem.split("_")[0],
        "manifest_key": manifest_key or stem.split("_")[0],
    }


def figure1() -> dict[str, object]:
    data = read_tsv(SOURCE["figure1"])
    summary = {r["metric"]: r for r in data if r["record_type"] == "collection_summary"}
    decisions = [r for r in data if r["record_type"] == "panel_decision"]
    workflow = sorted([r for r in data if r["record_type"] == "workflow"], key=lambda r: int(r["workflow_order"]))
    total = int(summary["cells"]["value"])
    labeled = int(summary["labeled_cells"]["value"])
    unknown = int(summary["unknown_or_missing_cells"]["value"])

    fig = plt.figure(figsize=(WIDTH, LAYOUT["figure1"]["height_inches"]), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.08])

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.barh([0], [labeled], color=C["labeled"], height=0.42, label="Author-labelled")
    ax.barh([0], [unknown], left=[labeled], color=C["unknown"], height=0.42, label="Unknown / missing")
    ax.set_xlim(0, total)
    ax.set_yticks([])
    ax.set_xlabel("Audited cells", fontsize=FS["axis_label"])
    ax.set_title("Resource audit", fontsize=FS["title"], loc="left")
    ax.text(labeled / 2, 0, f"{labeled:,}\n68.9% labelled", ha="center", va="center", fontsize=FS["annotation"], color="white", fontweight="bold")
    ax.text(labeled + unknown / 2, 0, f"{unknown:,}\n31.1% unknown", ha="center", va="center", fontsize=FS["annotation"], color=C["text"])
    ax.text(0.0, 0.93, f"121 readable H5AD  |  {total:,} cells  |  33 species", transform=ax.transAxes, fontsize=FS["axis_label"], va="top")
    clean_axes(ax, "x")

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Hierarchical label resolution", fontsize=FS["title"], loc="left")
    ax.axis("off")
    ax.text(0.20, 0.84, "Coarse system", transform=ax.transAxes, ha="center", fontsize=FS["axis_label"], fontweight="bold", color=C["resource"])
    ax.text(0.77, 0.84, "Resolved lineage", transform=ax.transAxes, ha="center", fontsize=FS["axis_label"], fontweight="bold", color=C["complete_family"])
    for y, coarse in [(0.62, "Vascular"), (0.34, "Stele")]:
        box = FancyBboxPatch((0.04, y - 0.075), 0.32, 0.15, boxstyle="round,pad=0.02,rounding_size=0.02", transform=ax.transAxes, facecolor=C["resource"], edgecolor="none", alpha=0.95)
        ax.add_patch(box)
        ax.text(0.20, y, coarse, transform=ax.transAxes, ha="center", va="center", color="white", fontsize=FS["axis_label"], fontweight="bold")
    for y, lineage in [(0.62, "Xylem"), (0.34, "Phloem")]:
        box = FancyBboxPatch((0.61, y - 0.075), 0.32, 0.15, boxstyle="round,pad=0.02,rounding_size=0.02", transform=ax.transAxes, facecolor=C["complete_family"], edgecolor="none", alpha=0.95)
        ax.add_patch(box)
        ax.text(0.77, y, lineage, transform=ax.transAxes, ha="center", va="center", color="white", fontsize=FS["axis_label"], fontweight="bold")
    for y0 in (0.62, 0.34):
        ax.plot([0.37, 0.48], [y0, y0], transform=ax.transAxes, color="#777777", lw=0.8, alpha=0.8)
    ax.plot([0.48, 0.48], [0.34, 0.62], transform=ax.transAxes, color="#777777", lw=0.8, alpha=0.8)
    for y1 in (0.62, 0.34):
        ax.annotate("", xy=(0.60, y1), xytext=(0.48, y1), xycoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>", color="#777777", lw=0.8, alpha=0.8))
    ax.text(0.04, 0.08, "Lineage resolution is label-aware; endpoint status remains prespecified.\n† Wheat xylem is sensitivity-only, not a primary pass.", transform=ax.transAxes, fontsize=FS["annotation"], color=C["text"], va="bottom")

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "C")
    panels = ["root", "leaf", "vascular"]
    included = [sum(r["accepted"].lower() == "true" and r["panel"] == p for r in decisions) for p in panels]
    excluded = [sum(r["accepted"].lower() != "true" and r["panel"] == p for r in decisions) for p in panels]
    x = np.arange(len(panels))
    ax.bar(x, included, color=C["primary_pass"], width=0.62, label="Candidate include", zorder=2)
    ax.bar(x, excluded, bottom=included, color=C["unknown"], width=0.62, label="Excluded", zorder=2)
    for i, (inc, exc) in enumerate(zip(included, excluded)):
        ax.text(i, inc / 2, str(inc), ha="center", va="center", fontsize=FS["annotation"], color="white", fontweight="bold")
        ax.text(i, inc + exc / 2, str(exc), ha="center", va="center", fontsize=FS["annotation"])
    ax.set_xticks(x, [p.capitalize() for p in panels])
    ax.set_ylabel("Dataset-panel decisions", fontsize=FS["axis_label"])
    ax.set_title("Frozen inclusion is panel-specific", fontsize=FS["title"], loc="left")
    ax.legend(frameon=False, fontsize=FS["legend"], ncol=2, loc="upper center")
    clean_axes(ax)

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "D")
    ax.set_title("Discovery → freeze → independent validation", fontsize=FS["title"], loc="left")
    ax.axis("off")
    colors = [C["resource"], C["single_gene"], C["primary_pass"], C["complete_family"]]
    short = [
        "Discovery\nresource audit",
        "Pre-score freeze\ndatasets · labels · programs",
        "Independent validation\nrice · Moricandia · wheat",
        "Prospective xylem\nEsa · Sir; Spa excluded",
    ]
    y = np.linspace(0.86, 0.14, len(workflow))
    for i, (yy, label, color) in enumerate(zip(y, short, colors)):
        box = FancyBboxPatch((0.08, yy - 0.075), 0.78, 0.13, boxstyle="round,pad=0.02,rounding_size=0.02", transform=ax.transAxes, facecolor=color, edgecolor="none", alpha=0.95)
        ax.add_patch(box)
        ax.text(0.47, yy - 0.01, label, transform=ax.transAxes, ha="center", va="center", color="white" if i != 1 else C["text"], fontsize=FS["axis_label"], fontweight="bold")
        if i < len(y) - 1:
            ax.annotate("", xy=(0.47, y[i + 1] + 0.075), xytext=(0.47, yy - 0.085), xycoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.2))
    return save(fig, "figure1_resource_to_validation", source_keys=["figure1", "figure2"])


def figure2() -> dict[str, object]:
    data = read_tsv(SOURCE["figure2"])
    pooled = [r for r in data if r["scope"] == "pooled"]
    replicates = [r for r in data if r["scope"] == "replicate"]
    order = ["moricandia_leaf_photosynthetic_v1", "moricandia_leaf_epidermal_v1", "wheat_root_stele_v1", "wheat_root_xylem_v1"]
    label = {order[0]: "Moricandia\nphotosynthetic", order[1]: "Moricandia\nepidermal", order[2]: "Wheat\nstele", order[3]: "Wheat\nxylem†"}
    status_color = {"primary_pass": C["primary_pass"], "primary_fail": C["primary_fail"], "sensitivity_support_not_primary": C["sensitivity"]}
    pidx = {r["contrast_id"]: r for r in pooled}

    fig, axs = plt.subplots(2, 2, figsize=(WIDTH, LAYOUT["figure2"]["height_inches"]), constrained_layout=False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.10, wspace=0.28, hspace=0.48)
    ax = axs[0, 0]
    panel_label(ax, "A")
    x = np.arange(4)
    vals = [fl(pidx[k]["auroc"]) for k in order]
    colors = [status_color[pidx[k]["endpoint_status"]] for k in order]
    bars = ax.bar(x, vals, color=colors, width=0.66, zorder=2)
    bars[-1].set_hatch("//")
    add_chance(ax)
    ax.set_ylim(0.45, 1.02)
    ax.set_xticks(x, [label[k] for k in order])
    ax.set_ylabel("Pooled AUROC", fontsize=FS["axis_label"])
    ax.set_title("Pooled endpoint outcomes", fontsize=FS["title"], loc="left")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.018, f"{v:.3f}", ha="center", fontsize=FS["annotation"])
    clean_axes(ax)

    ax = axs[0, 1]
    panel_label(ax, "B")
    for i, cid in enumerate(order):
        rr = [r for r in replicates if r["contrast_id"] == cid]
        ys = [fl(r["auroc"]) for r in rr]
        xx = np.linspace(i - 0.13, i + 0.13, len(ys))
        ax.plot(xx, ys, color=colors[i], linewidth=1, alpha=0.8)
        ax.scatter(xx, ys, color=colors[i], edgecolor="white", linewidth=0.6, s=34, zorder=3, marker="D" if i == 3 else "o")
    add_chance(ax)
    ax.set_ylim(0.45, 1.02)
    ax.set_xticks(x, [label[k] for k in order])
    ax.set_ylabel("Replicate AUROC", fontsize=FS["axis_label"])
    ax.set_title("Replicate consistency", fontsize=FS["title"], loc="left")
    clean_axes(ax)

    ax = axs[1, 0]
    panel_label(ax, "C")
    deltas = [fl(pidx[k]["family_minus_single_auroc"]) for k in order]
    ax.axhline(0, color="#777777", linewidth=0.9, linestyle="--")
    ax.bar(x, deltas, color=colors, width=0.66, zorder=2)
    ax.set_xticks(x, [label[k] for k in order])
    ax.set_ylabel("Family − single-gene AUROC", fontsize=FS["axis_label"])
    ax.set_title("Family − single-gene AUROC", fontsize=FS["title"], loc="left")
    ax.set_ylim(-0.02, 0.43)
    for i, v in enumerate(deltas):
        ax.text(i, v + 0.012, f"{v:+.3f}", ha="center", fontsize=FS["annotation"])
    clean_axes(ax)

    ax = axs[1, 1]
    panel_label(ax, "D")
    ax.axis("off")
    ax.set_title("Endpoint status and sample sizes", fontsize=FS["title"], loc="left")
    handles = [
        Patch(facecolor=C["primary_pass"], label="Primary pass"),
        Patch(facecolor=C["primary_fail"], label="Primary fail"),
        Patch(facecolor=C["sensitivity"], hatch="//", label="Sensitivity; not primary"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=FS["legend"], loc="upper left")
    lines = []
    for cid in order:
        r = pidx[cid]
        lines.append(f"{label[cid].replace(chr(10), ' ')}: n+={int(r['positive_cells']):,}, n−={int(r['negative_cells']):,}; AP={fl(r['average_precision']):.3f}; SMD={fl(r['smd']):.2f}; P={fl(r['empirical_p']):.4f}")
    ax.text(0.02, 0.56, "\n\n".join(textwrap.fill(t, 54) for t in lines), transform=ax.transAxes, va="top", fontsize=FS["annotation"])
    ax.text(0.02, 0.02, "† Wheat xylem was prospectively sensitivity-only because replicate cell counts were sparse.", transform=ax.transAxes, fontsize=FS["annotation"], color=C["sensitivity"])
    return save(fig, "figure2_independent_validation")


def figure3() -> dict[str, object]:
    data = read_tsv(SOURCE["figure3"])
    pooled = [r for r in data if r["scope"] == "pooled"]
    reps = [r for r in data if r["scope"] == "replicate" and r["species"] in ("esa", "sir")]
    all_rep = [r for r in data if r["scope"] == "replicate"]
    species_colors = {"esa": C["primary_pass"], "sir": C["complete_family"], "spa": C["excluded"]}

    fig = plt.figure(figsize=(WIDTH, LAYOUT["figure3"]["height_inches"]), constrained_layout=False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.10, wspace=0.32, hspace=0.48)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05])
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    for i, r in enumerate(pooled):
        y = fl(r["auroc"])
        ax.errorbar(i, y, yerr=[[y - fl(r["auroc_lower_95"])], [fl(r["auroc_upper_95"]) - y]], fmt="o", markersize=7, capsize=3, color=species_colors[r["species"]], linewidth=1.4)
        ax.text(i, y + 0.045, f"{y:.3f}", ha="center", fontsize=FS["annotation"])
    add_chance(ax)
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0.45, 1.02)
    ax.set_xticks([0, 1], ["Esa", "Sir"])
    ax.set_ylabel("Pooled AUROC", fontsize=FS["axis_label"])
    ax.set_title("Pooled replication", fontsize=FS["title"], loc="left")
    clean_axes(ax)

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    positions = {("esa", "R1"): 0, ("esa", "R2"): 1, ("sir", "R1"): 3, ("sir", "R2"): 4}
    for r in reps:
        x = positions[(r["species"], r["replicate"])]
        y = fl(r["auroc"])
        ax.errorbar(x, y, yerr=[[y - fl(r["auroc_lower_95"])], [fl(r["auroc_upper_95"]) - y]], fmt="o", markersize=5.5, capsize=2.5, color=species_colors[r["species"]], linewidth=1.1)
    add_chance(ax)
    ax.set_ylim(0.45, 1.02)
    ax.set_xticks([0, 1, 3, 4], ["Esa R1", "Esa R2", "Sir R1", "Sir R2"], rotation=20)
    ax.set_ylabel("Replicate AUROC", fontsize=FS["axis_label"])
    ax.set_title("Replicate consistency", fontsize=FS["title"], loc="left")
    clean_axes(ax)

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "C")
    order = [("esa", "R1"), ("esa", "R2"), ("sir", "R1"), ("sir", "R2"), ("spa", "R1"), ("spa", "R2")]
    lookup = {(r["species"], r["replicate"]): r for r in all_rep}
    counts = [int(lookup[k]["positive_cells"]) for k in order]
    colors = [species_colors[k[0]] for k in order]
    ax.bar(np.arange(6), counts, color=colors, width=0.68, zorder=2)
    ax.axhline(100, color=C["primary_fail"], linestyle="--", linewidth=1.1, label="Frozen minimum = 100")
    ax.set_xticks(np.arange(6), [f"{s.upper()}\n{r}" for s, r in order])
    ax.set_ylabel("Author-labelled xylem cells", fontsize=FS["axis_label"])
    ax.set_title("Pre-score eligibility", fontsize=FS["title"], loc="left")
    for i, n in enumerate(counts):
        ax.text(i, n + 4, str(n), ha="center", fontsize=FS["annotation"])
    ax.legend(frameon=False, fontsize=FS["legend"], loc="upper left")
    clean_axes(ax)

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "D")
    ax.axis("off")
    ax.set_title("Interpretation boundary", fontsize=FS["title"], loc="left")
    ax.text(0.03, 0.78, "Esa and Sir", color=C["primary_pass"], fontsize=10, fontweight="bold", transform=ax.transAxes)
    ax.text(0.03, 0.64, "Eligible and formally scored; both passed the frozen species-level rule.", fontsize=FS["axis_label"], transform=ax.transAxes, wrap=True)
    ax.text(0.03, 0.42, "Spa", color=C["excluded"], fontsize=10, fontweight="bold", transform=ax.transAxes)
    ax.text(0.03, 0.28, "Excluded before scoring because R1 contained 57 xylem cells (<100). No pooled rescue or relabelling was allowed.", fontsize=FS["axis_label"], transform=ax.transAxes, wrap=True)
    ax.text(0.03, 0.04, "Intervals: 2,000 stratified cell-bootstrap resamples; conditional on the observed biological replicates.", fontsize=FS["annotation"], fontweight="bold", transform=ax.transAxes)
    return save(fig, "figure3_prospective_xylem")


def figure4() -> dict[str, object]:
    data = read_tsv(SOURCE["figure4"])
    pooled_methods = [r for r in data if r["scope"] == "pooled" and r["record_type"] == "method"]
    deltas = [r for r in data if r["scope"] == "pooled" and r["record_type"] == "paired_difference" and r["estimand"] == "delta_auc"]
    random = [r for r in data if r["record_type"] == "random_null"]
    idx = {(r["species"], r["comparator"]): r for r in pooled_methods}
    species = ["esa", "sir"]
    sx = np.arange(2)

    fig, axs = plt.subplots(2, 2, figsize=(WIDTH, LAYOUT["figure4"]["height_inches"]), constrained_layout=False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.10, wspace=0.30, hspace=0.60)
    ax = axs[0, 0]
    panel_label(ax, "A", x=-0.17, y=1.08)
    methods = ["complete_family", "anti_circular", "phloem_off_target"]
    labels = ["Mapped members", "Marker-overlap exclusion", "Phloem off-target"]
    colors = [C["complete_family"], C["anti_circular"], C["negative_control"]]
    width = 0.22
    for j, (m, lab, col) in enumerate(zip(methods, labels, colors)):
        ax.bar(sx + (j - 1) * width, [fl(idx[(s, m)]["auroc"]) for s in species], width=width, color=col, label=lab, zorder=2)
    add_chance(ax)
    ax.set_xticks(sx, ["Esa", "Sir"])
    ax.set_ylim(0.35, 1.02)
    ax.set_ylabel("Pooled AUROC", fontsize=FS["axis_label"])
    ax.set_title("Negative controls", fontsize=FS["title"], loc="left")
    ax.legend(frameon=False, fontsize=FS["legend"], ncol=1, loc="upper center")
    clean_axes(ax)

    ax = axs[0, 1]
    panel_label(ax, "B", x=-0.17, y=1.08)
    methods = ["complete_family", "equal_weight", "single_gene"]
    labels = ["Mapped members", "Equal weight", "Single gene"]
    colors = [C["complete_family"], C["equal_weight"], C["single_gene"]]
    for j, (m, lab, col) in enumerate(zip(methods, labels, colors)):
        ax.bar(sx + (j - 1) * width, [fl(idx[(s, m)]["auroc"]) for s in species], width=width, color=col, label=lab, zorder=2)
    add_chance(ax)
    ax.set_xticks(sx, ["Esa", "Sir"])
    ax.set_ylim(0.75, 1.0)
    ax.set_ylabel("Pooled AUROC", fontsize=FS["axis_label"])
    ax.set_title("Method comparisons", fontsize=FS["title"], loc="left")
    ax.legend(frameon=False, fontsize=FS["legend"], loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, columnspacing=0.8, handletextpad=0.4)
    clean_axes(ax)

    ax = axs[1, 0]
    panel_label(ax, "C", x=-0.17, y=1.08)
    for i, r in enumerate(sorted(deltas, key=lambda x: x["species"])):
        est, lo, hi = fl(r["estimate"]), fl(r["lower_95"]), fl(r["upper_95"])
        ax.errorbar(i, est, yerr=[[est - lo], [hi - est]], fmt="o", markersize=7, capsize=4, color=C["complete_family"], linewidth=1.4)
        ax.text(i + 0.07, est, f"{est:+.3f}", ha="left", va="center", fontsize=FS["annotation"])
    ax.axhline(0, color="#777777", linewidth=0.9, linestyle="--")
    ax.set_xticks([0, 1], [r["species"].capitalize() for r in sorted(deltas, key=lambda x: x["species"])])
    ax.set_ylabel("Mapped-member − single-gene AUROC", fontsize=FS["axis_label"])
    ax.set_title("Mapped-member versus single-gene", fontsize=FS["title"], loc="left")
    ax.set_ylim(-0.018, 0.068)
    clean_axes(ax)

    ax = axs[1, 1]
    panel_label(ax, "D", x=-0.17, y=1.08)
    for i, r in enumerate(sorted(random, key=lambda x: x["species"])):
        median, q95 = fl(r["random_auc_median"]), fl(r["random_auc_q95"])
        observed = fl(idx[(r["species"], "complete_family")]["auroc"])
        ax.plot([i, i], [median, q95], color=C["unknown"], linewidth=7, solid_capstyle="round", zorder=1)
        ax.scatter(i, median, color="#555555", s=24, label="Random median" if i == 0 else None, zorder=3)
        ax.scatter(i, q95, color=C["primary_fail"], marker="_", s=120, linewidth=2, label="Random 95th percentile" if i == 0 else None, zorder=3)
        ax.scatter(i, observed, color=C["primary_pass"], marker="*", s=90, label="Observed mapped-member score" if i == 0 else None, zorder=3)
        ax.text(i, observed + 0.025, f"P={fl(r['empirical_p']):.4f}", ha="center", fontsize=FS["annotation"])
    ax.set_xticks([0, 1], [r["species"].capitalize() for r in sorted(random, key=lambda x: x["species"])])
    ax.set_ylim(0.44, 1.02)
    ax.set_ylabel("AUROC", fontsize=FS["axis_label"])
    ax.set_title("Random-programme null", fontsize=FS["title"], loc="left")
    ax.legend(frameon=False, fontsize=FS["legend"], loc="center right")
    clean_axes(ax)
    return save(fig, "figure4_controls_and_comparators")


def figure5() -> dict[str, object]:
    full = sorted(read_tsv(SOURCE["figure5"]), key=lambda r: int(r["frozen_rank"]))
    data = [r for r in full if r["function_descriptions"]]
    n = len(data)
    fig = plt.figure(figsize=(WIDTH, LAYOUT["figure5"]["height_inches"]), constrained_layout=False)
    fig.suptitle("Annotated orthogroups in the fixed xylem programme", fontsize=FS["title"] + 1, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.87, bottom=0.14, wspace=0.34)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.08, 0.34, 2.55])
    labels = [f"{r['frozen_rank']:>2}. {r['orthogroup_id']}" for r in data]
    species_fields = ["arabidopsis_members", "esa_members", "sir_members", "spa_members"]
    matrix = np.array([[1 if r[f] else 0 for f in species_fields] for r in data])

    ax = fig.add_subplot(gs[0, 0])
    ax.text(-0.18, 1.06, "A", transform=ax.transAxes, fontsize=FS["panel_label"], fontweight="bold", va="top")
    ax.imshow(matrix, cmap=matplotlib.colors.ListedColormap(["#F0F0F0", C["primary_pass"]]), aspect="auto", vmin=0, vmax=1)
    ax.set_yticks(np.arange(n), labels)
    ax.set_xticks(np.arange(4), ["Ath", "Esa", "Sir", "Spa"])
    ax.set_title("Gene mapping", fontsize=FS["title"], loc="left")
    ax.tick_params(labelsize=FS["tick"], length=0)
    for edge in np.arange(-0.5, n + 0.5, 1):
        ax.axhline(edge, color="white", linewidth=0.6)
    for edge in np.arange(-0.5, 4.5, 1):
        ax.axvline(edge, color="white", linewidth=0.6)
    mapped = sum(r["mapped_in_all_three_validation_species"].lower() == "true" for r in data)
    ax.text(0.0, -0.095, f"{mapped}/{n} annotated orthogroups map in all three validation species", transform=ax.transAxes, fontsize=FS["annotation"], fontweight="bold")

    ax = fig.add_subplot(gs[0, 1], sharey=ax)
    ax.text(-0.38, 1.06, "B", transform=ax.transAxes, fontsize=FS["panel_label"], fontweight="bold", va="top")
    tiers = np.ones((n, 1))
    ax.imshow(tiers, cmap=matplotlib.colors.ListedColormap([C["unknown"], C["primary_pass"]]), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks([0], ["Annotation\navailable"])
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=FS["tick"], length=0)
    ax.set_title("Class", fontsize=FS["title"], loc="center")
    for edge in np.arange(-0.5, n + 0.5, 1):
        ax.axhline(edge, color="white", linewidth=0.6)

    ax = fig.add_subplot(gs[0, 2], sharey=ax)
    panel_label(ax, "C")
    ax.set_xlim(0, 1)
    ax.set_ylim(n - 0.5, -0.5)
    ax.axis("off")
    ax.set_title("Descriptive functional annotation", fontsize=FS["title"], loc="left")
    for i, r in enumerate(data):
        desc = r["function_descriptions"]
        short = textwrap.shorten(desc.replace(" | ", "; "), width=72, placeholder="…")
        ax.text(0.0, i, short, va="center", fontsize=6.2, color=C["text"])
        ax.axhline(i + 0.5, color=C["grid"], linewidth=0.5)
    ax.text(0, n + 0.05, "Annotated subset only; descriptions are non-causal and do not establish paralog replacement.", fontsize=FS["annotation"], fontweight="bold", color=C["negative_control"])
    return save(fig, "figure5_xylem_program")


def supplementary_s1() -> dict[str, object]:
    data = read_tsv(SOURCE["figure1"])
    datasets = [r for r in data if r["record_type"] == "dataset"]
    species_cells: dict[str, int] = defaultdict(int)
    for r in datasets:
        species_cells[r["species"]] += int(r["cells"])
    top = sorted(species_cells.items(), key=lambda x: x[1], reverse=True)[:10][::-1]
    fig, ax = plt.subplots(figsize=(WIDTH, LAYOUT["supplementary_s1"]["height_inches"]), constrained_layout=True)
    ax.barh([x[0] for x in top], [x[1] for x in top], color=C["resource"], zorder=2)
    ax.set_xlabel("Cells", fontsize=FS["axis_label"])
    ax.set_title("Supplementary Figure S1 | Largest species collections", fontsize=FS["title"], loc="left")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(6, 6))
    clean_axes(ax, "x")
    return save(fig, "supplementary_figure_s1_species_collections", source_keys=["figure1"], layout_key="supplementary_s1", manifest_key="supplementary_s1")


def supplementary_s2() -> dict[str, object]:
    data = sorted(read_tsv(SOURCE["figure5"]), key=lambda r: int(r["frozen_rank"]))
    n = len(data)
    fig = plt.figure(figsize=(WIDTH, LAYOUT["supplementary_s2"]["height_inches"]), constrained_layout=False)
    fig.suptitle("Supplementary Figure S2 | Complete xylem mapping matrix", fontsize=FS["title"] + 1, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.22, right=0.90, top=0.91, bottom=0.10, wspace=0.35)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.24])
    labels = [f"{r['frozen_rank']:>2}. {r['orthogroup_id']}" for r in data]
    species_fields = ["arabidopsis_members", "esa_members", "sir_members", "spa_members"]
    matrix = np.array([[1 if r[f] else 0 for f in species_fields] for r in data])

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.imshow(matrix, cmap=matplotlib.colors.ListedColormap(["#F0F0F0", C["primary_pass"]]), aspect="auto", vmin=0, vmax=1)
    ax.set_yticks(np.arange(n), labels)
    ax.set_xticks(np.arange(4), ["Ath", "Esa", "Sir", "Spa"])
    ax.set_title("All 18 fixed orthogroups", fontsize=FS["title"], loc="left")
    ax.tick_params(labelsize=FS["tick"], length=0)
    for edge in np.arange(-0.5, n + 0.5, 1):
        ax.axhline(edge, color="white", linewidth=0.6)
    for edge in np.arange(-0.5, 4.5, 1):
        ax.axvline(edge, color="white", linewidth=0.6)
    ax.text(0.0, -0.07, "16/18 orthogroups mapped in all three validation species", transform=ax.transAxes, fontsize=FS["annotation"], fontweight="bold")

    ax = fig.add_subplot(gs[0, 1], sharey=ax)
    panel_label(ax, "B")
    classes = np.array([[1 if r["function_descriptions"] else 0] for r in data])
    ax.imshow(classes, cmap=matplotlib.colors.ListedColormap([C["unknown"], C["primary_pass"]]), aspect="auto", vmin=0, vmax=1)
    ax.set_xticks([0], ["Annotation\nstatus"])
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", labelsize=FS["tick"], length=0)
    ax.set_title("Class", fontsize=FS["title"], loc="center")
    for edge in np.arange(-0.5, n + 0.5, 1):
        ax.axhline(edge, color="white", linewidth=0.6)
    ax.legend(handles=[Patch(facecolor=C["primary_pass"], label="Annotation available"), Patch(facecolor=C["unknown"], label="Mapping only")], frameon=False, fontsize=FS["legend"], loc="upper left", bbox_to_anchor=(1.05, 1.0))
    return save(fig, "supplementary_figure_s2_full_xylem_program", source_keys=["figure5"], layout_key="supplementary_s2", manifest_key="supplementary_s2")


CAPTIONS = """# Figure captions and concise alt text — manuscript_figures_v1_1

## Figure 1. From audited resource collection to frozen validation

**Caption.** (A) The audited collection contained 121 readable H5AD files and 3,120,007 cells across 33 inferred species; 970,645 cells (31.1%) had unknown or missing labels. (B) Hierarchical label resolution separates coarse vascular or stele system labels from resolved xylem or phloem lineages. This hierarchy does not change prespecified endpoint status: wheat xylem remains sensitivity-only and is not a primary pass. (C) Dataset eligibility was frozen separately for root, leaf and vascular panels; green segments denote candidate inclusion and grey segments denote exclusion under the prespecified governance rules. (D) Analysis sequence from resource audit and discovery through pre-score freezing, independent validation and prospective GSE268881 xylem validation.

**Alt text.** Four-panel figure showing 3.12 million audited cells, hierarchical resolution from vascular or stele labels to xylem or phloem lineages, panel-specific frozen inclusion, and the sequence from discovery to validation; wheat xylem is explicitly sensitivity-only.

## Figure 2. Independent validation is endpoint-dependent

**Caption.** (A) Pooled AUROC for the two passing Moricandia primary endpoints, the failed wheat stele primary endpoint and the wheat xylem sensitivity endpoint. (B) Replicate-level AUROC. (C) Difference between complete-family and deterministic single-gene AUROC. (D) Endpoint sample sizes, average precision, standardized mean difference and empirical random-program P values. Wheat xylem remained sensitivity-only because prespecified replicate-level cell-count criteria were not met; it does not count as a primary pass.

**Alt text.** Pooled and replicate AUROCs show two primary Moricandia passes, a failed wheat stele primary endpoint and a strong but sensitivity-only wheat xylem result; family-versus-single-gene differences vary by endpoint.

## Figure 3. Prospective replication of the frozen xylem program

**Caption.** (A) Pooled AUROC and 95% conditional intervals for eligible *Eutrema salsugineum* (Esa) and *Sisymbrium irio* (Sir). (B) Replicate-level AUROC and intervals. (C) Author-labelled xylem cell counts used for expression-blind eligibility, with the frozen minimum of 100 cells per biological replicate. (D) Prespecified interpretation boundary. *Schrenkiella parvula* (Spa) was excluded before scoring because replicate R1 contained 57 xylem cells; pooling, relabelling or stress-condition substitution was prohibited. Intervals are percentile intervals from 2,000 stratified cell-bootstrap resamples, conditional on the observed biological replicates, and are not population-level biological-replicate confidence intervals.

**Alt text.** Esa and Sir have pooled AUROCs near 0.90 and positive results in both replicates; Spa is shown only as a prespecified eligibility failure because one replicate had 57 xylem cells, below the threshold of 100.

## Figure 4. Controls and method comparisons define the interpretation boundary

**Caption.** (A) The anti-circularity program retained high AUROC, whereas the phloem off-target control was near or below chance. (B) Complete-family, equal-weight orthogroup-presence and deterministic single-gene comparators. (C) Paired complete-family minus single-gene AUROC differences with 95% conditional intervals. The family advantage was positive for Esa but not supported for Sir, precluding a claim of universal family superiority. (D) Observed complete-family AUROCs relative to the median and 95th percentile of 500 frozen random programs per species. Conditional intervals are defined with respect to the observed biological replicates.

**Alt text.** Negative controls behave as expected, but family, equal-weight and single-gene rankings differ between Esa and Sir; paired intervals support a family advantage only in Esa, while both observed scores exceed random-program nulls.

## Figure 5. Annotated subset of the frozen xylem program

**Caption.** This main-text panel is the annotated subset of the frozen xylem program; the complete 18-orthogroup matrix is shown in Supplementary Figure S2. (A) Presence of Arabidopsis discovery members and mapped genes in Esa, Sir and Spa for orthogroups with annotation evidence. (B) All displayed rows are conservatively classified as annotation available; A/B program-specific priority tiers are not transferred across cell programs. (C) Available descriptive functional annotations. These annotations provide biological context but do not establish causal gene function, evolutionary replacement or paralog substitution.

**Alt text.** An annotated subset of the frozen xylem orthogroups is shown with cross-species mapping and descriptive functions; no causal, evolutionary replacement or paralog-substitution claim is made.

## Supplementary Figure S1. Largest audited species collections

**Caption.** Cell totals for the ten largest species collections in the audited resource. This resource-scale view is descriptive and is not an endpoint-performance result.

**Alt text.** Horizontal bars compare audited cell totals for the ten largest plant species collections.

## Supplementary Figure S2. Complete 18-orthogroup frozen xylem mapping matrix

**Caption.** (A) Presence of Arabidopsis discovery members and mapped genes in Esa, Sir and Spa for all 18 frozen orthogroups; 16 orthogroups (88.9%) were mapped in all three validation species. (B) Conservative evidence class distinguishing annotation-available from mapping-only orthogroups. A/B program-specific priority tiers are not transferred across cell programs, and no new xylem-specific annotation review was performed. Functional annotations remain descriptive and non-causal.

**Alt text.** The complete 18-row frozen xylem matrix shows mapping across Arabidopsis, Esa, Sir and Spa, with a side column distinguishing annotation-available and mapping-only entries.
"""


def verify_inputs() -> list[str]:
    return [f"Missing figure input: {path}" for path in SOURCE.values() if not path.is_file()]


def write_manifest(figures: list[dict[str, object]]) -> Path:
    path = OUT / "figure_source_manifest.tsv"
    fields = ["figure", "source_data_absolute_path", "source_data_sha256", "plot_script_absolute_path", "plot_script_sha256", "design_freeze_absolute_path", "design_freeze_sha256", "png_absolute_path", "png_sha256", "png_bytes", "png_width_px", "png_height_px", "dpi", "pdf_absolute_path", "pdf_sha256", "pdf_bytes", "width_inches", "height_inches"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for meta in figures:
            source_paths = [SOURCE[key] for key in meta["source_keys"]]
            png, pdf = meta["png"], meta["pdf"]
            writer.writerow({
                "figure": meta["manifest_key"],
                "source_data_absolute_path": ";".join(str(p) for p in source_paths),
                "source_data_sha256": ";".join(sha256(p) for p in source_paths),
                "plot_script_absolute_path": SCRIPT, "plot_script_sha256": sha256(SCRIPT),
                "design_freeze_absolute_path": DESIGN, "design_freeze_sha256": sha256(DESIGN),
                "png_absolute_path": png, "png_sha256": sha256(png), "png_bytes": png.stat().st_size,
                "png_width_px": meta["png_width_px"], "png_height_px": meta["png_height_px"], "dpi": DPI,
                "pdf_absolute_path": pdf, "pdf_sha256": sha256(pdf), "pdf_bytes": pdf.stat().st_size,
                "width_inches": WIDTH, "height_inches": LAYOUT[meta["layout_key"]]["height_inches"],
            })
    return path


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    global CFG, C, FS, LAYOUT, WIDTH, DPI
    CFG = json.loads(DESIGN.read_text(encoding="utf-8"))
    C, FS, LAYOUT = CFG["palette"], CFG["font_sizes"], CFG["layouts"]
    WIDTH, DPI = CFG["width_inches"], CFG["dpi"]
    plt.rcParams.update({
        "font.family": CFG["font_family"], "font.size": FS["tick"], "axes.labelcolor": C["text"],
        "axes.titlecolor": C["text"], "text.color": C["text"], "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "axes.linewidth": 0.7,
    })
    conflicts = verify_inputs()
    if conflicts:
        (OUT / "audit.json").write_text(json.dumps({"status": "blocked_input_conflict", "conflicts": conflicts}, indent=2) + "\n", encoding="utf-8")
        (OUT / "STATUS").write_text("BLOCKED_INPUT_CONFLICT\n", encoding="utf-8")
        return 2
    # Figures 1-3 use compact submission layouts generated by
    # plot_figure1_v1_2.py and build_main_figures_v1_2.py. This script
    # regenerates the remaining main and supplementary figures.
    figures = [figure4(), figure5(), supplementary_s1(), supplementary_s2()]
    captions = OUT / "figure_captions_and_alt_text.md"
    captions.write_text(CAPTIONS, encoding="utf-8")
    manifest = write_manifest(figures)
    generated = [p for f in figures for p in (f["png"], f["pdf"])] + [captions, manifest, DESIGN, SCRIPT]
    audit = {
        "status": "complete", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": f"{sys.executable} {SCRIPT}",
        "environment": {"python": sys.version, "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pillow": Image.__version__, "platform": platform.platform()},
        "input_directory": str(EVIDENCE), "input_status": "COMPLETE", "input_conflicts": 0,
        "input_files": [{"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} for p in SOURCE.values()],
        "design_freeze": {"path": str(DESIGN), "sha256": sha256(DESIGN), "frozen_before_rendering": True},
        "generated_files": [{"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} for p in generated],
        "conflicts": [],
        "manuscript_decisions_applied": [
            "Moved the top-ten species collection panel from Figure 1 to Supplementary Figure S1.",
            "Replaced Figure 1 panel B with hierarchical label resolution from vascular/stele to xylem/phloem and explicitly retained wheat xylem as sensitivity-only.",
            "Restricted main Figure 5 to the annotation-available subset and moved the complete 18-orthogroup matrix to Supplementary Figure S2.",
            "Retained annotation available versus mapping only; did not use A/B tiers or initiate a new xylem-specific annotation review.",
            "Moved the Figure 4 panel B method legend outside the plotting area."
        ],
        "decisions_for_manuscript_agent": [],
        "interpretation_guards": [
            "Wheat xylem is sensitivity-only and is not counted as a primary pass.",
            "Spa is eligibility-excluded and has no score.",
            "Family-versus-single-gene performance is not universally superior.",
            "Functional annotations are descriptive and non-causal.",
            "Bootstrap intervals are conditional on the observed biological replicates."
        ],
        "safety": {"new_data_downloaded": False, "training_started": False, "tuning_started": False, "thresholds_changed": False, "tex_modified": False, "authority_outputs_modified": False},
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUT / "STATUS").write_text("COMPLETE\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
