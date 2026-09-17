#!/usr/bin/env python3
"""Render PlantHOP Figure 1 v1.2 from frozen local evidence."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from reportlab.lib.colors import HexColor, white
from reportlab.pdfgen import canvas


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
RESOURCE = PROJECT / "manuscript_evidence_v1" / "figure1_resource_flow.tsv"
CANDIDATES = PROJECT / "conserved_program_candidates_v2_audit.json"
OUT_PDF = HERE / "figure1_resource_to_validation.pdf"
AUDIT = HERE / "figure1_v1_2_audit.json"

PAGE_W, PAGE_H = 7.2 * 72, 7.25 * 72
C = {
    "navy": HexColor("#315C78"), "blue": HexColor("#4C78A8"),
    "green": HexColor("#59A14F"), "purple": HexColor("#8E6C9F"),
    "gold": HexColor("#E2B04A"), "unknown": HexColor("#D9D9D9"),
    "text": HexColor("#24323D"), "grid": HexColor("#E8ECEF"),
    "pale": HexColor("#F4F6F8"), "grey": HexColor("#777777"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def text(c, x, y, value, size=8, bold=False, color=None, anchor="start"):
    c.setFillColor(color or C["text"])
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    if anchor == "middle":
        c.drawCentredString(x, y, value)
    elif anchor == "end":
        c.drawRightString(x, y, value)
    else:
        c.drawString(x, y, value)


def title(c, x, y, panel, value):
    text(c, x - 15, y + 2, panel, 12, True)
    text(c, x, y, value, 10.4, True)


def rounded_box(c, x, y, w, h, fill, stroke=None, radius=6):
    c.setFillColor(fill)
    c.setStrokeColor(stroke or fill)
    c.setLineWidth(0.9 if stroke else 0)
    c.roundRect(x, y, w, h, radius, fill=1, stroke=1 if stroke else 0)


def arrow_down(c, x, y1, y2):
    c.setStrokeColor(C["grey"]); c.setFillColor(C["grey"]); c.setLineWidth(0.9)
    c.line(x, y1, x, y2 + 4)
    p = c.beginPath(); p.moveTo(x - 3, y2 + 5); p.lineTo(x + 3, y2 + 5); p.lineTo(x, y2); p.close()
    c.drawPath(p, fill=1, stroke=0)


with RESOURCE.open("r", encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
summary = {r["metric"]: r for r in rows if r["record_type"] == "collection_summary"}
with CANDIDATES.open("r", encoding="utf-8") as handle:
    candidate_audit = json.load(handle)

total = int(summary["cells"]["value"])
labeled = int(summary["labeled_cells"]["value"])
unknown = int(summary["unknown_or_missing_cells"]["value"])
counts = candidate_audit["candidate_counts_by_panel_and_class"]
assert (total, labeled, unknown) == (3120007, 2149362, 970645)
assert candidate_audit["conserved_candidate_rows"] == 102
assert sum(counts.values()) == 102

c = canvas.Canvas(str(OUT_PDF), pagesize=(PAGE_W, PAGE_H))
c.setTitle("PlantHOP Figure 1 v1.2")
c.setAuthor("PlantHOP")
margin, gap = 30, 24
panel_w = (PAGE_W - 2 * margin - gap) / 2
x1, x2 = margin, margin + panel_w + gap
y_top, y_bottom = PAGE_H - 35, 246

# A: resource overview
title(c, x1, y_top, "A", "Resource overview")
text(c, x1, y_top - 20, f"121 readable H5AD  |  {total:,} cells  |  33 species", 7.8)
bar_y, bar_h, bar_w = y_top - 82, 32, panel_w
lab_w = bar_w * labeled / total
c.setFillColor(C["navy"]); c.rect(x1, bar_y, lab_w, bar_h, fill=1, stroke=0)
c.setFillColor(C["unknown"]); c.rect(x1 + lab_w, bar_y, bar_w - lab_w, bar_h, fill=1, stroke=0)
text(c, x1 + lab_w / 2, bar_y + 18, f"{labeled:,}", 7.4, True, white, "middle")
text(c, x1 + lab_w / 2, bar_y + 8, "68.9% labelled", 6.7, False, white, "middle")
text(c, x1 + lab_w + (bar_w - lab_w) / 2, bar_y + 18, f"{unknown:,}", 7.4, True, C["text"], "middle")
text(c, x1 + lab_w + (bar_w - lab_w) / 2, bar_y + 8, "31.1% unknown", 6.7, False, C["text"], "middle")
c.setStrokeColor(C["text"]); c.setLineWidth(0.7); c.line(x1, bar_y - 12, x1 + panel_w, bar_y - 12)
for frac, lab in [(0, "0"), (0.5, "1.56 M"), (1, "3.12 M")]:
    xx = x1 + frac * panel_w
    c.line(xx, bar_y - 12, xx, bar_y - 8)
    text(c, xx, bar_y - 23, lab, 6.7, anchor="middle")
text(c, x1 + panel_w / 2, bar_y - 37, "Reviewed cells", 7.4, anchor="middle")

# B: harmonized hierarchy
title(c, x2, y_top, "B", "Harmonized label hierarchy")
text(c, x2, y_top - 20, "Tissue panel", 7.2, True)
text(c, x2 + 105, y_top - 20, "Harmonized cell systems", 7.2, True)
hierarchy = [
    ("Leaf", "photosynthetic / ground  |  epidermal", C["blue"]),
    ("Root", "epidermal  |  ground  |  root cap  |  stele", C["green"]),
    ("Vascular", "xylem  |  phloem", C["purple"]),
]
for i, (panel, systems, color) in enumerate(hierarchy):
    yy = y_top - 65 - i * 48
    rounded_box(c, x2, yy, 70, 26, color)
    text(c, x2 + 35, yy + 9, panel, 7.4, True, white, "middle")
    c.setStrokeColor(C["grey"]); c.setLineWidth(0.8); c.line(x2 + 72, yy + 13, x2 + 96, yy + 13)
    rounded_box(c, x2 + 98, yy, panel_w - 98, 26, C["pale"], color)
    text(c, x2 + 98 + (panel_w - 98) / 2, yy + 9, systems, 6.5, False, C["text"], "middle")

# C: recurrent programme counts
title(c, x1, y_bottom, "C", "Recurrent programme candidates (n = 102)")
items = [
    ("Photosynthetic / ground [leaf]", 30, C["blue"]),
    ("Epidermal [leaf]", 9, C["blue"]),
    ("Epidermal [root]", 5, C["green"]),
    ("Ground tissue [root]", 8, C["green"]),
    ("Root cap [root]", 1, C["green"]),
    ("Stele [root]", 14, C["green"]),
    ("Phloem [vascular]", 17, C["purple"]),
    ("Xylem [vascular]", 18, C["purple"]),
]
label_w, max_bar = 116, panel_w - 136
start_y, row_h = y_bottom - 27, 20
for tick in [0, 10, 20, 30]:
    xx = x1 + label_w + max_bar * tick / 30
    c.setStrokeColor(C["grid"]); c.setLineWidth(0.55); c.line(xx, start_y + 15, xx, start_y - 7 * row_h - 7)
    text(c, xx, start_y - 7 * row_h - 19, str(tick), 6.4, anchor="middle")
for i, (label, value, color) in enumerate(items):
    yy = start_y - i * row_h
    text(c, x1 + label_w - 4, yy + 3, label, 6.5, anchor="end")
    c.setFillColor(color); c.rect(x1 + label_w, yy, max_bar * value / 30, 11, fill=1, stroke=0)
    text(c, x1 + label_w + max_bar * value / 30 + 4, yy + 2.5, str(value), 6.7, True)
text(c, x1 + label_w + max_bar / 2, start_y - 7 * row_h - 31, "Candidate orthogroups", 7.1, anchor="middle")
text(c, x1, 30, "Top 100 positive differences per species-class; recurrence in >=3 species", 6.2)

# D: analysis sequence
title(c, x2, y_bottom, "D", "Discovery and independent validation")
steps = [
    ("Programme discovery", "dataset review | hierarchy | orthogroups", C["navy"]),
    ("Prespecified design", "datasets | labels | programmes | thresholds", C["gold"]),
    ("Independent tests", "rice | Moricandia | wheat", C["green"]),
    ("Held-out xylem test", "Esa | Sir; Spa excluded", C["purple"]),
]
box_x, box_w, box_h = x2 + 15, panel_w - 30, 34
for i, (head, detail, color) in enumerate(steps):
    yy = y_bottom - 54 - i * 51
    rounded_box(c, box_x, yy, box_w, box_h, color)
    txt_color = C["text"] if color == C["gold"] else white
    text(c, box_x + box_w / 2, yy + 20, head, 7.8, True, txt_color, "middle")
    text(c, box_x + box_w / 2, yy + 8, detail, 6.6, False, txt_color, "middle")
    if i < len(steps) - 1:
        arrow_down(c, box_x + box_w / 2, yy - 2, yy - 20)

c.showPage(); c.save()

audit = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "status": "manuscript_figure1_v1_2",
    "inputs": {str(RESOURCE): sha256(RESOURCE), str(CANDIDATES): sha256(CANDIDATES)},
    "checks": {
        "total_cells": total, "unknown_cells": unknown,
        "candidate_programs": sum(counts.values()),
        "candidate_counts_by_panel_and_class": counts,
    },
    "outputs": {str(OUT_PDF): sha256(OUT_PDF)},
}
with AUDIT.open("w", encoding="utf-8") as handle:
    json.dump(audit, handle, indent=2, ensure_ascii=False)
print(json.dumps(audit, indent=2, ensure_ascii=False))
