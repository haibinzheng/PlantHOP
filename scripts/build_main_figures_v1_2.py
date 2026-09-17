#!/usr/bin/env python3
"""Build typography-optimized PlantHOP main Figures 2-5 as vector PDFs."""

from __future__ import annotations

import csv
import hashlib
import json
import textwrap
from datetime import datetime, timezone
from pathlib import Path

from reportlab.lib.colors import HexColor, white
from reportlab.pdfgen import canvas


REPOSITORY = Path(__file__).resolve().parents[1]
HERE = REPOSITORY / "results" / "figures"
EVIDENCE = REPOSITORY / "results" / "figure_data"
SOURCES = {i: EVIDENCE / f"figure{i}_{name}.tsv" for i, name in {
    2: "independent_validation", 3: "prospective_xylem",
    4: "controls_comparators", 5: "xylem_program",
}.items()}
W = 7.2 * 72
COL = {
    "pass": HexColor("#228833"), "fail": HexColor("#EE7733"),
    "sens": HexColor("#4477AA"), "complete": HexColor("#332288"),
    "equal": HexColor("#55BFE2"), "single": HexColor("#CCBB44"),
    "anti": HexColor("#44AA99"), "off": HexColor("#AA3377"),
    "grey": HexColor("#777777"), "light": HexColor("#E4E7EA"),
    "pale": HexColor("#F4F6F8"), "text": HexColor("#222222"),
    "excluded": HexColor("#999999"),
}


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def txt(c, x, y, s, size=8.2, bold=False, color=None, align="left"):
    c.setFillColor(color or COL["text"])
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    if align == "center": c.drawCentredString(x, y, s)
    elif align == "right": c.drawRightString(x, y, s)
    else: c.drawString(x, y, s)


def panel(c, x, y, label, title):
    txt(c, x - 18, y, label, 12, True)
    txt(c, x, y, title, 10.4, True)


def wrap(c, x, y, s, width, size=7.7, leading=10, bold=False, color=None):
    for i, line in enumerate(textwrap.wrap(s, width=width)):
        txt(c, x, y - i * leading, line, size, bold, color)


def axes(c, x, y, w, h, lo, hi, ylabel, ticks):
    c.setStrokeColor(COL["text"]); c.setLineWidth(0.75)
    c.line(x, y, x, y + h); c.line(x, y, x + w, y)
    for t in ticks:
        yy = y + h * (t - lo) / (hi - lo)
        c.setStrokeColor(COL["light"]); c.setLineWidth(0.55); c.line(x, yy, x + w, yy)
        txt(c, x - 6, yy - 2.5, f"{t:.2f}" if hi <= 1 else f"{t:g}", 7.2, align="right")
    c.saveState(); c.translate(x - 33, y + h / 2); c.rotate(90)
    txt(c, 0, 0, ylabel, 8.2, align="center"); c.restoreState()


def dot(c, x, y, r, color, stroke=None):
    c.setFillColor(color); c.setStrokeColor(stroke or color); c.circle(x, y, r, fill=1, stroke=1)


def bar(c, x, y, w, h, color, hatch=False):
    c.setFillColor(color); c.setStrokeColor(color); c.rect(x, y, w, h, fill=1, stroke=0)
    if hatch:
        c.saveState(); p = c.beginPath(); p.rect(x, y, w, h); c.clipPath(p, stroke=0)
        c.setStrokeColor(HexColor("#1D3557")); c.setLineWidth(0.8)
        for xx in range(int(x - h), int(x + w + h), 8): c.line(xx, y, xx + h, y + h)
        c.restoreState()


def finish(c):
    c.showPage(); c.save()


def figure2(rows):
    path = HERE / "figure2_independent_validation.pdf"
    c = canvas.Canvas(str(path), pagesize=(W, 5.75 * 72))
    pooled = [r for r in rows if r["scope"] == "pooled"]
    reps = [r for r in rows if r["scope"] == "replicate"]
    ids = ["moricandia_leaf_photosynthetic_v1", "moricandia_leaf_epidermal_v1", "wheat_root_stele_v1", "wheat_root_xylem_v1"]
    names = ["M.\nphoto", "M.\nepidermal", "Wheat\nstele", "Wheat\nxylem"]
    idx = {r["contrast_id"]: r for r in pooled}
    colors = [COL["pass"], COL["pass"], COL["fail"], COL["sens"]]
    xs = [38, 282]; top = 5.75 * 72 - 31; pw = 198; ph = 122
    # A
    panel(c, xs[0], top, "A", "Pooled endpoint outcomes")
    axx, ayy, aw, ah = xs[0] + 26, top - 150, 166, 112
    axes(c, axx, ayy, aw, ah, 0.5, 1.0, "Pooled AUROC", [0.5, .6, .7, .8, .9, 1.0])
    for i, cid in enumerate(ids):
        v = float(idx[cid]["auroc"]); bw = 29; bx = axx + 10 + i * 39
        bar(c, bx, ayy, bw, ah * (v - .5) / .5, colors[i], i == 3)
        txt(c, bx + bw/2, ayy + ah*(v-.5)/.5 + 5, f"{v:.3f}", 7.5, True, align="center")
        for j, line in enumerate(names[i].split("\n")): txt(c, bx + bw/2, ayy - 13 - j*9, line, 7.1, align="center")
    # B
    panel(c, xs[1], top, "B", "Replicate consistency")
    axx, ayy = xs[1] + 26, top - 150
    axes(c, axx, ayy, aw, ah, 0.5, 1.0, "Replicate AUROC", [0.5, .6, .7, .8, .9, 1.0])
    for i, cid in enumerate(ids):
        rr = [r for r in reps if r["contrast_id"] == cid]
        base = axx + 18 + i * 42
        pts = []
        for j, r in enumerate(rr):
            xx = base + (j - (len(rr)-1)/2)*7; v = float(r["auroc"]); yy = ayy + ah*(v-.5)/.5
            pts.append((xx,yy)); dot(c, xx, yy, 3.1, colors[i])
        if len(pts) > 1:
            c.setStrokeColor(colors[i]); c.setLineWidth(1); c.line(pts[0][0], pts[0][1], pts[-1][0], pts[-1][1])
        txt(c, base, ayy - 13, ["M-photo","M-epi","W-stele","W-xylem"][i], 7.1, align="center")
    # C
    bottom = 184
    panel(c, xs[0], bottom, "C", "Family minus single-gene AUROC")
    axx, ayy, aw, ah = xs[0] + 26, 38, 166, 112
    axes(c, axx, ayy, aw, ah, 0, .42, "AUROC difference", [0,.1,.2,.3,.4])
    for i, cid in enumerate(ids):
        v = float(idx[cid]["family_minus_single_auroc"]); bx = axx + 10 + i*39; bw=29
        bar(c, bx, ayy, bw, ah*v/.42, colors[i])
        txt(c, bx+bw/2, ayy+ah*v/.42+5, f"{v:+.3f}", 7.5, True, align="center")
        txt(c, bx+bw/2, ayy-13, ["M-photo","M-epi","W-stele","W-xylem"][i], 7.1, align="center")
    # D: readable endpoint table
    panel(c, xs[1], bottom, "D", "Endpoint metrics and status")
    x, y, w, rh = xs[1], 40, 205, 25
    headers = ["Endpoint", "AUROC", "AP", "SMD", "Status"]
    widths = [78, 31, 28, 32, 38]
    c.setFillColor(COL["pale"]); c.rect(x, y+4*rh, sum(widths), 21, fill=1, stroke=0)
    xx=x
    for h,ww in zip(headers,widths): txt(c,xx+3,y+4*rh+7,h,7.2,True); xx+=ww
    short=["M. photo","M. epidermal","Wheat stele","Wheat xylem"]
    status=["Pass","Pass","Fail","Sensitivity"]
    for i,cid in enumerate(ids):
        r=idx[cid]; yy=y+(3-i)*rh; c.setStrokeColor(COL["light"]); c.line(x,yy,x+sum(widths),yy)
        vals=[short[i],f"{float(r['auroc']):.3f}",f"{float(r['average_precision']):.3f}",f"{float(r['smd']):.2f}",status[i]]
        xx=x
        for j,(v,ww) in enumerate(zip(vals,widths)):
            txt(c,xx+3,yy+8,v,7.0,j==4,colors[i] if j==4 else None); xx+=ww
    finish(c); return path


def figure3(rows):
    path = HERE / "figure3_prospective_xylem.pdf"
    c = canvas.Canvas(str(path), pagesize=(W, 5.45*72)); top=5.45*72-30; xs=[38,282]; aw=174; ah=108
    pooled={r["species"]:r for r in rows if r["scope"]=="pooled"}; reps=[r for r in rows if r["scope"]=="replicate"]
    # A
    panel(c,xs[0],top,"A","Pooled held-out performance")
    axx,ayy=xs[0]+30,top-144; axes(c,axx,ayy,aw,ah,.5,1,"Pooled AUROC",[.5,.6,.7,.8,.9,1])
    for i,sp in enumerate(["esa","sir"]):
        r=pooled[sp]; v=float(r["auroc"]); lo=float(r["auroc_lower_95"]); hi=float(r["auroc_upper_95"]); xx=axx+52+i*72
        y=lambda z: ayy+ah*(z-.5)/.5
        c.setStrokeColor(COL["pass"] if sp=="esa" else COL["complete"]); c.setLineWidth(1.2); c.line(xx,y(lo),xx,y(hi)); c.line(xx-4,y(lo),xx+4,y(lo)); c.line(xx-4,y(hi),xx+4,y(hi)); dot(c,xx,y(v),4,COL["pass"] if sp=="esa" else COL["complete"])
        txt(c,xx,y(hi)+7,f"{v:.3f}",7.6,True,align="center"); txt(c,xx,ayy-14,sp.title(),7.6,align="center")
    # B
    panel(c,xs[1],top,"B","Replicate consistency")
    axx,ayy=xs[1]+30,top-144; axes(c,axx,ayy,aw,ah,.5,1,"Replicate AUROC",[.5,.6,.7,.8,.9,1])
    eligible_reps = [r for r in reps if r["auroc"]]
    for i,r in enumerate(eligible_reps):
        v=float(r["auroc"]); lo=float(r["auroc_lower_95"]); hi=float(r["auroc_upper_95"]); xx=axx+24+i*40; color=COL["pass"] if r["species"]=="esa" else COL["complete"]
        y=lambda z: ayy+ah*(z-.5)/.5
        c.setStrokeColor(color); c.line(xx,y(lo),xx,y(hi)); c.line(xx-3,y(lo),xx+3,y(lo)); c.line(xx-3,y(hi),xx+3,y(hi)); dot(c,xx,y(v),3.6,color)
        txt(c,xx,ayy-13,r["species"].title(),7.0,align="center")
        txt(c,xx,ayy-22,r["replicate"],7.0,align="center")
    # C
    by=168; panel(c,xs[0],by,"C","Pre-score eligibility")
    axx,ayy=xs[0]+30,34; axes(c,axx,ayy,aw,100,0,130,"Xylem cells",[0,50,100])
    ordered=[r for r in reps if r["species"] in ("esa","sir","spa")]
    for i,r in enumerate(ordered):
        n=int(r["positive_cells"]); xx=axx+8+i*27; color=COL["pass"] if r["species"]=="esa" else COL["complete"] if r["species"]=="sir" else COL["excluded"]
        bar(c,xx,ayy,18,100*n/130,color); txt(c,xx+9,ayy+100*n/130+5,str(n),7.2,True,align="center")
        txt(c,xx+9,ayy-13,r["species"].upper(),6.7,align="center"); txt(c,xx+9,ayy-22,r["replicate"],6.7,align="center")
    c.setStrokeColor(COL["fail"]); c.setDash(4,2); c.line(axx,ayy+100*100/130,axx+aw,ayy+100*100/130); c.setDash()
    txt(c,axx+aw,ayy+100*100/130+4,"minimum = 100",7.0,True,COL["fail"],"right")
    # D
    panel(c,xs[1],by,"D","Interpretation boundary")
    rounded=lambda yy,color: (c.setFillColor(COL["pale"]),c.setStrokeColor(color),c.roundRect(xs[1],yy,205,48,6,fill=1,stroke=1))
    rounded(88,COL["pass"]); txt(c,xs[1]+10,121,"Esa and Sir",8.8,True,COL["pass"]); wrap(c,xs[1]+10,107,"Eligible and formally scored; both passed the frozen species-level rule.",48,7.4,9)
    rounded(30,COL["excluded"]); txt(c,xs[1]+10,63,"Spa",8.8,True,COL["excluded"]); wrap(c,xs[1]+10,49,"Excluded before scoring: R1 contained 57 xylem cells (<100). No pooled rescue or relabelling.",48,7.4,9)
    finish(c); return path


def figure4(rows):
    path=HERE/"figure4_controls_and_comparators.pdf"; c=canvas.Canvas(str(path),pagesize=(W,5.65*72)); top=5.65*72-30; xs=[38,282]; aw=174; ah=105
    methods={(r["species"],r["comparator"]):r for r in rows if r["scope"]=="pooled" and r["record_type"]=="method"}
    # A controls
    panel(c,xs[0],top,"A","Negative controls")
    axx,ayy=xs[0]+30,top-140; axes(c,axx,ayy,aw,ah,.4,1,"Pooled AUROC",[.4,.5,.6,.7,.8,.9,1])
    comps=[("complete_family",COL["complete"],"Mapped members"),("anti_circular",COL["anti"],"Marker-overlap exclusion"),("phloem_off_target",COL["off"],"Phloem off-target")]
    for gi,sp in enumerate(["esa","sir"]):
        for j,(comp,color,_) in enumerate(comps):
            v=float(methods[(sp,comp)]["auroc"]); bx=axx+15+gi*83+j*18; bar(c,bx,ayy,15,ah*(v-.4)/.6,color)
        txt(c,axx+34+gi*83,ayy-14,sp.title(),7.6,align="center")
    for j,(_,color,label) in enumerate(comps): bar(c,xs[0]+j*66,top-25,10,7,color); txt(c,xs[0]+14+j*66,top-24,label,6.8)
    # B method comparison
    panel(c,xs[1],top,"B","Method comparisons")
    axx,ayy=xs[1]+30,top-140; axes(c,axx,ayy,aw,ah,.82,.96,"Pooled AUROC",[.82,.86,.90,.94])
    comps2=[("complete_family",COL["complete"],"Mapped members"),("equal_weight",COL["equal"],"Equal weight"),("single_gene",COL["single"],"Single gene")]
    for gi,sp in enumerate(["esa","sir"]):
        for j,(comp,color,_) in enumerate(comps2):
            v=float(methods[(sp,comp)]["auroc"]); bx=axx+15+gi*83+j*18; bar(c,bx,ayy,15,ah*(v-.82)/.14,color)
        txt(c,axx+34+gi*83,ayy-14,sp.title(),7.6,align="center")
    for j,(_,color,label) in enumerate(comps2): bar(c,xs[1]+j*66,top-25,10,7,color); txt(c,xs[1]+14+j*66,top-24,label,6.8)
    # C paired difference
    by=170; panel(c,xs[0],by,"C","Mapped-member - single-gene difference")
    paired={r["species"]:r for r in rows if r["scope"]=="pooled" and r["record_type"]=="paired_difference" and r["estimand"]=="delta_auc"}
    axx,ayy=xs[0]+30,34; axes(c,axx,ayy,aw,100,-.02,.07,"AUROC difference",[-.02,0,.02,.04,.06])
    y=lambda z: ayy+100*(z+.02)/.09
    c.setStrokeColor(COL["grey"]); c.setDash(4,2); c.line(axx,y(0),axx+aw,y(0)); c.setDash()
    for i,sp in enumerate(["esa","sir"]):
        r=paired[sp]; v=float(r["estimate"]); lo=float(r["lower_95"]); hi=float(r["upper_95"]); xx=axx+52+i*72
        c.setStrokeColor(COL["complete"]); c.setLineWidth(1.2); c.line(xx,y(lo),xx,y(hi)); c.line(xx-4,y(lo),xx+4,y(lo)); c.line(xx-4,y(hi),xx+4,y(hi)); dot(c,xx,y(v),4,COL["complete"]); txt(c,xx+8,y(v)-2,f"{v:+.3f}",7.5,True); txt(c,xx,ayy-14,sp.title(),7.5,align="center")
    # D random null
    panel(c,xs[1],by,"D","Size-matched random-programme null")
    null={r["species"]:r for r in rows if r["record_type"]=="random_null"}
    axx,ayy=xs[1]+30,34; axes(c,axx,ayy,aw,100,.45,.98,"AUROC",[.5,.6,.7,.8,.9])
    y2=lambda z: ayy+100*(z-.45)/.53
    for i,sp in enumerate(["esa","sir"]):
        r=null[sp]; xx=axx+52+i*72; med=float(r["random_auc_median"]); q=float(r["random_auc_q95"]); obs=float(methods[(sp,"complete_family")]["auroc"])
        c.setStrokeColor(COL["excluded"]); c.setLineWidth(5); c.line(xx,y2(med),xx,y2(q)); dot(c,xx,y2(med),3,COL["grey"]); txt(c,xx,y2(obs)-4,"*",15,True,COL["pass"],"center"); txt(c,xx,y2(obs)+10,"P=0.002",7.0,True,align="center"); txt(c,xx,ayy-14,sp.title(),7.5,align="center")
    finish(c); return path


def figure5(rows):
    path=HERE/"figure5_xylem_program.pdf"; c=canvas.Canvas(str(path),pagesize=(W,4.75*72)); H=4.75*72
    annotated=[r for r in rows if r["evidence_tier"] in ("A","B")]
    txt(c,W/2,H-25,"Annotated subset of the fixed xylem programme",12,True,align="center")
    x0,y0=65,48; row_h=27; heat_x=112; cell=24; class_x=229; desc_x=285
    panel(c,heat_x,H-52,"A","Frozen mapping"); panel(c,class_x,H-52,"B","Class"); panel(c,desc_x,H-52,"C","Descriptive functional annotation")
    species=["Ath","Esa","Sir","Spa"]
    for j,s in enumerate(species): txt(c,heat_x+j*cell+cell/2,y0-15,s,7.4,True,align="center")
    for i,r in enumerate(reversed(annotated)):
        yy=y0+i*row_h; label=f"{r['frozen_rank']}. {r['orthogroup_id']}"; txt(c,heat_x-7,yy+8,label,7.1,align="right")
        presence=[bool(r["arabidopsis_members"]),bool(r["esa_members"]),bool(r["sir_members"]),bool(r["spa_members"])]
        for j,p in enumerate(presence):
            c.setFillColor(COL["pass"] if p else COL["pale"]); c.setStrokeColor(white); c.rect(heat_x+j*cell,yy,cell,row_h,fill=1,stroke=1)
        c.setFillColor(COL["pass"]); c.setStrokeColor(white); c.rect(class_x,yy,29,row_h,fill=1,stroke=1)
        desc=r["function_descriptions"].replace(" | ","; ")
        lines=textwrap.wrap(desc,width=52)
        if len(lines)>2: lines=lines[:2]; lines[-1]=lines[-1].rstrip(".; ")+"..."
        if len(lines)==1: txt(c,desc_x,yy+8,lines[0],7.0)
        else:
            txt(c,desc_x,yy+14,lines[0],6.8); txt(c,desc_x,yy+5,lines[1],6.8)
        c.setStrokeColor(COL["light"]); c.line(desc_x,yy,500,yy)
    txt(c,class_x+14.5,y0-15,"Annotation",7.2,True,align="center")
    txt(c,class_x+14.5,y0-24,"available",7.2,True,align="center")
    txt(c,heat_x,y0+len(annotated)*row_h+8,"6/8 annotated orthogroups map in all three validation species",7.6,True)
    txt(c,desc_x,y0-15,"Annotations are descriptive and non-causal;",7.0,True,COL["off"])
    txt(c,desc_x,y0-24,"they do not establish paralog replacement.",7.0,True,COL["off"])
    finish(c); return path


def main():
    HERE.mkdir(parents=True,exist_ok=True)
    # Figures 4-5 and Supplementary Figures S1-S2 use the larger
    # matplotlib submission layouts generated by plot_manuscript_figures_v1_1.py.
    outputs=[figure2(read_tsv(SOURCES[2])),figure3(read_tsv(SOURCES[3]))]
    audit={"status":"complete","generated_at_utc":datetime.now(timezone.utc).isoformat(),"design":{"width_inches":7.2,"minimum_body_font_pt":7.0,"title_font_pt":10.4,"panel_label_pt":12,"format":"vector_pdf"},"inputs":{str(p):sha256(p) for p in SOURCES.values()},"outputs":{str(p):sha256(p) for p in outputs},"guards":["data values read from frozen manuscript evidence TSV files","no thresholds or endpoint status changed","typography and layout only"]}
    (HERE/"audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(json.dumps(audit,indent=2))


if __name__=="__main__": main()
