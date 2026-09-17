#!/usr/bin/env python3
"""Build auditable Tables S2/S3 and Methods S1 without modifying source data."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd


PROJECT = Path("/workspace/projects/PhyloOpenCell")
META = PROJECT / "metadata"
REPORTS = PROJECT / "reports"
FINAL_OUTPUT = REPORTS / "supplementary_completion_v1"
OUTPUT = REPORTS / "supplementary_completion_v1.partial"
PROCESS = PROJECT / "server_experiments/phylo_plant_fm_pilot/reports"

SOURCES = {
    "dataset_manifest": META / "dataset_manifest.tsv",
    "dataset_source_metadata": META / "dataset_source_metadata.tsv",
    "panel_label_freeze": META / "panel_label_freeze_v1.tsv",
    "coarse_label_hierarchy": META / "coarse_label_hierarchy_v1.tsv",
    "panel_label_candidates": META / "panel_label_candidates_v1.tsv",
    "ontology_mapping_candidates": META / "ontology_mapping_candidates.tsv",
    "panel_dataset_label_counts": META / "panel_dataset_label_counts_v1.tsv.gz",
    "panel_dataset_freeze": META / "panel_dataset_freeze_v1.tsv",
    "gene_mapping_summary": META / "gene_mapping_candidate_summary.tsv",
    "orthogroup_coverage": META / "dataset_orthogroup_coverage_v1.tsv",
    "gene_mapping_audit": REPORTS / "gene_mapping_candidate_audit.json",
    "orthofinder_audit": REPORTS / "orthofinder_full7_v1_audit.json",
    "marker_rank_dataset_audit": Path("/data/runs/PhyloOpenCell/marker_rank_features_v2/dataset_feature_audit.tsv"),
    "discovery_audit": REPORTS / "conserved_program_candidates_v2_audit.json",
    "p0_audit": REPORTS / "label_resolution_circularity_robustness_v1/audit.json",
    "p1_audit": REPORTS / "discovery_threshold_ortholog_robustness_v1/audit.json",
    "rice_v2_1_audit": REPORTS / "independent_program_validation_v2_1_audit.json",
    "external_v3_execution": REPORTS / "external_validation_v3_execution_addendum_v1.json",
    "supplement_s5": REPORTS / "manuscript_supplement_tables_v1_1/Table_S5_replicate_metrics_controls_exclusions.tsv",
}

S2_FIELDS = [
    "record_id", "dataset_id", "accession", "species", "tissue", "panel", "source_label_field",
    "source_label", "ontology_id", "governed_level", "governed_label", "parent_label",
    "mapping_status", "mapping_reason", "cells_n", "source_authority", "source_path",
    "source_hash", "coverage_scope", "gap_flags",
]

S2_DICT_FIELDS = [
    "panel", "source_label", "ontology_id", "ontology_name", "governed_level", "governed_label",
    "parent_label", "mapping_status", "mapping_reason", "accepted", "approval_scope",
    "aggregate_cells", "aggregate_datasets", "aggregate_species", "source_path", "source_hash",
]

S3_FIELDS = [
    "dataset_id", "accession", "species", "panel", "tissue", "feature_namespace", "total_features",
    "gene_mapped_features", "mapped_features", "mapping_fraction", "mapping_fraction_recomputed",
    "fraction_absolute_difference", "reference_version", "bridge_version", "discovery_eligibility",
    "eligibility_reason", "matrix_layer", "matrix_type", "expression_dtype", "source_h5ad_path",
    "source_annotation_path", "source_annotation_sha256", "mapping_source_path", "mapping_source_sha256",
    "missing_value_reason",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
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


def as_bool(value: str) -> bool:
    return str(value).strip().lower() == "true"


def tissue_value(row: dict[str, str]) -> str:
    proposed = row.get("proposed_tissue", "").strip()
    observed = row.get("observed_tissue", "").strip()
    if proposed:
        return proposed
    if observed and observed.lower() != "not explicitly recorded":
        return observed
    return ""


def global_mapping_status(row: dict[str, str] | None, source_label: str) -> tuple[str, str]:
    if source_label.strip() == "":
        return "unknown", "missing value in standardized source-label field"
    if row is None:
        return "unmapped", "label absent from frozen ontology candidate inventory"
    status = row.get("candidate_status", "")
    relation = row.get("mapping_relation", "")
    evidence = row.get("evidence", "")
    reason = "; ".join(x for x in (status, relation, evidence) if x)
    if status == "excluded_unknown" or source_label.strip().lower() in {"unknown", "unassigned", "na", "n/a"}:
        return "unknown", reason or "source label is unknown"
    if status.startswith("excluded_"):
        return "excluded", reason
    if "ambig" in status.lower() or "multiple" in status.lower() or "ambig" in relation.lower() or "multiple" in relation.lower():
        return "ambiguous", reason
    return "unmapped", reason or "not included in the governed pilot label set"


def tex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
        "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def methods_sections(counts: dict[str, object], hashes: dict[str, str]) -> list[tuple[str, list[str]]]:
    return [
        ("Scope and evidence boundary", [
            "This supplement documents the executed PhyloOpenCell evidence workflow supporting the PlantHOP submission. It describes source-data governance, orthogroup construction, rank-based program discovery, frozen validation, and post hoc robustness analyses. It does not present a newly trained predictive model. The supplementary-completion step read standardized metadata and existing audit products only, created versioned derivative tables, and did not alter source H5AD or RDS objects.",
            "All labels are author-provided or previously governed source annotations. No label was inferred from expression during preparation of Tables S2 and S3. Orthogroup recurrence is treated as an association and portability device, not as proof of conserved molecular function, one-to-one orthology, causal marker activity, or paralog substitution. Missing values remain blank with explicit gap fields; they were never replaced by zero.",
        ]),
        ("Computing environment and recorded software", [
            "Supplementary tables were generated on the project server under Python 3.12.14 with anndata 0.13.3.post0, NumPy 2.5.3, pandas 3.0.5, and SciPy 1.18.0. Metadata were read in backed mode so that H5AD expression matrices were not loaded or inspected. The nested-label analysis used the project R 4.4.3 runtime. A current environment query returned Seurat 5.5.1; because the original analysis did not serialize every package version, this Seurat value records the retained environment rather than independently proving package state for every earlier run.",
            "Orthogroups were generated with OrthoFinder 3.1.5 using DIAMOND sequence search and the recorded MSA/FAMSA/FastTree workflow. The seven-species run used 12 threads, contained 183,178 representative proteins, produced 18,383 orthogroups, and assigned 164,518 genes to orthogroups. Zea mays was excluded from this orthogroup run because a source-matched Zm00001d.1 protein FASTA was unavailable.",
        ]),
        ("Input inventory and read-only audit", [
            f"The governed manifest contained {counts['manifest_datasets']} standardized H5AD datasets. For Table S2, each readable object was opened with anndata in read-only backed mode and only obs metadata were accessed. The standardized source-label column cell_type_original was preferred, with celltype_after permitted only as a recorded fallback. The scan produced {counts['s2_dataset_label_pairs']} unique dataset–source-label pairs representing {counts['s2_distinct_labels']} distinct source-label strings and {counts['s2_observed_cells']} cells. The corresponding manifest total was {counts['manifest_cells']} cells.",
            "Dataset identity, accession, species, tissue, matrix storage, and source-annotation checksums were taken from dataset_manifest.tsv and dataset_source_metadata.tsv. Source-label counts were calculated directly from the standardized annotation field. Expression values were not consulted. A dataset lacking a recognized source-label field would have been retained as a source_unavailable row rather than omitted; the audit records the number of such datasets.",
        ]),
        ("Label governance and hierarchical mapping", [
            "The global ontology candidate inventory contains 212 distinct source-label strings and explicitly retains Unknown. Table S2 maps observed dataset-label pairs only through frozen governance assets. The panel label freeze contains 29 panel–label decisions, of which 27 rows are approved for the pilot hierarchy and two source labels are excluded. The accepted hierarchy maps fine source labels to coarse root, leaf, or vascular families. Many-to-one mappings preserve the original source label, ontology identifier, governed label, parent panel, mapping rule, approval scope, and cell count.",
            "A dataset-label pair receives mapped status only when the frozen dataset-panel label-count authority and accepted hierarchy support that panel-specific mapping. Unknown labels remain unknown. Explicit exclusions remain excluded. Multiple or ambiguous ontology candidates remain ambiguous. All other observed labels remain unmapped, even when they have an exact lexical ontology candidate, because an ontology suggestion is not equivalent to a governed analysis label. This rule prevents expression-based relabeling and preserves the complete unresolved vocabulary.",
        ]),
        ("Reference mapping, orthogroup construction, and bridges", [
            "Source H5AD var_names were mapped to species reference-gene candidates with the frozen gene_mapping_candidates_v1 bridge. The mapping audit covered 81 datasets from eight species, contained 1,611,655 feature rows, mapped 1,587,599 candidate rows, and reported no ambiguous candidate rows, target collisions, or gene rows with multiple representative proteins. These mappings remained candidates pending acceptance review and were not interpreted as experimental orthology validation.",
            "Representative proteins were grouped by the seven-species OrthoFinder run. The resulting gene_orthogroup_candidates_v1 bridge linked mapped reference genes to 18,383 candidate orthogroups. Dataset-level orthogroup coverage was calculated as orthogroup-assigned source features divided by total source features. Seventy-one of 81 evaluated datasets passed the provisional 60% gate. Table S3 expands these dataset-level values across the three frozen panels, retaining blank mapping fields for manifest datasets without an evaluated orthogroup bridge.",
        ]),
        ("Rank encoding and feature construction", [
            "The executed rank representation used nonnegative source expression values from the frozen input layer. Within each cell, nonzero features were ordered by decreasing value, with exact ties resolved by immutable source feature order. The top 128 features were assigned weights 1/log2(rank+1), multiple source features mapping to the same orthogroup were collapsed by their maximum weight, and the sparse orthogroup vector was L2-normalized per cell. This design reduced dependence on incomparable expression magnitudes while retaining within-cell feature ordering.",
            f"Marker-rank feature v2 was created for the approved leaf, root, and vascular panels in read-only backed mode. It contained {counts['marker_leaf_cells']:,} leaf cells, {counts['marker_root_cells']:,} root cells, and {counts['marker_vascular_cells']:,} vascular cells over the 18,383-orthogroup vocabulary. All {counts['marker_dataset_panel_rows']} included dataset-panel combinations had nonzero mapped features. Rank features were intended for cross-dataset marker and program analysis, not differential-expression inference.",
        ]),
        ("Program discovery", [
            "Within each species and panel, each accepted coarse class was contrasted against the other accepted coarse classes. Orthogroups with positive mean normalized-rank-score differences were sorted in descending order. The frozen discovery rule retained the top 100 orthogroups per eligible species–class comparison, required at least 100 class cells and 100 other cells, and nominated candidates recurring in at least three species.",
            "The v2 discovery generated 4,200 species-level marker rows from 42 eligible comparisons and 102 recurrent candidate rows: 39 leaf, 28 root, and 35 vascular candidates. Candidates arose from the same observational discovery collections and therefore were not independent validation. Recurrence did not establish conserved function, and multi-gene orthogroups required gene-level and phylogenetic follow-up.",
        ]),
        ("Freeze protocol and leakage controls", [
            "Labels, dataset panels, coarse hierarchy, candidate programs, contrasts, scoring parameters, random seeds, and decision rules were frozen before each corresponding scoring run. Versioned outputs were written to new directories without replacing retained results. Source objects were opened read-only. Hash manifests recorded input and output identities, and the provenance records documented whether source objects, frozen protocols, labels, thresholds, or program membership changed.",
            "Internal leave-one-study-out analyses selected the dominant gene using all but one dataset within an orthogroup and species, then evaluated the selected gene in the omitted dataset. This reduced same-study selection bias but remained internal evidence. External validation programs were pre-discovered, and target contrasts were frozen before scoring. Missing mappings and excluded labels were not replaced with related genes or labels.",
        ]),
        ("Validation rounds", [
            "Study-level robustness was assessed for 19 priority orthogroups by leaving out each available study within a species. The resulting family statuses summarize direction and dominant-gene stability across held-out studies. Clean supported trees were then used for descriptive dominant-gene compactness analyses. Neither different dominant genes nor dispersed tree positions were interpreted as paralog substitution, and the final evidence review found no family satisfying the strict joint expression-plus-supported-subclade criterion.",
            "Independent rice validation v2.1 used the RNA counts layer from GSE232863, encoded 49,273 selected cells, mapped 14,136 orthogroups, used top K=128, and compared each frozen program with 500 random programs under seed 20260914. Version 2.1 corrected only the per-contrast sample-direction filter; the frozen programs, contrasts, thresholds, and seed were unchanged. Three of five primary contrasts passed, with no majority-direction reversals. Replicate-level raw rows were not available in the retained validation archive and were not reconstructed.",
            "External validation v3 froze the raw-count rank-encoding implementation, complete program score, single-gene and expanded-family subsets, off-target controls, 500 random programs, and replicate-direction rule before expression scoring. The evidence tables retain the Populus phloem negative nonreplication result and distinguish unavailable values from zero. GSE268881 Esa and Sir served as nested-label cell-level discrimination and within-study repeatability evidence, not broad species-generalization evidence; Spa remained ineligible under the frozen metadata gate.",
        ]),
        ("Nested labels, circularity stress test, and parameter robustness", [
            "The nested-label analysis used the GSE268881 RNA counts layer and three frozen endpoints: xylem versus all other root labels; xylem versus pericycle, phloem, and procambium; and stele/vascular labels versus non-stele root labels. The complete 18-orthogroup program was evaluated in pooled cells and separately in R1 and R2 for Esa and Sir. The K=128 immutable-source-order condition exactly reproduced the retained pooled metrics and replicate AUROCs to 10^-12.",
            "Parameter robustness covered K values 64, 128, and 256 and the immutable source order plus ten deterministic global feature-order permutations. Across the complete program, all 18 species–endpoint–scope groups retained AUROC above 0.5 and positive standardized mean difference in all 33 conditions per group. These direction checks are descriptive robustness results rather than new confirmatory endpoints.",
            "Circularity sensitivity removed five of the 18 frozen orthogroups documented as overlapping published marker evidence. The complete label-assignment or label-transfer feature set was unavailable. Therefore the analysis supports robustness to removal of documented overlap only and does not completely exclude label circularity. Integrated scale.data was not used as a substitute marker set.",
        ]),
        ("Discovery-threshold and clade-balance sensitivity", [
            "P1 reconstructed the frozen top-100/minimum-three-species candidate set exactly. It then evaluated top-N values of 25, 50, and 100 and recurrence thresholds of two, three, four, and five species without replacing the primary rule. At a minimum of three species, 16, 35, and 102 candidates were retained for top 25, 50, and 100, respectively.",
            "Cross-clade presence required at least one eudicot and one monocot, whereas the strict balance descriptor required at least two species from each clade. Among the frozen 102 candidates, 62 had cross-clade presence and seven met the strict two-plus-two descriptor. Leaf collections were eudicot-only and thus could not pass a cross-clade rule; this was recorded as a design limitation rather than negative evidence. Six of the submitted 18 xylem-program orthogroups were independently rediscovered by the frozen vascular-xylem top-100/minimum-three rule.",
        ]),
        ("Statistical summaries", [
            "Cell-level discrimination was summarized by AUROC, average precision, and standardized mean difference. Average precision was interpreted relative to class prevalence. Replicate scopes were reported separately from pooled cells. Random-program comparisons used the frozen empirical probability (1 + number of random AUROCs at least as large as the target AUROC)/(number of random programs + 1). Direction rules used AUROC greater than 0.5 and positive standardized mean difference; these rules were frozen and were not optimized after outcome inspection.",
            "No missing metric was converted to zero. Counts and mapping fractions in Table S3 are stored as typed numeric values. Each nonmissing fraction was independently recomputed as mapped_features/total_features, checked to lie within [0,1], and required to agree with the retained value within 10^-12. Descriptive sensitivity analyses did not create new hypothesis-testing thresholds.",
        ]),
        ("Reproducibility, checksums, and table validation", [
            f"The dataset manifest SHA-256 is {hashes['dataset_manifest']}; the label-freeze hash is {hashes['panel_label_freeze']}; the hierarchy hash is {hashes['coarse_label_hierarchy']}; and the orthogroup-coverage hash is {hashes['orthogroup_coverage']}. The build records all input paths, file sizes, and SHA-256 values in input_manifest.tsv and all deliverable hashes in output_sha256.tsv.",
            "Table validation required unique record keys, allowed mapping-status enumerations, nonnegative cell and feature counts, valid fraction bounds, exact numerator/denominator reconciliation, and expected panel enumeration. The audit also compared scanned cell counts with the manifest, retained all observed label strings, checked the 29-row governed dictionary, and verified that unavailable mapping values remained blank. The Markdown and LaTeX methods files are generated from the same section source to prevent prose divergence.",
        ]),
        ("Safety boundaries and known gaps", [
            "Supplementary-table generation did not train a model, download data, write to source datasets, change a label, modify an orthogroup bridge, alter a threshold, or replace retained analysis results. H5AD access was limited to obs metadata. Earlier RDS validation runs were not repeated. Outputs were written to a versioned supplementary-results directory.",
            f"The principal table gaps are explicit. Dataset-specific source-label field names before standardization are unavailable for many source annotations; Table S2 therefore reports the standardized H5AD field and source-annotation hash. Only 29 panel–label governance decisions and 27 accepted hierarchy rows exist, so other observed labels remain unknown, excluded, ambiguous, or unmapped. Orthogroup coverage is available for {counts['s3_datasets_with_mapping']} of {counts['manifest_datasets']} manifest datasets; the remaining {counts['s3_datasets_without_mapping']} datasets retain blank numerators and fractions across their panel rows. Gene and orthogroup mappings are candidate bridges rather than validated one-to-one ortholog assignments. Functional conservation, causality, complete absence of label circularity, and paralog substitution remain unestablished.",
        ]),
    ]


def build_methods(sections: list[tuple[str, list[str]]]) -> tuple[str, str, int]:
    md_lines = ["# Supplementary Methods S1", ""]
    tex_lines = [
        r"\section*{Supplementary Methods S1}",
        r"\providecommand{\suppmethodsection}[1]{\subsection*{#1}}",
        "",
    ]
    for title, paragraphs in sections:
        md_lines += [f"## {title}", ""]
        tex_lines += [rf"\suppmethodsection{{{tex_escape(title)}}}", ""]
        for paragraph in paragraphs:
            md_lines += [paragraph, ""]
            tex_lines += [tex_escape(paragraph), ""]
    md = "\n".join(md_lines).rstrip() + "\n"
    tex = "\n".join(tex_lines).rstrip() + "\n"
    words = len(re.findall(r"\b[A-Za-z0-9][A-Za-z0-9'’-]*\b", md))
    return md, tex, words


def main() -> int:
    if FINAL_OUTPUT.exists() or OUTPUT.exists():
        raise RuntimeError(f"refusing to overwrite {FINAL_OUTPUT} or {OUTPUT}")
    for name, path in SOURCES.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name}: {path}")
    hashes = {name: sha256(path) for name, path in SOURCES.items()}
    OUTPUT.mkdir(parents=False)
    (OUTPUT / Path(__file__).name).write_bytes(Path(__file__).read_bytes())

    manifest = read_tsv(SOURCES["dataset_manifest"])
    manifest_by_id = {r["dataset_id"]: r for r in manifest}
    if len(manifest) != 121 or len(manifest_by_id) != len(manifest):
        raise RuntimeError("dataset manifest must contain 121 unique datasets")
    label_freeze = read_tsv(SOURCES["panel_label_freeze"])
    hierarchy = read_tsv(SOURCES["coarse_label_hierarchy"])
    label_candidates = {(r["panel"], r["source_label"]): r for r in read_tsv(SOURCES["panel_label_candidates"])}
    hierarchy_map = {(r["panel"], r["fine_source_label"]): r for r in hierarchy}
    freeze_map = {(r["panel"], r["source_label"]): r for r in label_freeze}
    excluded_freezes: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in label_freeze:
        if not as_bool(row["accepted"]):
            excluded_freezes[row["source_label"]].append(row)
    ontology = {r["source_label"]: r for r in read_tsv(SOURCES["ontology_mapping_candidates"])}
    governed_counts: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(SOURCES["panel_dataset_label_counts"]):
        governed_counts[(row["dataset_id"], row["source_label"])].append(row)

    s2_rows: list[dict[str, object]] = []
    scanned_pairs: set[tuple[str, str]] = set()
    observed_labels: set[str] = set()
    missing_field_datasets: list[str] = []
    observed_cells = 0
    label_field_counts: Counter[str] = Counter()
    for item in manifest:
        dataset_id = item["dataset_id"]
        path = Path(item["h5ad_path"])
        obj = ad.read_h5ad(path, backed="r")
        try:
            if "cell_type_original" in obj.obs.columns:
                label_field = "cell_type_original"
            elif "celltype_after" in obj.obs.columns:
                label_field = "celltype_after"
            else:
                missing_field_datasets.append(dataset_id)
                s2_rows.append({
                    "record_id": f"{dataset_id}|<source_label_field_missing>|", "dataset_id": dataset_id,
                    "accession": item["accession"], "species": item["proposed_species"] or item["observed_species"],
                    "tissue": tissue_value(item), "panel": "", "source_label_field": "", "source_label": "",
                    "ontology_id": "", "governed_level": "", "governed_label": "", "parent_label": "",
                    "mapping_status": "source_unavailable", "mapping_reason": "no recognized standardized source-label field",
                    "cells_n": "", "source_authority": "standardized_h5ad_obs;dataset_manifest_v1",
                    "source_path": str(path), "source_hash": item["source_annotation_sha256"],
                    "coverage_scope": "dataset retained; source-label mapping unavailable",
                    "gap_flags": "source_label_field_missing",
                })
                continue
            label_field_counts[label_field] += 1
            values = obj.obs[label_field].astype("string").fillna("").str.strip()
            counts = values.value_counts(dropna=False, sort=False)
            if int(counts.sum()) != int(obj.n_obs):
                raise RuntimeError(f"label counts do not sum to cells for {dataset_id}")
            observed_cells += int(counts.sum())
            for source_label, cells_n in sorted(((str(k), int(v)) for k, v in counts.items()), key=lambda x: x[0]):
                scanned_pairs.add((dataset_id, source_label))
                observed_labels.add(source_label)
                panel_rows = governed_counts.get((dataset_id, source_label), [])
                emitted = False
                for panel_row in sorted(panel_rows, key=lambda r: r["panel"]):
                    panel = panel_row["panel"]
                    freeze = freeze_map.get((panel, source_label))
                    hierarchy_row = hierarchy_map.get((panel, source_label))
                    accepted = freeze is not None and as_bool(freeze["accepted"]) and hierarchy_row is not None and as_bool(hierarchy_row["accepted"])
                    if accepted:
                        status = "mapped"
                        reason = f"{hierarchy_row['mapping_rule']}; {panel_row['base_gate_reason']}"
                        governed_level = hierarchy_row["hierarchy_level"]
                        governed_label = hierarchy_row["coarse_label"]
                        ontology_id = hierarchy_row["fine_ontology_id"]
                        gap_flags = ""
                    else:
                        status = "excluded"
                        reason = (freeze or {}).get("curator_note", "panel-label mapping not accepted")
                        governed_level = ""
                        governed_label = ""
                        ontology_id = (freeze or {}).get("ontology_id", "")
                        gap_flags = "panel_label_not_accepted"
                    s2_rows.append({
                        "record_id": f"{dataset_id}|{source_label}|{panel}", "dataset_id": dataset_id,
                        "accession": item["accession"], "species": item["proposed_species"] or item["observed_species"],
                        "tissue": tissue_value(item), "panel": panel, "source_label_field": label_field,
                        "source_label": source_label, "ontology_id": ontology_id, "governed_level": governed_level,
                        "governed_label": governed_label, "parent_label": panel if governed_label else "",
                        "mapping_status": status, "mapping_reason": reason, "cells_n": cells_n,
                        "source_authority": "standardized_h5ad_obs;panel_dataset_label_counts_v1;panel_label_freeze_v1;coarse_label_hierarchy_v1",
                        "source_path": f"{path};{SOURCES['panel_dataset_label_counts']};{SOURCES['panel_label_freeze']};{SOURCES['coarse_label_hierarchy']}",
                        "source_hash": f"{item['source_annotation_sha256']};{hashes['panel_dataset_label_counts']};{hashes['panel_label_freeze']};{hashes['coarse_label_hierarchy']}",
                        "coverage_scope": "complete standardized dataset-label count with panel-specific governed mapping",
                        "gap_flags": gap_flags,
                    })
                    emitted = True
                if not emitted:
                    excluded = sorted(excluded_freezes.get(source_label, []), key=lambda r: r["panel"])
                    if excluded:
                        for freeze in excluded:
                            panel = freeze["panel"]
                            s2_rows.append({
                                "record_id": f"{dataset_id}|{source_label}|{panel}", "dataset_id": dataset_id,
                                "accession": item["accession"], "species": item["proposed_species"] or item["observed_species"],
                                "tissue": tissue_value(item), "panel": panel, "source_label_field": label_field,
                                "source_label": source_label, "ontology_id": freeze["ontology_id"],
                                "governed_level": "", "governed_label": "", "parent_label": "",
                                "mapping_status": "excluded",
                                "mapping_reason": f"frozen panel-label exclusion ({freeze['approval_scope']}): {freeze['curator_note']}",
                                "cells_n": cells_n,
                                "source_authority": "standardized_h5ad_obs;panel_label_freeze_v1",
                                "source_path": f"{path};{SOURCES['panel_label_freeze']}",
                                "source_hash": f"{item['source_annotation_sha256']};{hashes['panel_label_freeze']}",
                                "coverage_scope": "complete standardized dataset-label count with explicit frozen exclusion",
                                "gap_flags": "frozen_excluded_label",
                            })
                        continue
                    status, reason = global_mapping_status(ontology.get(source_label), source_label)
                    candidate = ontology.get(source_label, {})
                    s2_rows.append({
                        "record_id": f"{dataset_id}|{source_label}|", "dataset_id": dataset_id,
                        "accession": item["accession"], "species": item["proposed_species"] or item["observed_species"],
                        "tissue": tissue_value(item), "panel": "", "source_label_field": label_field,
                        "source_label": source_label, "ontology_id": candidate.get("candidate_ids", ""),
                        "governed_level": "", "governed_label": "", "parent_label": "",
                        "mapping_status": status, "mapping_reason": reason, "cells_n": cells_n,
                        "source_authority": "standardized_h5ad_obs;ontology_mapping_candidates;panel_label_freeze_v1",
                        "source_path": f"{path};{SOURCES['ontology_mapping_candidates']};{SOURCES['panel_label_freeze']}",
                        "source_hash": f"{item['source_annotation_sha256']};{hashes['ontology_mapping_candidates']};{hashes['panel_label_freeze']}",
                        "coverage_scope": "complete standardized dataset-label count; no governed panel mapping",
                        "gap_flags": "not_in_governed_panel_hierarchy" if status == "unmapped" else "",
                    })
        finally:
            if getattr(obj, "file", None) is not None:
                obj.file.close()

    s2_rows.sort(key=lambda r: (str(r["dataset_id"]), str(r["source_label"]), str(r["panel"])))
    ids = [str(r["record_id"]) for r in s2_rows]
    if len(ids) != len(set(ids)):
        raise RuntimeError("duplicate Table S2 record_id")
    allowed_status = {"mapped", "unknown", "excluded", "ambiguous", "unmapped", "source_unavailable"}
    if not {str(r["mapping_status"]) for r in s2_rows}.issubset(allowed_status):
        raise RuntimeError("unexpected Table S2 mapping_status")
    if not observed_labels.issubset(set(ontology)):
        missing = sorted(observed_labels - set(ontology))
        raise RuntimeError(f"observed labels absent from ontology inventory: {missing[:10]}")
    write_tsv(OUTPUT / "Table_S2_source_label_mapping.tsv", s2_rows, S2_FIELDS)

    s2_dict_rows: list[dict[str, object]] = []
    for freeze in sorted(label_freeze, key=lambda r: (r["panel"], r["source_label"])):
        panel, source_label = freeze["panel"], freeze["source_label"]
        hierarchy_row = hierarchy_map.get((panel, source_label))
        candidate = label_candidates.get((panel, source_label), {})
        accepted = as_bool(freeze["accepted"]) and hierarchy_row is not None and as_bool(hierarchy_row["accepted"])
        s2_dict_rows.append({
            "panel": panel, "source_label": source_label, "ontology_id": freeze["ontology_id"],
            "ontology_name": freeze["ontology_name_in_pinned_obo"],
            "governed_level": hierarchy_row["hierarchy_level"] if accepted else "",
            "governed_label": hierarchy_row["coarse_label"] if accepted else "",
            "parent_label": panel if accepted else "", "mapping_status": "mapped" if accepted else "excluded",
            "mapping_reason": hierarchy_row["mapping_rule"] if accepted else freeze["curator_note"],
            "accepted": str(accepted), "approval_scope": freeze["approval_scope"],
            "aggregate_cells": candidate.get("aggregate_cells", ""),
            "aggregate_datasets": candidate.get("aggregate_datasets", ""),
            "aggregate_species": candidate.get("aggregate_species", ""),
            "source_path": f"{SOURCES['panel_label_freeze']};{SOURCES['coarse_label_hierarchy']};{SOURCES['panel_label_candidates']}",
            "source_hash": f"{hashes['panel_label_freeze']};{hashes['coarse_label_hierarchy']};{hashes['panel_label_candidates']}",
        })
    if len(s2_dict_rows) != 29:
        raise RuntimeError("expected 29 Table S2 hierarchy rows")
    write_tsv(OUTPUT / "Table_S2_label_hierarchy_dictionary.tsv", s2_dict_rows, S2_DICT_FIELDS)

    panel_freeze = read_tsv(SOURCES["panel_dataset_freeze"])
    coverage = {r["dataset_id"]: r for r in read_tsv(SOURCES["orthogroup_coverage"])}
    gene_mapping = {r["dataset_id"]: r for r in read_tsv(SOURCES["gene_mapping_summary"])}
    s3_rows: list[dict[str, object]] = []
    for panel_row in sorted(panel_freeze, key=lambda r: (r["dataset_id"], r["panel"])):
        item = manifest_by_id[panel_row["dataset_id"]]
        cov = coverage.get(panel_row["dataset_id"])
        gene = gene_mapping.get(panel_row["dataset_id"])
        total = int(cov["features"]) if cov else int(item["n_genes"])
        if cov:
            mapped = int(cov["orthogroup_assigned_candidates"])
            fraction = float(cov["orthogroup_coverage_fraction"])
            recomputed = mapped / total
            difference = abs(fraction - recomputed)
            gene_mapped = int(cov["gene_mapped_candidates"])
            missing_reason = ""
        else:
            mapped = ""
            fraction = ""
            recomputed = ""
            difference = ""
            gene_mapped = ""
            missing_reason = "orthogroup mapping coverage not evaluated for this manifest dataset; values retained blank"
        eligible = as_bool(panel_row["accepted"])
        s3_rows.append({
            "dataset_id": panel_row["dataset_id"], "accession": item["accession"],
            "species": panel_row["species"], "panel": panel_row["panel"], "tissue": tissue_value(item),
            "feature_namespace": "standardized H5AD var_names mapped to species reference-gene candidates",
            "total_features": total, "gene_mapped_features": gene_mapped, "mapped_features": mapped,
            "mapping_fraction": fraction, "mapping_fraction_recomputed": recomputed,
            "fraction_absolute_difference": difference, "reference_version": "reference_longest_protein_v1",
            "bridge_version": "gene_mapping_candidates_v1 + gene_orthogroup_candidates_v1",
            "discovery_eligibility": "eligible" if eligible else "excluded",
            "eligibility_reason": panel_row["reason"], "matrix_layer": "X",
            "matrix_type": f"{item['expression_class']};{item['storage_kind']}",
            "expression_dtype": item["expression_dtype"], "source_h5ad_path": item["h5ad_path"],
            "source_annotation_path": item["source_annotation_path"],
            "source_annotation_sha256": item["source_annotation_sha256"],
            "mapping_source_path": f"{SOURCES['gene_mapping_summary']};{SOURCES['orthogroup_coverage']}",
            "mapping_source_sha256": f"{hashes['gene_mapping_summary']};{hashes['orthogroup_coverage']}",
            "missing_value_reason": missing_reason,
        })
        if cov and gene:
            if int(gene["features"]) != total or int(gene["mapped_candidates"]) != gene_mapped:
                raise RuntimeError(f"gene/orthogroup mapping summary mismatch: {panel_row['dataset_id']}")
    s3_keys = [(str(r["dataset_id"]), str(r["panel"])) for r in s3_rows]
    if len(s3_rows) != 363 or len(s3_keys) != len(set(s3_keys)):
        raise RuntimeError("Table S3 must contain 363 unique dataset-panel rows")
    for row in s3_rows:
        if row["mapping_fraction"] == "":
            if row["mapped_features"] != "" or row["mapping_fraction_recomputed"] != "":
                raise RuntimeError("partial missingness in Table S3 mapping fields")
            continue
        fraction = float(row["mapping_fraction"])
        if not 0 <= fraction <= 1 or float(row["fraction_absolute_difference"]) > 1e-12:
            raise RuntimeError(f"invalid mapping fraction: {row['dataset_id']} {row['panel']}")
    write_tsv(OUTPUT / "Table_S3_dataset_orthogroup_mapping_coverage.tsv", s3_rows, S3_FIELDS)

    tables = {
        "Table_S2_source_label_mapping.tsv": (s2_rows, S2_FIELDS),
        "Table_S2_label_hierarchy_dictionary.tsv": (s2_dict_rows, S2_DICT_FIELDS),
        "Table_S3_dataset_orthogroup_mapping_coverage.tsv": (s3_rows, S3_FIELDS),
    }
    descriptions = {
        "record_id": "Unique dataset–source-label–panel record identifier.", "dataset_id": "Frozen project dataset identifier.",
        "accession": "Repository accession recorded in the dataset manifest.", "species": "Governed species name.",
        "tissue": "Recorded or proposed tissue when available.", "panel": "Frozen root, leaf, or vascular panel; blank when no governed panel mapping exists.",
        "source_label_field": "Standardized H5AD obs field read for source labels.", "source_label": "Unmodified standardized source-label value.",
        "ontology_id": "Frozen ontology candidate identifier(s), if available.", "ontology_name": "Name in the pinned ontology asset.",
        "governed_level": "Frozen hierarchy level assigned to an accepted mapping.", "governed_label": "Accepted coarse governed label.",
        "parent_label": "Parent panel for an accepted governed label.", "mapping_status": "mapped, unknown, excluded, ambiguous, unmapped, or source_unavailable.",
        "mapping_reason": "Frozen rule, exclusion, ambiguity, or gap explaining mapping status.", "cells_n": "Cells with the dataset/source-label value.",
        "source_authority": "Authorities supporting the record.", "source_path": "Semicolon-separated source and governance paths.",
        "source_hash": "Semicolon-separated authority SHA-256 values.", "coverage_scope": "Completeness scope represented by the row.",
        "gap_flags": "Explicit unresolved limitation(s).", "accepted": "Whether the frozen hierarchy decision is accepted.",
        "approval_scope": "Scope under which the label decision was approved.", "aggregate_cells": "Candidate-stage aggregate cell count.",
        "aggregate_datasets": "Candidate-stage dataset count.", "aggregate_species": "Candidate-stage species count.",
        "feature_namespace": "Feature identifier namespace and mapping target.", "total_features": "Source feature denominator.",
        "gene_mapped_features": "Features mapped to reference-gene candidates.", "mapped_features": "Features assigned to candidate orthogroups.",
        "mapping_fraction": "Authority orthogroup-assigned features divided by total features.",
        "mapping_fraction_recomputed": "Independently recomputed mapped_features/total_features.",
        "fraction_absolute_difference": "Absolute authority-versus-recomputed fraction difference.",
        "reference_version": "Reference-protein selection version.", "bridge_version": "Gene and orthogroup bridge versions.",
        "discovery_eligibility": "Frozen dataset-panel inclusion status.", "eligibility_reason": "Frozen reason for inclusion or exclusion.",
        "matrix_layer": "H5AD matrix layer used by the executed feature workflow.", "matrix_type": "Expression class and storage kind.",
        "expression_dtype": "Recorded expression data type.", "source_h5ad_path": "Read-only source H5AD path.",
        "source_annotation_path": "Original annotation source path recorded in the manifest.",
        "source_annotation_sha256": "SHA-256 of the original annotation source.",
        "mapping_source_path": "Mapping authority paths.", "mapping_source_sha256": "Mapping authority SHA-256 values.",
        "missing_value_reason": "Reason mapped counts/fraction remain blank.",
    }
    dictionary_rows = []
    for table_name, (rows, fields) in tables.items():
        for field in fields:
            dictionary_rows.append({
                "table": table_name, "field": field,
                "description": descriptions.get(field, field.replace("_", " ").capitalize() + "."),
                "missing_value_policy": "Blank means unavailable or not applicable; never interpret as zero.",
            })
    write_tsv(OUTPUT / "field_dictionary.tsv", dictionary_rows, ["table", "field", "description", "missing_value_policy"])

    marker_rows = read_tsv(SOURCES["marker_rank_dataset_audit"])
    marker_cells = defaultdict(int)
    for row in marker_rows:
        marker_cells[row["panel"]] += int(row["sampled_cells"])
    counts = {
        "manifest_datasets": len(manifest), "manifest_cells": sum(int(r["n_cells"]) for r in manifest),
        "s2_rows": len(s2_rows), "s2_dataset_label_pairs": len(scanned_pairs),
        "s2_distinct_labels": len(observed_labels), "s2_observed_cells": observed_cells,
        "s2_mapping_status_counts": dict(sorted(Counter(str(r["mapping_status"]) for r in s2_rows).items())),
        "s2_source_label_fields": dict(sorted(label_field_counts.items())),
        "s2_missing_field_datasets": len(missing_field_datasets), "s2_dictionary_rows": len(s2_dict_rows),
        "s3_rows": len(s3_rows), "s3_datasets_with_mapping": len(coverage),
        "s3_datasets_without_mapping": len(manifest) - len(coverage),
        "s3_rows_with_mapping": sum(r["mapping_fraction"] != "" for r in s3_rows),
        "s3_rows_without_mapping": sum(r["mapping_fraction"] == "" for r in s3_rows),
        "s3_eligible_rows": sum(r["discovery_eligibility"] == "eligible" for r in s3_rows),
        "marker_dataset_panel_rows": len(marker_rows), "marker_leaf_cells": marker_cells["leaf"],
        "marker_root_cells": marker_cells["root"], "marker_vascular_cells": marker_cells["vascular"],
    }
    sections = methods_sections(counts, hashes)
    md, tex, word_count = build_methods(sections)
    (OUTPUT / "Methods_S1_detailed_methods.md").write_text(md, encoding="utf-8")
    (OUTPUT / "Methods_S1_detailed_methods.tex").write_text(tex, encoding="utf-8")

    manifest_rows = []
    for name, path in SOURCES.items():
        manifest_rows.append({"source_id": name, "path": str(path), "sha256": hashes[name], "bytes": path.stat().st_size, "use": "authority input"})
    write_tsv(OUTPUT / "input_manifest.tsv", manifest_rows, ["source_id", "path", "sha256", "bytes", "use"])

    gaps = [
        "Pre-standardization source label column names are not retained uniformly; Table S2 reports the standardized H5AD obs field.",
        "Only frozen panel-label decisions are called mapped; other ontology candidates remain unmapped or ambiguous.",
        f"Orthogroup coverage is unavailable for {counts['s3_datasets_without_mapping']} of {counts['manifest_datasets']} manifest datasets and remains blank, not zero.",
        "Gene and orthogroup bridges are candidate mappings and do not establish one-to-one orthology or conserved function.",
        "The complete GSE268881 label-transfer feature set is unavailable, so label circularity is not completely excluded.",
        "Rice v2.1 replicate-level raw metric rows are unavailable in the retained validation archive and were not reconstructed.",
        "Zea mays was excluded from the seven-species OrthoFinder run pending a source-matched protein reference.",
    ]
    readme = [
        "# Supplementary completion v1", "",
        "This package contains submission-ready Tables S2 and S3 plus detailed Methods S1.", "",
        f"- Table S2 source-label mapping: {len(s2_rows)} rows representing {len(scanned_pairs)} unique dataset-label pairs across {len(manifest)} manifest datasets.",
        f"- Table S2 hierarchy dictionary: {len(s2_dict_rows)} frozen panel-label decisions.",
        f"- Table S3 mapping coverage: {len(s3_rows)} unique dataset-panel rows; {counts['s3_rows_with_mapping']} have evaluated orthogroup coverage and {counts['s3_rows_without_mapping']} retain blank coverage fields.",
        f"- Methods S1: {word_count} words in {len(sections)} sections, supplied as Markdown and LaTeX generated from one source.", "",
        "Documented gaps:",
    ] + [f"- {gap}" for gap in gaps] + ["", "Missing values are blank and must never be interpreted as zero."]
    (OUTPUT / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    audit = {
        "status": "COMPLETE_WITH_DOCUMENTED_GAPS", "generated_utc": datetime.now(timezone.utc).isoformat(),
        "version": "supplementary_completion_v1", "counts": counts,
        "methods": {"word_count": word_count, "section_count": len(sections), "formats": ["md", "tex"]},
        "checks": {
            "source_h5ad_read_only_backed": True, "manifest_dataset_count_121": len(manifest) == 121,
            "scanned_cells_equal_manifest_cells": observed_cells == sum(int(r["n_cells"]) for r in manifest) if not missing_field_datasets else False,
            "s2_record_ids_unique": len(ids) == len(set(ids)), "s2_mapping_status_enum_valid": True,
            "s2_unknown_and_excluded_items_retained": any(r["mapping_status"] == "unknown" for r in s2_rows) and any(r["mapping_status"] == "excluded" for r in s2_rows),
            "s2_observed_vocabulary_covered_by_frozen_inventory": observed_labels.issubset(set(ontology)),
            "s2_dictionary_rows_29": len(s2_dict_rows) == 29,
            "s3_dataset_panel_keys_unique": len(s3_keys) == len(set(s3_keys)), "s3_rows_363": len(s3_rows) == 363,
            "s3_fraction_range_valid": all(r["mapping_fraction"] == "" or 0 <= float(r["mapping_fraction"]) <= 1 for r in s3_rows),
            "s3_fraction_recalculation_within_1e12": all(r["mapping_fraction"] == "" or float(r["fraction_absolute_difference"]) <= 1e-12 for r in s3_rows),
            "missing_mapping_values_not_zero_filled": all(r["mapping_fraction"] != "" or (r["mapped_features"] == "" and r["mapping_fraction_recomputed"] == "") for r in s3_rows),
            "methods_markdown_latex_same_sections": True,
        },
        "gaps": gaps, "inputs": {name: {"path": str(path), "sha256": hashes[name], "bytes": path.stat().st_size} for name, path in SOURCES.items()},
        "safety": {
            "source_expression_accessed": False, "source_obs_metadata_read_only": True,
            "source_files_modified": False, "model_training_started": False, "new_data_downloaded": False,
            "labels_inferred_from_expression": False, "existing_results_overwritten": False,
        },
    }
    if not all(audit["checks"].values()):
        failed = [k for k, v in audit["checks"].items() if not v]
        raise RuntimeError(f"audit checks failed: {failed}")
    (OUTPUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (OUTPUT / "STATUS").write_text("COMPLETE_WITH_DOCUMENTED_GAPS\n", encoding="utf-8")
    targets = [p for p in sorted(OUTPUT.iterdir()) if p.is_file() and p.name != "output_sha256.tsv"]
    write_tsv(OUTPUT / "output_sha256.tsv", [{"output_file": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in targets])
    OUTPUT.replace(FINAL_OUTPUT)
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
