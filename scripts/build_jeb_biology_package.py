#!/usr/bin/env python3
"""Build a conservative JEB-facing evidence matrix and biology figure plan.

This script only summarizes existing PhyloOpenCell governance and biological
evidence tables. It does not read source H5AD files, train models, or access
external services.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("/workspace/projects/PhyloOpenCell")
META = ROOT / "metadata"
REPORTS = ROOT / "reports"
DOCS = ROOT / "docs"


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    governance = read_json(REPORTS / "governance_summary.json")
    source = read_json(REPORTS / "source_metadata_summary.json")
    freeze = read_json(REPORTS / "pilot_freeze_v1_audit.json")
    marker = read_json(REPORTS / "conserved_program_candidates_v2_audit.json")
    function = read_json(REPORTS / "conserved_program_function_evidence_v2_audit.json")
    paralog = read_json(REPORTS / "priority_paralog_evidence_v2_audit.json")
    tree = read_json(REPORTS / "priority_paralog_gene_tree_evidence_v2_audit.json")
    candidates = read_tsv(META / "conserved_program_function_evidence_v2.tsv")
    tree_rows = read_tsv(META / "priority_paralog_gene_tree_evidence_v2.tsv")
    gene_effects = read_tsv(META / "priority_paralog_gene_effects_v2.tsv")
    dataset_checks = read_tsv(META / "priority_paralog_dataset_checks_v2.tsv")

    priority = [r for r in candidates if r["priority_tier"] == "A"]
    annotated = [r for r in candidates if float(r["annotation_coverage_fraction"] or 0) > 0]
    zero_eligible = sum(int(r["eligible_orthogroup_label_tests"]) == 0 for r in dataset_checks)
    patterns = tree["patterns"]
    recurrence = {}
    for row in candidates:
        count = row["species_count"]
        recurrence[count] = recurrence.get(count, 0) + 1

    evidence_rows = [
        {
            "claim_id": "C1",
            "proposed_claim": "The governed atlas provides broad but heterogeneous plant single-cell coverage suitable for a cross-species resource.",
            "current_status": "supported_with_boundary",
            "supporting_evidence": f"{governance['readable_files']}/{governance['h5ad_files']} H5AD readable; {governance['cells']} cells; {governance['proposed_species']} proposed species; {source['datasets_with_official_source']} datasets with official source metadata.",
            "claim_boundary": f"Unknown-label fraction is {governance['unknown_fraction']:.3f}; expression representations are heterogeneous; three source species matches are complex or mismatched.",
            "primary_source": "reports/governance_summary.json; reports/source_metadata_summary.json",
            "recommended_placement": "Figure 1; Table 1; Results 1",
            "next_required_evidence": "Freeze the manuscript cohort and report dataset-level inclusion/exclusion reasons.",
        },
        {
            "claim_id": "C2",
            "proposed_claim": "A traceable hierarchical vocabulary reconciles source labels without rewriting the raw files.",
            "current_status": "supported_for_pilot_scope",
            "supporting_evidence": f"{freeze['counts']['approved_unique_labels']} unique labels, {freeze['counts']['approved_panel_label_rows']} panel-label rows, {freeze['counts']['approved_panel_dataset_rows']} panel-dataset rows, and {freeze['counts']['approved_split_rows']} split rows were frozen with hashes.",
            "claim_boundary": "The approved freeze covers a low-cost pilot subset, not every source label or every tissue.",
            "primary_source": "metadata/panel_label_freeze_v1.tsv; metadata/coarse_label_hierarchy_v1.tsv; reports/pilot_freeze_v1_audit.json",
            "recommended_placement": "Figure 2; Supplementary Tables S2-S4; Results 2",
            "next_required_evidence": "Expert review of medium-confidence collapses and an explicit unresolved-label accounting table.",
        },
        {
            "claim_id": "C3",
            "proposed_claim": "Orthogroup-rank analysis recovers recurrent coarse cell-system programs across several plant species.",
            "current_status": "supported_discovery",
            "supporting_evidence": f"{marker['conserved_candidate_rows']} recurrent candidates; {recurrence.get('3', 0)} in 3 species, {recurrence.get('4', 0)} in 4 species, and {recurrence.get('5', 0)} in 5 species; {len(priority)} Priority A candidates.",
            "claim_boundary": "Recurrence is based on frozen sampled cells and marker-ranking criteria; it is not yet independent validation.",
            "primary_source": "metadata/conserved_program_candidates_v2.tsv; reports/conserved_program_candidates_v2_audit.json",
            "recommended_placement": "Figure 3; Results 3",
            "next_required_evidence": "Repeat in an independent held-out dataset/species and report study-level uncertainty.",
        },
        {
            "claim_id": "C4",
            "proposed_claim": "Recurrent programs include biologically coherent photosynthetic, vascular, phloem, epidermal, and root-stele gene families.",
            "current_status": "supported_by_local_annotation",
            "supporting_evidence": f"{function['candidate_rows_with_any_description']}/{function['candidate_rows']} candidate orthogroups have a local sequence-header description; {len(annotated)} candidate rows have non-zero annotation coverage.",
            "claim_boundary": "Local FASTA-header annotation is incomplete and does not by itself establish conserved function or cell-type specificity.",
            "primary_source": "metadata/conserved_program_function_evidence_v2.tsv; reports/conserved_program_function_evidence_v2_audit.json",
            "recommended_placement": "Figure 3; Supplementary Table S5; Results 3",
            "next_required_evidence": "Curate literature/ontology support for the short list and test independent marker enrichment.",
        },
        {
            "claim_id": "C5",
            "proposed_claim": "Gene-level expression identifies species-specific dominant members within expanded conserved families.",
            "current_status": "supported_exploratory",
            "supporting_evidence": f"{paralog['priority_orthogroups']} priority orthogroups, {paralog['dataset_gene_effect_rows']} dataset-gene effects, {paralog['aggregated_gene_rows']} aggregated gene effects, and {paralog['species_dominance_rows']} species-dominance records; {paralog['orthogroups_ready_for_gene_tree_mapping']} orthogroups were tree-ready.",
            "claim_boundary": f"{zero_eligible} sampled dataset-panel rows had zero eligible orthogroup-label tests; dominance is not evidence of functional replacement.",
            "primary_source": "metadata/priority_paralog_gene_effects_v2.tsv; metadata/priority_paralog_species_dominance_v2.tsv; reports/priority_paralog_evidence_v2_audit.json",
            "recommended_placement": "Figure 4; Supplementary Table S6; Results 4",
            "next_required_evidence": "Require consistency across independent studies and quantify effect uncertainty for each selected family.",
        },
        {
            "claim_id": "C6",
            "proposed_claim": "AAP and LTP family drivers show exploratory conserved-subclade signals across species.",
            "current_status": "supported_exploratory",
            "supporting_evidence": f"{patterns.get('conserved_subclade_signal', 0)} orthogroups met the prespecified exploratory compactness tail: OG0000378 (AAP) and OG0000114 (LTP).",
            "claim_boundary": "The signal depends on the resolved OrthoFinder tree and random one-gene-per-species null; sequence-tree uncertainty remains.",
            "primary_source": "metadata/priority_paralog_gene_tree_evidence_v2.tsv; reports/priority_paralog_gene_tree_evidence_v2_audit.json",
            "recommended_placement": "Figure 5A-B; Results 4",
            "next_required_evidence": "Manually inspect alignments/tree support and validate expression in independent datasets before biological interpretation.",
        },
        {
            "claim_id": "C7",
            "proposed_claim": "PSBO is a candidate family for lineage-associated paralog switching.",
            "current_status": "hypothesis_only",
            "supporting_evidence": f"{patterns.get('dispersed_drivers_possible_paralog_switching', 0)} family (OG0003912) fell in the high-dispersion tail of the exploratory tree-distance null.",
            "claim_boundary": "Paralog substitution is not established; dispersion may reflect tree uncertainty, annotation error, duplication, or incomplete sampling.",
            "primary_source": "metadata/priority_paralog_gene_tree_evidence_v2.tsv; metadata/priority_paralog_gene_effects_v2.tsv",
            "recommended_placement": "Figure 5C only if independently validated; otherwise Supplementary candidate dossier",
            "next_required_evidence": "Alignment/tree-support review, orthology curation, reciprocal expression evidence, and external biological support.",
        },
    ]
    matrix_path = META / "jeb_biology_evidence_matrix_v1.tsv"
    write_tsv(matrix_path, evidence_rows, list(evidence_rows[0]))

    selected_ids = {"OG0000378", "OG0000114", "OG0003912"}
    candidate_lookup = {r["orthogroup_id"]: r for r in candidates}
    dossier_rows = []
    for tr in tree_rows:
        og = tr["orthogroup_id"]
        if og not in selected_ids:
            continue
        cr = candidate_lookup[og]
        effects = [r for r in gene_effects if r["orthogroup_id"] == og]
        positive = [r for r in effects if float(r["median_standardized_mean_difference"]) > 0]
        dossier_rows.append({
            "orthogroup_id": og,
            "panel": tr["panel"],
            "coarse_label": tr["coarse_label"],
            "local_symbols": cr["local_symbols_all_species"],
            "local_descriptions": cr["local_descriptions"],
            "species_recurrence": cr["species_count"],
            "tree_pattern": tr["gene_tree_pattern"],
            "distance_percentile": tr["distance_percentile_fraction_random_le_observed"],
            "dominant_genes_by_species": tr["dominant_genes_by_species"],
            "positive_gene_effect_rows": str(len(positive)),
            "paralog_substitution_claim": tr["paralog_substitution_claim"],
            "manuscript_decision": "retain_as_conserved_family_example" if tr["gene_tree_pattern"] == "conserved_subclade_signal" else "supplementary_hypothesis_until_independent_validation",
        })
    dossier_path = META / "jeb_priority_family_dossier_v1.tsv"
    write_tsv(dossier_path, dossier_rows, list(dossier_rows[0]))

    figure_plan = f"""# JEB biology-first figure plan (v1)

This plan is for the annotation/resource and biological-discovery strand of PhyloOpenCell. It does not claim results from the separate model benchmark.

## Proposed manuscript message

A governed cross-species plant single-cell resource and orthogroup-aware analysis reveal recurrent coarse cell-system programs, while gene-level and gene-tree checks sharply limit claims of paralog replacement.

## Main figures

### Figure 1 — Governed resource landscape

- Species/tissue/dataset flow from raw inventory to frozen analysis panels.
- Dataset-level cells, label completeness, expression representation, provenance status, and exclusions.
- Headline denominators: {governance['readable_files']} readable H5AD files, {governance['cells']:,} cells, {governance['proposed_species']} proposed species, unknown-label fraction {governance['unknown_fraction']:.1%}.
- Boundary: distinguish the full governed resource from the smaller frozen discovery cohort.

### Figure 2 — Hierarchical label reconciliation

- Sankey or alluvial map from source labels to tissue system and coarse lineage.
- Confidence and unresolved/excluded categories shown explicitly.
- Dataset-by-label coverage heat map for leaf, root, and vascular panels.
- Do not imply that pilot-approved labels cover the complete resource.

### Figure 3 — Cross-species conserved programs

- Orthogroup-by-species heat map for the {len(priority)} Priority A candidates.
- Separate panels for photosynthetic ground tissue, epidermal system, root stele, xylem, and phloem.
- Show per-species rank-score differences rather than pooled cells.
- Annotate coherent families (photosystem/carbon fixation, AAP, LTP, UXS/IRX) using local descriptions; distinguish curated names from uncharacterized groups.

### Figure 4 — From orthogroups to expressed family members

- For selected families, show study-level gene effects within each species.
- Use OG0000378 (AAP) and OG0000114 (LTP) as conserved-family examples.
- Display study consistency and missing/zero-eligible datasets; do not show only the strongest dataset.
- Caption must state that dominant expression is not functional replacement.

### Figure 5 — Gene-tree evidence and limits of paralog interpretation

- A-B: pruned, support-annotated gene trees for OG0000378 and OG0000114 with dominant genes highlighted.
- C: OG0003912 (PSBO) only as a falsifiable switching hypothesis if alignment/tree and independent-expression checks pass.
- D: null-distribution position for all 19 priority families ({patterns.get('conserved_subclade_signal', 0)} compact, {patterns.get('dispersed_drivers_possible_paralog_switching', 0)} dispersed, {patterns.get('indeterminate_tree_pattern', 0)} indeterminate).
- Main conclusion: gene-tree screening narrows rather than proves evolutionary mechanism.

### Figure 6 — Robustness and boundary conditions

- Leave-one-study/species recurrence of each conserved program, mapping coverage, label confidence, and input-expression class.
- Sensitivity to marker cut-off and orthogroup family size.
- Negative results and non-replicating families remain visible.
- This figure is pending independent validation and should not be replaced by model results from the separate benchmark task.

## Tables and supplementary items

- Table 1: frozen cohort, provenance, tissue, cells, labels, expression class, mapping coverage, and exclusion reason.
- Table 2: candidate programs with species recurrence, study recurrence, effect summary, annotation evidence, and claim tier.
- Supplementary Figure S1: governance and label audit.
- Supplementary Figure S2: all 102 recurrent candidates.
- Supplementary Figure S3: all 19 priority gene-tree null results.
- Supplementary Table S1: source manifest and checksums.
- Supplementary Tables S2-S4: vocabulary, mapping decisions, and frozen panels/splits.
- Supplementary Tables S5-S7: candidate annotations, gene-level effects, and family dossiers.

## Claim gates

1. Use “conserved program” only after independent study/species recurrence is quantified.
2. Use “paralog switching/substitution” only with curated orthology, supported gene tree, reciprocal expression evidence, and independent biological support.
3. Keep OG0003912 in the supplement unless every gate above passes.
4. Report study/species as the replication unit; never use millions of cells as independent biological replicates.
5. Keep model-performance claims outside this manuscript strand until the separate benchmark protocol is complete and reconciled.
"""
    plan_path = DOCS / "jeb_biology_figure_plan_v1.md"
    plan_path.write_text(figure_plan, encoding="utf-8")

    scope_note = """# Cross-task scope comparison (v1)

## Finding

No direct frozen-label conflict was established by the current read-only comparison, but the two tasks use intentionally different scopes and must not be silently merged.

- PhyloOpenCell governance/discovery uses leaf, root, and vascular panels and a 25-label pilot freeze.
- The independent `phylo_plant_fm_pilot` benchmark manifest is a 10-dataset root-focused cohort with its own broad root-family mapping.
- Several model-task labels are broader collapses (for example companion cell into phloem and stele states into vascular-undifferentiated). These are appropriate only for that benchmark protocol and do not replace the PhyloOpenCell frozen vocabulary.
- Dataset roles also differ: a source may be in the broader governance panel while having train, study-holdout, or species-holdout status in the model benchmark.

## Safety decision

No label, dataset role, split, or source file was changed. Any future shared analysis must join by exact dataset ID and retain both task-specific label columns and protocol version.
"""
    scope_path = REPORTS / "cross_task_scope_comparison_v1.md"
    scope_path.write_text(scope_note, encoding="utf-8")

    outputs = [matrix_path, dossier_path, plan_path, scope_path]
    audit = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "jeb_biology_evidence_package_complete",
        "evidence_claims": len(evidence_rows),
        "claim_status_counts": {},
        "priority_family_dossiers": len(dossier_rows),
        "gene_tree_pattern_counts": patterns,
        "cross_task_scope_finding": "different_scopes_no_direct_frozen_label_conflict_established",
        "outputs": {str(p.relative_to(ROOT)): {"sha256": sha256(p)} for p in outputs},
        "safety": {
            "source_h5ad_accessed": False,
            "model_training_started": False,
            "external_queries": False,
            "other_task_written": False,
            "existing_v1_run_overwritten": False,
        },
        "remaining_gates": [
            "independent study/species validation of conserved programs",
            "manual alignment and branch-support review for selected gene families",
            "curated external biological evidence before any paralog-switching claim",
        ],
    }
    for row in evidence_rows:
        status = row["current_status"]
        audit["claim_status_counts"][status] = audit["claim_status_counts"].get(status, 0) + 1
    audit_path = REPORTS / "jeb_biology_evidence_package_v1_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
