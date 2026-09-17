import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const base = ".";
const output = path.join(base, "PhyloOpenCell_conserved_program_candidates_v2.xlsx");
const previewDir = path.join(base, "workbook_previews_v2");
await fs.mkdir(previewDir, { recursive: true });

function parseTsv(text) {
  const lines = text.replace(/\r/g, "").trimEnd().split("\n");
  const headers = lines[0].split("\t");
  return lines.slice(1).filter(Boolean).map((line) => {
    const values = line.split("\t");
    const row = {};
    headers.forEach((h, i) => row[h] = values[i] ?? "");
    return row;
  });
}

const representation = parseTsv(await fs.readFile(path.join(base, "marker_feature_representation_audit_v2.tsv"), "utf8"));
const markers = parseTsv(await fs.readFile(path.join(base, "coarse_marker_orthogroups_v2.tsv"), "utf8"));
const candidates = parseTsv(await fs.readFile(path.join(base, "conserved_program_candidates_v2.tsv"), "utf8"));
const audit = JSON.parse(await fs.readFile(path.join(base, "conserved_program_candidates_v2_audit.json"), "utf8"));
const featureAudit = JSON.parse(await fs.readFile(path.join(base, "marker_rank_features_v2_summary.json"), "utf8"));
const functionEvidence = parseTsv(await fs.readFile(path.join(base, "conserved_program_function_evidence_v2.tsv"), "utf8"));
const functionAudit = JSON.parse(await fs.readFile(path.join(base, "conserved_program_function_evidence_v2_audit.json"), "utf8"));
const paralogSummary = parseTsv(await fs.readFile(path.join(base, "priority_paralog_orthogroup_summary_v2.tsv"), "utf8"));
const treeEvidence = parseTsv(await fs.readFile(path.join(base, "priority_paralog_gene_tree_evidence_v2.tsv"), "utf8"));
const paralogAudit = JSON.parse(await fs.readFile(path.join(base, "priority_paralog_evidence_v2_audit.json"), "utf8"));

const wb = Workbook.create();
const summary = wb.worksheets.add("Summary");
const candidate = wb.worksheets.add("Conserved Candidates");
const functionSheet = wb.worksheets.add("Function Evidence");
const paralogSheet = wb.worksheets.add("Paralog Evidence");
const repr = wb.worksheets.add("Representation");
const marker = wb.worksheets.add("Species Markers");
const methods = wb.worksheets.add("Methods & Sources");
for (const ws of [summary, candidate, functionSheet, paralogSheet, repr, marker, methods]) ws.showGridLines = false;
summary.tabColor = "#1F4E78";
methods.tabColor = "#7F8C8D";

const navy = "#1F4E78", blue = "#D9EAF7", red = "#FCE8E6", amber = "#FFF2CC", green = "#E2F0D9";
const font = "Arial";
function title(ws, address, text) {
  ws.getRange(address).values = [[text]];
  ws.getRange(address).format.font = { name: font, size: 14, bold: true, color: "#1F1F1F" };
}
function header(range) {
  range.format.fill = navy;
  range.format.font = { name: font, size: 10, bold: true, color: "#FFFFFF" };
  range.format.horizontalAlignment = "center";
  range.format.verticalAlignment = "center";
  range.format.borders = { preset: "inside", style: "thin", color: "#FFFFFF" };
}

title(summary, "A2", "PhyloOpenCell conserved-program candidates v2");
summary.getRange("A3:F3").format.borders = { bottom: { style: "thin", color: "#9EADBA" } };
summary.getRange("A5:B9").values = [
  ["Audit status", "Complete exploratory screen"],
  ["Panel × species combinations", null],
  ["Usable combinations", null],
  ["Zero-feature combinations", null],
  ["Conserved candidates (≥3 species)", audit.conserved_candidate_rows],
];
summary.getRange("B6").formulas = [["=COUNTA(Representation!A6:A19)"]];
summary.getRange("B7").formulas = [["=COUNTIF(Representation!G6:G19,\"Usable\")"]];
summary.getRange("B8").formulas = [["=B6-B7"]];
summary.getRange("A5:A9").format.fill = blue;
summary.getRange("A5:A9").format.font = { name: font, size: 10, bold: true };
summary.getRange("B5:B9").format.font = { name: font, size: 10 };
summary.getRange("D5:E8").values = [["Functional evidence", "Count"], ["Priority A candidates", functionAudit.priority_A_rows], ["Candidates with local descriptions", functionAudit.candidate_rows_with_any_description], ["Ready for gene-tree review", paralogAudit.orthogroups_ready_for_gene_tree_mapping]];
header(summary.getRange("D5:E5"));
summary.getRange("D6:D8").format.fill = blue;
summary.getRange("D6:D8").format.font = { name: font, size: 10, bold: true };
summary.getRange("E6:E8").format.font = { name: font, size: 10 };
summary.getRange("B5").format.fill = green;
summary.getRange("B5").format.font = { name: font, size: 10, bold: true, color: "#375623" };
summary.getRange("A11:F11").values = [["Decision", "Use the 102 orthogroups as candidates for biological annotation and independent validation. Do not yet call them conserved functions or paralog substitutions.", "", "", "", ""]];
summary.getRange("A11:F11").format.fill = amber;
summary.getRange("A11:F11").format.font = { name: font, size: 10, bold: true, color: "#7F6000" };
summary.getRange("A11:F11").format.wrapText = true;
summary.getRange("A13:F16").values = [
  ["Finding", "Evidence", "Interpretation", "Impact", "Action", "Status"],
  ["All 14 panel × species combinations are usable", "145,089 of 145,089 cells have nonzero features", "The mapping-filter defect is corrected", "Cross-species recurrence is now testable", "Retain v1 only as audit history", "RESOLVED"],
  ["102 orthogroups recur in at least 3 species", "83 recur in 3 species; 18 in 4; 1 in 5", "A focused annotation set is available", "Supports marker-program follow-up", "Prioritize by recurrence and rank", "CANDIDATES"],
  ["97 candidates contain a multi-gene family", "Orthogroup aggregation removes gene-level attribution", "Family expansion is a hint, not substitution evidence", "Paralog claims remain untested", "Return to gene-level expression", "LIMITATION"],
];
header(summary.getRange("A13:F13"));
summary.getRange("A14:F16").format.font = { name: font, size: 10 };
summary.getRange("A14:F16").format.wrapText = true;
summary.getRange("A14:F16").format.verticalAlignment = "top";
summary.getRange("F14:F15").format.fill = green;
summary.getRange("F16").format.fill = amber;
summary.getRange("F14:F16").format.font = { name: font, size: 10, bold: true };
summary.getRange("H3:J6").values = [["Panel", "Candidates", "Share"], ["leaf", 39, 39/102], ["root", 28, 28/102], ["vascular", 35, 35/102]];
header(summary.getRange("H3:J3"));
const chart = summary.charts.add("bar", summary.getRange("H3:I6"));
chart.title = "Conserved candidates by panel";
chart.titleTextStyle.fontSize = 12;
chart.titleTextStyle.typeface = font;
chart.legend = { position: "top", textStyle: { typeface: font } };
chart.xAxis = { axisType: "textAxis", textStyle: { typeface: font, fontSize: 10 } };
chart.yAxis = { numberFormatCode: "0", numberFormatSourceLinked: false, textStyle: { typeface: font, fontSize: 10 } };
chart.setPosition("H8", "N20");
summary.getRange("A1:N20").format.verticalAlignment = "center";
summary.getRange("A:A").format.columnWidth = 30;
summary.getRange("B:B").format.columnWidth = 28;
summary.getRange("C:F").format.columnWidth = 24;
summary.getRange("D:D").format.columnWidth = 34;
summary.getRange("E:E").format.columnWidth = 18;
summary.getRange("H:J").format.columnWidth = 14;
summary.getRange("J4:J6").format.numberFormat = "0.0%";
summary.getRange("14:16").format.rowHeight = 64;

title(candidate, "A2", "Cross-species coarse-label orthogroup candidates");
candidate.getRange("A3:L3").values = [["Exploratory screen: top-100 within-species markers recurring in at least three species. Functional conservation and paralog substitution require independent gene-level evidence.", "", "", "", "", "", "", "", "", "", "", ""]];
candidate.getRange("A3:L3").format.fill = amber;
candidate.getRange("A3:L3").format.font = { name: font, size: 10, italic: true, color: "#7F6000" };
const candidateHeaders = ["panel", "coarse_label", "orthogroup_id", "species_count", "species", "median_within_species_rank", "mean_rank_score_difference", "minimum_rank_score_difference", "species_with_multi_gene_family", "representative_genes_by_species", "coarse_labels_in_panel", "low_coarse_specificity_flag"];
const candidateRows = candidates.map(r => [r.panel, r.coarse_label, r.orthogroup_id, Number(r.species_count), r.species, Number(r.median_within_species_rank), Number(r.mean_rank_score_difference), Number(r.minimum_rank_score_difference), Number(r.species_with_multi_gene_family), r.representative_genes_by_species, Number(r.coarse_labels_in_panel), r.low_coarse_specificity_flag === "True" ? "Review" : "No"]);
candidate.getRangeByIndexes(4, 0, candidateRows.length + 1, candidateHeaders.length).values = [candidateHeaders, ...candidateRows];
header(candidate.getRange("A5:L5"));
candidate.getRange("A5:L5").format.wrapText = true;
candidate.getRange("5:5").format.rowHeight = 40;
const candidateTable = candidate.tables.add(`A5:L${5 + candidateRows.length}`, true, "ConservedCandidateTable");
candidateTable.style = "TableStyleMedium2";
candidate.getRange(`G6:H${5 + candidateRows.length}`).format.numberFormat = "0.0000";
candidateRows.forEach((row, i) => { if (row[11] === "Review") candidate.getRange(`L${6 + i}`).format.fill = amber; });
candidate.getRange("A:L").format.font = { name: font, size: 9 };
candidate.getRange("A:A").format.columnWidth = 11;
candidate.getRange("B:B").format.columnWidth = 28;
candidate.getRange("C:D").format.columnWidth = 16;
candidate.getRange("E:E").format.columnWidth = 55;
candidate.getRange("F:I").format.columnWidth = 21;
candidate.getRange("J:J").format.columnWidth = 80;
candidate.getRange("K:L").format.columnWidth = 20;
candidate.freezePanes.freezeRows(5);
candidate.freezePanes.freezeColumns(3);

title(functionSheet, "A2", "Local functional evidence for conserved candidates");
functionSheet.getRange("A3:L3").values = [["Descriptions come only from local Ensembl Plants and UniProt FASTA headers. They support prioritization but do not validate cell-type specificity or mechanism.", "", "", "", "", "", "", "", "", "", "", ""]];
functionSheet.getRange("A3:L3").format.fill = amber;
functionSheet.getRange("A3:L3").format.font = { name: font, size: 10, italic: true, color: "#7F6000" };
const functionHeaders = ["priority_tier", "panel", "coarse_label", "orthogroup_id", "species_count", "median_within_species_rank", "mean_rank_score_difference", "arabidopsis_genes", "arabidopsis_symbols", "local_descriptions", "annotation_coverage_fraction", "paralog_interpretation"];
const functionRows = functionEvidence.map(r => [r.priority_tier, r.panel, r.coarse_label, r.orthogroup_id, Number(r.species_count), Number(r.median_within_species_rank), Number(r.mean_rank_score_difference), r.arabidopsis_genes, r.arabidopsis_symbols, r.local_descriptions, Number(r.annotation_coverage_fraction), r.paralog_interpretation]);
functionSheet.getRangeByIndexes(4, 0, functionRows.length + 1, functionHeaders.length).values = [functionHeaders, ...functionRows];
header(functionSheet.getRange("A5:L5"));
functionSheet.getRange("A5:L5").format.wrapText = true;
functionSheet.getRange("5:5").format.rowHeight = 40;
const functionTable = functionSheet.tables.add(`A5:L${5 + functionRows.length}`, true, "FunctionEvidenceTable");
functionTable.style = "TableStyleMedium2";
functionSheet.getRange(`G6:G${5 + functionRows.length}`).format.numberFormat = "0.0000";
functionSheet.getRange(`K6:K${5 + functionRows.length}`).format.numberFormat = "0.0%";
functionRows.forEach((row, i) => { if (row[0] === "A") functionSheet.getRange(`A${6 + i}`).format.fill = green; });
functionSheet.getRange("A:L").format.font = { name: font, size: 9 };
functionSheet.getRange("A:B").format.columnWidth = 12;
functionSheet.getRange("C:C").format.columnWidth = 28;
functionSheet.getRange("D:G").format.columnWidth = 19;
functionSheet.getRange("H:I").format.columnWidth = 27;
functionSheet.getRange("J:J").format.columnWidth = 90;
functionSheet.getRange("K:K").format.columnWidth = 18;
functionSheet.getRange("L:L").format.columnWidth = 42;
functionSheet.freezePanes.freezeRows(5);
functionSheet.freezePanes.freezeColumns(4);

const treeByKey = new Map(treeEvidence.map(r => [`${r.panel}|${r.coarse_label}|${r.orthogroup_id}`, r]));
title(paralogSheet, "A2", "Gene-level dominance and gene-tree screening");
paralogSheet.getRange("A3:N3").values = [["Exploratory evidence only. A dispersed tree pattern prioritizes manual review but does not establish paralog substitution.", "", "", "", "", "", "", "", "", "", "", "", "", ""]];
paralogSheet.getRange("A3:N3").format.fill = amber;
paralogSheet.getRange("A3:N3").format.font = { name: font, size: 10, italic: true, color: "#7F6000" };
const paralogHeaders = ["panel", "coarse_label", "orthogroup_id", "species_with_positive_gene_level_marker", "species_with_multi_gene_family_and_positive_driver", "species_with_distributed_paralog_signal", "dominant_genes_by_species", "gene_level_followup_status", "gene_tree_pattern", "distance_percentile", "observed_mean_pairwise_tree_distance", "random_distance_median", "missing_tree_leaves", "claim_status"];
const paralogRows = paralogSummary.map(r => {
  const t = treeByKey.get(`${r.panel}|${r.coarse_label}|${r.orthogroup_id}`) ?? {};
  return [r.panel, r.coarse_label, r.orthogroup_id, Number(r.species_with_positive_gene_level_marker), Number(r.species_with_multi_gene_family_and_positive_driver), Number(r.species_with_distributed_paralog_signal), r.dominant_genes_by_species, r.gene_level_followup_status, t.gene_tree_pattern ?? "", Number(t.distance_percentile_fraction_random_le_observed ?? 0), Number(t.observed_mean_pairwise_tree_distance ?? 0), Number(t.random_mean_pairwise_distance_median ?? 0), t.missing_dominant_tree_leaves ?? "", "not established"];
});
paralogSheet.getRangeByIndexes(4, 0, paralogRows.length + 1, paralogHeaders.length).values = [paralogHeaders, ...paralogRows];
header(paralogSheet.getRange("A5:N5"));
paralogSheet.getRange("A5:N5").format.wrapText = true;
paralogSheet.getRange("5:5").format.rowHeight = 46;
const paralogTable = paralogSheet.tables.add(`A5:N${5 + paralogRows.length}`, true, "ParalogEvidenceTable");
paralogTable.style = "TableStyleMedium2";
paralogSheet.getRange(`J6:J${5 + paralogRows.length}`).format.numberFormat = "0.0%";
paralogSheet.getRange(`K6:L${5 + paralogRows.length}`).format.numberFormat = "0.0000";
paralogRows.forEach((row, i) => {
  if (row[8] === "conserved_subclade_signal") paralogSheet.getRange(`I${6 + i}`).format.fill = green;
  if (row[8] === "dispersed_drivers_possible_paralog_switching") paralogSheet.getRange(`I${6 + i}`).format.fill = amber;
});
paralogSheet.getRange("A:N").format.font = { name: font, size: 9 };
paralogSheet.getRange("A:A").format.columnWidth = 11;
paralogSheet.getRange("B:B").format.columnWidth = 28;
paralogSheet.getRange("C:F").format.columnWidth = 19;
paralogSheet.getRange("G:G").format.columnWidth = 95;
paralogSheet.getRange("H:I").format.columnWidth = 38;
paralogSheet.getRange("J:N").format.columnWidth = 20;
paralogSheet.freezePanes.freezeRows(5);
paralogSheet.freezePanes.freezeColumns(3);

title(repr, "A2", "Feature representation audit by panel and species");
const reprHeaders = ["panel", "species", "cells", "nonzero_feature_cells", "nonzero_cell_fraction", "matrix_nnz", "usable_for_marker_analysis"];
const reprRows = representation.map(r => [r.panel, r.species, Number(r.cells), Number(r.nonzero_feature_cells), Number(r.nonzero_cell_fraction), Number(r.matrix_nnz), r.usable_for_marker_analysis === "True" ? "Usable" : "Blocked"]);
repr.getRangeByIndexes(4, 0, reprRows.length + 1, reprHeaders.length).values = [reprHeaders, ...reprRows];
header(repr.getRange("A5:G5"));
repr.getRange("A5:G5").format.wrapText = true;
repr.getRange("5:5").format.rowHeight = 34;
const reprTable = repr.tables.add(`A5:G${5 + reprRows.length}`, true, "RepresentationTable");
reprTable.style = "TableStyleMedium2";
repr.getRange(`E6:E${5 + reprRows.length}`).format.numberFormat = "0.0%";
reprRows.forEach((row, i) => {
  if (row[6] === "Blocked") {
    repr.getRange(`G${6 + i}`).format.fill = red;
    repr.getRange(`G${6 + i}`).format.font = { name: font, size: 10, bold: true, color: "#9C0006" };
  }
});
repr.getRange("A:G").format.font = { name: font, size: 10 };
repr.getRange("A:A").format.columnWidth = 12;
repr.getRange("B:B").format.columnWidth = 26;
repr.getRange("C:F").format.columnWidth = 20;
repr.getRange("G:G").format.columnWidth = 28;
repr.freezePanes.freezeRows(5);

title(marker, "A2", "Exploratory within-species coarse-label marker orthogroups");
marker.getRange("A3:J3").values = [["Within-species marker rankings supporting the cross-species candidate table. These are rank-score contrasts, not differential-expression tests.", "", "", "", "", "", "", "", "", ""]];
marker.getRange("A3:J3").format.fill = amber;
marker.getRange("A3:J3").format.font = { name: font, size: 10, italic: true, color: "#7F6000" };
const markerHeaders = ["panel", "coarse_label", "species", "orthogroup_id", "within_species_rank", "mean_rank_score_difference", "class_cells", "other_cells", "representative_genes", "multi_gene_family_in_species"];
const markerRows = markers.map(r => [r.panel, r.coarse_label, r.species, r.orthogroup_id, Number(r.within_species_rank), Number(r.mean_rank_score_difference), Number(r.class_cells), Number(r.other_cells), r.representative_genes, r.multi_gene_family_in_species === "True"]);
marker.getRangeByIndexes(4, 0, markerRows.length + 1, markerHeaders.length).values = [markerHeaders, ...markerRows];
header(marker.getRange("A5:J5"));
marker.getRange("A5:J5").format.wrapText = true;
marker.getRange("5:5").format.rowHeight = 40;
const markerTable = marker.tables.add(`A5:J${5 + markerRows.length}`, true, "SpeciesMarkerTable");
markerTable.style = "TableStyleMedium2";
marker.getRange(`F6:F${5 + markerRows.length}`).format.numberFormat = "0.0000";
marker.getRange("A:J").format.font = { name: font, size: 9 };
marker.getRange("A:A").format.columnWidth = 11;
marker.getRange("B:B").format.columnWidth = 29;
marker.getRange("C:C").format.columnWidth = 24;
marker.getRange("D:D").format.columnWidth = 16;
marker.getRange("E:H").format.columnWidth = 18;
marker.getRange("I:I").format.columnWidth = 34;
marker.getRange("J:J").format.columnWidth = 22;
marker.freezePanes.freezeRows(5);
marker.freezePanes.freezeColumns(2);

title(methods, "A2", "Methods, provenance and interpretation limits");
const methodRows = [
  ["Item", "Value"],
  ["Purpose", "Audit whether frozen rank-orthogroup features can support conserved coarse-label marker discovery."],
  ["Comparison", audit.method.comparison],
  ["Marker ranking", audit.method.ranking],
  ["Cell thresholds", `Class ≥${audit.method.minimum_class_cells}; comparator ≥${audit.method.minimum_other_cells}`],
  ["Conservation threshold", `Top ${audit.method.top_n_per_species_class} in at least ${audit.method.minimum_species_recurrence} species`],
  ["Representation check", "All 145,089 sampled cells have nonzero v2 features across 14 panel × species combinations."],
  ["V1 correction", "V2 includes deterministic_transform_candidate mappings; v1 remains unchanged as audit history."],
  ["Safety", "Source H5AD files were read-only; no model training was run in this audit."],
  ["Source", "data/runs/PhyloOpenCell/marker_rank_features_v2/dataset_feature_audit.tsv"],
  ["Source", "data/PhyloOpenCell/mappings/gene_mapping_candidates_v1.tsv.gz"],
  ["Source", "./metadata/coarse_label_hierarchy_v1.tsv"],
  ["Audit JSON", "./reports/conserved_program_candidates_v2_audit.json"],
  ["Function evidence", "./metadata/conserved_program_function_evidence_v2.tsv"],
  ["Paralog evidence", "./metadata/priority_paralog_gene_tree_evidence_v2.tsv"],
  ["Interpretation limit", "Orthogroup recurrence is candidate evidence only. It does not establish conserved function, causality or gene-level paralog substitution."],
];
methods.getRangeByIndexes(4, 0, methodRows.length, 2).values = methodRows;
header(methods.getRange("A5:B5"));
methods.getRange("A6:A20").format.fill = blue;
methods.getRange("A6:A20").format.font = { name: font, size: 10, bold: true };
methods.getRange("B6:B20").format.font = { name: font, size: 10 };
methods.getRange("B6:B20").format.wrapText = true;
methods.getRange("A:A").format.columnWidth = 25;
methods.getRange("B:B").format.columnWidth = 90;
methods.getRange("6:20").format.rowHeight = 32;

wb.recalculate();
const checks = {
  summary: await wb.inspect({ kind: "region", sheetId: "Summary", range: "A1:J16", maxChars: 6000 }),
  candidates: await wb.inspect({ kind: "region", sheetId: "Conserved Candidates", range: "A1:L20", maxChars: 6000 }),
  functionEvidence: await wb.inspect({ kind: "region", sheetId: "Function Evidence", range: "A1:L20", maxChars: 6000 }),
  paralogEvidence: await wb.inspect({ kind: "region", sheetId: "Paralog Evidence", range: "A1:N24", maxChars: 6000 }),
  representation: await wb.inspect({ kind: "region", sheetId: "Representation", range: "A1:G19", maxChars: 6000 }),
  formulas: await wb.inspect({ kind: "formula", sheetId: "Summary", range: "A1:J20", maxChars: 3000, options: { maxResults: 20 } }),
  errors: await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, maxChars: 3000 }),
};
await fs.writeFile(path.join(previewDir, "inspection.json"), JSON.stringify(checks, null, 2));
for (const [sheetName, fileName, range] of [["Summary", "summary.png", "A1:N20"], ["Conserved Candidates", "conserved_candidates.png", "A1:L30"], ["Function Evidence", "function_evidence.png", "A1:L30"], ["Paralog Evidence", "paralog_evidence.png", "A1:N24"], ["Representation", "representation.png", "A1:G19"], ["Species Markers", "species_markers.png", "A1:J30"], ["Methods & Sources", "methods_sources.png", "A1:B20"]]) {
  const preview = await wb.render({ sheetName, range, scale: 1, format: "png" });
  await fs.writeFile(path.join(previewDir, fileName), new Uint8Array(await preview.arrayBuffer()));
}
const xlsx = await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(output);
console.log(JSON.stringify({ output, sheets: 7, candidateRows: candidateRows.length, functionRows: functionRows.length, paralogRows: paralogRows.length, markerRows: markerRows.length, representationRows: reprRows.length, auditStatus: audit.status }));
