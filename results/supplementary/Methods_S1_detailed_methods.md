# Supplementary Methods S1

## Scope and evidence boundary

This supplement documents the executed PhyloOpenCell evidence workflow supporting the PlantHOP submission. It describes source-data governance, orthogroup construction, rank-based program discovery, frozen validation, and post hoc robustness analyses. It does not present a newly trained predictive model. The supplementary-completion step read standardized metadata and existing audit products only, created versioned derivative tables, and did not alter source H5AD or RDS objects.

All labels are author-provided or previously governed source annotations. No label was inferred from expression during preparation of Tables S2 and S3. Orthogroup recurrence is treated as an association and portability device, not as proof of conserved molecular function, one-to-one orthology, causal marker activity, or paralog substitution. Missing values remain blank with explicit gap fields; they were never replaced by zero.

## Computing environment and recorded software

The supplementary-completion builder ran on the project server under Python 3.12.14 with anndata 0.13.3.post0, NumPy 2.5.3, pandas 3.0.5, and SciPy 1.18.0. Metadata were read in backed mode so that H5AD expression matrices were not loaded or inspected. The P0 nested-label analysis used the project R 4.4.3 runtime. A current environment query returned Seurat 5.5.1; because the original P0 audit did not serialize every package version, this Seurat value records the retained environment rather than independently proving package state for every earlier run.

Orthogroups were generated with OrthoFinder 3.1.5 using DIAMOND sequence search and the recorded MSA/FAMSA/FastTree workflow. The seven-species run used 12 threads, contained 183,178 representative proteins, produced 18,383 orthogroups, and assigned 164,518 genes to orthogroups. Zea mays was excluded from this orthogroup run because a source-matched Zm00001d.1 protein FASTA was unavailable.

## Input inventory and read-only audit

The governed manifest contained 121 standardized H5AD datasets. For Table S2, each readable object was opened with anndata in read-only backed mode and only obs metadata were accessed. The standardized source-label column cell_type_original was preferred, with celltype_after permitted only as a recorded fallback. The scan produced 924 unique dataset–source-label pairs representing 212 distinct source-label strings and 3120007 cells. The corresponding manifest total was 3120007 cells.

Dataset identity, accession, species, tissue, matrix storage, and source-annotation checksums were taken from dataset_manifest.tsv and dataset_source_metadata.tsv. Source-label counts were calculated directly from the standardized annotation field. Expression values were not consulted. A dataset lacking a recognized source-label field would have been retained as a source_unavailable row rather than omitted; the audit records the number of such datasets.

## Label governance and hierarchical mapping

The global ontology candidate inventory contains 212 distinct source-label strings and explicitly retains Unknown. Table S2 maps observed dataset-label pairs only through frozen governance assets. The panel label freeze contains 29 panel–label decisions, of which 27 rows are approved for the pilot hierarchy and two source labels are excluded. The accepted hierarchy maps fine source labels to coarse root, leaf, or vascular families. Many-to-one mappings preserve the original source label, ontology identifier, governed label, parent panel, mapping rule, approval scope, and cell count.

A dataset-label pair receives mapped status only when the frozen dataset-panel label-count authority and accepted hierarchy support that panel-specific mapping. Unknown labels remain unknown. Explicit exclusions remain excluded. Multiple or ambiguous ontology candidates remain ambiguous. All other observed labels remain unmapped, even when they have an exact lexical ontology candidate, because an ontology suggestion is not equivalent to a governed analysis label. This rule prevents expression-based relabeling and preserves the complete unresolved vocabulary.

## Reference mapping, orthogroup construction, and bridges

Source H5AD var_names were mapped to species reference-gene candidates with the frozen gene_mapping_candidates_v1 bridge. The mapping audit covered 81 datasets from eight species, contained 1,611,655 feature rows, mapped 1,587,599 candidate rows, and reported no ambiguous candidate rows, target collisions, or gene rows with multiple representative proteins. These mappings remained candidates pending acceptance review and were not interpreted as experimental orthology validation.

Representative proteins were grouped by the seven-species OrthoFinder run. The resulting gene_orthogroup_candidates_v1 bridge linked mapped reference genes to 18,383 candidate orthogroups. Dataset-level orthogroup coverage was calculated as orthogroup-assigned source features divided by total source features. Seventy-one of 81 evaluated datasets passed the provisional 60% gate. Table S3 expands these dataset-level values across the three frozen panels, retaining blank mapping fields for manifest datasets without an evaluated orthogroup bridge.

## Rank encoding and feature construction

The executed rank representation used nonnegative source expression values from the frozen input layer. Within each cell, nonzero features were ordered by decreasing value, with exact ties resolved by immutable source feature order. The top 128 features were assigned weights 1/log2(rank+1), multiple source features mapping to the same orthogroup were collapsed by their maximum weight, and the sparse orthogroup vector was L2-normalized per cell. This design reduced dependence on incomparable expression magnitudes while retaining within-cell feature ordering.

Marker-rank feature v2 was created for the approved leaf, root, and vascular panels in read-only backed mode. It contained 40,417 leaf cells, 74,821 root cells, and 29,851 vascular cells over the 18,383-orthogroup vocabulary. All 65 included dataset-panel combinations had nonzero mapped features. Rank features were intended for cross-dataset marker and program analysis, not differential-expression inference.

## Program discovery

Within each species and panel, each accepted coarse class was contrasted against the other accepted coarse classes. Orthogroups with positive mean normalized-rank-score differences were sorted in descending order. The frozen discovery rule retained the top 100 orthogroups per eligible species–class comparison, required at least 100 class cells and 100 other cells, and nominated candidates recurring in at least three species.

The v2 discovery generated 4,200 species-level marker rows from 42 eligible comparisons and 102 recurrent candidate rows: 39 leaf, 28 root, and 35 vascular candidates. Candidates arose from the same observational discovery collections and therefore were not independent validation. Recurrence did not establish conserved function, and multi-gene orthogroups required gene-level and phylogenetic follow-up.

## Freeze protocol and leakage controls

Labels, dataset panels, coarse hierarchy, candidate programs, contrasts, scoring parameters, random seeds, and decision rules were frozen before each corresponding scoring run. Versioned outputs were written to new directories with partial-file or partial-directory staging and refusal to overwrite authoritative results. Source objects were opened read-only. Hash manifests recorded input and output identities, and each authority audit documented whether source objects, frozen protocols, labels, thresholds, or program membership changed.

Internal leave-one-study-out analyses selected the dominant gene using all but one dataset within an orthogroup and species, then evaluated the selected gene in the omitted dataset. This reduced same-study selection bias but remained internal evidence. External validation programs were pre-discovered, and target contrasts were frozen before scoring. Missing mappings and excluded labels were not replaced with related genes or labels.

## Validation rounds

Study-level robustness was assessed for 19 priority orthogroups by leaving out each available study within a species. The resulting family statuses summarize direction and dominant-gene stability across held-out studies. Clean supported trees were then used for descriptive dominant-gene compactness analyses. Neither different dominant genes nor dispersed tree positions were interpreted as paralog substitution, and the final P1 audit found no family satisfying the strict joint expression-plus-supported-subclade criterion.

Independent rice validation v2.1 used the RNA counts layer from GSE232863, encoded 49,273 selected cells, mapped 14,136 orthogroups, used top K=128, and compared each frozen program with 500 random programs under seed 20260914. Version 2.1 corrected only the per-contrast sample-direction filter; the frozen programs, contrasts, thresholds, and seed were unchanged. Three of five primary contrasts passed, with no majority-direction reversals. Replicate-level raw rows were not available in the frozen authority package and were not reconstructed.

External validation v3 froze the raw-count rank-encoding implementation, complete program score, single-gene and expanded-family subsets, off-target controls, 500 random programs, and replicate-direction rule before expression scoring. The evidence tables retain the Populus phloem negative nonreplication result and distinguish unavailable values from zero. GSE268881 Esa and Sir served as nested-label cell-level discrimination and within-study repeatability evidence, not broad species-generalization evidence; Spa remained ineligible under the frozen metadata gate.

## Nested labels, circularity stress test, and parameter robustness

The P0 analysis used the GSE268881 RNA counts layer and three frozen nested endpoints: xylem versus all other root labels; xylem versus pericycle, phloem, and procambium; and stele/vascular labels versus non-stele root labels. The complete 18-orthogroup program was evaluated in pooled cells and separately in R1 and R2 for Esa and Sir. The K=128 immutable-source-order condition exactly reproduced the prior authority pooled metrics and replicate AUROCs to 10^-12.

Parameter robustness covered K values 64, 128, and 256 and the immutable source order plus ten deterministic global feature-order permutations. Across the complete program, all 18 species–endpoint–scope groups retained AUROC above 0.5 and positive standardized mean difference in all 33 conditions per group. These direction checks are descriptive robustness results rather than new confirmatory endpoints.

Circularity sensitivity removed five of the 18 frozen orthogroups documented as overlapping published marker evidence. The complete label-assignment or label-transfer feature set was unavailable. Therefore the analysis supports robustness to removal of documented overlap only and does not completely exclude label circularity. Integrated scale.data was not used as a substitute marker set.

## Discovery-threshold and clade-balance sensitivity

P1 reconstructed the frozen top-100/minimum-three-species candidate set exactly. It then evaluated top-N values of 25, 50, and 100 and recurrence thresholds of two, three, four, and five species without replacing the primary rule. At a minimum of three species, 16, 35, and 102 candidates were retained for top 25, 50, and 100, respectively.

Cross-clade presence required at least one eudicot and one monocot, whereas the strict balance descriptor required at least two species from each clade. Among the frozen 102 candidates, 62 had cross-clade presence and seven met the strict two-plus-two descriptor. Leaf collections were eudicot-only and thus could not pass a cross-clade rule; this was recorded as a design limitation rather than negative evidence. Six of the submitted 18 xylem-program orthogroups were independently rediscovered by the frozen vascular-xylem top-100/minimum-three rule.

## Statistical summaries

Cell-level discrimination was summarized by AUROC, average precision, and standardized mean difference. Average precision was interpreted relative to class prevalence. Replicate scopes were reported separately from pooled cells. Random-program comparisons used the frozen empirical probability (1 + number of random AUROCs at least as large as the target AUROC)/(number of random programs + 1). Direction rules used AUROC greater than 0.5 and positive standardized mean difference; these rules were frozen and were not optimized after outcome inspection.

No missing metric was converted to zero. Counts and mapping fractions in Table S3 are stored as typed numeric values. Each nonmissing fraction was independently recomputed as mapped_features/total_features, checked to lie within [0,1], and required to agree with the authority value within 10^-12. Descriptive sensitivity analyses did not create new hypothesis-testing thresholds.

## Reproducibility, checksums, and table validation

The dataset manifest SHA-256 is 74218fbeb25ed470fd9ee6802a7c08ab3c5a4124abb092fd7d34f3cf7f683b35; the label-freeze hash is 78a020ef314246281f19dc5eec25dbb28c50c5c376ee9df5a6b5090a0ec7d321; the hierarchy hash is e139f4d944b153210afd8aef98174679bf7caead33207fadc3a2e8c4cb619abb; and the orthogroup-coverage hash is 5bf556f16a0a3dc2bc78cea55cd5a1e852aef9f302ee79aee0e1baf98f8104b3. The build records all input paths, file sizes, and SHA-256 values in input_manifest.tsv and all deliverable hashes in output_sha256.tsv.

Table validation required unique record keys, allowed mapping-status enumerations, nonnegative cell and feature counts, valid fraction bounds, exact numerator/denominator reconciliation, and expected panel enumeration. The audit also compared scanned cell counts with the manifest, retained all observed label strings, checked the 29-row governed dictionary, and verified that unavailable mapping values remained blank. The Markdown and LaTeX methods files are generated from the same section source to prevent prose divergence.

## Safety boundaries and known gaps

The completion run did not train a model, download data, write to source datasets, change a label, modify an orthogroup bridge, alter a threshold, or overwrite an authority directory. H5AD access was limited to obs metadata. Earlier RDS validation runs were not repeated. Every output was written beneath a new supplementary_completion_v1 directory.

The principal table gaps are explicit. Dataset-specific source-label field names before standardization are unavailable for many source annotations; Table S2 therefore reports the standardized H5AD field and source-annotation hash. Only 29 panel–label governance decisions and 27 accepted hierarchy rows exist, so other observed labels remain unknown, excluded, ambiguous, or unmapped. Orthogroup coverage is available for 81 of 121 manifest datasets; the remaining 40 datasets retain blank numerators and fractions across their panel rows. Gene and orthogroup mappings are candidate bridges rather than validated one-to-one ortholog assignments. Functional conservation, causality, complete absence of label circularity, and paralog substitution remain unestablished.
