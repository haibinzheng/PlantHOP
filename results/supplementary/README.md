# PlantHOP supplementary-data submission workspace

Frozen machine-readable tables attached to the current manuscript package:

- `Table_S1_dataset_inventory.tsv`: 127 input records (121 audited resources and six additional validation inputs) covering repository accessions, resolved publication and download metadata, integrity scope, eligibility and cell counts.
- `Table_S1b_discovery_object_checksums.tsv`: filenames, byte sizes and direct SHA-256 checksums for the 44 standardized H5AD objects used in discovery.
- `Table_S2_label_hierarchy_dictionary.tsv`: 29 frozen panel--label governance decisions.
- `Table_S2_source_label_mapping.tsv`: 945 governed dataset--source-label--panel records, covering 924 unique dataset--label pairs and all 3,120,007 audited cells.
- `Table_S3_dataset_orthogroup_mapping_coverage.tsv`: 363 dataset--panel records; 243 rows contain evaluated coverage and 120 retain blank mapping fields because coverage was unavailable.
- `Table_S4_xylem_orthogroup_program.tsv`: all 18 Arabidopsis-derived xylem orthogroups, member-level selection provenance, recurrence status and available deterministic single-gene representatives.
- `Table_S5_replicate_metrics_controls_exclusions.tsv`: 83 pooled or replicate-level records covering validation metrics, controls, paired comparisons, exclusions and explicit source-data gaps; the retained Populus phloem failure is included.
- `Table_S6_label_resolution_circularity_parameter_robustness.tsv`: 54 records covering the three nested label contrasts, documented-marker exclusion and parameter-direction robustness.
- `Table_S7a_discovery_threshold_sensitivity.tsv`: 48 discovery-threshold sensitivity records.
- `Table_S7b_clade_balance_eligibility.tsv`: 11 clade-balance eligibility records.
- `Table_S7c_xylem_program_rediscovery.tsv`: 18 submitted-xylem-program rediscovery records.
- `Table_S7d_single_gene_ortholog_evidence_gaps.tsv`: 19 priority-family evidence-gap records; one-to-one orthology and paralogue substitution remain not established.
- `Methods_S1_detailed_methods.tex`: detailed preprocessing, programme specification, validation and robustness methods used to generate the supplementary PDF.

Supporting dictionaries:

- `field_dictionary.tsv`: field definitions and missing-value policy.
Tables S1--S7d are versioned derivatives of the retained analysis outputs. Blanks are not imputed or converted to zero. Table S4 does not transfer legacy A/B priority tiers into the xylem programme. Legacy internal identifiers containing `prospective` are retained solely for provenance and do not denote public preregistration. Build logs, source manifests and internal integrity records are retained with the project process files rather than included in the reviewer-facing supplement.

The version-controlled reproducibility repository is publicly available at `https://github.com/haibinzheng/PlantHOP`. It includes the portable Conda environment, fixed mappings, evidence tables, checksums and analysis scripts. A versioned archival release with a persistent identifier will be completed no later than publication.
