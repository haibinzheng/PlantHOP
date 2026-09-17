# Reproducibility guide

## Scope of the current snapshot

This snapshot preserves the executed analysis scripts and the compact evidence package used for submission. It supports method inspection, table verification and staged reruns after public inputs are placed locally.

## Executed software record

The supplementary completion step used Python 3.12.14, anndata 0.13.3.post0, NumPy 2.5.3, pandas 3.0.5, and SciPy 1.18.0. Nested-label analyses used R 4.4.3; the retained environment reported Seurat 5.5.1. Orthogroups were generated with OrthoFinder 3.1.5 using DIAMOND, FAMSA, and FastTree in the recorded workflow.

`environment/environment.yml` is the portable environment specification. `environment/conda-linux-64-explicit.txt` is the exact Linux package export retained for the orthogroup environment.

## Workflow map

| Stage | Representative entry points | Main frozen products |
|---|---|---|
| Collection audit | `audit_h5ad_collection.py` | Dataset inventory and integrity audit |
| Label governance | `phylo_open_cell_governance.py`, `phylo_open_cell_governance_candidates.py` | Label decisions and hierarchy |
| Gene/reference bridge | `phylo_open_cell_build_gene_mapping_candidates.py`, `phylo_open_cell_build_gene_protein_candidates.py` | Candidate gene mappings |
| Orthogroups | `phylo_open_cell_prepare_orthofinder_full.py`, `phylo_open_cell_run_orthofinder.py`, `phylo_open_cell_build_orthogroup_candidates.py` | Orthogroup bridge and coverage |
| Rank features | `build_marker_rank_features_v2.py` | Sparse rank-weighted orthogroup features |
| Programme discovery | `build_conserved_orthogroup_markers.py` | Frozen recurrent candidates |
| External validation | `freeze_external_validation_v3.py`, `score_external_validation_v3.py` | Frozen contrasts, metrics, and controls |
| Label resolution | `run_label_resolution_circularity_robustness_v1.R`, `build_label_resolution_circularity_robustness_v1.py` | Nested-label and marker-overlap checks |
| Threshold robustness | `build_discovery_threshold_ortholog_robustness_v1.py` | Tables S7a-S7d |
| Supplement build | `build_supplementary_completion_v1.py`, `build_manuscript_supplement_table_s6_v1.py`, `build_manuscript_supplement_table_s7_v1.py` | Methods S1 and Tables S1-S7 |

## Path configuration

Run scripts from the repository root. Relative paths refer to repository folders and to a local `data/` directory that is intentionally excluded from version control. Scripts with command-line path options may instead be directed to another read-only input location.

## Safety and leakage controls

- Treat source H5AD/RDS objects as read-only.
- Write each run to a new versioned output directory.
- Freeze labels, programs, contrasts, controls, seeds, and decision rules before scoring.
- Preserve missing values; do not convert unavailable evidence to zero.
- Validate input and output SHA-256 values.
- Do not infer or repair labels from expression during metadata governance.

## Known evidence gaps

- Orthogroup coverage is unavailable for 40 of 121 manifest datasets.
- Candidate bridges do not establish one-to-one orthology, functional conservation, causality, or paralog substitution.
- Complete label-transfer features are unavailable for GSE268881, so marker-removal sensitivity does not completely exclude circularity.
- Replicate-level raw rows for the retained rice v2.1 validation archive were not reconstructed.
- Zea mays was excluded from the seven-species OrthoFinder run because a source-matched Zm00001d.1 protein FASTA was unavailable.

## Release checklist

Before the archival release:

1. add a small synthetic smoke-test fixture;
2. run syntax and checksum validation in continuous integration;
3. verify all source accessions and third-party licenses;
4. tag the release and archive it with Zenodo.
