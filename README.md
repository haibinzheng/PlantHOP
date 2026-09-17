# PlantHOP

PlantHOP is a hierarchy-aware workflow for testing fixed orthogroup-level plant cell programmes across single-cell transcriptomic datasets and species. The repository accompanies the manuscript **"Label resolution is associated with cross-species transfer of orthogroup cell programmes in plant single-cell atlases."**

> Repository status: public submission snapshot. The manuscript is not yet represented here as a published article.

## What is included

- `scripts/`: executed Python, R, and JavaScript analysis scripts.
- `metadata/`: fixed label-governance and conserved-programme records.
- `results/supplementary/`: machine-readable Tables S1, S1b, S2--S6 and S7a--S7d, field dictionaries and detailed methods.
- `results/figures/`: the current manuscript and supplementary figures distributed with the reproducibility snapshot.
- `environment/`: a portable analysis environment and the exact Linux Conda package export used for the orthogroup workflow.
- `docs/`: data-access and reproduction guidance.

## What is not included

Raw or standardized H5AD/RDS objects, third-party reference proteomes, pretrained model weights, private server details, and manuscript working files are intentionally excluded. Public source accessions and the expected local layout are described in `docs/DATA_AVAILABILITY.md`. Publication and direct-download metadata are included where resolved; unresolved fields remain blank and carry explicit gap flags.

## Quick start

Create the portable analysis environment:

```bash
conda env create -f environment/environment.yml
conda activate planthop
```

Inspect the available command-line scripts:

```bash
python scripts/audit_h5ad_collection.py --help
python scripts/phylo_open_cell_governance.py --help
python scripts/build_marker_rank_features_v2.py --help
python scripts/build_conserved_orthogroup_markers.py --help
```

Run scripts from the repository root and place third-party inputs under a local `data/` directory, or override their command-line paths where available. See `docs/REPRODUCIBILITY.md` for the workflow map and known evidence gaps.

## Main analysis stages

1. Audit standardized single-cell datasets without modifying source objects.
2. Govern author-provided labels through frozen cell-system and lineage mappings.
3. Map source features to reference genes and candidate orthogroups.
4. Encode within-cell expression ranks as sparse orthogroup activity.
5. Discover recurrent programmes under fixed thresholds.
6. Specify programmes, contrasts, controls, seeds and decision rules before held-out scoring.
7. Evaluate external and nested-label transfer without retraining.
8. Build machine-readable supplementary evidence and checksums.

## Reuse and citation

Code in this repository is released under the MIT License. Third-party datasets and software remain governed by their original licenses. Citation metadata are provided in `CITATION.cff`; a DOI will be added when a public archived release is created.

## Contact

For repository questions, use the GitHub issue tracker. Scientific correspondence for the manuscript should follow the corresponding-author information in the article.
