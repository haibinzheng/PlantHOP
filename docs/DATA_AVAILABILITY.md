# Data availability

## Repository policy

This Git repository contains code, frozen metadata, compact evidence tables, checksums, and figures. It does not contain the approximately 400 GB source collection or large standardized objects.

## Public source datasets

The analysis uses publicly released plant single-cell datasets. Dataset identifiers, species, provenance, integrity status, eligibility decisions, and cell totals are recorded in:

- `results/supplementary/Table_S1_dataset_inventory.tsv`
- `results/supplementary/input_source_manifest.tsv`

External validation resources include the public accessions described by the manuscript and tables, including GSE232863, GSE268881, GSE270342, and PRJNA1186371. Users should obtain the original files from their authoritative repositories and comply with the corresponding source terms.

## Expected local data classes

Reproduction requires local paths for:

- standardized H5AD objects;
- external H5AD or Seurat/RDS validation objects;
- reference annotation and protein FASTA files;
- OrthoFinder outputs or the inputs required to rebuild them.

Large files are intentionally matched by `.gitignore`. Do not commit controlled, personal, or credential-bearing files.

## Integrity

The supplementary package includes SHA-256 manifests and audit JSON files. Missing orthogroup coverage is represented as blank rather than zero. Candidate gene and orthogroup bridges are not claims of experimentally validated one-to-one orthology.

## Planned archive

After the code and path-portability review is complete, a tagged GitHub release will be archived with Zenodo. The resulting DOI and final article links will replace this provisional note.

