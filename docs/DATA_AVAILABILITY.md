# Data availability

## Repository policy

This Git repository contains code, frozen metadata, compact evidence tables, checksums, and figures. It does not contain the approximately 400 GB source collection or large standardized objects.

## Public source datasets

The analysis uses publicly released plant single-cell datasets. Dataset identifiers, repository accessions, species, provenance, integrity status, eligibility decisions and cell totals are recorded in:

- `results/supplementary/Table_S1_dataset_inventory.tsv`

Source-publication and direct-download metadata are reported where resolved. Blank fields are not claims that the information does not exist; they identify metadata that was not resolved in the submitted inventory.

External validation resources include the public accessions described by the manuscript and tables, including GSE232863, GSE268881, GSE270342, and PRJNA1186371. Users should obtain the original files from their authoritative repositories and comply with the corresponding source terms.

## Expected local data classes

Reproduction requires local paths for:

- standardized H5AD objects;
- external H5AD or Seurat/RDS validation objects;
- reference annotation and protein FASTA files;
- OrthoFinder outputs or the inputs required to rebuild them.

Large files are intentionally matched by `.gitignore`. Do not commit controlled, personal, or credential-bearing files.

## Integrity

Supplementary tables identify the scope of available SHA-256 values and distinguish direct source-object checksums from metadata or derived-record hashes. Missing orthogroup coverage is represented as blank rather than zero. Candidate gene and orthogroup bridges are not claims of experimentally validated one-to-one orthology. Internal source paths are labelled as not publicly portable.

## Planned archive

After the code and path-portability review is complete, a tagged GitHub release will be archived with Zenodo. The resulting DOI and final article links will replace this provisional note.
