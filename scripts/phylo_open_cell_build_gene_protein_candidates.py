#!/usr/bin/env python3
"""Build auditable gene-to-protein candidates and longest-isoform FASTA inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


PROTEOMES = {
    "Arabidopsis thaliana": ("Arabidopsis_thaliana.TAIR10.pep.all.fa.gz", "ensembl"),
    "Brassica rapa": ("Brassica_rapa.Brapa_1.0.pep.all.fa.gz", "ensembl"),
    "Glycine max": ("Glycine_max.Glycine_max_v2.1.pep.all.fa.gz", "ensembl"),
    "Oryza sativa": ("Oryza_sativa.IRGSP-1.0.pep.all.fa.gz", "ensembl"),
    "Populus trichocarpa": ("Populus_trichocarpa.Pop_tri_v4.pep.all.fa.gz", "ensembl"),
    "Sorghum bicolor": ("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz", "ensembl"),
    "Catharanthus roseus": ("UP001060085.uniprot.fasta.gz", "uniprot"),
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fasta_records(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = None
        seq = []
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)


def identifiers(header: str, kind: str):
    protein = header.split()[0]
    if kind == "ensembl":
        match = re.search(r"(?:^|\s)gene:([^\s]+)", header)
        return (match.group(1) if match else None), protein
    match = re.search(r"(?:^|\s)GN=([^\s]+)", header)
    if "|" in protein:
        parts = protein.split("|")
        protein = parts[1] if len(parts) > 1 else protein
    return (match.group(1) if match else None), protein


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def atomic_gzip_text(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    return tmp, gzip.open(tmp, "wt", encoding="utf-8", newline="")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--reference-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--fasta-dir", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--audit", type=Path, required=True)
    args = p.parse_args()

    requested = defaultdict(set)
    with gzip.open(args.mapping, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["candidate_reference_gene_id"]
            if row["species"] in PROTEOMES and gene:
                requested[row["species"]].add(gene)

    output_rows = []
    summaries = []
    args.fasta_dir.mkdir(parents=True, exist_ok=True)
    for species, (filename, kind) in PROTEOMES.items():
        subdir = "uniprot" if kind == "uniprot" else "ensembl_plants_release_63"
        path = args.reference_root / subdir / filename
        source_hash = digest(path)
        proteins = defaultdict(list)
        no_gene_header = 0
        for header, sequence in fasta_records(path):
            gene, protein = identifiers(header, kind)
            if gene:
                proteins[gene].append((protein, sequence))
            else:
                no_gene_header += 1

        wanted = requested[species]
        matched = 0
        one_to_many = 0
        representative_count = 0
        fasta_path = args.fasta_dir / f"{slug(species)}.representative_proteins.fa.gz"
        fasta_tmp, fasta_handle = atomic_gzip_text(fasta_path)
        with fasta_handle:
            for gene in sorted(wanted):
                candidates = proteins.get(gene, [])
                if candidates:
                    matched += 1
                    one_to_many += int(len(candidates) > 1)
                    representative = max(candidates, key=lambda item: (len(item[1]), item[0]))
                    protein_ids = sorted(item[0] for item in candidates)
                    fasta_handle.write(f">{gene}|{representative[0]} species={species} source_sha256={source_hash}\n")
                    sequence = representative[1]
                    for start in range(0, len(sequence), 80):
                        fasta_handle.write(sequence[start:start + 80] + "\n")
                    representative_count += 1
                    status = "mapped_candidate"
                else:
                    representative = ("", "")
                    protein_ids = []
                    status = "no_protein_candidate"
                output_rows.append({
                    "species": species,
                    "candidate_reference_gene_id": gene,
                    "candidate_protein_accessions": ";".join(protein_ids),
                    "protein_cardinality": len(protein_ids),
                    "representative_protein_accession": representative[0],
                    "representative_length_aa": len(representative[1]),
                    "mapping_status": status,
                    "selection_rule": "longest_sequence_then_accession" if candidates else "none",
                    "proteome_sha256": source_hash,
                    "accepted": "False",
                })
        os.replace(fasta_tmp, fasta_path)
        summaries.append({
            "species": species,
            "candidate_genes": len(wanted),
            "genes_with_protein": matched,
            "coverage_pct": round(100 * matched / len(wanted), 4) if wanted else 0.0,
            "one_to_many_genes": one_to_many,
            "representative_sequences": representative_count,
            "proteome_records_with_no_gene_tag": no_gene_header,
            "reference_file": str(path),
            "reference_sha256": source_hash,
            "accepted": "False",
        })

    fields = list(output_rows[0])
    tmp, handle = atomic_gzip_text(args.output)
    with handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    os.replace(tmp, args.output)

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    tmp_summary = args.summary.with_suffix(args.summary.suffix + ".partial")
    with tmp_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    os.replace(tmp_summary, args.summary)

    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_mapping": str(args.mapping),
        "input_mapping_sha256": digest(args.mapping),
        "output": str(args.output),
        "output_sha256": digest(args.output),
        "accepted": False,
        "selection_rule": "longest protein sequence; accession lexical order breaks ties",
        "species": summaries,
        "totals": {
            "candidate_genes": sum(x["candidate_genes"] for x in summaries),
            "genes_with_protein": sum(x["genes_with_protein"] for x in summaries),
            "representative_sequences": sum(x["representative_sequences"] for x in summaries),
            "one_to_many_genes": sum(x["one_to_many_genes"] for x in summaries),
        },
        "excluded": {
            "Zea mays": "source-matched Zm00001d.1 protein FASTA was not found in the official MaizeGDB directory; d2 was not substituted",
        },
        "source_h5ad_modified": False,
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    tmp_audit = args.audit.with_suffix(args.audit.suffix + ".partial")
    tmp_audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp_audit, args.audit)
    print(json.dumps(audit["totals"], ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
