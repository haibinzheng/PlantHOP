import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def project_gene_id(token):
    return token.strip().split("|", 1)[0]


def read_project_gene_map(path):
    gene_to_ogs = defaultdict(set)
    gene_to_species = {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reference_columns = {
            "ath": "arabidopsis_thaliana.representative_proteins",
            "osa": "oryza_sativa.representative_proteins",
        }
        missing = set(reference_columns.values()) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing project orthogroup columns: {sorted(missing)}")
        for row in reader:
            orthogroup = row["Orthogroup"]
            for species, column in reference_columns.items():
                for token in row[column].split(","):
                    gene = project_gene_id(token)
                    if gene:
                        gene_to_ogs[gene].add(orthogroup)
                        gene_to_species[gene] = species
    return gene_to_ogs, gene_to_species


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--author-map", required=True)
    parser.add_argument("--project-orthogroups", required=True)
    parser.add_argument("--wheat-features", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--output-audit", required=True)
    args = parser.parse_args()

    project_gene_to_ogs, project_gene_to_species = read_project_gene_map(args.project_orthogroups)

    author_rows = []
    author_counts = Counter()
    author_unique_genes = defaultdict(set)
    gf_to_reference_matches = defaultdict(list)
    wheat_gene_to_gfs = defaultdict(set)
    with open(args.author_map, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene_id", "species", "gf_id"]:
            raise ValueError(f"Unexpected author-map columns: {reader.fieldnames}")
        for row in reader:
            gene = row["gene_id"].strip()
            species = row["species"].strip()
            gf_id = row["gf_id"].strip()
            author_rows.append((gene, species, gf_id))
            author_counts[species] += 1
            author_unique_genes[species].add(gene)
            if species == "tae":
                wheat_gene_to_gfs[gene].add(gf_id)
            elif species in {"ath", "osa"} and gene in project_gene_to_ogs:
                for orthogroup in project_gene_to_ogs[gene]:
                    gf_to_reference_matches[gf_id].append(
                        (species, gene, orthogroup)
                    )

    gf_crosswalk = {}
    gf_status = {}
    for gf_id, matches in gf_to_reference_matches.items():
        orthogroups = sorted({match[2] for match in matches})
        species = sorted({match[0] for match in matches})
        if len(orthogroups) == 1:
            gf_crosswalk[gf_id] = orthogroups[0]
            gf_status[gf_id] = (
                "unanimous_two_reference_species"
                if len(species) == 2
                else "unanimous_one_reference_species"
            )
        else:
            gf_status[gf_id] = "ambiguous_multiple_project_orthogroups"

    with open(args.wheat_features, encoding="utf-8-sig") as handle:
        wheat_features = [line.strip() for line in handle if line.strip()]
    if wheat_features and wheat_features[0] == "gene_id":
        wheat_features = wheat_features[1:]
    if len(wheat_features) != len(set(wheat_features)):
        raise ValueError("Wheat feature IDs are not unique")

    output_rows = []
    mapped_project_ogs = set()
    status_counts = Counter()
    for gene in wheat_features:
        author_gfs = sorted(wheat_gene_to_gfs.get(gene, set()))
        accepted_ogs = sorted({gf_crosswalk[gf] for gf in author_gfs if gf in gf_crosswalk})
        if not author_gfs:
            status = "unmapped_author_table"
        elif len(author_gfs) > 1:
            status = "ambiguous_multiple_author_families"
        elif len(accepted_ogs) == 1:
            status = gf_status[author_gfs[0]]
            mapped_project_ogs.add(accepted_ogs[0])
        elif len(accepted_ogs) > 1:
            status = "ambiguous_multiple_project_orthogroups"
        else:
            status = gf_status.get(author_gfs[0], "unresolved_no_shared_reference_gene")
        status_counts[status] += 1
        output_rows.append(
            {
                "wheat_gene_id": gene,
                "author_gf_id": ";".join(author_gfs),
                "project_orthogroup_id": ";".join(accepted_ogs),
                "mapping_status": status,
            }
        )

    output_path = Path(args.output_map)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    candidate_ogs = set()
    candidate_programs = defaultdict(set)
    with open(args.candidates, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            orthogroup = row["orthogroup_id"]
            candidate_ogs.add(orthogroup)
            candidate_programs[(row["panel"], row["coarse_label"])].add(orthogroup)

    mapped_count = sum(
        count for status, count in status_counts.items() if status.startswith("unanimous_")
    )
    audit = {
        "inputs": {
            "author_map": args.author_map,
            "author_map_sha256": sha256(args.author_map),
            "project_orthogroups": args.project_orthogroups,
            "project_orthogroups_sha256": sha256(args.project_orthogroups),
            "wheat_features": args.wheat_features,
            "wheat_features_sha256": sha256(args.wheat_features),
            "candidates": args.candidates,
            "candidates_sha256": sha256(args.candidates),
        },
        "author_mapping": {
            "rows_by_species": dict(sorted(author_counts.items())),
            "unique_genes_by_species": {
                species: len(genes) for species, genes in sorted(author_unique_genes.items())
            },
            "wheat_author_families": len(set().union(*wheat_gene_to_gfs.values())),
            "author_families_with_shared_reference_gene": len(gf_to_reference_matches),
            "author_families_unanimously_crosswalked": len(gf_crosswalk),
            "crosswalk_status_counts": dict(sorted(Counter(gf_status.values()).items())),
        },
        "wheat_feature_mapping": {
            "total_features": len(wheat_features),
            "status_counts": dict(sorted(status_counts.items())),
            "mapped_features": mapped_count,
            "mapped_fraction": mapped_count / len(wheat_features),
            "mapped_project_orthogroups": len(mapped_project_ogs),
        },
        "candidate_program_coverage": {
            "candidate_orthogroups": len(candidate_ogs),
            "covered_candidate_orthogroups": len(candidate_ogs & mapped_project_ogs),
            "coverage_fraction": len(candidate_ogs & mapped_project_ogs) / len(candidate_ogs),
            "by_program": {
                f"{panel}:{label}": {
                    "candidate_orthogroups": len(orthogroups),
                    "covered_candidate_orthogroups": len(orthogroups & mapped_project_ogs),
                    "coverage_fraction": len(orthogroups & mapped_project_ogs) / len(orthogroups),
                }
                for (panel, label), orthogroups in sorted(candidate_programs.items())
            },
        },
        "rule": {
            "description": "Bridge the paper-author family ID to the frozen project orthogroup namespace only when all exact shared Arabidopsis/rice gene IDs point to one project orthogroup; otherwise leave unresolved.",
            "expression_blind": True,
            "accepted_statuses": [
                "unanimous_two_reference_species",
                "unanimous_one_reference_species",
            ],
        },
        "output_map": str(output_path),
        "output_map_sha256": sha256(output_path),
    }
    audit_path = Path(args.output_audit)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
