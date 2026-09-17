import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--project-orthogroups", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--output-audit", required=True)
    args = parser.parse_args()

    ath_gene_to_og = {}
    with open(args.project_orthogroups, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        column = "arabidopsis_thaliana.representative_proteins"
        for row in reader:
            for token in row[column].split(","):
                gene = token.strip().split("|", 1)[0]
                if gene:
                    ath_gene_to_og[gene] = row["Orthogroup"]

    adata = ad.read_h5ad(args.h5ad, backed="r")
    feature_rows = []
    mapped_features_by_og = defaultdict(list)
    status_counts = Counter()
    for position, (reference_gene, row) in enumerate(adata.var.iterrows()):
        reference_gene = str(reference_gene)
        source_gene = str(row.get("gene_ids", ""))
        orthogroup = ath_gene_to_og.get(reference_gene, "")
        status = "exact_arabidopsis_reference_gene_to_project_orthogroup" if orthogroup else "unmapped_reference_gene"
        status_counts[status] += 1
        feature_rows.append(
            {
                "feature_position": position,
                "moricandia_source_gene_id": source_gene,
                "arabidopsis_reference_gene_id": reference_gene,
                "project_orthogroup_id": orthogroup,
                "mapping_status": status,
            }
        )
        if orthogroup:
            mapped_features_by_og[orthogroup].append(source_gene)

    output_path = Path(args.output_map)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(feature_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(feature_rows)

    programs = defaultdict(set)
    with open(args.candidates, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            programs[(row["panel"], row["coarse_label"])].add(row["orthogroup_id"])
    candidate_ogs = set().union(*programs.values())
    mapped_ogs = set(mapped_features_by_og)
    audit = {
        "inputs": {
            "h5ad": args.h5ad,
            "h5ad_sha256": sha256(args.h5ad),
            "project_orthogroups": args.project_orthogroups,
            "project_orthogroups_sha256": sha256(args.project_orthogroups),
            "candidates": args.candidates,
            "candidates_sha256": sha256(args.candidates),
        },
        "feature_mapping": {
            "total_features": len(feature_rows),
            "unique_arabidopsis_reference_gene_ids": len(set(map(str, adata.var_names))),
            "duplicate_reference_gene_rows": len(feature_rows) - len(set(map(str, adata.var_names))),
            "status_counts": dict(sorted(status_counts.items())),
            "mapped_project_orthogroups": len(mapped_ogs),
        },
        "candidate_program_coverage": {
            "candidate_orthogroups": len(candidate_ogs),
            "covered_candidate_orthogroups": len(candidate_ogs & mapped_ogs),
            "coverage_fraction": len(candidate_ogs & mapped_ogs) / len(candidate_ogs),
            "by_program": {
                f"{panel}:{label}": {
                    "candidate_orthogroups": len(orthogroups),
                    "covered_candidate_orthogroups": len(orthogroups & mapped_ogs),
                    "coverage_fraction": len(orthogroups & mapped_ogs) / len(orthogroups),
                    "mapped_moricandia_features": sum(len(mapped_features_by_og[og]) for og in orthogroups),
                }
                for (panel, label), orthogroups in sorted(programs.items())
            },
        },
        "rule": {
            "description": "Use the author-provided Arabidopsis reference gene stored as the H5AD var index and map it exactly into the frozen project orthogroup namespace; preserve duplicate Moricandia contigs as family members.",
            "expression_blind": True,
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
