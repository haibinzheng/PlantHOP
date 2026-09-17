#!/usr/bin/env python3
"""Rebuild priority-family trees from named sequences in a new immutable run.

Uses MAFFT and IQ-TREE with SH-aLRT and ultrafast bootstrap support. The script
refuses to overwrite either partial or final outputs and never touches the
original OrthoFinder v1 run.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path("/workspace/projects/PhyloOpenCell")
SOURCE = Path("/data/runs/PhyloOpenCell/orthofinder_full7_v1/results/Results_Sep10/Orthogroup_Sequences")
RUN_ROOT = Path("/data/runs/PhyloOpenCell/priority_gene_trees_v2")
PARTIAL = RUN_ROOT.with_name(RUN_ROOT.name + ".partial")
MAFFT = Path("/workspace/venvs/phyloopencell/bin/mafft")
IQTREE = Path("/workspace/venvs/phyloopencell/bin/iqtree3")
SEED = 20260911


def now():
    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_headers(path: Path):
    headers = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].split()[0])
    return headers


def run_version(command):
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return (result.stdout + result.stderr).strip().splitlines()[0]


def main():
    if RUN_ROOT.exists() or PARTIAL.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {RUN_ROOT} or {PARTIAL}")
    if not MAFFT.exists() or not IQTREE.exists():
        raise RuntimeError("MAFFT or IQ-TREE executable is unavailable")

    with (PROJECT / "metadata/priority_paralog_gene_tree_evidence_v2.tsv").open("r", encoding="utf-8", newline="") as handle:
        targets = sorted({row["orthogroup_id"] for row in csv.DictReader(handle, delimiter="\t")})

    (PARTIAL / "sequences").mkdir(parents=True)
    (PARTIAL / "alignments").mkdir()
    (PARTIAL / "trees").mkdir()
    (PARTIAL / "logs").mkdir()
    manifest = []
    status_path = PROJECT / "reports/priority_gene_trees_v2_status.json"
    status = {
        "started_utc": now(),
        "state": "running",
        "run_root": str(RUN_ROOT),
        "source": str(SOURCE),
        "orthogroups": len(targets),
        "seed": SEED,
        "mafft_version": run_version([str(MAFFT), "--version"]),
        "iqtree_version": run_version([str(IQTREE), "--version"]),
        "source_h5ad_accessed": False,
        "model_training_started": False,
        "original_v1_modified": False,
    }
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

    try:
        for index, og in enumerate(targets, start=1):
            source = SOURCE / f"{og}.fa"
            if not source.exists():
                raise FileNotFoundError(source)
            sequence_copy = PARTIAL / "sequences" / source.name
            shutil.copy2(source, sequence_copy)
            headers = fasta_headers(sequence_copy)
            if len(headers) != len(set(headers)):
                raise RuntimeError(f"duplicate sequence header in {og}")

            alignment = PARTIAL / "alignments" / f"{og}.mafft.fa"
            mafft_log = PARTIAL / "logs" / f"{og}.mafft.log"
            with alignment.open("w", encoding="utf-8") as out, mafft_log.open("w", encoding="utf-8") as err:
                subprocess.run([str(MAFFT), "--auto", "--thread", "4", str(sequence_copy)], stdout=out, stderr=err, text=True, check=True)
            if set(fasta_headers(alignment)) != set(headers):
                raise RuntimeError(f"alignment membership mismatch for {og}")

            prefix = PARTIAL / "trees" / og
            iq_log = PARTIAL / "logs" / f"{og}.iqtree.stdout.log"
            command = [
                str(IQTREE), "-s", str(alignment), "-m", "MFP", "-B", "1000",
                "--alrt", "1000", "-T", "4", "-seed", str(SEED),
                "--prefix", str(prefix), "--quiet",
            ]
            with iq_log.open("w", encoding="utf-8") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=True)
            treefile = Path(str(prefix) + ".treefile")
            if not treefile.exists() or treefile.stat().st_size == 0:
                raise RuntimeError(f"missing IQ-TREE treefile for {og}")
            tree_text = treefile.read_text(encoding="utf-8")
            tree_headers = set()
            for header in headers:
                if header in tree_text:
                    tree_headers.add(header)
            if tree_headers != set(headers):
                raise RuntimeError(f"tree membership mismatch for {og}")

            manifest.append({
                "orthogroup_id": og,
                "sequence_count": len(headers),
                "source_sha256": file_sha256(sequence_copy),
                "alignment_sha256": file_sha256(alignment),
                "tree_sha256": file_sha256(treefile),
                "tree_membership_match": True,
                "model_selection": "MFP",
                "ultrafast_bootstrap_replicates": 1000,
                "sh_alrt_replicates": 1000,
                "seed": SEED,
            })
            status.update({"completed_orthogroups": index, "last_completed": og, "updated_utc": now()})
            status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

        manifest_path = PARTIAL / "manifest.tsv"
        with manifest_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest)
        run_audit = {
            **status,
            "state": "complete",
            "completed_utc": now(),
            "completed_orthogroups": len(manifest),
            "manifest_sha256": file_sha256(manifest_path),
            "claim_status": "clean_supported_trees_built_topology_interpretation_pending",
        }
        (PARTIAL / "audit.json").write_text(json.dumps(run_audit, indent=2) + "\n", encoding="utf-8")
        os.replace(PARTIAL, RUN_ROOT)
        status_path.write_text(json.dumps(run_audit, indent=2) + "\n", encoding="utf-8")
        (PROJECT / "reports/priority_gene_trees_v2_audit.json").write_text(json.dumps(run_audit, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(run_audit, indent=2))
    except Exception as exc:
        status.update({"state": "failed", "failed_utc": now(), "error": repr(exc), "partial_retained": PARTIAL.exists()})
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    main()
