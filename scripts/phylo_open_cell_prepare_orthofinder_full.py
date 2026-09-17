#!/usr/bin/env python3
"""Materialize and validate full candidate FASTAs for an auditable OrthoFinder run."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


ALLOWED = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-?.")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source-dir", type=Path, required=True)
    p.add_argument("--run-root", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--readiness", type=Path, required=True)
    args = p.parse_args()

    input_dir = args.run_root / "input"
    partial_dir = args.run_root / "input.partial"
    if input_dir.exists() or partial_dir.exists():
        raise RuntimeError("refusing to overwrite existing full-run input or partial input")
    partial_dir.mkdir(parents=True)
    rows = []
    for source in sorted(args.source_dir.glob("*.fa.gz")):
        target = partial_dir / source.name.removesuffix(".gz")
        record_count = 0
        residues = 0
        empty_sequences = 0
        invalid_characters = set()
        duplicate_ids = 0
        identifiers = set()
        with gzip.open(source, "rt", encoding="utf-8") as inp, target.open("w", encoding="utf-8", newline="\n") as out:
            current_has_sequence = False
            for raw in inp:
                line = raw.rstrip("\r\n")
                if line.startswith(">"):
                    if record_count and not current_has_sequence:
                        empty_sequences += 1
                    identifier = line[1:].split()[0]
                    if identifier in identifiers:
                        duplicate_ids += 1
                    identifiers.add(identifier)
                    record_count += 1
                    current_has_sequence = False
                    out.write(line + "\n")
                else:
                    sequence = line.strip().upper()
                    if sequence:
                        current_has_sequence = True
                        residues += len(sequence)
                        invalid_characters.update(set(sequence) - ALLOWED)
                        out.write(sequence + "\n")
            if record_count and not current_has_sequence:
                empty_sequences += 1
        rows.append({
                "species_file": target.name,
                "records": record_count,
                "residues": residues,
                "empty_sequences": empty_sequences,
                "duplicate_sequence_ids": duplicate_ids,
                "invalid_characters": "".join(sorted(invalid_characters)),
                "source_path": str(source),
                "source_sha256": digest(source),
                "input_sha256": digest(target),
                "accepted": "False",
        })
    if len(rows) != 7:
        raise RuntimeError(f"expected 7 species FASTAs, found {len(rows)}")
    problems = [r for r in rows if r["records"] == 0 or r["empty_sequences"] or r["duplicate_sequence_ids"] or r["invalid_characters"]]
    if problems:
        raise RuntimeError(f"FASTA validation failed; partial input retained for inspection: {problems}")
    os.replace(partial_dir, input_dir)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest_tmp = args.manifest.with_suffix(args.manifest.suffix + ".partial")
    with manifest_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(manifest_tmp, args.manifest)

    readiness = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "ready": True,
        "accepted": False,
        "input_dir": str(input_dir),
        "species": len(rows),
        "proteins": sum(r["records"] for r in rows),
        "residues": sum(r["residues"] for r in rows),
        "manifest": str(args.manifest),
        "manifest_sha256": digest(args.manifest),
        "validation": {
            "empty_sequences": sum(r["empty_sequences"] for r in rows),
            "duplicate_sequence_ids_within_species": sum(r["duplicate_sequence_ids"] for r in rows),
            "invalid_characters": sorted(set().union(*(set(r["invalid_characters"]) for r in rows))),
        },
        "scope": "seven species with version-supported protein candidates; Zea mays remains excluded pending Zm00001d.1 proteins",
        "source_h5ad_modified": False,
    }
    args.readiness.parent.mkdir(parents=True, exist_ok=True)
    readiness_tmp = args.readiness.with_suffix(args.readiness.suffix + ".partial")
    readiness_tmp.write_text(json.dumps(readiness, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(readiness_tmp, args.readiness)
    print(json.dumps(readiness, ensure_ascii=False))


if __name__ == "__main__":
    main()
