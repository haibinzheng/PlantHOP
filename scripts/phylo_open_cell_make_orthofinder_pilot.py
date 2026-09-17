#!/usr/bin/env python3
"""Create deterministic, atomic small OrthoFinder pilot inputs from candidate FASTAs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def records(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = None
        sequence = []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line.rstrip("\n"), []
            else:
                sequence.append(line.strip())
        if header is not None:
            yield header, "".join(sequence)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--per-species", type=int, default=200)
    p.add_argument("--manifest", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in sorted(args.source_dir.glob("*.fa.gz")):
        target = args.output_dir / source.name.removesuffix(".gz")
        partial = target.with_suffix(target.suffix + ".partial")
        count = 0
        with partial.open("w", encoding="utf-8", newline="\n") as out:
            for header, sequence in records(source):
                if count >= args.per_species:
                    break
                if not sequence:
                    continue
                out.write(header + "\n")
                for start in range(0, len(sequence), 80):
                    out.write(sequence[start:start + 80] + "\n")
                count += 1
        if count != args.per_species:
            partial.unlink(missing_ok=True)
            raise RuntimeError(f"{source}: expected {args.per_species}, got {count}")
        os.replace(partial, target)
        rows.append({
            "species_file": target.name,
            "records": count,
            "source_path": str(source),
            "source_sha256": digest(source),
            "pilot_sha256": digest(target),
            "selection": "first records from gene-sorted candidate FASTA",
            "accepted": "False",
        })
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    partial_manifest = args.manifest.with_suffix(args.manifest.suffix + ".partial")
    with partial_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(partial_manifest, args.manifest)
    print(f"species={len(rows)} records={sum(r['records'] for r in rows)}")


if __name__ == "__main__":
    main()
