#!/usr/bin/env python3
"""Run the full OrthoFinder candidate analysis with status checkpoints and atomic promotion."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_status(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, path)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--run-root", type=Path, required=True)
    p.add_argument("--environment-bin", type=Path, required=True)
    p.add_argument("--threads", type=int, default=12)
    p.add_argument("--status", type=Path, required=True)
    p.add_argument("--log", type=Path, required=True)
    args = p.parse_args()

    final_dir = args.run_root / "results"
    partial_dir = args.run_root / "results.partial"
    if final_dir.exists() or partial_dir.exists():
        raise RuntimeError("refusing to overwrite existing final or partial results")
    args.run_root.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    started = now()
    base_status = {
        "run": "orthofinder_full7_v1",
        "started_utc": started,
        "wrapper_pid": os.getpid(),
        "input_dir": str(args.input_dir),
        "partial_results": str(partial_dir),
        "final_results": str(final_dir),
        "threads": args.threads,
        "accepted": False,
        "source_h5ad_modified": False,
    }
    env = os.environ.copy()
    env["PATH"] = str(args.environment_bin) + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = str(args.threads)
    command = [
        "/usr/bin/time", "-v",
        str(args.environment_bin / "orthofinder"),
        "-f", str(args.input_dir),
        "-o", str(partial_dir),
        "-t", str(args.threads),
        "-a", str(args.threads),
    ]
    os.nice(10)
    with args.log.open("w", encoding="utf-8", buffering=1) as log_handle:
        process = subprocess.Popen(command, stdout=log_handle, stderr=subprocess.STDOUT, env=env, text=True)
        write_status(args.status, {**base_status, "state": "running", "child_pid": process.pid, "command": command})
        return_code = process.wait()
    if return_code == 0:
        os.replace(partial_dir, final_dir)
        write_status(args.status, {**base_status, "state": "complete", "completed_utc": now(), "return_code": 0, "log": str(args.log)})
    else:
        write_status(args.status, {**base_status, "state": "failed", "failed_utc": now(), "return_code": return_code, "log": str(args.log), "partial_retained": partial_dir.exists()})
        raise SystemExit(return_code)


if __name__ == "__main__":
    main()
