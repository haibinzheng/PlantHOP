#!/usr/bin/env python3
"""Correct the runtime field in an existing OrthoFinder audit from GNU time output."""

from __future__ import annotations

import argparse
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--audit", type=Path, required=True)
    p.add_argument("--log", type=Path, required=True)
    args = p.parse_args()
    payload = json.loads(args.audit.read_text(encoding="utf-8"))
    clean_log = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", args.log.read_text(encoding="utf-8", errors="replace"))
    elapsed = re.search(r"Elapsed \(wall clock\) time.*\):\s*(\S+)", clean_log)
    if not elapsed:
        raise RuntimeError("elapsed wall-clock line not found")
    previous = payload["runtime"].get("elapsed_wall_clock")
    payload["runtime"]["elapsed_wall_clock"] = elapsed.group(1)
    payload.setdefault("audit_corrections", []).append({
        "corrected_utc": datetime.now(timezone.utc).isoformat(),
        "field": "runtime.elapsed_wall_clock",
        "previous": previous,
        "corrected": elapsed.group(1),
        "reason": "initial parser captured the time-format label rather than the value",
    })
    partial = args.audit.with_suffix(args.audit.suffix + ".partial")
    partial.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, args.audit)
    print(elapsed.group(1))


if __name__ == "__main__":
    main()
