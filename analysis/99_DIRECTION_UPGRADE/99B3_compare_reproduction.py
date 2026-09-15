#!/usr/bin/env python3
"""Compare deterministic stage-B artifacts between primary and clean rerun outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary", required=True)
    parser.add_argument("--reproduction", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    primary = Path(args.primary).resolve()
    reproduction = Path(args.reproduction).resolve()
    output = Path(args.output).resolve()

    excluded = {
        "_input_hashes_before.csv",
        "_protected_submission_tree_before.csv",
        "input_hashes_before_after.csv",
        "protected_submission_tree_before_after.csv",
        "session_info.txt",
        "output_manifest.csv",
    }
    candidates = sorted(
        path.relative_to(primary)
        for path in primary.rglob("*")
        if path.is_file()
        and path.name not in excluded
        and path.suffix.lower() != ".pdf"
        and not path.name.startswith("_RUN_")
    )
    rows: list[dict[str, object]] = []
    for relative in candidates:
        left = primary / relative
        right = reproduction / relative
        exists = right.is_file()
        left_hash = digest(left)
        right_hash = digest(right) if exists else ""
        rows.append(
            {
                "relative_path": str(relative),
                "primary_sha256": left_hash,
                "reproduction_sha256": right_hash,
                "byte_identical": exists and left_hash == right_hash,
            }
        )
    if not rows:
        raise RuntimeError("No deterministic artifacts found")
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite: {output}")
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    failed = [row for row in rows if not row["byte_identical"]]
    print(f"Compared {len(rows)} deterministic artifacts; mismatches={len(failed)}")
    return 2 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
