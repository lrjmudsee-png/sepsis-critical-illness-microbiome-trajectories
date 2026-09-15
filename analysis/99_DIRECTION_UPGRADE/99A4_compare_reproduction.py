#!/usr/bin/env python3
"""Compare deterministic A-stage outputs from two independent fresh runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


DETERMINISTIC_OUTPUTS = [
    "A03_supplement_sheet_index.csv",
    "A03B_supplement_table2_sample_ids.csv",
    "A04_patient_time_join_audit.csv",
    "A05_pair_eligibility.csv",
    "A06_population_counts.csv",
    "A07_taxonomy_mapping_proposal.csv",
    "A08_feature_and_zero_coverage.csv",
    "A09_healthy_reference_feasibility.csv",
    "A10_host_data_availability.csv",
    "A11_existing_analysis_overlap.md",
    "A12_feasibility_decision.md",
    "A13_questions_for_scientific_review.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    reference = Path(args.reference).resolve()
    candidate = Path(args.candidate).resolve()
    output = Path(args.output).resolve()
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite comparison output: {output}")

    rows = []
    all_equal = True
    for relative in DETERMINISTIC_OUTPUTS:
        left = reference / relative
        right = candidate / relative
        left_hash = sha256(left) if left.exists() else ""
        right_hash = sha256(right) if right.exists() else ""
        equal = left.exists() and right.exists() and left_hash == right_hash
        all_equal = all_equal and equal
        rows.append(
            {
                "relative_path": relative,
                "reference_sha256": left_hash,
                "reproduction_sha256": right_hash,
                "byte_identical": equal,
                "status": "PASS" if equal else "FAIL",
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Deterministic outputs identical: {all_equal} ({sum(row['byte_identical'] for row in rows)}/{len(rows)})")
    return 0 if all_equal else 2


if __name__ == "__main__":
    raise SystemExit(main())

