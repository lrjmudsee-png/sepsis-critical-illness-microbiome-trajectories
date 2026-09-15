#!/usr/bin/env python3
"""Initialize stage-B provenance guards in a fresh output directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"E:\sepsis_project")
INPUTS = [
    ROOT / r"code\03_data_processing\99_DIRECTION_UPGRADE\B00_FIXED_SECONDARY_ANALYSIS_SPECIFICATION.md",
    ROOT / r"results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY\A04_patient_time_join_audit.csv",
    ROOT / r"results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY\A05_pair_eligibility.csv",
    ROOT / r"results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY\A07_taxonomy_mapping_proposal.csv",
    ROOT / r"results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY\A08_feature_and_zero_coverage.csv",
    ROOT / r"results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY\A09_healthy_reference_feasibility.csv",
    ROOT / r"results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY\A12_feasibility_decision.md",
    ROOT / r"data\_V2_ANALYSIS_READY\01_PROJECTS\PRJNA691455\05_work\step87B_analysis_object\PRJNA691455_analysis_object_step87B.rds",
    ROOT / r"data\_V2_ANALYSIS_READY\01_PROJECTS\PRJNA851469\05_work\step87B_analysis_object\PRJNA851469_analysis_object_step87B.rds",
    ROOT / r"data\_V2_ANALYSIS_READY\01_PROJECTS\PRJNA516701\05_work\step87B_analysis_object\PRJNA516701_analysis_object_step87B.rds",
    ROOT / r"data\PRJNA1125274\03_dada2\PRJNA1125274_analysis_object_external_validation.rds",
    ROOT / r"results\V2_UPGRADE_20260907\98A_AITCHISON_LONGITUDINAL_ROBUSTNESS\01_SAMPLE_DISPLACEMENT\PRJNA691455_within_patient_displacement_bray_aitchison.csv",
    ROOT / r"results\V2_UPGRADE_20260907\98A_AITCHISON_LONGITUDINAL_ROBUSTNESS\01_SAMPLE_DISPLACEMENT\PRJNA851469_within_patient_displacement_bray_aitchison.csv",
    ROOT / r"results\V2_UPGRADE_20260907\98A_AITCHISON_LONGITUDINAL_ROBUSTNESS\01_SAMPLE_DISPLACEMENT\PRJNA516701_within_patient_displacement_bray_aitchison.csv",
    ROOT / r"results\V2_UPGRADE_20260907\98E_PRJNA1125274_EXTERNAL_VALIDATION\PRJNA1125274_external_validation_sample_displacement.csv",
    ROOT / r"results\V2_UPGRADE_20260910_FINAL_FREEZE\PRJNA1125274_FINAL_population_membership.csv",
    ROOT / r"results\V2_36G2_PRJNA851469_DAY3_CLINICAL_OUTCOME_MODELS_FIXED\10_MANUSCRIPT_SAFE_SUMMARY.txt",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    missing = [str(path) for path in INPUTS if not path.exists()]
    if missing:
        raise RuntimeError("Missing required stage-B inputs: " + " | ".join(missing))

    rows = []
    for path in INPUTS:
        rows.append(
            {
                "path": str(path),
                "exists_before": True,
                "bytes_before": path.stat().st_size,
                "mtime_utc_before": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
                "sha256_before": sha256_file(path),
            }
        )
    write_csv(output / "_input_hashes_before.csv", rows)

    submission = ROOT / "V2-SUBMISSION"
    submission_rows = []
    for path in sorted(item for item in submission.rglob("*") if item.is_file()):
        submission_rows.append(
            {
                "protected_root": str(submission),
                "relative_path": str(path.relative_to(submission)),
                "bytes_before": path.stat().st_size,
                "sha256_before": sha256_file(path),
            }
        )
    if not submission_rows:
        raise RuntimeError("Protected V2-SUBMISSION tree is missing or empty")
    write_csv(output / "_protected_submission_tree_before.csv", submission_rows)
    print(f"Stage-B provenance initialized with {len(rows)} named inputs and {len(submission_rows)} protected submission files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

