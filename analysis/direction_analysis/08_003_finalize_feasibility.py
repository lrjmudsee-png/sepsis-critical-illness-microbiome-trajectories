#!/usr/bin/env python3
"""Finalize A-stage integrity guards, QA, and output manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_PRE_FINAL = [
    "A01_input_inventory.csv",
    "A02_source_retrieval_log.csv",
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
    "session_info.txt",
    "python_source_inspection_session.txt",
    "_input_hashes_before.csv",
    "_protected_submission_tree_before.csv",
    "sources/PRJNA851469_Supplementary_Tables_2_17_OFFICIAL.xlsx",
]
EXPECTED_ELIGIBLE_PAIRS = {
    "PRJNA691455": 9,
    "PRJNA851469": 14,
    "PRJNA516701": 15,
    "PRJNA1125274": 24,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv_once(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite existing audit output: {path}")
    if not rows:
        raise RuntimeError(f"No rows for required output: {path}")
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text_once(path: Path, text: str) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite existing audit output: {path}")
    path.write_text(text, encoding="utf-8")


def bool_value(value: str) -> bool:
    return value.strip().upper() in {"TRUE", "T", "1", "YES", "Y"}


def finalize_input_hashes(output: Path) -> tuple[list[dict[str, object]], bool]:
    before_rows = read_csv(output / "_input_hashes_before.csv")
    final_rows: list[dict[str, object]] = []
    all_unchanged = True
    for row in before_rows:
        path = Path(row["path"])
        exists_after = path.exists()
        sha_after = sha256_file(path) if exists_after else ""
        bytes_after = path.stat().st_size if exists_after else ""
        unchanged = (
            bool_value(row["exists_before"]) == exists_after
            and row["sha256_before"] == sha_after
            and str(row["bytes_before"]) == str(bytes_after)
        )
        all_unchanged = all_unchanged and unchanged
        final_rows.append(
            {
                "guard_type": "named_input",
                "path": str(path),
                "exists_before": row["exists_before"],
                "exists_after": exists_after,
                "bytes_before": row["bytes_before"],
                "bytes_after": bytes_after,
                "sha256_before": row["sha256_before"],
                "sha256_after": sha_after,
                "status": "UNCHANGED" if unchanged else "CHANGED_OR_MISSING",
            }
        )

    tree_rows = read_csv(output / "_protected_submission_tree_before.csv")
    if tree_rows and tree_rows[0]["relative_path"]:
        root = Path(tree_rows[0]["protected_root"])
        before_by_relative = {row["relative_path"]: row for row in tree_rows}
        current_files = {
            str(path.relative_to(root)): path
            for path in sorted(item for item in root.rglob("*") if item.is_file())
        } if root.exists() else {}
        for relative in sorted(set(before_by_relative) | set(current_files)):
            before = before_by_relative.get(relative)
            current = current_files.get(relative)
            sha_before = before["sha256_before"] if before else ""
            sha_after = sha256_file(current) if current else ""
            bytes_before = before["bytes_before"] if before else ""
            bytes_after = current.stat().st_size if current else ""
            unchanged = before is not None and current is not None and sha_before == sha_after and str(bytes_before) == str(bytes_after)
            all_unchanged = all_unchanged and unchanged
            final_rows.append(
                {
                    "guard_type": "protected_submission_tree",
                    "path": str(root / relative),
                    "exists_before": before is not None,
                    "exists_after": current is not None,
                    "bytes_before": bytes_before,
                    "bytes_after": bytes_after,
                    "sha256_before": sha_before,
                    "sha256_after": sha_after,
                    "status": "UNCHANGED" if unchanged else "CHANGED_OR_MISSING",
                }
            )
    write_csv_once(output / "input_hashes_before_after.csv", final_rows)
    return final_rows, all_unchanged


def run_qa(output: Path, guards_unchanged: bool) -> list[dict[str, object]]:
    qa: list[dict[str, object]] = []
    for relative in REQUIRED_PRE_FINAL:
        path = output / relative
        passed = path.exists() and path.stat().st_size > 0
        qa.append(
            {
                "check": f"required_output:{relative}",
                "status": "PASS" if passed else "FAIL",
                "observed": path.stat().st_size if path.exists() else "missing",
                "expected": "exists_and_nonempty",
            }
        )

    pair_rows = read_csv(output / "A05_pair_eligibility.csv")
    observed_pairs: dict[str, int] = {}
    for project in EXPECTED_ELIGIBLE_PAIRS:
        observed = sum(
            1
            for row in pair_rows
            if row["project"] == project and bool_value(row["A_stage_eligible"])
        )
        observed_pairs[project] = observed
        qa.append(
            {
                "check": f"frozen_pair_count:{project}",
                "status": "PASS" if observed == EXPECTED_ELIGIBLE_PAIRS[project] else "FAIL",
                "observed": observed,
                "expected": EXPECTED_ELIGIBLE_PAIRS[project],
            }
        )

    sample_rows = read_csv(output / "A04_patient_time_join_audit.csv")
    nonexact = sum(row["audit_status"] != "EXACT_SAMPLE_COUNT_JOIN" for row in sample_rows)
    qa.append(
        {
            "check": "sample_to_count_exact_join",
            "status": "PASS" if nonexact == 0 else "FAIL",
            "observed": nonexact,
            "expected": 0,
        }
    )

    coverage_rows = read_csv(output / "A08_feature_and_zero_coverage.csv")
    unmatched = sum(int(float(row["n_feature_keys_unmatched"])) for row in coverage_rows)
    duplicate_keys = sum(int(float(row["n_duplicated_taxonomy_keys"])) for row in coverage_rows)
    no_candidate_group = [
        row["project"]
        for row in coverage_rows
        if int(float(row["n_P_features"])) == 0 or int(float(row["n_C_features"])) == 0
    ]
    qa.extend(
        [
            {
                "check": "taxonomy_feature_key_unmatched_total",
                "status": "PASS" if unmatched == 0 else "REVIEW",
                "observed": unmatched,
                "expected": 0,
            },
            {
                "check": "duplicated_taxonomy_key_total",
                "status": "PASS" if duplicate_keys == 0 else "REVIEW",
                "observed": duplicate_keys,
                "expected": 0,
            },
            {
                "check": "P_and_C_features_present_all_cohorts",
                "status": "PASS" if not no_candidate_group else "FAIL",
                "observed": "none_missing" if not no_candidate_group else "|".join(no_candidate_group),
                "expected": "P>0_and_C>0_each_cohort",
            },
            {
                "check": "protected_inputs_and_submission_tree_unchanged",
                "status": "PASS" if guards_unchanged else "FAIL",
                "observed": guards_unchanged,
                "expected": True,
            },
        ]
    )
    for row in coverage_rows:
        total_reads = float(row["P_total_reads"]) + float(row["C_total_reads"]) + float(row["OTHER_total_reads"])
        fractions_valid = (
            0 <= float(row["P_read_fraction"]) <= 1
            and 0 <= float(row["C_read_fraction"]) <= 1
            and 0 <= float(row["unresolved_family_read_fraction"]) <= 1
            and total_reads > 0
        )
        qa.append(
            {
                "check": f"P_C_OTHER_partition:{row['project']}",
                "status": "PASS" if fractions_valid else "FAIL",
                "observed": total_reads,
                "expected": "positive_total_and_valid_P_C_fractions",
            }
        )

    mapping_rows = read_csv(output / "A07_taxonomy_mapping_proposal.csv")
    mapping_review = [row for row in mapping_rows if bool_value(row["needs_manual_review"])]
    unexpected_candidate_like = [
        row for row in mapping_review
        if row["normalized_family_label"] != "<UNRESOLVED>"
    ]
    qa.append(
        {
            "check": "taxonomy_mapping_manual_review_rows",
            "status": "FAIL" if unexpected_candidate_like else ("REVIEW" if mapping_review else "PASS"),
            "observed": (
                f"total={len(mapping_review)};unexpected_candidate_like={len(unexpected_candidate_like)};"
                + "labels=" + "|".join(sorted({row['normalized_family_label'] for row in mapping_review}))
            ),
            "expected": "unresolved buckets retained in OTHER; no unapproved candidate-like aliases",
        }
    )

    population_rows = read_csv(output / "A06_population_counts.csv")
    ext_all = [
        row for row in population_rows
        if row["project"] == "PRJNA1125274" and row["center"] == "ALL_HOSPITALS"
    ]
    ext_patient_count = int(float(ext_all[0]["n_patients_object"])) if len(ext_all) == 1 else -1
    qa.append(
        {
            "check": "PRJNA1125274_object_patient_alias_groups",
            "status": "PASS" if ext_patient_count == 134 else "FAIL",
            "observed": ext_patient_count,
            "expected": 134,
        }
    )

    healthy_rows = read_csv(output / "A09_healthy_reference_feasibility.csv")
    healthy_ok = (
        len(healthy_rows) == 1
        and int(float(healthy_rows[0]["frozen_object_healthy_like_samples"])) == 13
        and int(float(healthy_rows[0]["official_table2_healthy_like_sample_ids"])) == 15
        and int(float(healthy_rows[0]["official_table2_healthy_ids_not_in_frozen_object"])) == 2
    )
    qa.append(
        {
            "check": "PRJNA851469_healthy_reference_roster",
            "status": "PASS" if healthy_ok else "FAIL",
            "observed": (
                f"frozen={healthy_rows[0]['frozen_object_healthy_like_samples']};"
                f"official={healthy_rows[0]['official_table2_healthy_like_sample_ids']};"
                f"official_not_frozen={healthy_rows[0]['official_table2_healthy_ids_not_in_frozen_object']}"
                if healthy_rows else "missing"
            ),
            "expected": "frozen=13;official=15;official_not_frozen=2",
        }
    )

    workbook_rows = read_csv(output / "A02_source_retrieval_log.csv")
    official = [row for row in workbook_rows if row["source_type"] == "official_publisher_supplement"]
    workbook_ok = len(official) == 1 and bool_value(official[0]["container_crc_valid"])
    qa.append(
        {
            "check": "official_supplement_crc",
            "status": "PASS" if workbook_ok else "FAIL",
            "observed": official[0]["detail"] if official else "missing",
            "expected": "valid_xlsx_zip_crc",
        }
    )
    return qa


def write_manifest(output: Path) -> None:
    manifest_path = output / "output_manifest.csv"
    if manifest_path.exists():
        raise RuntimeError(f"Refusing to overwrite existing audit output: {manifest_path}")
    rows: list[dict[str, object]] = []
    for path in sorted(item for item in output.rglob("*") if item.is_file() and item != manifest_path):
        rows.append(
            {
                "relative_path": str(path.relative_to(output)).replace("\\", "/"),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_csv_once(manifest_path, rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()
    if not output.exists():
        raise RuntimeError(f"Output directory does not exist: {output}")

    missing = [relative for relative in REQUIRED_PRE_FINAL if not (output / relative).exists()]
    if missing:
        raise RuntimeError("Missing required pre-final outputs: " + ", ".join(missing))

    _, guards_unchanged = finalize_input_hashes(output)
    qa = run_qa(output, guards_unchanged)
    write_csv_once(output / "A14_quality_control_checks.csv", qa)
    failures = [row for row in qa if row["status"] == "FAIL"]
    review_items = [row for row in qa if row["status"] == "REVIEW"]
    summary = "\n".join(
        [
            "# A14 质量控制摘要",
            "",
            f"- 生成时间（UTC）：{datetime.now(timezone.utc).isoformat()}",
            f"- 强制检查：{len(qa)}项",
            f"- FAIL：{len(failures)}项",
            f"- REVIEW：{len(review_items)}项",
            f"- 输入与 V2-SUBMISSION 保护树未改变：{'是' if guards_unchanged else '否'}",
            "- A阶段范围：未计算新的效应量、置信区间或p值；未重跑DADA2/FASTQ/98D。",
            "",
            "## 失败项",
            "",
            *([f"- {row['check']}：observed={row['observed']}，expected={row['expected']}" for row in failures] or ["- 无。"]),
            "",
            "## 需复核项",
            "",
            *([f"- {row['check']}：observed={row['observed']}，expected={row['expected']}" for row in review_items] or ["- 无。"]),
            "",
        ]
    )
    write_text_once(output / "A14_quality_control_summary.md", summary)

    with (output / "session_info.txt").open("a", encoding="utf-8") as handle:
        handle.write(f"\nfinalizer_timestamp_utc={datetime.now(timezone.utc).isoformat()}\n")
        handle.write(f"finalizer_python={sys.version}\n")
        handle.write(f"finalizer_platform={platform.platform()}\n")

    write_manifest(output)
    if failures:
        print("A-stage finalization completed with QA failures:")
        for row in failures:
            print(f"- {row['check']}: observed={row['observed']} expected={row['expected']}")
        return 2
    print(f"A-stage finalization passed {len(qa)} checks; review items={len(review_items)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
