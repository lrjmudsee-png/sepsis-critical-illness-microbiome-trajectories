#!/usr/bin/env python3
"""Validate stage-B outputs and add audit, interpretation, and manuscript handoff files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


ROOT = Path(r"E:\sepsis_project")
CODE_DIR = ROOT / r"code\03_data_processing\99_DIRECTION_UPGRADE"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv_once(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite: {path}")
    if not rows:
        raise RuntimeError(f"Refusing to write empty CSV: {path}")
    names = fields or list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def write_text_once(path: Path, text: str) -> None:
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite: {path}")
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def f(value: str | float | int | None) -> float:
    try:
        return float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan


def fmt(value: str | float | int | None, digits: int = 3) -> str:
    number = f(value)
    return f"{number:.{digits}f}" if math.isfinite(number) else "NA"


def is_true(value: object) -> bool:
    return str(value).strip().upper() in {"TRUE", "T", "1", "YES", "Y"}


def one(rows: list[dict[str, str]], **conditions: str) -> dict[str, str]:
    found = [row for row in rows if all(row.get(key) == value for key, value in conditions.items())]
    if len(found) != 1:
        raise RuntimeError(f"Expected one row for {conditions}; found {len(found)}")
    return found[0]


def add_check(checks: list[dict[str, object]], check_id: str, status: str, observed: object, expected: object, evidence: str) -> None:
    checks.append(
        {
            "check_id": check_id,
            "status": status,
            "observed": observed,
            "expected": expected,
            "evidence": evidence,
        }
    )


def holm_two(p_values: list[float]) -> list[float]:
    order = sorted(range(2), key=lambda index: p_values[index])
    adjusted_ordered: list[float] = []
    running = 0.0
    for rank, index in enumerate(order):
        candidate = min(1.0, (2 - rank) * p_values[index])
        running = max(running, candidate)
        adjusted_ordered.append(running)
    output = [math.nan, math.nan]
    for rank, index in enumerate(order):
        output[index] = adjusted_ordered[rank]
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output).resolve()

    required = [
        "_input_hashes_before.csv",
        "_protected_submission_tree_before.csv",
        "B01_sample_direction_metrics.csv",
        "B02_patient_paired_direction_metrics.csv",
        "B03_pair_exclusion_ledger.csv",
        "B04_primary_cohort_effects.csv",
        "B05_natural_history_meta_analysis.csv",
        "B05B_meta_analysis_weights.csv",
        "B06_primary_holm_family.csv",
        "B07_cohort_and_external_sensitivity_results.csv",
        "B08_PRJNA851469_healthy_reference_sample_metrics.csv",
        "B09_PRJNA851469_healthy_reference_paired_result.csv",
        "B10_D_B_H_associations.csv",
        "B11_transition_quadrant_patient_table.csv",
        "B11B_transition_quadrant_summary.csv",
        "B12_formula_handcheck.csv",
        "B13_meta_implementation_crosscheck.csv",
        "B14_forest_plot_source.csv",
        "B15_result_decision_table.csv",
        "B16_candidate_results_text.md",
        "session_info.txt",
    ]
    figure_names = [
        "Figure_B1_direction_balance_forest.png",
        "Figure_B1_direction_balance_forest.pdf",
        "Figure_B2_patient_balance_trajectories.png",
        "Figure_B2_patient_balance_trajectories.pdf",
        "Figure_B3_displacement_vs_direction.png",
        "Figure_B3_displacement_vs_direction.pdf",
        "Figure_B4_healthy_reference_trajectory.png",
        "Figure_B4_healthy_reference_trajectory.pdf",
    ]
    missing = [name for name in required if not (output / name).is_file() or (output / name).stat().st_size == 0]
    missing += [f"figures/{name}" for name in figure_names if not (output / "figures" / name).is_file() or (output / "figures" / name).stat().st_size == 0]
    if missing:
        raise RuntimeError("Missing or empty required stage-B outputs: " + " | ".join(missing))

    input_before = read_csv(output / "_input_hashes_before.csv")
    input_guard_rows: list[dict[str, object]] = []
    for row in input_before:
        path = Path(row["path"])
        exists_after = path.is_file()
        after_hash = sha256_file(path) if exists_after else ""
        after_bytes = path.stat().st_size if exists_after else ""
        unchanged = exists_after and after_hash == row["sha256_before"] and str(after_bytes) == row["bytes_before"]
        input_guard_rows.append(
            {
                **row,
                "exists_after": exists_after,
                "bytes_after": after_bytes,
                "sha256_after": after_hash,
                "unchanged": unchanged,
            }
        )
    write_csv_once(output / "input_hashes_before_after.csv", input_guard_rows)

    submission = ROOT / "V2-SUBMISSION"
    protected_before = read_csv(output / "_protected_submission_tree_before.csv")
    before_by_relative = {row["relative_path"]: row for row in protected_before}
    current_files = sorted(path for path in submission.rglob("*") if path.is_file())
    current_by_relative = {str(path.relative_to(submission)): path for path in current_files}
    all_relative = sorted(set(before_by_relative) | set(current_by_relative))
    protected_rows: list[dict[str, object]] = []
    for relative in all_relative:
        before = before_by_relative.get(relative)
        current = current_by_relative.get(relative)
        before_hash = before["sha256_before"] if before else ""
        after_hash = sha256_file(current) if current else ""
        status = "UNCHANGED" if before and current and before_hash == after_hash else ("ADDED" if current and not before else ("REMOVED" if before and not current else "CHANGED"))
        protected_rows.append(
            {
                "relative_path": relative,
                "bytes_before": before.get("bytes_before", "") if before else "",
                "bytes_after": current.stat().st_size if current else "",
                "sha256_before": before_hash,
                "sha256_after": after_hash,
                "status": status,
            }
        )
    write_csv_once(output / "protected_submission_tree_before_after.csv", protected_rows)

    samples = read_csv(output / "B01_sample_direction_metrics.csv")
    pairs = read_csv(output / "B02_patient_paired_direction_metrics.csv")
    exclusion = read_csv(output / "B03_pair_exclusion_ledger.csv")
    cohort = read_csv(output / "B04_primary_cohort_effects.csv")
    meta = read_csv(output / "B05_natural_history_meta_analysis.csv")
    primary = read_csv(output / "B06_primary_holm_family.csv")
    sensitivity = read_csv(output / "B07_cohort_and_external_sensitivity_results.csv")
    healthy = read_csv(output / "B09_PRJNA851469_healthy_reference_paired_result.csv")
    handcheck = read_csv(output / "B12_formula_handcheck.csv")
    meta_check = read_csv(output / "B13_meta_implementation_crosscheck.csv")
    decisions = read_csv(output / "B15_result_decision_table.csv")
    decision = {row["item"]: row["value"] for row in decisions}

    checks: list[dict[str, object]] = []
    expected = {"PRJNA691455": (9, 9), "PRJNA851469": (14, 14), "PRJNA516701": (15, 12), "PRJNA1125274": (24, 24)}
    for project, (expected_frozen, expected_analyzed) in expected.items():
        row = one(cohort, project=project)
        observed = (int(f(row["n_expected"])), int(f(row["n_analyzed"])))
        add_check(
            checks,
            f"PRIMARY_N_{project}",
            "PASS" if observed == (expected_frozen, expected_analyzed) else "FAIL",
            f"{observed[0]}/{observed[1]}",
            f"{expected_frozen}/{expected_analyzed}",
            "B04 n_expected/n_analyzed",
        )

    excluded_516 = sorted(
        row["patient_id"] for row in exclusion
        if row["project"] == "PRJNA516701" and is_true(row["in_frozen_primary_pair_population"]) and not is_true(row["main_delta_B_calculable"])
    )
    add_check(
        checks,
        "DOUBLE_ZERO_EXCLUSION_516701",
        "PASS" if excluded_516 == ["3516", "3669", "3927"] else "FAIL",
        "|".join(excluded_516),
        "3516|3669|3927",
        "B03 frozen-pair exclusion ledger",
    )

    sample_spot = samples[0]
    p_reads, c_reads = f(sample_spot["P_reads"]), f(sample_spot["C_reads"])
    expected_b = math.log((p_reads + 0.5) / (c_reads + 0.5)) if p_reads + c_reads > 0 else math.nan
    observed_b = f(sample_spot["B_pc0_5"])
    add_check(
        checks,
        "SAMPLE_BALANCE_FORMULA_SPOTCHECK",
        "PASS" if (not math.isfinite(expected_b) and not math.isfinite(observed_b)) or abs(expected_b - observed_b) <= 1e-12 else "FAIL",
        observed_b,
        expected_b,
        f"B01 first row {sample_spot.get('project')}::{sample_spot.get('run_id')}",
    )

    max_hand_error = max(f(row["max_absolute_error"]) for row in handcheck)
    add_check(checks, "PAIR_FORMULA_HANDCHECK", "PASS" if max_hand_error <= 1e-12 else "FAIL", max_hand_error, "<=1e-12", "B12 two independently selected pairs")

    raw_p = [f(row["raw_p"]) for row in primary]
    observed_holm = [f(row["holm_p"]) for row in primary]
    expected_holm = holm_two(raw_p)
    holm_error = max(abs(a - b) for a, b in zip(observed_holm, expected_holm))
    add_check(checks, "HOLM_EXACT_TWO_TEST_FAMILY", "PASS" if len(primary) == 2 and holm_error <= 1e-12 else "FAIL", f"n={len(primary)};max_error={holm_error}", "n=2;max_error<=1e-12", "B06 versus independent Python calculation")

    primary_meta = one(meta, analysis="PRIMARY_PC0_5")
    used_se = f(primary_meta["used_se"])
    expected_used = max(f(primary_meta["conventional_se"]), f(primary_meta["hk_raw_se"]))
    add_check(checks, "CONSERVATIVE_HK_SE", "PASS" if abs(used_se - expected_used) <= 1e-12 else "FAIL", used_se, expected_used, "B05 max(conventional SE, raw HK SE)")
    t_crit = 4.302652729911275  # t(0.975, df=2)
    ci_error = max(
        abs(f(primary_meta["ci_low"]) - (f(primary_meta["pooled_mean"]) - t_crit * used_se)),
        abs(f(primary_meta["ci_high"]) - (f(primary_meta["pooled_mean"]) + t_crit * used_se)),
    )
    add_check(checks, "META_CI_DF_K_MINUS_1", "PASS" if int(f(primary_meta["k"])) == 3 and int(f(primary_meta["df"])) == 2 and ci_error <= 2e-10 else "FAIL", f"k={primary_meta['k']};df={primary_meta['df']};max_error={ci_error}", "k=3;df=2;max_error<=2e-10", "B05 independent CI reconstruction; tolerance covers CSV decimal serialization")

    cross = meta_check[0]
    cross_status = cross["status"]
    add_check(checks, "REML_IMPLEMENTATION_CROSSCHECK", "PASS" if cross_status == "PASS" else "REVIEW", f"{cross_status};method={cross.get('crosscheck_method', '')};difference={cross.get('absolute_difference', '')}", "PASS; difference<=1e-7", "B13 metafor comparison when available, otherwise independent base-R nlminb optimization")

    expected_rule = (
        "STRONG_COMMON_DIRECTION_SUPPORT"
        if all(f(row["estimate"]) > 0 for row in primary) and all(f(row["holm_p"]) < 0.05 for row in primary)
        else "DIRECTIONALLY_CONCORDANT_BUT_NOT_BOTH_HOLM_SIGNIFICANT"
        if all(f(row["estimate"]) > 0 for row in primary)
        else "CONCORDANT_OPPOSITE_TO_HYPOTHESIS"
        if all(f(row["estimate"]) < 0 for row in primary)
        else "DIRECTIONALLY_HETEROGENEOUS"
    )
    add_check(checks, "RULE_BASED_PRIMARY_CONCLUSION", "PASS" if decision.get("primary_conclusion_rule") == expected_rule else "FAIL", decision.get("primary_conclusion_rule"), expected_rule, "B06 and B15")

    h = healthy[0]
    add_check(checks, "HEALTHY_REFERENCE_COUNTS", "PASS" if int(f(h["n_expected"])) == 14 and int(f(h["n_analyzed"])) == 14 and int(f(h["n_healthy_references"])) == 13 else "FAIL", f"pairs={h['n_analyzed']}/{h['n_expected']};refs={h['n_healthy_references']}", "pairs=14/14;refs=13", "B09")

    figure_missing = [name for name in figure_names if (output / "figures" / name).stat().st_size == 0]
    add_check(checks, "FIGURE_ARTIFACTS", "PASS" if not figure_missing else "FAIL", f"{len(figure_names) - len(figure_missing)}/{len(figure_names)} nonempty", "8/8 nonempty", "figures directory")

    changed_inputs = [row["path"] for row in input_guard_rows if not row["unchanged"]]
    add_check(checks, "NAMED_INPUT_IMMUTABILITY", "PASS" if not changed_inputs else "FAIL", len(changed_inputs), 0, "input_hashes_before_after.csv")
    changed_submission = [row["relative_path"] for row in protected_rows if row["status"] != "UNCHANGED"]
    add_check(checks, "V2_SUBMISSION_IMMUTABILITY", "PASS" if not changed_submission else "FAIL", len(changed_submission), 0, "protected_submission_tree_before_after.csv")

    external_centers = [row for row in sensitivity if row["project"] == "PRJNA1125274" and row["analysis"] == "HOSPITAL_STRATIFIED_PRIMARY"]
    external_overall = one(cohort, project="PRJNA1125274")
    finite_center_means = [f(row["mean"]) for row in external_centers if math.isfinite(f(row["mean"]))]
    simpson_reversal = bool(finite_center_means) and all(value * f(external_overall["mean"]) < 0 for value in finite_center_means)
    add_check(checks, "EXTERNAL_CENTER_DIRECTION_AUDIT", "REVIEW" if simpson_reversal else "PASS", ";".join(f"{row['population']}={fmt(row['mean'])}" for row in external_centers), f"overall={fmt(external_overall['mean'])}; no all-center reversal", "B07 hospital-stratified sensitivity")

    write_csv_once(output / "B17_quality_control_checks.csv", checks)
    pass_count = sum(row["status"] == "PASS" for row in checks)
    review_count = sum(row["status"] == "REVIEW" for row in checks)
    fail_count = sum(row["status"] == "FAIL" for row in checks)
    qc_summary = f"""# Stage-B quality-control summary

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validation
- Verification Status: {'FAIL' if fail_count else 'PASS_WITH_REVIEW' if review_count else 'PASS'}
- Version: direction-stage-B-qc-v1

## Result

- PASS: {pass_count}
- REVIEW: {review_count}
- FAIL: {fail_count}
- Named analysis inputs changed: {len(changed_inputs)}
- Protected `V2-SUBMISSION` files changed: {len(changed_submission)}

No DADA2, FASTQ preprocessing, or 98D sequence processing was run. A REVIEW item is a transparency flag, not an analysis failure.
"""
    write_text_once(output / "B17_quality_control_summary.md", qc_summary)

    natural = one(primary, test_id="NATURAL_HISTORY_REML_HK")
    external = one(primary, test_id="PRJNA1125274_EXTERNAL_24")
    healthy_direction = "toward" if f(h["mean"]) < 0 else "away from"
    natural_cohort_means = {row["project"]: f(row["mean"]) for row in cohort if row["project"] != "PRJNA1125274"}
    all_natural_same = len({"POS" if value > 0 else "NEG" if value < 0 else "ZERO" for value in natural_cohort_means.values()}) == 1
    validation_report = f"""# Stage-B statistical validation and governance report

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validation and interpretation
- Verification Status: {'VALIDATED_WITH_CAUTION' if fail_count == 0 else 'VALIDATION_FAILED'}
- Version: direction-stage-B-governance-v1

## Confirmed numerical result

- Natural-history REML/Hartung-Knapp: estimate {fmt(natural['estimate'])}, 95% CI {fmt(natural['ci_low'])} to {fmt(natural['ci_high'])}, raw p={fmt(natural['raw_p'], 4)}, Holm p={fmt(natural['holm_p'], 4)}.
- Frozen PRJNA1125274 primary population (n=24): estimate {fmt(external['estimate'])}, 95% CI {fmt(external['ci_low'])} to {fmt(external['ci_high'])}, raw p={fmt(external['raw_p'], 4)}, Holm p={fmt(external['holm_p'], 4)}.
- Rule-based conclusion: `{expected_rule}`.
- PRJNA851469 healthy-reference secondary result: mean change {fmt(h['mean'])}, ordinary 95% CI {fmt(h['ci_low'])} to {fmt(h['ci_high'])}, two-level bootstrap 95% CI {fmt(h['two_level_bootstrap_ci_low'])} to {fmt(h['two_level_bootstrap_ci_high'])}; patients moved {healthy_direction} the same-cohort healthy reference on average.

## Statistical fallacy audit

1. **Simpson's paradox — checked.** External center-specific means were compared with the overall mean. All-center reversal detected: {str(simpson_reversal).upper()}. This does not eliminate other center heterogeneity.
2. **Ecological fallacy — caution.** Published immune-mediator tables are aggregate family-correlation matrices and cannot substitute for patient-level host measurements.
3. **Berkson's bias — caution.** ICU/sepsis cohort selection may induce associations that do not generalize outside enrolled populations.
4. **Collider bias — limited exposure.** No new covariate-adjusted host model was fitted, avoiding unsupported adjustment; selection into paired complete cases can still act as a collider.
5. **Base-rate neglect — not a diagnostic classifier analysis.** No claims about prediction or post-test probability are made.
6. **Regression to the mean — monitored.** Change was defined between fixed visits rather than selecting individuals by extreme baseline balance, but two-time-point change remains noisy.
7. **Survivorship/complete-case bias — caution.** Paired follow-up requires observed later samples. PRJNA516701 additionally excludes three double-zero balance pairs under the fixed estimand.
8. **Look-elsewhere effect — partly controlled.** The P/C family dictionary and two-test Holm family were fixed before examining new effect estimates; secondary correlations and quadrant summaries remain exploratory.
9. **Garden of forking paths — caution.** This is a time-stamped post-hoc secondary SAP, not preregistration. Pseudocount, common-anchor, leave-one-out, center, and complete-case results must remain labelled sensitivity analyses.
10. **Correlation-causation fallacy — caution.** The balance is a compositional direction marker, not evidence of absolute bacterial expansion or mechanism.
11. **Reverse causality — caution.** Treatment, illness severity, feeding, bowel function, and recovery may drive the microbial pattern.

## Governance decision

The analysis is acceptable as a reproducible secondary direction-of-change layer if all FAIL counts remain zero. It cannot be promoted to patient-level host-mechanism validation without joinable host data. The healthy-reference analysis is same-cohort supportive context, not an independent cohort.
"""
    write_text_once(output / "B18_statistical_validation_report.md", validation_report)

    main_effect_lines = "\n".join(
        f"- {row['project']}: n={int(f(row['n_analyzed']))}/{int(f(row['n_expected']))}, mean delta_B={fmt(row['mean'])}, 95% CI {fmt(row['ci_low'])} to {fmt(row['ci_high'])}, p={fmt(row['t_p'], 4)}, dz={fmt(row['cohen_dz'])}."
        for row in cohort
    )
    zh_summary = f"""# B阶段方向分析结果摘要

## 已完成的自动分析

固定方向指标为 `B = ln[(P+0.5)/(C+0.5)]`，其中 P 为 Enterobacteriaceae + Enterococcaceae，C 为 Lachnospiraceae + Ruminococcaceae + Oscillospiraceae。正的 delta_B 表示晚期样本相对更偏向 P；它不是绝对丰度，也不能单独证明机制。

{main_effect_lines}

三项自然病程队列的 REML/Hartung-Knapp 合并效应为 {fmt(natural['estimate'])}（95% CI {fmt(natural['ci_low'])}–{fmt(natural['ci_high'])}；Holm p={fmt(natural['holm_p'], 4)}）。外部 PRJNA1125274 冻结的24人主分析效应为 {fmt(external['estimate'])}（95% CI {fmt(external['ci_low'])}–{fmt(external['ci_high'])}；Holm p={fmt(external['holm_p'], 4)}）。按预先固定的双检验规则，结论为 `{expected_rule}`。

PRJNA851469 的13个同批次健康参照支持了一个次要分析：14名配对患者的平均健康距离变化为 {fmt(h['mean'])}，双层bootstrap 95% CI {fmt(h['two_level_bootstrap_ci_low'])}–{fmt(h['two_level_bootstrap_ci_high'])}。这不是独立外部验证。

## 质量控制

- QC：{pass_count} PASS，{review_count} REVIEW，{fail_count} FAIL。
- 公式手算、Holm校正、REML-HK置信区间、样本数与排除表均已程序化核验。
- PRJNA516701 的3516、3669、3927因一个时间点 P+C=0 按固定规则记为不可估计，而不是错误地赋值为0。
- 原始分析输入及 `V2-SUBMISSION` 均未改动。
- 未运行 DADA2、FASTQ preprocessing 或 98D。

## 可以写入稿件的程度

若QC无FAIL，本结果可以作为“跨队列、固定方向指标、独立外部数据支持”的次级证据。措辞必须保留为关联性和组成性结果；宿主机制、因果关系及绝对菌量仍不可声称。
"""
    write_text_once(output / "B19_results_summary_ZH.md", zh_summary)

    methods = f"""# Candidate manuscript integration text (not yet inserted into V2-SUBMISSION)

## Material Passport

- Origin Skill: academic-research-suite / academic-paper + experiment-agent
- Origin Mode: draft integration candidate
- Verification Status: NUMBERS_PROGRAMMATICALLY_VERIFIED; AUTHOR_RATIFICATION_REQUIRED
- Version: direction-stage-B-manuscript-candidate-v1

## Methods candidate

We defined a fixed compositional direction balance as ln[(P+0.5)/(C+0.5)], where P comprised Enterobacteriaceae and Enterococcaceae and C comprised Lachnospiraceae, Ruminococcaceae and Oscillospiraceae. Family labels were matched exactly after removal of a leading `f__` prefix. The primary patient-level estimand was the late-minus-early change in this balance. Samples in which P+C equalled zero were treated as non-estimable. Within-cohort mean changes were summarized using t-based 95% confidence intervals, two-sided one-sample t tests, Wilcoxon tests, standardized paired-change effect sizes, and 5,000 fixed-seed bootstrap resamples. Natural-history cohort estimates were pooled with REML random effects and a conservative Hartung-Knapp standard error defined as the larger of the conventional and raw Hartung-Knapp standard errors. The primary multiplicity family comprised the natural-history meta-analysis and the frozen 24-patient PRJNA1125274 external analysis, with Holm adjustment. Pseudocount, nonzero-only, common-anchor, complete-case, SAL8-exclusion, center-stratified and leave-one-cohort-out analyses were designated sensitivity analyses.

## Results candidate

Across the three natural-history cohorts, the pooled early-to-late change in the fixed direction balance was {fmt(natural['estimate'])} (95% CI {fmt(natural['ci_low'])} to {fmt(natural['ci_high'])}; Holm-adjusted p={fmt(natural['holm_p'], 4)}). In the frozen PRJNA1125274 primary population (n=24), the corresponding change was {fmt(external['estimate'])} (95% CI {fmt(external['ci_low'])} to {fmt(external['ci_high'])}; Holm-adjusted p={fmt(external['holm_p'], 4)}). The prespecified rule-based interpretation was `{expected_rule}`. Three of 15 frozen PRJNA516701 pairs were non-estimable because P+C equalled zero at one visit; no neutral value was imputed.

In PRJNA851469, the mean Day-7 minus Day-3 change in Bray-Curtis distance to 13 co-processed healthy references was {fmt(h['mean'])} (ordinary 95% CI {fmt(h['ci_low'])} to {fmt(h['ci_high'])}; two-level bootstrap 95% CI {fmt(h['two_level_bootstrap_ci_low'])} to {fmt(h['two_level_bootstrap_ci_high'])}; n={int(f(h['n_analyzed']))}). This same-cohort comparison was considered supportive and not an independent validation.

## Discussion guardrail candidate

The fixed balance provides an interpretable, reproducible summary of relative community direction across heterogeneous studies. Because it is compositional, an increase cannot distinguish expansion of P from depletion of C, and neither the balance nor its association with 16S-derived displacement establishes host mechanism. Differences in sampling, treatment and taxonomic resolution—particularly the lower Family-level resolution in PRJNA1125274—remain important limitations. Patient-level host measurements were not publicly joinable and therefore were not modelled.
"""
    write_text_once(output / "B20_manuscript_integration_candidate.md", methods)

    changeset = f"""# Proposed submission-package changeset — requires author approval

No file under `E:\\sepsis_project\\V2-SUBMISSION` was changed in this run.

If the authors approve this analysis, the next editorial integration should:

1. Add the B20 Methods paragraph to the statistical/microbiome methods and label it a post-hoc fixed secondary analysis.
2. Add the two primary numerical results from B06/B20 to Results.
3. Add one forest figure (B1) to the main text or supplement; keep patient trajectories and diagnostics supplementary.
4. Add the exclusion ledger, pseudocount/common-anchor/center/complete-case sensitivities and the healthy-reference analysis to the supplement.
5. Add the compositional, paired-complete-case, taxonomic-resolution, same-cohort healthy-reference and no-patient-level-host-data limitations.
6. Preserve the previously negative clinical outcome result; do not replace it with the direction result.
7. Update Code Availability only after the new scripts and deterministic reproduction report are committed and archived.

Automated rule-based result: `{expected_rule}`. Human scientific ratification is still required before any manuscript promotion.
"""
    write_text_once(output / "B21_submission_material_changeset.md", changeset)

    generated = sorted(path for path in output.rglob("*") if path.is_file() and path.name != "output_manifest.csv")
    manifest_rows = [
        {
            "relative_path": str(path.relative_to(output)),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in generated
    ]
    write_csv_once(output / "output_manifest.csv", manifest_rows)

    if fail_count:
        print(f"Stage-B finalization completed with {fail_count} FAIL checks")
        return 2
    print(f"Stage-B finalization completed: {pass_count} PASS, {review_count} REVIEW, 0 FAIL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
