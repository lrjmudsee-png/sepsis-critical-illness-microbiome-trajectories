# Direction upgrade — A-stage feasibility audit

This folder contains the reproducible, read-only A-stage audit for the planned
clinical/mechanistic deepening of the sepsis manuscript.

## Scope

- retrieves and CRC-checks the official PRJNA851469 Supplementary Tables 2–17;
- inventories the four frozen analysis objects;
- audits strict sample/patient/time joins and frozen pair eligibility;
- applies the predeclared P/C Family dictionary for coverage and zero audits only;
- assesses healthy-reference and patient-level host-data feasibility;
- records overlap with existing Bray/Aitchison/EII and clinical analyses;
- does **not** calculate the proposed new effects or p values;
- does **not** rerun DADA2, FASTQ preprocessing, or Step98D;
- does **not** modify raw data, frozen results, or `V2-SUBMISSION`.

## Run

From PowerShell:

```powershell
& 'E:\sepsis_project\code\03_data_processing\99_DIRECTION_UPGRADE\RUN_A_FEASIBILITY.ps1'
```

For an independent reproduction, pass a fresh directory:

```powershell
& 'E:\sepsis_project\code\03_data_processing\99_DIRECTION_UPGRADE\RUN_A_FEASIBILITY.ps1' `
  -OutputDir 'E:\sepsis_project\results\V2_UPGRADE_20260913_DIRECTION\A_FEASIBILITY_REPRO'
```

The runner refuses to write into a non-empty output directory. This prevents
silent replacement of an earlier audit.

## Stage B direction analysis

Stage B implements the frozen direction-of-change secondary analysis without
rerunning sequence processing:

- `B00_FIXED_SECONDARY_ANALYSIS_SPECIFICATION.md`: time-stamped fixed rules.
- `99B0_initialize.py`: input and submission-tree immutability guards.
- `99B1_direction_analysis.R`: patient balance, meta-analysis, sensitivities,
  healthy-reference analysis, and diagnostics.
- `99B2_finalize.py`: independent checks, governance report, Chinese summary,
  and a manuscript integration candidate.
- `99B3_compare_reproduction.py`: byte-level deterministic rerun comparison.
- `RUN_B_DIRECTION.ps1`: fresh-output runner that retains failed runs.

Example:

```powershell
& 'E:\sepsis_project\code\03_data_processing\99_DIRECTION_UPGRADE\RUN_B_DIRECTION.ps1' `
  -OutputDirectory 'E:\sepsis_project\results\V2_UPGRADE_20260913_DIRECTION\B_DIRECTION_ANALYSIS'
```

## Main decision files

- `A12_feasibility_decision.md`
- `A13_questions_for_scientific_review.md`
- `A14_quality_control_summary.md`
- `input_hashes_before_after.csv`
- `output_manifest.csv`
