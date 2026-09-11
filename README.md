# Longitudinal gut-microbiome trajectories in sepsis and critical illness

Code-availability repository for the study of longitudinal gut-microbiome ecological displacement and heterogeneous taxonomic trajectories across sepsis and critical-illness cohorts.

## Repository contents

- `analysis/`: 342 source scripts (196 R, 139 PowerShell, and 7 Python) preserving the complete analysis history through the PRJNA1125274 validation and governance audit.
- `docs/RUN_ORDER.csv`: ordered map from metadata recovery and sequence processing to statistics and manuscript assembly.
- `docs/DATA_ACCESS.md`: public cohort accessions and local checkpoint definitions.
- `docs/ENVIRONMENT.md`: verified software and package versions.
- `evidence/`: records from the independently isolated terminal reproduction run performed on 2026-09-04, plus the downstream PRJNA1125274 validation audit completed on 2026-09-11.
- `tools/V2_REPRODUCE_TERMINAL_CHAIN.ps1`: runner for the 11-step frozen-input-to-manuscript terminal chain.

Raw sequencing reads, patient-level data, analysis workspaces, generated result trees, caches, and local directory links are not included in this repository.

## Public data

The study uses publicly available cohorts from NCBI SRA/BioProject, ENA, and NGDC-GSA:

`PRJNA691455`, `PRJNA516701`, `PRJNA851469`, `PRJNA578267`, `PRJNA430161`, `PRJNA1166732`, `PRJNA978257`, `PRJNA1010969`, `PRJEB82425`, and `CRA002354`.

Use the frozen sample/run inclusion manifests rather than downloading every run under each project accession. See `docs/DATA_ACCESS.md`.

## Verified terminal reproduction

The frozen-input terminal chain was rerun in an isolated result root on Windows 11 with R 4.4.0.

- 11/11 analysis and manuscript-integration steps exited successfully.
- 132/132 generated CSV files matched the frozen results after run-root normalization; 112 were identical without normalization.
- 10/10 PDFs were identical after rendering and SHA256 comparison.
- 8/8 PNG files were SHA256-identical.
- 42/42 text and completion-marker files matched after normalizing paths and timestamps.

See `evidence/REPRODUCTION_REPORT_20260904.md` and the machine-readable comparison tables in `evidence/`.

## PRJNA1125274 downstream validation

The Step98 series adds the external longitudinal cohort, patient-metadata reconciliation, frozen-result governance, non-target-feature and ASV-threshold sensitivities, complete-case selection audit, statistical-assumption checks, and reproducibility verification. The final report and machine-readable summaries are in `evidence/PRJNA1125274_20260911/`.

These downstream checks did not rerun FASTQ preprocessing or DADA2. The author-adjudicated 132-person clinical roster and clinical covariates needed for a definitive 132-versus-134 resolution and clinical attrition/confounding analysis are not available in this repository.

## Reproduction levels

### 1. Frozen inputs to manuscript outputs

This is the tested reviewer-facing route. It starts from frozen metadata and analysis objects and runs the following terminal chain:

1. trajectory ecology robustness;
2. longitudinal evidence freeze;
3. external validation evidence freeze;
4. infection-source evidence freeze;
5. cross-cohort taxonomic reproducibility;
6. global evidence synthesis;
7. manuscript source-pack assembly;
8. manuscript architecture and main tables;
9. Table 1 population correction and Results assembly;
10. final Results/Table 1 sentence audit;
11. supplementary tables and Discussion assembly.

The PowerShell runner expects a prepared Windows run root containing `_scripts/`, `data/`, `metadata/`, `processed_data/`, and the required frozen upstream `results/` directories.

```powershell
powershell -ExecutionPolicy Bypass -File tools/V2_REPRODUCE_TERMINAL_CHAIN.ps1 -RunRoot E:\path\to\prepared_run_root
```

### 2. Raw public reads to manuscript outputs

The complete historical workflow is documented in `docs/RUN_ORDER.csv`. It includes metadata reconstruction, DADA2 processing, the VSEARCH 97% OTU route for CRA002354, SILVA 138.2 taxonomy assignment, analysis-object freezing, longitudinal models, robustness analyses, and manuscript integration.

This full route has not been rerun after packaging and remains partly Windows-path-bound. Several historical scripts also install packages automatically. A future archival release should parameterize the remaining paths and lock the environment with `renv` or a container.

## Important provenance note

The repository preserves superseded and repair scripts because they document how the final frozen workflow was reached. Use `docs/RUN_ORDER.csv`, the `FINAL`, `FIXED`, and freeze-labelled scripts, and the 11-step terminal chain above when identifying authoritative outputs.

## Code availability statement

Analysis code, run order, environment information, public data accessions, and machine-readable terminal reproduction checks for this study are available at <https://github.com/lrjmudsee-png/sepsis-critical-illness-microbiome-trajectories>. Raw sequence data remain available from their originating public archives under the accessions listed above.
